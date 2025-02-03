### Set up ####
# Load in packages
library(readr)
library(tidyverse)
library(MCMCglmm)
library(MuMIn)
library(ape)


### Load in sex ratio data
setwd("~/Desktop/Mariel/Paper/Revisions/July2024 - workstudy data/FINAL CODE + DATA - to publish") 
df <- read_csv("sex ratio data final.csv")

# conservative data with juvenile data
conservative_full <- filter(df, Conservative == "Y")

# conservative data without juvenile data
conservative <- df %>% 
  filter(Conservative == "Y",
         Life.Stage %in% c("birth.sex.ratio", "adult.sex.ratio"))

# sda (sexual dimorphic assumed data added) data with juvenile data
sda_full <- df

# sda without juvenile data
sda <- filter(df, Life.Stage %in% c("birth.sex.ratio", "adult.sex.ratio"))

### Load in phylogeny data ####
setwd("~/Desktop/Mariel/Paper/Revisions/July2024 - workstudy data")
## phylogeny for sda models
phylo=read.tree("ultimate phylogeny.phy") 
inv.phylo = inverseA(phylo, nodes="ALL", scale=TRUE)

## phylogeny for conservative models with juvenile data
phylo_conservative <- read.tree("conservative full phylogeny.phy")
inv.phylo.conservative <- inverseA(phylo_conservative, nodes = "ALL", scale = TRUE)

## phylogeny for conservative models without juvenile data
phylo_conservativeBA <- read.tree("conservative BA phylogeny.phy")
inv.phylo.conservativeBA <- inverseA(phylo_conservativeBA, nodes = "ALL", scale = TRUE)

### Create a binary dataset ####
binary_ultimate <- function(data){
  
  # filter data so it only contains values with sex ratio (double check)
  temp <- filter(data, !is.na(Sex.Ratio))
  
  # create phylo column that is identical to species column (need for model)
  temp2 <- mutate(temp, phylo = species)
  
  # get column for number of males and females, remove those w/o sample size
  temp3 <- temp2 %>% 
    mutate(Nmales = round(Sex.Ratio * N)) %>% 
    mutate(Nfemales = round((1-Sex.Ratio)*N))
  temp4 <- filter(temp3, !is.na(N))
  
  # create empty data frame to fill in
  binary_data <- data.frame()
  
  # for each row in data; convert to binary
  for(i in 1:nrow(temp4)){
    row_fill <- temp4[i, ] # get basic information
    row_fillM <- cbind(row_fill, Sex  = c(1)) # fill in with column added if male
    row_fillF <- cbind(row_fill, Sex = c(0)) # fill in with column added if female
    
    # Get number of each sex reported
    NMales <- temp4[i, "Nmales"] # number of males (rows needed)
    NFemales <- temp4[i, "Nfemales"] # number of females(rows needed)
    
    # Create new data frame
    df <- rbind(row_fillM, row_fillF) # get data frame with both ones
    
    # Contingency in case a population is reported as 100% male or female
    if(NMales == 0|NFemales == 0){
      # If males are 0
      if(NMales == 0){  # remove male row if no males in this population reported
        df <- df[-1, ]
        df_temporary <- as.data.frame(lapply(df, rep, NFemales))
      }
      
      # If Females are 0
      if(NFemales == 0){df <- df[-2,]
      df_temporary <- as.data.frame(lapply(df, rep, NMales))
      } # remove if no females in this population reported
    }
    # Create temporary data set with replicated number of males and females
    else{df_temporary <- as.data.frame(lapply(df, rep, c(NMales, NFemales)))}
    
    # Add to total data frame
    binary_data <- rbind(binary_data, df_temporary)
    
  }
  # fix population id column name
  names(binary_data)[names(binary_data) == '?..population.ID'] <- 'population.ID'
  
  # return data
  return(binary_data)
  
}

# Making data into binary format
## for conservative data without juvenile data
conservative_binary <- binary_ultimate(conservative) 
## for conservative data with juvenile data
conservative_binary_full <- binary_ultimate(conservative_full)

## for sda data without juvenile data
sda_binary <- binary_ultimate(sda)
## for sda data with juvenile data
sda_binary_full <- binary_ultimate(sda_full)

### Model 1: single mechanism hypothesis #######
prior_SDM <-  list(G = list(G1 = list(V =1, nu = 0.02), G2 = list(V = 1, nu = 0.02)), R = list (V = 1, fix = 1))

# conservative model 1
model_conservative1.1<- MCMCglmm(Sex ~ Life.Stage * sex.determination , random=~sex.determination:species + species:population.ID, rcov=~sex.determination:units, family="categorical", ginverse=list(phylo=inv.phylo.conservativeBA$Ainv), prior=prior_SDM, data=conservative_binary,nitt=50000, burnin=3000, thin=50, verbose=F)


# sda model 1
model_sda1.1<- MCMCglmm(Sex ~ Life.Stage * sex.determination , random=~sex.determination:species + species:population.ID, rcov=~sex.determination:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_SDM, data= sda_binary,nitt=50000, burnin=3000, thin=50, verbose=F)


### Model 2: evolutionary history model ####
prior_phyloLH <- list(G = list(G1 = list(V =1, nu = 0.02), G2 = list(V = 1, nu = 0.02)), R = list (V = 1, fix = 1))

#conservative model 2
model_conservative1.2 <- MCMCglmm(Sex ~ Life.Stage * SuperTaxa * sex.determination , random= ~SuperTaxa:species + species:population.ID, rcov=~SuperTaxa:units, family="categorical", ginverse=list(phylo=inv.phylo.conservativeBA$Ainv), prior=prior_phyloLH, data=conservative_binary, nitt=50000, burnin=3000, thin=50, verbose=F)


# sda model 2
model_sda1.2 <- MCMCglmm(Sex ~ Life.Stage * SuperTaxa * sex.determination , random= ~SuperTaxa:species + species:population.ID, rcov=~SuperTaxa:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_phyloLH, data=sda_binary, nitt=50000, burnin=3000, thin=50, verbose=F)


### Model 3: TSD type model ####
prior_TSDtype<- list(G = list(G1 = list(V =1, nu = 0.02), G2=list(V= 1,nu=0.02)), R = list (V = 1, fix = 1)) 

# conservative model 3
model_conservative1.3 <-  MCMCglmm(Sex ~ Life.Stage * SDM.Type, random=~SDM.Type:species + species:population.ID, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo.conservativeBA$Ainv), prior=prior_TSDtype, data=conservative_binary,nitt=50000, burnin=3000, thin=50, verbose=F)

# sda model 3
model_sda1.3 <-  MCMCglmm(Sex ~ Life.Stage * SDM.Type, random=~SDM.Type:species + species:population.ID, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_TSDtype, data=sda_binary,nitt=50000, burnin=3000, thin=50, verbose=F)


### Model 4: Taxonomic Model 
prior_TSDtype<- list(G = list(G1 = list(V =1, nu = 0.02), G2=list(V= 1,nu=0.02)), R = list (V = 1, fix = 1)) 

# conservative model 4
model_conservative1.4 <-  MCMCglmm(Sex ~ Life.Stage * model4, random=~model4:species + species:population.ID, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo.conservativeBA$Ainv), prior=prior_TSDtype, data=conservative_binary,nitt=50000, burnin=3000, thin=50, verbose=F)


# sda model 4
model_sda1.4 <-  MCMCglmm(Sex ~ Life.Stage * model4, random=~model4:species + species:population.ID, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_TSDtype, data=sda_binary,nitt=50000, burnin=3000, thin=50, verbose=F)



### Identifying Best Model ####
# conservative models
output_conservative <- model.sel(model_conservative1.1, model_conservative1.2, model_conservative1.3, model_conservative1.4, rank = "DIC")

output_conservative # best one is model1.4

# sda models
output_sda <- model.sel(model_sda1.1, model_sda1.2, model_sda1.3, model_sda1.4,rank = "DIC")
output_sda # best one is model1.4


### Best Model with Full Data ####

# conservative model 4 - with juvenile data
model4_conservative_full <-  MCMCglmm(Sex ~ Life.Stage * model4, random=~model4:species + species:population.ID, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo.conservative$Ainv), prior=prior_TSDtype, data=conservative_binary_full,nitt=50000, burnin=3000, thin=50, verbose=F)


# sda model 4 - with juvenile data
model4_sda_full<-  MCMCglmm(Sex ~ Life.Stage * model4, random=~model4:species + species:population.ID, rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv), prior=prior_TSDtype, data=sda_binary_full,nitt=50000, burnin=3000, thin=50, verbose=F)




