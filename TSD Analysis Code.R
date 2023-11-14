### Set up ####

# Load in packages
library(readr)
library(MCMCglmm)
library(tidyverse)
library(MuMIn)

# Load in sex ratio data
sex_ratio_full <- read_csv("sex_ratio_data.csv")
sex_ratio <- filter(sex_ratio_full, Life.Stage %in% c("birth.sex.ratio", "adult.sex.ratio"))

# Load in phylogeny data
phylo=read.tree("final_tree.phy") 
inv.phylo = inverseA(phylo, nodes="ALL", scale=TRUE)


### Creating a binary data set #####
# The following code creates a function that transforms the sex ratio data sets into a binary format to use in the MCMCglmm models
# Each data is rewritten from a sex ratio record to individuals from that record (male = 1, female = 0) 

binary <- function(data){
  
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
df_binary <- binary(sex_ratio) # to use in series 1 of models of just birth vs adults
df_binary_full <- binary(sex_ratio_full)# to use in  model 2 of best of series of model 1

### Models ####
##### Model 1.1 Single Mechanism Hypothesis
# prior
prior_SDM <-  list(G = list(G1 = list(V =1, nu = 0.02), G2 = list(V = 1, nu = 0.02)), R = list (V = 1, fix = 1))

#model
model1.1<- MCMCglmm(Sex ~ Life.Stage * sex.determination , random=~sex.determination:species +
                      species:population.ID, rcov=~sex.determination:units, family="categorical", 
                    ginverse=list(phylo=inv.phylo$Ainv), prior=prior_SDM, data=df_binary,
                    nitt=50000, burnin=300, thin=50, verbose=F)


##### Model 1.2 Evolutionary History Hypothesis
# prior
prior_phylo <- list(G = list(G1 = list(V =1, nu = 0.02), G2 = list(V = 1, nu = 0.02)), R = list (V = 1, fix = 1))

# model
model1.2 <- MCMCglmm(Sex ~ Life.Stage * SuperTaxa * sex.determination , random= ~SuperTaxa:species +
                       species:population.ID, rcov=~SuperTaxa:units, family="categorical", 
                     ginverse=list(phylo=inv.phylo$Ainv), prior=prior_phylo, data=df_binary,
                     nitt=50000, burnin=300, thin=50, verbose=F)


##### Model 1.3 Life History Hypothesis
#prior
prior_LH <- list(G = list(G1 = list(V =1, nu = 0.02), G2 = list(V = 1, nu = 0.02)), R = list (V = 1, fix = 1)) 

# model
model1.3 <- MCMCglmm(Sex ~ Life.Stage * LifeHistory * sex.determination , random=~LifeHistory:species +
                       species:population.ID, rcov=~LifeHistory:units, family="categorical", 
                     ginverse=list(phylo=inv.phylo$Ainv), prior=prior_LH, data=df_binary,
                     nitt=50000, burnin=300, thin=50, verbose=F) 


##### Model 1.4 Sex Determining Mechanism Hypothesis
# prior
prior_TSDtype<- list(G = list(G1 = list(V =1, nu = 0.02), G2=list(V= 1,nu=0.02)), R = list (V = 1, fix = 1)) 

# model
model1.4 <-  MCMCglmm(Sex ~ Life.Stage * SDM.Type, random=~SDM.Type:species + species:population.ID,
                      rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv),
                      prior=prior_TSDtype, data=df_binary,nitt=50000, burnin=300, thin=50, verbose=F)

#### Identifying and running the best model ####
# Best fit model with DIC
output <- model.sel(model1.1, model1.2, model1.3, model1.4, rank = "DIC")
output # model1.4 is best

# running model 4 with juvenile data
model2.1 <- MCMCglmm(Sex ~ Life.Stage * SDM.Type, random=~SDM.Type:species + species:population.ID,
                     rcov=~SDM.Type:units, family="categorical", ginverse=list(phylo=inv.phylo$Ainv),
                     prior=prior_TSDtype, data=df_binary_full,nitt=500, burnin=300, thin=5, verbose=F)