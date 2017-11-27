rm(list=ls())

library(likelihood)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(forcats)

load("./data/Data.Rdata")
load("./Data/Climate_Data.Rdata")

targets <- cbind(targets, clima)
targets <- targets %>%
      mutate(DG = rowMeans(dplyr::select(.,DG_16:DG_14),na.rm=TRUE),
             D = DB_t - (DG_16 + DG_15 + DG_14)) %>%
      dplyr::select(ID_plot:Y_UTM,Age, DG, D,prec, temp)

neighbours <- neighbours %>%
      mutate(Species2 = fct_lump(Species, n=4))

# Create the datasets for computing neighborhood indices ------------------

neighbours$vecinos <- sequence(tabulate(neighbours$ID_pine+1))

dbhs <- dcast(data=neighbours, ID_pine ~ vecinos, value.var= "DB_n")[,-1]
distances <- dcast(data=neighbours, ID_pine ~ vecinos, value.var= "Dist")[,-1]
species <- dcast(data=neighbours, ID_pine ~ vecinos, value.var= "Species2")[,-1]

using.targets <- complete.cases(targets)
wk.targets <- targets[using.targets,]
wk.dbhs <- dbhs[using.targets,]
wk.distances <- distances[using.targets,]
wk.species <- species[using.targets,]

nci_function <- function(j) {
      
      all.spp <- levels(neighbours$Species2)
      neigh.lambda.indexes <- apply(wk.species, c(1,2),
                                    function(x) {
                                          if (is.na(x)) {
                                                NA
                                          } else {
                                                which(all.spp == x) 
                                          }})
model_functions= list (
      ## Null model
      
      
      Model.01 <- function (PotGrowth) {
            PotGrowth},
      
      ## Size effect
      Model.02 <- function (PotGrowth, sizeX0, sizeXb) {
            size.effect <- exp(-0.5*((log(wk.targets$D/sizeX0)/sizeXb)^2))
            PotGrowth * size.effect
      },
      
      ## Size effect + competition effect
      Model.03 <- function(PotGrowth, sizeX0, sizeXb, alpha, beta, C,D) {
            size.effect <- exp(-0.5*((log(wk.targets$D/sizeX0)/sizeXb)^2))
            nci <- rowSums((wk.dbhs^alpha)/(wk.distances^beta), na.rm=T)
            competition.effect <- exp(-(C/100) * nci^D)
            PotGrowth * size.effect * competition.effect 
      },
      
      ## Size effect + competition effect (with tree size)
      Model.04 <- function(PotGrowth, sizeX0, sizeXb, alpha, beta,
                           Cprim, D,gamma) {
            size.effect <- exp(-0.5*((log(wk.targets$D/sizeX0)/sizeXb)^2))
            nci <- rowSums((wk.dbhs^alpha)/(wk.distances^beta), na.rm=T)
            C <- Cprim/100 * wk.targets$D^gamma 
            competition.effect <- exp(-(C) * (nci)^D)
            PotGrowth * size.effect * competition.effect 
      },

        ## Size effect + competition effect (with lambda)
      Model.05 <- function(PotGrowth, sizeX0, sizeXb, lambdaVec, 
                           alpha, beta, Cprim, D, gamma) {
            size.effect <- exp(-0.5*((log(wk.targets$D/sizeX0)/sizeXb)^2))
            lambda.vals <- lambdaVec[neigh.lambda.indexes]
            dim(lambda.vals) <- dim(neigh.lambda.indexes)
            nci <- rowSums(lambda.vals *(wk.dbhs ^ alpha)/(wk.distances ^ beta), na.rm=T)
            C <- Cprim/100 * wk.targets$D^gamma
            competition.effect <- exp(-(C) * (nci)^D)
            PotGrowth * size.effect * competition.effect
      },
      
      ## Size effect + competition effect (with lambda) + climate effect
      Model.06 <- function(PotGrowth, sizeX0, sizeXb, lambdaVec, 
                           alpha, beta, Cprim, D, gamma, 
                           tempX0, tempXb, precX0, precXb) {
            size.effect <- exp(-0.5*((log(wk.targets$D/sizeX0)/sizeXb)^2))
            lambda.vals <- lambdaVec[neigh.lambda.indexes]
            dim(lambda.vals) <- dim(neigh.lambda.indexes)
            nci <- rowSums(lambda.vals *(wk.dbhs ^ alpha)/(wk.distances ^ beta), na.rm=T)
            C <- Cprim/100 * wk.targets$D^gamma
            competition.effect <- exp(-(C) * (nci)^D)
            temp.effect <- exp(-0.5*((wk.targets$temp-tempX0)/tempXb)^2)
            prec.effect <- exp(-0.5*((wk.targets$prec-precX0)/precXb)^2)
            
            PotGrowth * size.effect * competition.effect * temp.effect * prec.effect
      },
      
      Model.07 <- function(PotGrowth, sizeX0, sizeXb, alpha, beta, C,D,
                              tempX0, tempXb, precX0, precXb) {
            size.effect <- exp(-0.5*((log(wk.targets$D/sizeX0)/sizeXb)^2))
            nci <- rowSums((wk.dbhs^alpha)/(wk.distances^beta), na.rm=T)
            competition.effect <- exp(-(C/100) * nci^D)
            temp.effect <- exp(-0.5*((wk.targets$temp-tempX0)/tempXb)^2)
            prec.effect <- exp(-0.5*((wk.targets$prec-precX0)/precXb)^2)
            PotGrowth * size.effect * competition.effect * temp.effect * prec.effect
      }
      )


      # Get parameters ---------------------------------------------------------
            # Initial parameters ------------------------------------------------------
            par_model <- list(par_Model.01 = list(PotGrowth = 0.5,
                                      sd = 2),
                  
                              par_Model.02 = list(PotGrowth = 0.5, 
                                      sizeX0= 5, 
                                      sizeXb= 1,
                                      sd = 2),
                  
                              par_Model.03 = list(PotGrowth = 0.5, 
                                      sizeX0= 5, 
                                      sizeXb= 1,
                                      alpha = 2,
                                      beta = 1,
                                      C = 0.0001,
                                      D =1,
                                      sd = 2),
                  
                              par_Model.04 = list(PotGrowth = 0.5, 
                                      sizeX0= 5, 
                                      sizeXb= 1,
                                      alpha = 2,
                                      beta = 1,
                                      Cprim = 0.001,
                                      D =1,
                                      gamma = 0.5,
                                      sd = 2),
                              
                              par_Model.05 = list(PotGrowth = 0.5, 
                                                  sizeX0= 5, 
                                                  sizeXb= 1,
                                                  alpha = 2,
                                                  beta = 1,
                                                  lambdaVec=rep(0.5,5),
                                                  Cprim = 0.001,
                                                  D =1,
                                                  gamma = 0.5,
                                                  sd = 2),
                              par_Model.06 = list(PotGrowth = 0.5, 
                                                  sizeX0= 5, 
                                                  sizeXb= 1,
                                                  alpha = 2,
                                                  beta = 1,
                                                  lambdaVec=rep(0.5,5),
                                                  Cprim = 0.001,
                                                  D =1,
                                                  gamma = 0.5,
                                                  tempX0 = 12,
                                                  tempXb = 5,
                                                  precX0 = 600,
                                                  precXb = 100,
                                                  sd = 2),
                              par_Model.07 = list(PotGrowth = 0.5, 
                                                  sizeX0= 5, 
                                                  sizeXb= 1,
                                                  alpha = 2,
                                                  beta = 1,
                                                  C = 0.0001,
                                                  D =1,
                                                  tempX0 = 12,
                                                  tempXb = 5,
                                                  precX0 = 600,
                                                  precXb = 100,
                                                  sd = 2)
                              )

            # Low parameters ------------------------------------------------------
                  par_low <- list(par_lo_Model.01 = list(PotGrowth = 0.1,
                                       sd = 0.0001),
                
                              par_lo_Model.02 = list(PotGrowth = 0.1,
                                       sizeX0= 1,
                                       sizeXb= 0.5,
                                       sd = 0.0001),
                
                              par_lo_Model.03 = list(PotGrowth = 0.1,
                                       sizeX0= 1,
                                       sizeXb= 0.5,
                                       alpha = 0.5,
                                       beta = 0.5,
                                       C = 1e-10,
                                       D =0.05,
                                       sd = 0.0001),
                
                              par_lo_Model.04 = list(PotGrowth = 0.1,
                                       sizeX0= 1,
                                       sizeXb= 0.5,
                                       alpha = 0.5,
                                       beta = 0.5,
                                       Cprim =1e-10,
                                       D =0.05,
                                       gamma = -2,
                                       sd = 0.0001),
                              
                              par_lo_Model.05 = list(PotGrowth = 0.1, 
                                                  sizeX0= 1, 
                                                  sizeXb= 0.5,
                                                  alpha = 0.5,
                                                  beta = 0.5,
                                                  lambdaVec=rep(0,5),
                                                  Cprim = 1e-10,
                                                  D =0.05,
                                                  gamma = -2,
                                                  sd = 2),
                              
                              par_lo_Model.06 = list(PotGrowth = 0.1, 
                                                     sizeX0= 1, 
                                                     sizeXb= 0.5,
                                                     alpha = 0.5,
                                                     beta = 0.5,
                                                     lambdaVec=rep(0,5),
                                                     Cprim = 1e-10,
                                                     D =0.05,
                                                     gamma = -2,
                                                     tempX0 = 1,
                                                     tempXb = 0.01,
                                                     precX0 = 1,
                                                     precXb = 0.01,
                                                     sd = 2),
                              par_lo_Model.07 = list(PotGrowth = 0.1,
                                                     sizeX0= 1,
                                                     sizeXb= 0.5,
                                                     alpha = 0.5,
                                                     beta = 0.5,
                                                     C = 1e-10,
                                                     D =0.05,
                                                     tempX0 = 1,
                                                     tempXb = 0.01,
                                                     precX0 = 1,
                                                     precXb = 0.01,
                                                     sd = 0.0001)
                              )




            # High parameters ------------------------------------------------------
                  par_high <- list(par_hi_Model.01 = list(PotGrowth = 1,
                                        sd = 10),
                 
                                    par_hi_Model.02 = list(PotGrowth = 1, 
                                        sizeX0= 10, 
                                        sizeXb= 2,
                                        sd = 10),
                 
                                    par_hi_Model.03 = list(PotGrowth = 1, 
                                        sizeX0= 10, 
                                        sizeXb= 2,
                                        alpha = 2.5,
                                        beta = 2.5,
                                        C = 1000,
                                        D = 3,
                                        sd = 10),
                 
                                    par_hi_Model.04 = list(PotGrowth = 1, 
                                        sizeX0= 10, 
                                        sizeXb= 2,
                                        alpha = 3,
                                        beta = 3,
                                        Cprim = 1000,
                                        D = 3,
                                        gamma = 4,
                                        sd = 10),
                                   
                                   par_hi_Model.05 = list(PotGrowth = 1, 
                                                       sizeX0= 10, 
                                                       sizeXb= 2,
                                                       alpha = 3,
                                                       beta = 3,
                                                       lambdaVec=rep(1,5),
                                                       Cprim = 1000,
                                                       D =3,
                                                       gamma = 4,
                                                       sd = 2),
                                   
                                   par_hi_Model.06 = list(PotGrowth = 1, 
                                                          sizeX0= 10, 
                                                          sizeXb= 2,
                                                          alpha = 3,
                                                          beta = 3,
                                                          lambdaVec=rep(1,5),
                                                          Cprim = 1000,
                                                          D =3,
                                                          gamma = 4,
                                                          tempX0 = 20,
                                                          tempXb = 100,
                                                          precX0 = 10,
                                                          precXb = 1000,
                                                          sd = 2),
                                   par_hi_Model.07 = list(PotGrowth = 1, 
                                                          sizeX0= 10, 
                                                          sizeXb= 2,
                                                          alpha = 2.5,
                                                          beta = 2.5,
                                                          C = 1000,
                                                          D = 3,
                                                          tempX0 = 20,
                                                          tempXb = 100,
                                                          precX0 = 400,
                                                          precXb = 1000,
                                                          sd = 10)
                                   )


# Annealing algorithm -----------------------------------------------------
# Define the model to use
modelname<- c('Model.01','Model.02','Model.03','Model.04', 
              'Model.05','Model.06', 'Model.07')

model.flo = model_functions[[j]]

# Parameter list (for effect of stage)
par.flo = par_model[[j]]
par_hi.flo = par_high[[j]]
par_lo.flo  = par_low[[j]]

#  set up the list of variables to use in annealing
pdf.flo = dnorm
dep.var = "DG"
var.flo<- list(mean = "predicted", x = "DG", log=T)

#  Call the annealing algorithm, specifying the proper model
results=anneal(model = model.flo, var = var.flo,
               source_data =wk.targets,dep_var="DG",
               par = par.flo, par_lo=par_lo.flo, par_hi=par_hi.flo,
               pdf=pdf.flo, hessian = FALSE, max_iter=35000)

return(results)
# write_results(results,paste("./Results/NCI_",modelname[j],".txt", sep=""))

}

results_NCImodels <- map(1:7, nci_function)
modelname<- c('Model.01','Model.02','Model.03','Model.04', 
              'Model.05','Model.06', 'Model.07')
names(results_NCImodels) <- modelname [1:7]

modelsAIC <- map(results_NCImodels, "aic_corr")
best_pars <- map_df(results_NCImodels[7], "best_pars")
best_pars2 <- map_df(results_NCImodels[3], "best_pars")

save(best_pars, file= "./data/best_pars.Rdata")



# Check results of NCI models ---------------------------------------------

load("./data/best_pars.Rdata")

seqTemp = seq(6,15,by=0.2)
seqPrec = seq(0,3000,by=1)
seqDiam <- seq (0.1, 50,0.01)
seqDist <- seq(0.01, 25, 0.1)
seqNCI <- seq (0, 2000, 1)

## For adults...

pini_temp.effect <- exp(-0.5*((seqTemp-0.47)/11.59)^2)
pini_prec.effect <- exp(-0.5*((seqPrec-2569.78)/2131.21)^2)
pini_size.effect <- exp(-0.5*((log(seqDiam/18.68)/1.05)^2))
pini_competition.effect <- exp(-(0.01609) * 20^(-1.23) * seqNCI)
pini_efecto_diam <- seqDiam^1.96
pini_efecto_dist <- 1/seqDist^1.11

pisy_temp.effect <- exp(-0.5*((seqTemp-1.09)/11.75)^2)
pisy_prec.effect <- exp(-0.5*((seqPrec-2386.2)/2660.83)^2)
pisy_size.effect <- exp(-0.5*((log(seqDiam/19.92)/1.11)^2))
pisy_competition.effect <- exp(-(0.03506) * 20^(-1.11) * seqNCI)
pisy_efecto_diam <- seqDiam^1.82
pisy_efecto_dist <- 1/seqDist^1.45

plot(seqTemp, pisy_temp.effect, type ="l", col = "dark orange", ylim=c(0,1))
lines(seqTemp, pini_temp.effect, type ="l", col= "dark red")

plot(seqPrec, pisy_prec.effect, type ="l", col = "dark orange",  ylim=c(0,1))
lines(seqPrec, pini_prec.effect, type ="l", col= "dark red")

plot(seqDiam, pisy_size.effect, type ="l", col = "dark orange")
lines(seqDiam, pini_size.effect, type ="l", col= "dark red")

plot(seqDiam, pisy_efecto_diam, type ="l", col = "dark orange")
lines(seqDiam, pini_efecto_diam, type ="l", col= "dark red")

plot(seqDist, pisy_efecto_dist, type ="l", col = "dark orange")
lines(seqDist, pini_efecto_dist, type ="l", col= "dark red")

plot(seqNCI, pisy_competition.effect, type ="l", col = "dark orange")
lines(seqNCI, pini_competition.effect, type ="l", col= "dark red")



# Michaelis-Menten --------------------------------------------------------

DBH <- 6
GLI <- seq(0.01,100,0.02)
MM <- function(GLI,DBH,A,S, D){
      DBH^D * (A*GLI)/((A/S)+GLI)
}

pisy <- MM(GLI,DBH,0.804,0.008,1)
pini <- MM(GLI,DBH, 0.30,0.029,1)
abal <- MM(GLI,DBH,0.39,0.044,0.447)
piun <- MM(GLI,DBH,1.304,0.0074,1.041)

plot(GLI, pisy, type="l", col = "dark orange", lwd=2)
lines(GLI, piun, type="l", col = "black", lwd=2)
lines(GLI, abal, type="l", col = "dark green", lwd=2)
lines(GLI, pini, type="l", col = "dark red", lwd=2)

