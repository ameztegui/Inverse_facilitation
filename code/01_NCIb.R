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
            temp.effect <- exp(-0.5*((wk.targets$temp-tempX0/tempXb)^2))
            prec.effect <- exp(-0.5*((wk.targets$prec-precX0/precXb)^2))
            
            PotGrowth * size.effect * competition.effect * temp.effect * prec.effect
      })

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
                                                     sd = 2)
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
                                                          sd = 2)
                                   )


# Annealing algorithm -----------------------------------------------------
# Define the model to use
modelname<- c('Model.01','Model.02','Model.03','Model.04', 
              'Model.05','Model.06' )
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
               pdf=pdf.flo, hessian = FALSE, max_iter=25000,
               min_change=0.1,min_drops=2.5)
# 
write_results(results,paste("./Results/NCI_",modelname[j],".txt", sep=""))

}

nci_function(1)
nci_function(2)
nci_function(3)
nci_function(4)
nci_function(5)
nci_function(6)


DB <- seq (0.1, 10,0.01)
size.effect2 <-  exp(-0.5*((log(DB/10)/1.088)^2))
size.effect3 <-  exp(-0.5*((log(DB/10)/1.12)^2))
size.effect4 <-  exp(-0.5*((log(DB/10)/1.1053)^2))

ggplot()+
       geom_line(aes(DB, size.effect2)) +
       geom_line(aes(DB, size.effect3), color="steelblue") +
      geom_line(aes(DB, size.effect4), color="red", size=1.5) +
      ylab("Size effect") +
      ylim(0,1) +
      theme_bw()

wk.targets$nci3 <- rowSums((wk.dbhs^2.435)/(wk.distances^0.56), na.rm=T)
wk.targets$nci4 <- rowSums((wk.dbhs^2.667)/(wk.distances^0.702), na.rm=T)

test <-wk.targets %>% 
      mutate(size.effect4 = exp(-0.5*((log(((DB)/10)/1.088))^2)),
             comp.effect3 = exp(-(36.603/100) * nci3^0.0828),
            comp.effect4 = exp(-(0.02937/100) * nci4^0.8064*DB^-0.274),
            predRG3 = 1*size.effect4*comp.effect3,
             predRG4 = 0.58*size.effect4*comp.effect4)

ggplot() +
      geom_point(data=test, aes(RG, predRG3)) +
     geom_point(data=test, aes(RG, predRG4), col="steel blue") +
      xlim(0,0.75) + ylim(0,0.75) +
      geom_abline(slope =1, intercept = 0, color="gray") +
      theme_bw()

hist(wk.targets$nci3)
hist(wk.targets$nci4)



nci3 <- seq (1,4000,10)
nci4 <- seq (1,8000,10)

competition.effect3 <- exp(-(36.603/100) * nci3^0.0828)
competition.effect4 <- exp(-(0.02937/100) * nci4^0.8064*3.56^-0.274)

ggplot()+
      geom_line(aes(nci3, competition.effect3)) +
      geom_line(aes(nci4, competition.effect4), color="red") +
      ylim(0,1)


ggplot() +
      geom_point(data=test, aes(nci4,comp.effect4))+
      geom_line(aes(nci4, competition.effect4), col="red", size =1.5) +
      xlab("NCI") + ylab("Competition effect") + ylim(0.5,1)+
      theme_bw()



