rm(list=ls())

library(likelihood)
library(ggplot2)
library(tidyverse)
library(reshape2)
library(forcats)
library(broom)

load("./data/Data.Rdata")


targets <- targets %>%
      mutate(DG = rowMeans(dplyr::select(.,DG_16:DG_14),na.rm=TRUE),
             D = DB_t - (DG_16 + DG_15 + DG_14)) %>%
      dplyr::select(ID_plot:Y_UTM,Age, DG, D,prec, temp)

neighbours <- neighbours %>%
      mutate(Species2 = fct_lump(Species, n=4))

species_ne <- neighbours %>%
      filter(DBH_n >0) %>%
      group_by(Species2) %>%
      nest()


#  DBH as a function of diameter10 ----------------------------------------

dbh_d10_model <- function(df) {
      lm(DBH_n~ DB_n, data = df)
}

dbh_d10_models <- species_ne %>%
      mutate(mod = map(data,dbh_d10_model ),
             glance = mod %>% map(broom::glance),
             #AIC = glance %>% map_dbl("AIC"),
             R2 = glance %>% map_dbl("r.squared"),
             augment = mod %>% map(augment),
             tidy = mod %>% map(tidy))

dbh_d10_model_coefs <- dbh_d10_models %>%
      unnest(tidy) %>%
      dplyr::select(-(std.error: p.value)) %>%
      spread(term, estimate)

dbh_d10_model_predictions <- dbh_d10_models %>%
      unnest(augment)

ggplot() +
      geom_point(data= filter(neighbours, DBH_n >0) ,aes(DB_n, DBH_n), color="dark grey", alpha =0.5) +    
      geom_line(data=dbh_d10_model_predictions, aes(DB_n, .fitted), color="red", size=2) +
      facet_wrap(~Species2) +
      theme_bw()


# Height as a function of d10 ---------------------------------------------

height_d10_model <- function(df) {
      nls(Height_n ~ 0.1 + 30*(1 - exp(-B*DB_n)), data = df, 
          start = list( B=0.05), nls.control(maxiter=100))
}
height_d10_models <- species_ne %>%
      mutate(mod = map(data,height_d10_model ),
             glance = mod %>% map(broom::glance),
             #AIC = glance %>% map_dbl("AIC"),
             #R2 = glance %>% map_dbl("r.squared"),
             augment = mod %>% map(augment),
             tidy = mod %>% map(tidy))

height_d10_model_coefs <- height_d10_models %>%
      unnest(tidy) %>%
      dplyr::select(-(std.error: p.value)) %>%
      spread(term, estimate)

height_d10_model_predictions <- height_d10_models %>%
      unnest(augment)

ggplot() +
      geom_point(data= filter(neighbours, DB_n >0) ,aes(DB_n, Height_n), color="dark grey", alpha =0.5) +    
      geom_line(data=height_d10_model_predictions, aes(DB_n, .fitted), color="red", size=2) +
      facet_wrap(~Species2) +
      theme_bw()
