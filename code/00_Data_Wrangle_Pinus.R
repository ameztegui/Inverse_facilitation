rm(list=ls())

library(reshape2)
library(tidyverse)
library(forcats)
library(broom)


# Load and wrangle the data -----------------------------------------------------------

targets <- read_tsv("./data/Pinos.txt")
competitors <- read_tsv("./data/Competidores.txt")
crecimientos <- read_tsv("./data/Crecimientos.txt")

##### Wrangle target dataset

      targets <- targets %>%
            mutate(Met = ((Met1 + Met2 + Met3 + Met4 + Met5)/5),      # average height growth
                   RelHG = Met/Height,
                   DB = DB/10,                 # in cm
                   DBH = DBH/10,               # in cm
                   eq_DBH = (0.41257*DB - 0.3467))   # allometry


##### Wrangle competitors dataset
      
      # Replace 0s by NAs
      competitors$eqDBH <- map_dbl(competitors$eqDBH, function(x){replace(x, x == 0, NA)})
      competitors$eqDB <- map_dbl(competitors$eqDB, function(x){replace(x, x == 0, NA)})
      
      # Parse "species" to factor
      competitors$Species <- factor(competitors$Species)
      
      # Correct aberrant data
      competitors$eqDB[competitors$ID_compet == 82701] <- 25
      
      # Change unities and lump species into three categories
      competitors <- competitors %>%
            mutate(eqDBH = eqDBH/10,                     # in cm
                   eqDB = eqDB/10,                       # in cm 
                   Species2 = fct_lump(Species, n=4))    # lump species into 5 categories
      
      competitors$Species2[is.na(competitors$Species2)] <- "Other"

      
##### Wrangle crecimientos dataset

      # Substitute 0s by NAs
      crecimientos$a <- map_dbl(crecimientos$a, function(x){replace(x, x == 0, NA)})
      crecimientos$b <- map_dbl(crecimientos$b, function(x){replace(x, x == 0, NA)})
      crecimientos$Mean_growth <- map_dbl(crecimientos$Mean_growth, function(x){replace(x, x == 0, NA)})
      crecimientos$SD_growth <- map_dbl(crecimientos$SD_growth, function(x){replace(x, x == 0, NA)})



# Compute missing diameter values -----------------------------------------

# We want to determine DBH based on DB, and viceversa. However, we will do this on a species-specific way,
#       determining the allometry between both variables for each species

# Firs, let's have a look at the DBH - DB relationship
ggplot(competitors) +
      geom_point(aes(eqDB, eqDBH, color=Species2))

# And the height- DB relationship
ggplot(competitors) +
      geom_point(aes(eqDB, Height, color=Species2))


# define model to apply (lineal model)
species_model <- function (df) {
      lm(eqDBH ~ eqDB, data = df)
}

# nest data into a list of dataframes (by species)
# and apply the models
models <- competitors %>%
      group_by (Species2) %>%
      nest() %>%
      mutate(mod= map(data, species_model),
             glance = mod %>% map(broom::glance),    # R2
             tidy = mod %>% map(tidy),               # coefficients
             augment = mod %>% map(augment))         # predictions 

# extract the coefficients
model_coefs <- models %>%
      unnest(tidy) %>%
      dplyr::select(-(std.error: p.value)) %>%
      spread(term, estimate)

competitors <- competitors %>%
      mutate(eqDBH = ifelse (is.na(eqDBH) & Species2 == "Pinus nigra", -1.5908+0.8204*eqDB,
                              ifelse (is.na(eqDBH) & Species2 == "Quercus faginea", -1.0052+0.7913*eqDB,
                              ifelse (is.na(eqDBH) & Species2 == "Quercus ilex", -0.4908+0.7991*eqDB,
                              ifelse (is.na(eqDBH) & Species2 == "Quercus pubescens", -0.77784+0.7385*eqDB,
                              ifelse (is.na(eqDBH) & Species2 == "Other", -0.3873+0.7603*eqDB,
                                                              eqDBH)))))) %>%
      mutate(eqDB = ifelse (is.na(eqDB) & Species2 == "Pinus nigra", (eqDBH + 1.5908) / 0.8204,
                              ifelse (is.na(eqDB) & Species2 == "Quercus faginea", (eqDBH +1.0052) /0.7913,
                              ifelse (is.na(eqDB) & Species2 == "Quercus ilex", (eqDBH + 0.4908) / 0.7991,
                              ifelse (is.na(eqDB) & Species2 == "Quercus pubescens", (eqDBH + 0.77784) / 0.7385,
                              ifelse (is.na(eqDB) & Species2 == "Other", (eqDBH+0.3873) /0.7603,
                                                              eqDB))))))



# Compute radial growth for analyses --------------------------------------

dendro <- crecimientos %>%
      select(ID_pine,Year, Mean_growth) %>%
      mutate(year1="y") %>%
      unite(yearN, year1, Year) %>%
      spread(yearN, Mean_growth)

# make growth to be in diameter cm per year )divide by 1000 and multiply by 2
dendro[-1] <- dendro[-1]/500 

dendro <- dendro %>%
      mutate(RG=(y_2016 + y_2015 + y_2014 + y_2013 + y_2012)/5)
targets <- left_join(targets, dendro) %>%
      mutate(DB_2015 = DB - y_2016,
             DB_2014 = DB_2015 - y_2015,
             DB_2013 = DB_2014 - y_2014,
             DB_2012 = DB_2013 - y_2013,
             DB_2011 = DB_2012 - y_2012)


# Create the datasets for computing neighborhood indices ------------------

competitors$neighbours <- sequence(tabulate(competitors$ID_pine+1))
dbhs <- dcast(data=competitors, ID_pine ~ neighbours, value.var= "eqDB")[,-1]
distances <- dcast(data=competitors, ID_pine ~ neighbours, value.var= "Dist_pine")[,-1]
species <- dcast(data=competitors, ID_pine ~ neighbours, value.var= "Species2")[,-1]



save(targets, dbhs, distances, species, file="./data/NCI.Rdata")
