rm(list=ls())
library(meteoland)
library(tidyverse)
library(sp)
library(raster)

# Import data, including climatic data
load("./data/Data.Rdata")
load("./data/climate/SMC_interpolator.Rda")

# Import cartography and generate elevation and aspect
solsones <- shapefile("S:/Cartografia/Referencia/Limites_Admin/Catalunya/Comarcas/SOL_00LIMIT.shp")
elevation <- raster("S:/Cartografia/Referencia/MDT/MDT30CATb.tif")
slope <- terrain(elevation, opt = "slope",unit = "degrees")
aspect <- terrain(elevation, opt = "aspect",unit = "degrees")

# Transform points in SpatialPointsDataFrame object (S4)
points <- targets
coordinates(points) <- c("X_UTM", "Y_UTM")
proj4string(points) <- CRS(" +proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +towgs84=0,0,0 " )

# Plot points
plot(elevation)
lines(solsones)
points(points)

# Extract point values from rasters
points_elevation <- extract(elevation, points)
points_slope <- extract(slope, points)
points_aspect <- extract(aspect, points)

# Create spatialPointsTopography Object
points_topo <- SpatialPointsTopography(points, points_elevation, points_slope, points_aspect)

# Define dates of interest
hist_dates = seq.Date(as.Date("2008-01-01"),as.Date("2016-12-31"), by="day")

# Interpolate climatic data
mp = interpolationpoints(SMC_interpolator, points_topo, dates = hist_dates)

# Create monthly summaries
monthly_temp <- summarypoints(mp,var = "MeanTemperature", fun = mean, freq= "month")
monthly_mtemp <- summarypoints(mp,var = "MinTemperature", fun = mean, freq= "month")
monthly_Mtemp <- summarypoints(mp,var = "MaxTemperature", fun = mean, freq= "month")
monthly_prec <- summarypoints(mp,var = "Precipitation", fun = sum, freq= "month")

# Create yearly summaries
yearly_temp <- summarypoints(mp,var = "MeanTemperature", fun = mean, freq= "year")
yearly_prec <- summarypoints(mp,var = "Precipitation", fun = sum, freq= "year")

# Plot results
meteoplot(mp, index = 2, var = "MeanTemperature", freq = "year")
meteoplot(mp, index = 2, var = "Precipitation", freq = "year")

clima <- yearly_prec@data %>%
      mutate(prec= rowMeans(., na.rm=T),
             temp = rowMeans(yearly_temp@data, na.rm =T)) 
clima <- clima[, c("prec", "temp")]
clima$prec <- as.numeric(clima$prec)
clima$temp <- as.numeric(clima$temp)
save(clima, file ="./Data/Climate_Data.Rdata")
      
