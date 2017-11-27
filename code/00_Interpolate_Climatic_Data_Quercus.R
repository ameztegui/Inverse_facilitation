rm(list=ls())

library(meteoland)
library(tidyverse)
library(sp)
library(raster)

# Import data, including climatic data
qc_coords <- read_csv2("./data/quercus_plots_coords.csv")
load("./data/climate/SMC_interpolator.Rda")

# Import cartography and generate elevation and aspect
solsones <- shapefile("S:/Cartografia/Referencia/Limites_Admin/Catalunya/Comarcas/SOL_00LIMIT.shp")
iiff98 <- shapefile("S:/Cartografia/Tematica/INCENDIO_98/Limits98_Diss.shp")

elevation <- raster("S:/Cartografia/Tematica/INCENDIO_98/incend_mdt5_WGS84_HU31.tif")
hill <- raster("S:/Cartografia/Tematica/INCENDIO_98/incend_hillshade_WGS84_HU31.tif")
slope <- terrain(elevation, opt = "slope",unit = "degrees")
aspect <- terrain(elevation, opt = "aspect",unit = "degrees")


# Transform points in SpatialPointsDataFrame object (S4)
coordinates(qc_coords) <- c("X", "Y")
proj4string(qc_coords) <- CRS(" +proj=utm +zone=31 +ellps=WGS84 +datum=WGS84 +units=m +towgs84=0,0,0 " )

# Plot points

plot(hill,
     col=grey(1:100/100),  # create a color ramp of grey colors for hillshade
     legend=FALSE,         # no legend, we don't care about the grey of the hillshade
     axes=FALSE,          # makes for a cleaner plot, if the coordinates aren't necessary
     box =F)
plot(elevation,
     axes=FALSE,
     alpha=0.5,   # sets how transparent the object will be (0=transparent, 1=not transparent)
     add=T,col= terrain.colors(20))
# lines(iiff98, col = "orange", lwd = 2)
points(qc_coords, col = "dark red", pch = 20, cex = 6)



# Calculate meteo data ----------------------------------------------------


# Extract point values from rasters
points_elevation <- extract(elevation, qc_coords)
points_slope <- extract(slope, qc_coords)
points_aspect <- extract(aspect, qc_coords)

# Create spatialPointsTopography Object
points_topo <- SpatialPointsTopography(qc_coords, points_elevation, points_slope, points_aspect)

# Define dates of interest
hist_dates = seq.Date(as.Date("1998-01-01"),as.Date("2016-12-31"), by="day")

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

# Export results
temp_quercus <- write.table(yearly_temp@data, file = "./Data/qc_temp.txt", dec=",", row.names = F, sep = "\t")
prec_quercus <- write.table(yearly_prec@data, file = "./Data/qc_prec.txt", dec=",", row.names = F, sep = "\t")

