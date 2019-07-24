#!/usr/bin/r
# READ MODIS ATMOSPHERE DATA AND RUN FORTRAN CODE get_ed0 TO PRODUCE DAILY
# IRRADIANCE FILES AT DESIRED TEMPORAL RESOLUTION
# get_ed READS AND INTERPOLATES 5-DIMENSION LUT TABLES CREATED WITH SBDART BY
# SIMON BELANGER WITH DIMENSIONS: 83 WAVELENGTH, 19 SZA, 8 COT, 10 TO3, 7 SURF ALBEDO
# Marti Gali Tapias 2017-05-08

# ------------------------------------- EDIT -------------------------------------

## Choose Ice Camp or Green Edge cruise configuration (backwards search done only for cruise)

# # Ice camp configuration, only daily data searched
# stnlist <- read.table("~/Desktop/GreenEdge/Irradiance/input.noclim.IceCamp20152016.txt")
doydecr <- 0;

# # GreenEdge cruise configuration, backward search at all stations
# stnlist <- read.table("~/Desktop/GreenEdge/Irradiance/input.noclim.GE2016.txt")
# doydecr <- c(-2,-1,0); 

# Above (1) or below (0) water surface?
above = 1;

# What output frequency do we want?
hperiod = 1; # integer

## Where are files?
# mainpath <- "/Volumes/taku-njall/MODISA/L3BIN"
mainpath <- "~/Desktop/GreenEdge/Irradiance/Ed_RT_MODISA"; # for tests on MBP

# Where is the code that reads and interpolates the LUTs?
codepath <- "~/Desktop/GreenEdge/Irradiance/fortran_src_MGT"

# --------------------------------- END EDIT -------------------------------------

# Load libraries
library(rgdal)
library(gdalUtils)
library(sp)

# Choose above surface (Ed0plus LUT) or below surface (Ed0moins LUT)
if (above == 1){
  lutpath <- "/Users/martigalitapias/Documents/SBDART/Ed0plus_LUT_5nm_v2.dat"
  outpath <- "~/Desktop/GreenEdge/Irradiance/Ed0plus_MODISA_LUT_SurfAlb"
}else{
  lutpath <- "/Users/martigalitapias/Documents/SBDART/Ed0moins_LUT_5nm_v2.dat"
  outpath <- "~/Desktop/GreenEdge/Irradiance/Ed0moins_MODISA_LUT_SurfAlb"
}

# Define periods within day
hours <- 1:hperiod:23

# get_ed0 run in loop for stations, search days (backwards) at each station, and hours within day
# for (stn in 1:dim(stnlist)[1]){
for (stn in 1){
  
  # stn$LAT <- stnlist$V1[stn]
  # stn$LON <- stnlist$V2[stn]
  # stn$Y <- stnlist$V5[stn]
  # stn$DOY <- stnlist$V6[stn]
  
  # For tests
  stn$LAT <- 67.1
  stn$LON <- -57.1
  stn$Y <- 2003
  stn$DOY <- 181
  
  for (dd in 1:length(doydecr)){
    
    searchY <- stn$Y
    searchDOY <- stn$DOY + doydecr[dd]
    if (searchDOY < 1){
      # Case where changing day changes year
      searchY <- stn$Y - 1
      if (!(searchY %% 4) == 0){
        # Leap year
        searchDOY <- 366
      }else{
        # Normal year
        searchDOY <- 365
      }
    }
    
    # Preallocate output
    OUT <- matrix(NaN,83,length(hours))
    
    # Build outfile path
    outfilepath <- sprintf("%s/Ed0_%04i_%03i_%0.3f_%0.3f.txt", outpath, searchY, searchDOY, stn$LAT, stn$LON)
    
    # Build file path
    genericname <- sprintf("MYD08_D3.A%04i%03i.051.*", searchY, searchDOY)
    searchpath <- sprintf("%s/%04i/%03i/%s", mainpath, searchY, searchDOY, genericname)
    
    # Find the complete filename
    getfilepath <- paste("ls", searchpath, sep = " ")
    flist <- system(getfilepath, intern = TRUE)
    
    # Test if file found and either go ahead or put filename in list of files not found
    if (length(flist) == 0){
      
      write.table(stnlist[stn,], file = "files_not_found.txt", append = TRUE)
      
    }else{
      
      # Here we assume there is only one such file!
      filepath <- flist[1]
      
      # Read variables
      sds <- get_subdatasets(filepath)
      
      # Find pixel match and convert to right units
      
      
      # Define albedo as a function of sea ice cover
      stn$ALB <- 0.25
      
      # Call get_ed0 in hours loop
      system("rm TMP.txt")
      for (ih in hours){
        cmd <- sprintf("./get_ed0 %i %i %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f %0.2f >> TMP.txt",
                       above,searchDOY,ih,stn$LAT,stn$LON,stn.COT,stn.TO3,stn.CF,stn$ALB);
        system(cmd)
      }
      
      # Read output for daily data
      TMP <- read.table("TMP.txt", header = F, sep = " ")
      ncols <- dim(TMP)[2]
      TMP <- TMP[,(ncols-82):ncols]
      OUT <- t(TMP)
      
    }
  }  
}
