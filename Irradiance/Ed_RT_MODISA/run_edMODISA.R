#!/usr/bin/r
# Read list of station dates (coordinates, year, julian day), edit call_ed.R and run it for each station
# Martí Galí Tàpias 20170509

# WORK IN PROGRESS!

# ------------------------------------- EDIT -------------------------------------

## Define DOY decrements to output Ed for previous days (comment either case 0 or 1)

# # Case 0: normal configuration, only daily data searched
# stnlist <- read.table("~/Desktop/GreenEdge/Irradiance/input.noclim.IceCamp20152016.txt")
# doydecr <- 0; 

# Case 1: GreenEdge cruise, backward search at all stations
stnlist <- read.table("~/Desktop/GreenEdge/Irradiance/input.noclim.GE2016.txt")
doydecr <- c(-2,-1,0); 


## Where are files?
mainpath <- "/Volumes/taku-njall/MODISA/L3BIN"


## Where is the LUT? Where do we output files? Uncomment 2 of 3

# # Case Ed0moins A: LUT on taku-njall (Thetas(19)*Ozone(8)*TauCld(8), 83 lambdas)
# lutpath <- "/Volumes/taku-njall/LUTS/Ed0moins_LUT.dat"
# outpath <- "~/Desktop/GreenEdge/Irradiance/Ed0moins_MODISA_noclim_LUT_A"

# # Case Ed0moins B: LUT given by Srikanth (Thetas(19)*Ozone(10)*TauCld(8), 83 lambdas)
# lutpath <- "/Users/martigalitapias/Desktop/GreenEdge/Irradiance/LUT_SBDART/Ed0moins_LUT.dat"
# outpath <- "~/Desktop/GreenEdge/Irradiance/Ed0moins_MODISA_noclim_LUT_B"

# Case Ed0Plus B: LUT given by Srikanth (Thetas(19)*Ozone(10)*TauCld(8), 83 lambdas)
lutpath <- "/Users/martigalitapias/Desktop/GreenEdge/Irradiance/LUT_SBDART/Ed0Plus_LUT.dat"
outpath <- "~/Desktop/GreenEdge/Irradiance/Ed0Plus_MODISA_noclim_LUT_B"


# --------------------------------- END EDIT -------------------------------------

for (stn in 1:dim(stnlist)[1]){
  
  stnLAT <- stnlist$V1[stn]
  stnLON <- stnlist$V2[stn]
  stnY <- stnlist$V5[stn]
  stnDOY <- stnlist$V6[stn]
  
  for (dd in 1:length(doydecr)){
    
    searchY <- stnY
    searchDOY <- stnDOY + doydecr[dd]
    if (searchDOY < 1){
      # Case where changing day changes year
      searchY <- stnY - 1
      if (!(searchY %% 4) == 0){
        # Leap year
        searchDOY <- 366
      }else{
        # Normal year
        searchDOY <- 365
      }
    }
    
    # Build file path
    genericname <- sprintf("MYD08_D3.A%04i%03i.051.*.hdf", searchY, searchDOY)
    searchpath <- sprintf("%s/%04i/%03i/%s", mainpath, searchY, searchDOY, genericname)
    
    # Find the complete filename
    getfilepath <- paste("ls", searchpath, sep = " ")
    flist <- system(getfilepath, intern = TRUE)
    
    # Test if file found and either go ahead or put filename in list of files not found
    if (length(flist) == 0){
      
      write.table(stnlist[stn,], file = "files_not_found.txt", append = TRUE)
      
    }else{
      
      filepath <- flist[1] # Here we assume there is only one such file!
      
      # Build outfilepath
      outfilepath <- sprintf("%s/Ed0_%04i_%03i_%0.3f_%0.3f.txt", outpath, searchY, searchDOY, stnLAT, stnLON)
      
      # Call edMODISA
      cmd <- paste("~/svn/Takuvik/Teledetection/Util/Ed/edMODISA", filepath, lutpath, stnLAT, stnLON, searchDOY,">", outfilepath,
                   sep = " ")
      system(cmd)
      
    }
  }  
}
