# Preprocess flow cytometry data

library(readxl)
library(dplyr)

# Load xls
df <- readxl::read_excel(path = "~/Desktop/GreenEdge/GE_all_Cytometry version 1.6.xlsx",
                         sheet = "All data")
df <- df[ , c("operation","year","station","ctd","depth_m","Pico_mL","Nano_mL","Crypto_mL","Bact_mL") ]
df$station <- as.numeric(gsub(pattern = "G", replacement = "", x = df$station))

# Subset stations
# df <- df[ df$operation=="Amundsen" & df$station >= 400 , ]
# df <- df[ df$operation=="Amundsen" & df$station %in% c(418,507,512,605,703,707,713,719) , ]
df <- df[ df$operation=="Amundsen" & df$ctd %in% c(110,123,134,155,175,183,192,201) , ]

# Subset surface and SCM
scm <- data.frame(station = c(418,507,512,605,703,707,713,719),
                  depth = c(30,10,15,10,20,20,20,12))
FC <- df[df$depth_m==0,]
for (j in 1:dim(scm)[1]) {
  irow <- which(df$station==scm$station[j] & df$depth_m==scm$depth[j])
  FC <- rbind(FC, df[irow,])
}

# Add scm column
FC$scm <- ifelse(FC$depth_m==0, "surface", "SCM")
  
# Convert mL-1 to L-1
FC[, c("Pico_mL","Nano_mL","Crypto_mL","Bact_mL")] <- df[, c("Pico_mL","Nano_mL","Crypto_mL","Bact_mL")] * 1000

# Rename
FC <- rename(FC,
             stn = station,
             cast = ctd,
             depth = depth_m,
             Pico_FC = Pico_mL,
             Nano_FC = Nano_mL,
             Crypto_FC = Crypto_mL,
             Bact_FC = Bact_mL)

# Save
write.csv(file = "~/Desktop/GreenEdge/MS_DMS_GE16_Elementa/Fig_phyto/flow_cyto_preprocessed.Rda", FC)

