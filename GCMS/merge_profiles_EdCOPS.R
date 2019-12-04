# Merge casts with irradiance data based on in situ/SBDART Ed and C-OPS

# Load data
fpath <- '~/Desktop/GreenEdge/GCMS/GE2016.profiles.ALL.OK.backup.csv'
opath <- '~/Desktop/GreenEdge/GCMS/GE2016.profiles.ALL.OK.csv'
prof.all <- read.csv(file = fpath, header = T)
Ed <- read.csv(file = '~/Desktop/GreenEdge/Irradiance/par_daily_averaged_cops_amundsen_2016_31012018.csv', header = T)

# Merge
Ed$stn <- as.numeric(gsub(pattern = "G", replacement = "", x = Ed$station))
Ed$depth <- Ed$depth_m
m <- merge(x = prof.all, y = Ed[ ,c('stn','depth','par_d_p24h_ein_m_2_day_1')], by = c('stn','depth'), all.x = T, all.y = F)

# Save
prof.all <- m
write.csv(x = prof.all, file = opath, row.names = F)