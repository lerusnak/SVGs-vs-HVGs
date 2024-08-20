###################
#    Human CRC    #
# SpatialDE2 Prep #
###################

## Packages
library(here)
library(zellkonverter)

# Load filtered spatial experiment object

fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC_16um/outputs/spe_humanCRC16_lowFilt.rds")
spe_humanCRC <- readRDS(fn)
spe_humanCRC

# Add spatial information to colData

colData(spe_humanCRC)$spatial0 <- spatialCoords(spe_humanCRC)[,1]
colData(spe_humanCRC)$spatial1 <- spatialCoords(spe_humanCRC)[,2]

head(colData(spe_humanCRC))

spe_humanCRC


outs_dir <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC_16um/outputs")

# convert spatial experiment object int h5ad file using zellkonverter package

out_path <- tempfile("humanCRC16um_", 
                     tmpdir = outs_dir,
                     fileext = ".h5ad")

writeH5AD(spe_humanCRC, file = out_path)

