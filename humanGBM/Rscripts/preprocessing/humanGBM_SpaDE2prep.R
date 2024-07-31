###################
#    Human GBM    #
# SpatialDE2 Prep #
###################

## Packages
library(here)
library(zellkonverter)

# Load spatial experiment object

fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs/spe_humanGBM_lowFilt.rds")
spe_humanGBM <- readRDS(fn)
spe_humanGBM

outs_dir <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs")

# convert spatial experiment object int h5ad file using zellkonverter package

out_path <- tempfile("humanGBM", 
                     tmpdir = outs_dir,
                     fileext = ".h5ad")

writeH5AD(spe_humanGBM, file = out_path)

