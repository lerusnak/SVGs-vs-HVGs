###################
#   Human DLPFC   #
# SpatialDE2 Prep #
###################

## Packages
library(here)
library(zellkonverter)

# Load filtered spatial experiment object

fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/humanDLPFC_lowFilt.rds")
spe_humanDLPFC <- readRDS(fn)
spe_humanDLPFC

outs_dir <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs")

# convert spatial experiment object int h5ad file using zellkonverter package

out_path <- tempfile("humanDLPFC", 
                     tmpdir = outs_dir,
                     fileext = ".h5ad")
                     
writeH5AD(spe_humanDLPFC, file = out_path)

