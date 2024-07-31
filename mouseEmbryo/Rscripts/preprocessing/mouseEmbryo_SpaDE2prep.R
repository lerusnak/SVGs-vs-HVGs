###################
#  mouse Embryo   #
# SpatialDE2 Prep #
###################

## Packages
library(here)
library(zellkonverter)

# Load spatial experiment object

fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs/spe_mouseEmbryo_lowFilt.rds")
spe_mouseEmbryo <- readRDS(fn)
spe_mouseEmbryo

outs_dir <- here("/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs")

# convert spatial experiment object int h5ad file using zellkonverter package

out_path <- tempfile("mouseEmbryo", 
                     tmpdir = outs_dir,
                     fileext = ".h5ad")

writeH5AD(spe_mouseEmbryo, file = out_path)

