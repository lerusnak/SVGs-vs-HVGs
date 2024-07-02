######################
#   Human GBM Data   #
#  Sample: MGH 258   #
######################


# Packages

library(SpatialExperiment)


## load data with read10XVisium 

mgh258.outs <- "/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/data/outs"

spe <- read10xVisium(samples = mgh258.outs,
                     type = "sparse",
                     data = "filtered")

spe

saveRDS(spe, file = "/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs/spe_humanGBM.rds")
