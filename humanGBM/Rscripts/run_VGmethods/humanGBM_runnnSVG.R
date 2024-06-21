###############
#  human GBM  #
#  run nnSVG  #
###############

## Packages
library(SpatialExperiment)
library(nnSVG)
library(here)


#############
# load data #
#############

fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs/spe_humanGBM_lowFilt.rds")
spe <- readRDS(fn)


################
#  Run Method  #
################

# run method and save results, runtime, peak memory usage

# run nnSVG with default parameters (except increasing cores to 10 (n_threads))

# NOTE: NO ADDITIONAL FILTERING WAS PERFORMED WHEN RUNNING THIS METHOD #

set.seed(123)
runtime <- system.time({
  spe <- nnSVG(spe, n_threads = 10)
})

# store runtime in object
metadata(spe) <- list(
  runtime = runtime[["elapsed"]]
)


###############
# save object #
###############

file <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs", "spe_humanGBM_nnSVG.rds")
saveRDS(spe, file = file)