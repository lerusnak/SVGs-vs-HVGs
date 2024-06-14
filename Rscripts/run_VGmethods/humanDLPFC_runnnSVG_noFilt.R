#####################
#     Run Method    #
#   Method: nnSVG   #
# Data: Human DLPFC #
#####################

# Packages
library(SpatialExperiment)
library(nnSVG)
library(here)


################
#  Load Data   #
################

fn <- here("/projectnb/weber-lr/lerusnak/outputs/humanDLPFC_noFilt.rds")
spe <- readRDS(fn)

assayNames(spe)


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

file <- here("/projectnb/weber-lr/lerusnak/outputs", "spe_humanDLPFC_nnSVG_noFilt.rds")
saveRDS(spe, file = file)