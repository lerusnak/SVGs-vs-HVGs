###############
#  human CRC  #
#   run HVG   #
###############

## Packages
library(SpatialExperiment)
library(scran)
library(here)


#############
# load data #
#############

fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC_16um/outputs/spe_humanCRC16_lowFilt.rds")
spe <- readRDS(fn)


################
#  Run Method  #
################

# run method and save results, runtime, peak memory usage

# run HVGs

set.seed(123)
runtime <- system.time({
  dec <- modelGeneVar(spe)
})

# store in object
stopifnot(all(rownames(dec) == rowData(spe)$gene_id))
rowData(spe) <- cbind(rowData(spe), dec)

# calculate ranks
rowData(spe)$rank <- rank(-1 * rowData(spe)$bio, ties.method = "first")

# store runtime in object
metadata(spe) <- list(
  runtime = runtime[["elapsed"]]
)


###############
# save object #
###############

file <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC_16um/outputs", "spe_humanCRC16_HVGs.rds")
saveRDS(spe, file = file)
