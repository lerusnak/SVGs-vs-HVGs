#####################
#     Run Method    #
#    Method: HVGs   #
# Data: Human DLPFC #
#####################

# Packages
library(SpatialExperiment)
library(scran)
library(here)


################
#  Load Data   #
################

fn <- here("/projectnb/weber-lr/lerusnak/outputs/spe_humanDLPFC_preprocessed.rds")
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

file <- here("/projectnb/weber-lr/lerusnak/outputs", "spe_humanDLPFC_HVGs_noFilt.rds")
saveRDS(spe, file = file)