###############
#  human CRC  #
# run SPARK-X #
###############

library(SpatialExperiment)
library(SPARK)
library(here)


# ---------
# load data
# ---------

# load data object with preprocessing from previous script

fn <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs/spe_humanCRC_lowFilt.rds")
spe <- readRDS(fn)

dim(spe)


# ----------
# run method
# ----------

# run method and save results, runtime, peak memory usage

# run SPARK-X
set.seed(123)
runtime <- system.time({
  sparkx_out <- sparkx(
    count_in = counts(spe), 
    locus_in = spatialCoords(spe), 
    X_in = NULL, 
    numCores = 1, 
    option = "mixture", 
    verbose = TRUE
  )
})

# results for individual kernels
head(sparkx_out$stats)
head(sparkx_out$res_stest)
# results for combined kernels
head(sparkx_out$res_mtest)

# store results in SPE object
stopifnot(all(rownames(sparkx_out$res_mtest) == rowData(spe)$gene_id))

rowData(spe) <- cbind(rowData(spe), sparkx_out$res_mtest)

# calculate ranks
rowData(spe)$rank <- rank(rowData(spe)$combinedPval, ties.method = "first")

# store runtime in object
metadata(spe) <- list(
  runtime = runtime[["elapsed"]]
)


# -----------
# save object
# -----------

file <- here("/projectnb/weber-lr/SVGs-vs-HVGs/humanCRC/outputs", "spe_humanCRC_SPARKX.rds")
saveRDS(spe, file = file)