######################
#   Human GBM Data   #
######################


# Packages
library(GEOquery)
library(Biobase)
library(Seurat)
library(SpatialExperiment)
library(S4Vectors)


###########################

# Load Data - seurat?

gbm_mgh258 <- Load10X_Spatial(data.dir = "/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/data/outs",
                            filename = "filtered_feature_bc_matrix.h5",
                            assay = "Spatial",
                            slice = "detected_tissue_image.jpg",
                            filter.matrix = TRUE,
                            to.upper = FALSE)


# coerce into summarized experiment object

sce <- as.SingleCellExperiment(gbm_mgh258)

spatialCoords <- as.matrix(
  seu@images[[1]]@coordinates[, c("imagecol", "imagerow")])


# function
seurat_to_spe <- function(seu, sample_id, img_id) {
  ## Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  ## Extract spatial coordinates
  spatialCoords <- as.matrix(
    seu@images[[img_id]]@coordinates[, c("imagecol", "imagerow")])
  
  ## Extract and process image data
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))
  
  imgData <- DataFrame(
    sample_id = sample_id,
    image_id = img_id,
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
  
  # Convert to SpatialExperiment
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    sample_id = sample_id,
    spatialCoords = spatialCoords,
    imgData = imgData
  )
  # indicate all spots are on the tissue
  spe$in_tissue <- 1
  spe$sample_id <- sample_id
  # Return Spatial Experiment object
  spe
}


#####################################

seurat_to_spe2 <- function(seu, img_id) {
  ## Convert to SCE
  sce <- Seurat::as.SingleCellExperiment(seu)
  
  ## Extract spatial coordinates
  spatialCoords <- as.matrix(
    seu@images[[img_id]]@coordinates[, c("imagecol", "imagerow")])
  
  ## Extract and process image data
  img <- SpatialExperiment::SpatialImage(
    x = as.raster(seu@images[[img_id]]@image))
  
  imgData <- DataFrame(
    image_id = img_id,
    data = I(list(img)),
    scaleFactor = seu@images[[img_id]]@scale.factors$lowres)
  
  # Convert to SpatialExperiment
  spe <- SpatialExperiment(
    assays = assays(sce),
    rowData = rowData(sce),
    colData = colData(sce),
    metadata = metadata(sce),
    reducedDims = reducedDims(sce),
    altExps = altExps(sce),
    spatialCoords = spatialCoords,
    imgData = imgData
  )
  # indicate all spots are on the tissue
  spe$in_tissue <- 1
  # Return Spatial Experiment object
  spe
}


# se_ls <- list("posterior1" = se1, "anterior1" = se2)
# spe_ls <- lapply(names(se_ls), function(i) {
 # seurat_to_spatialexperiment(seu = se_ls[[i]], sample_id = i, img_id = i)
#})



spe <- seurat_to_spe(gbm_mgh258, sample_id = "MGH258", img_id = 1)
#spe <- seurat_to_spe2(gbm_mgh258, img_id = 1)





############################

# Load data from geo 

tmp_dir <- tempdir()

gse <- getGEO(GEO = "GSE237183")
show(gse)

gsm <- getGEO(GEO = "GSM7596587") # destdir = "/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/data"
show(gsm)

###################################


getGEOSuppFiles(GEO = "GSE237183", baseDir = tmp_dir, fetch_files = T, makeDirectory = F)
untar(file.path(tmp_dir, "GSE237183_RAW.tar"), exdir = tmp_dir)
geofiles  <- list.files(tmp_dir, full.names = T)

for (cel in cels) {
  gunzip(cel, overwrite = T)
}

cels <- list.files(tmp_dir, pattern = "[.]CEL$", full.names = T)
cels

celpath <- tmp_dir

fns <- list.celfiles(path = celpath, full.names = T)
fns