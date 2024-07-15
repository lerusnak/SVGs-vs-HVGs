###################
#  mouse Embryo   #
#   clustering    #
#     plots       #
###################


##############
#  packages  #
##############

library(SpatialExperiment)
library(here)
library(mclust)
library(ggplot2)
library(ggspavis)
library(readr)


###############
#  load data  #
###############

# directories
clustering_dir <- here('/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs')
plots_dir <- here('/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/plots/clustering')

# load data from clustering
fn <- here(clustering_dir, 'res_downstream_kmeansclust_subset.rds')
res_out <- readRDS(fn)

coldata_out <- res_out$coldata_out
spatialcoords_out <- res_out$spatialcoords_out

names(coldata_out)
names(spatialcoords_out)


##################
# match clusters #
##################

# match clusters to ground truth layers for each method

#match_HVGs <- c(2, 4, 6, 3, 8, 7, 1, 5)
#coldata_out$HVGs$label <- factor(
#  coldata_out$HVGs$label, levels = match_HVGs)

#match_nnSVG <- c(6, 7, 5, 8, 2, 4, 3, 1)
#coldata_out$nnSVG$label <- factor(
#  coldata_out$nnSVG$label, levels = match_nnSVG)


###################
#  spatial plots  #
###################

# color palette
pal <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", "#FF7F00", "#1A1A1A", "#666666")


for (i in seq_along(coldata_out)) {
  
  df <- as.data.frame(cbind(coldata_out[[i]], spatialcoords_out[[i]]))
  
  ggplot(df, aes(x = X2, y = X3, 
                 color = label)) + 
    geom_point(size = 0.2, alpha = 0.5) + 
    coord_fixed() + 
    scale_y_reverse() + 
    #scale_color_manual(values = pal, name = "label") + 
    guides(color = guide_legend(override.aes = list(size = 2))) + 
    ggtitle("Clustering", 
            subtitle = paste0("Mouse Embryo: ", names(coldata_out[i]))) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(plots_dir, 
                  paste0("kmeansclust_subset_mouseEmbryo_", names(coldata_out[i])))
  ggsave(paste0(fn, ".pdf"), width = 4, height = 4.25)
  ggsave(paste0(fn, ".png"), width = 4, height = 4.25)
}

