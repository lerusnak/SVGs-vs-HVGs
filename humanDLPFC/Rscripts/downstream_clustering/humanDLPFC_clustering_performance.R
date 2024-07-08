###################
#   human DLPFC   #
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
clustering_dir <- here('/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs')
plots_dir <- here('/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/plots/clustering')

# load data from clustering
fn <- here(clustering_dir, 'res_downstream_clustering.rds')
res_out <- readRDS(fn)

coldata_out <- res_out$coldata_out
spatialcoords_out <- res_out$spatialcoords_out

names(coldata_out)
names(spatialcoords_out)


##################
# match clusters #
##################

# match clusters to ground truth layers for each method

match_HVGs <- c(2, 4, 6, 3, 8, 7, 1, 5)
coldata_out$HVGs$label <- factor(
  coldata_out$HVGs$label, levels = match_HVGs)

match_nnSVG <- c(6, 7, 5, 8, 2, 4, 3, 1)
coldata_out$nnSVG$label <- factor(
  coldata_out$nnSVG$label, levels = match_nnSVG)


###################
#  spatial plots  #
###################

# color palette
pal <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", "#FF7F00", "#1A1A1A", "#666666")


for (i in seq_along(coldata_out)) {
  
  df <- as.data.frame(cbind(coldata_out[[i]], spatialcoords_out[[i]]))
  
  ggplot(df, aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                 color = label)) + 
    geom_point(size = 0.3) + 
    coord_fixed() + 
    scale_y_reverse() + 
    scale_color_manual(values = pal, name = "label") + 
    guides(color = guide_legend(override.aes = list(size = 2))) + 
    ggtitle("Clustering", 
            subtitle = paste0("Human DLPFC: ", names(coldata_out[i]))) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  fn <- file.path(plots_dir, 
                  paste0("clustering_humanDLPFC_", names(coldata_out[i])))
  ggsave(paste0(fn, ".pdf"), width = 4, height = 4.25)
  ggsave(paste0(fn, ".png"), width = 4, height = 4.25)
}


## Ground Truth plot

spe <- readRDS("/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/outputs/humanDLPFC_lowFilt.rds")

# plot ground truth labels in spatial 
xyplot_groundtruth <- plotSpots(spe, annotate = "ground_truth", 
                                pal = "libd_layer_colors",
                                point_size = 0.7) 
ggsave("xyplot_groundtruth.png", path = "/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/plots/clustering")


###########
#   ARI   #
###########

# re-label factor levels combine clusters 7 and 8 into a single cluster
# representing white matter

 coldata_out$HVGs$label <- factor(
  coldata_out$HVGs$label, labels = c("WM", "Layer6", "Layer2", "Layer5", "WM", "Layer4", "Layer1", "Layer3"))

 coldata_out$nnSVG$label <- factor(
  coldata_out$nnSVG$label, labels = c("WM", "Layer3", "WM", "Layer1", "Layer4", "Layer6", "Layer5", "Layer2"))


# calculate adjusted rand index

ari_HVGs <- round(adjustedRandIndex(coldata_out$HVGs$label, 
                              coldata_out$HVGs$ground_truth), 4)

ari_nnSVG <- round(adjustedRandIndex(coldata_out$nnSVG$label, 
                               coldata_out$nnSVG$ground_truth), 4)


# plot adjusted Rand index

df <- data.frame(
  method = c("HVGs", "nnSVG"), 
  ARI = c(ari_HVGs, ari_nnSVG)
)

df

write_tsv(df, file = "/projectnb/weber-lr/SVGs-vs-HVGs/humanDLPFC/plots/clustering/ari_table.tsv")

pal_methods <- c("blue3", "deepskyblue2")


ggplot(df, aes(x = method, y = ARI, shape = method, color = method)) + 
  geom_point(stroke = 1.5, size = 2) + 
  scale_shape_manual(values = c(4, 3)) + 
  scale_color_manual(values = pal_methods) + 
  ylim(c(0, 1)) + 
  ggtitle("Downstream clustering performance") + 
  labs(y = "Adjusted Rand Index") + 
  theme_bw() + 
  theme(axis.title.x = element_blank())

fn <- file.path(plots_dir, "summary_clustering_performance_ARI")
ggsave(paste0(fn, ".pdf"), width = 4.75, height = 3.1)
ggsave(paste0(fn, ".png"), width = 4.75, height = 3.1)
