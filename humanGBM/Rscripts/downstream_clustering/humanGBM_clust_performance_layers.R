###################
#    humanGBM     #
#   clustering    #
#      plots      #
#     LAYERS      #
###################


##############
#  packages  #
##############

library(SpatialExperiment)
library(here)
library(tidyverse)
library(mclust)
library(ggspavis)


###############
#  load data  #
###############

# directories
clustering_dir <- here('/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/outputs')
plots_dir <- here('/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/plots/clustering')

# load data from clustering
fn <- here(clustering_dir, 'res_downstream_clustering_layers.rds')
res_out <- readRDS(fn)

coldata_out <- res_out$coldata_out
spatialcoords_out <- res_out$spatialcoords_out

names(coldata_out)
names(spatialcoords_out)



##################
# 'ground truth' #
#     labels     #
##################

visium_metadata <- read.csv("/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/data/visium_metadata.csv")

spot_labels <- visium_metadata %>%
  select(spot_id, sample, mp, layer) %>%
  filter(sample == "MGH258") 

dim(spot_labels)


## Create Results Dataframe

HVGs_df <- as.data.frame(cbind(coldata_out$HVGs, spatialcoords_out$HVGs)) %>% 
  rownames_to_column(var = "spot_id") %>% arrange(spot_id)
HVGs_df <- full_join(HVGs_df, spot_labels, by = "spot_id")

SVGs_df <- as.data.frame(cbind(coldata_out$nnSVG, spatialcoords_out$nnSVG)) %>% 
  rownames_to_column(var = "spot_id") %>% arrange(spot_id)
SVGs_df <- full_join(HVGs_df, spot_labels, by = "spot_id")

VGdf_list <- list(HVGs_df, SVGs_df)



##################
# match clusters #
##################

# match clusters to ground truth layers for each method

match_HVGs <- c(1, 3, 5, 4, 2)
coldata_out$HVGs$label <- factor(
  coldata_out$HVGs$label, levels = match_HVGs)

match_nnSVG <- c(5, 1, 2, 3, 4)
coldata_out$nnSVG$label <- factor(
  coldata_out$nnSVG$label, levels = match_nnSVG)


###################
#  spatial plots  #
###################

# color palette
pal <- c("#F0027F", "#377EB8", "#4DAF4A", "#984EA3", "#FFD700", "#FF7F00", "#1A1A1A", "#666666")
VGs_df <- list()

for (i in seq_along(coldata_out)) {
  
  # Combine coldata_out and spatialcoords_out and convert to data frame
  combined_df <- as.data.frame(cbind(coldata_out[[i]], spatialcoords_out[[i]])) %>% 
    rownames_to_column(var = "spot_id") %>% arrange(spot_id)
  
  # Full join with spot_labels
  VGs_df[[i]] <- full_join(combined_df, spot_labels, by = "spot_id")
  
  # Plotting
  p <- ggplot(VGs_df[[i]], aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, 
                               color = label)) + 
    geom_point(size = 0.3) + 
    coord_fixed() + 
    scale_y_reverse() + 
    scale_color_manual(values = pal, name = "label") + 
    guides(color = guide_legend(override.aes = list(size = 2))) + 
    ggtitle("Clustering", 
            subtitle = paste0("Human GBM Layers: ", names(coldata_out)[i])) + 
    theme_bw() + 
    theme(panel.grid = element_blank(), 
          axis.title = element_blank(), 
          axis.text = element_blank(), 
          axis.ticks = element_blank())
  
  # Saving plots
  fn <- file.path(plots_dir, paste0("clustering_humanGBM_layers_", names(coldata_out)[i]))
  ggsave(paste0(fn, ".pdf"), plot = p, width = 4, height = 4.25)
  ggsave(paste0(fn, ".png"), plot = p, width = 4, height = 4.25)
}



## Ground Truth plot

# plot ground truth layer labels in spatial coordinates
xyplot_layers <-   ggplot(VGs_df[[1]], aes(x = pxl_col_in_fullres, y = pxl_row_in_fullres, color = layer)) + 
                        geom_point(size = 1) + 
                        coord_fixed() + 
                        scale_y_reverse() + 
                        scale_color_manual(values = pal, name = "label") + 
                        guides(color = guide_legend(override.aes = list(size = 2))) + 
                        ggtitle("Clustering ", subtitle = "Human GBM Layers") + 
                        theme_bw() + 
                        theme(panel.grid = element_blank(), 
                              axis.title = element_blank(), 
                              axis.text = element_blank(), 
                              axis.ticks = element_blank())

ggsave("xyplot_truelayers.png", path = "/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/plots/clustering")


###########
#   ARI   #
###########


## LAYERS ##

# re-label factor levels combine clusters 7 and 8 into a single cluster
# representing white matter

#coldata_out$HVGs$label <- factor(
#  coldata_out$HVGs$label, labels = c("WM", "Layer6", "Layer2", "Layer5", "WM", "Layer4", "Layer1", "Layer3"))

#coldata_out$nnSVG$label <- factor(
#  coldata_out$nnSVG$label, labels = c("WM", "Layer3", "WM", "Layer1", "Layer4", "Layer6", "Layer5", "Layer2"))


# calculate adjusted rand index

ari_HVGs_layers <- round(adjustedRandIndex(VGs_df[[1]]$label, 
                                     VGs_df[[1]]$layer), 4)

ari_nnSVG_layers <-round(adjustedRandIndex(VGs_df[[2]]$label, 
                                      VGs_df[[2]]$layer), 4)


# plot adjusted Rand index

df_layers <- data.frame(
  method = c("HVGs", "nnSVG"), 
  ARI = c(ari_HVGs_layers, ari_nnSVG_layers)
)

df_layers


write_tsv(df_layers, file = "/projectnb/weber-lr/SVGs-vs-HVGs/humanGBM/plots/clustering/ari_layers.tsv")


