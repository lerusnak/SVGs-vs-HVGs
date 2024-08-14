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
library(tidyverse)
library(mclust)
library(ggspavis)
library(readr)


###############
#  load data  #
###############

# directories
clustering_dir <- here('/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/outputs')
plots_dir <- here('/projectnb/weber-lr/SVGs-vs-HVGs/mouseEmbryo/plots/clustering')

# load data from clustering
fn <- here(clustering_dir, 'res_downstream_kmeansclust.rds')
res_out <- readRDS(fn)

coldata_out <- res_out$coldata_out
spatialcoords_out <- res_out$spatialcoords_out

names(coldata_out)
names(spatialcoords_out)


##################
# match clusters #
##################

# match clusters to ground truth layers for each method

match_HVGs <- c(8, 1, 3, 4, 5, 6, 7, 2, 9, 10, 11)
coldata_out$HVGs$label <- factor(
  coldata_out$HVGs$label, levels = match_HVGs)

match_nnSVG <- c(1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11)
coldata_out$nnSVG$label <- factor(
  coldata_out$nnSVG$label, levels = match_nnSVG)

match_SPARKX <- c(1, 2, 10, 4, 5, 6, 7, 8, 9, 3, 11)
coldata_out$SPARKX$label <- factor(
  coldata_out$SPARKX$label, levels = match_SPARKX)

match_SDE2 <- c(1, 2, 10, 8, 4, 6, 9, 11, 7, 3, 5)
coldata_out$SpatialDE2$label <- factor(
  coldata_out$SpatialDE2$label, levels = match_SDE2)

###################
#  spatial plots  #
###################

for (i in seq_along(coldata_out)) {
  
  df <- as.data.frame(cbind(coldata_out[[i]], spatialcoords_out[[i]]))
  
  ggplot(df, aes(x = X2, y = X3, 
                 color = label)) + 
    geom_point(size = 0.1, alpha = 0.2) + 
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
                  paste0("kmeansclust_mouseEmbryo_", names(coldata_out[i])))
  ggsave(paste0(fn, ".png"), width = 4, height = 4.25)
}



################################################

## Proportions of beads assigned to each cluster comapering VG methods

# HVG Proportions
hvg_props <- as.data.frame(table(coldata_out[["HVGs"]]$label)/length(coldata_out[["HVGs"]]$label))
colnames(hvg_props) <- c("cluster", "proportion")
hvg_props$proportion <- round(hvg_props$proportion, 4)
hvg_props

method <- rep("HVGs", 11)

hvg_props <- cbind(hvg_props, method)
hvg_props

# SVG Proportions
# nnSVG
nnsvg_props <- as.data.frame(table(coldata_out[["nnSVG"]]$label)/length(coldata_out[["nnSVG"]]$label))
colnames(nnsvg_props) <- c("cluster", "proportion")
nnsvg_props$proportion <- round(nnsvg_props$proportion, 4)
nnsvg_props

method <- rep("nnSVG", 11)

nnsvg_props <- cbind(nnsvg_props, method)
nnsvg_props

# SparkX
sparkx_props <- as.data.frame(table(coldata_out[["SPARKX"]]$label)/length(coldata_out[["SPARKX"]]$label))
colnames(sparkx_props) <- c("cluster", "proportion")
sparkx_props$proportion <- round(sparkx_props$proportion, 4)
sparkx_props

method <- rep("SPARKX", 11)

sparkx_props <- cbind(sparkx_props, method)
sparkx_props

# Moran's I
mori_props <- as.data.frame(table(coldata_out[["MorI"]]$label)/length(coldata_out[["MorI"]]$label))
colnames(mori_props) <- c("cluster", "proportion")
mori_props$proportion <- round(mori_props$proportion, 4)
mori_props

method <- rep("MoransI", 11)

mori_props <- cbind(mori_props, method)
mori_props


# Moran's I
sde2_props <- as.data.frame(table(coldata_out[["SpatialDE2"]]$label)/length(coldata_out[["SpatialDE2"]]$label))
colnames(sde2_props) <- c("cluster", "proportion")
sde2_props$proportion <- round(sde2_props$proportion, 4)
sde2_props

method <- rep("SpatialDE2", 11)

sde2_props <- cbind(sde2_props, method)
sde2_props


# Combined VG proportions data frame
vg_props <- rbind(hvg_props, nnsvg_props, sparkx_props, mori_props, sde2_props)
vg_props

vg_props_wide <- vg_props %>% pivot_wider(names_from = cluster, values_from = proportion) %>% 
  select(method, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `9`, `10`, `11`)
vg_props_wide

write_tsv(vg_props_wide, file = here(plots_dir, 'prop_spots_percluster.tsv'))

ggplot(data = vg_props, aes(x = method, y = proportion, fill = cluster)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(title = "Proporiton of beads belonging to ")
ggsave(filename = here(plots_dir, 'beadsprop_percluster_plot.png'))



## Frequencies

# HVG Frequencies
hvg_freqs <- as.data.frame(table(coldata_out[["HVGs"]]$label))
colnames(hvg_freqs) <- c("cluster", "frequency")
hvg_freqs

method <- rep("HVGs", 11)

hvg_freqs <- cbind(hvg_freqs, method)
hvg_freqs

# SVG Frequencies
nnsvg_freqs <- as.data.frame(table(coldata_out[["nnSVG"]]$label))
colnames(nnsvg_freqs) <- c("cluster", "frequency")
nnsvg_freqs

method <- rep("nnSVG", 11)

nnsvg_freqs <- cbind(nnsvg_freqs, method)
nnsvg_freqs

# Spark-X
sparkx_freqs <- as.data.frame(table(coldata_out[["SPARKX"]]$label))
colnames(sparkx_freqs) <- c("cluster", "frequency")
sparkx_freqs

method <- rep("SPARKX", 11)

sparkx_freqs <- cbind(sparkx_freqs, method)
sparkx_freqs

# Moran's I
mori_freqs <- as.data.frame(table(coldata_out[["MorI"]]$label))
colnames(mori_freqs) <- c("cluster", "frequency")
mori_freqs

method <- rep("MorI", 11)

mori_freqs <- cbind(mori_freqs, method)
mori_freqs

# SpatialDE2
sde2_freqs <- as.data.frame(table(coldata_out[["SpatialDE2"]]$label))
colnames(sde2_freqs) <- c("cluster", "frequency")
sde2_freqs

method <- rep("SpatialDE2", 11)

sde2_freqs <- cbind(sde2_freqs, method)
sde2_freqs


# Combined VG frequencies data frame
vg_freqs <- rbind(hvg_freqs, nnsvg_freqs, sparkx_freqs, mori_freqs, sde2_freqs)
vg_freqs

vg_freqs_wide <- vg_freqs %>% pivot_wider(names_from = cluster, values_from = frequency) %>% 
  select(method, `1`, `2`, `3`, `4`, `5`, `6`, `7`, `8`, `9`, `10`, `11`)
vg_freqs_wide

write_tsv(vg_freqs_wide, file = here(plots_dir, 'beadsfreq_percluster.tsv'))

# Stacked Barplot
ggplot(data = vg_freqs, aes(x = method, y = frequency, fill = cluster)) +
  geom_bar(stat = "identity") +
  theme_bw() +
  labs(title = "frequency of beads belonging to ")
ggsave(filename = here(plots_dir, 'beadsfreq_percluster_plot.png'))


##############################################################

# Comparing HVG cluster labels to nnSVG cluster labels

nnsvg_vs_hvg_ari <- adjustedRandIndex(coldata_out[["nnSVG"]]$label,coldata_out[["HVGs"]]$label)
nnsvg_vs_hvg_ari # 0.558479
