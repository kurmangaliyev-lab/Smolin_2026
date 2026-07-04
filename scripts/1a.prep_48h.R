library(tidyverse)
library(Seurat)
library(cowplot)

options(warn = 1)
rm(list = ls())
setwd("~/projects1/projVR1/analysis1/")

# 0. Plot Themes
theme1 <- theme(
  axis.title = element_blank(), axis.text = element_blank(), plot.title = element_blank(), 
  legend.position = "bottom", legend.text = element_text(size = 8), legend.key.size = unit(6, "pt")
)
theme2 <- theme(
  axis.title = element_blank(), legend.title = element_blank(), legend.position = "right", 
  legend.text = element_text(size = 6), legend.key.size = unit(8, "pt"),
  axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 8)
)

# 1. Load and Prep
data0 <- CreateSeuratObject(counts = Read10X("../data1/48h/outs/filtered_feature_bc_matrix/"), project = "projVR1_48h")
data0$Perc_MT <- PercentageFeatureSet(data0, pattern = "mt:")

design0 <- read.table("../data1/meta/design1.txt", header = F, sep = "\t", col.names = c("genotype", "group"))
demux0 <- read.table("../data1/multiplex1/48h/demuxlet0.best", header = T, sep = "\t") %>% 
  select(cellID = BARCODE, is_sng = DROPLET.TYPE, genotype = SNG.BEST.GUESS) %>% 
  mutate(genotype = ifelse(is_sng == "SNG", genotype, is_sng)) %>%
  left_join(design0) %>%
  replace_na(replace = list(group = "X")) %>%
  mutate(group = factor(group, c("ctrl", "br.Z3", "br.Z4", "X"))) %>%
  mutate(genotype = factor(genotype, gtools::mixedsort(unique(genotype)))) %>%
  column_to_rownames("cellID")
table(demux0$genotype, demux0$group)

data0 <- AddMetaData(data0, demux0)
saveRDS(data0, file = "data/data48_1.rds")

# 2. QC/Filter
data0 <- readRDS("data/data48_1.rds")
preQC <- as.matrix(table(data0$genotype, data0$group))
plot_grid(
  VlnPlot(data0, "nCount_RNA", group.by = "group", y.max = 60000) + geom_hline(yintercept = c(10000, 50000)) + NoLegend(),
  VlnPlot(data0, "Perc_MT", group.by = "group", y.max = 50) + geom_hline(yintercept = 10) + NoLegend(),
  nrow = 1
)
data0 <- subset(data0, group != "X" & nCount_RNA > 10000 & nCount_RNA < 50000 & Perc_MT < 10)
data0@meta.data <- droplevels(data0@meta.data)
postQC <- as.matrix(table(data0$genotype, data0$group))
c(sum(preQC), sum(postQC))

# 3. Cluster/Annotate
data0 <- data0 %>%
  NormalizeData() %>%
  FindVariableFeatures() %>%
  ScaleData(vars.to.regress = "nCount_RNA") %>%
  RunPCA()
ElbowPlot(data0, ndims = 50)

dims0 <- 1:7
data0 <- data0 %>%
  RunTSNE(dims = dims0) %>%
  FindNeighbors(dims = dims0) %>%
  FindClusters(resolution = 0.1)

data0 <- RenameIdents(data0, "0"="LPLC2a", "1"="LPLC2b", "2"="LPLC2c", "3"="X4", "4"="X5")
data0$cluster1 <- Idents(data0)

# 4. Average Expression
## By Clusters
expr1 <- AverageExpression(data0, group.by = "cluster1")[[1]] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "cluster1", values_to = "expr") %>%
  select(cluster1, gene, expr)

## By Groups and Replicates (for LPLC2-clusters only)
expr2 <- AverageExpression(subset(data0, idents = c("LPLC2a","LPLC2b","LPLC2c")), group.by = "genotype")[[1]] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(-gene, names_to = "genotype", values_to = "expr") %>%
  mutate(genotype = str_replace(genotype, "line-", "line_")) %>%
  left_join(design0) %>%
  select(group, genotype, gene, expr)

# 5. Plot: cluster_48h.pdf
p1 <- data0@meta.data %>%
  group_by(group, cluster1) %>%
  count() %>%
  ungroup() %>%
  complete(group, cluster1, fill = list(n = 0)) %>%
  ggplot(aes(x = cluster1, y = group, label = n, fill = n)) +
  geom_tile(col = "grey40", lwd = 0.5) + geom_text(size = 3) + 
  scale_fill_distiller(palette = "Oranges", direction = 1) +
  theme2 + theme(legend.position = "none")

pgenes1 <- c("His2A-GFP", "salm", "acj6", "Lim1", "fkh", "br", "dati", "VGlut")
p2 <- expr1 %>%
  filter(gene %in% pgenes1) %>%
  mutate(expr2 = ifelse(expr >= 20, 20, expr)) %>%
  ggplot(aes(x = cluster1, y = gene, fill = log1p(expr2))) +
  geom_tile(col = "grey40", lwd = 0.5) +
  scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  scale_y_discrete(limits = pgenes1) +
  theme2 + theme(legend.position = "bottom")

plot_grid(
  DimPlot(data0, label = T) + theme1,
  DimPlot(data0, group.by = "group") + theme1,
  plot_grid(p1, p2, nrow = 2, align = "v", axis = "btrl", rel_heights = c(1,2)), 
  nrow = 1, align = "h", axis = "rl", rel_widths = c(3,3,2)
)

ggsave(filename = "results/48h_cluster1.svg", units = "in", width = 8, height = 3.5)

# 6. Save
saveRDS(data0, file = "data/data48_2.rds")
saveRDS(expr1, file = "data/expr48_a.rds")
saveRDS(expr2, file = "data/expr48_b.rds")