library(tidyverse)
library(cowplot)

options(warn = 1)
rm(list = ls())
setwd("~/projects1/projVR1/analysis1/")

# 0. Plot Themes
theme1 <- theme(
  legend.title = element_blank(), legend.position = "right", 
  legend.text = element_text(size = 6), legend.key.size = unit(8, "pt"),
  axis.title = element_blank(), 
  axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5),
  axis.text.y = element_text(size = 8)
)

# 1. Load and Prep
## atlas V1.1
expr0 <- readRDS("~/atlas0/atlas_V1.1/expr_V1.1a.rds")
expr0 <- expr0 %>%
  filter(type2 %in% c("LPLC1", "LPLC2"), n >= 5) %>%
  group_by(group = type2, time, gene) %>%
  summarize(expr = mean(expr), .groups = "drop")

## perturb-seq
expr48 <- readRDS("data/expr48_b.rds")
expr72 <- readRDS("data/expr72_b.rds")

# 2. Define DEGs
dlist0 <- readRDS("data/deg0.rds")

WT_48_1 <- dlist0[["deg48_WT"]] %>% filter(deg_stringent) %>% arrange(avg_log2FC) %>% pull(gene)
WT_48_2 <- dlist0[["deg48_WT"]] %>% filter(deg_relaxed) %>% arrange(avg_log2FC) %>% pull(gene)

Z3_48_1 <- dlist0[["deg48_Z3"]] %>% filter(deg_stringent) %>% arrange(avg_log2FC) %>% pull(gene)
Z3_48_2 <- dlist0[["deg48_Z3"]] %>% filter(deg_relaxed) %>% arrange(avg_log2FC) %>% pull(gene)
Z3_72_1 <- dlist0[["deg72_Z3"]] %>% filter(deg_stringent) %>% arrange(avg_log2FC) %>% pull(gene)
Z3_72_2 <- dlist0[["deg72_Z3"]] %>% filter(deg_relaxed) %>% arrange(avg_log2FC) %>% pull(gene)

Z4_48_1 <- dlist0[["deg48_Z4"]] %>% filter(deg_stringent) %>% arrange(avg_log2FC) %>% pull(gene)
Z4_48_2 <- dlist0[["deg48_Z4"]] %>% filter(deg_relaxed) %>% arrange(avg_log2FC) %>% pull(gene)
Z4_72_1 <- dlist0[["deg72_Z4"]] %>% filter(deg_stringent) %>% arrange(avg_log2FC) %>% pull(gene)
Z4_72_2 <- dlist0[["deg72_Z4"]] %>% filter(deg_relaxed) %>% arrange(avg_log2FC) %>% pull(gene)


# 2. Correlate DEGs: atlas vs perturb-seq
get_mat0 <- function(expr00, genes00){
  if("genotype" %in% colnames(expr00)){
    expr00 <- expr00 %>%
      mutate(group = paste0(group,":",genotype))
  }
  mat00 <- expr00 %>%
    filter(gene %in% genes00) %>%
    select(gene, group, expr) %>%
    pivot_wider(names_from = group, values_from = expr) %>%
    column_to_rownames("gene") %>% as.matrix()
  return(mat00[genes00,])
}

plot_cor0 <- function(mat1, mat2){
  cor00 <- cor(mat1, mat2, method = "pearson")
  cor00 %>%
    as.data.frame() %>%
    rownames_to_column("type2") %>%
    pivot_longer(cols = -type2, names_to = "group2", values_to = "r") %>%
    separate(group2, into = c("group", "genotype"), sep = ":") %>% 
    mutate(group = factor(group, c("ctrl", "br.Z3", "br.Z4"))) %>%
    ggplot(aes(x = type2, y = r)) +
    geom_boxplot() + geom_jitter(size = 2, width = 0.25, shape = 8, col = "red") + 
    coord_cartesian(ylim = c(-0.5,1)) + geom_hline(yintercept = 0, linetype = "dotted") +
    facet_grid(. ~ group) + 
    theme_bw() + ylab("Pearson's r") + 
    theme(axis.title.x = element_blank(), axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5))
  
}

## atlas DEGs at 48h APF
mat1 <- get_mat0(subset(expr0, time == "48h"), WT_48_1)
mat2 <- get_mat0(expr48, WT_48_1)
plot_cor0(log1p(mat1), log1p(mat2))
ggsave(filename = "results/48h_cor1.svg", units = "in", width = 4, height = 2.5)

mat1 <- get_mat0(subset(expr0, time == "72h"), WT_48_1)
mat2 <- get_mat0(expr72, WT_48_1)
plot_cor0(log1p(mat1), log1p(mat2))
ggsave(filename = "results/72h_cor1.svg", units = "in", width = 4, height = 2.5)

# 3. Plot DEGs
plot_heatmap0 = function(atlas00, perturb00, pgenes00, title00 = "", add.legend = FALSE, max00 = 20){
  p1 <- atlas00 %>%
    filter(gene %in% pgenes00) %>%
    mutate(expr2 = ifelse(expr > max00, max00, expr)) %>%
    ggplot(aes(x = gene, y = group, fill = log1p(expr2))) +
    geom_tile(col = "grey40", lwd = 0.1) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1) +
    scale_x_discrete(limits = pgenes00) + scale_y_discrete(limits = c("LPLC2", "LPLC1")) +
    ggtitle(label = title00) +
    theme1 + theme(legend.position = "none", axis.text.x = element_blank(), plot.title = element_text(size = 9))
  p2 <- perturb00 %>%
    filter(gene %in% pgenes00) %>%
    mutate(expr2 = ifelse(expr > max00, max00, expr)) %>%
    mutate(group = factor(group, rev(c("ctrl", "br.Z3", "br.Z4")))) %>%
    ggplot(aes(x = gene, y = genotype, fill = log1p(expr2))) +
    geom_tile(col = "grey40", lwd = 0.1) +
    scale_fill_distiller(palette = "RdYlBu", direction = -1) +
    scale_x_discrete(limits = pgenes00) +
    facet_grid(group ~ ., scales = "free_y") +
    theme1 + theme(
      axis.text.x = element_text(size = 6, angle = 90, hjust = 1, vjust = 0.5),
      panel.border = element_rect(colour = "black", fill = NA, linewidth = 1), 
      panel.spacing = unit(0.05, "cm"), strip.text = element_text(size = 12),
      legend.position = "none"
    )
  if(add.legend){p2 <- p2 + theme(legend.position = "bottom")}
  plot_grid(p1, p2, ncol = 1, align = "v", axis = "lr",  rel_heights = c(1.5,6))
}

## Atlas DEGs
plot_heatmap0(subset(expr0, time == "48h"), expr48, WT_48_1, "LPLC1 vs LPLC2 (48h)", add.legend = TRUE)
ggsave(filename = "results/48h_heatmap1.svg", units = "in", width = 12, height = 5)

plot_heatmap0(subset(expr0, time == "48h"), expr48, WT_48_2, "LPLC1 vs LPLC2 (48h)", add.legend = TRUE)
ggsave(filename = "results/48h_heatmap2.svg", units = "in", width = 15, height = 5)

plot_heatmap0(subset(expr0, time == "72h"), expr72, WT_48_1, "LPLC1 vs LPLC2 (72h)", add.legend = TRUE)
ggsave(filename = "results/72h_heatmap1.svg", units = "in", width = 12, height = 5)

plot_heatmap0(subset(expr0, time == "72h"), expr72, WT_48_2, "LPLC1 vs LPLC2 (72h)", add.legend = TRUE)
ggsave(filename = "results/72h_heatmap2.svg", units = "in", width = 15, height = 5)

## Perturb-seq DEGs
plot_grid(
  plot_heatmap0(subset(expr0, time == "48h"), expr48, Z3_48_1, "br-Z3 vs ctrl (48h)"),
  plot_heatmap0(subset(expr0, time == "48h"), expr48, Z4_48_1, "br-Z4 vs ctrl (48h)")
)
ggsave(filename = "results/48h_perturb1.svg", units = "in", width = 15, height = 5)

plot_grid(
  plot_heatmap0(subset(expr0, time == "72h"), expr72, Z3_72_1, "br-Z3 vs ctrl (72h)"),
  plot_heatmap0(subset(expr0, time == "72h"), expr72, Z4_72_1, "br-Z4 vs ctrl (72h)")
)
ggsave(filename = "results/72h_perturb1.svg", units = "in", width = 15, height = 5)

# 4. Compare DEGs
plot_fc1 = function(contrast1, contrast2, top0){
  pdata0 <- inner_join(
    dlist0[[contrast1]] %>% select(gene, avg_log2FC, deg_relaxed, deg_stringent),
    dlist0[[contrast2]] %>% select(gene, avg_log2FC, deg_relaxed, deg_stringent),
    by = "gene"
  ) %>% filter(deg_relaxed.x | deg_relaxed.y | gene %in% top0)

  ggplot(pdata0, aes(x = avg_log2FC.x, y = avg_log2FC.y, label = gene)) +
    geom_point(size = 0.5, alpha = 0.7) + 
    ggrepel::geom_label_repel(
      data = subset(pdata0, gene %in% top0), size = 3, fill = "white", alpha = 0.8,
      max.overlaps = Inf, min.segment.length = 0, box.padding = 0.5, direction = "both"
    ) + 
    geom_hline(yintercept = 0, linetype = "dashed") +
    geom_vline(xintercept = 0, linetype = "dashed") +
    xlab(contrast1) + ylab(contrast2) +
    theme_bw()
}

## Genes to highlight
top0 <- c('br', "Eip93F", 'CG12836', '5-HT1A', "beat-IIIb", 'Sulf1', 'AstC')

## 48h APF
plot_grid(
  plot_fc1("deg48_Z3", "deg48_Z4", top0) + 
    labs(x = "br-Z3 vs ctrl (48h)", y = "br-Z4 vs ctrl (48h)") +
    coord_cartesian(xlim = c(-15,15), ylim = c(-7,7)),
  plot_fc1("deg48_WT", "deg48_Z3", top0) + 
    labs(x = "LPLC1 vs LPLC2 (48h)", y = "br-Z3 vs ctrl (48h)") +
    coord_cartesian(xlim = c(-15,15), ylim = c(-7,7)),
  plot_fc1("deg48_WT", "deg48_Z4", top0) +
    labs(x = "LPLC1 vs LPLC2 (48h)", y = "br-Z4 vs ctrl (48h)") +
    coord_cartesian(xlim = c(-15,15), ylim = c(-7,7)),
  align = "hv", axis = "tbrl", nrow = 1
)
ggsave(filename = "results/48h_logfc1.svg", units = "in", width = 12, height = 4)

## 72h APF (LPLC1 vs LPLC2 (at 48h APF))
plot_grid(
  plot_fc1("deg72_Z3", "deg72_Z4", top0) + 
    labs(x = "br-Z3 vs ctrl (72h)", y = "br-Z4 vs ctrl (72h)") +
    coord_cartesian(xlim = c(-15,15), ylim = c(-7,7)),
  plot_fc1("deg48_WT", "deg72_Z3", top0) + 
    labs(x = "LPLC1 vs LPLC2 (48h)", y = "br-Z3 vs ctrl (72h)") +
    coord_cartesian(xlim = c(-15,15), ylim = c(-7,7)),
  plot_fc1("deg48_WT", "deg72_Z4", top0) +
    labs(x = "LPLC1 vs LPLC2 (48h)", y = "br-Z4 vs ctrl (72h)") +
    coord_cartesian(xlim = c(-15,15), ylim = c(-7,7)),
  align = "hv", axis = "tbrl", nrow = 1
)
ggsave(filename = "results/72h_logfc1.svg", units = "in", width = 12, height = 4)
