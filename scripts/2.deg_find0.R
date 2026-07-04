library(tidyverse)
library(Seurat)
library(openxlsx)

options(warn = 1)
rm(list = ls())
setwd("~/projects1/projVR1/analysis1/")

# 1. Load and Prep
## atlas V1.1
data0 <- readRDS("/work1/yerbol/atlas0/atlas_V1.1/data_V1.1a.rds")
data0 <- subset(data0, type2 %in% c("LPLC1", "LPLC2"))

expr0 <- readRDS("~/atlas0/atlas_V1.1/expr_V1.1a.rds")
expr0 <- expr0 %>%
  filter(type2 %in% c("LPLC1", "LPLC2"), n >= 5) %>%
  group_by(group = type2, time, gene) %>%
  summarize(expr = mean(expr), .groups = "drop")

## perturb-seq 48h
data48 <- readRDS("data/data48_2.rds")
data48 <- subset(data48, idents = c("LPLC2a","LPLC2b","LPLC2c"))
Idents(data48) <- data48$group
expr48 <- readRDS("data/expr48_b.rds")

## perturb-seq 72h
data72 <- readRDS("data/data72_2.rds")
data72 <- subset(data72, idents = c("LPLC2a","LPLC2b","LPLC2c"))
Idents(data72) <- data72$group
expr72 <- readRDS("data/expr72_b.rds")

# 2. Find DEGs
find_deg0 <- function(data00, expr00, group1, group2){
  deg00 <- FindMarkers(
    data00, ident.1 = group1, ident.2 = group2,
    min.pct = 0, logfc.threshold = 0, pseudocount.use = 0.1
  ) %>% 
    rownames_to_column("gene") %>%
    mutate(diff.pct = pct.1 - pct.2) %>%
    relocate(p_val, p_val_adj, .after = last_col())
  deg00 <- expr00 %>%
    group_by(group, gene) %>%
    summarize(expr = mean(expr), .groups = "drop") %>%
    pivot_wider(names_from = group, values_from = expr) %>%
    select(gene, !!group1, !!group2) %>%
    left_join(deg00) %>%
    arrange(p_val)
}

##
dlist0 <- list(
  "deg48_WT" = find_deg0(subset(data0, time == "48h"), subset(expr0, time == "48h"), "LPLC1", "LPLC2"),
  "deg48_Z3" = find_deg0(data48, expr48, "br.Z3", "ctrl"),
  "deg48_Z4" = find_deg0(data48, expr48, "br.Z4", "ctrl"),
  "deg72_Z3" = find_deg0(data72, expr72, "br.Z3", "ctrl"),
  "deg72_Z4" = find_deg0(data72, expr72, "br.Z4", "ctrl")
)

# 3. Filter atlas genes
dlist0 <- lapply(dlist0, function(deg00) deg00 %>% filter(gene %in% levels(expr0$gene)))

# 4. Significant DEGs
find_sig0 <- function(deg00){
  deg00 %>%
    mutate(
      deg_relaxed = p_val_adj < 0.01 & abs(avg_log2FC) > 1 & (pick(2)[[1]] >= 0.5 | pick(3)[[1]] >= 0.5),
      deg_stringent = p_val_adj < 0.01 & abs(avg_log2FC) > 2 & (pick(2)[[1]] >= 1 | pick(3)[[1]] >= 1)
    )
}
dlist0 <- lapply(dlist0, find_sig0)

# 5. Save DEGs
saveRDS(dlist0, file = "data/deg0.rds")

names(dlist0) <- c(
  "LPLC1_vs_LPLC2_48h",
  "Z3_vs_CTRL_48h", "Z4_vs_CTRL_48h",
  "Z3_vs_CTRL_72h", "Z4_vs_CTRL_72h"
)
write.xlsx(dlist0, file = "results/deg0.xlsx")
