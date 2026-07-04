library(tidyverse)
library(Seurat)
library(cowplot)

options(warn = 1)
rm(list = ls())
setwd("~/projects1/projVR1/analysis1/")

# 0. Plot Themes
theme1 <- theme(
  axis.title.x = element_blank(), legend.title = element_blank(), legend.position = "right", 
  legend.text = element_text(size = 6), legend.key.size = unit(8, "pt"),
  axis.text.x = element_text(size = 8, angle = 90, hjust = 1, vjust = 0.5), axis.text.y = element_text(size = 8)
)

# 1. Load and Prep
## atlas V1.1
expr0 <- readRDS("~/atlas0/atlas_V1.1/expr_V1.1a.rds")
expr0 <- expr0 %>%
  filter(class1 == "N", n >= 5) %>%
  group_by(type2, time, gene) %>%
  summarize(expr = mean(expr), .groups = "drop")

## DEGs
dlist0 <- readRDS("data/deg0.rds")

# 2. TFs
tf0 <- read.table("extra/tf0.txt", header = T)[,1]

# 2.1 All TFs
dgenes0 <- dlist0[["deg48_WT"]] %>% filter(deg_stringent) %>% pull(gene)
tf_expr0 <- expr0 %>%
  filter(gene %in% tf0, type2 %in% c("LPLC1", "LPLC2"), time %in% c("48h")) %>%
  group_by(gene, type2) %>% 
  summarize(expr = mean(expr), .groups = "drop") %>%
  pivot_wider(names_from = type2, values_from = expr) %>%
  mutate() %>%
  filter(LPLC1 > 1 | LPLC2 > 1) %>%
  mutate(is_deg = gene %in% dgenes0, fc = LPLC1/LPLC2, delta = LPLC1 - LPLC2) %>%
  arrange(is_deg, delta)
  
pdata1 <- expr0 %>%
  filter(gene %in% tf_expr0$gene, type2 %in% c("LPLC1", "LPLC2"), time %in% c("24h","48h","72h","96h")) %>%
  mutate(expr2 = ifelse(expr > 20, 20, expr)) %>%
  mutate(gene = factor(gene, tf_expr0$gene))

ggplot(pdata1, aes(x = gene, y = time, fill = log1p(expr2))) +
  geom_tile(col = "grey40", lwd = 0.1) +
  scale_fill_distiller(palette = "RdYlBu", direction = -1) +
  scale_x_discrete(limits = tf_expr0$gene) +
  facet_grid(type2 ~ .) +
  theme1 + theme(legend.position = "bottom")

ggsave(filename = "results/atlas_tf1.svg", units = "in", width = 10, height = 3)

## 2.2 DEG TFs
pgenes1 <- tf_expr0 %>% filter(is_deg) %>% arrange(desc(abs(delta))) %>% pull(gene)
pdata1 <- expr0 %>%
  filter(gene %in% pgenes1) %>%
  mutate(gene = factor(gene, pgenes1), is_LPLC = type2 %in% c("LPLC1","LPLC2"))

ymax1 <- pdata1 %>%
  group_by(gene) %>% 
  summarize(ymin = 0, ymax = max(expr), .groups = "drop") %>%
  mutate(ymax = ifelse(ymax < 5, 5, ymax), type2 = NA, time = 0, expr = 0)

ggplot(data = subset(pdata1, !is_LPLC), aes(x = time, y = expr, group = type2)) +
  geom_line(col = "grey40", linewidth = 0.2) +
  geom_line(data = subset(pdata1, is_LPLC), aes(col = type2), linewidth = 1) +
  geom_blank(data = ymax1, aes(y = ymax)) +
  facet_wrap( ~ gene, scales = "free_y", nrow = 2) +
  scale_x_discrete(labels = c("24h","","48h","","72h","","96h")) +
  scale_color_manual(values = c("blue2", "red3")) +
  ylab("TP10K") + theme_bw() +
  theme(legend.position = "bottom", legend.title = element_blank())

ggsave(filename = "results/atlas_tf2.svg", units = "in", width = 6, height = 3.5)

# 3. Side/Beat-Ia
porder1 <- expr0 %>%
  filter(gene %in% c("side", "beat-Ia")) %>%
  mutate(type2 = as.character(type2), gene = as.character(gene)) %>%
  group_by(type2, gene) %>% 
  summarize(expr = median(expr), .groups = "drop") %>%
  arrange(desc(expr)) 

plot_cam1 <- function(pgene00, max00 = 75){
  ptypes00 <- porder1 %>% filter(gene == pgene00) %>% slice(1:max00) %>% pull(type2)
  pal00 <- RColorBrewer::brewer.pal(12, "Paired")[c(12, 7,8, 5,6, 3,4, 1,2)]
  expr0 %>%
    filter(gene == pgene00, type2 %in% ptypes00) %>%
    ggplot(aes(x = type2, y = expr)) +
    geom_boxplot() + geom_point(aes(col = time), size = 1.3) +
    scale_x_discrete(limits = ptypes00) + 
    scale_color_manual(values = pal00) + 
    theme_bw() + theme1 + ggtitle(pgene00) + ylab("TP10K")
}

##
plot_grid(plot_cam1("side"), plot_cam1("beat-Ia"), ncol = 1)
ggsave(filename = "results/atlas_cam1.svg", units = "in", width = 10, height = 4)