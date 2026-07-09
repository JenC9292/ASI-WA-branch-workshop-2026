# =============== Plotting in R workshop ===============

## --- packages --------------------------------------------------------------

# If any of these are not installed, run once: install.packages(c("tidyverse", "UpSetR", "ggalluvial", "circlize"))

library(tidyverse)     # dplyr, tidyr, ggplot2, readr, stringr
library(UpSetR)        # UpSet plots
library(ggalluvial)    # Sankey / alluvial plots
library(circlize)      # Circos / chord diagrams
library(here)

set.seed(42)



# ================= LOAD & EXPLORE ==============================

meta <- readRDS(here("data", "GSE168944_metadata.rds"))

## --- extract-cdr3 ----------------------------------------------------------
extract_cdr3 <- function(x, chain = c("TRA", "TRB")) {
  chain <- match.arg(chain)
  seg <- str_extract(x, paste0(chain, ":[^;]+"))   # the CHAIN:... segment
  str_extract(seg, "[^_]+$")                        # CDR3 aa = last token
}

meta <- meta %>%
  mutate(
    paired   = AA_clonotype,                         # TRA + TRB together
    TRA_cdr3 = extract_cdr3(AA_clonotype, "TRA"),
    TRB_cdr3 = extract_cdr3(AA_clonotype, "TRB"),
    lineage  = case_when(
      str_detect(Cluster, "CD8")  ~ "CD8",
      str_detect(Cluster, "CD4")  ~ "CD4",
      str_detect(Cluster, "Treg") ~ "Treg",
      str_detect(Cluster, "NKT")  ~ "NKT",
      TRUE                        ~ "Other"))



# ---- How many unique paired clonotypes overall? --------------------------

## --- overall-counts --------------------------------------------------------
n_cells  <- nrow(meta)
n_clones <- n_distinct(meta$paired)
c(cells = n_cells, unique_paired_clonotypes = n_clones)   # ~64,449 / ~47,115



# A small summary bar showing cells vs unique clonotypes:

## --- glance-plot -----------------------------------------------------------
glance_tbl <- tibble(
  measure = c("Cells", "Unique paired clonotypes"),
  value   = c(n_cells, n_clones))

ggplot(glance_tbl, aes(measure, value, fill = measure)) +
  geom_col(width = 0.6, show.legend = FALSE) +
  geom_text(aes(label = scales::comma(value)), vjust = -0.4) +
  scale_fill_manual(values = c("#34688f", "#8fbcd4")) +
  labs(title = "Dataset at a glance", x = NULL, y = "count") +
  theme_bw()



# ---- Unique paired clonotypes per cluster --------------------------------

## --- per-cluster -----------------------------------------------------------
clones_per_cluster <- meta %>%
  filter(!is.na(paired)) %>%
  distinct(Cluster, paired) %>%
  count(Cluster, name = "n_clonotypes") %>%
  arrange(desc(n_clonotypes))

print(clones_per_cluster, n = Inf)

ggplot(clones_per_cluster,
       aes(x = n_clonotypes, y = reorder(Cluster, n_clonotypes))) +
  geom_col(fill = "#34688f") +
  labs(title = "Unique paired (TRA+TRB) clonotypes per cluster",
       x = "unique clonotypes", y = NULL) +
  theme_bw() +
  theme(axis.text.y = element_text(size = 7))



# ============== UPSET: SHARED CLONOTYPES ACROSS THE 4 MICE =======================

## --- build-mouse-list ------------------------------------------------------
build_mouse_list <- function(df, clone_col) {
  df %>%
    filter(!is.na(.data[[clone_col]])) %>%
    distinct(Mouse, .data[[clone_col]]) %>%
    group_by(Mouse) %>%
    summarise(clones = list(unique(.data[[clone_col]])), .groups = "drop") %>%
    { setNames(.$clones, .$Mouse) }
}

paired_list <- build_mouse_list(meta, "paired")
sapply(paired_list, length)        # distinct paired clonotypes per mouse

# ---- Basic UpSet ---------------------------------------------------------

## --- upset-basic -----------------------------------------------------------
UpSetR::upset(UpSetR::fromList(paired_list),
              nsets    = length(paired_list),
              order.by = "freq",
              mainbar.y.label = "shared paired clonotypes",
              sets.x.label    = "clonotypes per mouse",
              text.scale = 1.3)


# ---- Play around with colours + order by degree --------------------------

# order.by = "degree" with decreasing = TRUE puts the 4-mouse intersection
# first, then all 3-mouse, then 2-mouse, then the single-mouse (private) sets.

## --- upset-degree ----------------------------------------------------------
UpSetR::upset(UpSetR::fromList(paired_list),
              nsets       = length(paired_list),
              order.by    = "degree",
              decreasing  = TRUE,
              main.bar.color = "#2c5f8a",
              sets.bar.color = "#8fbcd4",
              matrix.color   = "#c0392b",
              shade.color    = "#f2dede",
              mainbar.y.label = "shared paired clonotypes",
              sets.x.label    = "clonotypes per mouse",
              text.scale = 1.3)



# If not installed: install.packages("ComplexUpset")

## --- complexupset-lib ------------------------------------------------------
library(ComplexUpset)



# Step 1 — one row per paired clonotype, with a TRUE/FALSE column per mouse (the
# membership table ComplexUpset wants; like fromList() but tidy).

## --- complexupset-membership -----------------------------------------------
mice <- sort(unique(meta$Mouse))

clone_membership <- meta %>%
  filter(!is.na(paired)) %>%
  distinct(paired, Mouse) %>%
  mutate(present = TRUE) %>%
  pivot_wider(names_from = Mouse, values_from = present, values_fill = FALSE)



# Step 2 — attach the cluster composition of each clone. For every clonotype,
# count how many CELLS it has in each cluster. We summarise by lineage
# (CD8/CD4/Treg/...) to keep the legend short; swap to Cluster for the full
# breakdown.

## --- complexupset-composition ----------------------------------------------
clone_composition <- meta %>%
  filter(!is.na(paired)) %>%
  count(paired, lineage, name = "cells")

# one tidy table: membership (per mouse) + per-clone cell counts by lineage
upset_df <- clone_membership %>%
  left_join(clone_composition, by = "paired")



# Step 3 — the plot. The main intersection bars still count CLONOTYPES (one row
# each); the added panel shows the PROPORTION of cells by lineage for the clones
# in each intersection (position = "fill" -> every bar sums to 1).

## --- complexupset-plot -----------------------------------------------------
ComplexUpset::upset(
  upset_df,
  intersect   = mice,
  name        = "mouse",
  sort_intersections_by = "degree",          # 4-mouse first, then 3, 2, 1
  sort_intersections    = "descending",
  min_size    = 1,
  base_annotations = list(
    "shared clonotypes" = intersection_size(
      counts = TRUE,
      mapping = aes(fill = "clonotypes")
    ) + scale_fill_manual(values = c(clonotypes = "#2c5f8a"), guide = "none")
  ),
  annotations = list(
    ## the NEW panel: stacked cell-composition by cluster lineage (PROPORTION)
    "cell proportion by lineage" = (
      ggplot(mapping = aes(fill = lineage)) +
        geom_bar(stat = "count", position = "fill") +    # "fill" -> 0..1 per bar
        scale_fill_brewer(palette = "Set2", name = "lineage") +
        scale_y_continuous(labels = scales::percent_format()) +
        ylab("proportion of cells") +
        theme(legend.position = "right")
    )
  ),
  set_sizes = upset_set_size() + ylab("clonotypes per mouse"),
  width_ratio = 0.2,
  themes = upset_modify_themes(list(
    default = theme_minimal(base_size = 11)
  ))
) +
  ggtitle("Shared paired clonotypes across mice, with cluster composition")


# ---- Task 1 --------------------------------------------------------------

# (a) Repeat the UpSet for the SINGLE chains: build paired_list-style lists with
# build_mouse_list(meta, "TRA_cdr3") and "TRB_cdr3", and plot each.

# (b) Play around with colours: change main.bar.color / matrix.color, and swap
# order.by between "freq" and "degree".





# =========== SANKEY / ALLUVIAL: TRACK SHARED CLONES BETWEEN CLUSTERS ======================

# ---- Which cluster pairs share the most clones? --------------------------

# For each clonotype, list the clusters it appears in, then count co-
# occurrences.

## --- pair-sharing ----------------------------------------------------------
clone_clusters <- meta %>%
  filter(!is.na(paired)) %>%
  distinct(paired, Cluster)

# number of distinct clones shared by each pair of clusters
pair_sharing <- clone_clusters %>%
  inner_join(clone_clusters, by = "paired",
             relationship = "many-to-many") %>%
  filter(Cluster.x < Cluster.y) %>%                 # unordered pairs, no self
  count(Cluster.x, Cluster.y, name = "shared_clones") %>%
  arrange(desc(shared_clones))

head(pair_sharing, 15)   # top-sharing cluster pairs (mostly CD8 exhausted set)


# ---- Helper: build alluvial data for a chosen set of clusters ------------

# Keeps clones shared by >= min_clusters of the chosen clusters, then keeps only
# the top_n most-expanded of those so each ribbon is a distinct, readable clone
# (showing all ~170 shared clones gives an unreadable "swirl" of tiny ribbons).
# Returns one row per clone × cluster with the cell count, plus a short clone ID
# so the ribbons can be coloured and told apart. prefix is the cluster-name
# prefix to strip for tidy axis labels (CD8 by default).

## --- make-alluvial ---------------------------------------------------------
make_alluvial <- function(df, clusters, min_clusters = 2, top_n = 12,
                          prefix = "^C\\d+_CD8_") {
  d <- df %>% filter(Cluster %in% clusters, !is.na(paired))
  # clones present in at least `min_clusters` of the chosen clusters
  shared <- d %>% distinct(paired, Cluster) %>%
    count(paired, name = "k") %>% filter(k >= min_clusters) %>% pull(paired)
  # of those, keep the top_n biggest (by total cells across the chosen clusters)
  top_clones <- d %>% filter(paired %in% shared) %>%
    count(paired, name = "size") %>%
    slice_max(size, n = top_n) %>% pull(paired)
  d %>%
    filter(paired %in% top_clones) %>%
    count(paired, Cluster, name = "cells") %>%
    mutate(
      Cluster  = factor(str_remove(Cluster, prefix),
                        levels = str_remove(clusters, prefix)),
      # short clone label = the TRB CDR3 (compact and unique enough to read)
      clone_id = str_extract(paired, "(?<=TRB:)[^;]+") %>% str_extract("[^_]+$"))
}



# ---- Start with the two most-shared CD8 clusters first -------------------

# Each coloured band is ONE clone. Its height in each column = number of cells
# of that clone in that cluster, so you can read the PROPORTION of a clone that
# sits in each phenotype (e.g. a clone that is mostly Tex_Cd244 with a thin tail
# in Tex_Ccr7).

## --- sankey-two ------------------------------------------------------------
two_clusters <- c("C08_CD8_Tex_Cd244", "C09_CD8_Tex_Ccr7")   # 172 shared clones
allu2 <- make_alluvial(meta, two_clusters, top_n = 12)

ggplot(allu2,
       aes(x = Cluster, stratum = clone_id, alluvium = paired,
           y = cells, fill = clone_id)) +
  geom_alluvium(colour = "white", linewidth = 0.2, alpha = 0.9) +  # the ribbons
  geom_stratum(colour = "white", linewidth = 0.2, width = 0.18) +  # the blocks
  scale_fill_viridis_d(option = "turbo", name = "clone (TRB CDR3)") +
  labs(title = "Clones shared between two CD8 clusters",
       subtitle = "top 12 shared clones; band height = cells of that clone",
       x = NULL, y = "cells") +
  theme_minimal() +
  theme(legend.text = element_text(size = 7), legend.key.size = unit(3, "mm"))



# ---- Add more clusters: four exhausted/proliferating CD8 T cells ---------

## --- sankey-four -----------------------------------------------------------
four_clusters <- c("C07_CD8_Tex_ISG", "C08_CD8_Tex_Cd244",
                   "C09_CD8_Tex_Ccr7", "C12_CD8_Mki67_E2F")
allu4 <- make_alluvial(meta, four_clusters, top_n = 12)

ggplot(allu4,
       aes(x = Cluster, stratum = clone_id, alluvium = paired,
           y = cells, fill = clone_id)) +
  geom_alluvium(colour = "white", linewidth = 0.2, alpha = 0.9) +
  geom_stratum(colour = "white", linewidth = 0.2, width = 0.18) +
  scale_fill_viridis_d(option = "turbo", name = "clone (TRB CDR3)") +
  labs(title = "Tracking clones across four CD8 clusters",
       subtitle = "each clone's band shows how its cells spread across phenotypes",
       x = NULL, y = "cells") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 20, hjust = 1),
        legend.text = element_text(size = 7), legend.key.size = unit(3, "mm"))



# ---- Task 2 --------------------------------------------------------------

# Build the same clone-tracking alluvial for the CD4 clusters.

#   - The strongest CD4 pair is C15_CD4_Tfh  C16_CD4_Itgb1 (start here).
#   - Then add C17_CD4_Th1 and C18_CD4_S1pr1 for a four-cluster version.
#   - Pass prefix = "^C\\d+_CD4_" to make_alluvial() so axis labels are tidy.
#   - Try changing top_n (e.g. 8 or 15) to see more/fewer clones.





# ========= CIRCOS / CHORD: CLONE SHARING BETWEEN TISSUES ================


## --- tissue-matrix ---------------------------------------------------------
inc <- meta %>%
  filter(!is.na(paired)) %>%
  distinct(paired, Tissue) %>%
  mutate(present = 1L) %>%
  pivot_wider(names_from = Tissue, values_from = present, values_fill = 0) %>%
  column_to_rownames("paired") %>%
  as.matrix()

co <- t(inc) %*% inc
diag(co) <- 0
co



## --- tissue-circos ---------------------------------------------------------
grid.col <- setNames(
  c("#E41A1C", "#377EB8", "#4DAF4A", "#984EA3")[seq_len(ncol(co))],
  colnames(co))

circos.clear()
chordDiagram(co, grid.col = grid.col, symmetric = TRUE,
             annotationTrack = c("name", "grid"))
circos.clear()





# ========== CIRCOS BY CLUSTER, FOR EACH MOUSE =======================


## --- cluster-chord-matrix --------------------------------------------------
cluster_chord_matrix <- function(df, mouse, lineage_keep = "CD8") {
  d <- df %>%
    filter(Mouse == mouse,
           str_detect(Cluster, lineage_keep),
           !is.na(paired))
  inc <- d %>%
    distinct(paired, Cluster) %>%
    mutate(present = 1L) %>%
    pivot_wider(names_from = Cluster, values_from = present,
                values_fill = 0) %>%
    column_to_rownames("paired") %>%
    as.matrix()
  m <- t(inc) %*% inc
  diag(m) <- 0
  # tidy labels: drop the "C07_CD8_" style prefix
  short <- str_remove(colnames(m), "^C\\d+_CD8_")
  dimnames(m) <- list(short, short)
  m
}




# ---- Worked example: mouse m01. *Run the whole block as one unit. --------

## --- m01-circos ------------------------------------------------------------
co_m01 <- cluster_chord_matrix(meta, "m01", "CD8")
co_m01

circos.clear()
chordDiagram(co_m01, symmetric = TRUE,
             annotationTrack = c("grid"),
             preAllocateTracks = list(track.height = 0.12))
## add rotated cluster labels so they don't overlap
circos.trackPlotRegion(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.3,
              CELL_META$sector.index, facing = "clockwise",
              niceFacing = TRUE, adj = c(0, 0.5), cex = 0.6)
}, bg.border = NA)
title("m01: shared clonotypes between CD8 clusters")
circos.clear()




# ---- All four mice on one page for comparison ----------------------------

# Resets the layout to a 2×2 grid; each panel is one mouse.

## --- all-mice-circos -------------------------------------------------------
op <- par(mfrow = c(2, 2), mar = c(1, 1, 2, 1))
for (mouse in c("m01", "m02", "m03", "m04")) {
  m <- cluster_chord_matrix(meta, mouse, "CD8")
  circos.clear()
  chordDiagram(m, symmetric = TRUE, annotationTrack = "grid")
  title(paste0(mouse, " - CD8 cluster sharing"))
}
par(op)
circos.clear()




# ---- Task 3 --------------------------------------------------------------

# Make a per-mouse cluster chord diagram for one mouse of your choice using a
# different lineage. For example, the Treg clusters:

# r m <- cluster_chord_matrix(meta, "m02", "Treg") chordDiagram(m, symmetric =
# TRUE)

# Tip: the prefix stripped in cluster_chord_matrix() is CD8-specific — change
# "^C\\d+_CD8_" to "^C\\d+_Treg_" or "^C\\d+_(CD4|CD8|Treg)_" to keep labels
# tidy for other lineages.

