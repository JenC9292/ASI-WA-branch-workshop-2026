# =============== ASI Workshop: Heatmaps, boxplots, violin plots, functions and loops ===============


# Participant version
# Run this script line-by-line during the workshop.

library(here)
library(tidyverse)




# =============== Activity 27: Introduction to Cluster-Based Visualisation ===============

## ---- Data Import ----

# read in the main two files we will use from before
counts_matrix <- readRDS(here("data", "counts_matrix.rds"))
GSE168944_metadata <- readRDS(here("data", "GSE168944_metadata.rds"))




# =============== Activity 28: From PCA to Published Clusters ===============

# Look at the Cluster Column
head(GSE168944_metadata[ , c('barcode', 'Cluster', 'Tissue', 'Batch', 'Sample')])

unique(GSE168944_metadata$Cluster)



## ---- Activity ----

# TODO: Identify how many cells are in each cluster.

# TODO: Identify how many different clusters there are.



# Convert Cluster to a Factor
GSE168944_metadata <- GSE168944_metadata %>%
  mutate(
    Cluster = factor(Cluster)
  )


## ---- Downsample - Metadata ----
GSE168944_metadata <- GSE168944_metadata[GSE168944_metadata$barcode %in% colnames(counts_matrix), ]

GSE168944_metadata$Cell_annotation <- gsub("^[^_]*_", "", GSE168944_metadata$Cluster)

# check it worked as we want
head(GSE168944_metadata$Cell_annotation)
tail(GSE168944_metadata$Cell_annotation)

# Subset the metadata file
metadata <- GSE168944_metadata[grepl("^CD8", GSE168944_metadata$Cell_annotation), ]

# How many cells do we have now? Are they only CD8s?
nrow(metadata)

metadata %>%
  count(
    Cell_annotation,
    sort = TRUE
  )

rm(GSE168944_metadata)


## ---- Lets downsample to keep only 100 cells of each cell type ----

# Note: reduce 100 to 50 if plots run slowly on your laptop.
set.seed(42)

metadata <- metadata %>%
  group_by(Cell_annotation) %>% # Split the metadata by cell annotation
  group_modify(~ slice_sample( # Randomly keep up to 100 cells from each group
    .x,
    n = min(nrow(.x), 100)
  )) %>%
  ungroup() # Combine the groups back into one data frame

nrow(metadata)

metadata %>%
  count(
    Cell_annotation,
    sort = TRUE
  )


## ---- Downsample to Match - Expression Matrix ----
counts_matrix[1:5, 1:5]

counts_matrix <- counts_matrix %>%
  tibble::column_to_rownames("gene") %>%
  as.matrix()

# check
counts_matrix[1:5, 1:5]

counts_matrix <- counts_matrix[ , metadata$barcode, drop = FALSE]



## ---- Check the Metadata to the Expression Matrix Match ----

all(
  colnames(counts_matrix) %in% metadata$barcode
)

all(
  metadata$barcode %in% colnames(counts_matrix)
)


## ---- Check the Cell Order is the Same in the Metadata as the Expression Matrix ----
all(
  metadata$barcode == colnames(counts_matrix)
)

# Discussion: 
    # Q1. What is the difference between the code to compare the cells in the metadata to the count 
    #     matrix and to check the order of the cells?
    # Q2. Why did we have to compare the cells with both the metadata first and then the counts_matrix 
    #     first, but only once for the order of cells?




# =============== Activity 29: Explore Published Clusters ===============

## ---- Plot Cluster Sizes ----
metadata %>%
  count(
    Cell_annotation,
    sort = TRUE
  ) %>%
  ggplot(
    aes(
      x = reorder(Cell_annotation, n),
      y = n
    )
  ) +
  geom_col() +
  coord_flip() +
  labs(
    title = "Number of cells per published Cell_annotation",
    x = "Cell_annotation",
    y = "Number of cells"
  )


## ---- Cell_annotation Composition by Tissue ----
metadata %>%
  count(
    Tissue,
    Cell_annotation
  ) %>%
  ggplot(
    aes(
      x = Cell_annotation,
      y = n,
      fill = Tissue
    )
  ) +
  geom_col() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    title = "Cell_annotation composition by tissue",
    x = "Cell_annotation",
    y = "Number of cells"
  )


#This can be hard to interpret - it can also be displayed as a proportion bar plot
Cell_annotation_summary <- metadata %>%
  count(Tissue, Cell_annotation) %>%
  group_by(Tissue) %>%
  mutate(prop = n / sum(n))

Cell_annotation_summary

ggplot(Cell_annotation_summary,
       aes(Tissue, n, fill = Cell_annotation)) +
    geom_col(position = "fill",
    colour = "white",
    linewidth = 0.2) +
  scale_y_continuous(
    labels = scales::percent
  ) +
  labs(
    title = "Cell_annotation composition by tissue",
    x = "Tissue",
    y = "Percentage of cells",
    fill = "Cell_annotation"
  ) +
  theme_bw() +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

rm(Cell_annotation_summary)




# =============== Activity 30: Normalising and Scaling Expression Data ===============

## ---- Normalise, log-transform and scale expression before plotting ----

# Normalise counts
dim(counts_matrix)

counts_matrix[1:5, 1:5]

# Count genes with no expression across all cells
sum(rowSums(counts_matrix) == 0)

# Check the dimensions of the counts matrix
dim(counts_matrix)

# Remove genes with no expression across all cells
counts_matrix <- counts_matrix[
  rowSums(counts_matrix) > 0,
]

# Check the dimensions of the filtered matrix
dim(counts_matrix)

# Check all cells have a minimum of 200 reads
all(colSums(counts_matrix) >= 200)

# Calculate library size per cell
library_size <- colSums(counts_matrix)

summary(library_size)

# Normalise counts per 10,000 reads
counts_norm <- sweep(
  counts_matrix,
  2,
  library_size,
  FUN = "/"
) * 10000 # Divide each cell by its total counts, then scale to counts per 10,000

# Log-transform the normalised counts
counts_log <- log1p(counts_norm)

counts_log[1:5, 1:5]

# Scale expression by gene
counts_scaled <- t(
  scale( 
    t(counts_log)
  )
)
# scale() works column-wise, so we transpose the matrix to scale each gene

dim(counts_scaled)

counts_scaled[1:5, 1:5]

# Confirm the metadata still matches the matrix
all(
  colnames(counts_scaled) %in% metadata$barcode
)

rm(counts_log, counts_norm, library_size, counts_matrix)




# =============== Activity 31: Heatmaps with ComplexHeatmap ===============

# Load Packages
library(tidyverse)
library(ComplexHeatmap)
library(circlize)

# Select marker genes for plotting
marker_genes <- c(
  "Cd3d",
  "Cd3e",
  "Cd4",
  "Cd8a",
  "Cd8b1",
  "Nkg7",
  "Gzmb",
  "Gzmk",
  "Ifng",
  "Pdcd1",
  "Tox",
  "Cxcr6",
  "Ccr7",
  "Sell",
  "Mki67"
)

marker_genes <- marker_genes[
  marker_genes %in% rownames(counts_scaled)
]

marker_genes

## ---- Option 1: Create a smaller matrix first ----

# create a new matrix that contains only the marker genes we want to plot
heatmap_matrix <- counts_scaled[
  marker_genes,
  ,
  drop = FALSE
]

dim(heatmap_matrix)

Heatmap(
  heatmap_matrix,
  name = "scaled\nexpression",
  show_column_names = FALSE
)

## ---- Option 2: Select the genes inside `Heatmap()` ----

# keep the full expression matrix unchanged and select the marker genes directly inside the heatmap call
Heatmap(
  counts_scaled[
    marker_genes,
    ,
    drop = FALSE
  ],
  name = "scaled\nexpression",
  show_column_names = FALSE
)

rm(heatmap_matrix)





# =============== Activity 32: Basic Heatmap ===============

## ---- Activity: create a heatmap using only the first 8 marker genes ----
#Hint:
  # subset the matrix by rows
  # keep all columns
  # use Heatmap()


# Your code here






# =============== Activity 33: Comparing Global Colour Scaling with Row-wise Min-Max Scaling ===============

## ---- Create a row-wise min-max scaled matrix ----

# Function to scale a vector between 0 and 1
min_max_scale <- function(x) {

  # Subtract the minimum value so the smallest value becomes 0,
  # then divide by the range so the largest value becomes 1
  (x - min(x)) / (max(x) - min(x))
}

# Function to scale a vector between 0 and 1
min_max_scale <- function(x) {

  # Check whether all values are identical
  # If they are, the range is zero and division would not be possible
  if (max(x) == min(x)) {

    # Return a vector of zeros with the same length as the input
    return(
      rep(
        0,
        length(x)
      )
    )
  }

  # Otherwise, perform min-max scaling
  # Smallest value becomes 0
  # Largest value becomes 1
  # All other values fall between 0 and 1
  (x - min(x)) / (max(x) - min(x))
}

heatmap_minmax <- t(
  apply(
    counts_scaled,
    1,
    min_max_scale
  )
)

dim(heatmap_minmax)
dim(counts_scaled)


# What do the min and max values look like now?
min(counts_scaled[1,])
max(counts_scaled[1,])

min(heatmap_minmax[1,])
max(heatmap_minmax[1,])



# Plot the original scaled values
Heatmap(
  counts_scaled[marker_genes[1:8], , drop = FALSE],
  name = "scaled\nexpression",
  show_column_names = FALSE,
  cluster_rows = FALSE
)

# Plot the row-wise min-max scaled values
Heatmap(
  heatmap_minmax[marker_genes[1:8], , drop = FALSE],
  name = "row-wise\nmin-max",
  show_column_names = FALSE,
  cluster_rows = FALSE
)

# Compare the two heatmaps side-by-side
Heatmap(
  counts_scaled[marker_genes[1:8], , drop = FALSE],
  name = "scaled\nexpression",
  show_column_names = FALSE,
  cluster_rows = FALSE
) +
  Heatmap(
    heatmap_minmax[marker_genes[1:8], , drop = FALSE],
    name = "row-wise\nmin-max",
    show_column_names = FALSE,
    cluster_rows = FALSE
)

rm(min_max_scale, heatmap_minmax)


## ---- Custom colour scale ----

### ---- Option 1: Manual colours ----

# We can control the colours using colorRamp2() from circlize

# Use the full range of the data to automatically set the colour scale
limits <- range(counts_scaled[marker_genes[1:8], , drop = FALSE])
# note that you need to set the range based on the genes you are plotting, not all genes

Heatmap(
  counts_scaled[marker_genes[1:8], , drop = FALSE],
  name = "scaled\nexpression",
  col = colorRamp2(
    c(limits[1], mean(limits), limits[2]),
    c("#542788", "#F7F7F7", "#E08214")),
  show_column_names = FALSE
)

rm(limits)

### ---- Option 2: Colour palettes ----

library(viridisLite)

limits <- seq(
  min(counts_scaled[marker_genes[1:8], , drop = FALSE]),
  max(counts_scaled[marker_genes[1:8], , drop = FALSE]),
  length.out = 100)
  
Heatmap(
  counts_scaled[marker_genes[1:8], , drop = FALSE],
  name = "scaled\nexpression",
  col = colorRamp2(
    limits,
    magma(100)
  ),
  show_column_names = FALSE
)

rm(limits)





# =============== Activity 34: Add Metadata to the Heatmap ===============

## ---- Create a Column Annotation ----

column_anno <- HeatmapAnnotation(
  Cell_annotation = metadata$Cell_annotation,
  Tissue = metadata$Tissue
)

Heatmap(
  counts_scaled[marker_genes, , drop = FALSE],
  name = "scaled\nexpression",
  top_annotation = column_anno,
  show_column_names = FALSE
)


## ---- Split Columns by Cell_annotation ----
Heatmap(
  counts_scaled[marker_genes, , drop = FALSE],
  name = "scaled\nexpression",
  top_annotation = column_anno,
  column_split = metadata$Tissue,
  show_column_names = FALSE
)


## ---- Activity ----

#Create a heatmap using counts_scaled with marker_genes.
    # 1. Include Mouse and Tissue as column annotations.
    # 2. Split the columns by Mouse and then Tissue (two seperate heatmaps) and plot them side by side.

# Your code here





# =============== Activity 35: Write a Heatmap Function ===============

# Write a Heatmap Function
make_gene_heatmap <- function(
    genes,
    expression_matrix,
    metadata,
    split_by = "Tissue",
    heatmap_title = "scaled\nexpression") {
  
  # Keep only genes found in the expression matrix
  genes <- intersect(
    genes,
    rownames(expression_matrix)
  )  

  # Create column annotations from metadata
  column_anno <- HeatmapAnnotation(
    Cell_annotation = metadata$Cell_annotation,
    Tissue = metadata$Tissue
  )

  Heatmap(
    expression_matrix[genes, , drop = FALSE], 
    name = heatmap_title,
    top_annotation = column_anno,
    column_split = metadata[[split_by]], # Use metadata[[split_by]] so the split variable can be chosen by the user
    show_column_names = FALSE
  )
}


# Test the Function
make_gene_heatmap(
  genes = marker_genes,
  expression_matrix = counts_scaled,
  metadata = metadata,
  split_by = "Tissue"
)

## ---- Activity ----

# Adjust the function so that it also includes Mouse in the column annotation and 
# instead split the columns by Mouse. Be sure to test your function.


# Your code here






# =============== Activity 36: Loop Over Gene Sets and Save PDFs ===============

# Now we will create several gene sets and save the outputs
gene_sets <- list(
  T_cell_core = c("Cd3d", "Cd3e", "Trac", "Lck"),
  Cytotoxic = c("Nkg7", "Gzmb", "Gzmk", "Prf1", "Ifng"),
  Exhaustion = c("Pdcd1", "Tox", "Lag3", "Havcr2", "Ctla4"),
  Naive_memory = c("Ccr7", "Sell", "Tcf7", "Il7r"),
  Proliferation = c("Mki67", "Top2a", "Stmn1")
)

# Create an Output Folder
dir.create(
  here("results", "heatmaps"),
  recursive = TRUE,
  showWarnings = FALSE
)

## ---- Save Heatmaps in a Loop using the function we created ----

# Loop through each gene set
for (gene_set_name in names(gene_sets)) {

  # Create the heatmap
  heatmap_plot <- make_gene_heatmap(
    genes = gene_sets[[gene_set_name]],
    expression_matrix = counts_scaled,
    metadata = metadata,
    split_by = "Tissue"
  )
  
  # Open a PDF device
  pdf(
    file = here("results", "heatmaps",
                paste0(gene_set_name, "_heatmap.pdf")),
    width = 10,
    height = 6
  )

  # Draw the ComplexHeatmap object
  draw(
    heatmap_plot
  )

  # Close the PDF device
  dev.off()
}

# Check the Files
list.files(
  here("results", "heatmaps"),
  full.names = FALSE
)

## ---- Activity ----

# Change the function and loop to match the below specifications:
      # 1.They are saved in a new folder called heatmaps_activity
      # 2. They plot marker_genes
      # 3. Includes Tissue, Cell_annotation, Mouse, and Phase as column annotations
      # 4. Loops each column splits
      # 5. Saves each heatmap as a PDF

# Your code here





# =============== Activity 37: Box Plots and Violin Plots ===============

# Basic Box Plot
ggplot(
  metadata,
  aes(
    x = Tissue,
    y = nGene
  )
) +
  geom_boxplot()

# Basic Violin Plot
ggplot(
  metadata,
  aes(
    x = Tissue,
    y = nGene
  )
) +
  geom_violin()

# Combine Violin and Box Plot
ggplot(
  metadata,
  aes(
    x = Tissue,
    y = nGene
  )
) +
  geom_violin() +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA
  )

## ---- Activity ----

# 1. Make a violin plot showing percent.mito by Tissue. Add a box plot on top.
# 2. Make a violin plot showing percent.mito by Mouse. Add a box plot on top.


# Your code here



## ---- Add Groupings to Plots ----
ggplot(
  metadata,
  aes(
    x = Tissue,
    y = nGene,
    fill = Batch
  )
) +
  geom_boxplot()

## ---- Activity ----
# Make a violin plot showing nGene by Tissue. 
# Add a box plot on top and fill by Mouse. 
# What happens to the layout?


# Your code here




## ---- Add Individual Cells to Violin and Box Plots ----
ggplot(
  metadata,
  aes(
    x = Tissue,
    y = nGene,
    fill = Tissue
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 0.8,
    alpha = 0.4
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA,
    alpha = 0.8
  ) +
  theme_bw()

# Important: Turn Off Box Plot Outliers When Adding Points
ggplot(
  metadata,
  aes(
    x = Tissue,
    y = percent.mito,
    fill = Tissue
  )
) +
  geom_violin(
    trim = FALSE,
    alpha = 0.7
  ) +
  geom_jitter(
    width = 0.15,
    size = 0.8,
    alpha = 0.4
  ) +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA
  ) +
  theme_bw()

ggplot(
  metadata,
  aes(
    x = Tissue,
    y = nGene,
    fill = Tissue
  )
) +
  geom_jitter(
    aes(colour = Tissue),
    width = 0.18,
    size = 0.8,
    alpha = 0.25,
    show.legend = FALSE
  ) +
  geom_violin(
    alpha = 0.55,
    trim = FALSE,
    width = 0.8
  ) +
  geom_boxplot(
    width = 0.12,
    outlier.shape = NA,
    fill = "white",
    alpha = 0.8
  ) +
  theme_classic() +
  labs(
    x = "Tissue",
    y = "nGene"
  )






# =============== Activity 38: Plot Gene Expression with Metadata ===============

plot_gene <- "Pdcd1"

gene_expression_df <- tibble(
  barcode = colnames(counts_scaled),
  expression = as.numeric(
    counts_scaled[plot_gene, ]
  )
) %>%
  left_join(
    metadata,
    by = "barcode"
  )

gene_expression_df %>%
  select(barcode, expression, Cell_annotation, Tissue) %>%
  head()

## ---- Plot One Gene by Cell_annotation ----
ggplot(
  gene_expression_df,
  aes(
    x = Cell_annotation,
    y = expression
  )
) +
  geom_violin() +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )

## ---- Add Tissue Grouping ----
ggplot(
  gene_expression_df,
  aes(
    x = Tissue,
    y = expression,
    fill = Tissue
  )
) +
  geom_violin() +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA
  ) +
  facet_wrap(
    ~ Cell_annotation
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )



## ---- Activity ----

# Repeat the gene expression plot using Cxcr6 by Tissue. Group by Mouse.

# Your code here







# =============== Activity 39: Reshape Data for Multiple Genes ===============

genes_to_plot <- c(
  "Cd8a",
  "Gzmb",
  "Pdcd1",
  "Tox",
  "Cxcr6"
)

genes_to_plot <- genes_to_plot[
  genes_to_plot %in% rownames(counts_scaled)
]

## ---- Convert Multiple Genes to Long Format ----

multi_gene_expression <- counts_scaled[
  genes_to_plot,
  metadata$barcode
] %>%
  as.data.frame() %>%
  rownames_to_column("gene") %>%
  pivot_longer(
    cols = -gene,
    names_to = "barcode",
    values_to = "expression"
  ) %>%
  left_join(
    metadata %>%
      select(barcode, Cell_annotation, Tissue, Batch),
    by = "barcode"
  )

multi_gene_expression %>%
  head()

## ---- Plot Multiple Genes ----
ggplot(
  multi_gene_expression,
  aes(
    x = Tissue,
    y = expression
  )
) +
  geom_violin() +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA
  ) +
  facet_wrap(
    ~ gene,
    scales = "free_y"
  )

## ---- Plot Multiple Genes by Cell_annotation ----
ggplot(
  multi_gene_expression,
  aes(
    x = Cell_annotation,
    y = expression
  )
) +
  geom_violin() +
  geom_boxplot(
    width = 0.15,
    outlier.shape = NA
  ) +
  facet_wrap(
    ~ gene,
    ncol = 1,
    scales = "free_y"
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  )






# =============== Activity 40: Write a Plotting Function ===============

# Now we will write a reusable function for plotting one gene.
make_gene_violin <- function(
    gene,
    expression_matrix,
    metadata,
    group_var = "Tissue") {

  if (!gene %in% rownames(expression_matrix)) {
    stop("Gene not found in expression matrix: ", gene)
  }

  plot_df <- tibble(
    barcode = colnames(expression_matrix),
    expression = as.numeric(
      expression_matrix[gene, ]
    )
  ) %>%
    left_join(
      metadata,
      by = "barcode"
    )

  ggplot(
    plot_df,
    aes(
      x = .data[[group_var]],
      y = expression
    )
  ) +
    geom_violin() +
    geom_boxplot(
      width = 0.15,
      outlier.shape = NA
    ) +
    labs(
      title = gene,
      x = group_var,
      y = "scaled expression"
    ) +
    theme(
      axis.text.x = element_text(
        angle = 45,
        hjust = 1
      )
    )
}

# Test the Function
make_gene_violin(
  gene = "Pdcd1",
  expression_matrix = counts_scaled,
  metadata = metadata,
  group_var = "Tissue"
)




## ---- Activity ----

# 1. Write a loop to plot a violin plot for the genes in marker genes. 
#    Save them in a folder called violin_plots

# 2. Adjust the function to group the plots by Mouse and write a loop for 
#    the genes in marker genes. Save them in a folder called violin_plots_mouse



# Your code here








# =============== Activity 41: Combine a Heatmap with Box Plots ===============

library(grid)

# Use the marker genes already defined earlier
heatmap_matrix <- counts_scaled[
  marker_genes,
  ,
  drop = FALSE
]

dim(heatmap_matrix)


# Create a box plot annotation for each row/gene
row_boxplot <- rowAnnotation(
  "Expression\nsummary" = anno_boxplot(
    heatmap_matrix,
    outline = FALSE,
    size = unit(2, "mm")
  ),
  width = unit(3, "cm")
)

# Create a column annotation showing cell metadata
column_anno <- HeatmapAnnotation(
  Tissue = metadata$Tissue,
  Cell_annotation = metadata$Cell_annotation,
  show_annotation_name = TRUE,
  show_legend = FALSE
)

# Draw a heatmap with an aligned boxplot for each gene
Heatmap(
  heatmap_matrix,
  name = "scaled\nexpression",
  top_annotation = column_anno,
  show_column_names = FALSE,
  cluster_columns = TRUE,
  cluster_rows = TRUE
) +
  row_boxplot








# =============== Activity 42: Other similar plots ===============

## ---- Ridgeline Plot ----
library(ggridges)

ggplot(
  metadata,
  aes(
    x = nGene,
    y = Cell_annotation,
    fill = Cell_annotation
  )
) +
  geom_density_ridges(
    alpha = 0.7,
    scale = 1.2,
    show.legend = FALSE
  ) +
  theme_bw() +
  labs(
    x = "nGene",
    y = "Cell annotation"
  )


## ---- Dot Plot ----
dotplot_df <- multi_gene_expression %>%
  group_by(gene, Cell_annotation) %>%
  summarise(
    mean_expression = mean(expression),
    percent_expressing = mean(expression > 0) * 100,
    .groups = "drop"
  )

dotplot_df

ggplot(
  dotplot_df,
  aes(
    x = Cell_annotation,
    y = gene
  )
) +
  geom_point(
    aes(
      size = percent_expressing,
      colour = mean_expression
    )
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45,
      hjust = 1
    )
  ) +
  labs(
    x = "Cell_annotation",
    y = "Gene",
    size = "% expressing",
    colour = "Mean scaled expression"
  )
