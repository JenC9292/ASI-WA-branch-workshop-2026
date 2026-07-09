# libraries

library(tidyverse)
library(ggside)
library(here)

# Data reading

GSE168944_matrix_path <- here("data", "GSE168944_cell.counts.txt.gz")
n <- 10000

GSE168944_metadata_path <- here("data", "GSE168944_cell.metadata.csv.gz")

theme_set(theme_minimal(base_size = 14))

Sys.setenv("VROOM_CONNECTION_SIZE" = 1e7)


header <- read_lines(GSE168944_matrix_path, n_max = 1)

cell_ids <- str_split(header, pattern = " ")[[1]]
cell_ids <- cell_ids[cell_ids != ""]

final_header <- c("gene", cell_ids)


set.seed(42)

sampled_cells <- sample(cell_ids, size = n)

selected_cols <- c("gene", sampled_cells)

counts_matrix <- read_delim(
  file = GSE168944_matrix_path,
  delim = " ",
  col_names = final_header,
  skip = 1,
  col_select = all_of(selected_cols)
)

counts_matrix <-
  counts_matrix %>%
  rename_with(~ gsub('"', '', .x))

GSE168944_metadata <- read_delim(
  file = GSE168944_metadata_path,
)

# save the main two files we will use
saveRDS(counts_matrix, here("data", "counts_matrix.rds"))
saveRDS(GSE168944_metadata, here("data", "GSE168944_metadata.rds"))

# remove unused objects from the environment
rm(cell_ids, final_header, GSE168944_matrix_path, GSE168944_metadata_path,
   header, n, sampled_cells, selected_cols)

class(counts_matrix)

class(GSE168944_metadata)

typeof(counts_matrix)

typeof(GSE168944_metadata)


#

## Activity

#Determine:
  
#-   Number of genes
#-   Number of cells
#-   Number of metadata variables

#

# create and expression matrix
expr_mat <- counts_matrix %>%
  column_to_rownames("gene") %>%
  as.matrix()

dim(expr_mat)

# Transpose the matrix
expr_cells <- t(expr_mat)

dim(expr_cells)

# What do rows and columns represent now?


## Activity

#Inspect the metadata. Please note glimpse will exhaustively print all columns - hence the subset.

## Summarise Categories

GSE168944_metadata %>%
  count(Cluster)

GSE168944_metadata %>%
  count(Tissue)

GSE168944_metadata %>%
  count(Phase)

# Activity 4: Matching Cells Between Datasets

#Metadata and expression matrices must have matching cell names.

  
## Match Metadata and Counts
  
metadata <- GSE168944_metadata %>%
  filter(barcode %in% rownames(expr_cells))

expr_cells <- expr_cells[
  metadata$barcode,
]

all(
  rownames(expr_cells) ==
    metadata$barcode
)

## Depending on how you need the data, a left or right join will also remove unneeded data.

expr_mat_meta <- as_tibble(expr_cells, rownames = "barcode") %>%
  left_join(GSE168944_metadata, by = "barcode")

expr_mat_meta %>% select((last_col() - 5):last_col())

# Activity 5: Working with Tibbles

## Activity

#Determine:
  
#-   The number of rows
#-   The number of columns
#-   The number of missing values

# Activity 6: Exploring QC Metrics

## Start

#Common quality metrics include:
  
#-   nGene
#-   nUMI
#-   percent.mito


## Calculate Summary Statistics
  

metadata %>%
  summarise(
    mean_genes = mean(nGene),
    mean_umi = mean(nUMI),
    mean_mito = mean(percent.mito)
  )
  
## Visualise QC Metrics
  

ggplot(
  metadata,
  aes(nGene)
) +
  geom_histogram(
    bins = 50
  )

ggplot(
  metadata,
  aes(percent.mito)
) +
  geom_histogram(
    bins = 50
  )

# Activity 7: Conditional Logic

## Start

#Conditional logic is the 'backbone' of many different analysis and coding challenges.

  
## Worked Example - if we have mito at 12% but 10% is our cutoff
  

mito <- 0.12

if(mito > 0.1){
  
  print("Low quality cell")
  
} else {
  
  print("High quality cell")
}

## Activity
  
#Create a QC classification variable. Note this is an already cleaned dataset so we are setting the percent mito very low.

# Activity 8: Vectorised Logic

## Activity

#Create categories for all cells simultaneously. This is basically replacing chaining if_else statements - reads cleaner.

#~denotes the formula where the left side is the logical condition and the rightside is the replacement.

# Activity 9: Functions

## Start

#Functions allow code reuse. This is important and much better than copy
#and pasting the same code a million times! Below is the formula for the
#percentage coefficient of variation (CV), you don't want to have to
#write this out every time!


## Worked Example


cv <- function(x){

  sd(x, na.rm = TRUE) /
    mean(x, na.rm = TRUE) *
    100

}

## Activity

#Calculate the CV for metadata variables. Here we call cv() with the x
#input being a numeric vector. We access this from the metadata by
#accessing the named list using tibble\$named_list.



# Activity 10: Functions Returning Tables

## Create Summary Function

#This can be quite useful for example, many times when you are reporting
#say a metric such as the difference between two groups, you might want
#to able to quickly summarise the data.

summary_stats <- function(x){
  
  tibble(
    mean = mean(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    min = min(x, na.rm = TRUE),
    max = max(x, na.rm = TRUE)
  )
  
}

## Test the Function
  

summary_stats(metadata$nGene)

summary_stats(metadata$nUMI)

# Activity 11: Loops
  
## Start
  
#Loops automate repetitive calculations. For example, you want to loop
#through and perform a complex calculation on each row of a dataframe,
#that can't be vectorised.


## Simple Loop

for(i in 1:5){

  print(i)

}


## Loop Through Variables

vars <- c(
  "nGene",
  "nUMI",
  "percent.mito"
)

for(v in vars){

  print(v)

  print(
    mean(
      metadata[[v]],
      na.rm = TRUE
    )
  )

}

# Activity 12: Storing Results from Loops

#Depending on what you need, most of the time you will need to store data from a loop.
#For example, we use a loop to calculate the mean of the variables.
#This can be achieved in a number of ways, but below we init an empty vector of length vars.
#The loop then iterates along vars, and stores the means the the vector slots.

results <- numeric(
  length(vars)
)

for(i in seq_along(vars)){
  
  results[i] <-
    mean(
      metadata[[vars[i]]],
      na.rm = TRUE
    )
  
}

names(results) <- vars

results

# Activity 13: Loops versus purrr

## purrr Solution

#Overall, I would say Purrr map_x supersedes writing a manual loop. But
#good to know how it works if you need more granular control, or you need
#to maximise memory.

#There are also sapply, lapply, vapply, and apply in base R which are also quite useful. Purrr is the tidy version of these.
#For example, lapply will apply a function to each element of a list. Note, dataframes are lists of columns.

map_dbl(
  metadata[vars],
  mean,
  na.rm = TRUE
)

# Activity 15: Select Highly Variable Genes

## Start

#Using all genes can introduce noise.

#We first identify highly variable genes.

  
## Calculate Gene Variance
  
gene_var <- apply(
  expr_cells,
  2,
  var
)

top_genes <- names(
  sort(
    gene_var,
    decreasing = TRUE
  )[1:2000]
)

  
## Subset Variable Genes
  
expr_hvg <- expr_cells[
  ,
  top_genes
]

dim(expr_hvg)

  
# Activity 16: Log Transformation
  
## Start
  
#Raw counts are highly skewed.

#A log transformation stabilises variance.


## Transform Counts

expr_log <- log1p(
  expr_hvg
)

  
# Activity 17: Running PCA
  
## Activity
  
#Run PCA.


### Solution

sc_pca <- prcomp(
  expr_log,
  center = TRUE,
  scale. = TRUE
)

summary(sc_pca)

  
# Activity 18: Variance Explained
## Start
  
#PC1 explains the greatest amount of variation.

#PC2 explains the next largest amount.

#PCn should be capturing different pattern not in other PCs.

  
## Interpretation Questions
  
#  1.  How much variance does PC1 explain?
#  2.  How much variance do PC1 and PC2 explain?
#  3.  How many PCs explain more than 80% of the variance?
  

# Activity 19: PCA Scores
  
## Start
  
#Scores describe where cells sit in PCA space. We can extract the
#specific cells PCA score (fun fact - you can use this in machine
#                          learning as a method to reduce colinearity of variables and to reduce
#                          dimensions).

  
 ## Extract Scores
  
scores <- as_tibble(
  sc_pca$x,
  rownames = "barcode"
)

scores <- scores %>%
  left_join(
    metadata,
    by = "barcode"
  )

head(scores)

saveRDS(scores, here("data", "PCA_scores.rds"))

  
# Activity 20: PCA by Cluster
  
ggplot(
  scores,
  aes(
    PC1,
    PC2,
    colour = Cluster
  )
) +
  geom_point(
    alpha = 0.7,
    size = 1
  )

  ## Interpretation Questions
  
#  1.  Do clusters separate?
#  2.  Are there outliers?
#  3.  Which clusters overlap?
  

  
# Activity 21: PCA by Tissue
  
ggplot(
  scores,
  aes(
    PC1,
    PC2,
    colour = Tissue
  )
) +
  geom_point(
    alpha = 0.7,
    size = 1
  )


# Activity 22: PCA by Batch
  
ggplot(
  scores,
  aes(
    PC1,
    PC2,
    colour = Batch
  )
) +
  geom_point(
    alpha = 0.7,
    size = 1
  )

  
# Activity 23: PCA Loadings
  
## Start
  
# Loadings identify which genes drive each principal component.


## Extract Loadings
  
loadings <- as_tibble(
  sc_pca$rotation,
  rownames = "gene"
)

loadings

  
# Activity 24: Genes Driving PC1
  
loadings %>%
  mutate(
    abs_PC1 = abs(PC1)
  ) %>%
  arrange(
    desc(abs_PC1)
  ) %>%
  slice(1:20) %>%
  select(gene, PC1, abs_PC1)

# Activity 25: Visualising Loadings
  
## Plot Contributions to PC1
  
loadings %>%
  mutate(
    abs_PC1 = abs(PC1)
  ) %>%
  slice_max(
    abs_PC1,
    n = 20
  ) %>%
  ggplot(
    aes(
      reorder(gene, PC1),
      PC1
    )
  ) +
  geom_col() +
  coord_flip()
  
# Activity 26: PCA with ggside 
  
#ggside is a useful tool for adding side plots to other graphs. This is useful,
#for example in PCA a KDE can more elegantly describe the distribution differences
#between PC1 and PC2.


## Plot KDE for the PCA highlighting tissue


ggplot(
  scores,
  aes(
    PC1,
    PC2,
    colour = Tissue,
    fill = Tissue
  )
) +
  geom_point(
    alpha = 0.7,
    size = 1
  ) +
  geom_xsidedensity(
    alpha = 0.3
  ) +
  geom_ysidedensity(
    alpha = 0.3
  ) +
  theme_bw() +
  theme(
    ggside.panel.scale = 0.3
  )


## Plot boxplot for the PCA highlighting tissue


ggplot(
  scores,
  aes(
    PC1,
    PC2,
    colour = Tissue,
    fill = Tissue
  )
) +
  geom_point(
    alpha = 0.7,
    size = 1
  ) +
  geom_xsideboxplot(
    aes(
      y = Tissue
    ),
    alpha = 0.5,
    orientation = "y"
  ) +
  geom_ysideboxplot(
    aes(
      x = Tissue
    ),
    alpha = 0.5,
    orientation = "x"
  ) +
  theme_bw() +
  theme(
    ggside.panel.scale = 0.3
  )

