#!/usr/bin/Rscript
## Author: Mano Ranaweera 
## mranawee@bu.edu
## BU BF591
## Assignment Week 2

#### Bioconductor ####
# it is standard among R packages to define libraries and packages at the 
# beginning of a script. Also note that a package should NOT be installed every 
# time a script runs.
# The bioconductor repository has installation instructions for biomaRt: 
# https://bioconductor.org/install/

if (!require("BiocManager", quietly = TRUE)){
  install.packages("BiocManager")
}
BiocManager::install(version = "3.14")

  
if (!require("biomaRt", quietly = TRUE)){
  install.packages("biomaRt", quietly = TRUE)
}
BiocManager::install("biomaRt", force = TRUE)

install.packages('reshape2')
# load tidyverse and your new bioconductor package
library(tidyverse)
library(BiocManager)
library(ggplot2)
library(dplyr)
library(reshape2)
#### Loading and processing data ####
#' Load Expression Data
#'
#' @param filepath A text string of the full filepath to the file to load.
#'
#' @return A tibble containing the data loaded from the CSV in `filepath`. 
#' 
#' @details Note that not all CSVs are created equal, and there are often cases where 
#' the data will not load in correctly on the first try. You may want to write this functon to 
#' adjust the CSV file being loaded into this assignment so that it can be formed into a 
#' tibble correctly.
#'
#'Adam's example
#'readr::read_delim(filename, delim= ' ')   #read_delim assumes there is a header row
#'    expr_mat <- tibble::as_tibble(t(read.table(
#'       filename,
#'       header = TRUE,
#'       sep = " "
#'       )),
#'       
#'
#'
#'
#' @examples 
#' `data <- load_expression('/project/bf528/project_1/data/example_intensity_data.csv')`
load_expression <- function(filepath) {
  read_data <- read.table('data/example_intensity_data/example_intensity_data.csv', sep = ' ', header = TRUE) 
  tib_data <- as_tibble(read_data)
  return (tib_data)
}
expr <- load_expression('data/example_intensity_data/example_intensity_data.csv')
head(expr)
#' Filter 15% of the gene expression values.
#'
#' @param tibble A tibble of expression values, rows by probe and columns by sample.
#'
#' @return A tibble of affymetrix probe names from the input expression data tibble. 
#' These names match the rows with 15% or more of the expression values about log2(15).
#' 
#' @details This is similar to the filters being implemented in BF528's project 1. 
#' We do not necessarily want to capture all parts of the assay in our analysis, so 
#' filters like this serve to reduce the noise and amount of data to examine.
#'
#' @examples `samples <- filter_15(data_tib)`
#' `> str(samples)`
#' `tibble [40,158 × 1] (S3: tbl_df/tbl/data.frame)`
#' `$ probeids: chr [1:40158] "1007_s_at" "1053_at" "117_at" "121_at" ...`
#' 
#' 
#' only rows with 15% expression values above log_2(15)
filter_15 <- function(tibble){
  tibble$sig_vals <- rowSums(tibble[-1] > (log2(15))) #values per row > log2(15)...counts of this condition per row recorded in extra column "sig_vals"
  filtered_tib <- subset(tibble, ((sig_vals / length(select_if(tibble, is.numeric))) >= 0.15)) #at least 15 % of values meet threshold per row
  return(filtered_tib[,1])
}
#length(select_if(expr, is.numeric))
filtered <- filter_15(expr) #actual filtered data
head(filtered)
#### Gene name conversion ####

#' Convert affymetrix array names into hgnc_symbol IDs using biomaRt. Inputs and 
#' outputs will likely not be the same size.
#'
#' @param affy_tib A single column tibble of strings containing array names.
#'
#' @return A 2 column tibble that contains affy IDs in the first column,
#' and their corresponding HGNC gene ID in the second column. Note that not all affy IDs 
#' will necessarily correspond with a gene ID, and one gene may have multiple affy IDs.
#' 
#' @details Connecting to ensembl via biomaRt can be...hit or miss...so you may 
#' want to check if data was correctly returned (or if it was just empty). The 
#' `getBM()` function may not accept a tibble, so you might need to convert your 
#' input into a flat vector.
#'
#' @examples 
#' `> affy_to_hgnc(tibble(c('202860_at', '1553551_s_at')))`
#' `affy_hg_u133_plus_2 hgnc_symbol`
#' `1        1553551_s_at      MT-ND1`
#' `2        1553551_s_at       MT-TI`
#' `3        1553551_s_at       MT-TM`
#' `4        1553551_s_at      MT-ND2`
#' `5           202860_at     DENND4B`
#' 
#' find HGNC gene ID based on our probe ID
affy_to_hgnc <- function(affy_vector) {
  affyID <- pull(affy_vector['probeids'])
  ensembl <- useEnsembl(biomart="ensembl", dataset='hsapiens_gene_ensembl')
  gene_map <- as_tibble(
    getBM(
      attributes=c('affy_hg_u133_plus_2', 'hgnc_symbol'),
      filters = 'affy_hg_u133_plus_2',
      values = affyID,
      mart= ensembl
    )
  )
  return(gene_map)
}

#affy_vector <- filtered$probeids
#str_probe <- toString(affy_vector) #param
vector <- unlist(filtered$probeids)
sample_names <- affy_to_hgnc(filtered)
sample_names_test <- affy_to_hgnc(filtered)

#sample_names2 <- sample_names %>% drop_na()
head(sample_names_test)
#### ggplot ####

#' Reduce a tibble of expression data to only the rows in good_genes or bad_genes.
#'
#' @param expr_tibble A tibble holding the expression data, each row corresponding to
#' one affymetrix probe ID and each column to a sample.
#' @param names_ids A two column tibble that associates affy IDs with HGNC gene IDs. 
#' Generated `with affy_to_hgnc()`.
#' @param good_genes A list of gene names stored as a vector of strings.
#' @param bad_genes A list of gene names stored as a vector of strings.
#'
#' @return A tibble with two additional columns added:
#' 1. HGNC gene IDs 
#' 2. Does the gene is this row fall into "good" or "bad" genes?
#' This tibble should be reduced to only rows present in good or bad genes. All
#' other rows can be discarded.
#' 
#' @details In order to plot only our genes of interest, we need to rearrange our 
#' data to include only the elements we want to see. We also want to add to columns, 
#' one that associates the probeids with the HGNC gene name, and one that says if 
#' that gene is in the good or bad sets of genes.
#'
#' @examples 
#' `plot_tibble <- reduce_data(expr_tibble = expr, names_ids = sample_names,`
#' `                           goodGenes, badGenes)`
#' `> head(plot_tibble)`
#' `A tibble: 6 × 38`
#' `  probeids    hgnc    gene_set    GSM972389 ...`
#' `  <chr>       <chr>   <chr>       <dbl>     ...`
#' `1 202860_at   DENND4B good        7.16      ...`
#' `2 204340_at   TMEM187 good        6.40      ...`

#good <- 

reduce_data <- function(names_ids, expr_tibble, good_genes, bad_genes){
  #hgnc_symbols <- names_ids
  names(names_ids) <- c('probeids', 'hgnc_symbol')
  reduced_tib <- names_ids %>%
    filter(hgnc_symbol %in% good_genes | hgnc_symbol %in% bad_genes)
  #gene_set column
  reduced_tib <- reduced_tib %>%
    mutate(gene_set = if_else(hgnc_symbol %in% good_genes, 'good', 'bad'))
    
  #if (names_ids$hgnc_symbol == good_genes){
  #  reduced_tib$gene_set = 'good'
  #} else {
  #  reduced_tib$gene_set = 'bad'
  #}
  
  sample_tib <- expr_tibble %>%
    filter(probeids %in% reduced_tib$probeids)
  
  
  return(as_tibble(merge(reduced_tib, sample_tib)))

}

bad <- c("TP53", "EGFR", "BRAF", "KRAS", "PIK3CA", "ERBB2", "MAPK1", "NRAS")
good <- c("PKD1", "NOS3", "AGTR1", "COL4A5", "ECE1", "MEN1", "OLR1", "F7")

plot_tibble <- reduce_data(sample_names, expr, good, bad)

head(plot_tibble)


#new_expr_tibble <- affy_to_hgnc(expr)
new_expr_tibble_tib <- merge(sample_names_test, expr)

names(new_expr_tibble_tib) <- c('probeids', 'hgnc_symbol')

new_expr_tibble_tib %>%
  filter(hgnc_symbol %in% bad | hgnc_symbol %in% good)



plot_tibble <- reduce_data(expr_tibble = expr,
                           names_ids = sample_names,
                           good,
                           bad)

head(plot_tibble)

#' Plot a boxplot of good and bad genes.
#'
#' @param tibble A reduced tibble of expression data, with information about
#' good and bad genes and gene names.
#'
#' @return A ggplot object which contains a boxplot of the genes and samples we 
#' are interested in.
#' 
#' @details This function performs one additional step before using `ggplot()`: 
#' converting the _wide_ format of the input tibble to a _long_ format.
#'
#' @examples `p <- plot_ggplot(plot_tibble)`
plot_ggplot <- function(tibble) {
 tibble %>% ggplot(aes(x=hgnc_symbol, y = (GSM972409))) +
    geom_boxplot() +
    theme_bw(base_size = 16)
    
}

plot_ggplot(plot_tibble)
