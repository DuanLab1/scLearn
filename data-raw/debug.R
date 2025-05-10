# build package
roxygen2::roxygenise(getwd())
devtools::check(document = FALSE)
devtools::build(binary = FALSE, manual = TRUE, quiet = FALSE)

# build pkgdown
pkgdown::build_site() # Run to build the website

library(usethis)
library(SingleCellExperiment)

RefCellData <- readRDS(system.file("extdata/baron-human.rds", package = "scLearn"))

set.seed(123)
cell_types <- colData(RefCellData)$cell_type1
stratified_sample <- function(cell_types, n_total = 1000, min_per_type = 10) {
  sampled <- lapply(split(seq_along(cell_types), cell_types), function(idx) {
    if(length(idx) > min_per_type) {
      sample(idx, max(min_per_type, round(n_total * length(idx)/length(cell_types))))
    } else {
      idx
    }
  })
   unlist(sampled)
}

sampled_cells <- stratified_sample(cell_types)
RefCellData <- RefCellData[, sampled_cells]
RefCellData

QueryCellData <- readRDS(system.file("extdata/muraro-human.rds", package = "scLearn"))
cell_types <- colData(QueryCellData)$cell_type1
sampled_cells <- stratified_sample(cell_types)
QueryCellData <- QueryCellData[, sampled_cells]

use_data(RefCellData, overwrite = TRUE)
use_data(QueryCellData, overwrite = TRUE)
