#' @title Feature Selection Using M3Drop
#'
#' @description
#' Identifies highly variable genes using the M3Drop package, which implements
#' Michaelis-Menten modeling of dropouts for feature selection in single-cell RNA-seq data.
#'
#' @details
#' This function wraps the M3DropFeatureSelection method which:
#' \itemize{
#'   \item Models the relationship between gene expression mean and variance
#'   \item Identifies genes with significantly more dropouts than expected
#'   \item Uses false discovery rate (FDR) correction for multiple testing
#' }
#' The function automatically handles both log-normalized and raw count data.
#'
#' @author Bin Duan (binduan\@sjtu.edu.cn)
#'
#' @param expression_profile Gene expression matrix (genes x cells)
#' @param log_normalized Whether the input matrix is log-normalized (default: TRUE)
#' @param threshold FDR threshold for significant features (default: 0.05)
#'
#' @return A character vector of selected gene names that pass the significance threshold.
#'
#' @importFrom M3Drop M3DropFeatureSelection
#'
#' @examples
#' \dontrun{
#' # Using example data from scLearn package
#' data(QueryCellData)
#'
#' # Using log-normalized data (default)
#' lognorm_data <- logcounts(QueryCellData)
#' selected_genes <- Feature_selection_M3Drop(
#'   expression_profile = lognorm_data,
#'   threshold = 0.01
#' )
#'
#' # Using raw counts
#' raw_data <- counts(QueryCellData)
#' selected_genes <- Feature_selection_M3Drop(
#'   expression_profile = raw_data,
#'   log_normalized = FALSE,
#'   threshold = 0.1
#' )
#'
#' # Examine selected features
#' head(selected_genes)
#' length(selected_genes)
#' }
#'
#' @references
#' Andrews TS, Hemberg M (2018). "M3Drop: dropout-based feature selection
#' for scRNASeq." Bioinformatics, 35(16), 2865-2867.
#'
#' @seealso
#' \code{\link[M3Drop]{M3DropFeatureSelection}} for the underlying implementation
#'
#' @export
#'
Feature_selection_M3Drop <- function(
    expression_profile,
    log_normalized = TRUE,
    threshold = 0.05) {

  if (log_normalized) {
    high_varGenes <- M3DropFeatureSelection(exp(expression_profile) - 1, mt_method = "fdr", mt_threshold = threshold)
    high_varGene_names <- row.names(high_varGenes)

    return(high_varGene_names)
  } else {
    high_varGenes <- M3DropFeatureSelection(expression_profile, mt_method = "fdr", mt_threshold = threshold)
    high_varGene_names <- row.names(high_varGenes)

    return(high_varGene_names)
  }
}
