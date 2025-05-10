#' @title Single-Cell RNA-Seq Data Quality Control
#'
#' @description
#' Performs quality control filtering on single-cell RNA-seq expression data based on:
#' \itemize{
#'   \item Gene counts per cell (detected features)
#'   \item UMI counts per cell (library size)
#'   \item Mitochondrial gene percentage
#' }
#' The function applies user-defined thresholds to filter low-quality cells and optionally
#' generates diagnostic plots of QC metrics distributions.
#'
#' @details
#' Key processing steps:
#' \itemize{
#'   \item Calculates mitochondrial gene percentage (species-aware)
#'   \item Normalizes counts to 10,000 counts per cell (CP10K)
#'   \item Applies logarithmic transformation (optional)
#'   \item Filters cells based on user-defined thresholds
#'   \item Removes mitochondrial genes from final output
#'   \item Preserves sample metadata when provided
#' }
#'
#' @author Bin Duan (binduan\@sjtu.edu.cn)
#'
#' @param expression_profile Raw count matrix (genes x cells)
#' @param sample_information_cellType Optional cell type annotations (named vector)
#' @param sample_information_timePoint Optional time point annotations (named vector)
#' @param species Species specification ("Hs" for human, "Mm" for mouse)
#' @param gene_low Minimum genes required per cell (default: 500)
#' @param gene_high Maximum genes allowed per cell (default: 10000)
#' @param mito_high Maximum mitochondrial percentage allowed (default: 0.1)
#' @param umi_low Minimum UMIs required per cell (default: 1500)
#' @param umi_high Maximum UMIs allowed per cell (default: Inf)
#' @param logNormalize Whether to apply log transformation (default: TRUE)
#' @param plot Whether to generate QC diagnostic plots (default: FALSE)
#' @param plot_path Output path for QC plots (default: "./quality_control.pdf")
#'
#' @return A list containing:
#' \itemize{
#'   \item expression_profile: Filtered and normalized expression matrix
#'   \item sample_information_cellType: Filtered cell type annotations (if provided)
#'   \item sample_information_timePoint: Filtered time point annotations (if provided)
#' }
#'
#' @importFrom stringr str_detect
#' @importFrom grDevices pdf dev.off
#' @importFrom graphics hist lines par
#' @importFrom stats density
#' @import SingleCellExperiment
#'
#' @examples
#' \dontrun{
#' # Load example data
#' library(SingleCellExperiment)
#' data(QueryCellData)
#'
#' # Extract expression matrix (assuming counts are in 'counts' assay)
#' counts <- assay(QueryCellData, "counts")
#'
#' # Run QC with default parameters
#' qc_results <- Cell_qc(
#'   expression_profile = counts,
#'   species = "Hs"
#' )
#'
#' # Run QC with custom thresholds and plotting
#' qc_results <- Cell_qc(
#'   expression_profile = counts,
#'   species = "Hs",
#'   gene_low = 600,
#'   mito_high = 0.2,
#'   plot = TRUE,
#'   plot_path = "qc_plots.pdf"
#' )
#'
#' # With cell type annotations
#' cell_types <- colData(QueryCellData)$cell_type1
#' qc_results <- Cell_qc(
#'   expression_profile = counts,
#'   sample_information_cellType = cell_types,
#'   species = "Hs"
#' )
#' }
#'
#' @export
#'
Cell_qc <- function(
    expression_profile,
    sample_information_cellType = NULL,
    sample_information_timePoint = NULL,
    species = "Hs",
    gene_low = 500,
    gene_high = 10000,
    mito_high = 0.1,
    umi_low = 1500,
    umi_high = Inf,
    logNormalize = TRUE,
    plot = FALSE,
    plot_path = "./quality_control.pdf") {


  if (species == "Hs") {
    mito.genes <- grep("^MT-", rownames(expression_profile), value = FALSE)
  } else if (species == "Mm") {
    mito.genes <- grep("^mt-", rownames(expression_profile), value = FALSE)
  } else {
    stop("species should be 'Mm' or 'Hs'")
  }
  mito_percent <- function(x, mito.genes) {
    return(sum(x[mito.genes]) / sum(x))
  }
  percent.mito <- apply(expression_profile, 2, mito_percent, mito.genes = mito.genes)
  nUMI <- apply(expression_profile, 2, sum)
  nGene <- apply(expression_profile, 2, function(x) {
    length(x[x > 0])
  })
  expression_profile <- apply(expression_profile, 2, function(x) {
    x / (sum(x) / 10000)
  })
  if (logNormalize == TRUE) {
    expression_profile <- log(expression_profile + 1)
  }
  if (species == "Hs") {
    expression_profile <- expression_profile[!str_detect(row.names(expression_profile), "^MT-"), ]
  } else if (species == "Mm") {
    expression_profile <- expression_profile[!str_detect(row.names(expression_profile), "^mt-"), ]
  }
  filter_retain <- rep("filter", ncol(expression_profile))
  for (i in 1:length(filter_retain)) {
    if (nGene[i] > gene_low & nGene[i] < gene_high & percent.mito[i] < mito_high & nUMI[i] > umi_low & nUMI[i] < umi_high) {
      filter_retain[i] <- "retain"
    }
  }
  if (plot) {
    pdf(file = plot_path)
    par(mfrow = c(1, 3))
    hist(nGene, breaks = 20, freq = FALSE, xlab = "Gene numbers", ylab = "Density", main = "Gene numbers distribution")
    lines(nGene, col = "red", lwd = 1)
    hist(nUMI, breaks = 20, freq = FALSE, xlab = "UMI numbers", ylab = "Density", main = "UMI numbers distribution")
    lines(density(nUMI), col = "red", lwd = 1)
    hist(percent.mito, breaks = 20, freq = FALSE, xlab = "Percent of mito", ylab = "Density", main = "Percent of mito distribution")
    lines(density(percent.mito), col = "red", lwd = 1)
    dev.off()
  }
  SQ_filter <- as.matrix(expression_profile[, which(filter_retain == "retain")])
  SQ_data_qc <- SQ_filter # have adopted total_expr=10000 and log
  if (is.null(sample_information_cellType)) {
    if (is.null(sample_information_timePoint)) {
      return(list("expression_profile" = SQ_data_qc))
    }
  }
  if (is.null(sample_information_timePoint)) {
    sample_information_qc <- sample_information_cellType[colnames(SQ_data_qc)]
    return(list("expression_profile" = SQ_data_qc, "sample_information_cellType" = sample_information_qc))
  } else {
    sample_information_qc <- gsub("_", "-", sample_information_cellType[colnames(SQ_data_qc)])
    sample_information2_qc <- gsub("_", "-", sample_information_timePoint[colnames(SQ_data_qc)])
    return(list("expression_profile" = SQ_data_qc, "sample_information_cellType" = sample_information_qc, "sample_information_timePoint" = sample_information2_qc))
  }
}
