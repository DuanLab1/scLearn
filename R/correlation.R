#' @title Calculate Sample Similarity Matrix
#'
#' @description
#' Computes pairwise similarity/distance matrix between samples (columns) using various
#' metrics with parallel computation support. Suitable for large-scale genomic data.
#'
#' @details
#' Available similarity/distance measures:
#' \itemize{
#'   \item Pearson correlation (linear relationship)
#'   \item Spearman correlation (rank-based relationship)
#'   \item Cosine similarity (angle between vectors)
#'   \item Euclidean distance (geometric distance)
#' }
#' The function utilizes parallel computation via \code{parallel} package to accelerate
#' calculations for large matrices.
#'
#' @author Bin Duan (binduan\@sjtu.edu.cn)
#'
#' @param matrix Numeric matrix where columns represent samples and rows represent features
#' @param method Similarity/distance measure: "pearson", "spearman", "cosin", or
#'        "euclidean" (default: "pearson")
#' @param cpu_num Number of CPU cores to use for parallel computation (default: 8)
#'
#' @return A symmetric similarity/distance matrix with dimensions ncol(matrix) x ncol(matrix),
#' where row and column names match the input matrix column names.
#'
#' @importFrom parallel makeCluster stopCluster parApply
#' @importFrom stats cor
#'
#' @examples
#' \dontrun{
#' # Using example data from scLearn package
#' data(QueryCellData)
#'
#' # Get normalized expression matrix
#' norm_expr <- logcounts(QueryCellData)
#'
#' # Calculate Pearson correlation
#' cor_mat <- correlation(
#'   matrix = norm_expr,
#'   method = "pearson",
#'   cpu_num = 4
#' )
#'
#' # Calculate cosine similarity
#' cos_mat <- correlation(
#'   matrix = norm_expr,
#'   method = "cosin",
#'   cpu_num = 2
#' )
#'
#' # Visualize results
#' heatmap(cor_mat, symm = TRUE)
#' }
#'
#' @seealso
#' \code{\link[stats]{cor}} for correlation calculations,
#' \code{\link[parallel]{makeCluster}} for parallel computation setup
#'
#' @export
#'
correlation <- function(
    matrix,
    method = c("pearson", "spearman", "cosin", "euclidean"),
    cpu_num = 8) {

  cosdist <- function(x1, x2) {
    n1 <- sqrt(sum(x1^2))
    n2 <- sqrt(sum(x2^2))
    d <- as.numeric(x1 %*% x2) / n1 / n2
    return(d)
  }

  euclidean <- function(x1, x2) {
    return(sqrt(t(x1 - x2) %*% (x1 - x2)))
  }

  method <- method[1]

  if (!(method %in% c("pearson", "spearman", "cosin", "euclidean"))) {
    stop("method must be 'pearson','spearman','cosin' or 'euclidean'.")
  }

  cpu_num_set <- makeCluster(getOption("cluster.cores", cpu_num))
  sample_num <- ncol(matrix)
  cor_matrix <- matrix(rep(1, sample_num^2), sample_num)
  colnames(cor_matrix) <- colnames(matrix)
  row.names(cor_matrix) <- colnames(matrix)
  simi_cor <- function(vec, matrix, method) {
    if (method == "cosin") {
      cor_matrix <- apply(matrix, 2, cosdist, x2 = vec)
    } else if (method == "euclidean") {
      cor_matrix <- apply(matrix, 2, euclidean, x2 = vec)
    } else {
      cor_matrix <- apply(matrix, 2, cor, y = vec, method = method)
    }
    return(cor_matrix)
  }

  cor_matrix <- parApply(cpu_num_set, matrix, 2, simi_cor, matrix = matrix, method = method)
  stopCluster(cpu_num_set)

  return(cor_matrix)
}
