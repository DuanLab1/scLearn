#' @title Visualize Clusters Using t-SNE or UMAP
#'
#' @description
#' Generates 2D visualization of high-dimensional data using either t-SNE or UMAP dimensionality
#' reduction, with options for cluster labeling and aesthetic customization.
#'
#' @details
#' Key features:
#' \itemize{
#'   \item Supports both t-SNE (t-Distributed Stochastic Neighbor Embedding) and UMAP (Uniform Manifold Approximation and Projection)
#'   \item Option to use pre-computed coordinates or calculate new embeddings
#'   \item Automatic cluster labeling at median positions
#'   \item Customizable plotting parameters and reproducible results via seed setting
#' }
#' The function returns both the plot object and the coordinates for further analysis.
#'
#' @author Bin Duan (binduan\@sjtu.edu.cn)
#'
#' @param data Input data matrix (features x samples) or pre-computed coordinates if calculated=FALSE
#' @param label Vector of cluster labels for each sample (optional)
#' @param point_size Size of points in the plot (default: 1)
#' @param method Dimensionality reduction method: "tsne" or "umap" (default: "tsne")
#' @param draw_cluster_text Whether to display cluster labels at median positions (default: TRUE)
#' @param calculated Whether to compute dimensionality reduction (TRUE) or use provided coordinates (FALSE)) (default: TRUE)
#' @param pca Whether to perform PCA preprocessing for t-SNE (default: TRUE)
#' @param perplexity t-SNE perplexity parameter (default: 100)
#' @param plot Whether to display the plot (default: TRUE)
#' @param seed Random seed for reproducibility (default: 1)
#'
#' @return A list containing:
#' \itemize{
#'   \item p: ggplot2 object of the visualization
#'   \item x: Data frame with coordinates (V1 and V2 columns)
#'   \item cell_group: Factor vector of cluster labels
#' }
#'
#' @importFrom ggplot2 ggplot aes geom_point xlab ylab annotate theme_bw theme element_text guide_legend guides
#' @importFrom dplyr group_by summarise
#' @importFrom stats median
#'
#' @examples
#' \dontrun{
#' # Using example data from scLearn package
#' data(QueryCellData)
#'
#' # Get normalized expression matrix
#' norm_expr <- logcounts(QueryCellData)
#' cell_types <- colData(QueryCellData)$cell_type1
#'
#' # t-SNE visualization
#' tsne_res <- DrawCluster(
#'   data = norm_expr,
#'   label = cell_types,
#'   method = "tsne",
#'   perplexity = 30,
#'   point_size = 1.5
#' )
#'
#' # UMAP visualization without cluster labels
#' umap_res <- DrawCluster(
#'   data = norm_expr,
#'   label = cell_types,
#'   method = "umap",
#'   draw_cluster_text = FALSE,
#'   point_size = 2
#' )
#'
#' # Using pre-computed coordinates
#' precomputed_coords <- data.frame(V1 = rnorm(ncol(norm_expr)),
#'                                V2 = rnorm(ncol(norm_expr)))
#' DrawCluster(
#'   data = precomputed_coords,
#'   label = cell_types,
#'   calculated = FALSE
#' )
#' }
#'
#' @seealso
#' \code{\link[Rtsne]{Rtsne}} for t-SNE implementation,
#' \code{\link[umap]{umap}} for UMAP implementation
#'
#' @export
#'
DrawCluster <- function(
    data,
    label = NULL,
    point_size = 1,
    method = c("tsne", "umap"),
    draw_cluster_text = TRUE,
    calculated = TRUE,
    pca = TRUE,
    perplexity = 100,
    plot = TRUE,
    seed = 1) {


  set.seed((seed))
  method <- method[1]

  if (calculated) {
    if (method == "tsne") {
      tsneresult2 <- Rtsne::Rtsne(t(data), perplexity = perplexity, pca = pca)
      X <- as.data.frame(tsneresult2$Y)
    } else if (method == "umap") {
      umapresult1 <- umap::umap(t(data))
      X <- as.data.frame(umapresult1$layout)
    } else {
      stop("method must be tsne or umap.")
    }
  } else {
    X <- data
  }
  if (length(label) == 0) {
    label <- array(1, dim(X)[1])
    labelname <- c(1)
  }
  labelname <- names(table(label))
  p <- ggplot(X, aes(x = X[, 1], y = X[, 2]))
  cell_group <- factor(label)
  if (method == "tsne") {
    p <- p + geom_point(aes(color = cell_group), size = point_size) + xlab("tSNE1") + ylab("tSNE2")
  } else {
    p <- p + geom_point(aes(color = cell_group), size = point_size) + xlab("umap1") + ylab("umap2")
  }
  if (draw_cluster_text) {
    Label_cal <- X
    Label_cal$cluster <- label
    cluster_x_y <- Label_cal %>%
      group_by(cluster) %>%
      summarise(x_median = median(V1), y_median = median(V2))
    p <- p + annotate("text", x = cluster_x_y$x_median, y = cluster_x_y$y_median, label = cluster_x_y$cluster)
  }
  mytheme <- theme_bw() +
    theme(
      plot.title = element_text(size = rel(1.5), hjust = 0.5),
      axis.title = element_text(size = rel(1)),
      axis.text = element_text(size = rel(1)),
      panel.grid.major = element_line(color = "white"),
      panel.grid.minor = element_line(color = "white"),
      legend.text = element_text(size = 10),
      legend.title = element_text(size = 15)
    )
  p <- p + mytheme + guides(colour = guide_legend(override.aes = list(size = 4)))
  if (plot) {
    print(p)
  }
  # print(p + mytheme)
  return(list("p" = p, "x" = X, "cell_group" = cell_group))
}
