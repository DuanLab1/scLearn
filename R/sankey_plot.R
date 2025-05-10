#' @title Generate Sankey Diagram for Cell Type Projection Results
#'
#' @description
#' Creates a Sankey diagram to visualize the correspondence between original and projected
#' cell type labels, typically used to show scRNA-seq cell type annotation results.
#'
#' @details
#' This function:
#' \itemize{
#'   \item Processes prediction results to include all reference cell types
#'   \item Generates a Sankey diagram showing label transitions
#'   \item Handles cases where some reference cell types are not predicted
#' }
#' The visualization helps evaluate how cell types are mapped between reference and query datasets.
#'
#' @author Bin Duan (binduan\@sjtu.edu.cn)
#'
#' @param predict_result A data frame with 3 columns:
#'   \itemize{
#'     \item Column 1: Original labels
#'     \item Column 2: Projected labels
#'     \item Column 3: Sample names
#'   }
#' @param sample_information_reference Named vector of reference cell type annotations
#' @param plot Whether to display the Sankey diagram (default: FALSE)
#'
#' @return The augmented prediction data frame including all reference cell types.
#' If plot=TRUE, also displays the Sankey diagram.
#'
#' @importFrom dplyr %>%
#'
#' @examples
#' \dontrun{
#' # Mock prediction results
#' predict_df <- data.frame(
#'   orign_label = c("Tcell", "Bcell", "Tcell", "Macrophage"),
#'   projection_label = c("CD4+", "CD19+", "CD8+", "Mono"),
#'   sample_name = rep("sample1", 4)
#' )
#'
#' # Reference cell types
#' ref_labels <- c("Tcell", "Bcell", "Macrophage", "Neutrophil")
#'
#' # Generate plot
#' sankey_plot(
#'   predict_result = predict_df,
#'   sample_information_reference = ref_labels,
#'   plot = TRUE
#' )
#'
#' # Get augmented results without plotting
#' augmented_results <- sankey_plot(
#'   predict_result = predict_df,
#'   sample_information_reference = ref_labels
#' )
#' }
#'
#' @seealso
#' \code{\link[scmap]{getSankey}} for the underlying Sankey diagram implementation
#'
#' @export
#'
sankey_plot <- function(
    predict_result,
    sample_information_reference,
    plot = FALSE) {

  cell_type_reference <- names(table(sample_information_reference))
  other_cell_type <- setdiff(cell_type_reference, intersect(cell_type_reference, unique(predict_result[, 2])))

  if (length(other_cell_type) > 0) {
    other_cell_type <- as.data.frame(other_cell_type)
    other_cell_type$orign_label <- "NULL"
    other_cell_type$sample_name <- "for_sankey_plot"
    colnames(other_cell_type) <- c("projection_label", "orign_label", "sample_name")
    predict_result <- rbind(predict_result, other_cell_type)
  }

  if (plot) {
    plot(scmap::getSankey(predict_result[, 1], predict_result[, 2], plot_width = 300, plot_height = 300))
  }

  return(predict_result)
}
