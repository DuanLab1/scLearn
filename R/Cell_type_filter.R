#' @title Filter Cell Types by Minimum Cell Count
#'
#' @description
#' Filters out cell types (and optionally time points) that don't meet the minimum
#' required number of cells. This is particularly useful for ensuring sufficient
#' representation of each cell type for downstream analyses like differential
#' expression or classification.
#'
#' @details
#' The function performs the following operations:
#' \itemize{
#'   \item Combines cell type and time point information (if provided) to create
#'         composite group identifiers
#'   \item Calculates cell counts per group
#'   \item Removes groups with cell counts below the specified threshold
#'   \item Returns filtered expression matrix and corresponding metadata
#' }
#' When time points are provided, filtering is performed on the combined
#' cell_type-timepoint groups.
#'
#' @author Bin Duan (binduan\@sjtu.edu.cn)
#'
#' @param expression_profile Gene expression matrix (genes x cells)
#' @param sample_information_cellType Named vector of cell type annotations
#'        (names should match column names of expression_profile)
#' @param sample_information_timePoint Optional named vector of time point
#'        annotations (default: NULL)
#' @param min_cell_number Minimum number of cells required per group (default: 10)
#'
#' @return A list containing:
#' \itemize{
#'   \item expression_profile: Filtered expression matrix
#'   \item sample_information_cellType: Filtered cell type annotations
#'   \item sample_information_timePoint: Filtered time point annotations (if provided)
#' }
#'
#' @examples
#' \dontrun{
#' # Using example data from scLearn package
#' data(QueryCellData)
#'
#' # Get expression matrix (assuming counts are in 'counts' assay)
#' counts <- assay(QueryCellData, "counts")
#' cell_types <- colData(QueryCellData)$cell_type1
#' names(cell_types) <- colnames(QueryCellData)
#'
#' # Basic filtering
#' filtered_data <- Cell_type_filter(
#'   expression_profile = counts,
#'   sample_information_cellType = cell_types,
#'   min_cell_number = 20
#' )
#'
#' # With time points
#' time_points <- rep(c("Day1", "Day2"), length.out = ncol(QueryCellData))
#' names(time_points) <- colnames(QueryCellData)
#'
#' filtered_data <- Cell_type_filter(
#'   expression_profile = counts,
#'   sample_information_cellType = cell_types,
#'   sample_information_timePoint = time_points,
#'   min_cell_number = 5
#' )
#' }
#'
#' @export
#'
Cell_type_filter <- function(
    expression_profile,
    sample_information_cellType,
    sample_information_timePoint = NULL,
    min_cell_number = 10) {

  if (is.null(sample_information_timePoint)) {
    sample_information <- sample_information_cellType
  } else {
    sample_information <- paste(sample_information_cellType, sample_information_timePoint, sep = "_")
  }

  names(sample_information) <- names(sample_information_cellType)
  cell_number_type <- table(sample_information)
  cell_type_filtered <- names(cell_number_type[cell_number_type < min_cell_number])
  sample_information_filtered <- sample_information[!(sample_information %in% cell_type_filtered)]
  expression_profile_filtered <- expression_profile[, names(sample_information_filtered)]
  sample_information_cellType <- sample_information_cellType[names(sample_information_filtered)]

  if (!is.null(sample_information_timePoint)) {
    sample_information_timePoint <- sample_information_timePoint[names(sample_information_filtered)]
    return(list("expression_profile" = expression_profile_filtered, "sample_information_cellType" = sample_information_cellType, "sample_information_timePoint" = sample_information_timePoint))
  } else {
    return(list("expression_profile" = expression_profile_filtered, "sample_information_cellType" = sample_information_cellType))
  }
}
