#' @title Cell Type Prediction Using scLearn Model
#'
#' @description
#' Predicts cell types for query single-cell RNA-seq data using a pre-trained scLearn model.
#' Implements a voting-based assignment system with multiple quality control checks.
#'
#' @details
#' The prediction process involves:
#' \itemize{
#'   \item Feature matching between query data and model features
#'   \item Multiple transformation matrices application (if available)
#'   \item Correlation-based cell type assignment
#'   \item Voting mechanism for consensus prediction
#'   \item Novel cell type detection and quality checks
#' }
#' Three possible outcomes for each cell:
#' \itemize{
#'   \item Specific cell type assignment
#'   \item "unassigned" (low confidence)
#'   \item Quality flags (Gene_Missing/Novel_Cell/Too similar)
#' }
#'
#' @author Bin Duan (binduan\@sjtu.edu.cn)
#'
#' @param scLearn_model_learning_result A trained scLearn model object containing:
#'   \itemize{
#'     \item high_varGene_names: Vector of high-variance genes
#'     \item trans_matrix_learned: Transformation matrix/matrices
#'     \item feature_matrix_learned: Reference feature matrix/matrices
#'     \item simi_threshold_learned: Correlation threshold(s)
#'   }
#' @param expression_profile_query Query expression matrix (genes x cells)
#' @param vote_rate Minimum vote proportion for consensus (default: 0.6)
#' @param diff Minimum correlation difference between top candidates (default: 0.05)
#' @param threshold_use Whether to use correlation threshold (default: FALSE)
#'
#' @return A data frame with prediction results containing:
#'   \itemize{
#'     \item Query_cell_id: Cell identifiers
#'     \item Predict_cell_type: Predicted cell type or "unassigned"
#'     \item Additional_information: Quality flags (when single matrix used)
#'   }
#'
#' @importFrom stats cor
#'
#' @examples
#' \dontrun{
#' # Load example scLearn model and query data
#' data(scLearn_model)
#' data(QueryCellData)
#'
#' # Get query expression matrix
#' query_data <- logcounts(QueryCellData)
#'
#' # Basic prediction
#' predictions <- scLearn_cell_assignment(
#'   scLearn_model_learning_result = scLearn_model,
#'   expression_profile_query = query_data
#' )
#'
#' # Strict prediction with threshold
#' strict_pred <- scLearn_cell_assignment(
#'   scLearn_model_learning_result = scLearn_model,
#'   expression_profile_query = query_data,
#'   vote_rate = 0.7,
#'   diff = 0.1,
#'   threshold_use = TRUE
#' )
#'
#' # Examine results
#' table(predictions$Predict_cell_type)
#' head(predictions)
#' }
#'
#' @seealso
#' \code{\link{scLearn_model_learning}} for model training function
#'
#' @export
#'
scLearn_cell_assignment <- function(
    scLearn_model_learning_result,
    expression_profile_query,
    vote_rate = 0.6,
    diff = 0.05,
    threshold_use = FALSE) {

  Vote_class <- function(vec, vote_rate = 0.6) {
    vec_len <- length(vec)
    num <- length(table(vec)[table(vec) == max(table(vec))])
    if (num > 1) {
      return("unassigned")
    } else if (max(table(vec)) > vec_len * vote_rate) {
      return(names(table(vec)[table(vec) == max(table(vec))]))
    } else {
      return("unassigned")
    }
  }

  Get_query_hvg <- function(expression_profile, high_varGene_names) {
    missing_num <- length(high_varGene_names) - length(intersect(
      row.names(expression_profile),
      high_varGene_names
    ))
    missing_features <- setdiff(high_varGene_names, intersect(
      row.names(expression_profile),
      high_varGene_names
    ))
    missing_rate <- missing_num / length(high_varGene_names)
    print(paste("The number of missing features in the query data is ",
      missing_num,
      seq = ""
    ))
    print(paste("The rate of missing features in the query data is ",
      missing_rate,
      seq = ""
    ))
    if (missing_num > 0) {
      missing_data <- matrix(0, missing_num, ncol(expression_profile))
      row.names(missing_data) <- missing_features
      expression_profile <- rbind(expression_profile, missing_data)
    }
    expression_profile_hvg <- expression_profile[high_varGene_names, ]

    return(expression_profile_hvg)
  }

  Assignment_result <- function(
      expression_profile_query_hvg,
      feature_matrix,
      threshold,
      diff = 0.05,
      threshold_use = TRUE) {
    feature_matrix <- t(feature_matrix)
    result <- data.frame(
      cluster_lab = 1:ncol(expression_profile_query_hvg),
      cluster_cor = rep(0, ncol(expression_profile_query_hvg))
    )
    row.names(result) <- colnames(expression_profile_query_hvg)
    if (threshold_use == FALSE) {
      threshold <- -2
    }
    for (i in 1:ncol(expression_profile_query_hvg)) {
      cor_result <- rep(0, ncol(feature_matrix))
      names(cor_result) <- colnames(feature_matrix)
      if (sum(expression_profile_query_hvg[, i]) == 0) {
        result[i, 1] <- "unassigned"
        result[i, 2] <- "Gene_Missing"
        next
      }
      for (j in 1:ncol(feature_matrix)) {
        cor_result[j] <- cor(expression_profile_query_hvg[
          ,
          i
        ], feature_matrix[, j])
      }
      cor_compare <- cor_result - threshold
      if (length(cor_compare[cor_compare > 0]) == 0) {
        result[i, 1] <- "unassigned"
        result[i, 2] <- "Novel_Cell"
      }
      if (length(cor_compare[cor_compare > 0]) == 1) {
        result[i, 1] <- names(cor_compare[cor_compare >
          0])
        result[i, 2] <- cor_result[result[i, 1]]
      }
      if (length(cor_compare[cor_compare > 0]) >= 2) {
        if (sort(cor_result[names(cor_compare[cor_compare >
          0])], decreasing = T)[1] - sort(cor_result[names(cor_compare[cor_compare >
          0])], decreasing = T)[2] >= diff) {
          result[i, 1] <- names(cor_result[names(cor_compare[cor_compare >
            0])])[which(cor_result[names(cor_compare[cor_compare >
            0])] == max(cor_result[names(cor_compare[cor_compare >
            0])]))][1]
          result[i, 2] <- cor_result[result[i, 1]]
        } else {
          result[i, 1] <- "unassigned"
          result[i, 2] <- "Too similar to multiple cell type"
        }
      }
    }
    return(result)
  }

  expression_profile_query_hvg <- Get_query_hvg(
    expression_profile_query,
    scLearn_model_learning_result$high_varGene_names
  )
  if (is.list(scLearn_model_learning_result$trans_matrix_learned)) {
    if (length(scLearn_model_learning_result$trans_matrix_learned) > 1) {
      predict_result <- matrix(
        0, ncol(expression_profile_query),
        length(scLearn_model_learning_result$trans_matrix_learned)
      )
      for (r in 1:length(scLearn_model_learning_result$trans_matrix_learned)) {
        expression_profile_query_hvg_ml <- scLearn_model_learning_result$trans_matrix_learned[[r]] %*%
          expression_profile_query_hvg
        assignment_result <- Assignment_result(expression_profile_query_hvg_ml,
          scLearn_model_learning_result$feature_matrix_learned[[r]],
          threshold = scLearn_model_learning_result$simi_threshold_learned[[r]],
          diff = diff, threshold_use = threshold_use
        )
        predict_result[, r] <- assignment_result[, 1]
      }
      predict_result_final <- apply(predict_result, 1, Vote_class,
        vote_rate = vote_rate
      )
      predict_result_final <- as.data.frame(predict_result_final)
      predict_result_final$sample <- colnames(expression_profile_query)
      predict_result_final$predict_result_final <- as.character(predict_result_final$predict_result_final)
      colnames(predict_result_final) <- c(
        "Predict_cell_type",
        "Query_cell_id"
      )
      predict_result_final <- predict_result_final[, c(
        "Query_cell_id",
        "Predict_cell_type"
      )]
      return(predict_result_final)
    } else {
      expression_profile_query_hvg_ml <- scLearn_model_learning_result$trans_matrix_learned[[1]] %*%
        expression_profile_query_hvg
      assignment_result <- Assignment_result(expression_profile_query_hvg_ml,
        scLearn_model_learning_result$feature_matrix_learned[[1]],
        threshold = scLearn_model_learning_result$simi_threshold_learned[[1]],
        diff = diff, threshold_use = threshold_use
      )
      assignment_result$Query_cell_id <- row.names(assignment_result)
      assignment_result$Predict_cell_type <- as.character(assignment_result$cluster_lab)
      assignment_result$Additional_information <- as.character(assignment_result$cluster_cor)
      assignment_result <- assignment_result[, c(
        "Query_cell_id",
        "Predict_cell_type", "Additional_information"
      )]
      return(assignment_result)
    }
  } else {
    expression_profile_query_hvg_ml <- scLearn_model_learning_result$trans_matrix_learned %*%
      expression_profile_query_hvg
    assignment_result <- Assignment_result(expression_profile_query_hvg_ml,
      scLearn_model_learning_result$feature_matrix_learned,
      threshold = scLearn_model_learning_result$simi_threshold_learned,
      diff = diff, threshold_use = threshold_use
    )
    assignment_result$Query_cell_id <- row.names(assignment_result)
    assignment_result$Predict_cell_type <- as.character(assignment_result$cluster_lab)
    assignment_result$Additional_information <- as.character(assignment_result$cluster_cor)
    assignment_result <- assignment_result[, c(
      "Query_cell_id",
      "Predict_cell_type", "Additional_information"
    )]
    return(assignment_result)
  }
}
