#' @title Train scLearn Reference Model with Selectable Metric Learning
#'
#' @description
#' Trains a reference model for cell type annotation using either DCA or LMNN metric learning.
#' Supports both standard mode (cell type only) and time-course aware mode (cell type + time point).
#'
#' @param high_varGene_names Character vector of high-variance gene names
#' @param expression_profile Expression matrix (genes x cells)
#' @param sample_information_cellType Named vector of cell type annotations
#' @param sample_information_timePoint Optional time point annotations (default: NULL)
#' @param method Metric learning method ("DCA" or "LMNN", default: "DCA")
#' @param bootstrap_times Bootstrap iterations for standard mode (default: 10)
#' @param cutoff Percentile cutoff for similarity threshold (default: 0.01)
#' @param dim_para Variance explained threshold for time-course mode (default: 0.999)
#' @param ... Additional parameters passed to metric learning functions
#'
#' @return A scLearn model object containing:
#' \itemize{
#'   \item high_varGene_names: High-variance genes used
#'   \item simi_threshold_learned: List of correlation thresholds
#'   \item feature_matrix_learned: Reference feature matrices
#'   \item trans_matrix_learned: Transformation matrices
#'   \item method: Used metric learning method
#' }
#'
#' @importFrom Matrix nearPD
#' @importFrom stats cor
#' @export
scLearn_model_learning <- function(
    high_varGene_names,
    expression_profile,
    sample_information_cellType,
    sample_information_timePoint = NULL,
    method = c("DCA", "LMNN"),
    bootstrap_times = 10,
    cutoff = 0.01,
    dim_para = 0.999,
    ...) {

  # Validate inputs
  method <- match.arg(method)
  if (!all(high_varGene_names %in% rownames(expression_profile))) {
    stop("Not all high_varGene_names found in expression_profile")
  }

  # Internal helper functions -------------------------------------------------

  # Cluster mean features
  Feature_cluster <- function(expression_profile, sample_information) {
    sample_information <- sort(sample_information)
    expression_profile <- expression_profile[, names(sample_information)]
    num_each_class <- table(as.character(sample_information))
    feature_matrix <- matrix(0, length(num_each_class), nrow(expression_profile))
    rownames(feature_matrix) <- names(num_each_class)

    a <- 1
    for (i in seq_along(num_each_class)) {
      end_idx <- a + num_each_class[i] - 1
      feature_matrix[i, ] <- rowMeans(expression_profile[, a:end_idx, drop = FALSE])
      a <- end_idx + 1
    }
    return(feature_matrix)
  }

  # Calculate similarity thresholds
  Threshold_similarity <- function(expression_profile, sample_information, cutoff = 0.01) {
    sample_information <- sort(sample_information)
    expression_profile <- expression_profile[, names(sample_information)]
    num_each_class <- table(as.character(sample_information))
    simi_mean_list <- list()
    thre <- numeric(length(num_each_class))

    a <- 1
    for (i in seq_along(num_each_class)) {
      end_idx <- a + num_each_class[i] - 1
      current_expr <- expression_profile[, a:end_idx, drop = FALSE]
      feature_vec <- rowMeans(current_expr)
      correlations <- apply(current_expr, 2, cor, y = feature_vec)
      simi_mean_list[[i]] <- correlations
      thre[i] <- sort(correlations, decreasing = FALSE)[
        ceiling(length(correlations) * cutoff) + 1]
      a <- end_idx + 1
    }
    names(thre) <- names(num_each_class)
    names(simi_mean_list) <- names(num_each_class)
    return(list(simi_mean_list = simi_mean_list, threshold = thre))
  }

  # Time-course processing
  process_time_course <- function(expr, high_varGenes, cell_types, time_points, dim_para) {
    # Create label matrix
    combined_labels <- paste(cell_types, time_points, sep = "_")
    unique_labels <- c(levels(factor(cell_types)), levels(factor(time_points)))
    label_matrix <- matrix(0, length(unique_labels), ncol(expr))
    rownames(label_matrix) <- unique_labels
    colnames(label_matrix) <- colnames(expr)

    # Fill label matrix
    for (i in seq_along(unique_labels)) {
      label_matrix[i, combined_labels == unique_labels[i]] <- 1
    }

# Compute kernel matrix
L_kernel <- crossprod(label_matrix)

# Multi-domain discriminant analysis
mddm_analysis <- function(X, L, dim_para) {
  # Center kernel matrix
  centered_L <- scale(L, scale = FALSE)
  centered_L <- t(scale(t(centered_L), scale = FALSE))

  # Compute scatter matrix
  S <- X %*% centered_L %*% t(X)

  # Eigen decomposition
  eig <- eigen(S)
  eigenvalues <- Re(eig$values)
  eigenvectors <- Re(eig$vectors)

  # Determine projection dimension
  cum_var <- cumsum(eigenvalues) / sum(eigenvalues)
  proper_dim <- which(cum_var >= dim_para)[1]

  list(
    Projection_matrix = eigenvectors[, 1:proper_dim, drop = FALSE],
    eigenvalues = eigenvalues[1:proper_dim]
  )
}

# Perform analysis
mddm_result <- mddm_analysis(expr[high_varGenes, ], L_kernel, dim_para)
trans_matrix <- t(mddm_result$Projection_matrix)

list(
  trans_matrix = trans_matrix,
  transformed_data = trans_matrix %*% expr[high_varGenes, ],
  combined_labels = combined_labels
)
  }

  # Main training logic ------------------------------------------------------

  if (is.null(sample_information_timePoint)) {
    # Standard mode (with bootstrapping)
    threshold_cluster_trans <- vector("list", bootstrap_times)
    feature_matrix_trans <- vector("list", bootstrap_times)
    trans_matrix <- vector("list", bootstrap_times)

    for (r in seq_len(bootstrap_times)) {
      message("Metric learning iteration ", r, "/", bootstrap_times)

      # Select metric learning method
      if (method == "DCA") {
        trans_result <- runDCA(
          high_varGenes = high_varGene_names,
          expression_profile = expression_profile,
          sample_information = sample_information_cellType,
          seed = r,
          ...
        )
      } else {
        trans_result <- runLMNN(
          high_varGenes = high_varGene_names,
          expression_profile = expression_profile,
          sample_information = sample_information_cellType,
          seed = r,
          ...
        )
      }

      # Store results
      trans_matrix[[r]] <- trans_result$trans_matrix
      thre_result <- Threshold_similarity(
        trans_result$expression_profile_trans,
        trans_result$sample_information,
        cutoff = cutoff
      )
      threshold_cluster_trans[[r]] <- thre_result$threshold
      feature_matrix_trans[[r]] <- Feature_cluster(
        trans_result$expression_profile_trans,
        trans_result$sample_information
      )
    }
  } else {
    # Time-course aware mode
    message("Running in time-course aware mode...")
    time_result <- process_time_course(
      expression_profile,
      high_varGene_names,
      sample_information_cellType,
      sample_information_timePoint,
      dim_para
    )

    # Apply selected metric learning to transformed data
    if (method == "DCA") {
      trans_result <- runDCA(
        high_varGenes = high_varGene_names,
        expression_profile = expression_profile,
        sample_information = time_result$combined_labels,
        ...
      )
    } else {
      trans_result <- runLMNN(
        high_varGenes = high_varGene_names,
        expression_profile = expression_profile,
        sample_information = time_result$combined_labels,
        ...
      )
    }

    # Combine with time-course transformation
    final_trans_matrix <- trans_result$trans_matrix %*% time_result$trans_matrix
    final_transformed <- final_trans_matrix %*% expression_profile[high_varGene_names, ]

    # Calculate thresholds and features
    thre_result <- Threshold_similarity(
      final_transformed,
      time_result$combined_labels,
      cutoff = cutoff
    )
    threshold_cluster_trans <- list(thre_result$threshold)
    feature_matrix_trans <- list(
      Feature_cluster(final_transformed, time_result$combined_labels)
    )
    trans_matrix <- list(final_trans_matrix)
  }

  # Return standardized model object
  structure(
    list(
      high_varGene_names = high_varGene_names,
      simi_threshold_learned = threshold_cluster_trans,
      feature_matrix_learned = feature_matrix_trans,
      trans_matrix_learned = trans_matrix,
      method = method
    ),
    class = "scLearn_model"
  )
}


# scLearn_model_learning <- function(
#     high_varGene_names,
#     expression_profile,
#     sample_information_cellType,
#     sample_information_timePoint = NULL,
#     bootstrap_times = 10,
#     cutoff = 0.01,
#     dim_para = 0.999) {
#
#   qiu_ji <- function(x, y) {
#     return(x %*% y)
#   }
#
#   Feature_cluster <- function(
#       expression_profile,
#       sample_information) {
#     sample_information <- sort(sample_information)
#     expression_profile <- expression_profile[, names(sample_information)]
#     num_each_class <- table(as.character(sample_information))
#     feature_matrix <- matrix(0, length(num_each_class), nrow(expression_profile))
#     row.names(feature_matrix) <- names(num_each_class)
#
#     a <- 1
#
#     for (i in 1:(length(num_each_class))) {
#       # print(i)
#       expression_profile_choose <- expression_profile[, a:(a + num_each_class[i] - 1)]
#       a <- a + num_each_class[i]
#       feature_matrix[i, ] <- apply(expression_profile_choose, 1, mean)
#     }
#     return(feature_matrix)
#   }
#
#   Threshold_similarity <- function(
#       expression_profile,
#       sample_information,
#       cutoff = 0.01) {
#     sample_information <- sort(sample_information)
#     expression_profile <- expression_profile[, names(sample_information)]
#     num_each_class <- table(as.character(sample_information))
#     simi_mean_list <- list()
#
#     km <- function(vec, k) {
#       km <- mean(sort(vec, decreasing = TRUE)[2:(k + 1)])
#       return(km)
#     }
#
#     a <- 1
#     thre <- c()
#
#     for (i in 1:(length(num_each_class))) {
#       # print(i)
#       expression_profile_choose <- expression_profile[, a:(a + num_each_class[i] - 1)]
#       a <- a + num_each_class[i]
#       feature_vecter <- apply(expression_profile_choose, 1, mean)
#       simi_mean_list[[i]] <- apply(expression_profile_choose, 2, cor, y = feature_vecter)
#       thre[i] <- sort(simi_mean_list[[i]], decreasing = F)[ceiling(length(simi_mean_list[[i]]) * cutoff) + 1]
#     }
#     names(thre) <- names(num_each_class)
#     names(simi_mean_list) <- names(num_each_class)
#     return(list("simi_mean_list" = simi_mean_list, "threshold" = thre))
#   }
#
#   runDCA <- function(
#       high_varGenes,
#       expression_profile,
#       sample_information,
#       strength = 0.1,
#       seed = 1) {
#
#     sample_information <- sort(sample_information)
#     expression_profile <- expression_profile[high_varGenes, names(sample_information)]
#     num_eachclass <- table(sample_information)
#     # print("building chunks ...")
#
#     chks <- list()
#     a <- 1
#     k <- 1
#
#     for (i in 1:(length(num_eachclass))) {
#       s <- ceiling(num_eachclass[i] * strength) + 1
#       set.seed(seed)
#       chks[[k]] <- sample(a:(a + num_eachclass[i] - 1), s)
#       k <- k + 1
#       set.seed(seed)
#       chks[[k]] <- sample(setdiff(a:(a + num_eachclass[i] - 1), chks[[k - 1]]), s)
#       k <- k + 1
#       a <- a + num_eachclass[i]
#     }
#
#     chunks <- rep(-1, sum(num_eachclass))
#     for (i in 1:(2 * (length(num_eachclass)))) {
#       for (j in chks[[i]]) {
#         chunks[j] <- i
#       }
#     }
#
#     # print("building negtive links ...")
#     neglinks <- matrix(
#       rep(1, 4 * (length(num_eachclass)) * (length(num_eachclass))),
#       ncol = 2 * (length(num_eachclass)),
#       byrow = TRUE)
#
#     for (i in seq(1, ncol(neglinks), 2)) {
#       neglinks[i, i] <- 0
#       neglinks[i, i + 1] <- 0
#       neglinks[i + 1, i] <- 0
#       neglinks[i + 1, i + 1] <- 0
#     }
#
#     # print("performing dca ...")
#     dca_result <- dml::dca(data = t(expression_profile), chunks = chunks, neglinks = neglinks)
#     newData <- dca_result$newData
#     colnames(newData) <- paste("dc", 1:ncol(newData), sep = "")
#     trans_matrix <- dca_result$DCA
#     ex <- t(newData)
#     # print("Done!")
#     return(list(
#       "expression_profile_trans" = ex,
#       "expression_profile_origin" = expression_profile,
#       "trans_matrix" = trans_matrix,
#       "sample_information" = sample_information
#     ))
#   }
#
#   if (is.null(sample_information_timePoint)) {
#     threshold_cluster_trans <- list()
#     feature_matrix_trans <- list()
#     trans_matrix <- list()
#
#     for (r in 1:bootstrap_times) {
#       print(paste("Bootstrapying", r))
#       trans_result <- runDCA(
#         high_varGene_names,
#         expression_profile,
#         sample_information_cellType,
#         strength = 0.1,
#         seed = r
#       )
#
#       trans_matrix[[r]] <- trans_result$trans_matrix
#       thre_result_trans_cluster <- Threshold_similarity(
#         trans_result$expression_profile_trans,
#         trans_result$sample_information,
#         cutoff = cutoff
#       )
#       threshold_cluster_trans[[r]] <- thre_result_trans_cluster$threshold
#       feature_matrix_trans[[r]] <- Feature_cluster(trans_result$expression_profile_trans, trans_result$sample_information)
#     }
#
#     return(list(
#       "high_varGene_names" = high_varGene_names,
#       "simi_threshold_learned" = threshold_cluster_trans,
#       "feature_matrix_learned" = feature_matrix_trans,
#       "trans_matrix_learned" = trans_matrix
#     ))
#   } else {
#     label_len <- length(table(sample_information_cellType)) + length(table(sample_information_timePoint))
#     label_matrix <- matrix(0, label_len, ncol(expression_profile))
#
#     row.names(label_matrix) <- c(names(table(sample_information_cellType)), names(table(sample_information_timePoint)))
#
#     colnames(label_matrix) <- colnames(expression_profile)
#
#     sample_information_all <- c(sample_information_cellType, sample_information_timePoint)
#     label_all <- c(names(table(sample_information_cellType)), names(table(sample_information_timePoint)))
#
#     for (i in 1:nrow(label_matrix)) {
#       label_matrix[i, names(sample_information_all[sample_information_all == label_all[i]])] <- 1
#     }
#
#     L_kernel <- matrix(0, ncol(label_matrix), ncol(label_matrix))
#     colnames(L_kernel) <- colnames(expression_profile)
#     row.names(L_kernel) <- colnames(expression_profile)
#     L_kernel <- t(label_matrix) %*% label_matrix
#
#     getProperDim <- function(lambda, dim_para = 0.999) {
#       thr <- dim_para
#       sum_lambda <- sum(lambda)
#       lambda_num <- length(lambda)
#       tmp_lambda <- 0
#       for (lind in 1:lambda_num) {
#         tmp_lambda <- tmp_lambda + lambda[lind]
#         if (tmp_lambda >= thr * sum(lambda)) {
#           proper_dim <- lind
#           return(proper_dim)
#         }
#       }
#     }
#
#     mddm_linear <- function(X, L, dim_para = 0.999) {
#
#       D <- nrow(X)
#       N <- ncol(X)
#
#       nor_mean <- function(vec) {
#         vec <- vec - mean(vec)
#         return(vec)
#       }
#       tmpL <- apply(L, 2, nor_mean)
#       HLH <- apply(tmpL, 1, nor_mean)
#       S <- X %*% HLH %*% t(X)
#       B <- diag(D)
#       B_tr <- B
#       eig_result <- eigen(B_tr %*% S)
#       tmp_lambda <- diag(eig_result$values)
#       tmp_P <- eig_result$vectors
#       tmp_P <- Re(tmp_P)
#       tmp_lambda <- Re(tmp_lambda[row(tmp_lambda) == col(tmp_lambda)])
#       lambda <- sort(tmp_lambda, decreasing = T)
#       order <- order(tmp_lambda, decreasing = T)
#       P <- tmp_P[, order]
#       proper_dim <- getProperDim(lambda, dim_para)
#       P <- P[, 1:proper_dim]
#       lambda <- lambda[1:proper_dim]
#       return(list("Projection_matrix" = P, "eigenvalues" = lambda))
#     }
#
#     mddm_result <- mddm_linear(expression_profile[high_varGene_names, ], L_kernel, dim_para = dim_para)
#     trans_result <- list()
#     trans_result$trans_matrix <- t(mddm_result$Projection_matrix)
#     trans_result$expression_profile_orign <- expression_profile[high_varGene_names, ]
#     trans_result$expression_profile_trans <- t(mddm_result$Projection_matrix) %*% expression_profile[high_varGene_names, ]
#     sample_information_cb <- paste(sample_information_cellType, sample_information_timePoint, sep = "_")
#     names(sample_information_cb) <- names(sample_information_cellType)
#     trans_result$sample_information_combine <- sample_information_cb
#
#     thre_result_trans_cluster <- Threshold_similarity(
#       trans_result$expression_profile_trans,
#       trans_result$sample_information_combine,
#       cutoff = cutoff
#     )
#
#     threshold_cluster_trans <- thre_result_trans_cluster$threshold
#     feature_matrix_trans <- Feature_cluster(
#       trans_result$expression_profile_trans,
#       trans_result$sample_information_combine
#     )
#
#     return(list(
#       "high_varGene_names" = high_varGene_names,
#       "simi_threshold_learned" = threshold_cluster_trans,
#       "feature_matrix_learned" = feature_matrix_trans,
#       "trans_matrix_learned" = trans_result$trans_matrix
#     ))
#   }
# }


#' @title Discriminative Component Analysis (DCA) Transformation
#'
#' @description
#' Performs DCA metric learning for single-cell data. This is the original metric learning
#' method used in scLearn, implementing the algorithm from the dml package.
#'
#' @param high_varGenes Character vector of high-variance genes
#' @param expression_profile Numeric matrix (genes x cells)
#' @param sample_information Named vector of cell type labels
#' @param strength Subsampling strength (default: 0.1)
#' @param seed Random seed (default: 1)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return List containing:
#' \itemize{
#'   \item expression_profile_trans: Transformed matrix (genes x cells)
#'   \item expression_profile_origin: Original matrix
#'   \item trans_matrix: Learned transformation matrix
#'   \item sample_information: Input cell labels
#' }
#'
#' @importFrom dml dca
#'
#' @export
#'
runDCA <- function(
    high_varGenes,
    expression_profile,
    sample_information,
    strength = 0.1,
    seed = 1,
    verbose = TRUE) {

  # Input validation
  if (!requireNamespace("dml", quietly = TRUE)) {
    stop("dml package required for DCA. Install with: remotes::install_github('dml')")
  }

  if (!all(names(sample_information) %in% colnames(expression_profile))) {
    stop("Cell names in sample_information don't match expression_profile columns")
  }

  # Sort by cell type for consistency
  sample_information <- sort(sample_information)
  expression_profile <- expression_profile[high_varGenes, names(sample_information)]

  num_eachclass <- table(sample_information)

  # Build chunks for DCA
  set.seed(seed)
  chks <- list()
  a <- 1
  k <- 1

  for (i in seq_along(num_eachclass)) {
    s <- ceiling(num_eachclass[i] * strength) + 1
    chks[[k]] <- sample(a:(a + num_eachclass[i] - 1), s)
    k <- k + 1
    chks[[k]] <- sample(setdiff(a:(a + num_eachclass[i] - 1), chks[[k-1]]), s)
    k <- k + 1
    a <- a + num_eachclass[i]
  }

  chunks <- rep(-1, sum(num_eachclass))
  for (i in 1:(2 * (length(num_eachclass)))) {
    for (j in chks[[i]]) {
      chunks[j] <- i
    }
  }

  # Build negative links matrix
  neglinks <- matrix(1, ncol = 2*length(num_eachclass), nrow = 2*length(num_eachclass))
  for (i in seq(1, ncol(neglinks), 2)) {
    neglinks[i, i] <- 0
    neglinks[i, i+1] <- 0
    neglinks[i+1, i] <- 0
    neglinks[i+1, i+1] <- 0
  }

  # Run DCA
  if (verbose) message("Running DCA...")
  dca_result <- dml::dca(
    data = t(expression_profile),
    chunks = chunks,
    neglinks = neglinks
  )

  # Prepare output
  list(
    expression_profile_trans = t(dca_result$newData),
    expression_profile_origin = expression_profile,
    trans_matrix = dca_result$DCA,
    sample_information = sample_information
  )
}


#' @title LMNN Transformation for scLearn
#'
#' @description
#' Implements Large Margin Nearest Neighbor (LMNN) metric learning for single-cell data.
#' Returns transformed data in scLearn-compatible format.
#'
#' @param high_varGenes Character vector of high-variance genes
#' @param expression_profile Numeric matrix (genes x cells)
#' @param sample_information Named vector of cell type labels
#' @param k Number of target neighbors (default: 3)
#' @param learn_rate Learning rate (default: 1e-4)
#' @param max_iter Maximum iterations (default: 10)
#' @param seed Random seed (default: 1)
#' @param verbose Print progress (default: TRUE)
#'
#' @return List containing:
#' \itemize{
#'   \item expression_profile_trans: Transformed matrix (genes x cells)
#'   \item expression_profile_origin: Original matrix
#'   \item trans_matrix: Learned transformation matrix
#'   \item sample_information: Input cell labels
#' }
#'
#' @importFrom Matrix nearPD
#' @importFrom Rcpp sourceCpp
#' @importFrom stats setNames
#'
#' @export
#'
runLMNN <- function(
    high_varGenes,
    expression_profile,
    sample_information,
    k = 20,
    learn_rate = 1e-4,
    max_iter = 10,
    seed = 1,
    verbose = TRUE) {

  # Input validation
  if (!requireNamespace("Rcpp", quietly = TRUE)) {
    stop("Rcpp package required for LMNN")
  }

  if (!all(high_varGenes %in% rownames(expression_profile))) {
    stop("high_varGenes not found in expression_profile")
  }

  # Sort by cell type for consistency
  sample_information <- sort(sample_information)
  expr <- expression_profile[high_varGenes, names(sample_information), drop = FALSE]

  # Transpose for LMNN (cells x genes)
  X <- t(as.matrix(expr))
  y <- as.character(sample_information)
  n_cells <- nrow(X)
  n_genes <- ncol(X)

  # Initialize metric matrix
  M <- diag(n_genes)

  # Find target neighbors (k-NN within classes)
  find_target_neighbors <- function(X, y, k) {
    targets <- matrix(0, nrow(X), k)
    for (i in seq_len(nrow(X))) {
      same_class <- which(y == y[i])
      same_class <- setdiff(same_class, i) # Exclude self
      if (length(same_class) > 0) {
        dists <- rowSums(sweep(X[same_class, , drop = FALSE], 2, X[i, ], "-")^2)
        targets[i, seq_len(min(k, length(same_class)))] <-
          same_class[order(dists)[seq_len(min(k, length(same_class)))]]
      }
    }
    return(targets)
  }

  set.seed(seed)
  targets <- find_target_neighbors(X, y, k)

  # Main optimization loop
  for (iter in seq_len(max_iter)) {
    if (verbose && iter %% 10 == 0) {
      message(sprintf("Iteration %d/%d", iter, max_iter))
    }

    # Compute pairwise distances
    # D <- mahalanobis_distance_cpp(X, M)
    D <- mahalanobis_distance_r(X, M)

    # Compute gradient
    grad <- matrix(0, n_genes, n_genes)

    # Pull gradient (attract neighbors)
    for (i in seq_len(n_cells)) {
      valid_targets <- targets[i, ][targets[i, ] != 0]
      for (j in valid_targets) {
        diff <- X[i, ] - X[j, ]
        grad <- grad + tcrossprod(diff)
      }
    }

    # Push gradient (repel impostors)
    for (i in seq_len(n_cells)) {
      valid_targets <- targets[i, ][targets[i, ] != 0]
      for (j in valid_targets) {
        impostors <- which((y != y[i]) & (D[i, ] - D[i, j] < 1))
        for (l in impostors) {
          diff_il <- X[i, ] - X[l, ]
          diff_ij <- X[i, ] - X[j, ]
          grad <- grad + (tcrossprod(diff_il) - tcrossprod(diff_ij))
        }
      }
    }

    # Update metric matrix
    M <- M - learn_rate * grad

    # Ensure symmetry and positive semi-definiteness
    M <- (M + t(M)) / 2
    eig <- eigen(M, symmetric = TRUE)
    M <- eig$vectors %*% diag(pmax(eig$values, 0)) %*% t(eig$vectors)
  }

  # Compute transformation matrix
  M_pd <- Matrix::nearPD(M)$mat
  L <- chol(M_pd) # M = L'L

  # Transform data (transpose back to genes x cells)
  X_trans <- t(X %*% t(L))

  # Return in scLearn-compatible format
  list(
    expression_profile_trans = X_trans,
    expression_profile_origin = expr,
    trans_matrix = t(L), # Transpose to match DCA convention
    sample_information = sample_information
  )
}

#' Compute Mahalanobis Distance Matrix
#'
#' Calculates the pairwise Mahalanobis distance between all rows of a matrix `X`
#' using the given precision matrix `M` (inverse covariance).
#'
#' @param X A numeric matrix where each row represents an observation.
#' @param M A symmetric positive-definite matrix (precision matrix) used to compute the distance.
#'
#' @return A symmetric matrix `D` where `D[i, j]` is the Mahalanobis distance
#'         between rows `i` and `j` of `X`.
#'
#' @keywords internal
#'
mahalanobis_distance_r <- function(X, M) {
  n <- nrow(X)
  D <- matrix(0, n, n)
  for (i in 1:n) {
    diffs <- sweep(X, 2, X[i, ])
    D[i, ] <- rowSums((diffs %*% M) * diffs)
  }
  return(D)
}
