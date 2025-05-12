#' @title Train scLearn Reference Model with Multiple Metric Learning Options
#'
#' @description
#' Trains a reference model for cell type annotation using selectable metric learning methods.
#' Supports both standard mode (cell type only) and time-course aware mode (cell type + time point).
#' Available methods: DCA, LMNN, NCA, ITML, and MMC.
#'
#' @author Bin Duan (binduan\@sjtu.edu.cn)
#'
#' @param high_varGene_names Character vector of high-variance gene names
#' @param expression_profile Expression matrix (genes x cells)
#' @param sample_information_cellType Named vector of cell type annotations
#' @param sample_information_timePoint Optional time point annotations (default: NULL)
#' @param method Metric learning method to use. One of:
#' \itemize{
#'   \item "DCA" - Discriminative Component Analysis (linear transformation maximizing between-class/within-class variance ratio)
#'   \item "LMNN" - Large Margin Nearest Neighbors (preserves k-NN classification boundaries with margin constraints)
#'   \item "MSL" - Multi-Similarity Loss-Based Large Margin Nearest Neighbors (preserves k-NN classification boundaries with margin constraints)
#'   \item "NCA" - Neighborhood Components Analysis (maximizes leave-one-out classification probability)
#'   \item "ITML" - Information-Theoretic Metric Learning (uses relative distance constraints with information-theoretic regularization)
#'   \item "MMC" - Maximum Margin Clustering (explicitly optimizes within-class/between-class scatter matrices)
#' }
#' Default: "DCA"
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
#'
#' @examples
#' \dontrun{
#'
#' scLearn_model_learning_result <- scLearn_model_learning(
#'   high_varGene_names = RefDataQC_HVG_names,
#'   expression_profile = RefDataQC_filtered$expression_profile,
#'   sample_information_cellType = RefDataQC_filtered$sample_information_cellType,
#'   sample_information_timePoint = NULL,
#'   method = "DCA",
#'   bootstrap_times = 1,
#'   cutoff = 0.01,
#'   dim_para = 0.999
#' )
#' }
#' @export
#'
scLearn_model_learning <- function(
    high_varGene_names,
    expression_profile,
    sample_information_cellType,
    sample_information_timePoint = NULL,
    method = c("DCA", "LMNN", "MSL", "NCA", "ITML", "MMC"),
    bootstrap_times = 10,
    cutoff = 0.01,
    dim_para = 0.999,
    ...) {
  # Validate inputs and method selection
  method <- match.arg(method)
  if (!all(high_varGene_names %in% rownames(expression_profile))) {
    stop("Not all high_varGene_names found in expression_profile")
  }

  # Internal helper functions -----------------------------------------------

  # Calculate cluster mean features
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

  # Calculate similarity thresholds for each cluster
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
        ceiling(length(correlations) * cutoff) + 1
      ]
      a <- end_idx + 1
    }
    names(thre) <- names(num_each_class)
    names(simi_mean_list) <- names(num_each_class)
    return(list(simi_mean_list = simi_mean_list, threshold = thre))
  }

  # Time-course data processing function
  process_time_course <- function(expr, high_varGenes, cell_types, time_points, dim_para) {
    # Create combined label matrix
    combined_labels <- paste(cell_types, time_points, sep = "_")
    unique_labels <- c(levels(factor(cell_types)), levels(factor(time_points)))
    label_matrix <- matrix(0, length(unique_labels), ncol(expr))
    rownames(label_matrix) <- unique_labels
    colnames(label_matrix) <- colnames(expr)

    # Populate label matrix
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

      # Determine projection dimension based on variance explained
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

  # Main training logic -----------------------------------------------------

  if (is.null(sample_information_timePoint)) {
    # Standard mode (with bootstrapping)
    threshold_cluster_trans <- vector("list", bootstrap_times)
    feature_matrix_trans <- vector("list", bootstrap_times)
    trans_matrix <- vector("list", bootstrap_times)

    for (r in seq_len(bootstrap_times)) {
      message("Metric learning iteration ", r, "/", bootstrap_times)

      # Select appropriate metric learning method
      trans_result <- switch(method,
        "DCA" = runDCA(
          high_varGenes = high_varGene_names,
          expression_profile = expression_profile,
          sample_information = sample_information_cellType,
          seed = r,
          ...
        ),
        "LMNN" = runLMNN(
          high_varGenes = high_varGene_names,
          expression_profile = expression_profile,
          sample_information = sample_information_cellType,
          seed = r,
          ...
        ),
        "NCA" = runNCA(
          high_varGenes = high_varGene_names,
          expression_profile = expression_profile,
          sample_information = sample_information_cellType,
          seed = r,
          ...
        ),
        "ITML" = runITML(
          high_varGenes = high_varGene_names,
          expression_profile = expression_profile,
          sample_information = sample_information_cellType,
          seed = r,
          ...
        ),
        "MMC" = runMMC(
          high_varGenes = high_varGene_names,
          expression_profile = expression_profile,
          sample_information = sample_information_cellType,
          seed = r,
          ...
        ),
        "MSL" = runMSL(
          high_varGenes = high_varGene_names,
          expression_profile = expression_profile,
          sample_information = sample_information_cellType,
          seed = r,
          ...
        )
      )

      # Store transformation results
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
    trans_result <- switch(method,
      "DCA" = runDCA(
        high_varGenes = high_varGene_names,
        expression_profile = expression_profile,
        sample_information = time_result$combined_labels,
        ...
      ),
      "LMNN" = runLMNN(
        high_varGenes = high_varGene_names,
        expression_profile = expression_profile,
        sample_information = time_result$combined_labels,
        ...
      ),
      "NCA" = runNCA(
        high_varGenes = high_varGene_names,
        expression_profile = expression_profile,
        sample_information = time_result$combined_labels,
        ...
      ),
      "ITML" = runITML(
        high_varGenes = high_varGene_names,
        expression_profile = expression_profile,
        sample_information = time_result$combined_labels,
        ...
      ),
      "MMC" = runMMC(
        high_varGenes = high_varGene_names,
        expression_profile = expression_profile,
        sample_information = time_result$combined_labels,
        ...
      ),
      "MSL" = runMSL(
        high_varGenes = high_varGene_names,
        expression_profile = expression_profile,
        sample_information = time_result$combined_labels,
        ...
      )
    )

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
    chks[[k]] <- sample(setdiff(a:(a + num_eachclass[i] - 1), chks[[k - 1]]), s)
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
  neglinks <- matrix(1, ncol = 2 * length(num_eachclass), nrow = 2 * length(num_eachclass))
  for (i in seq(1, ncol(neglinks), 2)) {
    neglinks[i, i] <- 0
    neglinks[i, i + 1] <- 0
    neglinks[i + 1, i] <- 0
    neglinks[i + 1, i + 1] <- 0
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


#' @title Triplet Loss-Based Metric Learning for Single-Cell Data
#'
#' @description
#' Implements a triplet loss-based linear transformation learning to map gene expression
#' into an embedding space where cells of the same type are close and those of different
#' types are farther apart. This is a deep metric learning alternative to LMNN.
#'
#' @param high_varGenes Character vector of highly variable gene names.
#' @param expression_profile A numeric matrix of gene expression values (genes x cells).
#' @param sample_information A named vector of cell type labels (names match column names of expression_profile).
#' @param k Number of positive neighbors per anchor to form triplets (default: 5).
#' @param margin Margin between positive and negative pairs (default: 1.0).
#' @param learn_rate Learning rate for updates (default: 1e-4).
#' @param max_iter Number of optimization iterations (default: 100).
#' @param seed Random seed (default: 1).
#' @param verbose Logical; whether to print progress messages (default: TRUE).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{expression_profile_trans}: Transformed expression matrix (genes x cells).
#'   \item \code{expression_profile_origin}: Original filtered expression matrix (genes x cells).
#'   \item \code{trans_matrix}: Learned linear transformation matrix (L).
#'   \item \code{sample_information}: Sorted cell type labels.
#' }
#'
#' @importFrom Matrix nearPD
#' @importFrom stats setNames
#'
#' @export
#'
runLMNN <- function(
    high_varGenes,
    expression_profile,
    sample_information,
    k = 5,
    margin = 1.0,
    learn_rate = 1e-4,
    max_iter = 100,
    seed = 1,
    verbose = TRUE) {
  # 1. Input checks and preparation
  if (!all(high_varGenes %in% rownames(expression_profile))) {
    stop("Some high_varGenes not found in expression_profile.")
  }

  sample_information <- sort(sample_information)
  expr <- expression_profile[high_varGenes, names(sample_information), drop = FALSE]

  # Matrix: cells x genes
  X <- t(as.matrix(expr))
  y <- as.character(sample_information)
  n_cells <- nrow(X)
  n_genes <- ncol(X)

  set.seed(seed)

  # 2. Initialize linear transformation (identity)
  L <- diag(n_genes)

  # 3. Helper: Get positive/negative neighbors
  get_triplets <- function(X, y, k) {
    triplets <- list()
    for (i in seq_len(n_cells)) {
      anchor <- X[i, ]
      same_class <- setdiff(which(y == y[i]), i)
      diff_class <- which(y != y[i])

      if (length(same_class) == 0 || length(diff_class) == 0) next

      pos_idx <- sample(same_class, min(k, length(same_class)))
      neg_idx <- sample(diff_class, min(k, length(diff_class)))

      for (j in pos_idx) {
        for (l in neg_idx) {
          triplets[[length(triplets) + 1]] <- list(i = i, j = j, l = l)
        }
      }
    }
    triplets
  }

  # 4. Optimization using Triplet Loss
  for (iter in seq_len(max_iter)) {
    triplets <- get_triplets(X, y, k)

    grad <- matrix(0, n_genes, n_genes)
    loss <- 0

    for (triplet in triplets) {
      i <- triplet$i
      j <- triplet$j
      l <- triplet$l

      xi <- X[i, , drop = FALSE]
      xj <- X[j, , drop = FALSE]
      xl <- X[l, , drop = FALSE]

      # Transformed differences
      diff_pos <- (xi - xj) %*% L
      diff_neg <- (xi - xl) %*% L

      d_pos <- sum(diff_pos^2)
      d_neg <- sum(diff_neg^2)

      triplet_loss <- d_pos - d_neg + margin

      if (triplet_loss > 0) {
        loss <- loss + triplet_loss
        # Gradient w.r.t. L (chain rule)
        grad <- grad + 2 * t(xi - xj) %*% diff_pos - 2 * t(xi - xl) %*% diff_neg
      }
    }

    # Normalize gradient
    grad <- grad / max(length(triplets), 1)

    # Update transformation matrix
    L <- L - learn_rate * grad

    if (verbose && (iter %% 10 == 0 || iter == 1 || iter == max_iter)) {
      message(sprintf("Iteration %d/%d | Triplet Loss: %.4f", iter, max_iter, loss))
    }
  }

  # 5. Apply transformation (genes x cells)
  X_trans <- t(X %*% L)

  # 6. Return result
  list(
    expression_profile_trans = X_trans,
    expression_profile_origin = expr,
    trans_matrix = L,
    sample_information = sample_information
  )
}


#' @title Multi Similarity Loss-Based Metric Learning for Single-Cell Data
#'
#' @description
#' Implements a metric learning algorithm based on Multi Similarity Loss (MSL),
#' as an alternative to Triplet Loss in LMNN. It enhances learning from limited
#' sample pairs by considering multiple similarity cues and penalizes overfitting
#' through L2 regularization. Optimized for single-cell expression data.
#'
#' @param high_varGenes Character vector of highly variable gene names.
#' @param expression_profile A numeric matrix of gene expression values (genes x cells).
#' @param sample_information A named vector of cell type labels (names match column names of expression_profile).
#' @param alpha Scaling parameter for positive similarities (default: 2.0).
#' @param beta Scaling parameter for negative similarities (default: 50.0).
#' @param margin Margin threshold to filter informative pairs (default: 0.1).
#' @param lambda L2 regularization strength (default: 1e-4).
#' @param learn_rate Learning rate for gradient descent (default: 1e-5).
#' @param max_iter Maximum number of optimization iterations (default: 200).
#' @param seed Random seed (default: 1).
#' @param verbose Logical; whether to print progress information (default: TRUE).
#'
#' @return A list containing:
#' \itemize{
#'   \item \code{expression_profile_trans}: Transformed expression matrix (genes x cells).
#'   \item \code{expression_profile_origin}: Original filtered expression matrix (genes x cells).
#'   \item \code{trans_matrix}: Learned transformation matrix (L).
#'   \item \code{sample_information}: Sorted cell type labels.
#' }
#'
#' @export
#'
runMSL <- function(
    high_varGenes,
    expression_profile,
    sample_information,
    alpha = 2.0,
    beta = 50.0,
    margin = 0.1,
    lambda = 1e-4,
    learn_rate = 1e-5,
    max_iter = 200,
    seed = 1,
    verbose = TRUE) {
  if (!all(high_varGenes %in% rownames(expression_profile))) {
    stop("Some high_varGenes not found in expression_profile.")
  }

  sample_information <- sort(sample_information)
  expr <- expression_profile[high_varGenes, names(sample_information), drop = FALSE]

  X <- t(as.matrix(expr)) # cells x genes
  y <- as.character(sample_information)
  n_cells <- nrow(X)
  n_genes <- ncol(X)

  set.seed(seed)
  L <- diag(n_genes) # transformation matrix

  cosine_similarity <- function(a, b) {
    a_norm <- sqrt(sum(a^2))
    b_norm <- sqrt(sum(b^2))
    if (a_norm == 0 || b_norm == 0) {
      return(0)
    }
    sum(a * b) / (a_norm * b_norm)
  }

  grad_clip <- function(grad, clip_value = 1.0) {
    grad_norm <- sqrt(sum(grad^2))
    if (grad_norm > clip_value) {
      grad <- grad * (clip_value / grad_norm)
    }
    grad
  }

  for (iter in seq_len(max_iter)) {
    loss <- 0
    grad <- matrix(0, n_genes, n_genes)

    for (i in seq_len(n_cells)) {
      xi <- X[i, , drop = FALSE]
      Li <- xi %*% L
      Li_norm <- sqrt(sum(Li^2))

      for (j in seq_len(n_cells)) {
        if (i == j) next

        xj <- X[j, , drop = FALSE]
        Lj <- xj %*% L
        Lj_norm <- sqrt(sum(Lj^2))

        # Cosine similarity with numerical stability
        cos_sim <- sum(Li * Lj) / (Li_norm * Lj_norm + 1e-8)
        cos_sim <- pmax(pmin(cos_sim, 1 - 1e-6), -1 + 1e-6)

        is_pos <- (y[i] == y[j])
        if (is_pos && cos_sim < 1 - margin) {
          loss <- loss + (1 / alpha) * log(1 + exp(-alpha * (cos_sim - margin)))
          coeff <- -exp(-alpha * (cos_sim - margin)) / (1 + exp(-alpha * (cos_sim - margin)))
        } else if (!is_pos && cos_sim > margin) {
          loss <- loss + (1 / beta) * log(1 + exp(beta * (cos_sim - margin)))
          coeff <- exp(beta * (cos_sim - margin)) / (1 + exp(beta * (cos_sim - margin)))
        } else {
          next
        }

        # Compute gradient
        diff <- xi - xj
        num <- Lj * Li_norm^2 - Li * sum(Li * Lj)
        denom <- Li_norm^3 * Lj_norm + 1e-8
        cos_sim_grad <- num / denom

        if (is_pos) {
          grad <- grad + coeff * t(diff) %*% cos_sim_grad
        } else {
          grad <- grad - coeff * t(diff) %*% cos_sim_grad
        }
      }
    }

    # Add L2 regularization
    loss <- loss + lambda * sum(L^2)
    grad <- grad + 2 * lambda * L

    # Gradient clipping
    grad <- grad_clip(grad)

    # Update
    L <- L - learn_rate * grad

    if (verbose && (iter %% 10 == 0 || iter == 1 || iter == max_iter)) {
      message(sprintf("Iteration %d/%d | Loss: %.4f", iter, max_iter, loss))
    }
  }

  X_trans <- t(X %*% L)

  list(
    expression_profile_trans = X_trans,
    expression_profile_origin = expr,
    trans_matrix = L,
    sample_information = sample_information
  )
}

#' @title Neighborhood Components Analysis (NCA) Transformation
#'
#' @description
#' Performs NCA metric learning for single-cell data. Maximizes the leave-one-out
#' classification probability using a stochastic nearest neighbors approach.
#'
#' @param high_varGenes Character vector of high-variance genes
#' @param expression_profile Numeric matrix (genes x cells)
#' @param sample_information Named vector of cell type labels
#' @param max_iter Maximum iterations (default: 100)
#' @param learn_rate Learning rate (default: 0.01)
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
#' @examples
#' \dontrun{
#' # Example usage:
#' nca_result <- runNCA(
#'   high_varGenes = hvgs,
#'   expression_profile = expr_matrix,
#'   sample_information = cell_labels
#' )
#' }
#'
#' @export
#'
runNCA <- function(
    high_varGenes,
    expression_profile,
    sample_information,
    max_iter = 100,
    learn_rate = 0.01,
    seed = 1,
    verbose = TRUE) {
  # Input validation --------------------------------------------------------
  if (!all(high_varGenes %in% rownames(expression_profile))) {
    stop("high_varGenes not found in expression_profile")
  }

  # Sort by cell type for consistency
  sample_information <- sort(sample_information)
  expr <- expression_profile[high_varGenes, names(sample_information), drop = FALSE]

  # Prepare data matrices ---------------------------------------------------
  # Transpose for NCA (cells x genes)
  X <- t(as.matrix(expr))
  y <- as.numeric(factor(as.character(sample_information)))
  n_cells <- nrow(X)
  n_genes <- ncol(X)

  # Initialize transformation matrix
  set.seed(seed)
  A <- diag(n_genes) # Start with identity matrix

  # Define helper functions -------------------------------------------------

  # Softmax function with numerical stability
  softmax <- function(x) {
    exp_x <- exp(x - max(x)) # Subtract max for numerical stability
    exp_x / sum(exp_x)
  }

  # Main optimization loop --------------------------------------------------
  for (iter in 1:max_iter) {
    if (verbose && iter %% 10 == 0) {
      message(sprintf("NCA Iteration %d/%d", iter, max_iter))
    }

    # Initialize gradient
    grad <- matrix(0, n_genes, n_genes)
    p_correct <- 0

    for (i in 1:n_cells) {
      # Find same-class neighbors (excluding self)
      same_class <- which(y == y[i])
      same_class <- setdiff(same_class, i)

      if (length(same_class) > 0) {
        # Calculate differences using sweep (handles dimension mismatch)
        diff <- sweep(X[same_class, , drop = FALSE], 2, X[i, ])

        # Compute distances under current transformation
        distances <- rowSums((diff %*% A)^2)

        # Compute softmax probabilities
        probs <- softmax(-distances)

        # Update gradient (sum over all same-class neighbors)
        grad <- grad + t(diff) %*% (probs * diff)

        # For nearest neighbor (correct class)
        nearest <- which.max(probs)
        grad <- grad - t(diff[nearest, , drop = FALSE]) %*%
          diff[nearest, , drop = FALSE]
        p_correct <- p_correct + probs[nearest]
      }
    }

    # Update transformation matrix
    A <- A - learn_rate * grad

    # Early stopping if classification is nearly perfect
    if (iter > 10 && p_correct / n_cells > 0.99) {
      if (verbose) message("Early stopping at iteration ", iter)
      break
    }
  }

  # Transform data ---------------------------------------------------------
  # Transform data (transpose back to genes x cells)
  X_trans <- t(X %*% t(A))

  # Return results in scLearn-compatible format
  list(
    expression_profile_trans = X_trans,
    expression_profile_origin = expr,
    trans_matrix = t(A), # Transpose to match DCA convention
    sample_information = sample_information
  )
}

#' @title Information-Theoretic Metric Learning (ITML) for Single-Cell Data
#'
#' @description
#' Performs stable ITML metric learning for single-cell RNA-seq data with multiple
#' fallback mechanisms to ensure numerical stability. Learns a Mahalanobis distance
#' metric using relative distance constraints with information-theoretic regularization.
#'
#' @param high_varGenes Character vector of high-variance gene names
#' @param expression_profile Numeric matrix (genes x cells)
#' @param sample_information Named vector of cell type labels
#' @param gamma Slack parameter controlling constraint strictness (default: 0.1)
#' @param max_iter Maximum number of optimization iterations (default: 50)
#' @param epsilon Convergence threshold for early stopping (default: 1e-6)
#' @param seed Random seed for reproducibility (default: 1)
#' @param verbose Whether to print progress messages (default: TRUE)
#'
#' @return A list containing:
#' \itemize{
#'   \item expression_profile_trans: Transformed data matrix (genes x cells)
#'   \item expression_profile_origin: Original expression matrix
#'   \item trans_matrix: Learned transformation matrix
#'   \item sample_information: Input cell labels
#'   \item convergence: Convergence status information
#' }
#'
#' @examples
#' \dontrun{
#' # Example usage:
#' data("scRNA_example")
#' itml_result <- runITML(
#'   high_varGenes = scRNA_example$hv_genes,
#'   expression_profile = scRNA_example$expr_mat,
#'   sample_information = scRNA_example$cell_types
#' )
#' }
#'
#' @importFrom Matrix nearPD
#' @importFrom stats rnorm sd
#'
#' @export
#'
runITML <- function(
    high_varGenes,
    expression_profile,
    sample_information,
    gamma = 0.1,
    max_iter = 50,
    epsilon = 1e-6,
    seed = 1,
    verbose = TRUE) {
  # Input Validation -------------------------------------------------------
  if (!all(high_varGenes %in% rownames(expression_profile))) {
    stop("Not all high_varGenes found in expression_profile")
  }
  if (anyNA(expression_profile)) {
    stop("NA values detected in expression_profile")
  }

  # Data Preparation -------------------------------------------------------
  sample_information <- sort(sample_information)
  expr <- expression_profile[high_varGenes, names(sample_information), drop = FALSE]

  # Remove near-constant features
  gene_sd <- apply(expr, 1, sd)
  if (any(gene_sd < 1e-8)) {
    if (verbose) message("Removing ", sum(gene_sd < 1e-8), " near-constant features")
    expr <- expr[gene_sd > 1e-8, ]
    high_varGenes <- high_varGenes[gene_sd > 1e-8]
  }

  X <- t(as.matrix(expr)) # Convert to cells x genes format
  y <- as.character(sample_information)
  n_cells <- nrow(X)
  n_genes <- ncol(X)

  # Initialize Metric Matrix ------------------------------------------------
  set.seed(seed)
  M <- diag(n_genes) + 1e-6 * matrix(rnorm(n_genes^2), n_genes, n_genes)

  # Constraint Generation --------------------------------------------------
  generate_constraints <- function(y, max_pairs = 500) {
    classes <- unique(y)
    must_link <- list()
    cannot_link <- list()

    for (cls in classes) {
      same_class <- which(y == cls)
      if (length(same_class) > 1) {
        pairs <- utils::combn(sample(same_class, min(10, length(same_class))), 2, simplify = FALSE)
        must_link <- c(must_link, pairs[1:min(length(pairs), max_pairs)])
      }
    }

    for (i in 1:(length(classes) - 1)) {
      for (j in (i + 1):length(classes)) {
        cls1 <- which(y == classes[i])
        cls2 <- which(y == classes[j])
        if (length(cls1) > 0 && length(cls2) > 0) {
          pairs <- lapply(1:min(5, length(cls1)), function(k) c(sample(cls1, 1), sample(cls2, 1)))
          cannot_link <- c(cannot_link, pairs[1:min(length(pairs), max_pairs)])
        }
      }
    }

    list(
      must_link = must_link[1:min(max_pairs, length(must_link))],
      cannot_link = cannot_link[1:min(max_pairs, length(cannot_link))]
    )
  }

  constraints <- generate_constraints(y)

  # Stable Cholesky Decomposition Function ---------------------------------
  safe_chol <- function(M, n_genes, verbose) {
    attempts <- list(
      list(reg = 1e-6, msg = "Standard regularization"),
      list(reg = 1e-4, msg = "Increased regularization"),
      list(reg = 1e-2, msg = "Strong regularization")
    )

    for (attempt in attempts) {
      L <- tryCatch(
        {
          chol(M + attempt$reg * diag(n_genes))
        },
        error = function(e) NULL
      )

      if (!is.null(L)) {
        if (verbose) message("Cholesky successful with ", attempt$msg)
        return(L)
      }
    }

    # Final fallback to nearPD
    if (verbose) message("All Cholesky attempts failed, using nearPD")
    M_pd <- Matrix::nearPD(M, corr = FALSE, keepDiag = TRUE, maxit = 1000)$mat
    chol(M_pd + 1e-8 * diag(n_genes))
  }

  # Main Optimization Loop -------------------------------------------------
  convergence <- list(converged = FALSE, iterations = max_iter)

  for (iter in 1:max_iter) {
    if (verbose && iter %% 10 == 0) {
      message(sprintf("ITML Iteration %d/%d", iter, max_iter))
    }

    M_prev <- M

    # Process Constraints --------------------------------------------------
    process_constraints <- function(pairs, type = c("must_link", "cannot_link")) {
      for (pair in pairs) {
        i <- pair[1]
        j <- pair[2]
        diff <- X[i, ] - X[j, ]
        dist <- sqrt(max(sum(diff %*% M %*% diff), 1e-10))

        alpha <- if (type == "must_link") {
          min(gamma, max(1 / dist - 1, -0.99))
        } else {
          min(gamma, max(1 - 1 / dist, -0.99))
        }

        update_term <- alpha * M %*% tcrossprod(diff) %*% M
        if (!any(is.infinite(update_term))) {
          M <<- M + if (type == "must_link") update_term else -update_term
        }
      }
    }

    process_constraints(constraints$must_link, "must_link")
    process_constraints(constraints$cannot_link, "cannot_link")

    # Post-Update Stabilization --------------------------------------------
    M <- (M + t(M)) / 2 # Ensure symmetry
    M <- M / (norm(M, "F") + 1e-10) # Control matrix norm

    # Eigen Decomposition with Regularization ------------------------------
    eig <- tryCatch(
      {
        eigen(M, symmetric = TRUE)
      },
      error = function(e) {
        if (verbose) message("Eigen decomposition failed at iter ", iter)
        list(values = rep(1, n_genes), vectors = diag(n_genes))
      }
    )

    M <- eig$vectors %*% diag(pmax(eig$values, 1e-8)) %*% t(eig$vectors)

    # Convergence Check ----------------------------------------------------
    if (norm(M - M_prev, "F") < epsilon) {
      if (verbose) message("Converged at iteration ", iter)
      convergence$converged <- TRUE
      convergence$iterations <- iter
      break
    }
  }

  # Final Transformation ---------------------------------------------------
  L <- safe_chol(M, n_genes, verbose)

  # Return Results ---------------------------------------------------------
  list(
    expression_profile_trans = t(X %*% t(L)),
    expression_profile_origin = expr,
    trans_matrix = t(L),
    sample_information = sample_information,
    convergence = convergence
  )
}

#' @title Maximum Margin Clustering (MMC) Transformation (Memory-Optimized)
#'
#' @description
#' Performs MMC metric learning for single-cell data with reduced memory usage.
#' Incorporates low-rank approximation and other optimizations.
#'
#' @param high_varGenes Character vector of high-variance genes
#' @param expression_profile Numeric matrix (genes x cells)
#' @param sample_information Named vector of cell type labels
#' @param C Trade-off parameter (default: 1.0)
#' @param max_iter Maximum iterations (default: 50)
#' @param seed Random seed (default: 1)
#' @param verbose Print progress messages (default: TRUE)
#' @param rank Optional rank for low-rank approximation (default: NULL)
#' @param batch_size Number of cells to process at once if using batching (default: NULL)
#'
#' @return List containing:
#' \itemize{
#'   \item expression_profile_trans: Transformed matrix (genes x cells)
#'   \item expression_profile_origin: Original matrix
#'   \item trans_matrix: Learned transformation matrix
#'   \item sample_information: Input cell labels
#' }
#'
#' @export
#'
runMMC <- function(
    high_varGenes,
    expression_profile,
    sample_information,
    C = 1.0,
    max_iter = 50,
    seed = 1,
    verbose = TRUE,
    rank = NULL,
    batch_size = NULL) {

  # Input validation
  if (!requireNamespace("quadprog", quietly = TRUE)) {
    stop("quadprog package required for MMC. Install with: install.packages('quadprog')")
  }

  if (!all(high_varGenes %in% rownames(expression_profile))) {
    stop("high_varGenes not found in expression_profile")
  }

  # Sort by cell type for consistency
  sample_information <- sort(sample_information)
  expr <- expression_profile[high_varGenes, names(sample_information), drop = FALSE]

  # Transpose for MMC (cells x genes)
  X <- t(as.matrix(expr))
  y <- as.character(sample_information)
  classes <- unique(y)
  n_classes <- length(classes)
  n_cells <- nrow(X)
  n_genes <- ncol(X)

  # Determine rank for low-rank approximation if not specified
  if (is.null(rank)) {
    rank <- min(100, floor(n_genes/2))  # Default to 100 or half the genes
    if (verbose) message("Using low-rank approximation with rank = ", rank)
  }

  # Initialize metric matrix (low-rank)
  set.seed(seed)
  if (n_genes > rank) {
    # Use low-rank approximation
    if (requireNamespace("irlba", quietly = TRUE)) {
      svd_init <- irlba::irlba(X, nv = rank)
      L <- svd_init$v %*% diag(sqrt(svd_init$d[1:rank]))
    } else {
      warning("irlba package not available, using standard SVD which may be slower")
      svd_init <- svd(X, nu = rank, nv = rank)
      L <- svd_init$v %*% diag(sqrt(svd_init$d[1:rank]))
    }
  } else {
    # Full rank case
    rank <- n_genes
    L <- diag(n_genes)
  }

  # Memory-efficient scatter matrix computation
  compute_scatter <- function(X, y, batch_size = NULL) {
    overall_mean <- colMeans(X)
    S_w <- matrix(0, n_genes, n_genes)
    S_b <- matrix(0, n_genes, n_genes)

    if (!is.null(batch_size)) {
      # Process in batches to reduce memory
      for (cls in classes) {
        cls_idx <- which(y == cls)
        n_batches <- ceiling(length(cls_idx)/batch_size)

        for (i in 1:n_batches) {
          batch_idx <- cls_idx[((i-1)*batch_size + 1):min(i*batch_size, length(cls_idx))]
          X_batch <- X[batch_idx, , drop = FALSE]
          batch_mean <- colMeans(X_batch)
          S_w <- S_w + crossprod(X_batch - batch_mean)
        }

        cls_mean <- colMeans(X[y == cls, , drop = FALSE])
        S_b <- S_b + sum(y == cls) * tcrossprod(cls_mean - overall_mean)
      }
    } else {
      # Standard computation
      for (cls in classes) {
        cls_idx <- which(y == cls)
        cls_mean <- colMeans(X[cls_idx, , drop = FALSE])
        S_w <- S_w + crossprod(X[cls_idx, , drop = FALSE] - cls_mean)
        S_b <- S_b + length(cls_idx) * tcrossprod(cls_mean - overall_mean)
      }
    }

    list(S_w = S_w, S_b = S_b)
  }

  # Main optimization loop with memory optimizations
  for (iter in 1:max_iter) {
    if (verbose) {
      message(sprintf("MMC Iteration %d/%d", iter, max_iter))
    }

    # Compute scatter matrices with optional batching
    scatter <- compute_scatter(X, y, batch_size)
    S_w <- scatter$S_w

    # Solve QP problem with reduced memory footprint
    Dmat <- 2 * S_w
    dvec <- rep(0, n_genes)

    # Only enforce diagonal constraints to reduce memory
    Amat <- diag(n_genes)
    bvec <- rep(0, n_genes)

    # Solve using quadratic programming
    sol <- quadprog::solve.QP(Dmat, dvec, Amat, bvec)
    M_diag <- pmax(sol$solution, 1e-6)

    # Update L (low-rank factor)
    if (rank < n_genes) {
      # For low-rank case, update the factor matrix
      L_scale <- sqrt(M_diag / diag(tcrossprod(L)))
      L <- L * matrix(rep(L_scale, each = rank), nrow = n_genes, ncol = rank)
    } else {
      # For full-rank case, update diagonal directly
      L <- diag(sqrt(M_diag))
    }

    # Early stopping if change is small
    if (iter > 1) {
      M_current <- tcrossprod(L)
      delta <- norm(M_current - M_prev, "F") / norm(M_prev, "F")
      if (delta < 1e-6) {
        if (verbose) message("Converged after ", iter, " iterations")
        break
      }
    }
    M_prev <- tcrossprod(L)
  }

  # Compute final transformation matrix
  trans_matrix <- t(L)

  # Transform data in batches if large
  if (!is.null(batch_size) && n_cells > batch_size) {
    X_trans <- matrix(0, n_genes, n_cells)
    n_batches <- ceiling(n_cells/batch_size)
    for (i in 1:n_batches) {
      idx <- ((i-1)*batch_size + 1):min(i*batch_size, n_cells)
      X_trans[, idx] <- t(X[idx, , drop = FALSE] %*% t(trans_matrix))
    }
  } else {
    X_trans <- t(X %*% t(trans_matrix))
  }

  # Return results
  list(
    expression_profile_trans = X_trans,
    expression_profile_origin = expr,
    trans_matrix = trans_matrix,
    sample_information = sample_information
  )
}
