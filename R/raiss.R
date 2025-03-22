#' Core RAISS implementation for a single LD matrix
#'
#' @param ref_panel A data frame containing 'chrom', 'pos', 'variant_id', 'A1', and 'A2'.
#' @param known_zscores A data frame containing 'chrom', 'pos', 'variant_id', 'A1', 'A2', and 'z' values.
#' @param LD_matrix A square matrix of dimension equal to the number of rows in ref_panel.
#' @param lamb Regularization term added to the diagonal of the LD_matrix.
#' @param rcond Threshold for filtering eigenvalues in the pseudo-inverse computation.
#' @param R2_threshold R square threshold below which SNPs are filtered from the output.
#' @param minimum_ld Minimum LD score threshold for SNP filtering.
#' @param verbose Logical indicating whether to print progress information.
#'
#' @return A list containing filtered and unfiltered results, and filtered LD matrix.
#' @importFrom dplyr arrange
#' @export
raiss_single_matrix <- function(ref_panel, known_zscores, LD_matrix, lamb = 0.01, rcond = 0.01,
                                R2_threshold = 0.6, minimum_ld = 5, verbose = TRUE) {
  # Check that ref_panel and known_zscores are both increasing in terms of pos
  if (is.unsorted(ref_panel$pos) || is.unsorted(known_zscores$pos)) {
    stop("ref_panel and known_zscores must be in increasing order of pos.")
  }

  # Convert LD_matrix to matrix if it's a data frame
  if (is.data.frame(LD_matrix)) {
    LD_matrix <- as.matrix(LD_matrix)
  }

  # Define knowns and unknowns
  knowns_id <- intersect(known_zscores$variant_id, ref_panel$variant_id)
  knowns <- which(ref_panel$variant_id %in% knowns_id)
  unknowns <- which(!ref_panel$variant_id %in% knowns_id)

  # Handle edge cases
  if (length(knowns) == 0) {
    if (verbose) message("No known variants found, cannot perform imputation.")
    return(NULL)
  }

  if (length(unknowns) == 0) {
    if (verbose) message("No unknown variants to impute, returning known variants.")
    return(list(
      result_nofilter = known_zscores,
      result_filter = known_zscores,
      LD_mat = LD_matrix
    ))
  }

  # Extract zt, sig_t, and sig_i_t
  zt <- known_zscores$z[knowns]
  sig_t <- LD_matrix[knowns, knowns, drop = FALSE]
  sig_i_t <- LD_matrix[unknowns, knowns, drop = FALSE]

  # Call raiss_model
  results <- raiss_model(zt, sig_t, sig_i_t, lamb, rcond)

  # Format the results
  results <- format_raiss_df(results, ref_panel, unknowns)

  # Filter output
  results <- filter_raiss_output(results, R2_threshold, minimum_ld, verbose)

  # Merge with known z-scores
  result_nofilter <- merge_raiss_df(results$zscores_nofilter, known_zscores) %>% arrange(pos)
  result_filter <- merge_raiss_df(results$zscores, known_zscores) %>% arrange(pos)

  # Filter out variants not included in the imputation result
  filtered_out_variant <- setdiff(ref_panel$variant_id, result_filter$variant_id)

  # Update the LD matrix excluding filtered variants
  LD_extract_filtered <- if (length(filtered_out_variant) > 0) {
    filtered_out_id <- match(filtered_out_variant, ref_panel$variant_id)
    as.matrix(LD_matrix)[-filtered_out_id, -filtered_out_id]
  } else {
    as.matrix(LD_matrix)
  }

  # Return results
  return(list(
    result_nofilter = result_nofilter,
    result_filter = result_filter,
    LD_mat = LD_extract_filtered
  ))
}

#' Robust and accurate imputation from summary statistics
#'
#' This function is a part of the statistical library for SNP imputation from:
#' https://gitlab.pasteur.fr/statistical-genetics/raiss/-/blob/master/raiss/stat_models.py
#' It is R implementation of the imputation model described in the paper by Bogdan Pasaniuc,
#' Noah Zaitlen, et al., titled "Fast and accurate imputation of summary
#' statistics enhances evidence of functional enrichment", published in
#' Bioinformatics in 2014.
#'
#' This function can process either a single LD matrix or a list of LD matrices for different blocks.
#' For a list of matrices, it processes each block separately and combines the results.
#'
#' @param ref_panel A data frame containing 'chrom', 'pos', 'variant_id', 'A1', and 'A2'.
#' @param known_zscores A data frame containing 'chrom', 'pos', 'variant_id', 'A1', 'A2', and 'z' values.
#' @param LD_matrix Either a square matrix or a list of matrices for LD blocks.
#' @param variant_indices Optional data frame mapping variant IDs to block IDs when LD_matrix is a list.
#' @param lamb Regularization term added to the diagonal of the LD_matrix.
#' @param rcond Threshold for filtering eigenvalues in the pseudo-inverse computation.
#' @param R2_threshold R square threshold below which SNPs are filtered from the output.
#' @param minimum_ld Minimum LD score threshold for SNP filtering.
#' @param verbose Logical indicating whether to print progress information.
#'
#' @return A list containing filtered and unfiltered results, and filtered LD matrix.
#' @importFrom dplyr arrange bind_rows
#' @export
raiss <- function(ref_panel, known_zscores, LD_matrix, variant_indices = NULL, lamb = 0.01, rcond = 0.01,
                  R2_threshold = 0.6, minimum_ld = 5, verbose = TRUE) {
  # Check if LD_matrix is a list or a single matrix
  is_list_of_matrices <- is.list(LD_matrix) && !is.data.frame(LD_matrix) && !is.matrix(LD_matrix)

  # For single matrix, call the core function directly
  if (!is_list_of_matrices) {
    if (verbose) message("Processing single LD matrix...")
    return(raiss_single_matrix(
      ref_panel, known_zscores, LD_matrix,
      lamb, rcond, R2_threshold, minimum_ld, verbose
    ))
  }

  # For list of matrices, process each block
  if (is.null(variant_indices)) {
    stop("variant_indices must be provided when LD_matrix is a list.")
  }

  if (verbose) message("Processing multiple LD blocks...")

  # Prepare list to collect results from each block
  results_list <- list()

  # Get unique block IDs
  block_ids <- unique(variant_indices$block_id)

  # Process each block
  for (block_id in block_ids) {
    if (verbose) message(paste("Processing block", block_id, "of", length(block_ids)))

    # Get variants in this block
    block_variant_ids <- variant_indices$variant_id[variant_indices$block_id == block_id]

    # Subset ref_panel and LD_matrix for this block
    block_indices <- match(block_variant_ids, ref_panel$variant_id)
    block_ref_panel <- ref_panel[block_indices, ]
    block_LD_matrix <- LD_matrix$ld_matrices[[block_id]]

    # Check dimensions match
    if (nrow(block_LD_matrix) != nrow(block_ref_panel)) {
      stop(paste("Block", block_id, ": LD matrix dimension does not match number of variants in reference panel"))
    }

    # Process the block using the core function
    block_result <- raiss_single_matrix(
      block_ref_panel, known_zscores, block_LD_matrix,
      lamb, rcond, R2_threshold, minimum_ld,
      verbose = FALSE
    )

    # Skip if block returned NULL (no known variants)
    if (!is.null(block_result)) {
      results_list[[block_id]] <- block_result
    }
  }

  # If no valid blocks were processed
  if (length(results_list) == 0) {
    if (verbose) message("No blocks could be processed. Check that known_zscores overlap with variants in the blocks.")
    return(NULL)
  }

  # Combine results from all blocks
  nofilter_results <- lapply(results_list, function(x) x$result_nofilter)
  filter_results <- lapply(results_list, function(x) x$result_filter)
  ld_filtered_list <- lapply(results_list, function(x) x$LD_mat)

  # Combine into data frames
  result_nofilter <- dplyr::bind_rows(nofilter_results) %>% arrange(pos)
  result_filter <- dplyr::bind_rows(filter_results) %>% arrange(pos)

  # Return combined results
  return(list(
    result_nofilter = result_nofilter,
    result_filter = result_filter,
    LD_mat = ld_filtered_list
  ))
}

#' @param zt Vector of known z scores.
#' @param sig_t Matrix of known linkage disequilibrium (LD) correlation.
#' @param sig_i_t Correlation matrix with rows corresponding to unknown SNPs (to impute)
#'               and columns to known SNPs.
#' @param lamb Regularization term added to the diagonal of the sig_t matrix.
#' @param rcond Threshold for filtering eigenvalues in the pseudo-inverse computation.
#' @param batch Boolean indicating whether batch processing is used.
#'
#' @return A list containing the variance 'var', estimation 'mu', LD score 'raiss_ld_score',
#'         condition number 'condition_number', and correctness of inversion
#'         'correct_inversion'.
#' @noRd
raiss_model <- function(zt, sig_t, sig_i_t, lamb = 0.01, rcond = 0.01, batch = TRUE, report_condition_number = FALSE) {
  sig_t_inv <- invert_mat_recursive(sig_t, lamb, rcond)
  if (!is.numeric(zt) || !is.numeric(sig_t) || !is.numeric(sig_i_t)) {
    stop("zt, sig_t, and sig_i_t must be numeric.")
  }
  if (batch) {
    condition_number <- if (report_condition_number) rep(kappa(sig_t, exact = T, norm = "2"), nrow(sig_i_t)) else NA
    correct_inversion <- rep(check_inversion(sig_t, sig_t_inv), nrow(sig_i_t))
  } else {
    condition_number <- if (report_condition_number) kappa(sig_t, exact = T, norm = "2") else NA
    correct_inversion <- check_inversion(sig_t, sig_t_inv)
  }

  var_raiss_ld_score <- compute_var(sig_i_t, sig_t_inv, lamb, batch)
  var <- var_raiss_ld_score$var
  raiss_ld_score <- var_raiss_ld_score$raiss_ld_score

  mu <- compute_mu(sig_i_t, sig_t_inv, zt)
  var_norm <- var_in_boundaries(var, lamb)

  R2 <- ((1 + lamb) - var_norm)
  mu <- mu / sqrt(R2)

  return(list(var = var_norm, mu = mu, raiss_ld_score = raiss_ld_score, condition_number = condition_number, correct_inversion = correct_inversion))
}

#' @param imp is the output of raiss_model()
#' @param ref_panel is a data frame with columns 'chrom', 'pos', 'variant_id', 'ref', and 'alt'.
#' @noRd
format_raiss_df <- function(imp, ref_panel, unknowns) {
  result_df <- data.frame(
    chrom = ref_panel[unknowns, "chrom"],
    pos = ref_panel[unknowns, "pos"],
    variant_id = ref_panel[unknowns, "variant_id"],
    A1 = ref_panel[unknowns, "A1"],
    A2 = ref_panel[unknowns, "A2"],
    z = imp$mu,
    Var = imp$var,
    raiss_ld_score = imp$raiss_ld_score,
    condition_number = imp$condition_number,
    correct_inversion = imp$correct_inversion
  )

  # Specify the column order
  column_order <- c(
    "chrom", "pos", "variant_id", "A1", "A2", "z", "Var", "raiss_ld_score", "condition_number",
    "correct_inversion"
  )

  # Reorder the columns
  result_df <- result_df[, column_order]
  return(result_df)
}

merge_raiss_df <- function(raiss_df, known_zscores) {
  # Merge the data frames
  merged_df <- merge(raiss_df, known_zscores, by = c("chrom", "pos", "variant_id", "A1", "A2"), all = TRUE)

  # Identify rows that came from known_zscores
  from_known <- !is.na(merged_df$z.y) & is.na(merged_df$z.x)

  # Set Var to -1 and raiss_ld_score to Inf for these rows
  merged_df$Var[from_known] <- -1
  merged_df$raiss_ld_score[from_known] <- Inf

  # If there are overlapping columns (e.g., z.x and z.y), resolve them
  # For example, use z from known_zscores where available, otherwise use z from raiss_df
  merged_df$z <- ifelse(from_known, merged_df$z.y, merged_df$z.x)

  # Remove the extra columns resulted from the merge (e.g., z.x, z.y)
  merged_df <- merged_df[, !colnames(merged_df) %in% c("z.x", "z.y")]
  merged_df <- arrange(merged_df, pos)
  # assign imputed variants beta, se as NA to avoid confusion, since they are not imputed
  merged_df$beta[merged_df$Var == -1] <- NA
  merged_df$se[merged_df$Var == -1] <- NA
  return(merged_df)
}

filter_raiss_output <- function(zscores, R2_threshold = 0.6, minimum_ld = 5, verbose = TRUE) {
  # Reset the index and subset the data frame
  zscores <- zscores[, c("chrom", "pos", "variant_id", "A1", "A2", "z", "Var", "raiss_ld_score")]
  zscores$raiss_R2 <- 1 - zscores$Var

  # Count statistics before filtering
  NSNPs_bf_filt <- nrow(zscores)
  NSNPs_initial <- sum(zscores$raiss_R2 == 2.0)
  NSNPs_imputed <- sum(zscores$raiss_R2 != 2.0)
  NSNPs_ld_filt <- sum(zscores$raiss_ld_score < minimum_ld)
  NSNPs_R2_filt <- sum(zscores$raiss_R2 < R2_threshold)

  # Apply filters
  zscores_nofilter <- zscores
  zscores <- zscores[zscores$raiss_R2 > R2_threshold & zscores$raiss_ld_score >= minimum_ld, ]
  NSNPs_af_filt <- nrow(zscores)

  # Print report
  if (verbose) {
    message("IMPUTATION REPORT\n")
    message("Number of SNPs:\n")
    message("before filter:", NSNPs_bf_filt, "\n")
    message("not imputed:", NSNPs_initial, "\n")
    message("imputed:", NSNPs_imputed, "\n")
    message("filtered because of LD:", NSNPs_ld_filt, "\n")
    message("filtered because of R2:", NSNPs_R2_filt, "\n")
    message("after filter:", NSNPs_af_filt, "\n")
  }
  return(zscore_list = list(zscores_nofilter = zscores_nofilter, zscores = zscores))
}

compute_mu <- function(sig_i_t, sig_t_inv, zt) {
  return(sig_i_t %*% (sig_t_inv %*% zt))
}

compute_var <- function(sig_i_t, sig_t_inv, lamb, batch = TRUE) {
  if (batch) {
    var <- (1 + lamb) - rowSums((sig_i_t %*% sig_t_inv) * sig_i_t)
    raiss_ld_score <- rowSums(sig_i_t^2)
  } else {
    var <- (1 + lamb) - (sig_i_t %*% (sig_t_inv %*% t(sig_i_t)))
    raiss_ld_score <- sum(sig_i_t^2)
  }
  return(list(var = var, raiss_ld_score = raiss_ld_score))
}

check_inversion <- function(sig_t, sig_t_inv) {
  return(all.equal(sig_t, sig_t %*% (sig_t_inv %*% sig_t), tolerance = 1e-5))
}

var_in_boundaries <- function(var, lamb) {
  var[var < 0] <- 0
  var[var > (0.99999 + lamb)] <- 1
  return(var)
}

invert_mat <- function(mat, lamb, rcond) {
  tryCatch(
    {
      # Modify the diagonal elements of mat
      diag(mat) <- 1 + lamb
      # Compute the pseudo-inverse
      mat_inv <- MASS::ginv(mat, tol = rcond)
      return(mat_inv)
    },
    error = function(e) {
      # Second attempt with updated lamb and rcond in case of an error
      diag(mat) <- 1 + lamb * 1.1
      mat_inv <- MASS::ginv(mat, tol = rcond * 1.1)
      return(mat_inv)
    }
  )
}

invert_mat_recursive <- function(mat, lamb, rcond) {
  tryCatch(
    {
      # Modify the diagonal elements of mat
      diag(mat) <- 1 + lamb
      # Compute the pseudo-inverse
      mat_inv <- MASS::ginv(mat, tol = rcond)
      return(mat_inv)
    },
    error = function(e) {
      # Recursive call with updated lamb and rcond in case of an error
      invert_mat(mat, lamb * 1.1, rcond * 1.1)
    }
  )
}

invert_mat_eigen <- function(mat, tol = 1e-3) {
  eigen_mat <- eigen(mat)
  L <- which(cumsum(eigen_mat$values) / sum(eigen_mat$values) > 1 - tol)[1]
  if (is.na(L)) {
    # all eigen values are extremely small
    stop("Cannot invert the input matrix because all its eigen values are negative or close to zero")
  }
  mat_inv <- eigen_mat$vectors[, 1:L] %*%
    diag(1 / eigen_mat$values[1:L]) %*%
    t(eigen_mat$vectors[, 1:L])

  return(mat_inv)
}
