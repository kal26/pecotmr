context("raiss")
library(tidyverse)
library(MASS)

generate_dummy_data <- function(seed=1, ref_panel_ordered=TRUE, known_zscores_ordered=TRUE) {
    set.seed(seed)

    n_variants <- 100
    ref_panel <- data.frame(
        chrom = rep(1, n_variants),
        pos = seq(1, n_variants * 10, 10),
        variant_id = paste0("rs", seq_len(n_variants)),
        A1 = sample(c("A", "T", "G", "C"), n_variants, replace = TRUE),
        A2 = sample(c("A", "T", "G", "C"), n_variants, replace = TRUE)
    )

    n_known <- 50
    known_zscores <- data.frame(
        chrom = rep(1, n_known),
        pos = sample(ref_panel$pos, n_known),
        variant_id = sample(ref_panel$variant_id, n_known),
        A1 = sample(c("A", "T", "G", "C"), n_known, replace = TRUE),
        A2 = sample(c("A", "T", "G", "C"), n_known, replace = TRUE),
        z = rnorm(n_known)
    )

    LD_matrix <- matrix(rnorm(n_variants^2), nrow = n_variants, ncol = n_variants)
    diag(LD_matrix) <- 1 
    known_zscores <- if (known_zscores_ordered) known_zscores[order(known_zscores$pos),] else known_zscores
    ref_panel <- if (ref_panel_ordered) ref_panel else ref_panel[order(ref_panel$pos, decreasing = TRUE),]
    return(list(ref_panel=ref_panel, known_zscores=known_zscores, LD_matrix=LD_matrix))
}

test_that("Input validation for raiss works correctly", {
    input_data_ref_panel_unordered <- generate_dummy_data(ref_panel_ordered=FALSE)
    input_data_zscores_unordered <- generate_dummy_data(known_zscores_ordered=FALSE)
    expect_error(raiss(input_data_ref_panel_unordered$ref_panel, input_data$known_zscores, input_data$LD_matrix))
    expect_error(raiss(input_data$ref_panel, input_data_zscores_unordered$known_zscores, input_data$LD_matrix))
})

test_that("Default parameters for raiss work correctly", {
    # TODO - ask Gao about merging on removed columns
    input_data <- generate_dummy_data()
    result <- raiss(input_data$ref_panel, input_data$known_zscores, input_data$LD_matrix)
    expect_true(is.list(result))
})

test_that("Test Default Parameters for raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  result <- raiss_model(zt, sig_t, sig_i_t)

  expect_is(result, "list")
  expect_true(all(names(result) %in% c("var", "mu", "raiss_ld_score", "condition_number", "correct_inversion")))
})

test_that("Test with Different lamb Values for raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  lamb_values <- c(0.01, 0.05, 0.1)
  for (lamb in lamb_values) {
    result <- raiss_model(zt, sig_t, sig_i_t, lamb)
    expect_is(result, "list")
  }
})

test_that("Report Condition Number in raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  result_with_cn <- raiss_model(zt, sig_t, sig_i_t, report_condition_number = TRUE)
  result_without_cn <- raiss_model(zt, sig_t, sig_i_t, report_condition_number = FALSE)

  expect_is(result_with_cn, "list")
  expect_is(result_without_cn, "list")
})

test_that("Input Validation of raiss_model", {

  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)
  zt_invalid <- "not a numeric vector"
  sig_t_invalid <- "not a matrix"
  sig_i_t_invalid <- "not a matrix"

  expect_error(raiss_model(zt_invalid, sig_t, sig_i_t))
  expect_error(raiss_model(zt, sig_t_invalid, sig_i_t))
  expect_error(raiss_model(zt, sig_t, sig_i_t_invalid))
})

test_that("Boundary Conditions of raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  zt_empty <- numeric(0)
  sig_t_empty <- matrix(numeric(0), nrow = 0)
  sig_i_t_empty <- matrix(numeric(0), nrow = 0)

  expect_error(raiss_model(zt_empty, sig_t, sig_i_t))
  expect_error(raiss_model(zt, sig_t_empty, sig_i_t))
  expect_error(raiss_model(zt, sig_t, sig_i_t_empty))
})

test_that("Test with Different rcond Values for raiss_model", {
  zt <- c(1.2, 0.5)
  sig_t <- matrix(c(1, 0.5, 0.5, 1), nrow = 2)
  sig_i_t <- matrix(c(0.5, 0.2, 0.2, 0.5), nrow = 2)

  rcond_values <- c(0.01, 0.05, 0.1)
  for (rcond in rcond_values) {
    result <- raiss_model(zt, sig_t, sig_i_t, lamb = 0.01, rcond = rcond)
    expect_is(result, "list")
    expect_true(all(names(result) %in% c("var", "mu", "raiss_ld_score", "condition_number", "correct_inversion")))
  }
})

test_that("format_raiss_df returns correctly formatted data frame", {
  imp <- list(
    mu = rnorm(5),
    var = runif(5),
    raiss_ld_score = rnorm(5),
    condition_number = runif(5),
    correct_inversion = sample(c(TRUE, FALSE), 5, replace = TRUE)
  )
  
  ref_panel <- data.frame(
    chrom = sample(1:22, 10, replace = TRUE),
    pos = sample(1:10000, 10),
    variant_id = paste0("rs", 1:10),
    A1 = sample(c("A", "T", "G", "C"), 10, replace = TRUE),
    A2 = sample(c("A", "T", "G", "C"), 10, replace = TRUE)
  )

  unknowns <- sample(1:nrow(ref_panel), 5)

  result <- format_raiss_df(imp, ref_panel, unknowns)

  expect_true(is.data.frame(result))
  expect_equal(ncol(result), 10)
  expect_equal(colnames(result), c('chrom', 'pos', 'variant_id', 'A1', 'A2', 'z', 'Var', 'raiss_ld_score', 'condition_number', 'correct_inversion'))

  for (col in c('chrom', 'pos', 'variant_id', 'A1', 'A2')) {
    expect_equal(setNames(unlist(result[col]), NULL), unlist(ref_panel[unknowns, col, drop = TRUE]))
  }
  for (col in c('z', 'Var', 'raiss_ld_score', 'condition_number', 'correct_inversion')) {
    expected_col <- if (col == "z") "mu" else if (col == "Var") "var" else col
    expect_equal(setNames(unlist(result[col]), NULL), setNames(unlist(imp[expected_col]), NULL))
  }
})

test_that("Merge operation is correct for merge_raiss_df", {
    raiss_df_example <- data.frame(
        chrom = c("chr21", "chr22"),
        pos = c(123, 456),
        variant_id = c("var1", "var2"),
        A1 = c("A", "T"),
        A2 = c("T", "A"),
        z = c(0.5, 1.5),
        Var = c(0.2, 0.3),
        raiss_ld_score = c(10, 20),
        raiss_R2 = c(0.8, 0.7))

    known_zscores_example <- data.frame(
        chrom = c("chr21", "chr22"),
        pos = c(123, 456),
        variant_id = c("var1", "var2"),
        A1 = c("A", "T"),
        A2 = c("T", "A"),
        z = c(0.5, 1.5))

    merged_df <- merge_raiss_df(raiss_df_example, known_zscores_example)
    expect_equal(nrow(merged_df), 2)
    expect_true(all(c("chr21", "chr22") %in% merged_df$chrom))
})

generate_fro_test_data <- function(seed=1) {
    set.seed(seed)
    return(data.frame(
        chrom = paste0("chr", rep(22, 10)),
        pos = seq(1, 100, 10),
        variant_id = 1:10,
        A1 = rep("A", 10),
        A2 = rep("T", 10),
        z = rnorm(10),
        Var = runif(10, 0, 1),
        raiss_ld_score = rnorm(10, 5, 2)
    ))
}

test_that("Correct columns are selected in filter_raiss_output", {
    test_data <- generate_fro_test_data()
    output <- filter_raiss_output(test_data)$zscores
    expect_true(all(c('variant_id', 'A1', 'A2', 'z', 'Var', 'raiss_ld_score') %in% names(output)))
})

test_that("raiss_R2 is calculated correctly in filter_raiss_output", {
    test_data <- generate_fro_test_data()
    output <- filter_raiss_output(test_data)$zscores
    expected_R2 <- 1 - test_data[which(test_data$raiss_ld_score >= 5),]$Var
    expect_equal(output$raiss_R2, expected_R2[which(expected_R2 > 0.6)])
})

test_that("Filtering is applied correctly in filter_raiss_output", {
    test_data <- generate_fro_test_data()
    R2_threshold <- 0.6
    minimum_ld <- 5
    output <- filter_raiss_output(test_data, R2_threshold, minimum_ld)$zscores

    expect_true(all(output$raiss_R2 > R2_threshold))
    expect_true(all(output$raiss_ld_score >= minimum_ld))
})

test_that("Function returns the correct subset in filter_raiss_output", {
    test_data <- generate_fro_test_data()
    test_data$raiss_R2 <- 1 - test_data$Var
    output <- filter_raiss_output(test_data)$zscores

    manual_filter <- test_data[test_data$raiss_R2 > 0.6 & test_data$raiss_ld_score >= 5, ]

    expect_equal(nrow(output), nrow(manual_filter))
    expect_equal(sum(output$variant_id != manual_filter$variant_id), 0)
})

test_that("compute_mu basic functionality", {
    sig_i_t <- matrix(c(1, 2, 3, 4), nrow = 2)
    sig_t_inv <- matrix(c(5, 6, 7, 8), nrow = 2)
    zt <- matrix(c(9, 10, 11, 12), nrow = 2)

    expected_result <- matrix(c(517, 766, 625, 926), nrow = 2)
    result <- compute_mu(sig_i_t, sig_t_inv, zt)
    expect_equal(result, expected_result)
})

generate_mock_data_for_compute_var <- function(seed=1) {
    return(
        list(
            sig_i_t_1 = matrix(c(1, 2, 3, 4), nrow = 2),
            sig_t_inv_1 = matrix(c(5, 6, 7, 8), nrow = 2),
            lamb_1 = 0.5))
}

test_that("compute_var returns correct output for batch = TRUE", {
    input_data <- generate_mock_data_for_compute_var()
    result <- compute_var(input_data$sig_i_t_1, input_data$sig_t_inv_1, input_data$lamb_1, batch = TRUE)
    expect_is(result, "list")
    expect_length(result, 2)
    expect_true("var" %in% names(result))
    expect_true("raiss_ld_score" %in% names(result))
})

test_that("compute_var returns correct output for batch = FALSE", {
    input_data <- generate_mock_data_for_compute_var()
    result <- compute_var(input_data$sig_i_t_1, input_data$sig_t_inv_1, input_data$lamb_1, batch = FALSE)
    expect_is(result, "list")
    expect_length(result, 2)
    expect_true("var" %in% names(result))
    expect_true("raiss_ld_score" %in% names(result))
})

test_that("check_inversion correctly identifies inverse matrices in", {
  sig_t <- matrix(c(1, 2, 3, 4), nrow=2, ncol=2)
  sig_t_inv <- solve(sig_t)  
  expect_true(check_inversion(sig_t, sig_t_inv))
})

test_that("var_in_boundaries sets boundaries correctly", {
  lamb_test <- 0.05
  var <- c(-1, 0, 0.5, 1.04, 1.05)  

  result <- var_in_boundaries(var, lamb_test)

  expect_equal(result[1], 0)                   # Value less than 0 should be set to 0
  expect_equal(result[2], 0)                   # Value within lower boundary should remain unchanged
  expect_equal(result[3], 0.5)                 # Value within boundaries should remain unchanged
  expect_equal(result[4], 1.04)                   # Value greater than 0.99999 + lamb should be set to 1
  expect_equal(result[5], 1)                   # Value greater than 0.99999 + lamb should be set to 1
})

test_that("invert_mat computes correct pseudo-inverse", {
  mat <- matrix(c(1, 2, 3, 4), nrow = 2)
  lamb <- 0.5
  rcond <- 1e-7
  result <- invert_mat(mat, lamb, rcond)
  expect_true(is.matrix(result))
})

test_that("invert_mat handles errors and retries", {
  mat <- matrix(c(0, 0, 0, 0), nrow = 2) 
  lamb <- 0.1
  rcond <- 1e-7
  result <- invert_mat(mat, lamb, rcond)
  expect_true(is.matrix(result))
})

test_that("invert_mat_recursive correctly inverts a valid square matrix", {
  mat <- matrix(c(2, -1, -1, 2), nrow = 2)
  lamb <- 0.5
  rcond <- 0.01
  result <- invert_mat_recursive(mat, lamb, rcond)
  expect_true(is.matrix(result))
  expect_equal(dim(result), dim(mat))
})

test_that("invert_mat_recursive handles non-square matrices appropriately", {
  mat <- matrix(1:6, nrow = 2)
  lamb <- 0.5
  rcond <- 0.01
  expect_silent(invert_mat_recursive(mat, lamb, rcond))
})

test_that("invert_mat_recursive handles errors and performs recursive call correctly", {
  mat <- "not a matrix"
  lamb <- 0.5
  rcond <- 0.01
  expect_error(invert_mat_recursive(mat, lamb, rcond))
})

# Test with Different Tolerance Levels
test_that("invert_mat_eigen behaves differently with varying tolerance levels", {
  mat <- matrix(c(1, 0, 0, 1e-4), nrow = 2)
  tol_high <- 1e-2
  tol_low <- 1e-6
  result_high_tol <- invert_mat_eigen(mat, tol_high)
  result_low_tol <- invert_mat_eigen(mat, tol_low)
  expect_true(!is.logical(all.equal(result_high_tol, result_low_tol)))
})

test_that("invert_mat_eigen handles non-square matrices", {
  mat <- matrix(1:6, nrow = 2)
  expect_error(invert_mat_eigen(mat))
})

test_that("invert_mat_eigen returns the same matrix for an identity matrix", {
    mat <- diag(2)
    expected <- mat
    actual <- invert_mat_eigen(mat)
    expect_equal(actual, expected)
})

test_that("invert_mat_eigen returns a zero matrix for a zero matrix input", {
    mat <- matrix(0, nrow = 2, ncol = 2)
    expected <- mat
    expect_error(invert_mat_eigen(mat),
      "Cannot invert the input matrix because all its eigen values are negative or close to zero")
})

test_that("invert_mat_eigen handles matrices with negative eigenvalues", {
    mat <- matrix(c(-2, 0, 0, -3), nrow = 2)
    expect_silent(invert_mat_eigen(mat))
})

# Block-Diagonal LD data generator for RAISS testing
generate_block_diagonal_test_data <- function(seed = 123, block_structure = "overlapping", n_variants = 30) {
  set.seed(seed)
  
  # Create reference panel with variants
  ref_panel <- data.frame(
    chrom = rep(1, n_variants),
    pos = seq(1, n_variants * 10, 10),
    variant_id = paste0("var", seq_len(n_variants)),
    A1 = sample(c("A", "T", "G", "C"), n_variants, replace = TRUE),
    A2 = sample(c("A", "T", "G", "C"), n_variants, replace = TRUE),
    stringsAsFactors = FALSE
  )
  
  # Create known z-scores for every other variant
  known_indices <- seq(1, n_variants, by = 2)
  known_zscores <- data.frame(
    chrom = rep(1, length(known_indices)),
    pos = ref_panel$pos[known_indices],
    variant_id = ref_panel$variant_id[known_indices],
    A1 = ref_panel$A1[known_indices],
    A2 = ref_panel$A2[known_indices],
    z = rnorm(length(known_indices)),
    stringsAsFactors = FALSE
  )
  
  # Define block boundaries based on requested structure
  if (block_structure == "overlapping") {
    # Create blocks with overlaps at boundaries
    block_boundaries <- list(
      c(1, 11),    # Block 1: variants 1-11
      c(11, 21),   # Block 2: variants 11-21 (overlaps with block 1 at var11)
      c(21, n_variants)  # Block 3: variants 21-30 (overlaps with block 2 at var21)
    )
  } else if (block_structure == "non_overlapping") {
    # Create blocks without overlaps
    block_boundaries <- list(
      c(1, 10),    # Block 1: variants 1-10
      c(11, 20),   # Block 2: variants 11-20
      c(21, n_variants)  # Block 3: variants 21-30
    )
  } else if (block_structure == "uneven") {
    # Create uneven blocks (small, large, medium)
    block_boundaries <- list(
      c(1, 5),                 # Small block
      c(6, 20),                # Large block
      c(21, n_variants)        # Medium block
    )
  } else if (block_structure == "many_small") {
    # Create many small blocks
    block_size <- 5  # 5 variants per block
    n_blocks <- ceiling(n_variants / block_size)
    block_boundaries <- list()
    
    for (i in 1:n_blocks) {
      start_idx <- (i-1) * block_size + 1
      end_idx <- min(i * block_size, n_variants)
      block_boundaries[[i]] <- c(start_idx, end_idx)
    }
  } else if (block_structure == "single_block") {
    # Just one block covering all variants
    block_boundaries <- list(c(1, n_variants))
  }
  
  # Initialize matrix with zeros
  LD_matrix <- matrix(0, nrow = n_variants, ncol = n_variants)
  rownames(LD_matrix) <- colnames(LD_matrix) <- ref_panel$variant_id
  
  # Fill in block-diagonal structure - IMPORTANT: zero correlation between blocks
  for (b in seq_along(block_boundaries)) {
    start_idx <- block_boundaries[[b]][1]
    end_idx <- block_boundaries[[b]][2]
    block_range <- start_idx:end_idx
    
    # Fill this block with correlations
    for (i in block_range) {
      for (j in block_range) {
        if (i == j) {
          LD_matrix[i, j] <- 1  # Diagonal is 1
        } else {
          # Decreasing correlation with distance, but contained within block
          LD_matrix[i, j] <- 0.95^abs(i-j)
        }
      }
    }
  }
  
  # Create variant_indices and block matrices
  variant_indices <- data.frame(
    variant_id = character(),
    block_id = integer(),
    stringsAsFactors = FALSE
  )
  
  block_matrices <- list()
  for (i in seq_along(block_boundaries)) {
    start_idx <- block_boundaries[[i]][1]
    end_idx <- block_boundaries[[i]][2]
    
    # Get variant IDs for this block
    block_variant_ids <- ref_panel$variant_id[start_idx:end_idx]
    
    # Create block matrix with proper row/col names
    block_matrix <- LD_matrix[block_variant_ids, block_variant_ids, drop = FALSE]
    block_matrices[[i]] <- block_matrix
    
    # Add to variant indices
    block_indices <- data.frame(
      variant_id = block_variant_ids,
      block_id = i,
      stringsAsFactors = FALSE
    )
    variant_indices <- rbind(variant_indices, block_indices)
  }
  
  # Create block metadata for partition_LD_matrix
  block_sizes <- sapply(block_boundaries, function(b) b[2] - b[1] + 1)
  block_metadata <- data.frame(
    block_id = seq_along(block_boundaries),
    chrom = rep(1, length(block_boundaries)),
    size = block_sizes,
    start_idx = sapply(block_boundaries, function(b) b[1]),
    end_idx = sapply(block_boundaries, function(b) b[2]),
    stringsAsFactors = FALSE
  )
  
  # Create the block list structure expected by raiss and partition_LD_matrix
  LD_matrix_blocks <- list(
    ld_matrices = block_matrices,
    variant_indices = variant_indices,
    block_metadata = block_metadata,
    combined_LD_variants = ref_panel$variant_id
  )
  
  return(list(
    ref_panel = ref_panel,
    known_zscores = known_zscores,
    LD_matrix_full = LD_matrix,
    LD_matrix_blocks = LD_matrix_blocks,
    variant_indices = variant_indices,
    block_boundaries = block_boundaries,
    block_metadata = block_metadata
  ))
}

# Verify that matrix has block-diagonal structure
verify_block_diagonal <- function(matrix, block_boundaries) {
  for (i in 1:(length(block_boundaries)-1)) {
    for (j in (i+1):length(block_boundaries)) {
      block_i_range <- block_boundaries[[i]][1]:block_boundaries[[i]][2]
      block_j_range <- block_boundaries[[j]][1]:block_boundaries[[j]][2]
      
      # Get cross-block submatrix (excluding shared boundary variants)
      non_overlapping_i <- setdiff(block_i_range, block_j_range)
      non_overlapping_j <- setdiff(block_j_range, block_i_range)
      
      if (length(non_overlapping_i) > 0 && length(non_overlapping_j) > 0) {
        cross_block <- matrix[non_overlapping_i, non_overlapping_j]
        
        # Check if any elements are non-zero
        if (any(abs(cross_block) > 1e-10)) {
          return(FALSE)
        }
      }
    }
  }
  return(TRUE)
}

# Test equivalence between full matrix and block processing
test_that("full matrix and block processing produce identical results", {
  # Test with different block structures
  block_structures <- c("overlapping", "non_overlapping", "uneven", "many_small", "single_block")
  
  for (structure in block_structures) {
    # Generate block-diagonal test data
    test_data <- generate_block_diagonal_test_data(seed = 123, block_structure = structure)
    
    # Verify the generated matrix has proper block-diagonal structure
    is_block_diagonal <- verify_block_diagonal(test_data$LD_matrix_full, test_data$block_boundaries)
    expect_true(is_block_diagonal, info = paste("Test data for", structure, "should have block-diagonal structure"))
    
    # Run with full matrix
    result_full <- raiss(
      ref_panel = test_data$ref_panel,
      known_zscores = test_data$known_zscores,
      LD_matrix = test_data$LD_matrix_full,
      lamb = 0.01,
      rcond = 0.01,
      R2_threshold = 0.3,
      minimum_ld = 1,
      verbose = FALSE
    )
    
    # Run with blocks
    result_blocks <- raiss(
      ref_panel = test_data$ref_panel,
      known_zscores = test_data$known_zscores,
      LD_matrix = test_data$LD_matrix_blocks,
      variant_indices = test_data$variant_indices,
      lamb = 0.01,
      rcond = 0.01,
      R2_threshold = 0.3,
      minimum_ld = 1,
      verbose = FALSE
    )
    
    # Check that both methods produce results for the same variants
    expect_equal(
      sort(result_full$result_nofilter$variant_id),
      sort(result_blocks$result_nofilter$variant_id),
      info = paste("Variant IDs should match for", structure, "blocks")
    )
    
    # Sort both results by variant_id for comparison
    result_full$result_nofilter <- result_full$result_nofilter %>% arrange(variant_id)
    result_blocks$result_nofilter <- result_blocks$result_nofilter %>% arrange(variant_id)
    
    # Check that Z-scores match with appropriate tolerance
    expect_equal(
      result_full$result_nofilter$z,
      result_blocks$result_nofilter$z,
      tolerance = 1e-4,
      info = paste("Z-scores should match for", structure, "blocks")
    )
    
    # Check filtered results if they exist
    if (!is.null(result_full$result_filter) && !is.null(result_blocks$result_filter) &&
        nrow(result_full$result_filter) > 0 && nrow(result_blocks$result_filter) > 0) {
      # Check that both methods filter the same variants
      expect_equal(
        sort(result_full$result_filter$variant_id),
        sort(result_blocks$result_filter$variant_id),
        info = paste("Filtered variant IDs should match for", structure, "blocks")
      )
      
      # Sort filtered results
      result_full$result_filter <- result_full$result_filter %>% arrange(variant_id)
      result_blocks$result_filter <- result_blocks$result_filter %>% arrange(variant_id)
      
      # Check that filtered Z-scores match
      expect_equal(
        result_full$result_filter$z,
        result_blocks$result_filter$z,
        tolerance = 1e-4,
        info = paste("Filtered Z-scores should match for", structure, "blocks")
      )
    }
  }
})

# Test partition_LD_matrix integration with proper block-diagonal structure
test_that("partition_LD_matrix and raiss integration works correctly", {
  # Generate block-diagonal test data
  test_data <- generate_block_diagonal_test_data(seed = 456, block_structure = "non_overlapping", n_variants = 30)
  
  # Get the full matrix
  full_matrix <- test_data$LD_matrix_full
  
  # Create LD data for partition_LD_matrix
  ld_data <- list(
    combined_LD_matrix = full_matrix,
    combined_LD_variants = rownames(full_matrix),
    block_metadata = test_data$block_metadata
  )
  
  # Partition the matrix
  partitioned <- partition_LD_matrix(
    ld_data, 
    merge_small_blocks = TRUE,
    min_merged_block_size = 10,
    max_merged_block_size = 30
  )
  
  # Run raiss with the partitioned data
  result_partitioned <- raiss(
    ref_panel = test_data$ref_panel,
    known_zscores = test_data$known_zscores,
    LD_matrix = partitioned,
    lamb = 0.01,
    rcond = 0.01,
    R2_threshold = 0.3,
    minimum_ld = 1,
    verbose = FALSE
  )
  
  # Run raiss with the full matrix for comparison
  result_full <- raiss(
    ref_panel = test_data$ref_panel,
    known_zscores = test_data$known_zscores,
    LD_matrix = full_matrix,
    lamb = 0.01,
    rcond = 0.01,
    R2_threshold = 0.3,
    minimum_ld = 1,
    verbose = FALSE
  )
  
  # Verify results match
  result_full$result_nofilter <- result_full$result_nofilter %>% arrange(variant_id)
  result_partitioned$result_nofilter <- result_partitioned$result_nofilter %>% arrange(variant_id)
  
  expect_equal(
    result_full$result_nofilter$variant_id,
    result_partitioned$result_nofilter$variant_id,
    info = "Variant IDs should match between full and partitioned results"
  )
  
  expect_equal(
    result_full$result_nofilter$z,
    result_partitioned$result_nofilter$z,
    tolerance = 1e-4,
    info = "Z-scores should match between full and partitioned results"
  )
})

# Test boundary overlap handling with block-diagonal matrices
test_that("boundary overlaps are correctly handled in block-diagonal matrices", {
  # Generate test data with overlapping blocks
  test_data <- generate_block_diagonal_test_data(seed = 789, block_structure = "overlapping")
  
  # Find boundary variants (those that appear in multiple blocks)
  variant_counts <- table(test_data$variant_indices$variant_id)
  boundary_vars <- names(variant_counts[variant_counts > 1])
  
  # Run raiss with blocks
  result_blocks <- raiss(
    ref_panel = test_data$ref_panel,
    known_zscores = test_data$known_zscores,
    LD_matrix = test_data$LD_matrix_blocks,
    variant_indices = test_data$variant_indices,
    lamb = 0.01,
    rcond = 0.01,
    R2_threshold = 0.1,
    minimum_ld = 1,
    verbose = FALSE
  )
  
  # Verify each boundary variant appears only once in the results
  for (var in boundary_vars) {
    expect_equal(
      sum(result_blocks$result_nofilter$variant_id == var),
      1,
      info = paste("Boundary variant", var, "should appear exactly once in results")
    )
  }
  
  # Verify no duplicates in the results
  expect_equal(
    nrow(result_blocks$result_nofilter),
    length(unique(result_blocks$result_nofilter$variant_id)),
    info = "Results should have no duplicate variants"
  )
})

# Test single-block list case with block-diagonal matrix
test_that("raiss handles list with single block correctly", {
  # Generate test data with a single block
  test_data <- generate_block_diagonal_test_data(seed = 202, block_structure = "single_block")
  
  # Run with full matrix
  result_full <- raiss(
    ref_panel = test_data$ref_panel,
    known_zscores = test_data$known_zscores,
    LD_matrix = test_data$LD_matrix_full,
    lamb = 0.01,
    rcond = 0.01,
    R2_threshold = 0.3,
    minimum_ld = 1,
    verbose = FALSE
  )
  
  # Run with single-block list (testing your single-block handling)
  result_single_block <- raiss(
    ref_panel = test_data$ref_panel,
    known_zscores = test_data$known_zscores,
    LD_matrix = test_data$LD_matrix_blocks,
    variant_indices = test_data$variant_indices,
    lamb = 0.01,
    rcond = 0.01,
    R2_threshold = 0.3,
    minimum_ld = 1,
    verbose = FALSE
  )
  
  # Sort results
  result_full$result_nofilter <- result_full$result_nofilter %>% arrange(variant_id)
  result_single_block$result_nofilter <- result_single_block$result_nofilter %>% arrange(variant_id)
  
  # Check z-scores match exactly
  expect_equal(
    result_full$result_nofilter$z,
    result_single_block$result_nofilter$z,
    tolerance = 1e-6,  # Stricter tolerance for this case
    info = "Z-scores should match exactly for single-block case"
  )
})