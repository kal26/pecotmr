suppressPackageStartupMessages({
  library(qvalue)
  library(dplyr)
  library(tidyr)
  library(readr)
})

# Utility function for NULL coalescing
`%||%` <- function(a, b) if (is.null(a)) b else a

# Standardize chromosome format - VECTORIZED VERSION
standardize_chrom <- function(chrom_vector) {
  chrom_vector <- as.character(chrom_vector)
  # Vectorized operation without if statement
  needs_prefix <- !grepl("^chr", chrom_vector)
  chrom_vector[needs_prefix] <- paste0("chr", chrom_vector[needs_prefix])
  return(chrom_vector)
}

safe_qvalue <- function(p, ...) {
  # First try default settings
  tryCatch(
    {
      return(qvalue(p, ...))
    },
    error = function(e) {
      # For "missing or infinite values" error
      if (grepl("missing or infinite", e$message)) {
        return(qvalue(p, lambda = 0, ...))
      }
      # For "pi0 <= 0" error
      else if (grepl("pi0 <= 0", e$message)) {
        max_p <- max(p)
        lambda_seq <- seq(0, min(0.9, max_p * 0.95), length.out = 10)
        return(qvalue(p, lambda = lambda_seq, pi0.method = "bootstrap", ...))
      }
      # If all else fails
      else {
        q <- p.adjust(p, method = "BH")
        result <- list(call = match.call(), pi0 = 1, qvalues = q, pvalues = p)
        class(result) <- "qvalue"
        return(result)
      }
    }
  )
}

find_common_prefix <- function(files) {
  if (length(files) == 0) {
    return("")
  }
  if (length(files) == 1) {
    return(files)
  }

  # Extract basenames if full paths were provided
  filenames <- basename(files)

  # Split by dots
  parts_list <- strsplit(filenames, "\\.")

  # Find minimum number of parts
  min_parts <- min(sapply(parts_list, length))

  # Compare parts position by position
  common_parts <- character(0)
  for (i in 1:min_parts) {
    current_parts <- sapply(parts_list, `[`, i)
    if (length(unique(current_parts)) == 1) {
      common_parts <- c(common_parts, current_parts[1])
    } else {
      break
    }
  }

  # Return the joined common parts
  if (length(common_parts) == 0) {
    return("")
  }
  return(paste(common_parts, collapse = "."))
}

read_and_combine_files <- function(files, ...) {
  if (length(files) == 0) stop("No files provided")

  # Check file existence
  missing <- files[!file.exists(files)]
  if (length(missing) > 0) stop("Missing files: ", paste(missing, collapse = ", "))

  # Determine if files are compressed
  is_compressed <- grepl("\\.(gz|bz2|xz)$", files)

  # Load first file
  if (all(is_compressed)) {
    data <- data.table::fread(cmd = paste("zcat", shQuote(files[1])), ...)
  } else {
    data <- data.table::fread(files[1], ...)
  }

  # Load remaining files and append
  if (length(files) > 1) {
    for (i in 2:length(files)) {
      if (is_compressed[i]) {
        next_data <- data.table::fread(
          cmd = paste("zcat", shQuote(files[i])),
          skip = 1, header = FALSE, ...
        )
      } else {
        next_data <- data.table::fread(files[i], skip = 1, header = FALSE, ...)
      }

      # Set column names to match first file
      data.table::setnames(next_data, names(data))

      # Append using rbindlist (faster than rbind)
      data <- data.table::rbindlist(list(data, next_data), use.names = TRUE, fill = TRUE)
    }
  }

  return(data)
}

# ===================================
# Data Loaders
# ===================================
load_regional_data <- function(params) {
  setwd(params$workdir)

  # Load permutation-based regional files
  regional_files <- vector()
  if (!is.null(params$regional_pattern)) {
    regional_files <- list.files(pattern = params$regional_pattern, full.names = TRUE)
  }
  data <- NULL
  has_permutation <- FALSE
  n_var_col <- NULL

  if (length(regional_files) > 0) {
    data <- read_and_combine_files(regional_files) %>%
      {
        if ("chrom" %in% colnames(.)) mutate(., chrom = standardize_chrom(chrom)) else .
      }
    has_emt_n_var <- any(grepl("tests_emt", colnames(data)))
    # eigenMT test as number of effective variants
    if (has_emt_n_var) {
      data <- data %>%
        rename_with(~"molecular_trait_object_id", matches("phenotype_id"))
      n_var_col <- "tests_emt"
      has_permutation <- FALSE
      message("Found 'tests_emt' column in regional data, converting to n_variants")
    } else {
      n_var_col <- "n_variants"
      has_permutation <- TRUE
    }
  }

  return(list(
    regional_summary = data,
    regional_data_files = regional_files,
    has_permutation = has_permutation,
    n_var_col = n_var_col
  ))
}

extract_column_names <- function(file_path, pvalue_pattern = "pvalue", qvalue_pattern = "qvalue") {
  header <- system(paste0("zcat ", file_path, " | head -1"), intern = TRUE)
  cols <- strsplit(header, "\t")[[1]]
  column_info <- list(
    all_columns = cols
  )

  # Determine p-value column using pattern
  p_cols <- grep(pvalue_pattern, cols, value = TRUE)
  column_info$p_col <- if (length(p_cols) > 0) p_cols[1] else "pvalue"

  # Find q-value column using pattern
  q_cols <- grep(qvalue_pattern, cols, value = TRUE)
  column_info$q_col <- if (length(q_cols) > 0) q_cols[1] else "qvalue"

  column_info$p_idx <- which(cols == column_info$p_col)
  if (length(column_info$p_idx) == 0) {
    stop(sprintf("P-value column '%s' not found", column_info$p_col))
  }

  return(column_info)
}

# Gene coordinate data loader
load_gene_coordinates <- function(params) {
  gene_coords <- data.table::fread(params$gene_coordinates) %>% mutate(chr = standardize_chrom(chr))
  message(sprintf("Loaded gene coordinates with %d entries", nrow(gene_coords)))
  return(gene_coords)
}

# Calculate feature positions for cis-window filtering
calculate_feature_positions <- function(qtl_data, cis_window, gene_coords, start_distance_col = "start_distance", end_distance_col = "end_distance") {
  # Extract ENSEMBL IDs from molecular_trait_object_id
  extract_ensembl <- function(ids) {
    pattern <- "^.*?(ENSG\\d+).*$"
    ensembl_ids <- gsub(pattern, "\\1", ids)
    not_matched <- ensembl_ids == ids & !grepl("ENSG\\d+", ids)
    if (any(not_matched)) {
      ensembl_ids[not_matched] <- NA
    }
    return(ensembl_ids)
  }

  unique_traits <- qtl_data %>%
    select(molecular_trait_object_id, chrom) %>%
    distinct() %>%
    mutate(
      ensembl_id = extract_ensembl(molecular_trait_object_id)
    )

  # Join with lookup to get TSS/TES positions
  merged_traits <- unique_traits %>%
    left_join(
      gene_coords %>%
        select(chrom = chr, gene_id, gene_start = start, gene_end = end),
      by = c("ensembl_id" = "gene_id", "chrom" = "chrom")
    )

  # For unmapped genes, calculate approximate positions from variant data
  if (sum(is.na(merged_traits$gene_start)) > 0) {
    fallback_positions <- qtl_data %>%
      filter(molecular_trait_object_id %in%
        merged_traits$molecular_trait_object_id[is.na(merged_traits$gene_start)]) %>%
      group_by(molecular_trait_object_id, chrom) %>%
      summarize(
        approx_tss = first(pos) - first(!!sym(start_distance_col)),
        approx_tes = first(pos) - first(!!sym(end_distance_col)),
        .groups = "drop"
      )

    merged_traits <- merged_traits %>%
      left_join(fallback_positions, by = c("molecular_trait_object_id", "chrom")) %>%
      mutate(
        feature_tss = coalesce(gene_start, approx_tss),
        feature_tes = coalesce(gene_end, approx_tes)
      )
  } else {
    merged_traits <- merged_traits %>%
      mutate(
        feature_tss = gene_start,
        feature_tes = gene_end
      )
  }

  # Add cis window range
  feature_positions <- merged_traits %>%
    mutate(
      cis_start = feature_tss - cis_window,
      cis_end = feature_tes + cis_window
    ) %>%
    select(molecular_trait_object_id, chrom, feature_tss, feature_tes, cis_start, cis_end)

  return(feature_positions)
}

calculate_filtered_variant_counts <- function(filename, params, gene_coords) {
  message(sprintf("Counting per event variants in %s", basename(filename)))
  all_cols <- extract_column_names(filename, params$pvalue_pattern, params$qvalue_pattern)$all_columns
  required_cols <- c("molecular_trait_object_id", "chrom", "pos", params$af_col)

  col_indices <- sapply(required_cols, function(col) which(all_cols == col))
  if (any(sapply(col_indices, length) == 0)) {
    stop(sprintf("Required columns missing in file: %s", filename))
  }
  awk_cols <- paste(paste0("$", col_indices), collapse = ", ")
  cmd <- sprintf(
    "zcat %s | awk 'NR==1 {print \"%s\"} NR>1 {print %s}' OFS='\\t'",
    shQuote(filename), paste(required_cols, collapse = "\t"), awk_cols
  )
  message("Extracting minimal columns from file...")
  qtl_data <- data.table::fread(cmd = cmd) %>% mutate(chrom = standardize_chrom(chrom))
  trait_chrom <- qtl_data %>%
    group_by(molecular_trait_object_id) %>%
    summarize(chrom = first(chrom))
  original_counts <- qtl_data %>%
    group_by(molecular_trait_object_id) %>%
    summarize(n_variants_original = n())
  message("Calculating feature positions and cis windows...")
  feature_positions <- calculate_feature_positions(
    qtl_data,
    params$cis_window,
    gene_coords
  )

  message("Applying MAF and cis-window filters...")
  filtered_data <- qtl_data %>%
    left_join(
      feature_positions %>% select(molecular_trait_object_id, cis_start, cis_end),
      by = "molecular_trait_object_id"
    ) %>%
    filter(
      pmin(!!sym(params$af_col), 1 - !!sym(params$af_col)) > params$maf_cutoff,
      pos >= cis_start & pos <= cis_end
    )

  # Count filtered variants
  filtered_counts <- filtered_data %>%
    group_by(molecular_trait_object_id) %>%
    summarize(n_variants_filtered = n())

  # Combine results
  results <- trait_chrom %>%
    left_join(original_counts, by = "molecular_trait_object_id") %>%
    left_join(filtered_counts, by = "molecular_trait_object_id") %>%
    mutate(n_variants_filtered = ifelse(is.na(n_variants_filtered), 0, n_variants_filtered)) %>%
    rename(n_variants = n_variants_original) %>%
    select(chrom, molecular_trait_object_id, n_variants, n_variants_filtered)

  message(sprintf("Processed %d traits in %s", nrow(results), basename(filename)))
  return(results)
}

load_n_variants_data <- function(params, gene_coords) {
  setwd(params$workdir)
  if (is.null(params$qtl_pattern)) {
    stop("params$qtl_pattern must be provided")
  }

  if (is.null(params$n_variants_suffix)) {
    stop("params$n_variants_suffix must be provided")
  }
  qtl_files <- list.files(pattern = params$qtl_pattern, full.names = TRUE)
  # Create corresponding n_variants filenames by replacing the pattern
  n_variants_files <- vector("character", length(qtl_files))
  for (i in seq_along(qtl_files)) {
    # Extract base name by removing the pattern (removing $ from pattern)
    cleaned_pattern <- sub("\\$$", "", params$qtl_pattern) # Remove $ from end if present
    cleaned_pattern <- sub("\\*", "", cleaned_pattern) # Remove * from beginning
    base_name <- sub(cleaned_pattern, "", qtl_files[i])
    n_variants_suffix <- sub("\\*", "", params$n_variants_suffix)
    n_variants_suffix <- sub("\\$$", "", n_variants_suffix)
    n_variants_suffix <- sprintf(
      "maf_%s_window_%s_%s",
      params$maf_cutoff,
      format(params$cis_window, scientific = FALSE),
      n_variants_suffix
    )
    n_variants_files[i] <- paste0(base_name, ".", n_variants_suffix)
  }

  # Check if each n_variants file exists, if not, calculate it
  for (i in seq_along(qtl_files)) {
    qtl_file <- qtl_files[i]
    n_variants_file <- n_variants_files[i]

    if (!file.exists(n_variants_file)) {
      message(sprintf("Calculating n_variants for %s", qtl_file))
      n_variants <- calculate_filtered_variant_counts(qtl_file, params, gene_coords)
      write.table(n_variants,
        file = gzfile(n_variants_file),
        quote = FALSE, sep = "\t", row.names = FALSE
      )
    }
  }

  n_variants_data <- NULL
  if (length(n_variants_files) > 0) {
    message("Loading n_variants statistics")
    n_variants_data <- read_and_combine_files(n_variants_files) %>%
      mutate(chrom = standardize_chrom(chrom))
  }

  return(n_variants_data)
}

load_qtl_data <- function(params, load_n_variants = FALSE) {
  setwd(params$workdir)
  gene_coords <- load_gene_coordinates(params)

  files <- list.files(pattern = params$qtl_pattern, full.names = TRUE)
  if (length(files) == 0) {
    stop("No pair files found")
  }

  # Get column names
  column_info <- extract_column_names(files[1], params$pvalue_pattern, params$qvalue_pattern)

  # Only filter if pvalue_cutoff is less than 1
  filter_by_p <- params$pvalue_cutoff < 1

  if (filter_by_p) {
    # Filter rows with p-value < threshold for efficiency
    files_str <- paste(shQuote(files), collapse = " ")
    awk_cmd <- sprintf(
      "awk 'NR==1 {print; next} $%d < %s'",
      column_info$p_idx, params$pvalue_cutoff
    )
    cmd <- sprintf("zcat %s | %s", files_str, awk_cmd)

    message(sprintf("Only loading QTL data with p-value < %g", params$pvalue_cutoff))
    data <- data.table::fread(cmd = cmd) %>% mutate(chrom = standardize_chrom(chrom))
  } else {
    # Load all data (may be huge)
    message("Loading all pair data (no p-value filtering)")
    data <- read_and_combine_files(files) %>% mutate(chrom = standardize_chrom(chrom))
  }

  # Load n_variants data if requested
  n_variants_data <- NULL
  if (load_n_variants) {
    n_variants_data <- load_n_variants_data(params, gene_coords)
  }

  return(list(
    data = data,
    file_prefix = find_common_prefix(files),
    gene_coords = gene_coords,
    qtl_files = files,
    column_info = column_info,
    n_variants_data = n_variants_data
  ))
}

annotate_qtl_with_regional <- function(qtl_data, regional_data, n_variants_data = NULL, use_filtered = FALSE) {
  if (nrow(qtl_data) == 0) {
    # Return empty dataframe with n_variants column
    return(qtl_data %>% mutate(n_variants = integer(0)))
  }

  # First priority: If use_filtered is TRUE and n_variants_filtered exists in n_variants_data
  if (use_filtered && !is.null(n_variants_data) && "n_variants_filtered" %in% names(n_variants_data)) {
    n_variants_info <- n_variants_data %>%
      select(molecular_trait_object_id, n_variants = n_variants_filtered) %>%
      distinct()
    # Second priority: Use regional_data$regional_summary if available
  } else if (!is.null(regional_data$regional_summary)) {
    n_variants_info <- regional_data$regional_summary %>%
      select(molecular_trait_object_id, n_variants = !!sym(regional_data$n_var_col)) %>%
      distinct()
    # Third priority: Fallback to n_variants in n_variants_data
  } else if (!is.null(n_variants_data)) {
    n_variants_info <- n_variants_data %>%
      select(molecular_trait_object_id, n_variants = n_variants) %>%
      distinct()
    # No n_variants information available
  } else {
    stop("No n_variants data found. Using counts from data which can be biased if only loading partial QTL data eg with nominal pvalue_cutoff < 1")
  }

  # Join to add n_variants
  if ("n_variants" %in% colnames(qtl_data)) {
    qtl_data <- qtl_data %>%
      select(-n_variants)
  }
  annotated_data <- qtl_data %>%
    left_join(n_variants_info, by = "molecular_trait_object_id")

  return(annotated_data)
}

prepare_local_qtl_data <- function(qtl_data, regional_data, params, should_filter = TRUE) {
  original_data <- NULL
  filtered_data <- NULL

  if (should_filter) {
    feature_positions <- calculate_feature_positions(qtl_data$data, params$cis_window, qtl_data$gene_coords)
    filtered_data <- qtl_data$data %>%
      left_join(
        feature_positions %>% select(molecular_trait_object_id, cis_start, cis_end),
        by = "molecular_trait_object_id"
      ) %>%
      filter(pmin(!!sym(params$af_col), 1 - !!sym(params$af_col)) > params$maf_cutoff) %>%
      filter(pos >= cis_start & pos <= cis_end)

    # Check if filtered data has rows
    if (nrow(filtered_data) == 0) {
      message("Warning: Filtered data has 0 rows after filtering")
      # Return empty dataframe with appropriate columns
      return(list(
        original_data = original_data,
        filtered_data = filtered_data %>% mutate(n_variants = integer(0))
      ))
    }
    filtered_data <- annotate_qtl_with_regional(filtered_data, regional_data, n_variants_data = qtl_data$n_variants_data, use_filtered = TRUE)
  } else {
    # Annotate data with n_variants
    original_data <- annotate_qtl_with_regional(qtl_data$data, regional_data, n_variants_data = qtl_data$n_variants_data, use_filtered = FALSE)
  }

  if (!is.null(original_data)) {
    avg_variants <- original_data %>%
      group_by(molecular_trait_object_id) %>%
      slice(1) %>%
      ungroup() %>%
      summarize(avg = mean(n_variants)) %>%
      pull(avg)
    message(
      "Original data: ", nrow(original_data), " rows (avg n_variants per event: ",
      round(avg_variants, 0), ")"
    )
  }

  if (!is.null(filtered_data)) {
    avg_variants <- filtered_data %>%
      group_by(molecular_trait_object_id) %>%
      slice(1) %>%
      ungroup() %>%
      summarize(avg = mean(n_variants)) %>%
      pull(avg)
    message(
      "Filtered data: ", nrow(filtered_data), " rows (avg n_variants per event: ",
      round(avg_variants, 0), ")"
    )
  }

  return(list(
    original_data = original_data,
    filtered_data = filtered_data
  ))
}

############################################
# Local Adjustment (Step 1)
############################################

permutation_local_adjustment <- function(data, params) {
  message("Loading permutation-based local adjustment")
  reg_data <- data$regional_data$regional_summary
  if (!("p_perm" %in% colnames(reg_data)) && !("p_beta" %in% colnames(reg_data))) {
    stop("Permutation local adjustment requires p_perm or p_beta columns in regional data")
  }

  event_level_pvalues <- reg_data %>%
    select(molecular_trait_object_id) %>%
    distinct()
  for (col in c("p_perm", "p_beta")) {
    if (col %in% colnames(reg_data)) {
      event_level_pvalues[[col]] <- reg_data[[col]]
    }
  }

  data$event_level_pvalues <- event_level_pvalues
  data$local_adjustment_info <- list(
    method = "permutation",
    p_value_columns = c(
      if ("p_perm" %in% colnames(event_level_pvalues)) "p_perm",
      if ("p_beta" %in% colnames(event_level_pvalues)) "p_beta"
    ),
    is_filtered = FALSE
  )
  return(data)
}

bonferroni_local_adjustment <- function(data, params, should_filter = FALSE) {
  should_filter <- (params$maf_cutoff > 0 || params$cis_window > 0) && should_filter
  message(sprintf(
    "Applying Bonferroni local adjustment (filter applied: %s)...",
    ifelse(should_filter, "Yes", "No")
  ))

  p_col <- data$qtl_data$column_info$p_col
  qtl_data_updated <- prepare_local_qtl_data(
    data$qtl_data,
    data$regional_data,
    params,
    should_filter
  )

  if (should_filter) {
    if (is.null(qtl_data_updated$filtered_data) || nrow(qtl_data_updated$filtered_data) == 0) {
      message("No data remains after filtering. Returning empty result.")
      # Create empty event level p-values
      event_level_pvalues <- data.frame(
        molecular_trait_object_id = character(0),
        p_bonferroni_min = numeric(0)
      )

      result <- data
      result$qtl_data$data <- data.frame()
      result$event_level_pvalues <- event_level_pvalues
      result$local_adjustment_info <- list(
        method = "bonferroni",
        p_value_columns = "p_bonferroni_min",
        is_filtered = should_filter
      )
      return(result)
    }

    adjusted_qtl_data <- qtl_data_updated$filtered_data %>%
      mutate(p_bonferroni_adj = pmin(1, !!sym(p_col) * n_variants))
  } else {
    adjusted_qtl_data <- qtl_data_updated$original_data %>%
      mutate(p_bonferroni_adj = pmin(1, !!sym(p_col) * n_variants))
  }

  # Calculate gene-level p-values (min p-value per gene)
  event_level_pvalues <- adjusted_qtl_data %>%
    group_by(molecular_trait_object_id) %>%
    summarize(p_bonferroni_min = min(p_bonferroni_adj)) %>%
    ungroup()

  result <- data
  result$qtl_data$data <- adjusted_qtl_data
  result$event_level_pvalues <- event_level_pvalues
  result$local_adjustment_info <- list(
    method = "bonferroni",
    p_value_columns = "p_bonferroni_min",
    is_filtered = should_filter
  )
  return(result)
}

############################################
# Global Adjustment (Step 2)
############################################
perform_global_adjustment <- function(data, params) {
  message("Applying both FDR and qvalue global adjustments...")

  if (is.null(data$event_level_pvalues) || is.null(data$local_adjustment_info)) {
    stop("Cannot perform global adjustment: no event-level p-values found")
  }

  # Handle empty event_level_pvalues
  if (nrow(data$event_level_pvalues) == 0) {
    message("No events to adjust. Returning empty result.")

    result <- data
    if (is.null(result$regional_data)) {
      result$regional_data <- list()
    }

    result$regional_data$regional_summary <- data.frame()
    result$global_adjustment_info <- list(
      event_adjusted = data.frame(),
      local_method = result$local_adjustment_info$method,
      p_value_columns = character(0),
      fdr_columns = character(0),
      q_columns = character(0),
      statistics = list()
    )

    return(result)
  }

  p_columns <- data$local_adjustment_info$p_value_columns
  event_adjusted <- data.frame(molecular_trait_object_id = data$event_level_pvalues$molecular_trait_object_id)
  fdr_columns <- c()
  q_columns <- c()

  statistics <- list()

  for (p_col in p_columns) {
    # First add the original p-value column
    event_adjusted[[p_col]] <- data$event_level_pvalues[[p_col]]

    # Apply FDR adjustment (Benjamini-Hochberg)
    fdr_col <- sub("^p_", "fdr_", p_col)
    event_adjusted[[fdr_col]] <- p.adjust(data$event_level_pvalues[[p_col]], method = "fdr")
    fdr_columns <- c(fdr_columns, fdr_col)

    # Apply qvalue adjustment (Storey's q-value)
    q_col <- sub("^p_", "q_", p_col)
    event_adjusted[[q_col]] <- safe_qvalue(data$event_level_pvalues[[p_col]])$qvalues
    q_columns <- c(q_columns, q_col)

    # Collect statistics for summary
    for (threshold in c(0.05, 0.01)) {
      for (adj_col in c(fdr_col, q_col)) {
        stat_name <- sprintf("%s_%.2f", adj_col, threshold)
        statistics[[stat_name]] <- sum(event_adjusted[[adj_col]] < threshold, na.rm = TRUE)
      }
    }
  }

  result <- data # Copy input data to modify and carry over
  all_adjustment_cols <- c(fdr_columns, q_columns)
  if (is.null(result$regional_data)) {
    result$regional_data <- list()
  }

  if (!is.null(result$regional_data$regional_summary)) {
    # If regional_summary exists, update it with the adjustments
    result$regional_data$regional_summary <- result$regional_data$regional_summary %>%
      left_join(event_adjusted, by = "molecular_trait_object_id", suffix = c("", ".new")) %>%
      mutate(across(matches("\\.new$"),
        ~ coalesce(., get(sub("\\.new$", "", cur_column()))),
        .names = "{.col}"
      )) %>%
      select(-ends_with(".new"))
  } else {
    # Create a new regional_summary from the event_adjusted dataframe
    result$regional_data$regional_summary <- event_adjusted
  }

  result$global_adjustment_info <- list(
    local_method = result$local_adjustment_info$method,
    p_value_columns = p_columns,
    fdr_columns = fdr_columns,
    q_columns = q_columns,
    statistics = statistics
  )

  return(result)
}

############################################
# SNP Identification (Step 3)
############################################

identify_permutation_snps <- function(data, params) {
  message("Identifying significant SNPs using permutation thresholds...")

  # Check if permutation data (with beta parameters) is available
  if (is.null(data$regional_data) || is.null(data$regional_data$regional_summary) ||
    !all(c("p_beta", "q_beta", "beta_shape1", "beta_shape2") %in% colnames(data$regional_data$regional_summary))) {
    warning("Permutation data (p_beta, q_beta, beta shape parameters) not available")
    # Return data with empty significant QTLs
    if (is.null(data$significant_qtls)) data$significant_qtls <- list()
    data$significant_qtls[["permutation_significant"]] <- data.frame()
    return(data)
  }

  # Calculate threshold using p_beta and q_beta specifically
  regional_data <- data$regional_data$regional_summary

  # Find significant and non-significant p_beta values
  lb <- regional_data %>%
    filter(q_beta <= params$fdr_threshold) %>%
    pull(p_beta) %>%
    sort()

  ub <- regional_data %>%
    filter(q_beta > params$fdr_threshold) %>%
    pull(p_beta) %>%
    sort()

  if (length(lb) > 0) {
    lb_val <- tail(lb, 1) # Max p_beta that passes threshold
    threshold <- if (length(ub) > 0) {
      (lb_val + head(ub, 1)) / 2 # Average with min p_beta that fails
    } else {
      lb_val # Just use max passing if no failures
    }

    message(sprintf(
      "Min p-value threshold @ q_beta %.2f: %g",
      params$fdr_threshold, threshold
    ))

    # Update regional summary with permutation threshold
    updated_regional_data <- regional_data %>%
      mutate(p_nominal_threshold = qbeta(threshold, beta_shape1, beta_shape2))

    data$regional_data$regional_summary <- updated_regional_data

    # Identify significant SNPs using permutation threshold
    significant_pairs <- data$qtl_data$data %>%
      left_join(
        updated_regional_data %>% select(molecular_trait_object_id, p_nominal_threshold),
        by = "molecular_trait_object_id"
      ) %>%
      filter(!!sym(data$qtl_data$column_info$p_col) < p_nominal_threshold)

    # Store significant SNPs
    if (is.null(data$significant_qtls)) data$significant_qtls <- list()
    data$significant_qtls[["permutation_significant"]] <- significant_pairs
    significant_events <- data$regional_data$regional_summary %>%
      filter(q_beta < params$fdr_threshold)
    message(sprintf(
      "Identified %d significant SNPs from %d events using permutation thresholds",
      nrow(significant_pairs), nrow(significant_events)
    ))
  } else {
    message("No significant events found using q_beta threshold")
    if (is.null(data$significant_qtls)) data$significant_qtls <- list()
    data$significant_qtls[["permutation_significant"]] <- data.frame()
  }

  return(data)
}

identify_bonferroni_fdr_snps <- function(data, params) {
  message("Identifying significant SNPs using Bonferroni adjusted p-value thresholds...")

  # Get FDR column for Bonferroni
  fdr_bonferroni_col <- "fdr_bonferroni_min"
  p_bonferroni_col <- "p_bonferroni_min"

  # Check if columns exist in regional_summary
  if (!all(c(fdr_bonferroni_col, p_bonferroni_col) %in%
    names(data$regional_data$regional_summary))) {
    warning("Required p-value and FDR columns not found in regional summary")
    if (is.null(data$significant_qtls)) data$significant_qtls <- list()
    data$significant_qtls[["bonferroni_fdr_original"]] <- data.frame()
    return(data)
  }

  # Identify significant events using FDR threshold
  significant_events <- data$regional_data$regional_summary %>%
    filter(!!sym(fdr_bonferroni_col) < params$fdr_threshold)

  if (nrow(significant_events) == 0) {
    message(sprintf(
      "No significant events identified at %s threshold %g",
      fdr_bonferroni_col, params$fdr_threshold
    ))
    if (is.null(data$significant_qtls)) data$significant_qtls <- list()
    data$significant_qtls[["bonferroni_fdr_original"]] <- data.frame()
    return(data)
  }

  # Calculate threshold as maximum p-value among significant events
  threshold <- max(significant_events[[p_bonferroni_col]], na.rm = TRUE)

  # Identify significant SNPs
  significant_pairs <- data$qtl_data$data %>%
    filter(p_bonferroni_adj <= threshold)

  # Create result name
  result_name <- sprintf(
    "bonferroni_fdr_%s",
    if (data$local_adjustment_info$is_filtered) "filtered" else "original"
  )

  # Store significant SNPs
  if (is.null(data$significant_qtls)) data$significant_qtls <- list()
  data$significant_qtls[[result_name]] <- significant_pairs

  message(sprintf(
    "Identified %d significant SNPs from %d events using Bonferroni adjusted p-value threshold %g",
    nrow(significant_pairs), nrow(significant_events), threshold
  ))
  return(data)
}

identify_qvalue_snps <- function(data, params, base_data = NULL) {
  message("Identifying significant SNPs using q-value per event method...")
  if (is.null(base_data)) base_data <- data
  base_data$significant_qtls <- list()

  regional_data <- base_data$regional_data$regional_summary
  # First try to find q_beta (permutation-based)
  q_col <- if ("q_beta" %in% names(regional_data)) {
    message("Using permutation-based q_beta for significant events in qvalue-based QTL identification")
    "q_beta"
  } else if ("q_bonferroni_min" %in% names(regional_data)) {
    message("Using Bonferroni-based q_bonferroni_min for significant events in qvalue-based QTL identification")
    "q_bonferroni_min"
  } else {
    stop("Neither q_beta nor q_bonferroni_min found in regional data")
  }

  result_name <- sprintf("%s_adjusted_events_qvalue", q_col)

  # Identify significant events
  significant_events <- regional_data %>%
    filter(!!sym(q_col) < params$fdr_threshold) %>%
    pull(molecular_trait_object_id)
  if (length(significant_events) == 0) {
    message(sprintf("No significant events found using %s threshold %g", q_col, params$fdr_threshold))
    base_data$significant_qtls[[result_name]] <- data.frame()
    return(base_data)
  }
  snp_data <- data$qtl_data$data %>%
    filter(molecular_trait_object_id %in% significant_events)

  # Handle empty snp_data
  if (nrow(snp_data) == 0) {
    message("No SNP data for significant events")
    base_data$significant_qtls[[result_name]] <- data.frame()
    return(base_data)
  }

  # Check if the q-value column exists in the data
  q_value_col <- data$qtl_data$column_info$q_col
  if (q_value_col %in% names(snp_data)) {
    message(sprintf("Using existing q-value column '%s' for qvalue-based QTL identification", q_value_col))
    significant_snps <- snp_data %>%
      filter(!!sym(q_value_col) < params$fdr_threshold)
  } else {
    message("Computing q-values for each event's SNPs")
    p_col <- data$qtl_data$column_info$p_col
    significant_snps <- snp_data %>%
      group_by(molecular_trait_object_id) %>%
      mutate(qvalue = safe_qvalue(!!sym(p_col))$qvalues) %>%
      filter(qvalue < params$fdr_threshold) %>%
      ungroup()
  }

  # Create result name based on method used
  base_data$method_name <- q_col

  # Store significant SNPs
  base_data$significant_qtls[[result_name]] <- significant_snps
  message(sprintf(
    "Identified %d significant SNPs from %d events using q-value method with %s as event level significance test",
    nrow(significant_snps), length(significant_events), q_col
  ))
  return(base_data)
}

############################################
# Main Application Controller
############################################
hierarchical_multiple_testing_correction <- function(params) {
  # Step 0: Load all data upfront
  data <- list()
  data$regional_data <- load_regional_data(params)
  recount_n_variants <- FALSE
  if ((params$maf_cutoff > 0 || params$cis_window > 0) || is.null(params$regional_pattern)) {
    recount_n_variants <- TRUE
  }

  data$qtl_data <- load_qtl_data(params, load_n_variants = recount_n_variants)
  regional_base <- paste0(data$qtl_data$file_prefix, ".cis_regional")
  qtl_base <- paste0(data$qtl_data$file_prefix, ".cis_pairs")
  results <- list()

  # Step 1A: Permutation-based local adjustment
  if (data$regional_data$has_permutation) {
    perm_results <- permutation_local_adjustment(data, params)
    perm_results <- perform_global_adjustment(perm_results, params)
    perm_results <- identify_permutation_snps(perm_results, params)
    results$permutation <- list(
      regional_data = list(
        regional_summary = perm_results$regional_data$regional_summary
      ),
      significant_qtls = perm_results$significant_qtls[[1]],
      output_metadata = list(
        regional = paste0(regional_base, ".permutation_adjusted.fdr.gz"),
        events = list(
          path = paste0(regional_base, ".significant_events.permutation_adjusted.tsv.gz"),
          filter_column = "q_beta",
          filter_threshold = params$fdr_threshold
        ),
        qtls = paste0(qtl_base, ".significant_qtl.permutation_adjusted.tsv.gz")
      ),
      global_adjustment_info = list(
        statistics = perm_results$global_adjustment_info$statistics
      )
    )
  }

  # Step 1B: Bonferroni local adjustment (original)
  bonf_orig_results <- bonferroni_local_adjustment(data, params, should_filter = FALSE)
  bonf_orig_results <- perform_global_adjustment(bonf_orig_results, params)
  bonf_orig_results <- identify_bonferroni_fdr_snps(bonf_orig_results, params)
  results$bonferroni_original <- list(
    regional_data = list(
      regional_summary = bonf_orig_results$regional_data$regional_summary
    ),
    significant_qtls = bonf_orig_results$significant_qtls[[1]],
    output_metadata = list(
      regional = paste0(regional_base, ".original_bonferroni_BH_adjusted.fdr.gz"),
      events = list(
        path = paste0(regional_base, ".significant_events.original_bonferroni_BH_adjusted.tsv.gz"),
        filter_column = "fdr_bonferroni_min",
        filter_threshold = params$fdr_threshold
      ),
      qtls = paste0(qtl_base, ".significant_qtl.original_bonferroni_BH_adjusted.tsv.gz")
    ),
    global_adjustment_info = list(
      statistics = bonf_orig_results$global_adjustment_info$statistics %||% list()
    )
  )

  # Step 1C: Bonferroni local adjustment (filtered)
  if (params$maf_cutoff > 0 || params$cis_window > 0) {
    bonf_filt_results <- bonferroni_local_adjustment(data, params, should_filter = TRUE)
    bonf_filt_results <- perform_global_adjustment(bonf_filt_results, params)
    bonf_filt_results <- identify_bonferroni_fdr_snps(bonf_filt_results, params)
    results$bonferroni_filtered <- list(
      regional_data = list(
        regional_summary = bonf_filt_results$regional_data$regional_summary
      ),
      significant_qtls = bonf_filt_results$significant_qtls[[1]],
      output_metadata = list(
        regional = paste0(regional_base, ".filtered_bonferroni_BH_adjusted.fdr.gz"),
        events = list(
          path = paste0(regional_base, ".significant_events.filtered_bonferroni_BH_adjusted.tsv.gz"),
          filter_column = "fdr_bonferroni_min",
          filter_threshold = params$fdr_threshold
        ),
        qtls = paste0(qtl_base, ".significant_qtl.filtered_bonferroni_BH_adjusted.tsv.gz")
      ),
      global_adjustment_info = list(
        statistics = bonf_filt_results$global_adjustment_info$statistics %||% list()
      )
    )
  }

  # Step 1D: Apply q-value per SNP methodology
  qvalue_base <- results$permutation %||% results$bonferroni_filtered %||% results$bonferroni_original
  qvalue_results <- identify_qvalue_snps(data, params, qvalue_base)
  results$qvalue <- list(
    regional_data = list(
      regional_summary = qvalue_results$regional_data$regional_summary
    ),
    significant_qtls = qvalue_results$significant_qtls[[1]],
    output_metadata = list(
      events = list(
        path = paste0(regional_base, ".significant_events.", names(qvalue_results$significant_qtls), ".tsv.gz"),
        filter_column = qvalue_results$method_name,
        filter_threshold = params$fdr_threshold
      ),
      qtls = paste0(qtl_base, ".significant_qtl.", names(qvalue_results$significant_qtls), ".tsv.gz")
    ),
    global_adjustment_info = list(
      statistics = qvalue_results$global_adjustment_info$statistics %||% list()
    )
  )

  # Add summary metadata
  results$summary_metadata <- paste0(regional_base, ".summary.txt")
  return(results)
}
############################################
# Output Management
############################################

write_results <- function(results, out_dir, work_dir, to_cwd = c("regional")) {
  setwd(work_dir)
  summary_stats <- list()
  if (!dir.exists(out_dir)) {
    dir.create(out_dir, recursive = TRUE)
    message(sprintf("Created directory: %s", out_dir))
  }

  # Process each analysis result
  for (result_name in names(results)) {
    if (result_name == "summary_metadata") next

    result <- results[[result_name]]

    # Extract statistics for summary
    if (!is.null(result$global_adjustment_info$statistics)) {
      for (stat_name in names(result$global_adjustment_info$statistics)) {
        full_name <- paste0(result_name, "_", stat_name)
        summary_stats[[full_name]] <- result$global_adjustment_info$statistics[[stat_name]]
      }
    }

    # Write output files based on metadata
    meta <- result$output_metadata

    # Write regional data
    if (!is.null(meta$regional)) {
      # Choose directory based on to_cwd parameter
      dir_path <- if ("regional" %in% to_cwd) "." else out_dir
      full_path <- file.path(dir_path, meta$regional)
      write_delim(result$regional_data$regional_summary, gzfile(full_path), delim = "\t")
    }

    # Write events data
    if (!is.null(meta$events)) {
      dir_path <- if ("events" %in% to_cwd) "." else out_dir
      events_data <- result$regional_data$regional_summary %>%
        filter(!!sym(meta$events$filter_column) < meta$events$filter_threshold)

      full_path <- file.path(dir_path, meta$events$path)
      write_delim(events_data, gzfile(full_path), delim = "\t")
    }

    # Write QTLs data
    if (!is.null(meta$qtls)) {
      dir_path <- if ("qtls" %in% to_cwd) "." else out_dir
      full_path <- file.path(dir_path, meta$qtls)
      write_delim(result$significant_qtls, gzfile(full_path), delim = "\t")
    }
  }

  # Write summary table
  if (length(summary_stats) == 0) {
    summary_tbl <- data.frame(note = "No statistics calculated with available data")
  } else {
    summary_tbl <- data.frame(
      statistic = names(summary_stats),
      value = unlist(summary_stats)
    )
  }

  # Determine if summary should go to CWD
  dir_path <- if ("summary" %in% to_cwd) "." else out_dir
  write.table(summary_tbl, file.path(dir_path, results$summary_metadata),
    sep = "\t", row.names = FALSE, quote = FALSE
  )
}

archive_files <- function(params) {
  archive_dir <- params$archive_dir
  workdir <- params$workdir

  original_wd <- getwd()
  setwd(workdir)

  files_to_archive <- c(
    list.files(pattern = "parquet$", full.names = TRUE),
    list.files(pattern = "regional\\.(fdr\\.gz|summary\\.txt)$", full.names = TRUE)
  )

  existing_files <- files_to_archive[file.exists(files_to_archive)]
  existing_files <- unique(existing_files)

  if (length(existing_files) > 0) {
    # Create archive directory
    if (!dir.exists(archive_dir)) {
      dir.create(archive_dir, recursive = TRUE)
      message(sprintf("Created archive directory: %s", archive_dir))
    }

    # Move files to archive
    moved_count <- 0
    for (file_path in existing_files) {
      target_path <- file.path(archive_dir, basename(file_path))
      if (file.exists(target_path)) {
        message(sprintf("Skipping already archived file: %s", basename(file_path)))
      } else {
        file.rename(file_path, target_path)
        moved_count <- moved_count + 1
      }
    }

    message(sprintf("Archived %d/%d files", moved_count, length(existing_files)))
  } else {
    message("No files found to archive")
  }

  # Restore original working directory
  setwd(original_wd)
}

# source("~/GIT/pecotmr/code/tensorqtl_postprocessor.R")
# results <- hierarchical_multiple_testing_correction(params)
# write_results(results, params$output_dir, params$workdir, to_cwd = "regional")
# archive_files(params)
# saveRDS(results, "test.rds")
