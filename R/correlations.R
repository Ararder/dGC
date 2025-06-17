utils::globalVariables(c("cell"))
#' Prepare dGA format by subsetting cells and genes
#'
#' @param obj dGA object
#' @param ct celltype, "beta cells"
#' @param ct_column column of celltypes
#' @param prop_cells proportion of cells that need to express a gene to be kept
#'
#' @returns a list
#' @export
#'
#' @examples \dontrun{
#' prep_cluster_counts(obj, ct = "beta cells", ct_column = "named_celltype", prop_cells = 0.9)
#' }
prep_cluster_counts <- function(obj, ct, ct_column = "named_celltype", prop_cells = 0.9) {
  # Input validation
  rlang::check_required(ct)
  rlang::check_required(ct_column)
  rlang::check_required(obj)


  if (!ct_column %in% colnames(obj[["obs"]])) {
    stop("Column '", ct_column, "' not found in obs data")
  }
  obs_df <- obj[["obs"]]
  var_df <- obj[["var"]]

  # Filter cells by cell type
  obs_df <- dplyr::filter(obs_df, .data[[ct_column]] == ct)
  sel_cells <- obs_df[["cell"]]
  n_cells_selected <- length(sel_cells)


  if (n_cells_selected == 0) {
    stop("No cells found for cell type: ", ct)
  }


  cli::cli_inform("Selected {.emph {n_cells_selected}} cells of type '{ct}'")

  # Subset matrix and calculate gene filtering threshold
  M <- obj[["count_matrix"]]
  sub_M <- M[sel_cells, , drop = FALSE]
  min_cells_threshold <- as.integer(n_cells_selected * prop_cells)

  # Filter genes based on expression frequency
  genes_expressed_count <- Matrix::colSums(sub_M > 0)
  gene_index <- genes_expressed_count >= min_cells_threshold
  n_genes_kept <- sum(gene_index)

  cli::cli_alert_success("Keeping {n_genes_kept} genes expressed in > {prop_cells*100}% of cells")

  # Apply gene filtering and log transformation
  sub_M <- sub_M[, gene_index, drop = FALSE]
  sub_M@x <- log2(sub_M@x + 1)
  var_df <- dplyr::tibble(gene = colnames(sub_M)) |> dplyr::semi_join(var_df, by = "gene")

  cli::cli_alert_info("Selected {.emph {length(rownames(sub_M))}} cells and {.emph {length(colnames(sub_M))}} genes")

  list(
    matrix = sub_M,
    obs = obs_df,
    var = var_df
  )
}


#' Calculate the difference in correlation between two conditions
#'
#' @param M matrix (cells x genes)
#' @param obs_df observation data frame with a column for conditions
#' @param by column name in obs_df to split by (default: "status")
#' @param method correlation method (default: "pearson")
#'
#' @returns a matrix of differential correlations
#' @export
#'
#' @examples \dontrun{
#' corr_diff(M, obs_df, by = "status", method = "spearman")
#' }
corr_diff <- function(M, obs_df, by = "status", method = c("pearson", "spearman")) {

  method <- rlang::arg_match(method)
  # Validate inputs
  if (!by %in% colnames(obs_df)) {
    stop("Column '", by, "' not found in obs_df")
  }

  # Filter cells present in count matrix
  cells <- dplyr::filter(obs_df, cell %in% rownames(M))

  if (nrow(cells) == 0) {
    stop("No matching cells found between count_matrix and obs_df")
  }

  # Split cells by condition
  split_cells <- split(cells, cells[[by]])
  condition_names <- names(split_cells)

  if (length(condition_names) != 2) {
    stop("Expected exactly 2 conditions, found: ", length(condition_names))
  }


  # Calculate correlations for each condition
  cells_cond1 <- split_cells[[1]]$cell
  cells_cond2 <- split_cells[[2]]$cell

  cor_cond1 <- stats::cor(as.matrix(M[cells_cond1, , drop = FALSE]), method = method)
  cor_cond2 <- stats::cor(as.matrix(M[cells_cond2, , drop = FALSE]), method = method)

  # Return differential correlation (condition2 - condition1)
  cli::cli_inform("subtracting: '{condition_names[2]}' - '{condition_names[1]}'")
  cor_cond2 - cor_cond1
}




#' Apply bootstrap threshhold
#'
#' @param bootstrap_res output of [corr_bootstrap_diff()]
#' @param true_diff output of [corr_diff()]
#'
#' @returns a matrix
#' @export
#'
#' @examples \dontrun{
#' apply_threshold(bootstrap_res, true_diff)
#' }
apply_threshold <- function(bootstrap_res, true_diff) {
  # Validate input
  if (!is.list(bootstrap_res) || !all(c("cor_max", "cor_min", "cor_abs") %in% names(bootstrap_res))) {
    stop("Invalid bootstrap results format")
  }
  # one mask instead of four
  mask_pass <- (true_diff <  bootstrap_res$cor_min) |
    (true_diff >  bootstrap_res$cor_max)

  # stats before you overwrite
  n_pass      <- sum(mask_pass)
  n_smaller   <- sum(true_diff <  bootstrap_res$cor_min)
  n_larger    <- sum(true_diff >  bootstrap_res$cor_max)

  # in-place overwrite (no extra numeric copy)
  true_diff[!mask_pass] <- 0          # <- modifies original matrix

  cli::cli_alert_info(
    "Out of {length(mask_pass)} gene-gene correlations:
   {round(n_smaller / length(mask_pass) * 100, 1)}% smaller than min,
   {round(n_larger  / length(mask_pass) * 100, 1)}% larger than max,
   {round(n_pass    / length(mask_pass) * 100, 1)}% pass either bound.
   Under the null, {(1 / dim(bootstrap_res$bootstrap_array)[3]) * 100}% should pass."
  )

  true_diff   # renamed output, but still the same memory block

}


generate_random_splits <- function(df_obs, n_reps = 20) {
  n <- nrow(df_obs)
  half <- floor(n / 2)


  # Pre-allocate matrix for efficiency
  random_splits <- matrix(nrow = n, ncol = n_reps)

  for (i in seq_len(n_reps)) {
    idx <- sample(n)
    split_label <- rep("case", n)
    split_label[idx[(half + 1):n]] <- "control"
    random_splits[, i] <- split_label
  }

  # Convert to data frame and add to original data
  colnames(random_splits) <- paste0("rep_", seq_len(n_reps))
  split_df <- dplyr::as_tibble(random_splits)
  dplyr::bind_cols(df_obs, split_df)
}





estimate_graphs <- function(M, normalise_edges = TRUE) {
  if(normalise_edges) {
    M <- M / max(abs(M))
  }
  dissTOM <- WGCNA::TOMdist(as.matrix(M), TOMType = "signed")

}



