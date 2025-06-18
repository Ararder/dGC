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
  M <- obj[["matrix"]]
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

