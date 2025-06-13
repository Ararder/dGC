read_data <- function(path, type = c("seurat","h5ad", "loom")) {
  type <- rlang::arg_match(type)
  rlang::check_required(path)


  if(type=="seurat") {
    obj <- readr::read_rds(path)

    obs_df <- obj@meta.data |> dplyr::as_tibble(rownames = "cell")
    var_df <- dplyr::tibble(gene = rownames(obj))
    count_matrix <- obj@assays$RNA$counts

  }
  list(
    count_matrix = Matrix::t(count_matrix),
    obs = obs_df,
    var = var_df
  )

}
