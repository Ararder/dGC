#' validate data format
#'
#' @param data input list
#'
#' @returns a list
#' @export
#'
#' @examples \dontrun{
#' data <- validate_data(data)
#' }
validate_data <- function(data) {
  # ───────────────────────────────────────── checks 0: list & names
  if (!is.list(data))
    rlang::abort("`data` must be a list.")

  req <- c("matrix", "obs", "var")
  miss <- setdiff(req, names(data))
  if (length(miss))
    rlang::abort(c("`data` is missing component(s):", paste(miss, collapse = ", ")))

  mat <- data$matrix
  obs <- data$obs
  var <- data$var

  # ───────────────────────────────────────── checks 1: types
  if (!(is.matrix(mat) || inherits(mat, "Matrix")))
    rlang::abort("`count_matrix` must be a base or sparse matrix.")
  if (!is.data.frame(obs))
    rlang::abort("`obs` must be a data.frame / tibble.")
  if (!is.data.frame(var))
    rlang::abort("`var` must be a data.frame / tibble.")

  # ───────────────────────────────────────── checks 2: mandatory columns
  if (!"cell"  %in% names(obs))
    rlang::abort("`obs` must contain a `cell` column.")
  if (!"gene"  %in% names(var))
    rlang::abort("`var` must contain a `gene` column.")

  # ───────────────────────────────────────── checks 3: row / column names
  if (is.null(rownames(mat)))
    rlang::abort("`count_matrix` must have rownames (cells).")
  if (is.null(colnames(mat)))
    rlang::abort("`count_matrix` must have colnames (genes).")

  # ───────────────────────────────────────── subset mat if obs/var were filtered
  if (!all(obs$cell %in% rownames(mat)))
    rlang::abort("`obs$cell` contains unknown cells not found in count_matrix.")
  if (!all(var$gene %in% colnames(mat)))
    rlang::abort("`var$gene` contains unknown genes not found in count_matrix.")

  mat <- mat[obs$cell, var$gene, drop = FALSE]

  # ───────────────────────────────────────── reorder obs / var if needed
  if (!identical(obs$cell, rownames(mat))) {
    obs <- obs[match(rownames(mat), obs$cell), , drop = FALSE]
  }
  if (!identical(var$gene, colnames(mat))) {
    var <- var[match(colnames(mat), var$gene), , drop = FALSE]
  }

  list(matrix = mat, obs = obs, var = var)
}
