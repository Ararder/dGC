
validate_data <- function(data) {
  # ───────────────────────────────────────── checks 0: list & names
  if (!is.list(data))
    rlang::abort("`data` must be a list.")

  req <- c("count_matrix", "obs", "var")
  miss <- setdiff(req, names(data))
  if (length(miss))
    rlang::abort(c("`data` is missing component(s):", paste(miss, collapse = ", ")))

  mat <- data$count_matrix
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

  # identical sets?
  if (!setequal(obs$cell, rownames(mat)))
    rlang::abort("`obs$cell` must contain exactly the same values as rownames(count_matrix).")
  if (!setequal(var$gene, colnames(mat)))
    rlang::abort("`var$gene` must contain exactly the same values as colnames(count_matrix).")

  # ───────────────────────────────────────── reorder obs / var if needed
  if (!identical(obs$cell, rownames(mat))) {
    obs <- obs[match(rownames(mat), obs$cell), , drop = FALSE]
  }
  if (!identical(var$gene, colnames(mat))) {
    var <- var[match(colnames(mat), var$gene), , drop = FALSE]
  }

  # return re-ordered (invisible so you can `data <- validate_data(data)`)
  invisible(list(count_matrix = mat, obs = obs, var = var))
}
