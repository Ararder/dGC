
utils::globalVariables(c("donor", "cell"))


#' Run differential gene-correlation analysis
#' results are saved in `dir`, which defaults to the working directory
#'
#' @param obj a list created by `prep_cluster_counts()`
#' @param split_by column in `obj$obs` to split samples by condition
#' @param n_iter number of permutations
#' @param fit_models which models to fit
#' @param method correlation method
#' @param formula formula for mixed effect model
#' @param min_cells_per_donor Minimum number of cells per donor
#' @param n_ctrl number of control donors to sample
#' @param n_case number of case donors to sample
#' @param replace sample with replacement
#' @param dir output directory
#'
#' @returns NULL
#' @export
#'
#' @examples \dontrun{
#' run_dgc(object, split_by = "disease", n_iter=10, dir = tempdir())
#' }
run_dgc <- function(
    obj,
    split_by = "condition",
    n_iter=100,
    fit_models = c("none", "blmer", "lmer","glmer"),
    method = c("pearson","spearman"),
    formula = stats::as.formula("expr ~ 1 + (1|donor)"),
    min_cells_per_donor = 10,
    n_ctrl=10,
    n_case =10,
    replace=FALSE,
    dir = NULL

  ) {

  # check args --------------------------------------------------------------
  method <- rlang::arg_match(method)
  fit_models <- rlang::arg_match(fit_models)
  if(is.null(dir)) {
    dir <- getwd()
  } else {
    file.exists(dir) || cli::cli_abort("Output directory: {dir} does not exist")
  }

  cli::cli_h1("Starting dGC pipeline")

  # min cells per donor


  all_donors <- dplyr::count(obj$obs, donor) |>
    dplyr::filter(.data[["n"]] >= min_cells_per_donor) |>
    dplyr::pull(donor)

  M <- obj$matrix
  obs <- obj$obs


  cli::cli_inform("Splitting by {split_by}")
  split_by %in% colnames(obs) || cli::cli_abort("Could not find column {split_by} in obs data.frame")
  cond <- split(obs, obs[[split_by]])


  permutation_labels <- purrr::map(1:n_iter, \(i)  sample_donors(all_donors, n_ctrl, n_case, replace))
  readr::write_rds(permutation_labels ,file = file.path(dir, "permutation_labels.rds"), compress = "gz")



  real_diff <- corr_diff(
    M = M,
    obs_1 = cond[[1]],
    obs_2 = cond[[2]],
    fit_models = fit_models,
    formula = formula,
    method = method,
    ncores = 1
  )
  readr::write_rds(real_diff, file = file.path(dir, "real_diff.rds"), compress = "gz")

  # minimise memory
  rm(real_diff)
  gc()

  run_permutations(
    n_iter = n_iter,
    permutation_labels = permutation_labels,
    M = M,
    obs = obs,
    fit_models = fit_models,
    formula = formula,
    method = method

  )




}

#' Generate permutations
#'
#' @param n_iter number of iterations
#' @param permutation_labels permutation list
#' @param M Matrix
#' @param obs observation df
#' @inheritParams  run_dgc
#'
#' @returns NULL
#' @export
#'
#' @examples \dontrun{
#' run_permutations(100, labels, M, obs)
#' }
#'
run_permutations <- function(n_iter, permutation_labels, M, obs, fit_models, formula, method) {
  purrr::walk(1:n_iter, purrr::in_parallel(function(i) {
    cli::cli_h1("Generating permutation {i}")

    d1 <- permutation_labels[[i]][[1]]
    d2 <- permutation_labels[[i]][[2]]
    diff <- corr_diff(
      M = M,
      obs_1 = dplyr::filter(obs, donor %in% d1),
      obs_2 = dplyr::filter(obs, donor %in% d2),
      fit_models = fit_models,
      formula = formula,
      method = method,
      ncores = 1
    )

    readr::write_rds(diff, file = file.path(dir, paste0("perm_diff_", i, ".rds")), compress = "gz")
    rm(diff)
    gc()

  }, .progress = list(type = "tasks", name = "computing permutations"),
  permutation_labels = permutation_labels,
  corr_diff = dGC::corr_diff,
  M = M,
  dir = dir,
  obs = obs,
  fit_models = fit_models,
  formula = formula,
  method = method
  ))
}




#' Calculate the difference in Gene-Gene correlation across two conditions
#'
#' @param M a matrix of gene expression values, with genes as columns and cells as rows
#' @param obs_1 a data frame containing cell identifiers and donor information for the first condition
#' @param obs_2 a data frame containing cell identifiers and donor information for the second condition
#' @inheritParams run_dgc
#' @param ncores number of cores to use for parallel processing, default is 1
#'
#' @returns a matrix of the difference in correlation between the two conditions
#' @export
#'
#' @examples \dontrun{
#' corr_diff(M, obs_1, obs_2, fit_models = "blmer", method = "pearson", ncores = 4)
#' }
corr_diff <- function(
    M,
    obs_1,
    obs_2,
    fit_models = c("none", "blmer", "lmer","glmer"),
    method = c("pearson", "spearman"),
    formula = stats::as.formula("expr ~ 1 + (1|donor)"),
    ncores=1
) {
  cli::cli_h2("Calculating a correlation difference matrix")
  method <- rlang::arg_match(method)
  fit_models <- rlang::arg_match(fit_models)

  if(fit_models == "none") {
    m1 <- Matrix::Matrix(M)[obs_1$cell, , drop = FALSE]
    m2 <- Matrix::Matrix(M)[obs_2$cell, , drop = FALSE]
  } else {
    cli::cli_alert_info("Residualizing expression values using a {fit_models} model")
    cli::cli_alert_info("Condition 1: {nrow(obs_1)} cells from {length(unique(obs_1$donor))} donors")
    m1 <- compute_residuals(matrix = M, cells = obs_1$cell, donor_vec = obs_1$donor, ncores = ncores, engine = fit_models, formula=formula)
    cli::cli_alert_info("Condition 2: {nrow(obs_2)} cells from {length(unique(obs_2$donor))} donors")
    m2 <- compute_residuals(matrix = M, cells = obs_2$cell, donor_vec = obs_2$donor, ncores = ncores, engine = fit_models, formula=formula)
  }

  cor_cond1 <- stats::cor(as.matrix(m1), method = method)
  cor_cond2 <- stats::cor(as.matrix(m2), method = method)

  cor_cond2 - cor_cond1
}




#' compute residuals of a single-cell data matrix
#'
#' @param matrix count matrix
#' @param engine method to fit the model, one of "blmer", "lmer", or "glmer"
#' @param formula a formula to use for the model, default is "expr ~ 1 + (1|donor)"
#' @param cells a vector of cell identifiers to use, if NULL all cells are used
#' @param donor_vec a vector of donor identifiers, must be the same length as cells
#' @param ncores number of cores to use for parallel processing, default is 1
#'
#' @returns a list()
#' @export
#'
#' @examples \dontrun{
#' compute_residuals(count_matrix)
#' }
compute_residuals <- function(matrix, engine = c("blmer", "lmer","glmer"), formula = stats::as.formula("expr ~ 1 + (1|donor)"), cells=NULL, donor_vec, ncores=1) {
  engine <- rlang::arg_match(engine)
  stopifnot(length(cells) == length(donor_vec))
  if(!is.null(cells)) {
    matrix <- Matrix::Matrix(matrix)[cells, ]
  }
  stopifnot(nrow(matrix) == length(cells))

  if(ncores > 1) {
    future::plan(future::multisession, workers = ncores)
    res <- furrr::future_map(1:ncol(matrix), \(idx) {
      fit_model(expr = Matrix::Matrix(matrix)[, idx], donor_vec = donor_vec, formula = formula, engine = engine)
    }, .progress = TRUE)

  } else {
    res <- purrr::map(1:ncol(matrix), \(idx) {
      fit_model(expr= Matrix::Matrix(matrix)[, idx], donor_vec = donor_vec, formula = formula, engine = engine)
    }, .progress = list(type = "tasks"))

  }

  M <- do.call(cbind, res)
  rownames(M) <- cells
  colnames(M) <- colnames(matrix)

  M
}

fit_model <- function(expr, donor_vec, formula, engine) {
  df <- dplyr::tibble(expr = expr, donor = donor_vec)
  if (engine == "blmer") {
    m <- blme::blmer(formula, data = df)
  } else if (engine == "lmer") {
    m <- lme4::lmer(formula, data = df)
  } else if (engine == "glmer") {
    m <- lme4::glmer(formula, data = df, family = stats::gaussian())
  }
  stats::residuals(m)
}

sample_donors <- function(all_donors, n_ctrl, n_case, replace) {
  group1 <- sample(all_donors, size = n_ctrl, replace = replace)
  group2 <- sample(setdiff(all_donors, group1), size = n_case, replace = replace)

  list(
    group1 = group1,
    group2 = group2
  )
}




mask_from_perm <- function(R, P_path) {


  n_perms <- length(P_path)
  mask <- matrix(FALSE, ncol = ncol(R),nrow = nrow(R))

  for(i in seq_along(P_path)) {
    perm <- readr::read_rds(P_path[[i]])
    mask <- mask + (abs(R) > abs(perm))
  }

  emp_p <- (mask + 1) / (n_perms + 1)
  n_edges <- length(mask)


  alpha <- 0.05
  n_edges_sig <- sum(emp_p < alpha)

  cli::cli_alert_info("{round(n_edges_sig / n_edges,4)} edges significant at alpha = {alpha}")

  emp_p

}

