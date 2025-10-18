



#' Run differential gene-correlation analysis
#' results are saved in `dir`, which defaults to the working directory
#'
#' @param obj a list created by `prep_cluster_counts()`
#' @param split_by column in `obj$obs` to split samples by condition
#' @param n_iter number of permutations
#' @param fit_models which models to fit
#' @param method correlation method
#' @param formula formula for mixed effect model
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


  all_donors <- dplyr::count(obj$obs, donor) |>
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
