utils::globalVariables(c("donor", "cell", "celltype"))




#' Setup the analysis for differential gene-correlation
#'
#' @param obj an object created by [read_data()]
#' @param n_iter number of iterations for bootstrap
#' @param fit_models regress out donor effects?
#' @param method pearson or spearman for correlation?
#' @param replace allow replacement in permutation?
#' @param dir directory to setup analysis in
#'
#' @returns NULL
#' @export
#'
#' @examples \dontrun{
#' setup_dgc(obj, "case")
#'}
#'
setup_dgc <- function(
    obj,
    n_iter=100,
    fit_models = c("none", "blmer", "lmer","glmer"),
    method = c("pearson","spearman"),
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

  cli::cli_h1("Setting up dGC pipeline")

  all_donors <- obj$obs$donor |> unique()

  cli::cli_alert_info("Found {length(all_donors)} unique donors across both conditions")

  n_per_cond <- dplyr::group_by(obj$obs, .data[["condition"]]) |>
    dplyr::summarise(n = dplyr::n_distinct(donor))

  n_control <- n_per_cond$n[1]
  n_case <- n_per_cond$n[2]

  #
  cli::cli_alert_info("Using {n_control} donors for condition 1 and {n_case} donors for condition 2")


  M <- obj$matrix
  obs <- obj$obs


  # -------------------------------------------------------------------------
  cli::cli_alert_info("Cell by Gene matrix dimensions: {dim(M)}")


  cli::cli_alert("Saving data in {dir}")
  permutation_labels <- purrr::map(1:n_iter, \(i)  sample_donors(all_donors, n_control, n_case, replace))
  readr::write_rds(permutation_labels ,file = file.path(dir, "permutation_labels.rds"), compress = "gz")
  readr::write_rds(M,file = file.path(dir, "M.rds"), compress = "gz")
  readr::write_rds(obs,file = file.path(dir, "obs.rds"), compress = "gz")
  readr::write_rds(list(fit_models = fit_models,method = method), file = file.path(dir, "arg.rds"))

}




#' Generate permutations
#'
#' @param dir directory where [setup_dgc()] has created data
#' @param ncores number of cores for residualizing
#' @returns NULL
#' @export
#'
#' @examples \dontrun{
#' run_permutations(100, labels, M, obs)
#' }
#'
run_permutations <- function(
    dir,
    ncores =1
    ) {


  params <- readr::read_rds(file.path(dir, "arg.rds"))

  # 2. Extract variables (and convert formula back)
  method <- params$method
  fit_models  <- params$fit_models


  permutation_labels <- readr::read_rds(file = file.path(dir, "permutation_labels.rds"))
  obs <- readr::read_rds(file = file.path(dir, "obs.rds"))
  n_iter <- length(permutation_labels)

  purrr::walk(1:n_iter, purrr::in_parallel(function(i) {
    cli::cli_h1("Generating permutation {i}")

    M <- readr::read_rds(file = file.path(dir, "M.rds"))
    permutation_labels <- readr::read_rds(file = file.path(dir, "permutation_labels.rds"))
    obs <- readr::read_rds(file = file.path(dir, "obs.rds"))

    d1 <- permutation_labels[[i]][[1]]
    d2 <- permutation_labels[[i]][[2]]
    diff <- dGC::corr_diff(
      M = M,
      obs_1 = dplyr::filter(obs, donor %in% d1),
      obs_2 = dplyr::filter(obs, donor %in% d2),
      fit_models = fit_models,
      method = method,
      ncores = ncores
    )

    readr::write_rds(diff, file = file.path(dir, paste0("perm_diff_", i, ".rds")), compress = "gz")
    rm(diff)
    gc()

  },
  dir = dir,
  fit_models = fit_models,
  method = method,
  ncores = ncores
  ))
}

#' Calculate the true correlation difference between conditions
#'
#' @param dir directory where [setup_dgc()] has created data
#' @param ncores number of cores for residualizing
#'
#' @returns NULL
#' @export
#'
#' @examples \dontrun{
#' run_real_diff(tempdir(), ncores=3)
#' }
run_real_diff <- function(dir, ncores=1) {
  params <- readr::read_rds(file.path(dir, "arg.rds"))
  obs <- readr::read_rds(file = file.path(dir, "obs.rds"))
  M <- readr::read_rds(file = file.path(dir, "M.rds"))
  # 2. Extract variables (and convert formula back)
  method <- params$method
  fit_models  <- params$fit_models

  conds <- split(obs, obs$condition)


  # -------------------------------------------------------------------------
  cli::cli_h1("Calculating the true correlation difference matrix")
  cli::cli_inform("Using: {ncores} cores, number of genes: {ncol(M)}")
  # and unique donors in each condition
  cli::cli_inform("Condition 1: {length(unique(conds[[1]]$donor))} unique donors ( {names(conds)[1]} ), {nrow(conds[[1]])} cells")
  cli::cli_inform("Condition 2: {length(unique(conds[[2]]$donor))} unique donors ( {names(conds)[2]} ), {nrow(conds[[2]])} cells)")

  real_diff <- corr_diff(
    M = M,
    obs_1 = conds[[1]],
    obs_2 = conds[[2]],
    fit_models = fit_models,
    method = method,
    ncores = ncores
  )

  readr::write_rds(real_diff, file = file.path(dir, "real_diff.rds"))


}




#' Calculate the difference in Gene-Gene correlation across two conditions
#'
#' @param M a matrix of gene expression values, with genes as columns and cells as rows
#' @param obs_1 a data frame containing cell identifiers and donor information for the first condition
#' @param obs_2 a data frame containing cell identifiers and donor information for the second condition
#' @param ncores number of cores to use for compute_residuals
#' @param fit_models model to residualise expression
#' @param method Pearson or spearman for correlations
#'
#' @returns a matrix of the difference in correlation between the two conditions
#' @export
#'
#' @examples \dontrun{
#' corr_diff(M, obs_1, obs_2, fit_models = "blmer", method = "pearson")
#' }
corr_diff <- function(
    M,
    obs_1,
    obs_2,
    fit_models = c("none", "blmer", "lmer","glmer"),
    method = c("pearson", "spearman"),
    ncores = 1
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
    m1 <- compute_residuals(matrix = M, cells = obs_1$cell, donor_vec = obs_1$donor, engine = fit_models, ncores=ncores)
    cli::cli_alert_info("Condition 2: {nrow(obs_2)} cells from {length(unique(obs_2$donor))} donors")
    m2 <- compute_residuals(matrix = M, cells = obs_2$cell, donor_vec = obs_2$donor, engine = fit_models, ncores=ncores)
  }

  cor_cond1 <- stats::cor(as.matrix(m1), method = method)
  cor_cond2 <- stats::cor(as.matrix(m2), method = method)

  cor_cond2 - cor_cond1
}




#' compute residuals of a single-cell data matrix
#'
#' @param matrix count matrix
#' @param engine method to fit the model, one of "blmer", "lmer", or "glmer"
#' @param cells a vector of cell identifiers to use, if NULL all cells are used
#' @param donor_vec a vector of donor identifiers, must be the same length as cells
#' @param ncores number of cores to use
#' @returns a list()
#' @export
#'
#' @examples \dontrun{
#' compute_residuals(count_matrix)
#' }
compute_residuals <- function(matrix, engine = c("blmer", "lmer","glmer"), ncores = 1, cells=NULL, donor_vec) {
  engine <- rlang::arg_match(engine)

  if(!is.null(cells)) {
    matrix <- Matrix::Matrix(matrix)[cells, ]
    stopifnot(length(cells) == length(donor_vec))
    stopifnot(nrow(matrix) == length(cells))
  }

  if(ncores == 1) {
    res <- purrr::map(1:ncol(matrix), \(idx) dGC::fit_model(expr= matrix[, idx], donor_vec = donor_vec,engine = engine), .progress = list(type = "tasks"))

  } else {
    future::plan(future::multisession, workers = ncores)


    res <- furrr::future_map(1:ncol(matrix), \(idx) {
      fit_model(expr = Matrix::Matrix(matrix)[, idx], donor_vec = donor_vec, engine = engine)
    }, .progress = TRUE)

  }

  M <- do.call(cbind, res)
  rownames(M) <- cells
  colnames(M) <- colnames(matrix)

  M
}

#' Fit a regression model to gene expression
#'
#' @param expr column of gene expression
#' @param donor_vec donor vector
#' @param engine engine to ues
#'
#' @returns a vector of residuals
#' @export
#'
#' @examples \dontrun{
#' fit_model(gene, donors, "blmer")
#' }
fit_model <- function(expr, donor_vec, engine) {

  df <- dplyr::tibble(expr = expr, donor = donor_vec)
  if (engine == "blmer") {
    m <- blme::blmer(expr ~ 1 + (1|donor), data = df)
  } else if (engine == "lmer") {
    m <- lme4::lmer(expr ~ 1 + (1|donor), data = df)
  } else if (engine == "glmer") {
    m <- lme4::glmer(expr ~ 1 + (1|donor), data = df, family = stats::gaussian())
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





#' Get the empirical p-value for each link
#'
#' @param dir directory where permutations are stored
#'
#' @returns NULL
#' @export
#'
#' @examples \dontrun{
#' get_empirical_p(tempdir())
#' }
get_empirical_p <- function(dir) {
  P_path <- fs::dir_ls(dir, glob = "*perm_diff_*.rds")
  R <- readr::read_rds(file.path(dir, "real_diff.rds"))

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

