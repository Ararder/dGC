
utils::globalVariables(c("donor", ""))




#' Compute permutations of differential gene correlation
#'
#' @param obj a list in the format of [read_data()]
#' @param n_iter number of iterations to perform
#' @param residualise whether to fit a linear-mixed-model to remove donor effects from each gene
#' @param method method for correlation, pearson or spearman
#' @param n_ctrl number of donors in "control" correlation group
#' @param n_case number of donors in "case" correlation group
#' @param min_cells minimum number of cells per donor to include in the analysis
#' @param replace whether to sample donors with replacement
#' @param ncores number of cores to use for parallel processing
#'
#' @returns a list of matrices, each matrix is the difference in correlation between the case and control groups for each iteration
#' @export
#'
#' @examples \dontrun{
#' permuts <- corr_permute(obj)
#' }
#'
corr_permute <- function(
    obj, n_iter=100, residualise = FALSE, method = c("pearson","spearman"),
    n_ctrl=10, n_case =10, min_cells = 20, replace=FALSE,ncores=1
    ) {
  method <- rlang::arg_match(method)

  all_donors <- dplyr::count(obj$obs, donor) |>
    dplyr::filter(n >= min_cells) |>
    dplyr::pull(donor)

  M <- obj$matrix
  obs <- obj$obs


  if(ncores == 1) {
    output <- purrr::map(1:n_iter, function(i) {
      .perm_iter(
        M = M,obs = obs,
        all_donors = all_donors, n_ctrl = n_ctrl, n_case = n_case,
        replace = replace, fit_models = residualise, method = method,ncores=1
      )

    }, .progress = list(type = "tasks", name = "computing permutations"))

  } else {

    future::plan(future::multisession, workers = ncores)
    output <- furrr::future_map(1:n_iter, \(i) {
      .perm_iter(
        M = M,obs = obs,
        all_donors = all_donors, n_ctrl = n_ctrl, n_case = n_case,
        replace = replace, fit_models = residualise, method = method,ncores=1
        )

    }, .progress =TRUE, .options = furrr::furrr_options(seed = TRUE))
  }

  output

}



.perm_iter <- function(M, obs, all_donors, n_ctrl, n_case, replace, fit_models, method,ncores=1) {
  group1 <- sample(all_donors, size = n_ctrl, replace = replace)
  group2 <- sample(setdiff(all_donors, group1), size = n_case, replace = replace)

  if(fit_models) {
    cell_info <- dplyr::filter(obs, donor %in% group1)
    cell_info2 <- dplyr::filter(obs, donor %in% group2)
    m1 <- compute_residuals(matrix = M,cells = cell_info$cell,donor_vec = cell_info$donor, ncores = ncores)
    m2 <- compute_residuals(matrix = M,cells = cell_info2$cell,donor_vec = cell_info2$donor, ncores = ncores)

  } else {
    m1 <- Matrix::Matrix(M)[obs$donor %in% group1, ]
    m2 <- Matrix::Matrix(M)[obs$donor %in% group2, ]
  }


  stats::cor(as.matrix(m2), method = method) - stats::cor(as.matrix(m1), method = method)

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





#' calculate a mask from a list of permutations
#'
#' @param P list of permutations, from [corr_permute()]
#' @param R a matrix of correlations, e.g. the output from [corr_diff()]
#'
#' @returns a matrix of the same dimensions as R, with values between 0 and 1 indicating the proportion of permutations that were greater than the observed value in R
#' @export
#'
#' @examples \dontrun{
#' mask_from_perm(perm, real)
#' }
mask_from_perm <- function(P, R) {
  # pos_mask <- matrix(0L, ncol = ncol(R),nrow = nrow(R))
  # neg_mask <- matrix(0L, ncol = ncol(R),nrow = nrow(R))
  # pos_edges <- R > 0
  # neg_edges <- R < 0
  mask <- matrix(0L, ncol = ncol(R),nrow = nrow(R))

  for(i in seq_along(P)) {
    perm <- P[[i]]
    mask <- mask + (abs(R) > abs(perm))

    # # positive edges
    # pos <- matrix(FALSE, ncol = ncol(R),nrow = nrow(R))
    # pos[pos_edges] <- R[pos_edges] > perm[pos_edges]
    # pos_mask <- pos_mask + pos
    #
    # #
    # neg <- matrix(FALSE, ncol = ncol(R),nrow = nrow(R))
    # neg[neg_edges] <- R[neg_edges] > perm[neg_edges]
    # neg_mask <- neg_mask + neg




  }

  emp_p <- 1 - ( mask  / (length(P)+1) )
  tot <- length(mask)
  alpha <- 0.01
  n_edges <- sum(emp_p < alpha)
  cli::cli_alert_info(
   "{round(n_edges / tot, 3)}% of edges are significant at alpha = {alpha}
    Detected {n_edges} edges, under which {round(alpha*tot)} are expected by chance.
   False-proportion rate is estimated to be {round(((alpha*tot) / n_edges),3)*100}%"
  )
  alpha <- 0.025
  n_edges <- sum(emp_p < alpha)
  cli::cli_alert_info(
  "{round(n_edges / tot, 3)}% of edges are significant at alpha = {alpha}
    Detected {n_edges} edges, under which {round(alpha*tot)} are expected by chance.
   False-proportion rate is estimated to be {round(((alpha*tot) / n_edges),3)*100}%"
  )


}






perm_gene_test <- function(permutations) {
  G = nrow(permutations[[1]])
  N = length(permutations)

  purrr::map(seq_along(permutations), \(idx){

    mask_x <- mask_from_perm(
      P = permutations[-idx],
      R = permutations[[idx]]
    )



  })


}





