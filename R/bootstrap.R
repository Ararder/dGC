
utils::globalVariables(c("donor", ""))




#' Compute permutations of differential gene correlation
#'
#' @param obj a list in the format of [read_data()]
#' @param n_iter number of iterations to perform
#' @param fit_models whether to fit a linear-mixed-model to remove donor effects from each gene
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
    obj, n_iter=100, fit_models = c("none", "blmer", "lmer","glmer"), method = c("pearson","spearman"),
    n_ctrl=10, n_case =10, replace=FALSE,ncores=1
    ) {
  method <- rlang::arg_match(method)

  all_donors <- dplyr::count(obj$obs, donor) |>
    dplyr::pull(donor)

  M <- obj$matrix
  obs <- obj$obs


  if(ncores == 1) {
    output <- purrr::map(1:n_iter, function(i) {
      corr_permute_diff(
        M = M,obs = obs,all_donors = all_donors, n_ctrl = n_ctrl, n_case = n_case,
        replace = replace, fit_models = fit_models, method = method,ncores=1
      )

    }, .progress = list(type = "tasks", name = "computing permutations"))

  } else {

    future::plan(future::multisession, workers = ncores)
    output <- furrr::future_map(1:n_iter, \(i) {
      corr_permute_diff(
        M = M,obs = obs,all_donors = all_donors, n_ctrl = n_ctrl, n_case = n_case,
        replace = replace, fit_models = fit_models, method = method,ncores=1
        )
    }, .progress =TRUE, .options = furrr::furrr_options(seed = TRUE))
  }

  output

}



corr_permute_diff <- function(M, obs, all_donors, n_ctrl, n_case, replace, fit_models, method,ncores=1) {
  # sample donors in the specified structure
  group1 <- sample(all_donors, size = n_ctrl, replace = replace)
  group2 <- sample(setdiff(all_donors, group1), size = n_case, replace = replace)

  corr_diff(
    M = M,
    obs_1 = dplyr::filter(obs, donor %in% group1),
    obs_2 = dplyr::filter(obs, donor %in% group2),
    fit_models = fit_models,
    method = method,
    ncores = ncores
  )


}


#' Calculate the difference in Gene-Gene correlation across two conditions
#'
#' @param M a matrix of gene expression values, with genes as columns and cells as rows
#' @param obs_1 a data frame containing cell identifiers and donor information for the first condition
#' @param obs_2 a data frame containing cell identifiers and donor information for the second condition
#' @param fit_models a character vector indicating whether to fit a linear mixed model to the data, one of "none", "blmer", "lmer", or "glmer"
#' @param method a character string indicating the method to use for correlation calculation, one of "pearson" or "spearman"
#' @param ncores number of cores to use for parallel processing, default is 1
#'
#' @returns a matrix of the difference in correlation between the two conditions
#' @export
#'
#' @examples \dontrun{
#' corr_diff(M, obs_1, obs_2, fit_models = "blmer", method = "pearson", ncores = 4)
#' }
corr_diff <- function(M, obs_1, obs_2, fit_models =c("none", "blmer", "lmer","glmer"), method = c("pearson", "spearman"), ncores=1) {

  method <- rlang::arg_match(method)
  fit_models <- rlang::arg_match(fit_models)

  if(fit_models == "none") {
    m1 <- Matrix::Matrix(M)[obs_1$cell, , drop = FALSE]
    m2 <- Matrix::Matrix(M)[obs_2$cell, , drop = FALSE]
  } else {
    m1 <- compute_residuals(matrix = M, cells = obs_1$cell, donor_vec = obs_1$donor, ncores = ncores, engine = fit_models)
    m2 <- compute_residuals(matrix = M, cells = obs_2$cell, donor_vec = obs_2$donor, ncores = ncores, engine = fit_models)
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





#' calculate a mask from a list of permutations
#'
#' @param P list of permutations, from [corr_permute()]
#' @param R a matrix of correlations, e.g. the output from [corr_diff()]
#' @param alpha significance level, default is 0.01
#' @returns a matrix of the same dimensions as R, with values between 0 and 1 indicating the proportion of permutations that were greater than the observed value in R
#' @export
#'
#' @examples \dontrun{
#' mask_from_perm(perm, real)
#' }
mask_from_perm <- function(P, R,alpha=0.01) {

  mask <- matrix(FALSE, ncol = ncol(R),nrow = nrow(R))
  mask_pos <- matrix(FALSE, ncol = ncol(R),nrow = nrow(R))
  pos_increment <- matrix(FALSE, ncol = ncol(R),nrow = nrow(R))
  mask_neg <- matrix(FALSE, ncol = ncol(R),nrow = nrow(R))
  neg_increment <- matrix(FALSE, ncol = ncol(R),nrow = nrow(R))

  pos_edges <- R > 0
  neg_edges <- R < 0

  for(i in seq_along(P)) {
    perm <- P[[i]]
    mask <- mask + (abs(R) > abs(perm))

    # for sign tested, test only against instances where the sign is the same
    pos_test_mask <- perm > 0 & pos_edges
    test_result <- R > perm
    mask_pos[pos_test_mask] <- mask_pos[pos_test_mask] + test_result[pos_test_mask]
    pos_increment[pos_test_mask] <- pos_increment[pos_test_mask] + 1

    neg_test_mask <- perm < 0 & neg_edges
    test_result <- R < perm
    mask_neg[neg_test_mask] <- mask_neg[neg_test_mask] + test_result[neg_test_mask]
    neg_increment[neg_test_mask] <- neg_increment[neg_test_mask] + 1


  }

  emp_p_neg <- 1 - ( mask_neg  / (neg_increment +1) )
  emp_p_pos <- 1 - ( mask_pos  / (pos_increment +1) )
  emp_p <- 1 - ( mask  / (length(P)+1) )

  sum(emp_p_neg < 0.05) / length(emp_p_neg)
  sum(emp_p_pos < 0.05) / length(emp_p_pos)
  sum(emp_p < 0.05) / length(emp_p_pos)









  emp_p <- 1 - ( mask  / (length(P)+1) )
  tot <- length(mask)
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
  emp_p


}






perm_gene_test <- function(P, R) {
  G = nrow(R)
  N = length(P)


  # future::plan(future::multisession, workers = ncores)
  gene_null_distrib <- purrr::map(seq_along(P[1:25]), \(idx){
    p_iter <- P[[idx]]

    mask_x <- mask_from_perm(
      P = P[-idx],
      R = p_iter
    )


    R[mask_x > 0.025] <- 0
    R[mask_x < 0.025] <- 1

    edges <- rowSums(R)

    as.matrix(edges)

  },.progress = list(type = "tasks"))


  # null_dist <- purrr::reduce(gene_null_distrib, cbind)
  # row_max <- apply(null_dist, 1, max)
  #
  # mask <- mask_from_perm(P, R)
  # R[mask > 0.025] <- 0
  # R[mask < 0.025] <- 1
  # edges <- rowSums(R)
  # obs <- dplyr::tibble(links = edges, genes = names(edges))
  # null_res <- dplyr::tibble(links = row_max, genes = names(row_max))
  # dplyr::inner_join(obs, null_res,by = "genes") |>
  #   dplyr::mutate(diff = links.x - links.y) |>
  #   dplyr::filter(links.x > links.y) |>
  #   print(n = 21)








}





