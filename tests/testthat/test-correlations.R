test_that("reduce matrix works", {
  skip()
  data <- read_data("~/projects/dcgna/workflow/t2d/qc_merged_seurat.rds")


  data
  ct <- "alpha cells"
  ct_column <- "named_celltype"

  reduced <- reduce_matrix(data, ct, ct_column = ct_column, prop_cells = 0.9)

  true_diff <- corr_diff(
    reduced$log2_matrix,
    reduced$obs,
    by = "dataset",
    method = "pearson"
  )

  generate_random_splits(reduced$obs, 50)


})





test_that("The full pipeline", {
  skip()

  data <- read_data("~/projects/dcgna/workflow/t2d/qc_merged_seurat.rds")
  data$count_matrix <- (data$count_matrix / rowSums(data$count_matrix))*10^6



  reduced <- prep_cluster_counts(data, "beta cells", ct_column = "named_celltype", prop_cells = 0.9)

  all <-
    reduced$obs |> dplyr::count(donor, status) |>
    dplyr::filter(n >= 20)

  reduced$obs <- dplyr::filter(reduced$obs, donor %in% all$donor)
  reduced <- validate_data(reduced)
  M <- compute_residuals(reduced$matrix, cells = reduced$obs$cell, donor_vec = reduced$obs$donor,ncores = 6)
  reduced$matrix <- M

  real_diff <- corr_diff(reduced$matrix, reduced$obs, by = "status", method = "pearson")
  permutations <- corr_permute(reduced, residualise = FALSE, n_case = 11, n_ctrl = 11, n_iter = 256, ncores=6)



  idx <- 26

  t_mask <- mask_from_perm(P = permutations[-idx], R = permutations[[idx]])
  sum(t_mask < 0.025) / length(t_mask)

  l_mask <- mask_from_perm(P = permutations, R = real_diff)
  sum(l_mask < 0.025) / length(l_mask)







})

test_that("plotting works", {
  skip()
  idx <- 25
  real_diff <- permutations[[idx]]
  mask <- mask_from_perm(P = permutations[-25], R = real_diff)

  M <- real_diff
  M[mask < 0.05] <- 0
  normalise_edges <- TRUE
  if(normalise_edges) {
    M <- M / max(abs(M))
  }

  dissTOM <- WGCNA::TOMdist(as.matrix(M), TOMType = "signed")
  hierTOM <- stats::hclust(stats::as.dist(dissTOM), method = "average")



  cutoff <- c(0.85, 0.981, 0.987, 0.992, 0.992)
  deep_split <- c(FALSE, FALSE, FALSE, FALSE, TRUE)
  module_list <- purrr::map2(cutoff, deep_split, \(cut_off, deep) dynamicTreeCut::cutreeDynamic(hierTOM, method = "tree", minClusterSize = 30, cutHeight = cut_off, deepSplit = deep) |> WGCNA::labels2colors())


  WGCNA::plotDendroAndColors(
    hierTOM,
    data.frame(module_list),
    c("Modules1", "Modules2", "Modules3", "Modules4", "Modules5"),
    dendroLabels = FALSE,
    hang = 0.03,
    addGuide = TRUE,
    guideHang = 0.05,
    main = "Hierarchical clustering",
    autoColorHeight = FALSE
  )



  ll <- dplyr::tibble(cluster = module_list[[1]], gene = rownames(M))
  sets <- purrr::map(setdiff(unique(ll$cluster), "grey"), \(clust) dplyr::filter(ll, cluster == clust) |> dplyr::pull(gene)) |>
    purrr::set_names(setdiff(unique(ll$cluster), "grey"))



  bg_ensgid <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = colnames(M), keytype = c("ALIAS"), column = "ENSEMBL")
  gs <- readr::read_rds("inst/extdata/go-terms.rds")

  all <- purrr::map(sets, \(set) {
    gsa(
      geneset=AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = set, keytype = c("ALIAS"), column = "ENSEMBL"),
      pathways = gs,
      backg = bg_ensgid,
    ) |>
      dplyr::filter(p_hyper < (0.05 / n_sets_tested))

  })






})
