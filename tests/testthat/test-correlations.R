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

  d1 <- dplyr::filter(data$obs, dataset == "sandberg")
  m <- data$count_matrix[d1$cell, ]
  data2 <- list(
    count_matrix = m,
    obs = d1,
    var = data$var
  )



  reduced <- prep_cluster_counts(data2, "beta cells", ct_column = "named_celltype", prop_cells = 0.7)
  # M <- compute_residuals(reduced)
  M <- reduced$log2_matrix

  reduced$obs[["random_split"]] <- sample(c("case","control"), size = 735, replace = TRUE)
  true_diff <- corr_diff(M, reduced$obs, by = "random_split",method = "pearson")


  bootstrap_res <- corr_bootstrap_diff(M, obs_df = reduced$obs, method = "pearson", n_boot = 20)

  M <- apply_threshold(bootstrap_res, true_diff)


  normalise_edges <- TRUE
  if(normalise_edges) {
    M <- M / max(abs(M))
  }



})

test_that("plotting works", {
  skip()

  dissTOM <- WGCNA::TOMdist(as.matrix(M), TOMType = "signed")
  hierTOM <- stats::hclust(stats::as.dist(dissTOM), method = "average")



  cutoff <- c(0.974, 0.981, 0.987, 0.992, 0.992)
  deep_split <- c(FALSE, FALSE, FALSE, FALSE, TRUE)
  module_list <- purrr::map2(cutoff, deep_split, \(cut_off, deep) cutreeDynamic(hierTOM, method = "tree", minClusterSize = 30, cutHeight = cut_off, deepSplit = deep) |> labels2colors())


  plotDendroAndColors(
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

  dplyr::tibble(cluster = module_list[[1]], gene = rownames(M)) |>
    dplyr::filter(cluster == "yellow") |> dplyr::pull(gene)





})
