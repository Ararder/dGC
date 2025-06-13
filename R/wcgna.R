# remove utils warning
dummy <- function() {
  utils::apropos()

}


run_wgcna <- function(M) {

  dissTOM <- WGCNA::TOMdist(as.matrix(M), TOMType = "signed")
  hierTOM <- stats::hclust(stats::as.dist(dissTOM), method = "average")
  cutoff <- c(0.974, 0.981, 0.987, 0.992, 0.992)
  deep_split <- c(FALSE, FALSE, FALSE, FALSE, TRUE)
  module_list <- purrr::map2(cutoff, deep_split, \(cut_off, deep) {
    dynamicTreeCut::cutreeDynamic(
      hierTOM,
      method = "tree",
      minClusterSize = 30,
      cutHeight = cut_off,
      deepSplit = deep
    ) |>
      WGCNA::labels2colors()
  })


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


}

