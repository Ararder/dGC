test_that("reduce matrix works", {
  skip()



  data <- read_data("~/Downloads/qc_merged_seurat.rds")
  reduced <- prep_cluster_counts(
    data,
    prop_cells = 0.5,
    min_cells_per_donor = 10,
    ct = "beta cells",
    ct_column = "named_celltype",
    donor_col = "donor",
    condition_col = "status"
  )


  setup_dgc(obj = reduced, n_iter = 6,fit_models = "blmer",dir = tempdir())
  run_real_diff(tempdir(), ncores = 6)
  mirai::daemons(3)
  run_permutations(tempdir(), ncores = 3)
  list.files(tempdir())






  # -------------------------------------------------------------------------
  tmp <- fs::dir_create(file.path(tempdir(), "dgc"))
  setup_dgc(
    obj = reduced,
    split_by = "status",
    n_iter = 6,
    fit_models = "blmer",
    dir = tmp
  )

  mirai::daemons(2)
  run_permutations(
    tmp,
    ncores =3
  )


  conds <- split(reduced$obs, reduced$obs$status)
  tictoc::tic()
  tt <- corr_diff(
    M = reduced$matrix,
    obs_1 = conds[[1]],
    obs_2 = conds[[2]],
    fit_models = "blmer",
    ncores = 6
  )
  tictoc::toc()
  # This compares both approaches in one table



  compute_residuals(reduced$matrix, cells = reduced$obs$cell, donor_vec = reduced$obs$donor, ncores =4)
  dir <- fs::dir_create(tempdir(), "dgc")
  setup_dgc(
    reduced,
    n_iter = 256,
    split_by ="status",
    fit_models = "blmer",
    dir = dir
  )

  mirai::daemons(2)
  run_permutations(dir)









  generate_random_splits(reduced$obs, 50)


})





test_that("The full pipeline", {
  skip()

  data <- read_data("~/projects/dcgna/workflow/t2d/qc_merged_seurat.rds")
  data$matrix <- (data$count_matrix / rowSums(data$count_matrix))*10^6



  reduced <- prep_cluster_counts(data, "alpha cells", ct_column = "named_celltype", prop_cells = 0.9)

  M <- compute_residuals(reduced$matrix, cells = reduced$obs$cell, donor_vec = reduced$obs$donor)
  reduced$matrix <- M
  reduced <- validate_data(reduced)

  conditions <- split(reduced$obs,reduced$obs$status)


  real_diff <- corr_diff(
    reduced$matrix,
    conditions[[1]],
    conditions[[2]],
    fit_models = "none",
    method = "pearson",
    ncores = 6
  )

  mask_from_perm

  permutations <- corr_permute(reduced, fit_models = "none", n_case = 12, n_ctrl = 16, n_iter = 12, ncores=6)
  # t_mask <- mask_from_perm(P = permutations, R = real_diff)

  mask <- mask_from_perm(P = permutations, R = real_diff, alpha = 0.05)



  idx <- 54
  f_mask <- mask_from_perm(P = permutations[-idx], R = permutations[[idx]])
  fake_diff <- permutations[[idx]]
  fake_diff[f_mask > 0.025] <- 0

  t_mask <- mask_from_perm(P = permutations, R = real_diff)
  real_diff[t_mask > 0.025] <- 0








})

test_that("plotting works", {
  skip()



  noise_minig <- purrr::map(sample(1:length(permutations), 10), \(idx) {
    f_mask <- mask_from_perm(P = permutations[-idx], R = permutations[[idx]])
    fake_diff <- permutations[[idx]]
    fake_diff[f_mask > 0.025] <- 0
    network_clustering(fake_diff)
  })

  noise_minig[[1]]
  nc_res <- network_clustering(real_diff)
  fake_res <- network_clustering(fake_diff)

  purrr::imap(fake_res$genes, \(fake_set, n2) {
    purrr::imap(nc_res$genes, \(real_set, n1) {
      ov <- length(intersect(fake_set, real_set)) / length(union(fake_set, real_set))
      cli::cli_inform("comparing {n2} with {n1}: {ov}")
    })
  })


  real_diff <- readr::read_rds("~/real_diff.rds")
  results <- network_clustering(real_diff)
  network_clustering <- function(real_diff) {
    real_diff <- real_diff / max(abs(real_diff))
    dissTOM <- WGCNA::TOMdist(as.matrix(real_diff), TOMType = "signed")
    hierTOM <- stats::hclust(stats::as.dist(dissTOM), method = "average")
    cutoff <- c(0.974, 0.981, 0.987, 0.992, 0.992)
    deep_split <- c(FALSE, FALSE, FALSE, FALSE, TRUE)
    module_list <- purrr::map2(cutoff, deep_split, \(cut_off, deep) dynamicTreeCut::cutreeDynamic(hierTOM, method = "tree", minClusterSize = 30, cutHeight = cut_off, deepSplit = deep) |> WGCNA::labels2colors())


    plot <- WGCNA::plotDendroAndColors(
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
    ll <- dplyr::tibble(cluster = module_list[[1]], gene = rownames(real_diff))
    sets <- purrr::map(setdiff(unique(ll$cluster), "grey"), \(clust) dplyr::filter(ll, cluster == clust) |> dplyr::pull(gene)) |>
      purrr::set_names(setdiff(unique(ll$cluster), "grey"))


    bg_ensgid <- AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = colnames(real_diff), keytype = c("ALIAS"), column = "ENSEMBL")
    gs <- readr::read_rds("inst/extdata/go-terms.rds")

    all <- purrr::map(sets, \(set) {
      gsa(
        geneset=AnnotationDbi::mapIds(org.Hs.eg.db::org.Hs.eg.db, keys = set, keytype = c("ALIAS"), column = "ENSEMBL"),
        pathways = gs,
        backg = bg_ensgid,
      ) |>
        dplyr::filter(p_hyper < (0.05 / n_sets_tested))

    })

    list(plot, all, genes = sets)
  }












})


test_that("mouse", {
  skip()
  data <- readr::read_rds("~/projects/dcgna/workflow/22q_mouse/layer2_3.rds")
  data$var <- dplyr::tibble(gene = data$var)
  data$matrix <- data$matrix |> Matrix::t()
  data$obs <- data$obs |> dplyr::rename(donor = mouseID, status = genotype) |>
    dplyr::filter(nFeature_RNA > 6000)

  data <- validate_data(data)


  reduced <- prep_cluster_counts(data, "L2/3 IT Stard8", ct_column = "celltype_annotation", prop_cells = 0.99)

  M <- compute_residuals(reduced$matrix, cells = reduced$obs$cell, donor_vec = reduced$obs$donor,ncores = 6)




  conditions <- split(reduced$obs,reduced$obs$status)

  cond <- list(conditions[[2]], conditions[[1]]) |> purrr::set_names("wt","q22")


  real_diff <- corr_diff(
    reduced$matrix,
    conditions[[1]],
    conditions[[2]],
    fit_models = "none",
    method = "pearson",
    ncores = 6
  )


  permuts <- corr_permute(reduced,n_ctrl = 8, n_case=8, ncores =6)

  t_mask <- mask_from_perm(P = permuts, R = real_diff)




})

test_that("parallel", {
  skip()
  # last version, apache arrow parallelisation
  dd <- fs::dir_create(fs::path(tempdir(), "dgc"))
  dplyr::as_tibble(as.matrix(matrix), rownames = "cell") |>
    dplyr::inner_join(dplyr::select(obs,donor,cell)) |>
    dplyr::select(-cell) |>
    arrow::write_parquet(fs::path(dd, "file", ext = "parquet"))

  dns <- obs$donor |> unique()
  dnrs <- dns[1:10]
  ds <- arrow::open_dataset(fs::path(dd, "file", ext = "parquet"))
  genes <- colnames(ds)[-length(colnames(ds))]


  path = fs::path(dd, "file", ext = "parquet")
  mirai::daemons(4)

  res <- mirai::mirai_map(genes, \(x)
                   fit_models_pq(gene =x, path = path, dns = dns, engine = engine),
                   path = path, dns = dns, engine = engine, fit_models_pq  =fit_models_pq
  )[.progress]


  fit_models_pq <- function(gene, dns = dns, engine, path) {
    dl <-
      arrow::read_parquet(path, col_select = dplyr::all_of(c("donor", gene))) |>
      dplyr::filter(donor %in% dns)
    dGC::fit_model(expr = dl[[2]], donor_vec = dl[[1]], engine = engine)

  }




  matrix <- reduced$matrix
  donor_vec <- reduced$obs$donor
  engine = "blmer"


  r <- parallel::mclapply(
    1:ncol(matrix)[1:10],
    \(idx) dGC::fit_model(expr= Matrix::Matrix(matrix)[, idx], donor_vec = donor_vec,engine = engine)
  )



})
