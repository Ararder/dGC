#' @importFrom dplyr .data
utils::globalVariables(c(
  "ensgid", "overlap_size","p", "geneset_size", "background_size",
  "pathway_size","n","p_hyper", ":="
))







#' gene-set enrichment analysis
#'
#' @param geneset a character vector of gene ids ENSGID
#' @param pathways a data.frame with columns ensgid, go, pathway
#' @param backg a character vector of gene ids (ENSGID)
#' @param set_name a character string to name the gene set
#' @param overlap_threshold an integer to filter pathways with overlap size greater than or equal to this value
#'
#' @return a [dplyr::tibble()]
#' @export
#'
#' @examples \dontrun{
#' gsa(geneset, pathways, backg)
#' }
gsa <- function(geneset, pathways, backg, set_name = "set1", overlap_threshold = 10) {
  rlang::check_required(geneset)
  stopifnot(rlang::is_character(geneset))
  rlang::check_required(backg)
  stopifnot(rlang::is_character(backg))
  stopifnot(rlang::is_scalar_character(set_name))
  stopifnot(rlang::is_scalar_integerish(overlap_threshold))

  # filter pathways and gs to only include those in background
  sub_pathways <- dplyr::filter(pathways, ensgid %in% backg)
  sub_gs <- intersect(geneset, backg)


  n_bg <- length(backg)
  n_gs <- length(geneset)
  n_gs_after_filter <- length(sub_gs)
  n_pathways <- nrow(sub_pathways)
  n_sub_pathways <- nrow(sub_pathways)




  # compute overlap
  dplyr::mutate(sub_pathways, {{ set_name}} := dplyr::if_else(ensgid %in% sub_gs, 1,0)) |>
    #
    dplyr::summarise(
      background_size = length(backg),
      geneset_size = length(sub_gs),
      pathway_size  = dplyr::n(),
      overlap_size  = sum(.data[[set_name]]),
      .by = c("go", "pathway")
    ) |>
    dplyr::filter(overlap_size >= overlap_threshold) |>
    dplyr::mutate(
      p_hyper = stats::phyper(overlap_size - 1, geneset_size, background_size - geneset_size, pathway_size, lower.tail=FALSE),
      set_name =  {{ set_name }},
      n_sets_tested = dplyr::n(),
      log_fold_enrich = calc_log_fold_enrichment(background_size, geneset_size, pathway_size, overlap_size),
      odds_ratio = calc_odds_ratio(background_size, geneset_size, pathway_size, overlap_size),
      relative_risk = calc_relative_risk(background_size, geneset_size, pathway_size, overlap_size)

    ) |>
    dplyr::arrange(p_hyper)


}


calc_odds_ratio <- function(background_size, geneset_size, pathway_size, overlap_size) {
  a <- overlap_size
  b <- geneset_size - overlap_size
  c <- pathway_size - overlap_size
  d <- background_size - geneset_size - c

  OR <- (a / b) / (c / d)
  return(OR)
}



# Log-Fold Enrichment (LFE)
calc_log_fold_enrichment <- function(background_size, geneset_size, pathway_size, overlap_size) {
  p1 <- overlap_size / geneset_size
  p2 <- pathway_size / background_size

  LFE <- log2(p1 / p2)
  return(LFE)
}


calc_relative_risk <- function(background_size, geneset_size, pathway_size, overlap_size) {
  p1 <- overlap_size / geneset_size
  p2 <- (pathway_size - overlap_size) / (background_size - geneset_size)

  RR <- p1 / p2
  return(RR)
}



# map_to <- function(geneset, current = c("ALIAS","ENSEMBL", "ENTREZID", "SYMBOL"), to = c("ENSEMBL", "ENTREZID", "SYMBOL")) {
#   rlang::check_required(geneset)
#   rlang::check_installed("org.Hs.eg.db")
#   rlang::check_installed("AnnotationDbi")
#   current <- rlang::arg_match(current)
#   to <- rlang::arg_match(to)
#
#   AnnotationDbi::mapIds(
#     org.Hs.eg.db::org.Hs.eg.db,
#     keys = stringr::str_remove(geneset, "##"),
#     keytype = current,
#     column = to
#   ) |>
#     unname()
#
# }
