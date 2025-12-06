







# perm_gene_test <- function(P, R) {
#   G = nrow(R)
#   N = length(P)
#
#
#   # future::plan(future::multisession, workers = ncores)
#   gene_null_distrib <- purrr::map(seq_along(P[1:25]), \(idx){
#     p_iter <- P[[idx]]
#
#     mask_x <- mask_from_perm(
#       P = P[-idx],
#       R = p_iter
#     )
#
#
#     R[mask_x > 0.025] <- 0
#     R[mask_x < 0.025] <- 1
#
#     edges <- rowSums(R)
#
#     as.matrix(edges)
#
#   },.progress = list(type = "tasks"))
#
#
#   # null_dist <- purrr::reduce(gene_null_distrib, cbind)
#   # row_max <- apply(null_dist, 1, max)
#   #
#   # mask <- mask_from_perm(P, R)
#   # R[mask > 0.025] <- 0
#   # R[mask < 0.025] <- 1
#   # edges <- rowSums(R)
#   # obs <- dplyr::tibble(links = edges, genes = names(edges))
#   # null_res <- dplyr::tibble(links = row_max, genes = names(row_max))
#   # dplyr::inner_join(obs, null_res,by = "genes") |>
#   #   dplyr::mutate(diff = links.x - links.y) |>
#   #   dplyr::filter(links.x > links.y) |>
#   #   print(n = 21)
#
# }





