


sample_donors <- function(obj, condition_col="status", donor_col="donor",n=3, replace = TRUE) {


  # check that only two levels exist
  n_levels <- obj[["obs"]][[condition_col]] |> unique() |> length()
  stopifnot(n_levels == 2)

  unique_donors <- split(obj[["obs"]], obj[["obs"]][[condition_col]]) |>
    purrr::map(\(x) dplyr::pull(x, .data[[donor_col]])) |>
    purrr::map(unique)



  A = sample(unique_donors[[1]], n, replace = replace)
  B = sample(unique_donors[[2]], n, replace = replace)



  sub_A <- dplyr::select(dplyr::filter(obj[["obs"]], .data[[donor_col]] %in% A) ,dplyr::all_of(c("cell", donor_col)))
  sub_B <- dplyr::select(dplyr::filter(obj[["obs"]], .data[[donor_col]] %in% B) ,dplyr::all_of(c("cell", donor_col)))

  list(
    A = sub_A,
    B = sub_B
  ) |>
    purrr::set_names(names(unique_donors))
}

