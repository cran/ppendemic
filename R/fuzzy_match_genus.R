#' Fuzzy Match Genus Name
#'
#' @description
#' This function performs a fuzzy match of genus names against the ppendemic database using fuzzyjoin::stringdist() to account for slight variations in spelling.
#'
#' @param df A tibble containing the genus names to be matched.
#' @param target_df A tibble representing the ppendemic database containing the reference list of endemic species.
#'
#' @return
#' A tibble with two additional columns:
#' - fuzzy_match_genus: A logical column indicating whether the genus was successfully matched (`TRUE`) or not (`FALSE`).
#' - fuzzy_genus_dist: A numeric column representing the distance for each match.
#'
#' @keywords internal
fuzzy_match_genus <- function(df, target_df = NULL){
  assertthat::assert_that(all(c('Orig.Genus',
                                'Orig.Species',
                                'Orig.Infraspecies') %in% colnames(df)))

  ## solve issue of empty input tibble, and needed to ensure compatilbility with sequential_matching: because there the columns already exists for the second backbone
  if(nrow(df) == 0){
    if(!all(c('fuzzy_match_genus', 'fuzzy_genus_dist') %in% colnames(df))){
      return(tibble::add_column(df, fuzzy_match_genus = NA, fuzzy_genus_dist = NA))
    }
    else{
      return(df)
    }
  }
  ## solve issue in second iteration of sequential_matching: necessary to remove fuzzy_species_dist column: otherwise 2 columns are generated 'fuzzy_species_dist...1, fuzzy_species_dist...2'
  if('fuzzy_genus_dist' %in% colnames(df)){
    df <- df |>
      dplyr::mutate(fuzzy_genus_dist = NULL)
  } ## TODO: can potentially be removed again????

  Tree.Genera <- target_df |>
    dplyr::distinct(Genus)
  # fuzzy match
  matched_temp <- df |>
    fuzzyjoin::stringdist_left_join(Tree.Genera,
                                    by = c('Orig.Genus' = 'Genus'),
                                    max_dist = 1,
                                    distance_col = 'fuzzy_genus_dist') |>
    # save matched Genus name to Matched.Genus
    dplyr::mutate(Matched.Genus = Genus) |>
    dplyr::select(-c('Genus')) |>
    dplyr::group_by(Orig.Genus, Orig.Species) |>
    dplyr::filter(fuzzy_genus_dist == min(fuzzy_genus_dist))


  ## If there are multiple matches for the same genus: raise warning and advise for manual checking
  if(matched_temp |>
     dplyr::filter(dplyr::n() > 1) |>
     nrow() > 0){
    message("Multiple fuzzy matches for genera with similar string distance:
            Please consider curating the ambiguous entries by hand and re-run the pipeline.
            The ambiguous matched genera were saved to 'treemendous_ambiguous_genera.csv' in the current working directory.
             The algorithm will choose the first genus to continue.")
    #Do you want save a list of the ambiguous matched genera current working directory in 'treemendous_ambiguous_genera.csv'?")
    ## Save ambiguous genera for manual curation:
    matched_temp |>
      dplyr::filter(dplyr::n() > 1) |>
      dplyr::select(Genus.x, Species, Matched.Genus) |>
      readr::write_csv(file = 'ambiguous_genera.csv') ##
    ## Alternative Idea: prompt the user to insert the correct name. Caution here however because this might cause trouble with unit testing
  }

  ## continue selecting first genus if more than one match
  matched <- matched_temp |>
    dplyr::group_modify(
      ~ifelse(nrow(.x) == 0,
              return(.x),
              return(dplyr::slice_head(.x,n = 1)))
      ## In cases of multiple matches: we choose first match.
      ## Alternatively could use something more sophisticated here:
      ## like for instance choosing the one with more support (present
      ## in more databases)
    ) |>
    dplyr::ungroup()


  unmatched <- df |>
    fuzzyjoin::stringdist_anti_join(Tree.Genera,
                                    by = c('Orig.Genus' = 'Genus'),
                                    max_dist = 1)
  #unmatched
  assertthat::assert_that(nrow(df) == (nrow(matched) + nrow(unmatched)))

  res <-  dplyr::bind_rows(matched, unmatched,
                           .id = 'fuzzy_match_genus') |>
    dplyr::mutate(fuzzy_match_genus = (fuzzy_match_genus == 1)) |>  ## convert to Boolean
    dplyr::arrange(Orig.Genus, Orig.Species) |>
    dplyr::relocate(c('Orig.Genus',
                      'Orig.Species',
                      'Orig.Infraspecies')) ## Genus & Species column at the beginning of tibble
  #res
  return(res)
}
