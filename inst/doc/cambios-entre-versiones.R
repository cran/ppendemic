## ----setup, include = FALSE---------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>",
  fig.width = 7,
  fig.height = 4.5
)

library(dplyr)
library(tibble)
library(tidyr)

data("ppendemic_tab14", package = "ppendemic")
data("ppendemic_tab15", package = "ppendemic")
data("ppendemic_tab16", package = "ppendemic")

## ----comparison-functions-----------------------------------------------------
candidate_replacements <- function(removed, added) {
  if (nrow(removed) == 0L || nrow(added) == 0L) {
    return(tibble::tibble())
  }

  same_val <- function(x, y) {
    !is.na(x) & !is.na(y) & x != "" & y != "" & x == y
  }

  name_sim <- function(x, y) {
    distance <- diag(adist(tolower(x), tolower(y)))
    lengths <- pmax(nchar(x), nchar(y))
    1 - distance / lengths
  }

  pairs <- dplyr::cross_join(
    removed %>% dplyr::rename_with(~ paste0(.x, "_old")),
    added %>% dplyr::rename_with(~ paste0(.x, "_new"))
  )

  candidates <- pairs %>%
    dplyr::mutate(
      similarity = name_sim(taxon_name_old, taxon_name_new),
      same_species = same_val(Species_old, Species_new),
      same_infraspecies = same_val(infraspecies_old, infraspecies_new),
      same_author = same_val(taxon_authors_old, taxon_authors_new),
      same_family = same_val(family_old, family_new),
      same_year = same_val(year_actual_old, year_actual_new),
      
      score = 0.40 * similarity +
        0.25 * same_species +
        0.10 * same_infraspecies +
        0.15 * same_author +
        0.05 * same_family +
        0.05 * same_year
    ) %>%
    dplyr::filter(score >= 0.50)

  if (nrow(candidates) == 0L) {
    return(tibble::tibble())
  }

  candidates %>%
    dplyr::mutate(
      interpretacion = dplyr::case_when(
        same_species & same_infraspecies & coalesce(infraspecific_rank_old != infraspecific_rank_new, FALSE) ~ "Posible cambio de rango",
        same_species & coalesce(Genus_old != Genus_new, FALSE) ~ "Posible transferencia de genero",
        similarity >= 0.85 & (same_author | same_year) ~ "Posible correccion ortografica",
        TRUE ~ "Posible reemplazo taxonomico"
      ),
      puntuacion = round(score, 3)
    ) %>%
    dplyr::select(
      exclusion_observada = taxon_name_old,
      inclusion_observada = taxon_name_new,
      familia = family_new,
      interpretacion,
      puntuacion
    ) %>%
    dplyr::arrange(dplyr::desc(puntuacion))
}

compare_versions <- function(old, new) {
  removed <- dplyr::anti_join(old, new, by = "taxon_name")
  added <- dplyr::anti_join(new, old, by = "taxon_name")
  candidates <- candidate_replacements(removed, added)

  linked_removed <- unique(candidates$exclusion_observada)
  linked_added <- unique(candidates$inclusion_observada)

  list(
    summary = tibble::tibble(
      version_anterior = unique(old$version),
      version_nueva = unique(new$version),
      registros_anteriores = nrow(old),
      registros_nuevos = nrow(new),
      inclusiones_observadas = nrow(added),
      exclusiones_observadas = nrow(removed),
      cambio_neto = registros_nuevos - registros_anteriores,
      posibles_reemplazos = nrow(candidates)
    ),
    candidates = candidates,
    probable_inclusions = added %>%
      dplyr::filter(!taxon_name %in% linked_added) %>%
      dplyr::select(taxon_name, family, year_actual),
    probable_exclusions = removed %>%
      dplyr::filter(!taxon_name %in% linked_removed) %>%
      dplyr::select(taxon_name, family, year_actual)
  )
}

comparison_14_15 <- compare_versions(ppendemic_tab14, ppendemic_tab15)
comparison_15_16 <- compare_versions(ppendemic_tab15, ppendemic_tab16)

## ----summary-table------------------------------------------------------------
change_summary <- dplyr::bind_rows(
  comparison_14_15$summary,
  comparison_15_16$summary
)

knitr::kable(
  change_summary,
  caption = "Cambios observados y posibles reemplazos entre versiones."
)

## ----change-plot--------------------------------------------------------------
plot_values <- rbind(
  inclusiones = change_summary$inclusiones_observadas,
  exclusiones = -change_summary$exclusiones_observadas,
  cambio_neto = change_summary$cambio_neto
)

barplot(
  plot_values,
  beside = TRUE,
  names.arg = paste(
    change_summary$version_anterior,
    change_summary$version_nueva,
    sep = " a "
  ),
  col = c("#2E8B57", "#B22222", "#4169E1"),
  ylab = "Numero de registros",
  legend.text = rownames(plot_values),
  args.legend = list(x = "topleft", bty = "n")
)
abline(h = 0, col = "grey40")

## ----candidates-14-15---------------------------------------------------------
knitr::kable(
  comparison_14_15$candidates,
  caption = "Posibles reemplazos entre V-14 y V-15."
)

## ----candidates-15-16---------------------------------------------------------
knitr::kable(
  comparison_15_16$candidates,
  caption = "Posibles reemplazos entre V-15 y V-16."
)

## ----inclusion-summary--------------------------------------------------------
inclusions_15_16 <- comparison_15_16$probable_inclusions

inclusion_families <- inclusions_15_16 %>%
  dplyr::count(family, name = "posibles_inclusiones") %>%
  dplyr::arrange(dplyr::desc(posibles_inclusiones)) %>%
  dplyr::slice_head(n = 15) %>%
  dplyr::rename(familia = family)

knitr::kable(
  inclusion_families,
  caption = "Familias con mas posibles inclusiones entre V-15 y V-16."
)

## ----recent-inclusions--------------------------------------------------------
recent_inclusions <- inclusions_15_16 %>%
  dplyr::filter(!is.na(year_actual), year_actual >= 2024) %>%
  dplyr::arrange(dplyr::desc(year_actual)) %>%
  dplyr::slice_head(n = 25)

knitr::kable(
  recent_inclusions,
  caption = "Ejemplos de posibles inclusiones publicadas desde 2024."
)

## ----exclusion-summary--------------------------------------------------------
exclusions_15_16 <- comparison_15_16$probable_exclusions

exclusion_families <- exclusions_15_16 %>%
  dplyr::count(family, name = "posibles_exclusiones") %>%
  dplyr::arrange(dplyr::desc(posibles_exclusiones)) %>%
  dplyr::rename(familia = family)

knitr::kable(
  exclusion_families,
  caption = "Familias de las posibles exclusiones entre V-15 y V-16."
)

## ----exclusion-list-----------------------------------------------------------
knitr::kable(
  exclusions_15_16 %>% dplyr::arrange(family),
  caption = "Posibles exclusiones entre V-15 y V-16."
)

## ----family-deltas------------------------------------------------------------
family_change <- dplyr::full_join(
  ppendemic_tab15 %>% dplyr::count(family, name = "v15"),
  ppendemic_tab16 %>% dplyr::count(family, name = "v16"),
  by = "family"
) %>%
  dplyr::mutate(
    v15 = tidyr::replace_na(v15, 0L),
    v16 = tidyr::replace_na(v16, 0L),
    change = v16 - v15
  ) %>%
  dplyr::arrange(dplyr::desc(abs(change)), family)

knitr::kable(
  head(family_change, 20),
  caption = "Mayores cambios absolutos por familia entre V-15 y V-16."
)

