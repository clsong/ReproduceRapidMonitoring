#' Reproduce Figure S11 in the main text
#'
#' @export
plot_figure_S11 <- function() {
  df <- tibble(level1 = names(empirical_data)) %>%
    mutate(level2 = map(level1, ~ names(pluck(empirical_data, .)))) %>%
    unnest(level2) %>%
    mutate(web = map2(level1, level2, ~ pluck(pluck(empirical_data, .x), .y))) %>%
    unnest(web) %>%
    mutate(web = map(web, function(web) {
      web[web > 0] <- 1
      t(web)
    })) %>%
    group_by(level1, level2) %>%
    mutate(level3 = row_number()) %>%
    ungroup() %>%
    select(level1, level2, level3, web) %>%
    mutate(species_number = map_dbl(web, nrow))

  ggthemr(layout = "clean")
  df %>%
    ggplot(aes(level3, species_number)) +
    geom_line(aes(group = level2)) +
    facet_wrap(~level2) +
    labs(
      x = "Month",
      y = "Species richness"
    )
  # ggsave('SI_Fig_empirical_richness.pdf',  width = 6, height = 4)
  df %>%
    rowwise() %>%
    mutate(interaction = sum(web)) %>%
    ungroup() %>%
    ggplot(aes(level3, interaction)) +
    geom_line(aes(group = level2)) +
    facet_wrap(~level2) +
    labs(
      x = "Month",
      y = "# of interactions"
    )
  # ggsave('SI_Fig_empirical_interaction.pdf', width = 6, height = 4)
}

#' Reproduce Figure S12 in the main text
#'
#' @export
plot_figure_S12 <- function() {
  df <- tibble(level1 = names(empirical_data)) %>%
    mutate(level2 = map(level1, ~ names(pluck(empirical_data, .)))) %>%
    unnest(level2) %>%
    mutate(web = map2(level1, level2, ~ pluck(pluck(empirical_data, .x), .y))) %>%
    unnest(web) %>%
    mutate(web = map(web, function(web) {
      web[web > 0] <- 1
      t(web)
    })) %>%
    group_by(level1, level2) %>%
    mutate(level3 = row_number()) %>%
    ungroup() %>%
    select(level1, level2, level3, web) %>%
    mutate(species_number = map_dbl(web, nrow))

  ggthemr(layout = "clean")

  df %>%
    rowwise() %>%
    mutate(interaction = sum(web)) %>%
    ungroup() %>%
    ggplot(aes(level3, interaction)) +
    geom_line(aes(group = level2)) +
    facet_wrap(~level2) +
    labs(
      x = "Month",
      y = "# of interactions"
    )
  # ggsave('SI_Fig_empirical_interaction.pdf', width = 6, height = 4)
}
