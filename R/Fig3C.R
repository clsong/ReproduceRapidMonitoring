library(tidyverse)
library(bmotif)
library(here)
library(furrr)
library(tidymodels)
library(ggrepel)
library(purrrgress)
source("r/toolbox_biomonitoring.r")

theme_set(jtools::theme_nice())

plan(multisession)

# minimum sampled points --------------------------------------------------

library(furrr)

#' Reproduce the data used in Figure 3C in the main text
#'
#' @export
generate_figure3C_data <- function() {
  plan(multisession, workers = 8)

  df_fig3C <- tibble(
    nrow = c(1:4) * 10,
    ncol = 10
  ) %>%
    expand_grid(
      interaction_type = c("antagonistic"),
      n_rep_network = 1:10
    ) %>%
    mutate(
      bipartite = map2(nrow, ncol, function(nrow, ncol) {
        bipartite <- matrix(rbinom(nrow * ncol, size = 1, prob = .5), nrow = nrow)
        while (sum(rowSums(bipartite) == 0) > 0 | sum(colSums(bipartite) == 0) > 0) {
          bipartite <- matrix(rbinom(nrow * ncol, size = 1, prob = .5), nrow = nrow)
        }
        bipartite
      }),
    ) %>%
    mutate(
      Inte_component = map2(
        bipartite, interaction_type,
        ~ construct_random_interaction_component(.x, interaction_type = .y)
      )
    ) %>%
    mutate(
      Inte_gradient = future_map2(Inte_component, nrow, function(Inte_component, nrow) {
        tibble(
          # mutualistic_strength = .1
          mutualistic_strength = seq(from = 0.01, to = 1, by = .01)
        ) %>%
          mutate(Inte = map(mutualistic_strength, function(delta) {
            rbind(
              cbind(Inte_component[[1]], delta * Inte_component[[2]]),
              cbind(delta * Inte_component[[3]], Inte_component[[4]])
            )
          })) %>%
          mutate(stability = map(Inte, function(Inte) {
            try(if_else(sum(Re(eigen(Inte)$values) > 0) == 0, "stable", "unstable"))
          })) %>%
          unnest(stability) %>%
          # filter(row_number() == 1)
          filter(row_number() == rle(stability)$lengths[1] - 2)
      }, .progress = T)
    ) %>%
    unnest(Inte_gradient) %>%
    mutate(r_set = map(Inte, function(Inte) {
      size <- nrow(Inte)
      Nsample <- 100
      r_unfeasible <- -Inte %*% runif(size, -1, 0) %>% as.vector()
      r_feasible <- -Inte %*% rep(1, size) %>% as.vector()
      tibble(
        lambda = 1:Nsample / Nsample
      ) %>%
        filter(lambda %in% c(1 / Nsample, 1)) %>%
        mutate(r = map(lambda, function(lambda) {
          (1 - lambda) * r_unfeasible + lambda * r_feasible
        }))
    })) %>%
    mutate(motif_information = future_pmap(list(Inte, r_set, bipartite),
                                           ~ get_motif_gradient(..1, ..2, ..3),
                                           .progress = T)) %>%
    select(-bipartite, -stability, -Inte_component,
           -mutualistic_strength, -Inte, -r_set)

  fic3C_df
}

#' Reproduce Figure 3C in the main text
#'
#' @export
generate_figure3C_plot <- function(){
  bootstrap_times <- 50

  df_convergence <- fic3C_df %>%
    mutate(feasibility_whole = map_chr(motif_information, ~ unique(.$feasibility_whole))) %>%
    mutate(compensation = ifelse(interaction_type == "mutualistic", .5, .5)) %>%
    mutate(motif_information = map(motif_information, ~ mutate(select(., -feasibility_whole), ratio = log(ratio)))) %>%
    mutate(motif_information = map2(motif_information, compensation, function(motif_information, compensation) {
      motif_information %>%
        mutate(ratio = ifelse(feasibility_motif == "feasible", ratio - log10(.5 / compensation), ratio + log10((1 - compensation) / .5)))
    })) %>%
    mutate(confidence_sample = map(motif_information, function(motif_information) {
      1:bootstrap_times %>%
        map(function(x) {
          motif_information$ratio %>%
            sample() %>%
            cumsum() %>%
            enframe(name = "sampled_number", value = "confidence") %>%
            bind_rows(
              tibble(
                sampled_number = 0,
                confidence = 0
              )
            )
        }) %>%
        bind_rows(.id = "bootstrap_id")
    })) %>%
    select(-motif_information) %>%
    mutate(richness = as.factor(nrow + ncol)) %>%
    unnest(confidence_sample) %>%
    mutate(feasibility_whole = as_factor(feasibility_whole)) %>%
    mutate(
      feasibility_whole = fct_recode(feasibility_whole,
                                     "Persistent whole network" = "feasible",
                                     "Non-persistent whole network" = "infeasible"
      )
    )

  df_stats %>%
    ungroup() %>%
    mutate(n = ifelse(is.na(n), 0, n)) %>%
    filter(decisive) %>%
    filter(sampled_number < 20) %>%
    ggplot(aes(sampled_number, n, group = richness, color = feasibility_whole, alpha = richness)) +
    geom_line(size = 2) +
    facet_wrap(~feasibility_whole) +
    # scale_color_grey() +
    labs(
      x = "Number of sampled subnetworks",
      y = "Proporotion with\ndecisive Bayes factor"
    ) +
    scale_colour_manual(
      values = c("dodgerblue2", "firebrick2"),
      guide = "none"
    ) +
    scale_alpha_discrete(range = c(.1, .5, 0.8, 1)) +
    scale_x_continuous(
      breaks = 0:15,
      limits = c(0, 15)
    ) +
    jtools::theme_nice() +
    theme(
      legend.title = element_blank(),
      legend.position = c(.9, .3),
      # aspect.ratio = 1,
      panel.spacing = unit(2, "lines"),
      text = element_text(size = 14),
      legend.key.size = unit(.8, "cm"),
      legend.text = element_text(size = 17)
    )
}
