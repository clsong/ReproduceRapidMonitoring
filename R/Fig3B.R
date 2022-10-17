#' Reproduce the data used in Figure 3B in the main text
#'
#' @export
generate_figure3B_data <- function() {
  plan(multisession, workers = 8)

  set.seed(1010)

  df_interaction_matrix <- expand_grid(
    nrow = 12,
    ncol = 12,
    interaction_type = c("mutualistic"),
    n_rep_network = 1:50
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
      Inte_gradient = map(Inte_component, function(Inte_component) {
        tibble(
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
          filter(row_number() == rle(stability)$lengths[1] - 2)
      })
    ) %>%
    unnest(Inte_gradient) %>%
    mutate(r_set = map(Inte, function(Inte) {
      size <- nrow(Inte)
      Nsample <- 100
      # set.seed(1024)
      r_unfeasible <- -Inte %*% runif(size, -1, 0) %>% as.vector()
      r_feasible <- -Inte %*% rep(1, size) %>% as.vector()
      tibble(
        lambda = 1:Nsample / Nsample
      ) %>%
        mutate(r = map(lambda, function(lambda) {
          (1 - lambda) * r_unfeasible + lambda * r_feasible
        }))
    })) %>%
    mutate(motif_information = future_pmap(
      list(Inte, r_set, bipartite),
      ~ get_motif_gradient(..1, ..2, ..3),
      .progress = T
    ))

  fig3B_df_summarized <-
    df_interaction_matrix %>%
    mutate(motif_information = map(motif_information, function(motif_information) {
      motif_information %>%
        select(-num_row, -num_col, -Var1, -Var2, -species) %>%
        mutate(motif_type = case_when(
          feasibility_motif == "feasible" & stability_motif == "stable" ~ "feasible & stable",
          feasibility_motif == "feasible" & stability_motif == "unstable" ~ "infeasible",
          feasibility_motif == "infeasible" ~ "infeasible"
        )) %>%
        group_by(lambda, feasibility_whole, motif_type) %>%
        count() %>%
        ungroup()
    }))
}

#' Reproduce Figure 3B in the main text
#'
#' @export
generate_figure3B_plot <- function() {
  df_plot <- fig3B_df_summarized %>%
    unnest(motif_information) %>%
    group_by(nrow, ncol, n_rep_network, lambda, interaction_type) %>%
    mutate(prop = n / sum(n)) %>%
    ungroup() %>%
    select(-bipartite, -Inte_component, -mutualistic_strength, -stability, -r_set, -Inte) %>%
    nest(data = -c(nrow, ncol, n_rep_network, interaction_type)) %>%
    mutate(critical = map(data, function(data) {
      point <- rle(data$feasibility_whole)$lengths[1]
      mean(c(data$lambda[point], data$lambda[point + 1]))
    })) %>%
    unnest(critical) %>%
    mutate(data = map2(data, critical, function(data, critical) {
      data %>%
        mutate(lambda = ifelse(lambda < critical,
          .5 / critical * lambda,
          .5 / (1 - critical) * lambda + 1 - .5 / (1 - critical)
        ))
    })) %>%
    select(-critical) %>%
    unnest(data) %>%
    mutate(motif_type = if_else(motif_type == "feasible & stable", "coexistent\nsubnetworks", "non-coexistent\nsubnetworks")) %>%
    mutate(lambda = 1 - lambda)

  labelInfo <- split(df_plot, df_plot$motif_type) %>%
    lapply(function(t) {
      data.frame(
        predAtMax = loess(prop ~ lambda, span = 0.8, data = t) %>%
          predict(newdata = data.frame(lambda = max(t$lambda))),
        max = max(t$lambda)
      )
    }) %>%
    bind_rows()
  labelInfo$label <- levels(factor(df_plot$motif_type))

  df_plot %>%
    ggplot(aes(lambda, prop, color = motif_type, group = interaction(n_rep_network, motif_type))) +
    geom_vline(xintercept = .5, size = 1) +
    geom_line(size = .5, alpha = .3) +
    scale_colour_manual(
      values = c("dodgerblue2", "firebrick2")
    ) +
    geom_smooth(aes(group = interaction(motif_type)), se = FALSE, size = 1.5) +
    annotate(geom = "label", x = .85, y = .15, label = "Isolation persistent\nsubnetworks", color = "dodgerblue2", size = 5) +
    annotate(geom = "label", x = .85, y = .85, label = "Isolation non-persistent\nsubnetworks", color = "firebrick2", size = 5) +
    labs(
      y = "Proportion of subnetworks",
      x = "Distance to the centroid"
    ) +
    theme_nice() +
    theme(
      text = element_text(size = 20),
      legend.position = "none"
    )
}
