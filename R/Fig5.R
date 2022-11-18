#' Reproduce the data used in Figure 5 in the main text
#'
#' @export
generate_figure_5_data <- function() {
  # clean the data ----------------------------------------------------------

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
    mutate(species_number = map_dbl(web, nrow)) %>%
    expand_grid(
      sample_number = 5:9,
      sample_id = 1:50
    ) %>%
    mutate(web = map2(web, sample_number, ~ .x[sample(1:nrow(.x), .y), ])) %>%
    mutate(motif = map(web, ~ mcount(., six_node = F, normalisation = T, mean_weight = F, standard_dev = F)))

  # randomization -----------------------------------------------------------
  get_random_motif <- function(web, Nrand = 2) {
    from_incidence_matrix_to_bipartite <- function(tmp){
      g <- graph_from_incidence_matrix(tmp)
      cl <- clusters(g, "strong")
      g2 <- induced.subgraph(g, which(cl$membership == which.max(cl$csize)))
      return(g2)
    }

    my.sample <- function(x, size=1) x[sample(length(x), size)]

    # create a random walk through the complete graph
    # add a link every time a new node is discovered
    # this will provide a skeleton for the graph:
    # a (random) spanning tree will connect all of the nodes
    random_spanning_tree <- function(num_rows, num_cols) {
      tree_B <- matrix(0, num_rows, num_cols)
      discovered_rows <- rep(0, num_rows)
      discovered_cols <- rep(0, num_cols)
      my_row <- sample(1:num_rows, 1)
      my_col <- sample(1:num_cols, 1)
      discovered_rows[my_row] <- 1
      discovered_cols[my_col] <- 1
      tree_B[my_row, my_col] <- 1
      while(sum(c(discovered_cols, discovered_rows) == 0) > 0){
        my_row <- sample(1:num_rows, 1)
        if (discovered_rows[my_row] == 0){
          discovered_rows[my_row] <- 1
          tree_B[my_row, my_col] <- 1
        }
        my_col <- sample(1:num_cols, 1)
        if (discovered_cols[my_col] == 0){
          discovered_cols[my_col] <- 1
          tree_B[my_row, my_col] <- 1
        }
      }
      return(tree_B)
    }

    sparse_erdos_renyi <- function(B) {
      ## create a random spanning tree to ensure connectedness
      ER_B <- random_spanning_tree(nrow(B), ncol(B))
      ## now add remaining links to make number of edges match original network
      ER_B <- as.vector(ER_B)
      ER_B[sample(which(ER_B == 0), sum(B) - sum(ER_B), replace=FALSE)] <- 1
      ER_B <- matrix(ER_B, nrow(B), ncol(B))
      dimnames(ER_B) <- dimnames(B)
      ## convert back to igraph object
      return(graph_from_incidence_matrix(ER_B))
    }

    erdos_renyi <- function(g) {
      MAX_ATTEMPTS <- 100

      ## get B from the igraph object
      B <- get.incidence(g)
      ## take some measures to ensure connectedness:
      connected <- FALSE
      n_attempts <- 0
      while (!connected) {
        if (n_attempts > MAX_ATTEMPTS) {
          cat("\nReached", MAX_ATTEMPTS, "ER randomizations without finding a connected one")
          return(NULL)
        }
        n_attempts <- n_attempts + 1
        if (ecount(g) > max(dim(B)) * log(vcount(g))) {
          ## then the ER graph is almost surely connected
          ER_g <- sample_bipartite(nrow(B), ncol(B), m=ecount(g), type="gnm")
        } else {
          ## ensure no empty rows before assigning remaining links
          ER_g <- sparse_erdos_renyi(B)
        }
        ## double check that it is connected
        connected <- is.connected(ER_g)
      }
      #cat("     ER")
      return(ER_g)
    }
    g_connected <- from_incidence_matrix_to_bipartite(web)
    1:Nrand %>%
      map_dfr(~ mutate(
        mcount(get.incidence(erdos_renyi(g_connected)),
          six_node = F,
          # six_node = T,
          normalisation = T, mean_weight = F, standard_dev = F
        ),
        num_rand = .
      ))
  }

  plan(multisession)
  df <- df %>%
    mutate(motif_random = future_map(web,
                                     ~ get_random_motif(., Nrand = 500),
                                     .progress = T))

  # compute the z-score ----------------------------------------------------------------
  scale_this <- function(x) {
    (x - min(x)) / (max(x) - min(x))
  }

  summarize_statistics <- function(motif, motif_random) {
    motif_random <- motif_random %>%
      as_tibble() %>%
      gather(key, value, -motif, -nodes, -num_rand) %>%
      group_by(motif, key) %>%
      summarise(
        mean = mean(value),
        sd = sd(value),
        .groups = "drop"
      )
    motif <- motif %>%
      as_tibble() %>%
      gather(key, value, -motif, -nodes) %>%
      select(-nodes)

    left_join(motif, motif_random, by = c("motif", "key")) %>%
      mutate(z_score = (value - mean) / sd)
  }

  df_statistics <- df %>%
    mutate(motif_statistics = pro_map2(motif, motif_random, ~ summarize_statistics(.x, .y), .progress = T)) %>%
    select(-web, -motif, -motif_random) %>%
    unnest(motif_statistics) %>%
    na.omit()

  # identify a motif with its motif type -----------------------------------------------------------------------
  identify_motif <- function(M) {
    motif_identification %>%
      filter(
        nrow == nrow(M),
        ncol == ncol(M),
        sum == sum(M),
        product == prod(diag(M %*% t(M))),
        product2 == sum(M %*% t(M))
      ) %>%
      pull(motif_type)
  }
  group <- tibble(x = double(), y = double())
  for (i in 1:44) {
    if (!i %in% group$y) {
      group <- group %>%
        bind_rows(tibble(x = i, y = identify_motif(t(motifs[[i]]))))
    }
  }
  group <- group %>%
    mutate(group = paste0(x, "-", y)) %>%
    gather(key, motif_number, -group) %>%
    select(-key) %>%
    mutate(motif_number = as.factor(motif_number)) %>%
    distinct()


  # structural stability --------------------------------------------------------------------
  df_statistics_rank <- df_statistics %>%
    select(level1, level2, level3, motif, z_score, sample_number, sample_id) %>%
    left_join(
      rank %>%
        ungroup() %>%
        select(-group, -value) %>%
        pivot_wider(names_from = key, values_from = rank) %>%
        select(-coef) %>%
        mutate(motif_number = as.numeric(motif_number)),
      by = c("motif" = "motif_number")
    )
}

#' Reproduce Figure 5 in the main text
#'
#' @export
plot_figure_5 <- function(){
 p1 <-
    df_empirical_1 %>%
    mutate(group = factor(group, levels =
                            c("2-3", "4-7", "8-17", "9-13", "10-14", "11-15", "12-16")
    )) %>%
    filter(group == "12-16") %>%
    filter(color_label == "internally more persistent motifs") %>%
    mutate(level1 = if_else(level1 == 'restored', 'Restored regions', 'Unrestored regions')) %>%
    ggplot(aes(level3, z_score, color = level1)) +
    geom_hline(yintercept = 2) +
    geom_hline(yintercept = -2) +
    geom_line(aes(group = interaction(level2, motif_number)), linewidth = .8, alpha = .1) +
    scale_color_manual(values=c('#2AB157', '#A16DAF'), guide = 'none') +
    geom_smooth(aes(group = interaction(level1, motif_number)), se = FALSE) +
    labs(
      y = 'Z scores of subnetworks\nthat persist in isolation',
      x = 'Month',
      title = 'A. Monitoring the whole network'
    ) +
    scale_x_continuous(breaks=1:8) +
    jtools::theme_nice() +
    theme(
      legend.title = element_blank(),
      text = element_text(size=15)
    )

  p2 <- df_empirical_2 %>%
    mutate(group = factor(group, levels =
                            c("2-3", "4-7", "8-17", "9-13", "10-14", "11-15", "12-16")
    )) %>%
    filter(sample_number == 6) %>%
    mutate(level1 = if_else(level1 == 'restored', 'Restored regions', 'Unrestored regions')) %>%
    mutate(sample_number = paste0(sample_number, ' monitored species')) %>%
    filter(color_label == "internally more persistent subnetworks") %>%
    ggplot(aes(level3, z_score, color = level1, fill = level1)) +
    geom_hline(yintercept = 2) +
    scale_color_manual(values=c('#2AB157', '#A16DAF')) +
    scale_fill_manual(values=c('#2AB157', '#A16DAF')) +
    stat_smooth(aes(group = level1), alpha=.3, se=T)+
    scale_linetype_manual(values=c("solid"))+
    labs(
      y = 'Z scores of subnetworks\nthat persist in isolation',
      x = 'Month',
      title = 'B. Monitoring 6 species'
    ) +
    scale_x_continuous(breaks=1:8) +
    jtools::theme_nice() +
    theme(
      legend.position = c(.85, .2),
      legend.title = element_blank(),
      text = element_text(size=15),
      strip.text.x = element_text(size = 10),
      panel.spacing = unit(1, "lines"),
      legend.spacing.y = unit(0, "cm"),
      axis.title.y=element_blank(),
    )


  p1 + p2 +
    plot_layout(guides = 'collect')  &
    theme(legend.position = 'bottom',
          plot.title.position = "plot",
          plot.title = element_text(size=18))
}
