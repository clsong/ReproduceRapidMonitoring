# library(deSolve)
# library(corrplot)
# library(here)
# library(combinat)
# library(bmotif)
# library(patchwork)
# library(tidyverse)
# library(tidygraph)
# library(ggraph)
# library(furrr)

#' identify the type of a motif
#' @param row_label
#' @param col_label
#' @param bipartite
#' @return motif type
identify_motif <- function(row_label, col_label, bipartite) {
  M <- bipartite[row_label, col_label]
  if (!is.matrix(M)) M <- matrix(M, ncol = length(col_label))
  if (sum(nrow(M) == 0) > 0 | sum(ncol(M) == 0) > 0) {
    return(NA)
  }

  motif_identification %>%
    filter(
      nrow == nrow(M),
      ncol == ncol(M),
      sum == sum(M),
      product == prod(diag(M %*% t(M))),
      product2 == sum(M %*% t(M))
    ) %>%
    pull(motif_label)
}

#' check if motif is feasible
#' @param alpha
#' @param r
#' @param label
#' @return motif feasibility
check_motif_feasibility <- function(alpha, r, label) {
  abundance_equilibrium <- solve(alpha, -r[label])
  feasibility <- if_else(sum(abundance_equilibrium < 0) == 0, "feasible", "infeasible")
  if (feasibility == "feasible") {
    stability <- if_else(sum(Re(eigen(alpha %*% diag(abundance_equilibrium))$values) > 0) == 0, "stable", "unstable")
  } else {
    stability <- "unstable"
  }
  tibble(
    feasibility_motif = feasibility,
    stability_motif = stability
  )
}

#' calculate structural stability of an interaction network
#' @param alpha parametrized network structure
#' @return structural stability
calculate_omega <- function(alpha) {
  S <- nrow(alpha)
  Sigma <- solve(t(alpha) %*% alpha)
  m <- matrix(0, S, 1)
  a <- matrix(0, S, 1)
  b <- matrix(Inf, S, 1)
  d <- pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
  d[1]^(1 / S)
}


#' calculate motif profiles with a list of intrinic
#' growth rates
#' @param Inte parametrized network structure
#' @param r_set
#' @param bipartite
#' @return a tibble
get_motif_gradient <- function(Inte, r_set, bipartite) {
  get_all_motifs <- function(num_row, num_col, bipartite, Inte) {
    numb <- expand.grid(
      combn(1:nrow(bipartite), num_row) %>% split(rep(1:ncol(.), each = nrow(.))),
      combn(1:ncol(bipartite), num_col) %>% split(rep(1:ncol(.), each = nrow(.)))
    ) %>%
      as_tibble()

    if (nrow(numb) > 100) {
      numb <- numb %>% sample_n(100)
    } else {
      (
        numb <- numb %>% sample_frac(.5)
      )
    }

    numb %>%
      mutate(species = map2(Var1, Var2, ~ c(.x, .y + nrow(bipartite)))) %>%
      mutate(alpha = map(species, ~ Inte[., .])) %>%
      mutate(motif = map2(Var1, Var2, ~ identify_motif(.x, .y, bipartite))) %>%
      unnest(motif)
  }

  extract_motifs <- function(bipartite, Inte, six_node = F) {
    if (six_node) {
      res <- tibble(
        num_row = c(2, 1, 1, 2, 3, 1, 2, 3, 4, 1, 4, 3, 2, 5),
        num_col = c(1, 2, 3, 2, 1, 4, 3, 2, 1, 5, 2, 3, 4, 1)
      )
    } else {
      res <- tibble(
        num_row = c(2, 1, 1, 2, 3, 1, 2, 3, 4),
        num_col = c(1, 2, 3, 2, 1, 4, 3, 2, 1)
      )
    }

    res %>%
      mutate(motif = map2(num_row, num_col, ~ get_all_motifs(.x, .y, bipartite, Inte))) %>%
      unnest(motif) %>%
      mutate(omega = map(alpha, ~ calculate_omega(.))) %>%
      unnest(omega)
  }

  motifs_all <- extract_motifs(bipartite, Inte)

  res <- r_set %>%
    mutate(feasibility_whole = map(r, function(r) {
      abundance_equilibrium <- solve(Inte, -r)
      feasibility <- if_else(sum(abundance_equilibrium < 0) == 0, "feasible", "infeasible")
      if (feasibility == "feasible") {
        stability <- if_else(sum(Re(eigen(Inte %*% diag(abundance_equilibrium))$values) > 0) == 0, "stable", "unstable")
      } else {
        stability <- "unstable"
      }
      tibble(
        feasibility_whole = feasibility,
        stability_whole = stability
      )
    })) %>%
    unnest(feasibility_whole) %>%
    mutate(feasibility_motif = map(r, function(r) {
      motifs_all %>%
        mutate(feasible_motif = map2(alpha, species, ~ check_motif_feasibility(.x, r, .y))) %>%
        unnest(feasible_motif)
    })) %>%
    unnest(feasibility_motif)
}

get_subnetwork_gradient <- function(Inte, r_set, bipartite) {
  subnetwork_all <-
    r_set %>%
    mutate(feasibility_whole = map(r, function(r) {
      abundance_equilibrium <- solve(Inte, -r)
      feasibility <- if_else(sum(abundance_equilibrium < 0) == 0, "feasible", "infeasible")
      if (feasibility == "feasible") {
        stability <- if_else(sum(Re(eigen(Inte %*% diag(abundance_equilibrium))$values) > 0) == 0, "stable", "unstable")
      } else {
        stability <- "unstable"
      }
      tibble(
        feasibility_whole = feasibility,
        stability_whole = stability
      )
    })) %>%
    unnest(feasibility_whole) %>%
    mutate(feasibility_motif = map(r, function(r) {
      subnetwork_all %>%
        mutate(feasible_motif = map2(alpha, species, ~ check_motif_feasibility(.x, r, .y))) %>%
        unnest(feasible_motif)
    })) %>%
    unnest(feasibility_motif)
}
