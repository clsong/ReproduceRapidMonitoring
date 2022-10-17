library(deSolve)
library(corrplot)
library(here)
library(combinat)
library(bmotif)
library(patchwork)
library(tidyverse)
library(tidygraph)
library(ggraph)
library(furrr)



construct_pool_matrix <- function(bipartite, mu, sigma) {
  SA <- nrow(bipartite)
  SP <- ncol(bipartite)
  size <- SA + SP

  rho <- 0
  alphaA <- matrix(rho, SA, SA) - (1 - rho) * diag(rep(1, SA))
  alphaP <- matrix(rho, SP, SP) - (1 - rho) * diag(rep(1, SP))

  gammaA <- runif(SA * SP, min = mu / size - sigma * sqrt(3 / size), max = mu / size + sigma * sqrt(3 / size)) %>% matrix(nrow = SA)
  gammaP <- runif(SP * SA, min = mu / size - sigma * sqrt(3 / size), max = mu / size + sigma * sqrt(3 / size)) %>% matrix(nrow = SP)

  gammaA <- gammaA * bipartite
  gammaP <- gammaP * t(bipartite)

  rbind(cbind(alphaA, gammaA), cbind(gammaP, alphaP))
}



calculate_embedded_omega <- function(alpha, alpha_embedded, species) {
  # alpha <- dff$alpha[[1]]
  # alpha_embedded <- dff$alpha_embedded[[1]]
  # species <- dff$species[[1]]
  # alpha
  # alpha_embedded

  # normalize <- function(x){x/sqrt(sum(x^2))}
  # alpha_embedded <- apply(alpha_embedded,2,normalize)

  # tic()
  index <- setdiff(1:ncol(alpha_embedded), species) %>%
    map(~ solve(alpha, alpha_embedded[, .])) %>%
    map_dbl(~ length(which(. < 0))) %>%
    {
      which(. > 0)
    } %>%
    c(species) %>%
    combn(length(species))

  index_max <- 1:ncol(index) %>%
    map_dbl(~ calculate_omega(alpha_embedded[, index[, .]])) %>%
    {
      which(. == max(., na.rm = TRUE))
    } %>%
    {
      index[, .]
    }
  # index_max
  calculate_omega(alpha_embedded[, index_max])
  # toc()

}

extract_motifs <- function(bipartite, Inte, six_node = F, subsampling) {
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

  if(nrow(bipartite) < 6 | ncol(bipartite) < 6){
    res <- res %>%
      filter(num_row + num_col < 5)
  }
    res %>%
      mutate(motif = map2(num_row, num_col, ~ get_all_motifs(.x, .y, bipartite, Inte))) %>%
      unnest(motif) %>%
      mutate(internal_omega = map_dbl(alpha, ~ calculate_omega(.)))

    # mutate(embedded_omega = pmap_dbl(list(alpha, alpha_embedded, species), calculate_embedded_omega))
}


get_all_motifs <- function(num_row, num_col, bipartite, Inte, subsampling = T) {
  numb <- expand.grid(
    combn(1:nrow(bipartite), num_row) %>% split(rep(1:ncol(.), each = nrow(.))),
    combn(1:ncol(bipartite), num_col) %>% split(rep(1:ncol(.), each = nrow(.)))
  ) %>%
    as_tibble()

  # if (subsampling) {
    if (nrow(numb) > 100) {
      numb <- numb %>% sample_n(100)
    }
   #else {
    #   (
    #     numb <- numb %>% sample_frac(.5)
    #   )
    # }
  # }

  numb %>%
    as_tibble() %>%
    mutate(species = map2(Var1, Var2, ~ c(.x, .y + nrow(bipartite)))) %>%
    mutate(alpha = map(species, ~ Inte[., .])) %>%
    mutate(alpha_embedded = map(species, ~ Inte[., ])) %>%
    mutate(motif = map2(Var1, Var2, ~ identify_motif(.x, .y, bipartite))) %>%
    unnest(motif)
}



get_motif_gradient <- function(Inte, r_set, bipartite, subsampling = T, cleaning = T) {
  # bipartite <- df$bipartite[[1]]
  # Inte <- df$Inte[[1]]
  # r_set <- df$r_set[[1]]

  motifs_all <- extract_motifs(bipartite, Inte, subsampling)

  res <- r_set %>%
    mutate(feasibility_whole = map_chr(r, function(r) {
      abundance_equilibrium <- solve(Inte, -r)
      if_else(sum(abundance_equilibrium < 0) == 0, "feasible", "infeasible")
    })) %>%
    mutate(feasibility_motif = map(r, function(r) {
      motifs_all %>%
        mutate(feasible_motif = map2(alpha, species, ~ check_motif_feasibility(.x, r, .y))) %>%
        unnest(feasible_motif)
    })) %>%
    unnest(feasibility_motif)

  if(cleaning == T){
    tryCatch({
      res <- res %>%
        mutate(internal_omega = internal_omega^(num_row + num_col)) %>%
        # mutate(internal_omega = internal_omega) %>%
        mutate(ratio = ifelse(feasibility_motif == "feasible",
                              .5/(internal_omega),
                              .5/(1-internal_omega)
        )) %>%
        select(feasibility_whole, ratio, feasibility_motif) %>%
        group_split(feasibility_whole)
    }, error = function(e) {
      res <- NA
    })
  }
  res
}

quadraticRoots <- function(a, b, c) {

  # print(paste0("You have chosen the quadratic equation ", a, "x^2 + ", b, "x + ", c, "."))

  discriminant <- (b^2) - (4*a*c)

    x_int_plus <- (-b + sqrt(discriminant)) / (2*a)
    # x_int_neg <- (-b - sqrt(discriminant)) / (2*a)

  x_int_plus
}

std <- function(x) sd(x)/sqrt(length(x))


# df_minimum
#
# df_minimum %>%
#   mutate(minimum_points_central = map2_dbl(motif_information, feasibility_whole, function(motif_information, feasibility_whole){
#     if(feasibility_whole != 'feasible'){
#       motif_information$ratio <- 1/motif_information$ratio
#     }
#     x <- motif_information%>%
#       mutate(ratio = log10(ratio)) %>%
#       summarize(
#         mean = mean(ratio),
#         sd = sd(ratio)
#       )
#     ceiling(quadraticRoots(x$mean, -1.96*x$sd, -2)^2)
#   }))  %>%
#   # mutate(id = row_number()) %>%
#   # filter(is.na(minimum_points_central))
#   na.omit() %>%
#   ggplot(aes(minimum_points,minimum_points_central)) +
#   geom_point() +
#   ggpmisc::stat_poly_eq(formula = y~x,
#                aes(label = paste(..eq.label.., ..rr.label.., sep = "~~~")),
#                parse = TRUE) +
#   geom_smooth(method = 'lm') +
#   geom_abline(slope = 1, intercept = 0)
#
# x <- df_minimum$motif_information[[17]] %>%
#   mutate(ratio = log10(ratio)) %>%
#   summarize(
#     mean = mean(ratio),
#     sd = sd(ratio)
#   )
#
# quadraticRoots(x$mean, -1.96*x$sd, -2)
# ggpubr::ggqqplot()
