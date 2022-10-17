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

calculate_gamma_hat <- function(web, rho, delta) {
  f_eig <- function(gamma_avg, web, rho, delta) {
    alpha <- construct_interaction_matrix(web, gamma_avg, rho, delta)
    out <- (min(Re(eigen(alpha)$values)))^2
    out
  }
  out <- optimize(f_eig, c(0, 1000), web = web, rho = rho, delta = delta)$minimum
  return(out)
}

construct_interaction_matrix <- function(web, gamma_avg, rho, delta) {
  SA <- nrow(web)
  SP <- ncol(web)
  alphaA <- matrix(rho, SA, SA) + (1 - rho) * diag(rep(1, SA))
  alphaP <- matrix(rho, SP, SP) + (1 - rho) * diag(rep(1, SP))
  gammaA <- diag(rowSums(web)^-delta) %*% web
  gammaP <- diag(colSums(web)^-delta) %*% t(web)
  f <- sum(gammaA[web == 1] + gammaP[t(web) == 1]) / (2 * sum(web == 1))
  gammaA <- gamma_avg / f * diag(rowSums(web)^-delta) %*% web
  gammaP <- gamma_avg / f * diag(colSums(web)^-delta) %*% t(web)
  alpha <- rbind(cbind(alphaA, -gammaA), cbind(-gammaP, alphaP))
  out <- list(alpha = alpha, alphaA = alphaA, alphaP = alphaP, gammaA = gammaA, gammaP = gammaP)

  -out$alpha
}

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

plot_matrix <- function(x) {
  corrplot(x, col = col3(200), cl.lim = c(-1, 1))
}


load(here("motif_list.RData"))
motif_identification <- motifs %>%
  map_dfr(~ tibble(nrow = nrow(.), ncol = ncol(.), sum = sum(.), product = prod(diag(. %*% t(.))), product2 = sum(. %*% t(.)))) %>%
  mutate(motif_label = row_number())
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

construct_random_interaction_component <- function(web, strength = 1, interaction_type) {
  SA <- nrow(web)
  SP <- ncol(web)
  alphaA <- web %*% t(web)
  alphaA <- scale(alphaA, center = FALSE, scale = colSums(alphaA))
  diag(alphaA) <- 1
  alphaA <- -alphaA
  alphaP <- t(web) %*% web
  alphaP <- scale(alphaP, center = FALSE, scale = colSums(alphaP))
  diag(alphaP) <- 1
  alphaP <- -alphaP
  if (interaction_type == "mutualistic") {
    gammaA <- matrix(rlnorm(prod(dim(web)), 0, strength), nrow = nrow(web)) * web
    gammaP <- t(matrix(rlnorm(prod(dim(web)), 0, strength), nrow = nrow(web)) * web)
  }
  if (interaction_type == "antagonistic") {
    gammaA <- matrix(rlnorm(prod(dim(web)), 0, strength), nrow = nrow(web)) * web
    gammaP <- -t(matrix(rlnorm(prod(dim(web)), 0, strength), nrow = nrow(web)) * web)
  }

  list(alphaA, gammaA, gammaP, alphaP)
}


check_motif_feasibility <- function(alpha, r, label) {
  abundance_equilibrium <- solve(alpha, -r[label])
  feasibility <- if_else(sum(abundance_equilibrium < 0) == 0, "feasible", "infeasible")
  tibble(
    feasibility_motif = feasibility,
  )
}

require(mvtnorm)
calculate_omega <- function(alpha) {
  # S <- nrow(alpha)
  # # Sigma <- solve(t(alpha) %*% alpha)
  # omega <- function(S, Sigma) {
  #   m <- matrix(0, S, 1)
  #   a <- matrix(0, S, 1)
  #   b <- matrix(Inf, S, 1)
  #   d <- tryCatch({
  #     pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
  #     },
  #     error = function(e){
  #       0
  #     })
  #   d[1]^(1 / S)
  # }
  #
  # f <- function(m) class(try(solve(t(m) %*% m), silent = T)) == "matrix"
  # if (f(alpha) == FALSE) {
  #   return(0)
  # }
  # else {
  #   Sigma <- solve(t(alpha) %*% alpha)
  #   return(omega(S, Sigma))
  # }

  S <- nrow(alpha)
  Sigma <- tryCatch(
    {solve(t(alpha) %*% alpha)},
    error = function(e){
      matrix(0, nrow = S, ncol = S)
    }
  )
  m <- matrix(0, S, 1)
  a <- matrix(0, S, 1)
  b <- matrix(Inf, S, 1)
  d <- tryCatch(
    {
      pmvnorm(lower = rep(0, S), upper = rep(Inf, S), mean = rep(0, S), sigma = Sigma)
    },
    error = function(e) {
      0
    }
  )
  d[1]^(1 / S)
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