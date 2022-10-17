#' Calcualte the maximum stable gamma
#' @param web binary network structure
#' @param rho mean field interspecific competition
#' @param delta trade-off
#' @return stability condition
calculate_gamma_hat <- function(web, rho, delta) {
  f_eig <- function(gamma_avg, web, rho, delta) {
    alpha <- construct_interaction_matrix(web, gamma_avg, rho, delta)
    out <- (min(Re(eigen(alpha)$values)))^2
    out
  }

  optimize(f_eig, c(0, 1000), web = web, rho = rho, delta = delta)$minimum
}

#' @param web binary network structure
#' @param gamma_avg average interaction strength
#' @param rho mean-field strength
#' @param delta trade-off
#' @return paramatrized interaction matrix
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

#' Construct a random matrix
#' @param web binary network structure
#' @param strength average interaction strength
#' @param interaction_type mutualistic or antagonistic
#' @return paramatrized interaction matrix
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

