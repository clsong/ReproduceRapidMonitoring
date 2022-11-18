#' Reproduce Figure 4 in the main text
#'
#' @export
plot_figure_4 <- function() {
  base = exp(1)
  beta <- .01
  result<- expand_grid(
    Nspecies = c(30, 50, 70),
    Nsample = 1:100,
    persistence = c(T, F),
    alpha = c(.5)
  ) %>%
    mutate(prob = ifelse(persistence, alpha, beta)) %>%
    mutate(bayes = pmap_dbl(list(Nsample, prob, alpha),
                            function(Nsample, prob, alpha){
                              Nsample * prob * log(alpha / beta, base = base) +
                                Nsample * (1-prob) * log((1 - alpha) / (1 - beta), base = base)
                            })) %>%
    mutate(bayes_sd = pmap_dbl(list(Nsample, prob, alpha),
                               function(Nsample, prob, alpha){
                                 sqrt(Nsample * prob * (1-prob)) *
                                   (log(alpha / beta, base = base) -
                                      log((1 - alpha) / (1 - beta), base = base))
                               })) %>%
    mutate(bayes = ifelse(persistence, bayes, -bayes)) %>%
    mutate(
      bayes = base^bayes,
      bayes_sd = base^bayes_sd
    ) %>%
    mutate(bayes = ifelse(persistence, 1/Nspecies * bayes, Nspecies*bayes),
           bayes_sd = ifelse(persistence, 1/Nspecies * bayes_sd, Nspecies*bayes_sd)) %>%
    mutate(persistence = ifelse(persistence,
                                "Persistent\nwhole network",
                                "Non-persistent\nwhole network"
    ))

  ggthemr(palette = 'fresh', layout = 'clean')
  result %>%
    filter(Nsample<=Nspecies/3) %>%
    expand_grid(
      evidence = c('Decisive', 'Strong',
                   'Substantial', 'not worth')
    ) %>%
    mutate(percent = case_when(
      evidence == 'Decisive' ~ pnorm((bayes - 10^2)/bayes_sd),
      evidence == 'Strong' ~ pnorm((bayes - 10)/bayes_sd) -pnorm((bayes - 100)/bayes_sd) ,
      evidence == 'Substantial' ~ pnorm((bayes - 10^(.5))/bayes_sd) - pnorm((bayes - 10)/bayes_sd),
      evidence == 'not worth' ~ 1-pnorm((bayes - 10^(.5))/bayes_sd)
    )) %>%
    mutate(evidence = ordered(evidence,
                              levels = c('Decisive', 'Strong',
                                         'Substantial', 'not worth'))) %>%
    mutate(Nsample = Nsample*3/Nspecies) %>%
    mutate(Nspecies = paste0(Nspecies, " species")) %>%
    filter(Nsample <= .6) %>%
    ggplot(aes(Nsample, percent,
               group = evidence,
               fill = evidence)) +
    geom_area(alpha=0.6 , size=1) +
    facet_grid(persistence~Nspecies) +
    scale_x_continuous(breaks= scales::pretty_breaks()) +
    labs(
      x = 'Percentage of monitored species',
      y = 'Proportion'
    ) +
    guides(fill = guide_legend(reverse = TRUE)) +
    theme(aspect.ratio = 1,
          legend.position = 'bottom',
          legend.title = element_blank())
}


