#
#
# # stat_bayes <- expand_grid(
# #   Nsample = 1:20,
# #   rand_label = 1:100,
# #   persistence = c(T, F),
# #   Nspecies = c(1:5)
# # ) %>%
# #   mutate(persistence_profile = map(
# #     persistence,
# #     function(persistence) {
# #       if (persistence) {
# #         rbinom(1000, 1, alpha)
# #       } else {
# #         rbinom(1000, 1, beta)
# #       }
# #     }
# #   )) %>%
# #   mutate(sampled = map2(
# #     Nsample, persistence_profile,
# #     ~ sample(.y, size = .x)
# #   )) %>%
# #   mutate(
# #     bayes = map_dbl(sampled, function(sampled) {
# #       (alpha / beta)^sum(sampled) *
# #         ((1 - alpha) / (1 - beta))^(sum(1 - sampled))
# #     })
# #   ) %>%
# #   mutate(bayes = map2_dbl(
# #     bayes, Nspecies,
# #     ~ .x / (2^.y)
# #   )) %>%
# #   select(Nsample, persistence, bayes, Nspecies) %>%
# #   mutate(bayes = ifelse(persistence, bayes, 1 / bayes))
# #
# # stat_bayes %>%
# #   mutate(Nspecies = as.factor(Nspecies)) %>%
# #   group_by(Nsample, persistence, Nspecies) %>%
# #   summarise(percent = sum(bayes > 1000) / n()) %>%
# #   ungroup() %>%
# #   mutate(persistence = ifelse(persistence,
# #     "Persistent whole network",
# #     "Non-persistent whole network"
# #   )) %>%
# #   ggplot(aes(Nsample, percent, color = Nspecies)) +
# #   geom_point() +
# #   geom_line() +
# #   facet_wrap(~persistence) +
# #   jtools::theme_nice()
#
# # stat_bayes %>%
# #   mutate(persistence = ifelse(persistence,
# #     "Persistent whole network",
# #     "Non-persistent whole network"
# #   )) %>%
# #   mutate(Nspecies = as.factor(Nspecies)) %>%
# #   filter(Nsample < 10) %>%
# #   ggplot(aes(Nsample, bayes, color = Nspecies)) +
# #   geom_boxplot(aes(group = interaction(Nsample, Nspecies))) +
# #   scale_y_log10(
# #     breaks = scales::trans_breaks("log10", function(x) 10^x),
# #     labels = scales::trans_format("log10", scales::math_format(10^.x))
# #   ) +
# #   geom_hline(yintercept = 100) +
# #   facet_grid(persistence ~ Nspecies, scales = "free") +
# #   jtools::theme_nice()
#
# # stats_bayes <-
# # expand_grid(
# #   Nspecies = c(1:5) * 10,
# #   rand_id = 1:1000
# # ) %>%
# #   mutate(alpha = map_dbl(Nspecies, function(n) {
# #     mat <- runif(n^2, 0, 1 / sqrt(n)) %>%
# #       matrix(nrow = n)
# #     diag(mat) <- -1
# #     feasoverlap::calculate_omega(mat[1:3, 1:3], method = "sphere")
# #   })) %>%
# #   group_by(Nspecies) %>%
# #   summarise(alpha = mean(alpha), .groups = 'drop') %>%
# #   mutate(beta = .01) %>%
# #   expand_grid(
# #     rand_label = 1:100,
# #     Nsample = 1:10,
# #     persistence = c(T, F),
# #   ) %>%
# #   mutate(persistence_profile = pmap(
# #     list(persistence, alpha, beta),
# #     function(persistence, alpha, beta) {
# #       if (persistence) {
# #         data <- rbinom(1000, 1, alpha)
# #       } else {
# #         data <- rbinom(1000, 1, beta)
# #       }
# #       data
# #     }
# #   )) %>%
# #   mutate(sampled = map2(
# #     Nsample, persistence_profile,
# #     ~ sample(.y, size = .x)
# #   )) %>%
# #   mutate(
# #     bayes = pmap_dbl(
# #       list(sampled, alpha, beta),
# #       function(sampled, alpha, beta) {
# #         (alpha / beta)^sum(sampled) *
# #           ((1 - alpha) / (1 - beta))^(sum(1 - sampled))
# #       }
# #     )
# #   ) %>%
# #   select(Nsample, persistence, bayes, Nspecies, rand_label) %>%
# #   mutate(bayes = ifelse(persistence, bayes, 1 / bayes))
#
# # stats_bayes %>%
# #   mutate(Nspecies = as.factor(Nspecies)) %>%
# #   group_by(Nsample, persistence, Nspecies) %>%
# #   summarise(percent = sum(bayes > 100) / n()) %>%
# #   ungroup() %>%
# #   mutate(persistence = ifelse(persistence,
# #                               "Persistent whole network",
# #                               "Non-persistent whole network"
# #   )) %>%
# #   ggplot(aes(Nsample, percent, color = Nspecies)) +
# #   geom_point() +
# #   geom_line() +
# #   facet_wrap(~persistence) +
# #   jtools::theme_nice()
# #
# # stats_bayes %>%
# #   mutate(persistence = ifelse(persistence,
# #                               "Persistent whole network",
# #                               "Non-persistent whole network"
# #   )) %>%
# #   filter(Nsample < 5) %>%
# #   ggplot(aes(Nsample, bayes, color = Nspecies,
# #              group = interaction(Nspecies, Nsample))) +
# #   geom_boxplot() +
# #   scale_y_log10(
# #     breaks = scales::trans_breaks("log10", function(x) 10^x),
# #     labels = scales::trans_format("log10", scales::math_format(10^.x))
# #   ) +
# #   facet_wrap(~persistence) +
# #   jtools::theme_nice()
#
# # alpha <- .5
# # beta <- .1
# # expand_grid(
# #   # Nspecies = c(1:5)*10,
# #   Nsample = 1:10,
# #   persistence = c(T, F)
# # ) %>%
# #   mutate(prob = ifelse(persistence, alpha, beta)) %>%
# #   mutate(bayes = map2_dbl(Nsample, prob, function(Nsample, prob){
# #     Nsample * prob * log(alpha / beta) +
# #       Nsample * (1-prob) * log((1 - alpha) / (1 - beta))
# #   })) %>%
# #   mutate(bayes_lower = map2_dbl(Nsample, prob, function(Nsample, prob){
# #     (Nsample * prob - 2*sqrt(Nsample*prob*(1-prob))) *
# #       log(alpha / beta, base = 10) +
# #       (Nsample * (1-prob) + 2*sqrt(Nsample*prob*(1-prob))) *
# #       log((1 - alpha) / (1 - beta), base = 10)
# #   })) %>%
# #   mutate(bayes_higher = map2_dbl(Nsample, prob, function(Nsample, prob){
# #     (Nsample * prob + 2*sqrt(Nsample*prob*(1-prob))) *
# #       log(alpha / beta, base = 10) +
# #       (Nsample * (1-prob) - 2*sqrt(Nsample*prob*(1-prob))) *
# #       log((1 - alpha) / (1 - beta), base = 10)
# #   })) %>%
# #   mutate(bayes = ifelse(persistence, bayes, -bayes)) %>%
# #   mutate(bayes_lower = ifelse(persistence, bayes_lower, -bayes_lower)) %>%
# #   mutate(bayes_higher = ifelse(persistence, bayes_higher, -bayes_higher)) %>%
# #   mutate(persistence = ifelse(persistence,
# #                               "Persistent whole network",
# #                               "Non-persistent whole network"
# #   )) %>%
# #   ggplot(aes(Nsample, bayes)) +
# #   geom_point() +
# #   geom_errorbar(aes(ymin = bayes_lower, ymax = bayes_higher)) +
# #   facet_wrap(~persistence) +
# #   jtools::theme_nice()
#
#
# # alpha <- .3
# base = exp(1)
# beta <- .01
# result<- expand_grid(
#   Nspecies = c(30, 50, 70),
#   Nsample = 1:100,
#   persistence = c(T, F),
#   alpha = c(.5)
# ) %>%
#   mutate(prob = ifelse(persistence, alpha, beta)) %>%
#   mutate(bayes = pmap_dbl(list(Nsample, prob, alpha),
#                           function(Nsample, prob, alpha){
#     Nsample * prob * log(alpha / beta, base = base) +
#       Nsample * (1-prob) * log((1 - alpha) / (1 - beta), base = base)
#   })) %>%
#   mutate(bayes_sd = pmap_dbl(list(Nsample, prob, alpha),
#                               function(Nsample, prob, alpha){
#     sqrt(Nsample * prob * (1-prob)) *
#       (log(alpha / beta, base = base) -
#          log((1 - alpha) / (1 - beta), base = base))
#   })) %>%
#   mutate(bayes = ifelse(persistence, bayes, -bayes)) %>%
#   mutate(
#     bayes = base^bayes,
#     bayes_sd = base^bayes_sd
#   ) %>%
#   mutate(bayes = ifelse(persistence, 1/Nspecies * bayes, Nspecies*bayes),
#          bayes_sd = ifelse(persistence, 1/Nspecies * bayes_sd, Nspecies*bayes_sd)) %>%
#   mutate(persistence = ifelse(persistence,
#                               "Persistent\nwhole network",
#                               "Non-persistent\nwhole network"
#   ))
#
# ggthemr::ggthemr(palette = 'fresh', layout = 'clean')
# result %>%
#   # filter(Nsample<=10) %>%
#   filter(Nsample<=Nspecies/3) %>%
#   expand_grid(
#     evidence = c('Decisive', 'Strong',
#                  'Substantial', 'not worth')
#   ) %>%
#   mutate(percent = case_when(
#     evidence == 'Decisive' ~ pnorm((bayes - 10^2)/bayes_sd),
#     evidence == 'Strong' ~ pnorm((bayes - 10)/bayes_sd) -pnorm((bayes - 100)/bayes_sd) ,
#     evidence == 'Substantial' ~ pnorm((bayes - 10^(.5))/bayes_sd) - pnorm((bayes - 10)/bayes_sd),
#     evidence == 'not worth' ~ 1-pnorm((bayes - 10^(.5))/bayes_sd)
#   )) %>%
#   mutate(evidence = ordered(evidence,
#                             levels = c('Decisive', 'Strong',
#                                        'Substantial', 'not worth'))) %>%
#   mutate(Nsample = Nsample*3/Nspecies) %>%
#   mutate(Nspecies = paste0(Nspecies, " species")) %>%
#   filter(Nsample <= .6) %>%
#   ggplot(aes(Nsample, percent,
#              group = evidence,
#              fill = evidence)) +
#   geom_area(alpha=0.6 , size=1) +
#   facet_grid(persistence~Nspecies) +
#   scale_x_continuous(breaks= scales::pretty_breaks()) +
#   # jtools::theme_nice() +
#   labs(
#     x = 'Percentage of monitored species',
#     y = 'Proportion'
#   ) +
#   guides(fill = guide_legend(reverse = TRUE)) +
#   theme(aspect.ratio = 1,
#         legend.position = 'bottom',
#         legend.title = element_blank())
# ggsave('Fig4.pdf', width = 6, height = 4.5, dpi = 300)
# # result %>%
# #   filter(Nsample<=15) %>%
# #   expand_grid(
# #     evidence = c('Decisive', 'Strong',
# #                  'Substantial', 'not worth')
# #   ) %>%
# #   mutate(percent = case_when(
# #     evidence == 'Decisive' ~ pnorm((bayes - 2)/bayes_sd),
# #     evidence == 'Strong' ~ pnorm((bayes - 1)/bayes_sd) -pnorm((bayes - 2)/bayes_sd) ,
# #     evidence == 'Substantial' ~ pnorm((bayes - .5)/bayes_sd) - pnorm((bayes - 1)/bayes_sd),
# #     evidence == 'not worth' ~ 1-pnorm((bayes - .5)/bayes_sd)
# #   )) %>%
# #   mutate(evidence = ordered(evidence,
# #                             levels = c('Decisive', 'Strong',
# #                                        'Substantial', 'not worth'))) %>%
# #   ggplot(aes(Nsample, percent,
# #              group = evidence,
# #              fill = evidence)) +
# #   geom_area(alpha=0.6 , size=1) +
# #   facet_grid(~persistence) +
# #   scale_x_continuous(breaks= scales::pretty_breaks()) +
# #   # jtools::theme_nice() +
# #   labs(
# #     x = 'Number of sampled subnetworks',
# #     y = 'Proportion'
# #   ) +
# #   guides(fill = guide_legend(reverse = TRUE)) +
# #   theme(aspect.ratio = 1.5,
# #         legend.position = 'bottom',
# #         legend.title = element_blank())
# # ggsave('Fig3C.pdf', width = 6, height = 5, dpi = 300)
# # ggsave('benno2.jpg', width = 10, height = 5, dpi = 300)
#
# # result %>%
# #   ggplot(aes(Nsample, bayes)) +
# #   geom_line() +
# #   geom_errorbar(aes(ymin = bayes - bayes_sd,
# #                     ymax = bayes + bayes_sd)) +
# #   facet_wrap(~persistence) +
# #   jtools::theme_nice() +
# #   theme(
# #     legend.position = 'bottom'
# #   )
# # ggsave('benno1.jpg', width = 6, height = 3, dpi = 300)
#
# ggthemr::ggthemr(palette = 'fresh', layout = 'clean')
# result %>%
#   # filter(Nsample<=10) %>%
#   filter(Nsample<=Nspecies/3) %>%
#   mutate(Nsample = Nsample*3/Nspecies) %>%
#   mutate(Nspecies = paste0(Nspecies, " species")) %>%
#   filter(Nsample <= .6) %>%
#   ggplot(aes(Nsample, bayes)) +
#   geom_point(alpha=0.6 , size=1) +
#   geom_errorbar(aes(ymin = bayes - 2*bayes_sd, ymax = bayes + 2*bayes_sd)) +
#   facet_grid(persistence~Nspecies, scales = 'free') +
#   # jtools::theme_nice() +
#   labs(
#     x = 'Percentage of monitored species',
#     y = 'Proportion'
#   ) +
#   # scale_x_log10(
#   #   breaks = scales::trans_breaks("log10", function(x) 10^x),
#   #   labels = scales::trans_format("log10", scales::math_format(10^.x))
#   # ) +
#   scale_y_log10(
#     breaks = scales::trans_breaks("log10", function(x) 10^x),
#     labels = scales::trans_format("log10", scales::math_format(10^.x))
#   ) +
#   guides(fill = guide_legend(reverse = TRUE)) +
#   theme(aspect.ratio = 1,
#         legend.position = 'bottom',
#         legend.title = element_blank())
# ggsave('Fig3C_alternatice.pdf', width = 6, height = 4.5, dpi = 300)
