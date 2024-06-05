### Author: Philip Berg
library(pacman)
p_load(dplyr, tidyr, stringr, purrr, mavis, ranger, ggplot2, forcats)

calc_lfc <- function(data, condi1, condi2) {
  data %>%
    transmute(
      !!paste0(condi1, "_vs_", condi2) := rowMeans(across(contains(condi1)), na.rm = T) - rowMeans(across(contains(condi2)), na.rm = T)
    )
}

human_str_replace <- c(
  'spike_prop_6' = '1:6',
  'spike_prop_12' = '1:12',
  'spike_prop_25' = '1:25'
)

human_design <- model.matrix(~ 0 + factor(rep(1:3, each = 23)))
colnames(human_design) <- paste0('spike_prop_', c(25, 12, 6))

human <- readxl::read_excel('../mavis-paper-code/diaWorkflowResults_allDilutions.xlsx', na = 'NA') %>%
  select(identifier = 1, everything(), -2:-24) %>%
  filter(!if_all(-identifier, is.na)) %>%
  rename_with(
    ~ str_replace_all(.x, setNames(c('', 'spike_prop_'), c('_.*', '1-'))) %>%
      paste0('_', 1:23),
    -identifier
  ) %>%
  mutate(
    across(where(is.numeric), ~ magrittr::raise_to_power(2, .x))
  ) %>%
  psrn('identifier')


human_filt <- human %>%
  drop_na()

human_mv_filt <- human_filt %>%
  calculate_mean_sd_trends(human_design) %>%
  select(mean, sd)

human_filt <- map2(colnames(human_design)[c(1, 1, 2)], colnames(human_design)[c(2, 3, 3)], ~ calc_lfc(human_filt, .x, .y)) %>%
  bind_cols(human_filt, human_mv_filt, .) %>%
  select(-matches("^spike_prop_[0-9]*_[0-9]*$")) %>%
  mutate(
    class = if_else(str_detect(identifier, "ECOLI"), "TP", "TN") %>%
      as.factor()
  )

human_filt %>%
  pivot_longer(-c(1:3, class)) %>%
  mutate(
    name = stringr::str_replace_all(name, human_str_replace) %>%
      stringr::str_replace_all("_", " "),
    class = stringr::str_replace_all(class, c("TP" = "True Positive", "TN" = "True Negative"))
  ) %>%
  arrange(class) %>%
  ggplot(aes(mean, value, color = class)) +
  geom_point(size = 1/10) +
  scale_color_viridis_d("Class", guide = guide_legend(override.aes = list(size=1)), end = 0, begin = 1) +
  theme_bw() +
  facet_grid(name ~ .) +
  theme(legend.position = c(.9, .75), legend.background = element_blank()) +
  labs(y = expression(Log[2]~Fold~Change), x = expression(bar(X)))
ggsave("no_imp_human.png", dpi = 600, width = 178, height = 200, units = "mm")


human_filt %>%
  count(str_detect(identifier, "ECOLI"))

w <- sum(as.logical(human_filt$class == "TP")) / nrow(human_filt)
rf <- ranger::ranger(class ~ ., human_filt[-1], class.weights = c(w, 1 - w))


human_mv <- human %>%
  calculate_mean_sd_trends(human_design) %>%
  select(mean, sd)

human2 <- map2(colnames(human_design)[c(1, 1, 2)], colnames(human_design)[c(2, 3, 3)], ~ calc_lfc(human, .x, .y)) %>%
  bind_cols(human, human_mv, .) %>%
  select(-matches("^spike_prop_[0-9]*_[0-9]*$")) %>%
  mutate(
    class = if_else(str_detect(identifier, "ECOLI"), "TP", "TN") %>%
      as.factor()
  )

human2 %>%
  pivot_longer(-c(1:3, class)) %>%
  mutate(
    name = stringr::str_replace_all(name, human_str_replace) %>%
      stringr::str_replace_all("_", " "),
    class = stringr::str_replace_all(class, c("TP" = "True Positive", "TN" = "True Negative"))
  ) %>%
  arrange(class) %>%
  ggplot(aes(mean, value, color = class)) +
  geom_point(size = 1/10) +
  scale_color_viridis_d("Class", guide = guide_legend(override.aes = list(size=1)), end = 0, begin = 1) +
  theme_bw() +
  facet_grid(name ~ .) +
  theme(legend.position = c(.9, .75), legend.background = element_blank()) +
  labs(y = expression(Log[2]~Fold~Change), x = expression(bar(X)))
ggsave("no_imp_ignore_na_human.png", dpi = 600, width = 178, height = 200, units = "mm")

human2 <- human2 %>%
  drop_na() %>%
  mutate(
    pred = predict(rf, data = drop_na(human2))$predictions,
    class = if_else(str_detect(identifier, "ECOLI"), "TP", "TN")
  )

human2 %$%
  table(class, pred)

human2 <- human2 %>%
  unite(cls_prd, class, pred, sep = '-') %>%
  select(1, cls_prd)

human2 %>%
  filter(cls_prd != "TP-TN") %>%
  select(1) %>%
  readr::write_csv("human_ids.csv")

human_mv_imp <- human %>%
  calculate_mean_sd_trends(human_design) %>%
  select(mean, sd)

human_imp <- map2(colnames(human_design)[c(1, 1, 2)], colnames(human_design)[c(2, 3, 3)], ~ calc_lfc(human, .x, .y)) %>%
  bind_cols(human, human_mv_imp, .) %>%
  select(-matches("^spike_prop_[0-9]*_[0-9]*$")) %>%
  mutate(
    class = if_else(str_detect(identifier, "ECOLI"), "TP", "TN") %>%
      as.factor()
  )

human_imp %>%
  pivot_longer(-c(1:3, class)) %>%
  mutate(
    name = stringr::str_replace_all(name, human_str_replace) %>%
      stringr::str_replace_all("_", " "),
    class = stringr::str_replace_all(class, c("TP" = "True Positive", "TN" = "True Negative"))
  ) %>%
  arrange(class) %>%
  ggplot(aes(mean, value, color = class)) +
  geom_point(size = 1/10) +
  scale_color_viridis_d("Class", guide = guide_legend(override.aes = list(size=1)), end = 0, begin = 1) +
  theme_bw() +
  facet_grid(name ~ .) +
  theme(legend.position = c(.9, .75), legend.background = element_blank()) +
  labs(y = expression(Log[2]~Fold~Changhuman_impe), x = expression(bar(X)))
ggsave("imp_human.png", dpi = 600, width = 178, height = 200, units = "mm")



human_imp %>%
  select(-class) %>%
  left_join(human2) %>%
  mutate(
    cls_prd = replace_na(cls_prd, "Discarded") %>%
      forcats::fct_infreq()
  ) %>%
  pivot_longer(-c(1:3, cls_prd)) %>%
  mutate(
    name = stringr::str_replace_all(name, human_str_replace) %>%
      stringr::str_replace_all("_", " ")
  ) %>%
  arrange(cls_prd) %>%
  ggplot(aes(mean, value, color = cls_prd)) +
  geom_point(size = 1/10) +
  scale_color_viridis_d("Data-Prediction", guide = guide_legend(override.aes = list(size=1)), end = .1, begin = .9, option = "H") +
  theme_bw() +
  facet_grid(name ~ .) +
  theme(legend.position = c(.9, .78), legend.background = element_blank()) +
  labs(y = expression(Log[2]~Fold~Change), x = expression(bar(X)))
ggsave("imp_human_pred.png", dpi = 600, width = 178, height = 200, units = "mm")
