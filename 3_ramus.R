#### Author: Philip Berg
#### Load packages ####
if (!require(pacman)) {
  install.packages("pacman")
}
library(pacman)
if (!require(limma)) {
  BiocManager::install("limma")
}
if (!require(mavis)) {
  devtools::install_github("PhilipBerg/mavis")
}

p_load(dplyr, tidyr, stringr, purrr, mavis)

#### Define functions ####
baldur_wrapper <- function(data, design, contrast, gam_model, workers) {
  uncertainty <- data %>%
    mavis::estimate_uncertainty('identifier', design, gam_model)
  gam_model %>%
    mavis::estimate_gamma_hyperparameters(data, design) %>%
    baldur::infer_data_and_decision_model(
      id_col = 'identifier',
      design_matrix = design,
      contrast_matrix = contrast,
      uncertainty_matrix = uncertainty,
      stan_model = baldur::empirical_bayes,
      clusters = workers
    )
}

make_pairwise_contrast <- function(x, conditions) {
  m <- matrix(0, conditions, 1)
  m[x] <- c(1, -1)
  return(m)
}

# Function for running t-test
ttest_wrapper <- function(contrast, data) {
  comp <- str_split(contrast, '-')[[1]]
  data %>%
    group_by(identifier) %>%
    mutate(
      comparison = str_flatten(comp, ' vs '),
      p_val = t.test(across(contains(comp[1])), across(contains(comp[2])), var.equal = T)$p.value,
      lfc   = t.test(across(contains(comp[1])), across(contains(comp[2])), var.equal = T)$estimate[1] -
        t.test(across(contains(comp[1])), across(contains(comp[2])), var.equal = T)$estimate[2]
    ) %>%
    group_by(comparison) %>%
    mutate(
      p_val = p.adjust(p_val, 'fdr')
    ) %>%
    select(identifier, comparison, p_val, lfc)
}

# Function for finding alpha used in creating ROC curves for multiple imputation
find_alpha <- function(results, data, regex, cluster) {
  tmp <- results$limma_results[[1]] %>%
    filter(comparison == .data$comparison[[1]])
  p <- sum(str_detect(tmp$identifier, regex))
  n <- nrow(tmp) - p
  long <- results %>%
    select(1, 3) %>%
    unnest(limma_results) %>%
    group_by(imputation, comparison) %>%
    transmute(
      identifier = identifier,
      p_val = p.adjust(p_val, "fdr")
    ) %>%
    ungroup()
  N <- nrow(results)
  multidplyr::cluster_copy(cl,
                           c("p", "n", "regex", "N")
  )
  alpha <- tmp %>%
    use_series(p_val) %>%
    unlist(T, F) %>%
    unique() %>%
    c(., 1)
  out <- list()
  if (length(alpha) > 2000) {
    alpha <- sort(alpha)
    alpha <- alpha[seq.int(1, length(alpha), length.out = 2000)]
  }
  multidplyr::cluster_assign_partition(cl, alpha = alpha)

  for (i in unique(long$comparison)) {
    tmp <- filter(long, comparison == i)
    multidplyr::cluster_copy(cl, c("tmp", "i"))
    out[[i]] <- multidplyr::cluster_call(cl,
                                         map(alpha,
                                             ~ tmp %>%
                                               group_by(identifier) %>%
                                               summarise(
                                                 binom_p_value = binom.test(sum(p_val < .x), N, alternative = "greater")$p.val
                                               ) %>%
                                               mutate(
                                                 tp = str_detect(identifier, regex), tn = !tp
                                               ) %>%
                                               calc_tpr_fpr(.05, ., p, n, "binom_p_value") %>%
                                               mutate(
                                                 comparison = i,
                                                 alpha = .x
                                               ) %>%
                                               bind_rows()
                                         )
    ) %>%
      bind_rows()
  }
  bind_rows(out)
}

# Saving tiff and pdf
ggsave_wrapper <- function(file_name, width, height = NA){
  ggsave(paste0(file_name, '.tiff'), width = width, height = height, units = 'mm', dpi = 320, compression = 'lzw')
  ggsave(paste0(file_name, '.pdf'), width = width, height = height, units = 'mm', dpi = 320)
}


### Performance
calc_tpr_fpr <- function(alpha, hits, p, n, p_val_col){
  p_val <- pull(hits, p_val_col)
  tp <- hits$tp
  tn <- hits$tn
  tibble(
    TP = sum((p_val<=alpha) & tp),
    TN = sum((p_val>=alpha) & tn),
    FP = sum((p_val<=alpha) & tn),
    FN = sum((p_val>=alpha) & tp)
  ) %>%
    mutate(
      TPR = TP/p,
      FPR = FP/n,
      NPV = TN/(TN+FN),
      precision = TP/(TP + FP),
      MCC = sqrt(precision)*sqrt(TPR)*sqrt(1-FPR)*sqrt(NPV) - sqrt(1-precision)*sqrt(1-TPR)*sqrt(FPR)*sqrt(1-NPV)
    )
}

roc_method_within <- function(p_val_col, data, regex, cl) {
  data <- data %>%
    split.data.frame(.$comparison)
  out <- map(data,
             pull, p_val_col
  ) %>%
    map(unique) %>%
    map(~c(-0.1, .)) %>%
    map(enframe, value = 'alpha', name = NULL)
  data <- data %>%
    map(
      mutate,
      tp = str_detect(identifier, regex),
      tn = !tp
    )
  p <- data %>%
    map(
      use_series, tp
    ) %>%
    map(sum)
  n <- data %>%
    map(
      use_series, tn
    ) %>%
    map(sum)
  if(!is.null(cl)){
    multidplyr::cluster_copy(cl,
                             c(
                               "calc_tpr_fpr",
                               "data",
                               "regex",
                               "p",
                               "n"
                             )
    )
  }
  pmap(list(out, data, p, n), create_roc_helper, cl, p_val_col) %>%
    bind_rows()
}

roc_method_between <- function(p_val_col, data, regex, cl) {
  p <- data$identifier %>%
    unique() %>%
    str_detect(regex)
  n <- length(p) - sum(p)
  p <- sum(p)
  if(!is.null(cl)){
    multidplyr::cluster_copy(cl,
                             c(
                               "calc_tpr_fpr",
                               "regex",
                               "p",
                               "n",
                               "between_helper",
                               "p_val_col"
                             )
    )
  }
  data %>%
    nest(data = -comparison) %>%
    multidplyr::partition(cl) %>%
    mutate(
      out = map(data, pull, p_val_col),
      out = map(out, enframe, name = NULL, value = 'alpha'),
      out = map(out, distinct),
      out = map(out, bind_rows, tibble(alpha = c(-.1, 0))),
      data = map(data,
                 ~ transmute(.x,
                             !!p_val_col := !!sym(p_val_col),
                             tp = str_detect(identifier, regex),
                             tn = !tp
                 )
      )
    ) %>%
    mutate(
      out = map2(out, data, between_helper)
    ) %>%
    collect() %>%
    select(comparison, out) %>%
    unnest(cols = out)
}

between_helper <- function(out, data) {
  out %>%
    mutate(
      results = map(alpha,
                    ~ calc_tpr_fpr(.x, data, p, n, p_val_col)
      )
    ) %>%
    unnest(cols = results)
}

create_roc <- function(p_val_col, data, regex, cl = NULL){
  if (!is.null(cl)) {
    multidplyr::cluster_library(cl,
                                c("dplyr",
                                  "stringr",
                                  "tidyr",
                                  "purrr",
                                  "tibble"
                                )
    )
  }
  comps <- data %>%
    use_series(comparison) %>%
    unique() %>%
    length()
  if(comps < 10) {
    roc_method_within( p_val_col, data, regex, cl)
  } else {
    roc_method_between(p_val_col, data, regex, cl)
  }
}

create_roc_helper <- function(out, data, p, n, cl, p_val_col) {
  if(!is.null(cl)){
    out <- out %>%
      multidplyr::partition(cl)
  }
  out <- out %>%
    mutate(
      results = map(alpha, calc_tpr_fpr, data, p, n, p_val_col),
      comparison = data$comparison[[1]]
    )
  if(!is.null(cl)){
    out <- out %>%
      collect()
  }
  out %>%
    unnest(results)
}

integrater <- function(x, y) {
  x_na <- is.na(x)
  y_na <- is.na(y)
  x <- x[!x_na&!y_na]
  y <- y[!x_na&!y_na]
  x_order <- order(x)
  x <- x[x_order]
  y <- y[x_order]
  dx <- diff(x)
  end <- length(y)
  my <- (y[1:(end - 1)] + y[2:end]) / 2
  sum(dx * my)
}

#### Setup variables ####
ramus_design <- model.matrix(~0+factor(rep(1:9, each = 3)))
colnames(ramus_design) <- paste0("condi", 1:9)
ramus_contrast <- combn(colnames(ramus_design), 2) %>%
  t() %>%
  apply(1, str_flatten, '-') %>%
  limma::makeContrasts(contrasts = ., levels = ramus_design)
ramus_cont <- combn(1:9, 2) %>%
  apply(2, make_pairwise_contrast, 9)

full_page <- 178
half_page <- 86
color_theme <- set_names(
  viridisLite::turbo(7, end = .9),
  c(
    'GCR-Baldur',
    'GR-Baldur',
    'GR-Limma',
    'GCR-Limma',
    'Limma-trend',
    't-test',
    "Mavis"
  )
)

ramus_str_replace <- c(
  'condi1' = '0.05',
  'condi2' = '0.125',
  'condi3' = '0.25',
  'condi4' = '0.5',
  'condi5' = '2.5',
  'condi6' = '5',
  'condi7' = '12.5',
  'condi8' = '25',
  'condi9' = '50'
)
ramus_plot_order <- ramus_str_replace %>%
  combn(2) %>%
  t() %>%
  apply(1, str_flatten, ' vs ')


workers <- round(parallel::detectCores()/2)

#### Normalize and partition ####
if (!file.exists('ramus_clean.csv')) {
  readr::read_csv('https://figshare.com/ndownloader/files/35592290?private_link=28e837bfe865e8f13479', show_col_types = FALSE) %>%
    janitor::clean_names() %>%
    mutate(
      across(where(is.numeric), ~na_if(.x, 0))
    ) %>%
    readr::write_csv('ramus_clean.csv')
}

ramus_norm <- readr::read_csv('ramus_clean.csv') %>%
  rename_with(~paste0('condi', rep(1:9, each = 3), '_', rep(1:3, lenght.out = 9*3)), where(is.numeric)) %>%
  psrn(load_info = F, id_col = 'identifier')

ramus_part <- ramus_norm %>%
  single_imputation(ramus_design, workers = workers) %>%
  calculate_mean_sd_trends(ramus_design)

ramus_part <- ramus_norm %>%
  single_imputation(ramus_design, workers = workers) %>%
  calculate_mean_sd_trends(ramus_design) 

ramus_opti <- ramus_part %>%
  grid_search(ramus_design, n_h1 = 20, n_h2 = 20, formula = c(sd ~ mean + c), g = "linear")


ramus_part <- ramus_part %>%
  iterater(ramus_design, sd ~ mean + c, c(ramus_opti$h1[1], ramus_opti$h2[1]), "linear", verbose = 1, lambda = 1, penalty = p, plot = T) %>% 
  use_series(data)

ramus_norm$c <- ramus_part$c

#### Fit regression models ####
gr  <- fit_gamma_regression(ramus_part)
gcr <- fit_gamma_regression(ramus_part, sd ~ mean + c)

#### Run Baldur and t-test ####
ramus_gr_baldur  <- baldur_wrapper(ramus_part, ramus_design, ramus_cont, gr,  workers)
ramus_gcr_baldur <- baldur_wrapper(ramus_part, ramus_design, ramus_cont, gcr, workers)
ramus_ttest      <- map(colnames(ramus_contrast), ttest_wrapper, ramus_part) %>%
  bind_rows()

#### Run limma ####
cl <- multidplyr::new_cluster(workers)
multidplyr::cluster_library(
  cl,
  c(
    'purrr',
    'dplyr',
    'magrittr',
    'mavis',
    "stringr"
  )
)
multidplyr::cluster_copy(cl,
                         "calc_tpr_fpr"
)
ramus_trend <- multiple_imputation_and_limma(
  ramus_norm, ramus_design, ramus_contrast, 1000, workers, "identifier", TRUE,
  FALSE, FALSE
) 

ramus_trend_roc <- ramus_trend %>%
  find_alpha(ramus_norm, "UPS", cl) %>%
  mutate(
    method = "Limma-trend"
  )
ramus_trend <- ramus_trend %>% 
  extract_results(ramus_norm, .05, 0, "fdr", "identifier")

ramus_gr_limma <- multiple_imputation_and_limma(
  ramus_norm, ramus_design, ramus_contrast, 1000, workers, "identifier", TRUE,
  TRUE, FALSE, sd_p ~ mean, sd ~ mean
) 

ramus_gr_limma_roc <- ramus_gr_limma %>%
  find_alpha(ramus_norm, "UPS", cl) %>%
  mutate(
    method = "GR-Limma"
  )
ramus_gr_limma <- ramus_gr_limma %>% 
  extract_results(ramus_norm, .05, 0, "fdr", "identifier")

ramus_gcr_limma <- multiple_imputation_and_limma(
  ramus_norm, ramus_design, ramus_contrast, 1000, workers, "identifier", TRUE,
  TRUE, FALSE, sd_p ~ mean, sd ~ mean*c
) 

ramus_gcr_limma_roc <- ramus_gcr_limma %>%
  find_alpha(ramus_norm, "UPS", cl) %>%
  mutate(
    method = "GCR-Limma"
  )
ramus_gcr_limma <- ramus_gcr_limma %>% 
  extract_results(ramus_norm, .05, 0, "fdr", "identifier")
rm(cl)
gc()

#### Run Mavis ####
stan_mod <- cmdstanr::cmdstan_model("~/Downloads/mavis.stan")

data_regex <- "ramus_.*(limma|baldur|ttest)$"

ramus_all_err <- mget(ls(pattern = data_regex)) %>% 
  map(
    select, 1, comparison, contains('err'), matches('(median_|^)p_val')
  ) %>% 
  imap(~ rename(.x, !!.y := 3)) %>% 
  reduce(left_join, by = c("identifier", "comparison")) %>% 
  split.data.frame(.$comparison) %>% 
  map(select, -comparison)

weights <- mget(ls(pattern = data_regex)) %>%
  map(ungroup) %>%
  map(select, identifier, comparison, any_of(c("err", "p_val", "median_p_val", "median_lfc", "lfc"))) %>%
  map(~ split.data.frame(.x, .x$comparison)) %>% 
  map_depth(2,
            ~ cor(abs(.x[[3]]), abs(.x[[4]]), method = "spear")
  ) %>% 
  map(unlist)

weights <- names(weights[[1]]) %>% 
  setNames(., .) %>% 
  map(
    ~ lapply(weights, `[[`, .x)
  ) %>% 
  map(unlist) %>% 
  map(
    ~ .x/sum(.x)
  )

samp <- ramus_all_err %>% 
  map(make_stan_input) %>% 
  map2(weights, add_weights) %>% 
  map(stan_mod$sample, parallel_chains = 4, iter_warmup = 500, iter_sampling = 500)

mu_smry <- samp %>% 
  map(~ .x$summary("mu"))

post_means <- mu_smry %>% 
  map(use_series, mean) %>% 
  map(multiply_by, -1) %>% 
  map(add, 1)

ramus_mavis <- map2(ramus_all_err, post_means,
                    ~ .x[1] %>% 
                      mutate(
                        p_val = .y
                      )
) %>% 
  imap(~ mutate(.x, comparison = .y)) %>% 
  bind_rows()


#### Performance ####
cl <- multidplyr::new_cluster(workers)
multidplyr::cluster_library(cl,
                            c("dplyr",
                              "stringr",
                              "tidyr",
                              "purrr",
                              "tibble",
                              "magrittr"
                            )
)
multidplyr::cluster_copy(cl,
                         "calc_tpr_fpr"
)

ramus_gr_baldur_roc  <- create_roc("err", ramus_gr_baldur,  "UPS", cl = cl) %>%
  mutate(
    method = 'GR-Baldur'
  )
ramus_gcr_baldur_roc <- create_roc("err", ramus_gcr_baldur, "UPS", cl = cl) %>%
  mutate(
    method = 'GCR-Baldur'
  )
ramus_ttest_roc      <- create_roc("p_val", ramus_ttest,    "UPS", cl = cl) %>%
  mutate(
    method = 't-test'
  )
ramus_mavis_roc <- create_roc('p_val', ramus_mavis, "UPS", cl) %>%
  mutate(
    method = 'Mavis'
  )
rm(cl)
gc()

# Concatenate results from all methods
ramus_roc <- mget(
  ls(pattern = "^ramus.*_roc")
) %>%
  bind_rows() %>%
  mutate(
    comparison = str_replace_all(comparison, ramus_str_replace),
    comparison = factor(comparison, ramus_plot_order)
  )

# Calculate auROC
ramus_auroc <- ramus_roc %>%
  group_by(comparison, method) %>%
  summarise(
    auROC = integrater(FPR, TPR)
  ) %>%
  arrange(desc(auROC)) %>%
  mutate(
    FPR = seq(0.1, .98, length.out = n()),
    TPR = 0
  )

dir.create("roc_data", FALSE)
save(ramus_roc, ramus_auroc, file = "roc_data/ramus_roc.RData")


ramus_roc %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_theme,
                     guide = guide_legend(
                       nrow = 1,
                       title.position = "top",
                       title.vjust = -2.5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0, 1, .5), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0, 1, .5), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(r = 2),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.spacing = unit(0, 'mm')
  ) +
  labs(
    y = 'True Positive Rate',
    x = 'False Positive Rate',
    shape = expression(alpha~phantom())
  ) +
  facet_wrap(. ~ comparison)
ggsave_wrapper('ramus_roc', full_page, full_page)

ramus_roc %>%
  filter(between(alpha, 0, .2)) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, MCC, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_theme,
                     guide = guide_legend(
                       nrow = 1,
                       title.position = "top",
                       title.vjust = -2.5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0, 1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0, 1,.1), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(r = 3, l = 1),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.spacing = unit(0, 'mm')
  ) +
  labs(
    x = expression(alpha),
    y = 'Matthews Correlation Coefficient'
  ) +
  facet_wrap(. ~ comparison)
ggsave_wrapper('ramus_mcc', full_page, full_page)

ramus_roc %>%
  filter(between(alpha, 0, .2)) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, TPR, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_theme,
                     guide = guide_legend(
                       nrow = 1,
                       title.position = "top",
                       title.vjust = -2.5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0, 1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0, 1,.1), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(r = 3, l = 1),
    legend.direction = 'horizontal',
    legend.key.size = unit(.25, "cm"),
    legend.background = element_blank()
  ) +
  labs(
    x = expression(alpha),
    y = 'True Positive Rate'
  ) +
  facet_wrap(. ~ comparison)
ggsave_wrapper('ramus_tpr', width = full_page, height = full_page)

ramus_roc %>%
  filter(between(alpha, 0, .2)) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, FPR, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_theme,
                     guide = guide_legend(
                       nrow = 1,
                       title.position = "top",
                       title.vjust = -2.5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0, 1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0, 1,.1), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(r = 3, l = 1),
    legend.direction = 'horizontal',
    legend.key.size = unit(.25, "cm"),
    legend.background = element_blank()
  ) +
  labs(
    x = expression(alpha),
    y = 'False Positive Rate'
  ) +
  facet_wrap(. ~ comparison)
ggsave_wrapper('ramus_fpr', width = full_page, height = full_page)


ramus_roc %>%
  filter(between(alpha, 0, .2)) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, precision, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_theme,
                     guide = guide_legend(
                       nrow = 1,
                       title.position = "top",
                       title.vjust = -2.5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0, 1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0, 1,.1), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = "bottom",
    plot.margin = margin(r = 3, l = 1),
    legend.direction = 'horizontal',
    legend.key.size = unit(.25, "cm"),
    legend.background = element_blank()
  ) +
  labs(
    x = expression(alpha),
    y = 'Precision'
  ) +
  facet_wrap(. ~ comparison)
ggsave_wrapper('ramus_precision', width = full_page, height = full_page)

