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

p_load(dplyr, tidyr, purrr, mavis, magrittr, ggplot2)

#### Define functions ####
make_stan_input <- function(data) {
  N <- nrow(data)
  data <- data %>% 
    select(where(is.numeric)) %>% 
    mutate(
      across(where(is.numeric), ~ if_else(. != 0, ., min(.[. != 0])))
    )
  M <- ncol(data)
  list(
    N = N,
    M = M,
    y = 
      matrix(
        unlist(data), N, M
      )
  )
}

add_weights <- function(stan_input, weights) {
  stan_input$tau <- weights
  stan_input
}

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
  p <- sum(str_detect(results$limma_results[[1]]$identifier, regex))
  n <- nrow(results$limma_results[[1]]) - p
  multidplyr::cluster_copy(cl,
                           c("p", "n", "regex")
  )
  alpha <- results$limma_results[[1]] %>%
    filter(comparison == .data$comparison[[1]]) %>%
    use_series(p_val) %>%
    unlist(T, F) %>%
    unique() %>%
    c(., 1)
  if (length(alpha) > 2000) {
    alpha <- sort(alpha)
    alpha <- alpha[seq.int(1, length(alpha), length.out = 2000)]
  }
  alpha %>%
    tibble::enframe(name = NULL, value = "alpha") %>%
    multidplyr::partition(cluster) %>%
    mutate(
      tmp = map(alpha,
                ~ extract_results(
                  results,
                  data,
                  alpha = .x,
                  abs_lfc = 0,
                  id_col = 'identifier',
                  pcor = 'fdr'
                ) %>%
                  mutate(
                    tp = str_detect(identifier, regex), tn = !tp
                  ) %>%
                  split(.$comparison) %>%
                  imap(
                    ~ calc_tpr_fpr(.05, .x, p, n, "binom_p_value") %>%
                      mutate(
                        comparison = .y
                      )
                  ) %>%
                  bind_rows()
      )
    ) %>%
    collect() %>%
    unnest(tmp)
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
    TP = sum((p_val<alpha) & tp),
    TN = sum((p_val>alpha) & tn),
    FP = sum((p_val<alpha) & tn),
    FN = sum((p_val>alpha) & tp)
  ) %>%
    mutate(
      TPR = TP/p,
      FPR = FP/n,
      NPV = TN/(TN+FN),
      precision = TP/(TP + FP),
      accuracy = (TP + TN)/(p + n),
      MCC = sqrt(precision)*sqrt(TPR)*sqrt(1-FPR)*sqrt(NPV) - sqrt(1-precision)*sqrt(1-TPR)*sqrt(FPR)*sqrt(1-NPV)
    )
}

create_roc <- function(p_val_col, data, regex, cl = NULL){
  data <- data %>%
    split.data.frame(.$comparison)
  out <- map(data,
             pull, p_val_col
  ) %>%
    map(unique) %>%
    map(c, 1.1) %>%
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
  pmap(list(out, data, p, n), create_roc_helper, cl, p_val_col) %>%
    bind_rows()
}

create_roc_helper <- function(out, data, p, n, cl, p_val_col) {
  if(!is.null(cl)){
    multidplyr::cluster_copy(cl,
                             c(
                               "data",
                               "regex",
                               "p",
                               "n"
                             )
    )
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
  sum(dx *my)
}

#### Setup variables ####
yeast_design <- model.matrix(~0+factor(rep(1:2, each = 3)))
colnames(yeast_design) <- paste0("ng", c(50, 100))
yeast_contrast <- limma::makeContrasts('ng50-ng100', levels = yeast_design)
yeast_cont <- matrix(c(1, -1), 2)

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

workers <- round(parallel::detectCores()/2)

#### Normalize and partition ####
yeast_norm <- mavis::yeast %>%
  psrn("identifier")

yeast_part <- yeast_norm %>%
  single_imputation(yeast_design, workers = workers) %>%
  calculate_mean_sd_trends(yeast_design) 

yeast_part <- yeast_part %>%
  grid_search(yeast_design, n_h1 = 30, n_h2 = 30, formula = c(sd ~ mean + c), workers = workers, h1_prop = c(12e-3, .01), h2_prop = c(.0025, .07))

yeast_part <- yeast_part$clustered_data[[1]]

yeast_part_na <- yeast_norm %>% 
  grid_search(yeast_design, n_h1 = 25, n_h2 = 25, formula = c(sd ~ mean + c), workers = workers)
  
yeast_part_na <- yeast_part_na$clustered_data[[1]]


#### Fit regression models ####
gr  <- fit_gamma_regression(yeast_part)
gcr <- fit_gamma_regression(yeast_part, sd ~ mean+c)

#### Run Baldur and t-test ####
yeast_gr_baldur  <- baldur_wrapper(yeast_part, yeast_design, yeast_cont, gr,  workers)
yeast_gcr_baldur <- baldur_wrapper(yeast_part, yeast_design, yeast_cont, gcr, workers)
yeast_ttest      <- ttest_wrapper(colnames(yeast_contrast), yeast_part)

#### Run limma ####
cl <- multidplyr::new_cluster(workers)
multidplyr::cluster_library(
  cl,
  c(
    'purrr',
    'dplyr',
    'magrittr',
    'mavis',
    'stringr'
  )
)
multidplyr::cluster_copy(cl,
                         "calc_tpr_fpr"
)
yeast_trend <- multiple_imputation_and_limma(
  yeast_norm, yeast_design, yeast_contrast, 1000, workers, "identifier", TRUE,
  FALSE, FALSE
) 

yeast_trend_roc <- yeast_trend %>%
  find_alpha(yeast_norm, "YEAST", cl) %>%
  mutate(
    method = "Limma-trend"
  )
yeast_trend <- yeast_trend %>% 
  extract_results(alpha = .05, abs_lfc = 0, pcor = "fdr", id_col = "identifier")

yeast_gr_limma <- multiple_imputation_and_limma(
  yeast_norm, yeast_design, yeast_contrast, 1000, workers, "identifier", TRUE,
  TRUE, FALSE, sd_p ~ mean, sd ~ mean
) 

yeast_gr_limma_roc <- yeast_gr_limma %>%
  find_alpha(yeast_norm, "YEAST", cl) %>%
  mutate(
    method = "GR-Limma"
  )

yeast_gr_limma <- yeast_gr_limma %>% 
  extract_results(alpha = .05, abs_lfc = 0, pcor = "fdr", id_col = "identifier")

yeast_gcr_limma <- multiple_imputation_and_limma(
  yeast_part_na, yeast_design, yeast_contrast, 1000, workers, "identifier", TRUE,
  TRUE, FALSE, sd_p ~ mean, sd ~ mean+c
) 

yeast_gcr_limma_roc <- yeast_gcr_limma %>%
  find_alpha(yeast_norm, "YEAST", cl) %>%
  mutate(
    method = "GCR-Limma"
  )
yeast_gcr_limma <- yeast_gcr_limma %>% 
  extract_results(alpha = .05, abs_lfc = 0, pcor = "fdr", id_col = "identifier")

rm(cl)
gc()

#### Run Mavis ####
data_regex <- "yeast_.*(limma|baldur|trend)$"


yeast_all_err <- mget(ls(pattern = data_regex)) %>%
  map(
    select, 1, comparison, contains('err'), matches('(median_|^)p_val')
  ) %>%
  imap(~ rename(.x, !!.y := 3)) %>%
  reduce(left_join, by = c("identifier", "comparison")) %>%
  split.data.frame(.$comparison) %>%
  map(select, -comparison)

yeast_mavis <- yeast_all_err$`ng50 vs ng100` %>% 
  mavis()

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

yeast_gr_baldur_roc  <- create_roc("err", yeast_gr_baldur,  "YEAST", cl = cl) %>%
  mutate(
    method = 'GR-Baldur'
  )
yeast_gcr_baldur_roc <- create_roc("err", yeast_gcr_baldur, "YEAST", cl = cl) %>%
  mutate(
    method = 'GCR-Baldur'
  )
yeast_ttest_roc      <- create_roc("p_val", yeast_ttest,    "YEAST", cl = cl) %>%
  mutate(
    method = 't-test'
  )
yeast_mavis_roc <- create_roc('mavis', yeast_mavis, "YEAST", cl) %>%
  mutate(
    method = 'Mavis'
  )
rm(cl)
gc()

yeast_roc <- mget(
  ls(pattern = "^yeast.*_roc")
) %>%
  bind_rows()

yeast_auroc <- yeast_roc %>%
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
save(yeast_roc, yeast_auroc, file = "roc_data/yeast_roc.RData")

yeast_roc %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  geom_text(data = yeast_auroc,
            aes(FPR, TPR, label = round(auROC, 3)),
            size = 3,
            show.legend = F
  ) +
  scale_color_manual('Method',
                     values = color_theme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.7, .3),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.key.size = unit(.2, "cm"),
    legend.spacing.y = unit(.01, 'cm')
  ) +
  labs(
    y = 'True Positive Rate',
    x = 'False Positive Rate',
    shape = expression(alpha~phantom())
  )
ggsave_wrapper('yeast_roc', half_page, half_page)

yeast_roc %>%
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
                       ncol = 1,
                       title.position = "top",
                       title.vjust = -2.5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0, 1, .2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0, 1, .05), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.795, .24),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.spacing = unit(0, 'mm')
  ) +
  labs(
    x = expression(alpha),
    y = 'Matthews Correlation Coefficient'
  )
ggsave_wrapper('yeast_mcc', width = half_page, height = half_page)

yeast_roc %>%
  filter(between(alpha, 0, .2)) %>%
  pivot_longer(c(TPR, FPR, precision)) %>%
  mutate(
    name = str_replace(name, 'pre', 'Pre'),
    name = str_replace(name, 'FPR', 'False Positive Rate'),
    name = str_replace(name, 'TPR', 'True Positive Rate')
  ) %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(alpha, value, color = method)) +
  geom_vline(xintercept = .05, linetype = 'dotted', color = 'red') +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  scale_color_manual('Method',
                     values = color_theme,
                     guide = guide_legend(
                       ncol = 1,
                       title.position = "top",
                       title.vjust = -2.5,
                       override.aes = list(linewidth = 1)
                     )
  ) +
  scale_y_continuous(breaks = seq(0, 1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0, 1, .05), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.9, .24),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.key.size = unit(.25, "cm"),
    legend.background = element_blank()
  ) +
  labs(
    x = expression(alpha),
    y = 'True-, False- Positive Rate, Precision'
  ) +
  facet_grid(.~name)
ggsave_wrapper('yeast_decomposed', width = full_page, height = 2/3*full_page)
