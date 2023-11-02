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

# Function for running t-test
ttest_wrapper <- function(contrast, data) {
  comp <- str_split(contrast, '-')[[1]]
  data %>%
    group_by(identifier) %>%
    mutate(
      comparison = str_flatten(comp, ' vs '),
      p_val = t.test(across(contains(comp[1])), across(contains(comp[2])), var.equal = T)$p.value
    ) %>%
    group_by(comparison) %>%
    mutate(
      p_val = p.adjust(p_val, 'fdr')
    ) %>%
    select(identifier, comparison, p_val)
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
  ggsave(paste0(file_name, '.tiff'), width = width, height = height, units = 'in', dpi = 320, compression = 'lzw')
  ggsave(paste0(file_name, '.pdf'), width = width, height = height, units = 'in', dpi = 320)
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
ups_design <- model.matrix(~0+factor(rep(1:3, each = 4)))
colnames(ups_design) <- paste0("fmol", c(25, 50, 100))
ups_contrast <- combn(colnames(ups_design), 2) %>%
  t() %>%
  apply(1, str_flatten, '-') %>%
  limma::makeContrasts(
    contrasts = .,
    levels = ups_design
  )
ups_cont <- matrix(
  c(
    1,   -1,    0,
    1,    0,   -1,
    0,    1,   -1
  ),, 3, byrow = T
) %>%
  t()

full_page <- 178
half_page <- 86
color_theme <- set_names(
  viridisLite::turbo(6, end = .9),
  c(
    'GCR-Baldur',
    'GR-Baldur',
    'GR-Limma',
    'GCR-Limma',
    'Limma-trend',
    't-test'
  )
)

workers <- round(parallel::detectCores()/2)

#### Normalize and partition ####
ups_norm <- mavis::ups %>%
  psrn("identifier")
ups_part <- ups_norm %>%
  single_imputation(ups_design, workers = 10) %>%
  calculate_mean_sd_trends(ups_design) %>%
  trend_partitioning(human_design, sd ~ mean * c, eps = c(.2, .9), eps_scale = 'linear')

#### Fit regression models ####
gr  <- fit_gamma_regression(ups_part)
gcr <- fit_gamma_regression(ups_part, sd ~ mean*c)

#### Run Baldur and t-test ####
ups_gr_baldur  <- baldur_wrapper(ups_part, ups_design, ups_cont, gr,  workers)
ups_gcr_baldur <- baldur_wrapper(ups_part, ups_design, ups_cont, gcr, workers)
ups_ttest      <- map(colnames(ups_contrast), ttest_wrapper, ups_part) %>%
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
ups_trend_roc <- multiple_imputation_and_limma(
  ups_norm, ups_design, ups_contrast, 1000, workers, "identifier", TRUE,
  FALSE, FALSE
) %>%
  find_alpha(ups_norm, "UPS", cl) %>%
  mutate(
    method = "Limma-trend"
  )
ups_gr_limma_roc <- multiple_imputation_and_limma(
  ups_norm, ups_design, ups_contrast, 1000, workers, "identifier", TRUE,
  TRUE, FALSE, sd_p ~ mean, sd ~ mean
) %>%
  find_alpha(ups_norm, "UPS", cl) %>%
  mutate(
    method = "GR-Limma"
  )
ups_gcr_limma_roc <- multiple_imputation_and_limma(
  ups_norm, ups_design, ups_contrast, 1000, workers, "identifier", TRUE,
  TRUE, FALSE, sd_p ~ mean, sd ~ mean*c
) %>%
  find_alpha(ups_norm, "UPS", cl) %>%
  mutate(
    method = "GCR-Limma"
  )
rm(cl)
gc()

#### Performance ####
cl <- multidplyr::new_cluster(workers)
multidplyr::cluster_library(cl,
                            c("dplyr",
                              "stringr",
                              "tidyr",
                              "purrr",
                              "tibble",
                              "magrittr",
                              "stringr"
                            )
)
multidplyr::cluster_copy(cl,
                         "calc_tpr_fpr"
)

ups_gr_baldur_roc  <- create_roc("err", ups_gr_baldur,  "UPS", cl = cl) %>%
  mutate(
    method = 'GR-Baldur'
  )
ups_gcr_baldur_roc <- create_roc("err", ups_gcr_baldur, "UPS", cl = cl) %>%
  mutate(
    method = 'GCR-Baldur'
  )
ups_ttest_roc      <- create_roc("p_val", ups_ttest,    "UPS", cl = cl) %>%
  mutate(
    method = 't-test'
  )
rm(cl)
gc()

# Concatenate results from all methods
ups_roc <- mget(
  ls(pattern = "^ups.*_roc")
) %>%
  bind_rows()

# Calculate auROC
ups_auroc <- ups_roc %>%
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
save(ups_roc, ups_auroc, file = "roc_data/ups_roc.RData")

ups_roc %>%
  group_by(comparison, method) %>%
  arrange(alpha) %>%
  ggplot(aes(FPR, TPR, color = method)) +
  geom_path(linewidth = 1/6) +
  theme_classic() +
  geom_text(data = ups_auroc,
            aes(FPR, TPR, label = round(auROC, 3)),
            size = 2.8,
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
    legend.position = c(.925, .2),
    plot.margin = margin(r = 2),
    legend.direction = 'horizontal',
    legend.key.size = unit(.2, "cm"),
    legend.spacing.y = unit(.01, 'cm')
  ) +
  labs(
    y = 'True Positive Rate',
    x = 'False Positive Rate',
    shape = expression(alpha~phantom())
  ) +
  facet_grid(. ~ comparison)
ggsave_wrapper('ups_roc', full_page, 1/2*full_page)

ups_roc %>%
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
  scale_y_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.925, .895),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.background = element_blank(),
    legend.key.size = unit(.25, "cm"),
    legend.spacing = unit(0, 'mm')
  ) +
  labs(
    x = expression(alpha),
    y = 'Matthews Correlation Coefficient'
  ) +
  facet_grid(. ~ comparison)
ggsave_wrapper('ups_mcc', full_page, 1/2*full_page)

ups_roc %>%
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
  scale_y_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  scale_x_continuous(breaks = seq(0,1,.2), labels = function(x) ifelse(x == 0, "0", x)) +
  theme(
    legend.position = c(.92, .075),
    plot.margin = margin(),
    legend.direction = 'horizontal',
    legend.key.size = unit(.25, "cm"),
    legend.background = element_blank()
  ) +
  labs(
    x = expression(alpha),
    y = 'True-, False- Positive Rate, Precision'
  ) +
  facet_grid(name ~ comparison)
ggsave_wrapper('ups_decomposed', width = full_page, height = full_page)
