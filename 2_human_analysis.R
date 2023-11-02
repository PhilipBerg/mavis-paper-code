### Author: Mohit Verma
### Packages ####
library(pacman)
p_load(dplyr, multidplyr, magrittr ,purrr, tidyr, ggplot2, mavis, rlang, stringr)

### Functions ####
data <- read.csv(file = "human_data_without_control.csv", sep= "\t", stringsAsFactors = F)%>%
  filter(
    identifier %in% readr::read_csv("human_ids.csv")$identifier # Put correct path to file here
  )
# Wrapper for baldur
baldur_wrapper <- function(data, design, contrast, gam_model, workers) {
  uncertainty <- data %>%
    mavis::estimate_uncertainty('identifier', design, gam_model)
  gam_model %>%
    estimate_gamma_hyperparameters(data, design) %>%
    infer_data_and_decision_model(
      id_col = 'identifier',
      design_matrix = design,
      contrast_matrix = contrast,
      uncertainty_matrix = uncertainty,
      clusters = workers
    )
}

# Function for finding alpha used in creating ROC curves for multiple imputation
find_alpha <- function(results, data, cluster) {
  alpha <- results$limma_results[[1]] %>%
    filter(comparison == .data$comparison[[1]]) %>%
    use_series(p_val) %>%
    unlist(T, F) %>%
    unique() %>%
    c(., 1)
  alpha %>%
    tibble::enframe(name = NULL) %>%
    multidplyr::partition(cluster) %>%
    mutate(
      tmp = map(value,
                ~ extract_results(
                  results,
                  data,
                  alpha = .x,
                  abs_lfc = 0,
                  id_col = 'identifier',
                  pcor = 'fdr'
                )
      )
    ) %>%
    collect() %>%
    unnest(tmp) %>%
    filter(binom_p_value < .05) %>%
    group_by(identifier, comparison) %>%
    summarise(
      min_alpha = min(value)
    ) %>%
    ungroup()
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
    select(identifier, comparison, p_val)
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
      MCC = sqrt(precision)*sqrt(TPR)*sqrt(1-FPR)*sqrt(NPV) - sqrt(1-precision)*sqrt(1-TPR)*sqrt(FPR)*sqrt(1-NPV),
      MCC2 = (TP*TN - FP*FN) / sqrt(
        (TP+FP) * (TP + FN) * (TN + FP) * (TN + FN)
      )
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

### Setup ####
## Figures
sinlge <- 3.33
double <- 7
### Setup ####
## Design matrix
human_design <- model.matrix(~ 0 + factor(rep(1:3, each = 23)))
colnames(human_design) <- paste0('A', c( 100, 150, 200))

## Limma contrast
human_contrast <- limma::makeContrasts(
  contrasts = paste0('A', c( 100, 100, 150), '-A', c( 150, 200, 200)),
  levels = human_design
)
## Baldur contrast
human_cont <- matrix(
  c(
    -1, 1, 0,
    -1, 0, 1,
    0, -1, 1
  ), nrow = 3
)

## Set number of workers for parallel computation
workers <- round(parallel::detectCores()/2)

conditions <- colSums(human_design) %>%
  imap(
    ~ paste0(.y, '_1:', .y, '_', .x)
  )

filter_all_miss <- parse_expr(
  paste0(
    '!if_all(',
    word(conditions[[1]], 1, sep = ':'),
    ':',
    word(conditions[[length(conditions)]], 2, sep = ':'),
    ', is.na)'
  )
)

### Pre-processing and Imputation ####
human <- data %>%
   ## Remove data where all values are missing
   filter(!!filter_all_miss) %>%
   ## Inverse-Log
    mutate(
        across(matches(str_flatten(colnames(human_design), "|")), ~ raise_to_power(2, .))
    ) %>%
    ## Normalize
    psrn("identifier") %>%
  ## M-V trends pre-imputation
  calculate_mean_sd_trends(human_design)


## Run single imputation with one trend
human_sin <- human %>%
  single_imputation(human_design, workers = workers) %>%
  ## Update M-V trends
  calculate_mean_sd_trends(human_design)%>%
  mavis::trend_partitioning(human_design, eps = c(.2, .9), eps_scale = 'linear')
  print("human_sin")
save(human_sin, file='human_sin.RData')

## Gamma regressions for single trend imputation
# Gamma regressions for one trend
gam_sin_sin <- human_sin %>%
  fit_gamma_regression(sd ~ mean)
print("gam_sin_sin")
save(gam_sin_sin, file='human_p_gam_sin_sin.RData')

# Gamma regressions for one trend
gam_sin_mix <- human_sin %>%
  fit_gamma_regression(sd ~ mean*c)
print("gam_sin_mix")
save(gam_sin_mix, file='human_p_gam_sin_mix.RData')

## Run baldur for single trend imputation
# One trend imputation and one trend prior
bald_sin_sin <- baldur_wrapper(human_sin, human_design, human_cont, gam_sin_sin, workers)
print("bald_sin_sin")
save(bald_sin_sin, file='human_p_bald_sin_sin.RData')

# One trend imputation and mixed trends prior
bald_sin_mix <- baldur_wrapper(human_sin, human_design, human_cont, gam_sin_mix, workers)
print("bald_sin_mix")
save(bald_sin_mix, file='human_p_bald_sin_mix.RData')

### Multiple imputation and limma
## Run mixture

limma_sin_mix <- multiple_imputation_and_limma(
  data = human,
  design = human_design,
  contrast_matrix = human_contrast,
  imputations = 1000,
  workers = workers,
  id_col = 'identifier',
  .robust = T,
  weights = T,
  plot_trend = F,
  formula_weights =  sd ~ mean * c
)
print("limma_sin_mix")
save(limma_sin_mix, file='human_p_limma_sin_mix.RData')

limma_sin_sin <- multiple_imputation_and_limma(
  data = human,
  design = human_design,
  contrast_matrix = human_contrast,
  imputations = 1000,
  workers = workers,
  id_col = 'identifier',
  .robust = T,
  weights = T,
  plot_trend = F,
  formula_weights =  sd ~ mean
)
print("limma_sin_sin")
save(limma_sin_sin, file='human_p_limma_sin_sin.RData')

## Run limma-trend
# Single imp
limma_trend_sin <- multiple_imputation_and_limma(
  data = human,
  design = human_design,
  contrast_matrix = human_contrast,
  imputations = 1000,
  workers = workers,
  id_col = 'identifier',
  .robust = T,
  weights = F,
  plot_trend = F,
)
print("limma_trend_sin")
save(limma_trend_sin, file='human_p_limma_trend_sin.RData')

## Find what alpha value peptides became significant at
## Lets do this with parallel computation to reduce the time
cl <- multidplyr::new_cluster(workers)
multidplyr::cluster_library(
  cl,
  c(
    'purrr',
    'dplyr',
    'magrittr',
    'mavis'
  )
)

#### Mult Imp Limma ####
limma_sin_mix_alpha <- find_alpha(limma_sin_mix, human_sin, cl)
print("limma_sin_mix_alpha")
save(limma_sin_mix_alpha, file='human_p_limma_sin_mix_alpha.RData')

## Single trend
limma_sin_sin_alpha <- find_alpha(limma_sin_sin, human_sin, cl)
print("limma_sin_sin_alpha")
save(limma_sin_sin_alpha, file='human_p_limma_sin_sin_alpha.RData')

#### Limma trend ####
# Single
limma_trend_sin_alpha <- find_alpha(limma_trend_sin, human_sin, cl)
print("limma_trend_sin_alpha")
save(limma_trend_sin_alpha, file='human_p_limma_trend_sin_alpha.RData')

# Lets clean up our parallel computation
rm(cl)
gc()


### t-test
## single
sin_t_test <- map(colnames(human_contrast), ttest_wrapper, human_sin) %>%
  bind_rows() %>%
  group_by(comparison) %>%
    mutate(
      p_val = p.adjust(p_val, 'fdr')
    ) %>%
  ungroup()

### Performance ####
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

#### Baldur Number:-3 ####
sin_sin_roc <- create_roc('err', bald_sin_sin, 'ECOLI', cl) %>%
  mutate(
    method = 'Single-Imputation+Gamma-Baldur'
  )
print("sin_sin_roc")
save(sin_sin_roc, file='human_p_sin_sin_roc.RData')

#### Baldur Number:-4 ####
sin_mix_roc <- create_roc('err', bald_sin_mix, 'ECOLI', cl) %>%
  mutate(
    method = 'Single-Imputation+Mix-Baldur'
  )
print("sin_mix_roc")
save(sin_mix_roc, file='human_p_sin_mix_roc.RData')

#### limma Number:-5 ####
sin_sin_gam_roc <- create_roc('min_alpha', limma_sin_sin_alpha, 'ECOLI', cl) %>%
  mutate(
    method = 'Single-Imputation+Gamma-Limma'
  )
print("sin_sin_gam_roc")
save(sin_sin_gam_roc, file='human_p_sin_sin_gam_roc.RData')

### limma Number:-7 ####
sin_mix_gam_roc <- create_roc('min_alpha', limma_sin_mix_alpha, 'ECOLI', cl) %>%
  mutate(
    method = 'Single-Imputation+Mix-Limma'
  )
print("sin_mix_gam_roc")
save(sin_mix_gam_roc, file='human_p_sin_mix_gam_roc.RData')

#### Limma trend Number:-9 ####
sin_trend_roc <- create_roc('min_alpha', limma_trend_sin_alpha, 'ECOLI', cl) %>%
  mutate(
    method = 'Single-Imputation+Limma-trend'
  )
print("sin_trend_roc")
save(sin_trend_roc, file='human_p_sin_trend_roc.RData')

#### T-test Number:-12 ####
sin_ttest_roc <- create_roc('p_val', sin_t_test, 'ECOLI', cl) %>%
  mutate(
    method = 'Single-Imputation+t-test'
  )
print("sin_ttest_roc")
save(sin_ttest_roc, file='human_p_sin_ttest_roc.RData')


rm(cl)
gc()

# Starting from the performance data
human_roc <- bind_rows(
  sin_sin_roc,
  sin_mix_roc,
  sin_mix_gam_roc,
  sin_sin_gam_roc,
  sin_trend_roc,
  sin_ttest_roc
) %>%
  mutate(
    comparison = str_replace(comparison, 'Vs', 'vs'),
    method = str_replace_all(method, set_names(c("GCR", "GR"),c("Mix", "Gamma")))
  ) %>%

  # Remove the ones with all values missing
  # for baldur without imputation
  filter(alpha < 2) %>%
  # Separate imputation and decision into two columns
  tidyr::separate(method, c('imputation', 'method'), '\\+') %>%
  # Split on the imputation methods
  split.data.frame(.$imputation)

human_roc<- mutate(human_roc$`Single-Imputation`, comparison=str_replace(comparison, 'A200 vs A150', 'A150 vs A200'))
human_roc<-mutate(human_roc, comparison=str_replace(comparison, 'A150 vs A100', 'A100 vs A150'))
human_roc<-mutate(human_roc, comparison=str_replace(comparison, 'A200 vs A100', 'A100 vs A200'))

human_roc<-list(human_roc=human_roc)


# Get the auroc
human_auroc <- human_roc %>%
  map(
    ~ .x %>%
      group_by(comparison, method) %>%
      summarise(
        auROC = integrater(FPR, TPR)
      ) %>%
      arrange(desc(auROC)) %>%
      mutate(
        FPR = seq(0.1, .98, length.out = n()),
        TPR = 0
      )
  )

### Plotting ####
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

# ROC curves
map2(human_roc, human_auroc,
          ~ .x %>%
            group_by(comparison, method) %>%
            arrange(alpha) %>%
            ggplot(aes(FPR, TPR, color = method)) +
            geom_abline(intercept = 0, slope = 1, linetype = 'dashed') +
            geom_path(linewidth = 1/4) +
            theme_classic() +
            geom_text(data = .y,
                      aes(FPR, TPR, label = round(auROC, 1.5)),
                      size = 3.5,
                      show.legend = F
            ) +
            scale_color_manual('Methods',
                               values = color_theme,
                               guide = guide_legend(
                                 ncol = 1,
                                 title.position = "top",
                                 title.hjust = .05,
                                 override.aes = list(linewidth=1)
                               )
            ) +
            scale_y_continuous(breaks = seq(0, 1, .1), limits = c(0,1)) +
            scale_x_continuous(labels = function(x) ifelse(x == 0, "0", x), breaks = seq(0, 1, .1)) +
            theme(
              legend.position = c(.830, .095),
              plot.margin = margin(l = 13.5),
              legend.direction = 'horizontal',
              legend.background = element_blank(),
              legend.key.size = unit(.25, "cm"),
              legend.spacing.y = unit(.01, 'cm'),
              axis.title.y = element_blank()
            ) +
            labs(
              x = 'False Positive Rate',
              y = 'True Positive Rate'
            ) +
            facet_wrap(comparison~., ncol = 1)
) %>%
  iwalk(
    ~ ggsave(paste0('human_prog_roc_', .y, '.png'), .x,  width = 2.99, units = "in")
  )


map(human_roc,
        ~ .x %>%
          pivot_longer(c(TPR, FPR, precision)) %>%
          mutate(
            intercept = if_else(name == 'precision', 448/2225, 0),
            slope = if_else(name == 'precision', 0, 1),
            name = str_replace(name, 'pre', 'Pre'),
            name = str_replace(name, 'FPR', 'False Positive Rate'),
            name = str_replace(name, 'TPR', 'True Positive Rate')
          ) %>%
          group_by(comparison, method) %>%
          arrange(alpha) %>%
          ggplot(aes(alpha, value, color = method)) +
          geom_vline(xintercept = .05, linetype = 'dashed', color = 'red') +
          geom_abline(aes(intercept = intercept, slope = slope), linetype = 'dashed') +
          geom_path(linewidth = 1/4) +
          theme_classic() +
          scale_color_manual('Methods',
                             values = color_theme,
                             guide = guide_legend(
                               ncol = 1,
                               title.position = "top",
                               title.hjust = .5,
                               override.aes = list(linewidth=1)
                             )
          ) +
          scale_y_continuous(breaks = seq(0,1,.1), limits = c(0,1)) +
          theme(
            legend.position = c(.915, .08),
            plot.margin = margin(),
            legend.direction = 'horizontal',
            legend.background = element_blank(),
            legend.key.size = unit(.25, "cm")
          ) +
          labs(
            x = expression(alpha),
            y = 'True-, False- Positive Rate, Precision'
          ) +
          facet_wrap(comparison~name, scales = 'free')
) %>%
  # Save the plots
  iwalk(
    ~ ggsave(paste0('human_prog_decomposed_', .y, '.png'), .x, width = 7, units = "in")
  )

map(human_roc,
        ~ .x %>%
          group_by(comparison, method) %>%
          arrange(alpha) %>%
          ggplot(aes(alpha, MCC, color = method)) +
          geom_vline(xintercept = .05, linetype = 'dashed') +
          geom_path(linewidth = 1/4) +
          theme_classic() +
          scale_color_manual('Methods',
                             values = color_theme,
                             guide = guide_legend(
                               ncol = 1,
                               title.position = "top",
                               title.hjust = .5,
                               override.aes = list(linewidth=1)
                             )
          ) +
          scale_y_continuous(breaks = seq(0, 1, .1), limits = c(0,1)) +
          theme(
            legend.position = c(.830, .90),
            plot.margin = margin(),
            legend.direction = 'horizontal',
            legend.background = element_blank(),
            legend.key.size = unit(.25, "cm")
          ) +
          labs(
            x = expression(alpha),
            y = 'Matthews Correlation Coefficient'
          ) +
          facet_wrap(comparison~.,ncol = 1)
) %>%
  # Save the plots
  iwalk(
    ~ ggsave(paste0('human_prog_mcc_', .y, '.png'), .x, width = 2.99, units = "in")
  )

