### Author: Philip Berg
p_load(dplyr, tidyr, ggplot2, stringr, purrr, mavis, readr, readxl, tibble)

ggsave_wrapper <- function(file_name, width, height = NA, dpi = 1200){
  ggsave(paste0(file_name, '.tiff'), width = width, height = height, units = 'mm', dpi = dpi, compression = 'lzw')
  ggsave(paste0(file_name, '.pdf'), width = width, height = height, units = 'mm', dpi = dpi)
}

ups_design <- model.matrix(~0+factor(rep(1:3, each = 4)))
colnames(ups_design) <- paste0("fmol", c(25, 50, 100))
human_design <- model.matrix(~0+factor(rep(1:3, each = 23)))
colnames(human_design) <- paste0("spike_prop_", c(25, 12, 6))
yeast_design <- model.matrix(~0+factor(rep(1:2, each = 3)))
colnames(yeast_design) <- paste0("ng", c(50, 100))
ramus_design <- model.matrix(~0+factor(rep(1:9, each = 3)))
colnames(ramus_design) <- paste0("condi", 1:9)

human_part <- readxl::read_excel('diaWorkflowResults_allDilutions.xlsx', na = 'NA') %>%
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
  psrn('identifier') %>%
  single_imputation(human_design, workers = 10) %>%
  calculate_mean_sd_trends(human_design) %>%
  trend_partitioning(human_design, sd ~ mean * c, eps = c(.2, .9), eps_scale = 'linear')

yeast_part <- mavis::yeast %>%
  psrn("identifier") %>%
  single_imputation(yeast_design, workers = 10) %>%
  calculate_mean_sd_trends(yeast_design) %>%
  trend_partitioning(yeast_design, sd ~ mean * c, eps = c(.1, .4), eps_scale = 'log-linear')

readr::read_csv('https://figshare.com/ndownloader/files/35592290?private_link=28e837bfe865e8f13479', show_col_types = FALSE) %>%
  janitor::clean_names() %>%
  mutate(
    across(where(is.numeric), ~na_if(.x, 0))
  ) %>%
  readr::write_csv('ramus_clean.csv')

ramus_part <- readr::read_csv('ramus_clean.csv') %>%
  rename_with(~paste0('condi', rep(1:9, each = 3), '_', rep(1:3, lenght.out = 9*3)), where(is.numeric)) %>%
  psrn(load_info = F, id_col = 'identifier') %>%
  single_imputation(ramus_design, workers = 10) %>%
  calculate_mean_sd_trends(ramus_design) %>%
  trend_partitioning(ramus_design, sd ~ mean * c, eps = c(.1, .4), eps_scale = 'log-linear')

ups_part <- mavis::ups %>%
  psrn("identifier") %>%
  single_imputation(ups_design, workers = 10) %>%
  calculate_mean_sd_trends(ups_design) %>%
  trend_partitioning(ups_design, sd ~ mean * c, eps = c(.1, .6), eps_scale = 'log-linear')


full_page <- 178
half_page <- 86

gr_plots <- mget(paste0(c('yeast', 'ups', 'ramus', 'human'), '_part')) %>%
  map(
    plot_gamma
  ) %>%
  map2(c('Yeast', 'UPS', 'Ramus', 'Human'),
    ~ .x +
    ggtitle(.y) +
      scale_shape_manual(values = 16) +
      scale_size_manual(values = 0)
  ) %>%
  map(
    ~ {.x$layers[[1]]$aes_params$shape <- '.'; .x$layers[[2]]$aes_params$size <- .3; .x}
  ) %>%
  plot_grid(plotlist = ., align = 'hv', nrow = 1)

eps <- list(
  yeast = c(.1, .4),
  ups = c(.1, .6),
  ramus = c(.1, .4),
  human = c(.2, .9)
)
eps_scale <- list(
  yeast = 'log-linear',
  ups = 'log-linear',
  ramus = 'log-linear',
  human = 'linear'
)

gmr_plots <- mget(paste0(c('yeast', 'ups', 'ramus', 'human'), '_part')) %>%
  map(
    plot_gamma_partition, formula = sd ~ mean*c
  ) %>%
  map(
       ~ .x +
         theme(
           plot.title = element_blank(),
           legend.position = c(.79, .79),
           legend.spacing.y = unit(0, 'mm')
           ) +
         guides(
           color = guide_legend(
             override.aes = list(size=1, shape = 16),
             keywidth = unit(.1, 'mm'),
             keyheight = unit(.1, 'mm')
             )
           )
  ) %>%
  map(
    ~ {
      .x$layers[[1]]$aes_params$shape <- '.';
      .x$layers[[2]]$aes_params$size <- .3;
      .x$layers[[3]]$aes_params$size <- .3;
      .x
    }
  ) %>%
  plot_grid(plotlist = ., align = 'hv', nrow = 1)


plot_grid(gr_plots, gmr_plots, align = 'hv', nrow = 2, labels = 'AUTO')
ggsave_wrapper('trends_plot', full_page, full_page/2)


mget(paste0(c('yeast', 'ups', 'ramus', 'human'), '_part')) %>%
  map(
    fit_gamma_regression, sd ~ mean
  ) %>%
  map_dbl(
    ~ 1 - .x$deviance/.x$null.deviance
  ) %>%
  tibble::enframe()

mget(paste0(c('yeast', 'ups', 'ramus', 'human'), '_part')) %>%
  map(
    fit_gamma_regression, sd ~ mean*c
  ) %>%
  map_dbl(
    ~ 1 - .x$deviance/.x$null.deviance
  ) %>%
  tibble::enframe()




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

all_data <- list.files("roc_data/", full.names = T) %>%
  map(~ {load(.x); mget(ls(pattern = '_roc'))}) %>%
  map(magrittr::extract2, 1) %>%
  imap(
    ~ mutate(.x, data = case_when(
      .y %in% 1:6 ~ "Human",
      .y == 7 ~ "Ramus",
      .y == 8 ~ "UPS",
      .y == 9 ~ "Yeast",
    ))
  ) %>%
  bind_rows() %>%
  select(-MCC2) %>%
  mutate(
    method = str_remove(method, "Single-Imputation\\+"),
    method = str_replace_all(method, c("Mix" = "GCR", "Gamma" = "GR")),
    method = str_replace(method, "limma", "Limma")
  )

all_data2 <- all_data %>%
  select(-c(accuracy, NPV)) %>%
  pivot_longer(TPR:MCC) %>%
  group_by(comparison, method, data, name) %>%
  filter(alpha <= .05) %>%
  slice_max(alpha) %>%
  mutate(
    method = factor(method,
                    levels = c("t-test", "Limma-trend", paste0(c("GR", "GCR"), "-", rep(c("Limma", "Baldur"), each = 2)))
    ),
    name = str_replace(name, "prec", "Prec")
  )


all_data2 %>%
  group_by(name, data, method) %>%
  summarise(
    mean = mean(value, na.rm = T)
  ) %>%
  ggplot(aes(name, mean, fill = method)) +
  geom_col(position = position_dodge(1)) +
  geom_point(data = all_data2, aes(name, value), position = position_dodge(1), color = 'grey', size = 1/10) +
  scale_fill_manual("Method", values = color_theme) +
  coord_flip() +
  theme_classic() +
  theme(
    legend.position = "top",
    plot.margin = margin()
  ) +
  facet_grid(. ~ data) +
  scale_y_continuous(breaks = seq(0,1,.25), labels = function(x) case_when(x == 0 ~ "0", x == 1 ~ "1", T ~ as.character(x))) +
  labs(x = "Metric", y = "Value")
ggsave_wrapper("summary_bar", full_page, 2/3*full_page, dpi = 1200)
