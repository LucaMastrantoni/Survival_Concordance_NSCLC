library(readxl)

m <- read_xlsx("ALEX_KM.xlsx")

fit_cox_bayes <- brm(
  formula = as.formula(paste0("time | cens(1 - status) ~ ", "treat")),
  data   = m,
  family = brmsfamily("cox"),
  prior  = prior(normal(0, 3), class = "b"),
  chains = 4, iter = 8000, warmup = 2000,
  cores  = 4,
  seed   = 123,
  backend = "cmdstanr"
)

post <- posterior::as_draws_df(fit_cox_bayes) |>
  dplyr::select(starts_with("b_treat"))

saveRDS(post, file = "posterior_b_treat.rds")