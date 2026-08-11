# Metodo IDF do SGB aplicado a df_imax_wilson.pqt
#
# Dados:
#   gauge_code = codigo da estacao
#   wateryear  = ano hidrologico
#   d          = duracao (horas)
#   imax       = intensidade maxima anual (mm/h)
#
# Etapas:
#   1) GEV por momentos-L para cada estacao e duracao;
#   2) quantis para os tempos de retorno do SGB;
#   3) ajuste da equacao IDF
#          i = k * T^gamma / (t + b)^c
#      com t em minutos.

library(arrow)
library(dplyr)
library(lmom)

T_sgb <- c(2, 5, 10, 15, 20, 25, 30, 40, 50, 60, 75, 100)

fit_gev <- function(x, T = T_sgb) {
  par <- pelgev(samlmu(x, nmom = 3))
  
  tibble(
    T = T,
    intensity = quagev(1 - 1 / T, par),
    gev_location = par[1],
    gev_scale = par[2],
    gev_shape = par[3]
  )
}

fit_idf_sgb <- function(data, objective = c("REQM", "DPMA"),
                        b_interval = c(0, 120)) {
  objective <- match.arg(objective)
  
  b0 <- optimize(
    function(b) {
      data %>%
        group_by(T) %>%
        summarise(
          r2 = summary(lm(log(intensity) ~ log(duration_min + b)))$r.squared,
          .groups = "drop"
        ) %>%
        summarise(r2 = mean(r2)) %>%
        pull(r2) %>%
        `-`()
    },
    interval = b_interval
  )$minimum
  
  initial_by_T <- data %>%
    group_by(T) %>%
    group_modify(~ {
      m <- lm(log(intensity) ~ log(duration_min + b0), data = .x)
      tibble(
        A = exp(unname(coef(m)[1])),
        c = -unname(coef(m)[2])
      )
    }) %>%
    ungroup()
  
  mA <- lm(log(A) ~ log(T), data = initial_by_T)
  
  theta0 <- c(
    log(exp(unname(coef(mA)[1]))),
    unname(coef(mA)[2]),
    log(max(b0, 1e-6)),
    log(mean(initial_by_T$c))
  )
  
  objective_function <- function(theta) {
    k <- exp(theta[1])
    gamma <- theta[2]
    b <- exp(theta[3])
    c <- exp(theta[4])
    
    calc <- with(data, k * T^gamma / (duration_min + b)^c)
    
    if (objective == "REQM") {
      sqrt(sum((data$intensity - calc)^2) / (nrow(data) - 1))
    } else {
      mean(abs(data$intensity - calc) / data$intensity) * 100
    }
  }
  
  opt <- optim(
    par = theta0,
    fn = objective_function,
    method = "Nelder-Mead",
    control = list(maxit = 50000, reltol = 1e-12)
  )
  
  pars <- c(
    k = exp(opt$par[1]),
    gamma = opt$par[2],
    b = exp(opt$par[3]),
    c = exp(opt$par[4])
  )
  
  fitted <- data %>%
    mutate(
      fitted = pars["k"] * T^pars["gamma"] /
        (duration_min + pars["b"])^pars["c"]
    )
  
  tibble(
    k = pars["k"],
    gamma = pars["gamma"],
    b = pars["b"],
    c = pars["c"],
    REQM = sqrt(sum((fitted$intensity - fitted$fitted)^2) /
                  (nrow(fitted) - 1)),
    DPMA = mean(abs(fitted$intensity - fitted$fitted) /
                  fitted$intensity) * 100,
    convergence = opt$convergence
  )
}

# -----------------------------------------------------------------------------
# Aplicacao ao banco de maximas anuais
# -----------------------------------------------------------------------------

df_imax <- read_parquet(
  "C:/Users/wilso/OneDrive/My papers/IDF_Dirceu/database/df_imax_wilson.pqt"
)

station_info <- df_imax %>%
  group_by(gauge_code) %>%
  summarise(
    n_years = first(n_years),
    first_year = min(wateryear),
    last_year = max(wateryear),
    .groups = "drop"
  )

# GEV por estacao e duracao.
# No arquivo, d varia de 1 a 24 horas; a equacao do SGB usa t em minutos.
idf_quantiles <- df_imax %>%
  mutate(duration_min = 60 * d) %>%
  group_by(gauge_code, d, duration_min) %>%
  group_modify(~ fit_gev(.x$imax)) %>%
  ungroup()

# Parametros GEV: uma linha por estacao e duracao.
gev_parameters <- idf_quantiles %>%
  distinct(
    gauge_code, d, duration_min,
    gev_location, gev_scale, gev_shape
  )

# Ajuste da IDF para cada estacao usando a REQM, como no teste de Itauna.
idf_results <- idf_quantiles %>%
  group_by(gauge_code) %>%
  group_modify(~ fit_idf_sgb(.x, objective = "REQM")) %>%
  ungroup() %>%
  left_join(station_info, by = "gauge_code") %>%
  mutate(DPMA_le_10 = DPMA <= 10) %>%
  select(
    gauge_code, n_years, first_year, last_year,
    k, gamma, b, c, REQM, DPMA, DPMA_le_10, convergence
  )

idf_results

write.csv(idf_results, "idf_sgb_results_reqm.csv", row.names = FALSE)
write.csv(gev_parameters, "gev_parameters_by_duration.csv", row.names = FALSE)