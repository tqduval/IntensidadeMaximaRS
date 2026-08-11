# Metodo IDF do SGB 
#
# Etapas:
#   1) ajuste da GEV por momentos-L (pacote lmom);
#   2) calculo dos quantis para os tempos de retorno;
#   3) ajuste da equacao IDF:
#
#        i = k * T^gamma / (t + b)^c
#
#      i = intensidade (mm/h)
#      T = tempo de retorno (anos)
#      t = duracao (minutos)
#
# O SGB usa o GRG2 do Excel Solver. Aqui a mesma funcao objetivo e ajustada
# com optim(..., method = "Nelder-Mead"), mantendo o codigo inteiramente em R.

library(dplyr)
library(tidyr)
library(lmom)

# -----------------------------------------------------------------------------
# 1. GEV por momentos-L
# -----------------------------------------------------------------------------

fit_gev <- function(x, T = c(2, 5, 10, 15, 20, 25, 30, 40, 50, 60, 75, 100)) {
  par <- pelgev(samlmu(x, nmom = 3))

  list(
    parameters = setNames(par, c("location", "scale", "shape")),
    quantiles = data.frame(
      T = T,
      quantile = quagev(1 - 1 / T, par)
    )
  )
}

# -----------------------------------------------------------------------------
# 2. Ajuste da equacao IDF do SGB
# -----------------------------------------------------------------------------

fit_idf_sgb <- function(data, objective = c("REQM", "DPMA"), b_interval = c(0, 120)) {
  objective <- match.arg(objective)

  data <- data %>%
    filter(
      is.finite(duration_min), duration_min > 0,
      is.finite(T), T > 1,
      is.finite(intensity), intensity > 0
    )

  # Valor inicial de b: procurar a melhor linearizacao de
  # log(i) versus log(t + b), usando o R2 medio entre os tempos de retorno.
  b0 <- optimize(
    function(b) {
      r2 <- data %>%
        group_by(T) %>%
        group_modify(~ {
          m <- lm(log(intensity) ~ log(duration_min + b), data = .x)
          data.frame(r2 = summary(m)$r.squared)
        }) %>%
        ungroup() %>%
        summarise(r2 = mean(r2)) %>%
        pull(r2)

      -r2
    },
    interval = b_interval
  )$minimum

  # Para cada T:
  # log(i) = log(A) - c * log(t + b)
  # A = k * T^gamma
  initial_by_T <- data %>%
    group_by(T) %>%
    group_modify(~ {
      m <- lm(log(intensity) ~ log(duration_min + b0), data = .x)
      data.frame(
        A = exp(unname(coef(m)[1])),
        c = -unname(coef(m)[2])
      )
    }) %>%
    ungroup()

  mA <- lm(log(A) ~ log(T), data = initial_by_T)

  k0 <- exp(unname(coef(mA)[1]))
  gamma0 <- unname(coef(mA)[2])
  c0 <- mean(initial_by_T$c)

  # Parametrizacao logaritmica para manter k, b e c positivos.
  theta0 <- c(log(k0), gamma0, log(max(b0, 1e-6)), log(c0))

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
        (duration_min + pars["b"])^pars["c"],
      error_pct = 100 * (fitted - intensity) / intensity
    )

  metrics <- fitted %>%
    summarise(
      REQM = sqrt(sum((intensity - fitted)^2) / (n() - 1)),
      DPMA = mean(abs(intensity - fitted) / intensity) * 100
    )

  list(
    coefficients = pars,
    metrics = metrics,
    fitted = fitted,
    convergence = opt$convergence
  )
}

# -----------------------------------------------------------------------------
# 3. Teste com toda a tabela publicada para Itauna/MG - SGB, 2020
#    18 duracoes x 12 tempos de retorno = 216 pontos
# -----------------------------------------------------------------------------

itauna <- data.frame(
  duration_min = c(5, 10, 15, 20, 30, 45, 60, 120, 180, 240, 300, 360,
                   420, 480, 720, 840, 1200, 1440),
  `2`   = c(134.1,115.4,101.6,91.0,75.7,60.9,51.2,32.2,23.9,19.2,16.1,14.0,12.4,11.1,8.0,7.1,5.3,4.5),
  `5`   = c(150.8,129.8,114.3,102.3,85.1,68.4,57.6,36.2,26.9,21.6,18.2,15.7,13.9,12.5,9.0,8.0,5.9,5.1),
  `10`  = c(164.7,141.8,124.8,111.8,92.9,74.8,62.9,39.5,29.3,23.6,19.8,17.2,15.2,13.7,9.9,8.7,6.5,5.6),
  `15`  = c(173.5,149.3,131.5,117.7,97.9,78.7,66.3,41.6,30.9,24.8,20.9,18.1,16.0,14.4,10.4,9.2,6.8,5.9),
  `20`  = c(180.0,154.9,136.4,122.1,101.5,81.7,68.8,43.2,32.1,25.8,21.7,18.8,16.6,14.9,10.8,9.5,7.1,6.1),
  `25`  = c(185.2,159.4,140.3,125.7,104.5,84.1,70.7,44.4,33.0,26.5,22.3,19.3,17.1,15.4,11.1,9.8,7.3,6.3),
  `30`  = c(189.5,163.1,143.6,128.6,106.9,86.0,72.4,45.5,33.8,27.1,22.8,19.8,17.5,15.7,11.3,10.0,7.5,6.4),
  `40`  = c(196.6,169.2,149.0,133.4,110.9,89.2,75.1,47.2,35.0,28.2,23.7,20.5,18.2,16.3,11.8,10.4,7.7,6.7),
  `50`  = c(202.3,174.1,153.3,137.3,114.1,91.8,77.3,48.5,36.0,29.0,24.4,21.1,18.7,16.8,12.1,10.7,8.0,6.9),
  `60`  = c(207.1,178.2,156.9,140.5,116.8,94.0,79.1,49.7,36.9,29.6,24.9,21.6,19.1,17.2,12.4,10.9,8.2,7.0),
  `75`  = c(213.1,183.4,161.5,144.6,120.2,96.7,81.4,51.1,38.0,30.5,25.7,22.2,19.7,17.7,12.7,11.2,8.4,7.2),
  `100` = c(221.0,190.2,167.5,150.0,124.7,100.3,84.4,53.0,39.4,31.6,26.6,23.1,20.4,18.3,13.2,11.7,8.7,7.5),
  check.names = FALSE
) %>%
  pivot_longer(
    cols = -duration_min,
    names_to = "T",
    values_to = "intensity"
  ) %>%
  mutate(T = as.numeric(T))

# Ajuste aos 216 valores publicados.
fit_itauna <- fit_idf_sgb(itauna, objective = "REQM")

fit_itauna$coefficients
fit_itauna$metrics

# Comparacao com a equacao publicada pelo SGB:
# i = 1829 * T^0.1277 / (t + 20.4)^0.8351


