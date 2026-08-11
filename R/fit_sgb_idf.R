
# Funções p/ ajuste da IDF-SGB baseado na metodologia
# do Atlas Pluviométrico do Brasil

# Códigos elaborados pelo Wilson

# Dados:
#   gauge_code = codigo da estacao
#   wateryear  = ano hidrologico
#   d          = duracao (horas)
#   imax       = intensidade maxima anual (mm/h)

# Etapas:
#   1) GEV por momentos-L para cada estacao e duracao;
#   2) quantis para os tempos de retorno do SGB;
#   3) ajuste da equacao IDF
#          i = k * T^gamma / (t + b)^c
#      com t em minutos.


# PACOTES ----------------------------------------------------------------

if(!require(pacman)) install.packages("pacman")
pacman::p_load(tidyverse, lmom, arrow)


# FUNÇÕES ----------------------------------------------------------------

# Ajustar GEV e estimar quantis c/ L-momentos
fgev.lmom <- function(x, rp){
  par <- pelgev(samlmu(x, nmom = 3))
  tibble(
    return_period = rp,
    intensity = quagev(1 - 1/rp, par),
    xi = par[1],
    alpha = par[2],
    kappa = par[3]
  )
}

fitsgb.idf <- function(data, objective = c("REQM", "DPMA"), b.interval = c(0, 120)){

  objective <- match.arg(objective, c("REQM", "DPMA"))

  # Otimizar `b`
  b0 <- optimize(
    function(b){
      r2 <- data |>                # minimizar R²
        group_by(return_period) |> # agrupar por tempo de retorno
        summarise(
          r2 = summary(lm(log(intensity) ~ log(duration_min + b)))$r.squared, # modelo log(intensidade) ~ log(duração + b)
          .groups = "drop"
        ) |> 
        summarise(r2 = mean(r2)) |> 
        pull(r2)
      -r2  # negativar R²
    },
    interval = b.interval
  )$minimum

  initial.by.rp <- data |> # estimativa inicial por tempo de retorno
    group_by(return_period) |>  # agrupar por tempo de retorno
    group_modify(~{
      m <- lm(log(intensity) ~ log(duration_min + b0), data = .x)
      tibble(
        A = exp(unname(coef(m)[1])),
        c = -unname(coef(m)[2])
      )
    }) |> 
    ungroup()

  mA <- lm(log(A) ~ log(return_period), data = initial.by.rp)

  # Parâmetros iniciais 
  par.start <- c(
    log(exp(unname(coef(mA)[1]))),
    unname(coef(mA)[2]),
    log(max(b0, 1e-6)), # não permitir `b` < 0
    log(mean(initial.by.rp$c))
  )

  # Função objetivo p/ ajuste final (REQM ou DPMA)
  fo <- function(par){

    k <- exp(par[1])
    gamma <- par[2]
    b <- exp(par[3])
    c <- exp(par[4])

    calc <- with(data, k*return_period^gamma/(duration_min + b)^c)

    switch(objective,
      "REQM" = sqrt(sum((data$intensity - calc)^2)/(nrow(data) - 1)),
      "DPMA" = mean(abs(data$intensity - calc)/data$intensity)*100
    )

  } # fim `fo`

  opt <- optim(
    par = par.start, # parâmetros iniciais
    fn = fo,         # função objetivo
    method = "Nelder-Mead",
    control = list(maxit = 50000, reltol = 1e-12)
  )

  # Parâmetros finais
  par <- c(
    k = exp(opt$par[1]),
    gamma = opt$par[2],
    b = exp(opt$par[3]),
    c = exp(opt$par[4])
  )

  fitted <- data |> mutate(fitted = par["k"]*return_period^par["gamma"]/(duration_min + par["b"])^par["c"])

  tibble(
    k = par["k"],
    gamma = par["gamma"],
    b = par["b"],
    c = par["c"],
    reqm = sqrt(sum((fitted$intensity - fitted$fitted)^2)/(nrow(fitted) - 1)),
    dpma = mean(abs(fitted$intensity - fitted$fitted)/fitted$intensity)*100,
    convergence = opt$convergence
  )

}


# CÁLCULO ----------------------------------------------------------------

# Ler dados
na.yr <- 0.2
min.durations <- 3
d.interval <- c(1, 24)
min.yr <- 8
imax.df <- read_parquet("base/gerados/df_imax.pqt") |> 
  filter(
    na_prct_yr <= na.yr,
    between(d, d.interval[1], d.interval[2])
  ) |> 
  group_by(gauge_code) |> 
  mutate(
    n_years = n_distinct(wateryear),
    gauge_label = sprintf("%s | %2g anos", gauge_code, n_years)
  ) |> 
  filter(
    n_distinct(d) >= min.durations,
    n_years >= min.yr
  ) |> 
  ungroup()

# Ajustar GEV por estação e duração
rp <- c(2, 5, 10, 15, 20, 25, 30, 40, 50, 60, 75, 100)
idf.df <- imax.df |> 
  mutate(duration_min = d*60) |> 
  group_by(gauge_code, d, duration_min) |> 
  group_modify(~fgev.lmom(.x$imax, rp)) |> 
  ungroup()

# Informações
station.info <- imax.df |> 
  group_by(gauge_code) |> 
  summarise(
    n_years = first(n_years),
    first_year = min(wateryear),
    last_year = max(wateryear),
    .groups = "drop"
  )

# Extrair tabela somente c/ parâmetros da GEV
gev.par.df <- idf.df |> 
  distinct(gauge_code, d, duration_min, xi, alpha, kappa)

# Ajuste da IDF p/ cada estação usando REQM
idf.sgb.df <- idf.df |> 
  group_by(gauge_code) |> 
  group_modify(~fitsgb.idf(.x, objective = "REQM")) |> 
  ungroup() |> 
  left_join(station.info, by = "gauge_code") |> 
  mutate(dpma_lim_10 = dpma <= 10) |> 
  select(gauge_code, n_years, first_year, last_year, k, gamma, b, c, reqm, dpma, dpma_lim_10, convergence)