rm(list = ls())
invisible(gc())

# Importar funções
source("R/funcoes/fdgev.R")

# Intensidades máximas anuais
df.imax <- arrow::read_parquet("base/gerados/df_imax.pqt")

# Aplicar critérios de qualidade
data <- df.imax %>%
  filter(na_prct_yr <= 0.2, between(d, 1, 24)) %>%
  group_by(gauge_code) %>%
  filter(n_distinct(d) >= 3, n_distinct(wateryear) >= 8) %>%
  ungroup() %>%
  select(gauge_code, wateryear, d, imax)

prior.info <- NULL
conf <- 0.95
n.grid <- 21
se.factor <- 4
data.gg <- data[data$gauge_code == "A812", ]
cols <- names(data)

fit <- fdgev(
  data.gg,
  prior.info,
  cols = names(data.gg),
  method = "L-BFGS-B",
  scale = "wide"
)

par <- fit$par
nll.min <- fit$value
se <- tryCatch(sqrt(diag(solve(fit$hessian))), error = function(e) rep(NA, 5))

y.id <- data.gg %>%
  pivot_wider(names_from = d, values_from = imax, names_prefix = "d_") %>%
  arrange(wateryear) %>%
  select(-c(gauge_code, wateryear)) %>%
  as.matrix() %>%
  na.omit()

durations <- unique(data.gg[[cols[3]]])

# Limites inferior e superior
lower <- c(-Inf, 1e-6, if (is.null(prior.info)) -Inf else -0.5 + 1e-6, 1e-6, 0)
upper <- c(
  Inf,
  Inf,
  if (is.null(prior.info)) Inf else 0.5 - 1e-6,
  1 - 1e-6,
  Inf
)

grad.ll.dgev <- function(param, durations, y.id, prior.info) {
  # Parameters
  xit <- param[1]
  alpha0 <- param[2]
  kappa <- param[3]
  H <- param[4]
  theta <- param[5]

  ds <- durations
  N <- nrow(y.id)
  D <- ncol(y.id)

  if (any(ds + theta <= 0) || alpha0 <= 0 || H < 0 || H > 1 - 1e-6) {
    return(rep(1e6, 5))
  }

  # Variáveis auxiliares
  A <- (ds + theta)^H
  yA <- sweep(x = y.id, MARGIN = 2, STATS = A, FUN = "*") # produto: y.id * Ad
  B <- yA / alpha0 - xit # y.id*A/alpha0 - xit
  z.id <- 1 - kappa * B

  if (any(z.id <= 0)) {
    return(rep(1e6, 5))
  }

  # Pre-calcular termos com z.id
  z.inv <- 1 / z.id
  z.pow <- z.id^(1 / kappa - 1)
  ln.z <- log(z.id)
  z.pow.k <- z.id^(1 / kappa)

  # Derivadas parciais
  # Posição (xit)
  dl.dxit <- (1 - kappa) * sum(z.inv) - sum(z.pow)

  # Escala (alpha0)
  dl.dalpha0 <- -(N * D) /
    alpha0 +
    (1 - kappa) * sum(yA * z.inv / (alpha0^2)) -
    sum(z.pow * yA / (alpha0^2))

  # Forma (kappa)
  if (!is.null(prior.info)) {
    mu.kappa <- prior.info[1]
    var.kappa <- prior.info[2]^2

    a <- 0.5
    b <- mu.kappa + a
    p <- b^2 * (1 - b) / var.kappa - b
    q <- p * (1 / b - 1)
    if (a - abs(kappa) <= 1e-6) {
      return(rep(1e6, 5))
    }

    dl.prior <- (p - 1) / (a + kappa) - (q - 1) / (a - kappa)
  } else {
    dl.prior <- 0
  }

  term1 <- (1 / kappa - 1) * B * z.inv + ln.z / (kappa^2)
  term2 <- (1 / kappa) * ((ln.z * z.pow.k) / kappa + B * z.pow)
  dl.dkappa <- -sum(term1) + sum(term2) + dl.prior

  # Expoente de invariância de escala (H)
  ln.dt <- log(ds + theta)
  aux1 <- sweep(yA / alpha0, MARGIN = 2, ln.dt, FUN = "*")
  dl.dH <- N * sum(ln.dt) + (kappa - 1) * sum(z.inv * aux1) + sum(z.pow * aux1)

  # Deslocamento (theta)
  dt.inv <- 1 / (ds + theta)
  aux2 <- sweep(yA / alpha0, MARGIN = 2, STATS = H * dt.inv, FUN = "*")
  dl.dtheta <- N *
    H *
    sum(dt.inv) +
    (kappa - 1) * sum(z.inv * aux2) +
    sum(z.pow * aux2)

  par.names <- c("xit", "alpha0", "kappa", "H", "theta")
  grad <- setNames(
    c(dl.dxit, dl.dalpha0, dl.dkappa, dl.dH, dl.dtheta),
    par.names
  )
  return(-grad)
}

nll.dgev <- function(param, durations, y.id, prior.info) {
  # Parâmetros iniciais
  xit <- param[1] # posição: xi~
  alpha0 <- param[2] # escala: alpha0
  kappa <- param[3] # forma
  H <- param[4] # invariância de escala (expoente de escala - H ou eta)
  theta <- param[5] # invariância de escala (parâmetro de posição - theta)

  N <- nrow(y.id)
  D <- ncol(y.id)

  # Limites teóricos
  if (any(durations + theta <= 0)) {
    return(1e6 + abs(min(durations + theta)) * 1e4)
  } # impedir offset negativo
  if (theta < 0) {
    return(1e6 + abs(theta) * 1e4)
  } # impedir offset negativo
  if (alpha0 <= 0) {
    return(1e6 + abs(alpha0) * 1e4)
  } # impedir escala GEV nulo
  if (H < 1e-6 || H > 1 - 1e-6) {
    return(1e6 + abs(H) * 1e4)
  } # restrição invariância de escala

  # Calcular termos z.id = 1 - kappa*(y.id*(d + theta)^H)/alpha - xi~)
  # Função sweep() realiza uma operação nas linhas ou colunas de euma matriz
  A <- (durations + theta)^H # termo de alpha(d)
  aux <- sweep(
    x = y.id / alpha0, # matriz c/ imax
    MARGIN = 2, # operação por colunas
    STATS = A, # argumentos (nesse caso oq vai ser multiplicado por y.id)
    FUN = "*"
  ) # multiplicação

  z.id <- 1 - kappa * (aux - xit) # variável reduzida vetorizada

  # Conferir limites teóricos
  # if(any(z.id <= 0)) return(1e8)
  if (any(z.id <= 1e-6)) {
    return(1e6 + abs(min(z.id)) * 1e4)
  }

  # Definir priori log-beta caso tenha sido informada
  if (!is.null(prior.info)) {
    mu.kappa <- prior.info[1]
    var.kappa <- prior.info[2]^2

    a <- 0.5
    b <- mu.kappa + a
    p <- b^2 * (1 - b) / var.kappa - b
    q <- p * (1 / b - 1)

    if (a - abs(kappa) <= 1e-6) {
      return(1e6 + abs(kappa) * 1e4)
    } # impedir log(<0)
    # if(a - abs(kappa) <= 1e-6) return(1e6) # impedir log(<0)
    l.prior <- (p - 1) * log(a + kappa) + (q - 1) * log(a - kappa)
  } else {
    l.prior <- 0
  }

  # Log-verossimilhança
  nll <- -(-N *
    D *
    log(alpha0) +
    H * N * sum(log(durations + theta)) +
    (1 / kappa - 1) * sum(log(z.id)) -
    sum(z.id^(1 / kappa)) +
    l.prior)

  return(nll) # retornando negativo p/ "minimização"
}

prof.nll <- function(
  par.free,
  which.par,
  par.value,
  durations,
  y.id,
  prior.info
) {
  i <- which.par # posição do parâmetro que vai ser analisado (de 1:length(param))
  full <- numeric(5)
  full[i] <- par.value # substituir parâmetro em análise
  full[-i] <- par.free # trazer demais parâmetros que serão otimizados
  nll.dgev(full, durations, y.id, prior.info)
}

prof.grad <- function(
  par.free,
  which.par,
  par.value,
  durations,
  y.id,
  prior.info
) {
  i <- which.par # posição do parâmetro que vai ser analisado (de 1:length(param))
  full <- numeric(5)
  full[i] <- par.value
  full[-i] <- par.free
  gr <- grad.ll.dgev(full, durations, y.id, prior.info)
  gr[-i]
}

# par.names <- c("xit", "alpha0", "kappa", "H", "theta")
par.names <- names(par)

# Limite baseado no quantil da qui-quadrado
lim <- nll.min + qchisq(conf, 1) / 2
result <- vector("list", 5)
names(result) <- par.names

for (i in length(par)) {
  # Criar uma grade ao longo da função de verossimilhança
  if (is.na(se[i]) || se[i] <= 0) {
    step <- ifelse(abs(par[i]) < 1e-6, 0.1, abs(par[i] * 0.1))
    grid <- seq(
      par[i] - se.factor * step,
      par[i] + se.factor * step,
      length.out = n.grid
    )
  } else {
    grid <- seq(
      par[i] - se.factor * se[i],
      par[i] + se.factor * se[i],
      length.out = n.grid
    )
  }

  # Conferir limites
  grid <- grid[grid >= lower[i] & grid <= upper[i]]

  # Iterar pelo perfil
  proflik.values <- numeric(length(grid))
  for (j in seq_along(grid)) {
    val <- grid[j] # valor da grid
    start <- par[-i] # remover o parâmetro em análise do vetor
    bounds <- list(lower = lower[i], upper = upper[-i])

    # Otimizar
    fit <- tryCatch(
      optim(
        start,
        prof.nll,
        gr = prof.grad,
        method = "L-BFGS-B",
        which.par = i,
        par.value = val,
        durations = durations,
        y.id = y.id,
        prior.info = prior.info,
        control = list(maxit = 1e8, factr = 1e9)
      ),
      error = function(e) NULL
    )
    if (is.null(fit) || fit$convergence != 0) {
      fit <- tryCatch(
        optim(
          start,
          prof.nll,
          gr = NULL,
          method = "Nelder-Mead",
          which.par = i,
          par.value = val,
          durations = durations,
          y.id = y.id,
          prior.info = prior.info,
          control = list(maxit = 1e8)
        ),
        error = function(e) NULL
      )
    }
    proflik.values[j] <- if (is.na(fit)) NA else fit$value
  }

  result[[i]] <- list(grid = grid, nll = proflik.values, lower = ci_l)
}

ci.proflik <- function() {}
