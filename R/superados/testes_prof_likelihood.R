rm(list = ls()); invisible(gc())
pacman::p_load(pacman, tidyverse)

# Set theme
suppressMessages(extrafont::loadfonts())
font <- "Consolas"
fontsize <- 8
textsize <- (fontsize - 2) / .pt
ggplot2::set_theme(
  theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box.margin = margin(-5, 2, 2, 2),
      legend.key.size = unit(fontsize, "pt"),
      text = element_text(color = "black", size = fontsize, family = font),
      axis.text = element_text(color = "black"),
      strip.text = element_text(hjust = 0, size = fontsize - 2, margin = margin(2, 2, 2, 2), color = "black"),
      panel.border = element_rect(color = "black"),
      panel.grid = element_blank()
    )
)

# Importar funções
source("R/funcoes/fdgev.R") # ajuste com invariância de escala
source("R/funcoes/fgev.R")  # ajuste por duração

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
# prior.info <- c(0.0178, 0.0627)
conf <- 0.95
data.gg <- data[data$gauge_code == "A812", ]
cols <- names(data)
maxit <- 1e8

y.id <- data.gg %>%
  pivot_wider(names_from = d, values_from = imax, names_prefix = "d_") %>%
  arrange(wateryear) %>%
  select(-c(gauge_code, wateryear)) %>%
  as.matrix() %>%
  na.omit()


# AJUSTE GEV -------------------------------------------------------------

# Função para otimização com theta_{-i}
prof.nll.gev <- function(par.free, which.par, par.value, y, prior.info){
  i <- which.par       # posição do parâmetro que vai ser analisado (de 1:length(param))
  full <- numeric(3)
  full[i] <- par.value # substituir parâmetro em análise
  full[-i] <- par.free # trazer demais parâmetros que serão otimizados
  nll.gev(full, y, prior.info)
}

# Função gradiente com theta_{-i}
prof.grad.gev <- function(par.free, which.par, par.value, y, prior.info){
  i <- which.par # posição do parâmetro que vai ser analisado (de 1:length(param))
  full <- numeric(3)
  full[i] <- par.value
  full[-i] <- par.free
  gr <- grad.ll.gev(full, y, prior.info)
  gr[-i]
}

# Limites inferior e superior
lower <- c(-Inf, 1e-6, if (is.null(prior.info)) -Inf else -0.5 + 1e-6)
upper <- c(Inf, Inf, if (is.null(prior.info)) Inf else 0.5 - 1e-6)

# Ajuste com série original
d <- 20 # no final usar apply() não vai precisar de índice
y <- y.id[,d]
fit <- fgev(y, prior.info)
par <- fit$par                     # original parameters
nll.min <- fit$value               # observed minimum

i <- 3

# Output matrix
fctr <- 6
del <- 0.0025
npar <- length(par)
par.range <- seq(-fctr*abs(par[i]), fctr*abs(par[i]), del)
prof.par <- matrix(NA_real_, nrow = length(par.range) + 1, ncol = length(par) + 4)
colnames(prof.par) <- c(names(par), "loglik", "failure_lbfgsb", "failure_nm", "jump")
prof.par[1,] <- c(par, -nll.min, 0, 0, NA_real_)
# j <- 2

for(j in 2:length(par.range)){

  flag.lbfgsb <- 1 # assumir que L-BFGS-B falhou
  flag.nm <- 1     # assumir que Nelder-Mead falhou
  val.i <- par.range[j]
  last.par <- prof.par[j-1, seq_len(npar)]
  last.ll <- prof.par[j-1, npar + 1]
  last.jump <- prof.par[j-1, ncol(prof.par)]
  fit.prof <- tryCatch(
    optim(
      last.par[-i], prof.nll.gev, gr = prof.grad.gev, method = "L-BFGS-B", lower = lower, upper = upper, 
      which.par = i, par.value = val.i, prior.info = prior.info, y = y,
      control = list(maxit = maxit, factr = 1e9, pgtol = 1e-8)
    ),
    error = function(e) NULL
  )

  if(!is.null(fit.prof) && fit.prof$convergence == 0){
    ll <- -fit.prof$value
    flag.lbfgsb <- 0
    flag.nm <- NA_real_
  } else {
    fit.prof <- tryCatch(
      optim(
        last.par[-i], prof.nll.gev, gr = NULL, method = "Nelder-Mead",
        which.par = i, par.value = val.i, prior.info = prior.info, y = y,
        control = list(maxit = maxit)
      ),
      error = function(e) NULL
    )
    if(!is.null(fit.prof) && fit.prof$convergence == 0){
      ll <- -fit.prof$value
      flag.nm <- 0
      if(ll < -1e6){
        ll <- NA_real_
        flag.nm <- 1
      }
    } else{
      ll <- NA_real_
    }
  }

  if(!is.na(ll) && !is.na(last.ll)){
    jump <- (ll - last.ll)/del
  } else{
    jump <- NA_real_
  }
  
  full <- numeric(npar)
  if(!is.null(fit.prof) && !is.na(ll)){
    full[-i] <- fit.prof$par # se não falhar, atualizar parâmetros
  } else{
    full[-i] <- last.par[i]  # se falhar, manter rodada anterior
  }
  full[i] <- val.i
  prof.par[j,] <- c(full, ll, flag.lbfgsb, flag.nm, jump)

}

prof.par <- as_tibble(prof.par)

lim <-  -nll.min - qchisq(conf, 1)/2 # threshold for (1 - alpha)*100 confidence
# lim <-  max(prof.par[["loglik"]], na.rm = TRUE) - qchisq(conf, 1)/2 # threshold for (1 - alpha)*100 confidence

ggplot() +
  geom_vline(xintercept = par[i], linetype = "dashed", color = "grey50") +
  geom_line(aes(x = prof.par[[i]], y = prof.par[["loglik"]]), linewidth = 0.8) +
  geom_hline(yintercept = lim, color = "red") +
  geom_text(
    aes(x = par[i], y = lim, label = sprintf("κ = %0.3f", par[i])), 
    family = font, size = textsize, angle = 90, vjust = -1.1, hjust = 1.2
  ) +
  labs(
    x = names(par)[i], y = "Perfil de verossimilhança", 
    subtitle = sprintf("Duração: %g h | Priori: %s | Delta: %0.03f", d, gsub("c", "", deparse(prior.info)), del)
  ) +
  theme(panel.grid = element_blank())

# AJUSTE d-GEV -----------------------------------------------------------

fit <- fdgev(y.id, prior.info, cols = names(data.gg), method = "L-BFGS-B", scale = "wide" )

par <- fit$par
nll.min <- fit$value
se <- tryCatch(sqrt(diag(solve(fit$hessian))), error = function(e) rep(NA, 5))

durations <- unique(data.gg[[cols[3]]])

# Limites inferior e superior
lower <- c(-Inf, 1e-6, if (is.null(prior.info)) -Inf else -0.5 + 1e-6, 1e-6, 0)
upper <- c(Inf, Inf, if (is.null(prior.info)) Inf else 0.5 - 1e-6, 1 - 1e-6, Inf)

prof.nll.dgev <- function(par.free, which.par, par.value, durations, y.id, prior.info){
  i <- which.par       # posição do parâmetro que vai ser analisado (de 1:length(param))
  full <- numeric(5)   # build back full set of parameters
  full[i] <- par.value # substituir parâmetro em análise
  full[-i] <- par.free # trazer demais parâmetros que serão otimizados
  nll.dgev(full, durations, y.id, prior.info)
}

prof.grad.dgev <- function(par.free, which.par, par.value, durations, y.id, prior.info){
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
lim <- nll.min + qchisq(conf, 1) / 2 # máximo observado - quantil qui-quadrado/2
result <- vector("list", 5)
names(result) <- par.names

# Teste com kappa
par.range <- seq(-0.5, 0.5, 0.005) # usando -0.5 e 0.5 como limites
i <- 3 # kappa
j <- 1

ls.prof.par <- pbapply::pblapply(par.range, function(val.i){

  fit <- tryCatch(
    optim(
      par[-i], prof.nll.dgev, gr = prof.grad.dgev, method = "L-BFGS-B", lower = lower, upper = upper, 
      which.par = i, par.value = val.i, durations = durations, prior.info = prior.info, y.id = y.id,
      control = list(maxit = maxit, factr = 1e9, pgtol = 1e-8)
    ),
    error = function(e) NULL
  )

  if(!is.null(fit) && fit$convergence == 0){
    ll <- -fit$value
  } else {
    fit <- tryCatch(
      optim(
        par[-i], prof.nll.dgev, gr = NULL, method = "Nelder-Mead",
        which.par = i, par.value = val.i, durations = durations, prior.info = prior.info, y.id = y.id,
        control = list(maxit = maxit)
      ),
      error = function(e) NULL
    )
    if(!is.null(fit) && fit$convergence == 0){
      ll <- -fit$value
      if(ll < -1e6) ll <- NA_real_
    } else{
      ll <- NA_real_
    }
  }

  return(data.frame(parname = names(par)[i], par = val.i, loglik = ll))

})

prof.par <- do.call(rbind, ls.prof.par)

ggplot() +
  geom_segment(aes(x = par[i], y = min(prof.par$loglik, na.rm = TRUE), yend = -nll.min), linetype = "dashed", color = "grey50") +
  geom_line(aes(x = par, y = loglik), prof.par, linewidth = 0.8) +
  geom_hline(yintercept = -lim) +
  labs(x = names(par)[i], y = "Perfil de verossimilhança") +
  # scale_y_continuous(limits = c(-1.05*lim, 0)) +
  theme(panel.grid = element_blank())

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
    bounds <- list(lower = lower[-i], upper = upper[-i])

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
