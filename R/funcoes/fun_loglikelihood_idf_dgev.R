# Construção da função de log-verossimilhança do modelo de equação de IDF baseado
# na distribuição dGEV (duration-dependent Generalized Extreme Value distribution)

# Provavelmente, vai precisar fazer uma organização dos dados igual na função de
# invariância de escala, organizando em cada coluna a AMS de uma duração diferente
# pq a função vai precisar iterar em um loop das durações e das intensidades
# então, entra na primeira coluna, extrai o primeiro 'd', e os 'y_i' pro somatório

rm(list = ls()); invisible(gc())

# Set theme
suppressMessages(extrafont::loadfonts())
font <- "Consolas"
fontsize <- 8
textsize <- (fontsize - 2)/.pt
set_theme(theme_minimal() +
            theme(text = element_text(color = "black", size = fontsize, family = font),
                  axis.text = element_text(color = "black"),
                  panel.border = element_rect(color = "black")))

pacman::p_load(pacman, tidyverse, arrow, numDeriv, evd)

df.imax <- arrow::read_parquet("base/gerados/df_imax.pqt")

data <- df.imax %>% 
  filter(na_prct_yr <= 0.2, na_prct_mon <= 0.1, between(d, 1, 24)) %>% 
  group_by(gauge_code) %>% 
  filter(n_distinct(d) >= 3, n_distinct(wateryear) >= 8) %>% 
  ungroup() %>% 
  select(gauge_code, wateryear, d, imax)

gauges <- unique(data$gauge_code)
gg <- "A801"
gg <- "430350901A"
gg <- "430510802A"
gg <- "SBPA"
gg <- "A884" # falha na otimização
gg <- "431020703A"
gg <- "A803"

data.gg <- data %>% 
  filter(gauge_code == gg)

cols <- c("wateryear", "d", "imax")

teste <- data.gg %>% 
  pivot_wider(names_from = cols[1],
              values_from = cols[2],
              names_prefix = "d_")

y.id <- data.gg %>% 
  pivot_wider(names_from = d,
              values_from = imax,
              names_prefix = "d_") %>% 
  arrange(wateryear) %>% 
  select(-c(gauge_code, wateryear)) %>% 
  as.matrix() %>% 
  na.omit()

durations <- as.numeric(stringr::str_remove_all(colnames(y.id), "d_"))


# FUNÇÃO ------------------------------------------------------------------


## nll.dgev ---------------------------------------------------------------

#' @title Função de log-verossimilhança do modelo de d-GEV 
#' @param param initial parameters: 1. location, 2. scale, 3. shape, 4. eta (power scale inv), 5. theta (offset scale inv)
#' @param durations durations: same as colnames(y.id) [hours]
#' @param y.id matrix with a column for each duration and annual intensity maximas on each row

nll.dgev <- function(param, durations, y.id){
  
  # Parâmetros iniciais
  xit <- param[1]    # posição: xi~
  alpha0 <- param[2] # escala: alpha0
  kappa <- param[3]  # forma
  H <- param[4]      # invariância de escala (expoente de escala - H ou eta)
  theta <- param[5]  # invariância de escala (parâmetro de posição - theta)
  
  N <- nrow(y.id)
  D <- ncol(y.id)
  
  # Limites teóricos
  if(any(durations + theta <= 0)) return(1e6) # impedir offset negativo
  if(theta <= 0) return(1e6)                  # impedir offset negativo
  if(alpha0 <= 0) return(1e6)                 # impedir escala GEV nulo
  if(H <= 0 | H > 1) return(1e6)              # restrição invariância de escala
  
  # Calcular termos z.id = 1 - kappa*(y.id*(d + theta)^H)/alpha - xi~)
  # Função sweep() realiza uma operação nas linhas ou colunas de euma matriz
  A <- (durations + theta)^H      # termo de alpha(d)
  aux <- sweep(x = y.id/alpha0,   # matriz c/ imax
               MARGIN  = 2,       # operação por colunas
               STATS = A,         # argumentos (nesse caso oq vai ser multiplicado por y.id)
               FUN = "*")         # multiplicação
  
  z.id <- 1 - kappa*(aux - xit) # variável reduzida vetorizada
  
  # Conferir limites teóricos
  if(any(z.id <= 0)) return(1e6)
  # return(1e10)
  
  # Log-verossimilhança
  nll <- -(-N*D*log(alpha0) + H*N*sum(log(durations + theta)) + (1/kappa - 1)*sum(log(z.id)) - sum(z.id^(1/kappa)))
  
  return(nll) # retornando negativo p/ "minimização"
  
}

# Função gradiente log-verossimilhança
grad.ll.dgev <- function(param, durations, y.id){
  
  # Parameters
  xit <- param[1]
  alpha0 <- param[2]
  kappa <- param[3]
  H <- param[4]
  theta <- param[5]
  
  ds <- durations
  N <- nrow(y.id)
  D <- ncol(y.id)
  
  if(any(ds + theta <= 0) || alpha0 <= 0 || H <= 0 || H > 1){
    return(rep(1e6, 5))
  }
  
  # Variáveis auxiliares
  A <- (ds + theta)^H
  yA <- sweep(x = y.id, MARGIN = 2, STATS = A, FUN = "*") # produto: y.id * Ad
  B <- yA/alpha0 - xit                                    # y.id*A/alpha0 - xit
  z.id <- 1 - kappa*B
  
  if(any(z.id <= 0)) return(rep(1e6, 5))
  
  # Pre-calcular termos com z.id
  z.inv <- 1/z.id
  z.pow <- z.id^(1/kappa - 1)
  ln.z <- log(z.id)
  z.pow.k <- z.id^(1/kappa)
  

  ## DERIVADAS PARCIAIS ---------------------------------------------------
  
  # Posição (xit)
  dl.dxit <- (1 - kappa)*sum(z.inv) - sum(z.pow)
  
  # Escala (alpha0)
  dl.dalpha0 <- -(N*D)/alpha0 + (1 - kappa)*sum(yA*z.inv/(alpha0^2)) - sum(z.pow*yA/(alpha0^2))
  
  # Forma (kappa)
  term1 <- (1/kappa - 1)*B*z.inv + ln.z/(kappa^2)
  term2 <- (1/kappa)*((ln.z*z.pow.k)/kappa + B*z.pow)
  dl.dkappa <- -sum(term1) + sum(term2)
  
  # Expoente de invariância de escala (H)
  ln.dt <- log(ds + theta)
  aux1 <- sweep(yA/alpha0, MARGIN = 2, ln.dt, FUN = "*")
  dl.dH <- N*sum(ln.dt) + (kappa - 1)*sum(z.inv*aux1) + sum(z.pow*aux1)
  
  # Deslocamento (theta)
  dt.inv <- 1/(ds + theta)
  aux2 <- sweep(yA/alpha0, MARGIN = 2, STATS = H*dt.inv, FUN = "*")
  dl.dtheta <- N*H*sum(dt.inv) + (kappa - 1)*sum(z.inv*aux2) + sum(z.pow*aux2)
  
  par.names <- c("xit", "alpha0", "kappa", "H", "theta")
  grad <- setNames(c(dl.dxit, dl.dalpha0, dl.dkappa, dl.dH, dl.dtheta), par.names)
  return(-grad)
  
}


# Estimar parâmetros iniciais
start.dgev <- function(data.gg, durations, scale.inv = c("wide", "scrict")){
  
  data.imax <- data.gg[,c("d", "imax", "wateryear")]
  
  if(!identical(unique(data.imax$d), durations)) stop("Durations in 'durations' vector  and 'data.gg$d' do not match.")
  scale.inv <- match.arg(scale.inv, c("wide", "scrict"))
  
  y.id <- data.imax %>% 
    pivot_wider(names_from = d,
                values_from = imax,
                names_prefix = "d_") %>%
    select(-wateryear) %>% 
    as.matrix() %>% 
    na.omit()
  
  # Expoente 'scale' de invariância de escala
  scale <- switch(
    scale.inv,
    wide = {
      mom <- apply(y.id, 2, function(imax) mean(imax))
      scale <- -lm(log(mom) ~ log(durations))$coefficient[[2]]
    },
    scrict = {
      scale.xi <- -lm(log(param.d["alpha",]) ~ log(durations))$coefficient[[2]]
      scale.alpha <- -lm(log(param.d["alpha",]) ~ log(durations))$coefficient[[2]]
      scale <- mean(scale.xi, scale.alpha)
    }
  )
  
  theta <- 0
  
  # Estimar parâmetros da GEV p/ cada duração
  par.d <- apply(y.id, 2, function(imax){
    sample.lmom <- lmom::samlmu(imax)
    par <- lmom::pelgev(sample.lmom)
  })
  
  alpha <- par.d[2,1]                          # estimate for first duration
  xit <- par.d[1,1]/alpha                      # rescaled xi
  alpha0 <- alpha*(durations[1] + theta)^scale # rescaled alpha
  kappa <- mean(par.d[3,])                     # average kappa
  
  names <- c("xit", "alpha0", "kappa", "scale", "theta")
  
  return(setNames(c(xit, alpha0, kappa, scale, theta), names))
  
}

# start <- start.dgev(data.gg, durations, scale.inv = "scrict")
start <- start.dgev(y.id, durations, scale.inv = "wide")

method <- "Nelder-Mead"
# method <- "L-BFGS-B"
maxit <- 1e9

# Otimização com funções gradiente
fit <- optim(par = start, fn = nll.dgev, hessian = TRUE, gr = grad.ll.dgev,
             method = method, control = list(maxit = maxit),
             y.id = y.id, durations = durations)

par <- fit$par %>% print()

convergence <- fit$convergence
message <- fit$message
if(convergence == 0){
  message("Otimização concluída!")
} else{
  message("Houve um erro durante o processo de otimização com método ", method)
  message(convergence, ": ", message)
}

ls.data <- split(data, data$gauge_code)
teste.nprior <- pbapply::pblapply(ls.data, fit.dgev, cols = names(data), method = "Nelder-Mead", scale.inv = "strict") # só L-BFGS-B está funcionando
teste.nprior <- pbapply::pblapply(ls.data, fit.dgev, cols = names(data), method = "L-BFGS-B", scale.inv = "strict") # só L-BFGS-B está funcionando
# teste <- pbapply::pblapply(ls.data, fit.dgev, 
teste.wprior <- pbapply::pblapply(ls.data, fit.dgev,
                           prior.info = c(0.0178, 0.0627),
                           cols = names(data), 
                           # method = "Nelder-Mead", 
                           method = "L-BFGS-B",
                           scale.inv = "strict") # só L-BFGS-B está funcionando

par.nprior.gg <- teste.nprior[[gg]]$par %>% print()
par.wprior.gg <- teste.wprior[[gg]]$par %>% print()

df.teste <- bind_rows(
  lapply(seq_along(teste), function(i){
    fit <- teste[[i]]
    gg <- names(teste[i])
    row <- tibble(gauge_code = gg, par_name = names(fit$par), par = fit$par, convergence = fit$convergence, message = fit$message,
                  value = fit$value, hessian = list(as_tibble(fit$hessian)))
  })
)


# QUANTIS -----------------------------------------------------------------

# fit <- teste[[gg]]
fit <- teste.nprior[[gg]]
fit <- teste.wprior[[gg]]
par <- fit$par %>% print()
# ds <- as.numeric(gsub("d+_","", colnames(y.id)))
durations <- c(1, 2, 4, 6, 8, 12, 24)

q.dgev <- function(param, duration, p.nexceed){
  
  param <- unname(param)
  xit <- param[1]
  alpha0 <- param[2]
  kappa <- param[3]
  H <- param[4]
  theta <- param[5]
  
  p <- p.nexceed; d <- duration
  
  # Calcular parâmetros convencionais da GEV
  alpha <- alpha0/(d + theta)^H
  xi <- xit*alpha
  
  yp <- xi + alpha/kappa*(1 - (-log(p))^kappa)
  
  return(yp)
  
}

data.gg <- data[data$gauge_code == gg,]
y.id.gg <- data.gg %>% 
  pivot_wider(names_from = cols[3],
              values_from = cols[4],
              names_prefix = "d_") %>% 
  arrange(wateryear) %>% 
  select(-cols[1:2]) %>% 
  na.omit() %>% 
  as.matrix()

# Posiçao de plotagem de Gringorten
freq.obs <- apply(y.id.gg, 2, function(ams){
  sorted <- sort(ams, decreasing = TRUE)
  n <- length(ams)
  i <- 1:length(ams)
  p <- (i - 0.44)/(n + 0.12)
})

rp <- c(5, 10, 50, 100, 200)
p <- 1 - 1/rp
rp <- 1/freq.obs[,1]
quantiles <- lapply(durations, q.dgev, p.nexceed = p, param = par)
quantiles <- lapply(durations, q.dgev, p.nexceed = 1 - freq.obs[,1], param = par) # passar prob. não excedência
idf.aux <- matrix(unlist(quantiles), length(freq.obs[,1]), length(durations), dimnames = list(rp, durations))
idf.aux <- matrix(unlist(quantiles), length(p), length(durations), dimnames = list(rp, durations))
idf <- bind_cols(gauge_code = gg, idf.aux) %>% 
  mutate(return_period = rp) %>% 
  pivot_longer(cols = -c(gauge_code, return_period),
               names_to = "d",
               names_transform = list(d = as.numeric),
               values_to = "q") %>% print()

# Frequencias observadas e valores
df.fq.obs <- y.id.gg %>% 
  as_tibble() %>% 
  pivot_longer(cols = where(is.numeric),
               names_prefix = "d_",
               names_to = "d",
               names_transform = list(d = as.numeric),
               values_to = "imax") %>% 
  group_by(d) %>% 
  arrange(desc(imax)) %>% 
  mutate(n = length(imax), i = row_number(),
         p_weibull = i/(n + 1), tr_weibull = 1/p_weibull,
         p_cunanne = (i - 0.4)/(n + 0.2), tr_cunanne = 1/p_cunanne,
         p_gring = (i - 0.44)/(n + 0.12), tr_gring = 1/p_gring) %>% 
  ungroup() %>% print()


# VISUALIZATION -----------------------------------------------------------

# Plotar curva de frequência
ggplot() +
  geom_line(aes(x = return_period, y = q, color = factor(d)), idf %>% mutate(return_period = as.numeric(return_period))) +
  geom_point(aes(x = tr_gring, y = imax, color = factor(d)), df.fq.obs %>% filter(d %in% durations)) +
  scale_x_log10() +
  labs(x = "T [years]", y = "Intensity [mm/h]", color = "Duration [h]", caption = sprintf("Subdaily gauge: %s", gg)) +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.key.height = unit(fontsize, "pt"),
        legend.title = element_text(size = fontsize - 2))

# Plot IDF
ggplot() +
  geom_line(aes(x = d, y = q, color = factor(return_period)), idf, linewidth = 0.5) +
  geom_point(aes(x = d, y = imax), data.gg[data.gg$d %in% durations,], alpha = 0.4, color = "grey10", size = 0.8) +
  scale_x_log10(breaks = durations) +
  # scale_y_log10() +
  labs(x = "Duration [h]", y = "Intensity [mm/h]", color = "T [years]", caption = sprintf("Subdaily gauge: %s", gg)) +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.key.height = unit(fontsize - 1, "pt"),
        legend.title = element_text(size = fontsize - 2))

q.nprior <- lapply(durations, q.dgev, p.nexceed = p, param = par.nprior.gg)
q.wprior <- lapply(durations, q.dgev, p.nexceed = p, param = par.wprior.gg)
idf.nprior.aux <- matrix(unlist(q.nprior), length(p), length(durations), dimnames = list(rp, durations))
idf.wprior.aux <- matrix(unlist(q.wprior), length(p), length(durations), dimnames = list(rp, durations))
idf.nprior <- bind_cols(gauge_code = gg, idf.nprior.aux) %>% 
  mutate(return_period = rp) %>% 
  pivot_longer(cols = -c(gauge_code, return_period),
               names_to = "d",
               names_transform = list(d = as.numeric),
               values_to = "q") %>% print()
idf.wprior <- bind_cols(gauge_code = gg, idf.wprior.aux) %>% 
  mutate(return_period = rp) %>% 
  pivot_longer(cols = -c(gauge_code, return_period),
               names_to = "d",
               names_transform = list(d = as.numeric),
               values_to = "q") %>% print()

ggplot() +
  geom_line(aes(x = d, y = q, color = factor(return_period)), idf.nprior, linewidth = 0.5, linetype = "dashed") +
  geom_line(aes(x = d, y = q, color = factor(return_period)), idf.wprior, linewidth = 0.5, linetype = "solid") +
  # geom_point(aes(x = d, y = imax), data.gg[data.gg$d %in% durations,], alpha = 0.4, color = "grey10", size = 0.8) +
  scale_x_log10(breaks = durations) +
  scale_y_log10() +
  labs(x = "Duration [h]", y = "Intensity [mm/h]", color = "T [years]", caption = sprintf("Subdaily gauge: %s", gg)) +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.key.height = unit(fontsize - 1, "pt"),
        legend.title = element_text(size = fontsize - 2))

ggsave("figuras/idf_compare.png", width = 12, height = 8, units = "cm", dpi = 300)

data.gg %>% 
  filter(d %in% c(1, 2, 6, 12, 24)) %>% 
  ggplot(aes(x = wateryear, y = imax, group = factor(d))) +
  geom_line() +
  labs(caption = sprintf("Gauge: %s", unique(data.gg$gauge_code))) +
  theme(legend.position = "bottom")
  


# VALIDAÇÃO ---------------------------------------------------------------

library(IDF)
IDF.agg(list(as.data.frame(df.data[df.data$gauge_code == gg,])), ds = durations, names = c("datetime", "rain_mm"), na.accept = 0.2)
xdat <- data.gg %>% rename(xdat = imax, ds = d)
fit.ulrich <- IDF::gev.d.fit(xdat = xdat$xdat, ds = xdat$ds)
par.ulrich <- IDF::gev.d.params(fit.ulrich) %>% print()
par.eu <- fit$par %>% print()
