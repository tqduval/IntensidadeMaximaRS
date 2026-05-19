# Construção da função de log-verossimilhança do modelo de equação de IDF baseado
# na distribuição dGEV (duration-dependent Generalized Extreme Value distribution)

# Provavelmente, vai precisar fazer uma organização dos dados igual na função de
# invariância de escala, organizando em cada coluna a AMS de uma duração diferente
# pq a função vai precisar iterar em um loop das durações e das intensidades
# então, entra na primeira coluna, extrai o primeiro 'd', e os 'y_i' pro somatório

rm(list = ls()); invisible(gc())


pacman::p_load(pacman, tidyverse, arrow, numDeriv, evd)

df.imax <- arrow::read_parquet("base/gerados/df_imax.pqt")

data <- df.imax %>% 
  filter(na_prct_yr <= 0.2, na_prct_mon <= 0.1, between(d, 1, 24)) %>% 
  group_by(gauge_code) %>% 
  filter(n_distinct(d) >= 3, n_distinct(wateryear) >= 8) %>% 
  ungroup() %>% 
  select(gauge_code, wateryear, d, imax)

gauges <- unique(data$gauge_code)
gauge <- "A802"

data.gauge <- data %>% 
  filter(gauge_code == gauge)

y.id <- data.gauge %>% 
  pivot_wider(names_from = d,
              values_from = imax,
              names_prefix = "d_") %>% 
  arrange(wateryear) %>% 
  select(-c(gauge_code, wateryear)) %>% 
  as.matrix() %>% 
  na.omit()

durations <- as.numeric(stringr::str_remove_all(colnames(y.id), "d_"))


# FUNÇÃO ------------------------------------------------------------------

nll.dgev <- function(param,     # vetor c/ parâmetros 1. posicao 2. escala 3. forma 4. eta 5. theta
                     durations, # vetor c/ durações (nomes de y.id)
                     y.id       # matrix c/ série de máximos anuais observada p/ cada duração em colunas
                     ){
  
  # Parâmetros iniciais
  xit <- param[1]     # posição: xi~
  alpha0 <- param[2]  # escala: alpha0
  kappa <- param[3]   # forma
  scale <- param[4]   # invariância de escala (expoente de escala - H ou eta)
  offset <- param[5]  # invariância de escala (parâmetro de posição - theta)
  
  N <- nrow(y.id)
  D <- ncol(y.id)
  
  # Limites teóricos
  if(any(durations + offset <= 0)) return(1e6) # impedir offset negativo
  if(alpha0 <= 0) return(1e6)                  # impedir escala GEV nulo
  if(scale <= 0 | scale > 1) return(1e6)       # restrição invariância de escala
  
  # Calcular termos z.id = 1 - kappa(y.id*(d + offset)^scale)/alpha - xi~
  # Função sweep() realiza uma operação nas linhas ou colunas de euma matriz
  A <- (durations + offset)^scale # termo de alpha(d)
  aux <- sweep(x = y.id/alpha0,   # matriz c/ imax
               MARGIN  = 2,       # operação por colunas
               STATS = A,         # argumentos (nesse caso oq vai ser multiplicado por y.id)
               FUN = "*")         # multiplicação
  
  z.id <- 1 - kappa*(aux - xit) # variável reduzida vetorizada
  
  # Conferir limites teóricos
  if(any(z.id <= 0)) return(1e6)
  # return(1e10)
  
  # Log-verossimilhança
  nll <- -(-N*D*log(alpha0) + scale*N*sum(log(durations + offset)) + (1/kappa - 1)*sum(log(z.id)) - sum(z.id^(1/kappa)))
  
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


# Estimar parâmetros da GEV p/ cada duração
param.d <- apply(X = y.id, MARGIN = 2, FUN = function(imax){
  sample.lmom <- lmom::samlmu(imax)
  gev.par <- lmom::pelgev(sample.lmom)
  
})


# Parâmetros iniciais
d1 <- durations[1]
xi1 <- param.d[1.1]
alpha1 <- param.d[2,1]
H0 <- 0.6
theta0 <- 0

xit0 <- xi1/alpha1
alpha0 <- alpha1*(d1 + theta0)^H0

# Estimar 

start <- c(xi1 = xit0,
           alpha0 = alpha0,
           kappa = param.d[3,1],
           scale = H0,
           offset = theta0)


# Otimizar
optim(par = start,
      fn = nll.dgev,
      hessian = TRUE,
      y.id = y.id,
      durations = durations,
      control = list(maxit = 1e6))

optim(par = start,
      fn = nll.dgev,
      hessian = TRUE,
      gr = grad.ll.dgev,
      y.id = y.id,
      durations = durations)

param <- start
ds <- as.numeric(gsub("d+_","", colnames(y.id)))

