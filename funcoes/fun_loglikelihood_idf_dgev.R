# Construção da função de log-verossimilhança do modelo de equação de IDF baseado
# na distribuição dGEV (duration-dependent Generalized Extreme Value distribution)

# Provavelmente, vai precisar fazer uma organização dos dados igual na função de
# invariância de escala, organizando em cada coluna a AMS de uma duração diferente
# pq a função vai precisar iterar em um loop das durações e das intensidades
# então, entra ana primeira coluna, extrai o primeiro 'd', e os 'y_i' pro somatório

rm(list = ls()); invisible(gc())


pacman::p_load(pacman, tidyverse, arrow)

df.imax <- arrow::read_parquet("C:/Users/tomas/OneDrive/1_Academico/Mestrado/Projeto IDF/Proj_ANA_UnB_comp1/base/gerados/df_imax.parquet")
data <- df.imax %>% 
  filter(if_else(year == 2024, TRUE, na_prct <= 0.2),
         between(d, 1, 24)) %>% 
  group_by(gauge_code) %>% 
  filter(n_distinct(d) >= 3,
         n_distinct(year) >= 8) %>% 
  ungroup() %>% 
  select(gauge_code, year, d, imax)

gauges <- unique(data$gauge_code)
gauge <- gauges[[2]]

data.gauge <- data %>% 
  filter(gauge_code == gauge)

y.id <- data.gauge %>% 
  pivot_wider(names_from = d,
              values_from = imax,
              names_prefix = "d_") %>% 
  arrange(year) %>% 
  select(-c(gauge_code, year)) %>% 
  as.matrix()

durations <- as.numeric(stringr::str_remove_all(colnames(y.id), "d_"))

# Parâmetros iniciais
xi <- param[1]     # posição: xi~
alpha <- param[2]  # escala: alpha0
kappa <- param[3]  # forma
scale <- param[4]  # invariância de escala (expoente de escala - H ou eta)
offset <- param[5] # invariância de escala (parâmetro de posição - theta)

N <- nrow(y.id)
D <- ncol(y.id)

# Limites teóricos
if(any(durations + offset <= 0)) return(1e10)
if(alpha <= 0) return(1e10)
if(scale <= 0 | scale > 1) return(1e10)

# Calcular termos z.id = 1 - kappa(y.id*(d + offset)^scale)/alpha - xi~
# Função sweep() realiza uma operação nas linhas ou colunas de euma matriz
aux1 <- (durations + offset)^scale # termo de alpha(d)
aux2 <- sweep(x = y.id,            # matriz c/ imax
             MARGIN  = 2,          # operação por colunas
             STATS = aux1,         # argumentos (nesse caso oq vai ser multiplicado por y.id)
             FUN = "*")            # multiplicação

z.id <- 1 - kappa*(aux2 - xi)

# Conferir limites teóricos
if(any(z.id <= 0)) return(1e10)

# Log-verossimilhança
l <- -N*D*log(alpha) + scale*N*sum(log(durations + offset)) + (1/kappa - 1)*sum(log(z.id)) - sum(z.id^(1/kappa))

return(-l) # retornando negativo p/ "minimização"


# FUNÇÃO ------------------------------------------------------------------

fun.l.dgev <- function(param,     # vetor c/ parâmetros 1. posicao 2. escala 3. forma 4. eta 5. theta
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
  alpha.d <- (durations + offset)^scale # termo de alpha(d)
  aux <- sweep(x = y.id,                # matriz c/ imax
               MARGIN  = 2,             # operação por colunas
               STATS = alpha.d,         # argumentos (nesse caso oq vai ser multiplicado por y.id)
               FUN = "*")               # multiplicação
  
  z.id <- 1 - kappa*(aux - xit) # variável reduzida vetorizada
  
  # Conferir limites teóricos
  if(any(z.id <= 0)) return(1e6)
  # return(1e10)
  
  # Log-verossimilhança
  l <- -N*D*log(alpha0) + scale*N*sum(log(durations + offset)) + (1/kappa - 1)*sum(log(z.id)) - sum(z.id^(1/kappa))
  
  return(-l) # retornando negativo p/ "minimização"
  
}

# Função gradiente log-verossimilhança


# Estimar parâmetros da GEV p/ cada duração
param.d <- apply(X = y.id, MARGIN = 2, FUN = function(imax){
  sample.lmom <- lmom::samlmu(imax)
  gev.par <- lmom::pelgev(sample.lmom)
  
})


# Parâmetros iniciais
start <- c(xit = param.d[1,1],
           alpha0 = param.d[2,1],
           kappa = param.d[3,1],
           scale = 0.6,
           offset = 0)


# Otimizar
optim(par = start,
      fn = fun.l.dgev,
      hessian = TRUE,
      y.id = y.id,
      durations = durations,
      control = list(maxit = 1e6,
                     trace = 1))

    
# Forma alternativa da função: usando for loop explícito
N <- nrow(y.id)
D <- ncol(y.id)

for(i in 1:N){
  
  for(d in 1:D){
    
    y <- y.id[i,d]
    z <- 
    
  }
  
}