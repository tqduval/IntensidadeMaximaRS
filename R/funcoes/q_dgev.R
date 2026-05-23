#' @title Estimar quantis de intensidade de chuva
#' @details
#' Esta função constrói as curvas de intensidade-duração-frequência a partir dos quantis estimados para diferentes durações e frequências
#' @param fit.obj o resultado da função `fit.dgev()`, contendo uma lista com um resultados padrão da função `optim()` para cada estação
#' @param durations um vetor com durações [h] para as quais serão calculadas as intensidades
#' @param return.period um vetor com períodos de retorno (inverso da probabildiade de excedência) [anos]
#' @returns um `data.frame` com intervalos de confiança estimados para quantil de cada duração, tempo de retorno e estação.
#' 
q.dgev <- function(fit.obj, durations, return.period){
  
  #### CHECAGENS ----
  
  if(!inherits(fit.obj, "list")) stop("fit.obj' deve conter uma lista de componentes do 'optim()' com ao menos: par, value e hessian")
  if(length(durations) <= 3) stop("Informe ao menos 3 durações em horas.")
  
  gauges <- names(fit.obj)
  
  #### FUNÇÕES ----
  
  # Calcular quantis
  get.q <- function(param, duration, p.nexceed){
    
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
  
  # Calcular quantis p/ todas as estações
  q.all <- function(fit, durations, p.nexceed){
    
    par <- fit.obj$par; p <- p.nexceed; rp <- round(1/(1 - p), 0)
    ls.q <- lapply(durations, get.q, p.nexceed = p, param = par)
    idf.aux <- matrix(unlist(ls.q), length(p), length(durations), dimnames = list(rp , durations))
    
    return(idf.aux)
    
  }
  
  #### CÁLCULOS ----
  
  # Converter return.period em probabilidade de não excedência
  p <- 1 - 1/return.period
  idf <- q.all(fit.obj, durations, p)
  
  return(idf)
  
}