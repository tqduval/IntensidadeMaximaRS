# source("R/funcoes/fit_dgev.R")
# source("R/funcoes/q_dgev.R")
# 
# 
# gg <- "A801"
# data.gg <- data[data$gauge_code == gg,]
# durations <- unique(data.gg$d)
# prior.info <- c(0.0178, 0.0627)
# method <-  "L-BFGS-B"
# conf <- 0.95

fitdelta.dgev <- function(data, prior.info = NULL,
                          durations, return.period,
                          cols = c("gauge_code", "year", "duration", "imax"),
                          method = c("Nelder-Mead", "L-BFGS-B"),
                          scale.inv = c("wide", "strict"),
                          conf = 0.95,
                          maxit = 1e6){
  
  #### CHECAGENS ----
  
  if(!inherits(data, "list")) stop("'data' deve ser uma lista com data.frames e não: ", class(data), ".")
  if(!length(prior.info) %in% c(0, 2)) stop("'prior.info' dever ser 'NULL' ou um vetor de comprimento 2.")
  
  #### FUNÇÕES ----
  
  # Ajustar d-GEV e estimar quantis
  source("R/funcoes/fit_dgev.R")
  
  # Função vetorizada p/ quantil
  q.dgev.vec <- function(param, d, p){
    
    xit <- param[1]; alpha0 <- param[2]; kappa <- param[3]; H <- param[4]; theta <- param[5]
    A <- alpha0/(d + theta)^H
    B <- xit + (1/kappa)*(1 - (-log(p))^kappa)
    
    yp <- A*B
    
    return(unname(yp))
    
  }
  
  # Função gradiente p/ quantil
  grad.q.dgev <- function(param, d, p){
    
    xit <- param[1]; alpha0 <- param[2]; kappa <- param[3]; H <- param[4]; theta <- param[5]
    logp <- -log(p)
    A <- alpha0/(d + theta)^H
    B <- xit + (1/kappa)*(1 - logp^kappa)
    yp <- A*B
    
    dy.dxit <- A
    dy.dalpha0 <- B/(d + theta)^H
    dy.dkappa <- -A/kappa^2*(1 + logp^kappa*(kappa*log(logp) - 1))
    dy.dH <- -yp*log(d + theta)
    dy.dtheta <- -yp*H/(d + theta)
    
    rp <- 1/(1 - p)
    grad.names <- paste0(names(param), "_d", rep(d, each = length(param)), "_rp", rp)
    return(setNames(c(dy.dxit, dy.dalpha0, dy.dkappa, dy.dH, dy.dtheta), grad.names))
    
  }
  
  #### CÁLCULOS ----
  
  # Rodar para cada estação
  message(paste0(rep("_", getOption("width")*0.5), collapse = ""), "\n")
  message("Ajustando d-GEV, estimando quantis.")
  message(">> Número de estações: ", length(data))
  message(">> Otimizador: ", method)
  message(">> Durações [h]: ", paste(durations, collapse = ", "))
  message(">> Tempos de retorno [anos]: ", paste(return.period, collapse = ", "))
  message(">> Nível de confiança para intervalo percentil com método delta: ", paste0(conf*100, "%"))
  
  probs <- 1 - 1/return.period
  
  ls.res <- pbapply::pblapply(data, function(data.gg){
    
    gg <- unique(data.gg$gauge_code)
    fit <- fit.dgev(data.gg, prior.info, names(data.gg), method, "wide")
    par <- fit$par
    hes <- fit$hessian
    
    # Intervalos de confiança parâmetros
    # Extrair elementos das diagonais
    cov.par <- try(solve(hes), silent = TRUE)
    if(inherits(cov.par, "try-error")) return(matrix(NA, length(par), length(par)))
    var.par <- diag(cov.par)
    sd.par <- sqrt(var.par)
    
    ci.par <- cbind(ci_lower = par - qnorm(1 - (1 - conf)/2)*sd.par, 
                    ci_upper = par + qnorm(1 - (1 - conf)/2)*sd.par)
    
    # Intervalos de confiança quantis
    grid <- expand.grid(d = durations, p = probs)
    ls.ci.q <- lapply(1:nrow(grid), function(i){
      
      d <- grid[i, "d"]
      p <- grid[i, "p"]
      
      q <- q.dgev.vec(par, d, p)
      grad <- grad.q.dgev(par, d, p)
      var.q <- t(grad) %*% cov.par %*% grad
      sd.q <- sqrt(var.q)
      ci.q.l <- q - qnorm(1 - (1 - conf)/2)*sd.q
      ci.q.u <- q + qnorm(1 - (1 - conf)/2)*sd.q
      
      out <- c(d = d, return_period = 1/(1 - p), q = q, sd = sd.q, ci_lower = ci.q.l, ci_upper = ci.q.u)
      
    }) # fim 'ls.ci.q'
    
    parameters <- bind_cols(gauge_code = gg, parname = names(par), par = par, ci.par)
    quantiles <- bind_cols(gauge_code = gg, bind_rows(ls.ci.q)) %>% arrange(d, return_period)
    
    out <- list(par = parameters,
                idf = quantiles)
    
  }) # fim 'ls.res'
  
  message(paste0(rep("_", getOption("width")*0.5), collapse = ""), "\n")
  
  return(ls.res)
  
}

