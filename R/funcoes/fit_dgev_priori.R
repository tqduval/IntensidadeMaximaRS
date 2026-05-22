#' @param data.gg a data.frame specific to a given gauge (gg) containing 'year', 'duration' and 'imax' columns
#' @param cols names(data.gg), following the same order specified above
#' @param method either Nelder-Mead or L-BFGS-B, the latter is tried first, if there is an issue with convergence, fallsback to NM
#' @param scale.inv 
#' @param maxit maximum number of iterations for 'optim()'

fit.dgev <- function(data.gg,
                     prior.info = NULL,
                     cols = c("gauge_code", "year", "duration", "imax"),
                     method = c("Nelder-Mead", "L-BFGS-B"),
                     scale.inv = c("wide", "strict"),
                     # offset = TRUE,
                     maxit = 1e6){
  
  
  #### CHECAGENS ----
  
  # Conferir argumentos
  if(!is.data.frame(data.gg)) stop("Argumento 'data.gg' deve ser 'tbl_df' e não ", class(data.gg))
  if(sum(is.element(cols, names(data.gg))) != length(cols)) stop("Colunas indicadas em 'cols' devem estar contidas em 'data.gg'\ncols: ", cols, "\nnames(data.gg): ", names(data.gg))
  
  method <- match.arg(method, c("Nelder-Mead", "L-BFGS-B"))
  scale.inv <- match.arg(scale.inv, c("wide", "strict"))
  
  #### FUNÇÕES ----
  
  # Função de log-verossimilhança do modelo d-GEV
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
    if(theta < 0) return(1e6)                  # impedir offset negativo
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
    if(any(z.id <= 0)) return(1e8)
    # return(1e10)
    
    # Log-beta
    if(!is.null(prior.info)){
      
      mu.kappa <- prior.info[1]
      var.kappa <- prior.info[2]^2
      
      a <- 0.5
      b <- mu.kappa+a
      p <- b^2*(1-b)/var.kappa-b
      q <- p*(1/b-1)
      
      l.prior <- (p - 1)*log(a + mu.kappa) + (q - 1)*log(a - mu.kappa)
      
    } else{
      
      l.prior <- 0
      
    }
    
    
    # Log-verossimilhança
    nll <- -(-N*D*log(alpha0) + H*N*sum(log(durations + theta)) + (1/kappa - 1)*sum(log(z.id)) - sum(z.id^(1/kappa)) + l.prior)
    
    return(nll) # retornando negativo p/ "minimização"
    
  }
  
  # Função de estimativas iniciais do modelo d-GEV
  start.dgev <- function(y.id, durations, scale.inv){
    
    # Estimar parâmetros da GEV p/ cada duração
    par.d <- apply(y.id, 2, function(imax){
      sample.lmom <- lmom::samlmu(imax)
      par <- lmom::pelgev(sample.lmom)
    })
    
    # Expoente 'scale' de invariância de escala
    scale <- switch(
      scale.inv,
      wide = {
        mom <- apply(y.id, 2, function(imax) mean(imax))
        scale <- -lm(log(mom) ~ log(durations))$coefficient[[2]]
      },
      strict = {
        scale.xi <- -lm(log(par.d["alpha",]) ~ log(durations))$coefficient[[2]]
        scale.alpha <- -lm(log(par.d["alpha",]) ~ log(durations))$coefficient[[2]]
        scale <- mean(scale.xi, scale.alpha)
      }
    )
    
    if(scale < 0) scale <- 0.5
    theta <- 0
    alpha <- par.d[2,1]                          # estimate for first duration
    xit <- par.d[1,1]/alpha                      # rescaled xi
    alpha0 <- alpha*(durations[1] + theta)^scale # rescaled alpha
    kappa <- mean(par.d[3,])                     # average kappa
    kappa <- ifelse(abs(kappa) >= 0.5, prior.info[1], kappa)
    # kappa <- min(0.5 - 1e-6, max(-0.5 + 1e-6, kappa))
    
    
    names <- c("xit", "alpha0", "kappa", "scale", "theta")
    
    return(setNames(c(xit, alpha0, kappa, scale, theta), names))
    
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
    
    # Derivadas parciais
    # Posição (xit)
    dl.dxit <- (1 - kappa)*sum(z.inv) - sum(z.pow)
    
    # Escala (alpha0)
    dl.dalpha0 <- -(N*D)/alpha0 + (1 - kappa)*sum(yA*z.inv/(alpha0^2)) - sum(z.pow*yA/(alpha0^2))
    
    # Forma (kappa)
    
    if(!is.null(prior.info)){
      
      mu.kappa <- prior.info[1]
      var.kappa <- prior.info[2]^2
      
      a <- 0.5
      b <- mu.kappa+a
      p <- b^2*(1-b)/var.kappa-b
      q <- p*(1/b-1)
      if(a - abs(kappa) <= 1e-8) return(1e6)
      
      dl.prior <- (p - 1)/(a + prior.info[1]) + (q - 1)/(a - prior.info[1])
      
    } else{
      
      dl.prior <- 0
      
    }
    
    term1 <- (1/kappa - 1)*B*z.inv + ln.z/(kappa^2)
    term2 <- (1/kappa)*((ln.z*z.pow.k)/kappa + B*z.pow)
    dl.dkappa <- -sum(term1) + sum(term2) + dl.prior
    
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
  
  #### CÁLCULOS ----
  
  # Configurar dados
  durations <- unique(data.gg[[cols[3]]])
  gg <- unique(data.gg[[cols[1]]])
  data.gg <- data.gg[, cols] # manter somente imax e d
  y.id <- data.gg %>% 
    pivot_wider(names_from = cols[3], # transformar durações em colunas
                values_from = cols[4],
                names_prefix = "d_") %>% 
    arrange(cols[2]) %>% 
    select(-cols[1:2]) %>% 
    na.omit() %>% as.matrix()
  
  # Definir limites para L-BFGS-B
  lower <- c(-Inf, 1e-6, ifelse(is.null(prior.info), -Inf, -0.5 + 1e-6), 1e-6, 0)
  upper <- c(Inf, Inf, ifelse(is.null(prior.info), -Inf, 0.5 - 1e-6), 1 - 1e-6, Inf)
  
  # Parâmetros iniciais (chute)
  start <- start.dgev(y.id, durations, scale.inv)
  
  # Ajustar d-GEV
  # Usando "L-BFGS-B"
  if(method == "L-BFGS-B"){
    fit <- tryCatch(
      optim(par = start, fn = nll.dgev, gr = grad.ll.dgev,
            method = "L-BFGS-B", lower = lower, upper = upper,
            hessian = TRUE,
            durations = durations, y.id = y.id,
            control = list(maxit = maxit)),
      error = function(e){
        warning("L-BFGS-B failed (", e$message, "). Otimizando com Nelder-Mead.")
        return(NULL)
      }
    )
    
    if(!is.null(fit) && fit$convergence == 0) return(fit)
  }
  
  # Usando Nelder-Mead
  fit <- optim(par = start, fn = nll.dgev, gr = grad.ll.dgev,
               method = "Nelder-Mead",
               durations = durations, y.id = y.id,
               control = list(maxit = maxit))
  
  if(fit$value == 1e6) warning("Não houve otimização na estação ", gg, "! Retornando parâmetros iniciais:\n", paste(round(start, 4), collapse = " "))
  
  return(fit)
  
}
