
# INTERNAL FUNCTIONS -----------------------------------------------------

# Fit GEV distribution with maximum likelihood or
# generalized maximum likelihood --- 'fgev()'
source("R/funcoes/fgev.R")

#' @title (n - 1) negative log-likelihood function
#' @details Takes the original negative log-likelihood and adapts it
#' to optimize n - 1 parameters, keeping one of them fixed
prof.nll.gev <- function(par.free, which.par, par.value, y, prior.info){
  i <- which.par                        # index in 'par' vector of which parameter to build CIs for
  full <- numeric(length(par.free) + 1) # build 'new' parameters vector
  full[i] <- par.value                  # fixed value for current parameter
  full[-i] <- par.free                  # remaining parameters to be optimized
  nll.gev(full, y, prior.info)          # standard GEV negative log-likelihood
}

#' @title (n - 1) log-likelihood gradient function
#' @details Takes the original negative log-likelihood gradient function and adapts it
#' to optimize n - 1 parameters, keeping one of them fixed
prof.grad.gev <- function(par.free, which.par, par.value, y, prior.info){
  i <- which.par
  full <- numeric(length(par.free) + 1)
  full[i] <- par.value
  full[-i] <- par.free
  gr <- grad.ll.gev(full, y, prior.info)
  gr[-i] # remove par[i] from output vector
}

#' @title Find the 'roots' of the profile log-likelihood function
#' @details This function optimizes nll.prof.gev() with respecto to a desired parameter
#' and estimates its Confidence Intervals to a given confidence level 'conf'.
#' It may happen that this function ends up finding 'better' optimal parameters
#' than the original fgev() function
roots.proflik <- function(step, conf, par, i, ll){
  
  # val.i <- par[[i]] + step      # iniciar contador
  val.i <- par[[i]]             # iniciar contador
  npar <- length(par)           # number of parameters
  last.par <- par               # salvar último conjunto de parâmetros
  lim <- ll - qchisq(conf, 1)/2 # limite de definição dos ICs
  result <- list()              # iniciar lista com resultados

  # Limites inferior e superior
  lower <- c(-Inf, 1e-6, if (is.null(prior.info)) -Inf else -0.5 + 1e-6)
  upper <- c(Inf, Inf, if (is.null(prior.info)) Inf else 0.5 - 1e-6)
  
  # Start of 'while search'
  # while(ll > lim){
  while(ll > lim){
    
    fit.prof <- tryCatch(
      optim(
        last.par[-i], fn = prof.nll.gev, gr = prof.grad.gev, method = "L-BFGS-B", lower = lower[-i], upper = upper[-i],
        which.par = i, par.value = val.i, prior.info = prior.info, y = y,
        control = list(maxit = maxit, factr = 1e9, pgtol = 1e-8)
      ),
      error = function(e) NULL
    )
    
    if(!is.null(fit.prof) && fit.prof$convergence == 0){
      ll <- -fit.prof$value
    } else{
      fit.prof <- tryCatch(
        optim(
          last.par[-i], fn = prof.nll.gev, gr = NULL, method = "Nelder-Mead",
          which.par = i, par.value = val.i, prior.info = prior.info, y = y,
          control = list(maxit = maxit)
        ),
        error = function(e) NULL
      )
      if(!is.null(fit.prof) && fit.prof$convergence == 0){
        ll <- -fit.prof$value
        if(ll < -1e6){
          ll <- NA_real_
        }
      } else{
        ll <- NA_real_
      }
    }
    
    full <- numeric(npar)
    if(!is.null(fit.prof) && !is.na(ll)){
      full[-i] <- fit.prof$par # se não falhar, atualizar parâmetros
    } else{
      full[-i] <- last.par[-i] # se falhar, usar rodada anterior
    }
    
    full[i] <- val.i      # 
    val.i <- val.i + step # atualizar contador
    last.par <- full      # atualizar conjunto de parâmetros
    current.iter <- length(result) + 1
    result[[current.iter]] <- c(full, ll) # incluir resultados na lista
    
    # message(sprintf("Intervalo: %-8s | Parâmetro: %s = %0.3f | Loglik: %0.3f\n", 
    #                 ifelse(step < 0, "inferior", "superior"), names(par)[i], val.i, ll))
    
  } # end of 'while search'
  
  prof.par <- do.call(rbind, result)
  colnames(prof.par) <- c(names(par), "loglik")
  
  return(prof.par)
  
}

# TEST FUNCTION ----------------------------------------------------------

rm(list = ls()); invisible(gc())

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
gg <- "A801"
# gg <- "A812"
data.gg <- data[data$gauge_code == gg, ]
cols <- names(data)
maxit <- 1e8

y.id <- data.gg %>%
  pivot_wider(names_from = d, values_from = imax, names_prefix = "d_") %>%
  arrange(wateryear) %>%
  select(-c(gauge_code, wateryear)) %>%
  as.matrix() %>%
  na.omit()

d <- 20
y <- y.id[,d]

# Original optimal parameters and maximum
fit <- fgev(y, prior.info)
par <- fit$par
npar <- length(par)
max.ll <- -fit$value

# Define confidence level and a search step
lltol <- 0.1   # stablish a tolerance for differences in fgev-ll and proflik-ll
i <- 1         # which parameter to derive CIs for
conf <- 0.95   # confidence level
step <- 0.0025 # search step

# Find values to the left and to the right of the original 'par[i]' estimate
left <- roots.proflik(-step, conf, par, i, max.ll)
right <- roots.proflik(step, conf, par, i, max.ll)
proflik <- rbind(left, right)
ci.lower <- min(proflik[[i]]) # lower CI limit
ci.upper <- max(proflik[[i]]) # upper CI limit

# Update optimal parameters if proflik finds a higher optimal point
ll.max.proflik <- max(proflik[,"loglik"])
ll.dif <- abs((ll.max.proflik - ll)/ll) <= lltol
if(isFALSE(ll.dif)){
  message <- 1 # updated parameters
  out.par <- unlist(proflik[which.max(proflik[,"loglik"]), seq_len(npar)])
} else{
  message <- 0 # original fgev() parameters
  out.par <- par
}

# Manter ajuste com série original, extração de 'par' e 'll' antes do loop
profile <- lapply(X = seq_along(par), FUN = function(i){

  # Find values to the left and to the right of the original 'par[i]' estimate
  left <- roots.proflik(-step, conf, par, i, max.ll)
  right <- roots.proflik(step, conf, par, i, max.ll)
  proflik <- rbind(left, right)
  ci.lower <- min(proflik[[i]]) # lower CI limit
  ci.upper <- max(proflik[[i]]) # upper CI limit

  # Update optimal parameters if proflik finds a higher optimal point
  ll.max.proflik <- max(proflik[,"loglik"])
  ll.dif <- abs((ll.max.proflik - ll)/ll) <= lltol
  if(isFALSE(ll.dif)){
    message <- 1 # updated parameters
    out.par <- unlist(proflik[which.max(proflik[,"loglik"]), seq_len(npar)])
  } else{
    message <- 0 # original fgev() parameters
    out.par <- par
  }

  parname <- names(par)[i]
  out <- tibble(
    parname = parname, 
    par = unname(par[i]), 
    ci_lower = ci.lower,
    ci_upper = ci.upper, 
    max_ll = max(ll, ll.max.proflik, na.rm = TRUE), 
    message = message, 
    conf = conf, 
    step = step
  )

}) # fim 'profile'

do.call(rbind, profile)