
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
#' @details This function optimizes `nll.prof.gev()` with respecto to a desired parameter
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
        control = list(maxit = maxit, factr = 1e9, pgtol = 1e-8, parscale = abs(last.par[-i]))
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
          control = list(maxit = maxit, parscale = abs(last.par[-i]))
        ),
        error = function(e) NULL
      )
      if(!is.null(fit.prof) && fit.prof$convergence %in% c(0, 10)){
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

#' @title Estimate profile likelihood confidence intervals for a parameter
#' @details This functions calls `roots.proflik()` and estimates confidence intervals
#' for a given parameter. It's called inside an `lapply()` loop and returns a list with:
#' 1) a row tbl_df with original parameters from the log-likelihood maximization and 
#' estimated confidence intervals;
#' 2) a matrix with all evaluated values of the parameters and respective profile values; and
#' 3) a vector with optimal profile log-likelihood parameters and its maximum.
ciprof.par <- function(par, max.ll, lltol, step, conf){

  # Loop through all GEV parameters
  ls.profile <- lapply(X = seq_along(par), FUN = function(i){

    # Find values to the left and to the right of the original 'par[i]' estimate
    left <- roots.proflik(-step, conf, par, i, max.ll)
    right <- roots.proflik(step, conf, par, i, max.ll)
    proflik <- rbind(left, right)
    ci.lower <- min(proflik[,i]) # lower CI limit
    ci.upper <- max(proflik[,i]) # upper CI limit

    # Update optimal parameters if proflik finds a higher optimal point
    ll.max.proflik <- max(proflik[,"loglik"])
    ll.dif <- abs((ll.max.proflik - max.ll)/max.ll) <= lltol
    if(isFALSE(ll.dif)){
      message <- 1 # updated parameters
      out.par <- unlist(proflik[which.max(proflik[,"loglik"]), seq_len(npar)])
    } else{
      message <- 0 # original fgev() parameters
      out.par <- par
    }

    parname <- names(par)[i]
    profile <- proflik[,c(i, ncol(proflik))]
    colnames(profile) <- c("par", "loglik")
    out <- list(
      ci = data.frame(
        gauge_code = gg,
        d = d,
        parname = parname,
        par = unname(par[i]),
        ci_lower = ci.lower,
        ci_upper = ci.upper,
        max = max.ll,
        prior = deparse(prior.info),
        message = message, # if message == 1, 'prof.pars' should be checked
        algo = algo,
        conf = conf,
        step = step
      ),
      prof.par = unlist(proflik[which.max(proflik[,"loglik"]), seq_len(npar + 1)]),
      profile = as.data.frame(cbind(parname = parname, profile))
    )

    return(out)

  }) # fim 'ls.profile'

  names(ls.profile) <- names(par) # name parameters

  results <- bind_rows(lapply(ls.profile, `[[`, "ci"))
  prof.par <- bind_rows(lapply(ls.profile, `[[`, "prof.par"))
  profile <- bind_rows(lapply(ls.profile, `[[`, "profile"))

  if(any(results$message == 1)){
    new.par <- unlist(prof.par[which.max(prof.par$loglik), seq_len(npar)])
    new.max <- max(prof.par$loglik)
    results$par <- new.par
    results$max <- new.max
  }

  return(list(ci = results, profile = profile))

}


# MAIN FUNCTION ----------------------------------------------------------

fitprof.gev <- function(
  data, prior.info = NULL,
  durations, return.period,
  cols = c("gauge_code", "year", "duration", "imax"),
  conf = 0.95, lltol = 0.1, step = 0.01
){}

# TEST FUNCTION ----------------------------------------------------------

rm(list = ls()); invisible(gc())
pacman::p_load(pacman, tidyverse, arrow)

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

# Define confidence level and a search step
conf <- 0.95
lltol <- 0.1   # stablish a tolerance for differences between fgev-ll and proflik-ll
step <- 0.0025 # search step


## PARÂMETROS ------------------------------------------------------------

# Function
# Original optimal parameters and maximum
fit <- fgev(y, prior.info)       # this is the original estimate with a maximum log-likelihood algorithm
gradient <- fit[["counts"]][[2]] # check if gradient is NA or number
if(is.na(gradient)){
  algo <- "Nelder-Mead"
} else{
  algo <- "L-BFGS-B"
}
par <- fit$par
npar <- length(par)
convergence <- fit$convergence
max.ll <- -fit$value

# Test run
ls.ci.par <- ciprof.par(par, max.ll, lltol, step, conf)
opt.par <- ls.ci.par[["ci"]][["par"]]
opt.ll <- max(ls.ci.par[["ci"]][["max"]])

# Run again with optimal 'new.par' parameters
rerun.ci.par <- ciprof.par(setNames(opt.par, names(par)), opt.ll, lltol, step, conf)
opt.par <- setNames(rerun.ci.par[["ci"]][["par"]], names(par))
opt.ll <- max(rerun.ci.par[["ci"]][["max"]])

## QUANTIS ---------------------------------------------------------------

# Define durations and return periods
# durations <- c(1, 6, 12, 24) # durations is already set as 'd'
return.period <- 10
p <- 1 - 1/return.period

# Estimate quantiles with optimal parameters
q <- qgev(param = opt.par, p.nexceed = p)

# Reparametrize quantile expression in terms of xi
alpha <- opt.par[[2]]; kappa <- opt.par[[3]]
xi.q <- q - alpha/kappa*(1 - (-log(p))^kappa) # estimate a new 'xi'
par.q <- opt.par # create a par vector to be optimized with respect to alpha and kappa
par.q[1] <- xi.q # replace existing 'xi' estimate with new 'xi.q'

# Now the new parameter to be fixed is xi.q, which is in reality 'q'
prof.nll.gev(par.free = par.q[-1], which.par = 1, par.value = xi.q, y, prior.info)
ll <- opt.ll
roots.proflik(step, conf, par = par.q, i = 1, ll = opt.ll)
# ALTERAR: isso não funcionaria pq que deve ser alterado é q e não xi.q

# rp must have only one value

roots.proflik.q <- function(step, conf, par, i, ll, rp){

  # Reparametrize quantile expression in terms of xi
  p <- 1 - 1/rp
  q <- qgev(par, p)
  alpha <- par[[2]]
  kappa <- par[[3]]
  xi.q <- q - alpha/kappa*(1 - (-log(p))^kappa) # estimate a new 'xi'
  par.q <- opt.par # create a par vector to be optimized with respect to alpha and kappa
  par.q[1] <- xi.q # replace existing 'xi' estimate with new 'xi.q'
  i <- 1

  # val.i <- par[[i]] + step      # iniciar contador
  val.i <- par[[1]]             # iniciar contador
  npar <- length(par.q)         # number of parameters
  last.par <- par.q             # salvar último conjunto de parâmetros
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
        control = list(maxit = maxit, factr = 1e9, pgtol = 1e-8, parscale = abs(last.par[-i]))
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
          control = list(maxit = maxit, parscale = abs(last.par[-i]))
        ),
        error = function(e) NULL
      )
      if(!is.null(fit.prof) && fit.prof$convergence %in% c(0, 10)){
        ll <- -fit.prof$value
        if(ll < -1e6){
          ll <- NA_real_
        }
      } else{
        ll <- NA_real_
      }
    }
    
    # FAZER: aqui deve ser feita a atualização de 'q'
    full <- numeric(npar + 1)
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

ciprof.q <- function(param, return.period){

  p <- 1 - 1/return.period
  q <- qgev(param, p)
  

}