
rm(list = ls()); invisible(gc())

# INTERNAL FUNCTIONS -----------------------------------------------------

# Fit GEV distribution with maximum likelihood or
# generalized maximum likelihood --- 'fgev()'
source("R/funcoes/fgev.R")

# @title (n - 1) negative log-likelihood function
# @details Takes the original negative log-likelihood and adapts it
# to optimize n - 1 parameters, keeping one of them fixed
prof.nll.gev <- function(par.free, which.par, par.value, y, prior.info){
  # i <- which.par                        # index in 'par' vector of which parameter to build CIs for
  full <- numeric(length(par.free) + 1) # build 'new' parameters vector
  full[which.par] <- par.value          # fixed value for current parameter
  full[-which.par] <- par.free          # remaining parameters to be optimized

  # if(full[1] < 0) return(1e9) # causam instabilidade nos intervalos
  # if(full[2] < 1e-3) return(1e9)
  # if(abs(kappa) >=)
  
  nll.gev(full, y, prior.info)          # standard GEV negative log-likelihood

}

# @title (n - 1) log-likelihood gradient function
# @details Takes the original negative log-likelihood gradient function and adapts it
# to optimize n - 1 parameters, keeping one of them fixed
prof.grad.gev <- function(par.free, which.par, par.value, y, prior.info){
  full <- numeric(length(par.free) + 1)
  full[which.par] <- par.value
  full[-which.par] <- par.free
  gr <- grad.ll.gev(full, y, prior.info)
  gr[-which.par] # remove par[i] from output vector
}

# @title Find the 'roots' of the profile log-likelihood function
# @details This function optimizes `nll.prof.gev()` with respecto to a desired parameter
# and estimates its Confidence Intervals to a given confidence level 'conf'.
# It may happen that this function ends up finding 'better' optimal parameters
# than the original fgev() function
roots.proflik <- function(step, conf, par, i, ll){
  
  val.i <- par[[i]]
  npar <- length(par)
  last.par <- par
  lim <- ll - qchisq(conf, 1)/2
  
  result <- list()
  
  lower <- c(-Inf, 1e-6, if (is.null(prior.info)) -Inf else -0.5 + 1e-6)
  upper <- c(Inf, Inf, if (is.null(prior.info)) Inf else 0.5 - 1e-6)
  
  while(ll > lim && !is.na(ll)){
    
    last.ll <- ll # keep last value of 'll'
    
    # Check 'last.par'
    if(i == 2){
      if(val.i <= lower[i]){
        ll <- NA_real_
        next
      }
    }
    
    fit.prof <- tryCatch(
      optim(
        par = last.par[-i],
        fn = function(internal.par){
          prof.nll.gev(
            par.free = internal.par,
            which.par = i,
            par.value = val.i,
            y = y,
            prior.info = prior.info
          )
        },
        gr = NULL,
        method = "L-BFGS-B",
        lower = lower[-i],
        upper = upper[-i],
        control = list(
          maxit = maxit,
          factr = 1e9,
          pgtol = 1e-8,
          parscale = abs(last.par[-i])
        )
      ),
      error = function(e) NULL
    )
    
    if(!is.null(fit.prof) && fit.prof$convergence == 0){
      
      ll <- -fit.prof$value
      
    } else{
      
      fit.prof <- tryCatch( # try again with Nelder-Mead
        optim(
          par = last.par[-i],
          fn = function(internal.par){
            prof.nll.gev(
              par.free = internal.par,
              which.par = i,
              par.value = val.i,
              y = y,
              prior.info = prior.info
            )
          },
          method = "Nelder-Mead",
          control = list(
            maxit = maxit,
            parscale = abs(last.par[-i])
          )
        ),
        error = function(e) NULL
      )
      
      if(!is.null(fit.prof) && fit.prof$convergence %in% c(0, 10)){
        ll <- -fit.prof$value
        if(ll < -1e4) ll <- NA_real_ # optim might reach convergence == 0, but converge to a break wall
      } else{
        ll <- NA_real_
      }
      
    }
    
    full <- numeric(npar)
    if(!is.null(fit.prof) && !is.na(ll)){
      full[-i] <- fit.prof$par
    } else{
      full[-i] <- last.par[-i]
    }
    
    full[i] <- val.i
    
    # Check GEV limits
    gev.lim <- full[1] + full[2]/full[3]  # xi + alpha/kappa
    
    if(full[3] < 0){
      if(any(y < gev.lim)) ll <- NA_real_ # domain violation for Fréchet law
    }
    if(full[3] > 0){
      if(any(y > gev.lim)) ll <- NA_real_ # domain violation for reverse Weibull law
    }
    
    val.i <- val.i + step
    slope <- (ll - last.ll)/step    # how is the algorithm jumping around
    last.par <- full
    current.iter <- length(result) + 1
    result[[current.iter]] <- c(full, ll, slope)
    
    message(
      sprintf(
        "ci %-5s | %s: %0.3f | ll: %0.3f | slope: %0.3f", 
        ifelse(step < 0, "lower", "upper"), 
        names(par)[i], val.i, ll, slope)
    )
    
  } # end of 'while' search
  
  prof.par <- do.call(rbind, result)
  colnames(prof.par) <- c(names(par), "loglik", "slope")
  
  return(prof.par)
  
}

# @title Estimate profile likelihood confidence intervals for a parameter
# @details This functions calls `roots.proflik()` and estimates confidence intervals
# for a given parameter. It's called inside an `lapply()` loop and returns a list with:
# 1) a row tbl_df with original parameters from the log-likelihood maximization and 
# estimated confidence intervals;
# 2) a matrix with all evaluated values of the parameters and respective profile values; and
# 3) a vector with optimal profile log-likelihood parameters and its maximum.
ciprof.par <- function(par, max.ll, step, conf, lltol, slopetol = 30, smoothtol = 10){

  # Loop through all GEV parameters
  ls.profile <- lapply(X = seq_along(par), FUN = function(i){

    # Find values to the left and to the right of the original 'par[i]' estimate
    left <- roots.proflik(-step, conf, par, i, max.ll)
    right <- roots.proflik(step, conf, par, i, max.ll)
    proflik <- rbind(left, right)
    proflik <- proflik[order(proflik[,i]),] # order profile
    slope <- proflik[,"slope"]           # extract profile slope
    smoothness <- abs(c(0, diff(slope))) # extract profile smoothness
    
    # Filter out any unstable ll estimates
    valid.idx <- which(abs(slope) < slopetol & smoothness < smoothtol)
    # proflik0 <- proflik
    proflik <- proflik[valid.idx,] 
    ci.lower <- min(proflik[,i]) # lower CI limit
    ci.upper <- max(proflik[,i]) # upper CI limit

    # TERMINAR DE ARRUMAR, A SEGUNDA RODADA ESTÁ RETORNANDO QUASE NENHUM VALOR

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
    profile <- proflik[,c(parname, "loglik", "slope")]
    # profile <- profile[order(profile[,1]),] # order by parameter value
    colnames(profile) <- c("par", "loglik", "slope")
    out <- list(
      ci = data.frame(
        gauge_code = gg,
        d = d,
        parname = parname,
        par = unname(par[i]),
        ci_lower = ci.lower,
        ci_upper = ci.upper,
        max = max.ll,
        lim = max.ll - qchisq(conf, 1)/2,
        prior = deparse(prior.info),
        message = message, # if message == 1, 'prof.pars' should be checked
        conf = conf,
        step = step
      ),
      prof.par = unlist(proflik[which.max(proflik[,"loglik"]), seq_len(npar + 1)]),
      profile = cbind(parname = parname, as.data.frame(profile))
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

# New negative log-likelihood funcion specific to quantile CI estimation
# @title (n - 1) negative log-likelihood function reparameterizzed to be a 
# function of GEV's quantile expression by replaxing $\xi$ with $\xi(q_p)$
# @details Takes the original negative log-likelihood and adapts it
# to optimize n - 1 parameters, keeping one of them fixed.
prof.nll.gev.q <- function(par.free, q0, p.nexceed, y, prior.info){

  p <- p.nexceed
  alpha <- par.free[[1]]
  kappa <- par.free[[2]]

  # if(abs(kappa) <= 1e-8) return(abs(kappa)*1e4)
  xi.q <- q0 - alpha/kappa*(1 - (-log(p))^kappa) # estimate a new 'xi' tied to 'q0'
  full <- c(xi.q, alpha, kappa)

  return(nll.gev(full, y, prior.info))

}

# Find roots of the profile log-likelihood to estimate quantile confidence intervals
# @title Find the 'roots' of the profile log-likelihood function
# @details This function optimizes `nll.prof.gev.q()` with respect to a disired quantile
# and estimates its Confidence Intervals to a given confidence level 'conf'.
roots.proflik.q <- function(step, conf, par, ll, rp, y, prior.info, lower, upper){

  p <- 1 - 1/rp                 # non-exceedance prob
  q0 <- qgev(unname(par), p)    # original quantile for 'p'
  last.par.free <- par[-1]      # initiate 'last.par.free'
  lim <- ll - qchisq(conf, 1)/2 # CI threshold
  result <- list()

  while(ll > lim){

    fit.prof <- tryCatch(
      optim(
        par = last.par.free, fn = function(internal.par){
          prof.nll.gev.q(par.free = internal.par, q0 = q0, p.nexceed = p, y = y, prior.info = prior.info)
        },
        gr = NULL, # allow L-BFGS-B to use its own estimated gradient
        method = "L-BFGS-B", lower = lower[-1], upper = upper[-1],
        control = list(maxit = maxit, factr = 1e9, pgtol = 1e-8, parscale = abs(par.free))
      ),
      error = function(e) NULL
    )

    if(!is.null(fit.prof) && fit.prof$convergence == 0){
      ll <- -fit.prof$value
    } else{
      fit.prof <- tryCatch(
        optim(
          par = last.par.free, fn = function(internal.par){
            prof.nll.gev.q(par.free = internal.par, q0 = q0, p.nexceed = p, y = y, prior.info = prior.info)
          },
          method = "Nelder-Mead", control = list(maxit = maxit, parscale = abs(last.par.free))
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
    gradient <- fit.prof[["counts"]][[2]]
    if(!is.null(fit.prof) && !is.na(ll)){
      last.par.free <- fit.prof$par
    }
  
    opt.alpha <- last.par.free[1] # extract optimal 'alpha'
    opt.kappa <- last.par.free[2] # extract optimal 'kappa'
    opt.xi <- q0 - opt.alpha/opt.kappa*(1 - (-log(p))^opt.kappa) # rebuild optimal 'xi'
    q0 <- q0 + step               # update quantile
    current.iter <- length(result) + 1 # update counter
    result[[current.iter]] <- c(q0, opt.xi, opt.alpha, opt.kappa, ll)
    message(sprintf(
      "q: %0.3f mm/h | ll: %0.3f | method: %s | convergence: %g",
      q0, ll, ifelse(is.na(gradient), "nm    ", "lbfgsb"), fit.prof[["convergence"]]
    ))

  }

  if (length(result) == 0) return(NULL)
  
  prof.par <- do.call(rbind, result)
  colnames(prof.par) <- c(paste0("q_", rp), "xi", "alpha", "kappa", "loglik")
  
  return(prof.par)

}

# 
ciprof.q <- function(par, return.period, max.ll, step, conf){

}

# MAIN FUNCTION ----------------------------------------------------------

fitprof.gev <- function(
  data, prior.info = NULL,
  durations, return.period,
  cols = c("gauge_code", "year", "duration", "imax"),
  conf = 0.95, lltol = 0.1, step = 0.01
){}

# TEST FUNCTION ----------------------------------------------------------

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
# gg <- "A801"
gg <- "A812"
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
# step <- 0.0025
step <- 0.005 # search step
# step <- 0.01


## PARÂMETROS ------------------------------------------------------------

# Function
# Original optimal parameters and maximum
fit <- fgev(y, prior.info)
par <- fit$par
npar <- length(par)
convergence <- fit$convergence
max.ll <- -fit$value

# First run
ls.ci.par <- ciprof.par(par, max.ll, step, conf, lltol, slopetol = 50, smoothtol = 15)
opt.par <- ls.ci.par[["ci"]][["par"]]
opt.ll <- max(ls.ci.par[["ci"]][["max"]])

# Run again with optimal 'new.par' parameters
rerun.ci.par <- ciprof.par(setNames(opt.par, names(par)), opt.ll, lltol, step, conf)
opt.par <- setNames(rerun.ci.par[["ci"]][["par"]], names(par))
opt.ll <- max(rerun.ci.par[["ci"]][["max"]])
lim.par <- opt.ll - qchisq(conf, 1)/2 # limite de definição dos ICs

## QUANTIS ---------------------------------------------------------------

# Define durations and return periods
# durations <- c(1, 6, 12, 24) # durations is already set as 'd'
return.period <- c(10, 20, 50)
probs <- 1 - 1/return.period

par.free <- opt.par[-1]
p <- probs[1]
rp <- return.period[1]
q0 <- qgev(opt.par, 1 - 1/rp)
ll <- opt.ll
# rp must have only one value

# Limites inferior e superior
lower <- c(-Inf, 1e-6, if (is.null(prior.info)) -Inf else -0.5 + 1e-6)
upper <- c(Inf, Inf, if (is.null(prior.info)) Inf else 0.5 - 1e-6)

# Setup an lapply loop to go through return periods
teste <- lapply(return.period, function(rp){
  
  # Find roots
  p <- 1 - 1/rp                 # later: this is being done twice: here and inside roots.proflik.q()
  q0 <- qgev(unname(opt.par), p)    # original quantile for 'p'
  left <- roots.proflik.q(-step, conf, opt.par, ll, rp, y, prior.info, lower, upper)
  right <- roots.proflik.q(step, conf, opt.par, ll, rp, y, prior.info, lower, upper)
  proflik <- rbind(left, right)
  
  q.name <- paste0("q_", rp)
  ci.lower <- min(proflik[,q.name]) # lower CI limit
  ci.upper <- max(proflik[,q.name]) # upper CI limit

  profile <- proflik[,c(q.name, "loglik")]
  colnames(profile) <- c("q", "loglik")
  out <- list(
    ci = data.frame(
      gauge_code = gg,
      d = d,
      return_period = rp,
      q = q0,
      ci_lower = ci.lower,
      ci_upper = ci.upper,
      max = max.ll,
      prior = deparse(prior.info),
      # algo = algo,
      conf = conf,
      step = step
    ),
    # prof.par = unlist(proflik[which.max(proflik[,"loglik"]), seq_len(npar + 1)]),
    profile = cbind(return_period = rp, as.data.frame(profile))
  )

  return(out)
  
})

proflik <- as_tibble(rbind(left, right))

ggplot(proflik, aes(x = q_10, y = loglik)) +
  geom_line() +
  # geom_hline()
  theme_bw()



# Figura
plot.par <- "alpha"
# profile <- rerun.ci.par$profile
profile <- ls.ci.par$profile
# ci.par <- rerun.ci.par$ci
ci.par <- ls.ci.par$ci
profile.par <- profile[profile$parname == plot.par,]
ggplot() +
  geom_vline(xintercept = ci.par[ci.par$parname == plot.par, "par"], linetype = "dashed", color = "grey50") +
  geom_line(aes(x = par, y = loglik), profile.par, linewidth = 0.3, color = "grey40") +
  geom_hline(yintercept = ci.par[ci.par$parname == plot.par, "lim"], color = "red") +
  geom_point(aes(x = par, y = loglik), profile.par, shape = 1, alpha = 0.3) +
  geom_point(aes(x = ci_lower, y = lim), ci.par[ci.par$parname == plot.par,], color = "red", size = 4, shape = 4) +
  geom_point(aes(x = ci_upper, y = lim), ci.par[ci.par$parname == plot.par,], color = "red", size = 4, shape = 4) +
  labs(
    x = plot.par, y = "loglik",
    subtitle = sprintf("Estação: %s | d: %g h", gg, d),
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Consolas"),
    panel.grid = element_blank()
  )
