
rm(list = ls()); gc()

pacman::p_load(pacman, tidyverse, arrow)

df.imax <- arrow::read_parquet("base/gerados/df_imax.pqt")

data <- df.imax %>%
  filter(na_prct_yr <= 0.2, between(d, 1, 24)) %>%
  group_by(gauge_code) %>%
  filter(n_distinct(d) >= 3, n_distinct(wateryear) >= 8) %>%
  ungroup() %>%
  select(gauge_code, wateryear, d, imax)

prior.info <- NULL
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

conf <- 0.95
lltol <- 0.1
step <- 0.005 # search step


fit <- fgev(y, prior.info)
par <- fit$par
parnames <- names(par)

# New 'roots.proflik()' with restrictions on fixed parameters
# to avoid 'while()' trying values that violate GEV limits
# i <- 1
# i <- 2 # alpha
# i <- 3

# This is inside `roots,prpolik()`
val.i <- par[[i]]
npar <- length(par)
last.par <- par
nll <- fit$value
ll <- -nll
lim <- ll - qchisq(conf, 1)/2

result <- list()

lower <- c(-Inf, 1e-6, if (is.null(prior.info)) -Inf else -0.5 + 1e-6)
upper <- c(Inf, Inf, if (is.null(prior.info)) Inf else 0.5 - 1e-6)

while(ll > lim && !is.na(ll)){

  last.ll <- ll
  
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
  delta.ll <- (ll - last.ll)/step
  last.par <- full
  current.iter <- length(result) + 1
  result[[current.iter]] <- c(full, ll, delta.ll)
  
  message(
    sprintf(
      "ci %-5s | %s: %0.3f | loglik: %0.3f", 
      ifelse(step < 0, "lower", "upper"), 
      names(par)[i], val.i, ll)
    )
  
} # end of 'while' search

do.call(rbind, result)


roots.proflik <- function(step, conf, par, i, ll){
  
  val.i <- par[[i]]
  npar <- length(par)
  last.par <- par
  lim <- ll - qchisq(conf, 1)/2
  
  result <- list()
  
  lower <- c(-Inf, 1e-6, if (is.null(prior.info)) -Inf else -0.5 + 1e-6)
  upper <- c(Inf, Inf, if (is.null(prior.info)) Inf else 0.5 - 1e-6)
  
  while(ll > lim && !is.na(ll)){
    
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
    last.par <- full
    current.iter <- length(result) + 1
    result[[current.iter]] <- c(full, ll)
    
    message(
      sprintf(
        "ci %-5s | %s: %0.3f | loglik: %0.3f", 
        ifelse(step < 0, "lower", "upper"), 
        names(par)[i], val.i, ll)
    )
    
  } # end of 'while' search
  
  prof.par <- do.call(rbind, result)
  colnames(prof.par) <- c(names(par), "loglik")
  
  return(prof.par)
  
}

max.ll <- -fit$value

i <- 3
left <- roots.proflik(-step, conf, par, i, max.ll)
right <- roots.proflik(step, conf, par, i, max.ll)
proflik <- rbind(left, right)

proflik.par <- proflik[,c(parnames[i], "loglik")]
ci.par <- data.frame(
  x = c(
    min(proflik.par[,parnames[i]]),
    max(proflik.par[,parnames[i]])
  ),
  y = rep(max.ll - qchisq(conf, 1) / 2, 2)
)

# Figura
ggplot() +
  geom_vline(xintercept = par[i], linetype = "dashed", color = "grey50") +
  geom_line(aes(x = proflik[, i], y = proflik[, "loglik"]), linewidth = 0.3, color = "grey40") +
  geom_hline(yintercept = max.ll - qchisq(conf, 1) / 2, color = "red") +
  geom_point(aes(x = proflik[, i], y = proflik[, "loglik"]), shape = 1, alpha = 0.3) +
  geom_point(aes(x = x, y = y), ci.par, color = "red", size = 4, shape = 4) +
  labs(
    x = parnames[i], y = "loglik",
    subtitle = sprintf("Estação: %s | d: %g h", gg, d),
  ) +
  theme_bw() +
  theme(
    text = element_text(family = "Consolas"),
    panel.grid = element_blank()
  )
