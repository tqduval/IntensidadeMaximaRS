rm(list = ls()); invisible(gc())

# INTERNAL FUNCTIONS -----------------------------------------------------

# Fit d-GEV distribution with maximum likelihood or generalized maximum likelihood
source("R/funcoes/fdgev.R")

# d-GEV negative log-likelihood function for all but one parameter
prof.nll.dgev <- function(par.free, which.par, par.value, y.id, ds, prior.info){

  full <- numeric(length(par.free) + 1) # build 'new' parameters vector
  full[which.par] <- par.value          # fixed value for current parameter
  full[-which.par] <- par.free          # remaining parameters to be optimized

  nll.dgev(full, ds, y.id, prior.info)

}

# Gradient function for d-GEV negative log-likelihood for all but one parameter
prof.grad.dgev <- function(par.free, which.par, par.value, y.id, durations, prior.info){

  full <- numeric(length(par.free) + 1) # build 'new' parameters vector
  full[which.par] <- par.value          # fixed value for current parameter
  full[-which.par] <- par.free          # remaining parameters to be optimized
  gr <- grad.ll.dgev(full, durations, y.id, prior.info)
  
  gr[-which.par]

}

# Find the 'roots' of the profile log-likelihood function
roots.proflik <- function(y.id, ds, step, conf, par, i, ll, lower, upper){

  val.i <- par[[i]]
  last.par <- par
  npar <- length(par)
  lim <- ll - qchisq(conf, 1)/2
  if(val.i != 0) step <- step*par[[i]] # use a dynamically defined step
  parnames <- c("xit", "alpha0", "kappa", "H", "theta")

  result <- list()

  while(ll > lim && !is.na(ll)){

    last.ll <- ll
    if(i == 2){
      if(val.i <= lower[i]){
        ll <- NA_real_
        next
      }
    }
    if(i == 5 && step < 0){
      if(val.i <= lower[i]){
        ll <- NA_real_
        next
      }
    }

    fit.prof <- tryCatch(
      optim(
        par = last.par[-i],
        fn = function(internal.par){
          prof.nll.dgev(
            par.free = internal.par, which.par = i, par.value = val.i,
            y.id = y.id, ds = ds, prior.info = prior.info
          )
        },
        gr = NULL, method = "L-BFGS-B", lower = lower[-i], upper = upper[-i],
        control = list(maxit = maxit, factr = 1e9, pgtol = 1e-8)
      ), error = function(e) NULL
    )

    if(!is.null(fit.prof) && fit.prof$convergence == 0){

      ll <- -fit.prof$value # save new profile value 'll'

    } else{

      fit.prof <- tryCatch(
        optim(
          par = last.par[-i],
          fn = function(internal.par){
            prof.nll.dgev(
              par.free = internal.par, which.par = i, par.value = val.i,
              y.id = y.id, ds = ds, prior.info = prior.info
            )
          },
          method = "Nelder-Mead", control = list(maxit = maxit)
        )
      )
    }

    full <- numeric(npar)
    if(!is.null(fit.prof) && !is.na(ll)){
      full[-i] <- fit.prof$par
    } else{
      full[-i] <- last.par[-i]
    }
    full[i] <- val.i

    # Check GEV limits
    xit <- full[1]
    alpha0 <- full[2]
    kappa <- full[3]
    H <- full[4]
    theta <- full[5]

    A <- (ds + theta)^H
    xi0 <- xit*alpha0
    y <- sweep(y.id, 2, A, "*")      # this is the Koutsoyiannis et al. (1998) rescaled variable
    zi <- 1 - kappa/alpha0*(y - xi0) # GEV internal variable
    if(any(zi < 1e-6)) ll <- NA_real_

    # val.i <- val.i + step.i             # update 'val.i' for next iteration
    val.i <- val.i + step               # update 'val.i' for next iteration
    slope <- (ll - last.ll)/step        # profile surface slope (check for inconsistency)
    last.par <- full                    # update 'last.par' with current results
    current.iter <- length(result) + 1  # update counter

    result[[current.iter]] <- c(full, ll, slope) # save current results
    # result[[current.iter]] <- data.frame(t(setNames(full, parnames)), ll = ll, slope = slope)

    # message(sprintf(
    #   "%s | ci %-3s | %s: %0.3f | ll: %0.3f | slope: %0.3f",
    #   gg, ifelse(step < 0, "low", "up"), parnames[i], val.i, ll, slope
    # ))

  } # end of 'while' search

  prof.par <- do.call(rbind, result)
  colnames(prof.par) <- c(parnames, "loglik", "slope")

  return(prof.par)

}

# Estimate profile likelihood confidence intervals for d-GEV parameters
ciprof.par <- function(y.id, ds, par, max.ll, conf, step, lltol){

  lower <- c(-Inf, 1e-6, ifelse(is.null(prior.info), -Inf, -0.5 + 1e-6), 1e-6, 0)
  upper <- c(Inf, Inf, ifelse(is.null(prior.info), Inf, 0.5 - 1e-6), 1 - 1e-6, Inf)
  parnames <- c("xit", "alpha0", "kappa", "H", "theta")
  
  
  ls.profile <- lapply(X = seq_along(par), FUN = function(i){
  
    # Search left (lower CI limit)
    if(par[i] > lower[i] + step){
      left <- roots.proflik(y.id, ds, -step, conf, par, i, max.ll, lower, upper)
    } else{
      message(sprintf("Optimal `%s` is at the edge of its domain`, returning `ci_lower == par[i]`", parnames[i]))
      left <- matrix(c(par, max.ll, 0), nrow = 1)
      left[1, i] <- lower[i]
      colnames(left) <- c(parnames, "loglik", "slope")
    }
    
    # Search right (upper CI limit)
    if(par[i] < upper[i] - step){
      right <- roots.proflik(y.id, ds, step, conf, par, i, max.ll, lower, upper)
    } else{
      message(sprintf("Optimal `%s` is at the edge of its domain`, returning `ci_upper == par[i]`", parnames[i]))
      right <- matrix(c(par, max.ll, 0), nrow = 1)
      right[1, i] <- upper[i]
      colnames(right) <- c(parnames, "loglik", "slope")
    }
    
    proflik <- rbind(left, right)
    proflik <- proflik[order(proflik[,i]),]
    ci.lower <- min(proflik[,i], na.rm = TRUE)
    ci.upper <- max(proflik[,i], na.rm = TRUE)
    ll.max.proflik <- max(proflik[,"loglik"])
    par.max.proflik <- proflik[which.max(proflik[,"loglik"]), seq_len(length(par))]
    if(par[i] != 0){
      if(abs(par.max.proflik[i] - par[i])/par[i] > lltol){
        message <- 1               # new max found
        out.ll <- ll.max.proflik   # update max ll
        out.par <- par.max.proflik # update parameters
      } else{
        message <- 0     # keeping same max
        out.ll <- max.ll # same max ll
        out.par <- par   # same parameters
      }
    } else{
      message <- 0     # keeping same max
      out.ll <- max.ll # same max ll
      out.par <- par   # same parameters
    }
    
    # plot(proflik[,c(parnames[i], "loglik")], type = "p")
    # segments(x0 = par[i], x1 = par[i], y0 = max.ll, y1 = max.ll - qchisq(conf, 1)/2, lty = "dashed", col = "grey50")
    # segments(min(proflik[,i], na.rm = TRUE),  
    #   max(proflik[,i], na.rm = TRUE), 
    #   max.ll - qchisq(conf, 1)/2,
    #   max.ll - qchisq(conf, 1)/2,
    #   col = "red"
    # )
  
    profile <- proflik[,c(parnames[i], "loglik", "slope")]
    colnames(profile) <- c("par", "loglik", "slope")
    out <- list(
      ci = data.frame(
        gauge_code = gg,
        parname = parnames[i],
        par = unname(par[i]),
        ci_lower = ci.lower,
        ci_upper = ci.upper,
        ci_type = "proflik",
        rel_width = (ci.upper - ci.lower)/abs(out.par[i]),
        max_ll = max.ll,                  # original max ll
        lim = max.ll - qchisq(conf, 1)/2, # original chi-squared limit
        found_max = out.ll, # max(max.ll, max(proflik[,"loglik"]))
        prior = deparse(prior.info),
        n_profile = nrow(proflik),
        message = message, # if message == 1, 'prof.pars' should be checked
        conf = conf,
        step = step
      ),
      prof.par = unlist(proflik[which.max(proflik[,"loglik"]), seq_len(length(par) + 1)]), # + 1 pega o loglik
      # prof.par = out.par,
      profile = cbind(
        gauge_code = gg, # keep as comment to avoid having to call it inside every function
        parname = parnames[i], as.data.frame(profile)
      )
    )
  
    return(out)
  
  })
  
  results <- bind_rows(lapply(ls.profile, `[[`, "ci"))
  prof.par <- bind_rows(lapply(ls.profile, `[[`, "prof.par"))
  profile <- bind_rows(lapply(ls.profile, `[[`, "profile"))
  
  if(any(results$message == 1)){
    new.par <- unlist(prof.par[which.max(prof.par$loglik), seq_len(length(par))])
    new.max <- max(profile[,"loglik"], na.rm = TRUE)
    results$par <- new.par
    results$max <- new.max
  }
  
  return(list(ci = results, profile = profile))

}

# d-GEV negative log-likelihood function specific to quantile CI estimation
prof.nll.dgev.q <- function(
  par.free, q0, p.nexceed, y.id, 
  ds, # vetor com durações (nome coluna y.id)
  duration, # escalar (1 valor somente) para calcular o quantil
  prior.info
){

  p <- p.nexceed; d <- duration
  alpha0 <- par.free[[1]]
  kappa <- par.free[[2]]
  H <- par.free[[3]]
  theta <- par.free[[4]]

  A <- (d + theta)^H
  alpha <- alpha0/A

  if(abs(kappa) < 1e-5){
    xi0 <- q0*A + alpha0*log(-log(p)) # calcular quantil Gumbel p/ evitar erro numérico
  } else{
    xi0 <- q0*A - alpha0/kappa*(1 - (-log(p))^kappa)
  }

  xit <- xi0/alpha0
  full <- c(xit, par.free)

  return(nll.dgev(full, ds, y.id, prior.info))

}

# Vectorized quantile function
q.dgev.vec <- function(param, d, p){
  
  xit <- param[1]; alpha0 <- param[2]; kappa <- param[3]; H <- param[4]; theta <- param[5]
  A <- alpha0/(d + theta)^H
  B <- xit + (1/kappa)*(1 - (-log(p))^kappa)
  
  yp <- A*B
  
  return(unname(yp))
  
}

roots.proflik.q <- function(step, q0, y.id, par, ds, d, rp, conf, ll, lower, upper){

  npar <- length(par)
  last.par.free <- par[-1]
  # step <- step*q0
  lim <- ll - qchisq(conf, 1)/2
  parnames <- c("xit", "alpha0", "kappa", "H", "theta")

  result <- list()
  
  while(ll > lim && !is.na(ll)){
  
    fit.prof <- tryCatch(
      optim(
        par = last.par.free, fn = function(internal.par){
          prof.nll.dgev.q(
            par.free = internal.par,
            q0 = q0, p.nexceed = p, y.id = y.id,
            ds = ds, duration = d,
            prior.info = prior.info
          )
        },
        gr = NULL, method = "L-BFGS-B", lower = lower[-1], upper = upper[-1],
        control = list(maxit = maxit, factr = 1e9, pgtol = 1e-8)
      ),
      error = function(e) NULL
    )
    
    if(!is.null(fit.prof) && fit.prof$convergence == 0){
      ll <- -fit.prof$value
    } else{
      fit.prof <- tryCatch(
        optim(
          par = last.par.free, fn = function(internal.par){
            prof.nll.dgev.q(
              par.free = internal.par,
              q0 = q0, p.nexceed = p, y.id = y.id,
              ds = ds, duration = d,
              prior.info = prior.info
            )
          },
          method = "Nelder-Mead", control = list(maxit = maxit)
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
    
    if(!is.null(fit.prof) && !is.na(ll)) last.par.free <- fit.prof$par
    
    # Check GEV limits
    alpha0 <- last.par.free[[1]]
    kappa <- last.par.free[[2]]
    H <- last.par.free[[3]]
    theta <- last.par.free[[4]]
    
    A <- (ds + theta)^H # `A` calculado com o vetor p/ conferir a distribuição
    if(abs(kappa) < 1e-5){
      xi0 <- q0*A + alpha0*log(-log(p)) # calcular quantil Gumbel p/ evitar erro numérico
    } else{
      xi0 <- q0*A - alpha0/kappa*(1 - (-log(p))^kappa)
    }
    y <- sweep(y.id, 2, A, "*")      # this is the Koutsoyiannis et al. (1998) rescaled variable
    zi <- 1 - kappa/alpha0*(y - xi0) # GEV internal variable
    if(any(zi < 1e-6)) ll <- NA_real_
  
    A <- (d + theta)^H # `A` calculado como escalar p/ atualizar o valor de `xit` no while
    if(abs(kappa) < 1e-5){
      xi0 <- q0*A + alpha0*log(-log(p)) # calcular quantil Gumbel p/ evitar erro numérico
    } else{
      xi0 <- q0*A - alpha0/kappa*(1 - (-log(p))^kappa)
    }
    xit <- xi0/alpha0 # agora como `xi0` é um escalar, `xit` também é
    
    current.iter <- length(result) + 1
    result[[current.iter]] <- c(d, q0, xit, last.par.free, ll)
    q0 <- q0 + step
  
    message(sprintf(
      "%s | rp: %g | d: %g | q: %0.3f mm/h | ll: %0.3f | conv: %g",
      gg, rp, d, q0, ll, fit.prof[["convergence"]]
    ))
  
  }
  
  if(length(result) == 0) return(NULL)
  prof.par <- do.call(rbind, result)
  colnames(prof.par) <- c("d", paste0("q_", rp), parnames, "loglik")
  
  return(prof.par)

}

# MAIN FUNCTION ----------------------------------------------------------


fitdelta.dgev <- function(
  data, prior.info = NULL,
  return.period,
  cols = c("gauge_code", "year", "duration", "imax"),
  conf = 0.95, lltol = 0.005, step = 0.005
){

  if(step <= 0) stop("A positive 'step' must be defined.")

}


# TESTING THE FUNCTION ---------------------------------------------------

pacman::p_load(pacman, tidyverse, arrow)
imax.df <- arrow::read_parquet("base/gerados/df_imax.pqt")

# Filter dada
na.yr <- 0.2
min.durations <- 3
d.interval <- c(1, 24)
min.yr <- 8

imax.df <- imax.df |> 
  filter(na_prct_yr <= na.yr,
         between(d, d.interval[1], d.interval[2])) |> 
  group_by(gauge_code) |> 
  filter(n_distinct(d) >= min.durations,
         n_distinct(wateryear) >= min.yr) |> 
  ungroup() |> 
  select(gauge_code, wateryear, d, imax)

imax.ls <- split(imax.df, imax.df$gauge_code)

# Define function arguments
which.prior <- list(NULL, c(0.0178, 0.0627))
prior.info <- which.prior[[1]]
cols <- names(imax.df)
method <- "L-BFGS-B"
scale.inv <- "wide"
durations <- c(1, 2, 4, 6, 8, 10, 12, 18, 24)
return.period <- c(1.1, 1.5, 2, 5, 10, 20, 50, 100, 200, 500)
conf <- 0.95
maxit <- 1e8
step <- 0.005
lltol <- 0.005

## INSIDE THE FUNCTION ---------------------------------------------------

# Everything below this is inside a 'pblapply()' loop
# data.gg <- imax.ls[[1]]
data.gg <- imax.ls[["A812"]]
gg <- unique(data.gg$gauge_code)

y.id <- data.gg |> 
  pivot_wider(
    names_from = cols[3],  # turn durations into columns
    values_from = cols[4], # fill with 'imax' values
    names_prefix = "d_"
  ) |> 
  arrange(cols[2]) |>   # arrange by 'year' (or 'wateryear')
  select(-cols[1:2]) |> # remove 'gauge_code' and 'year' columns
  na.omit() |> 
  as.matrix()

ds <- as.integer(gsub("d_", "", colnames(y.id))) # duration vector

# This case does not require a loop thorugh all durations as they
# are alll fitted in the same model
fit <- fdgev(
  y.id = y.id,
  prior.info = prior.info,
  cols = names(data.gg),
  method = method,
  scale.inv = scale.inv
)

par <- fit$par
npar <- length(par)
convergence <- fit$convergence
max.ll <- -fit$value

ls.ci.par <- ciprof.par(y.id, ds, par, max.ll, conf, step, lltol)

# início roots.proflik.q()
lower <- c(-Inf, 1e-6, ifelse(is.null(prior.info), -Inf, -0.5 + 1e-6), 1e-6, 0)
upper <- c(Inf, Inf, ifelse(is.null(prior.info), Inf, 0.5 - 1e-6), 1 - 1e-6, Inf)

if(is.null(durations)){
  grid <- expand.grid(d = ds, rp = return.period) # use durations from the data given
} else{
  grid <- expand.grid(d = durations, rp = return.period) # use user defined durations
}

# entra em um loop do comprimento de `grid`
ls.profile <- lapply(seq_len(nrow(grid)), function(i){
  
  # d <- grid[i, "d"]
  # rp <- grid[i, "rp"]
  d <- 20
  rp <- 20
  p <- 1 - 1/rp
  q0 <- q.dgev.vec(par, d, p)
  
  left <- roots.proflik.q(-step, q0, y.id, par, ds, d, rp, conf, max.ll, lower, upper)
  right <- roots.proflik.q(step, q0, y.id, par, ds, d, rp, conf, max.ll, lower, upper)
  proflik <- rbind(left, right)
  q.name <- paste0("q_", rp)
  ci.lower <- min(proflik[,q.name], na.rm = TRUE)
  ci.upper <- max(proflik[,q.name], na.rm = TRUE)
  neval <- nrow(proflik)
  
  if(neval <= 20){
    message <- 1 # possible error, few profile values
  } else{
    message <- 0 # might be ok, but check anyway
  }
  
  profile <- proflik[order(proflik[,q.name]),][,c(q.name, "loglik")]
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
      lim = max.ll - qchisq(conf, 1)/2,
      prior = deparse(prior.info),
      n_profile = nrow(proflik),
      flag_shape = if(abs(par[3]) > 0.5) 1 else 0, # flag '1' whenever |kappa| is larger than 0.
      conf = conf,
      step = step
    ),
    # prof.par = unlist(proflik[which.max(proflik[,"loglik"]), seq_len(npar + 1)]),
    profile = cbind(
      gauge_code = gg, d = d,
      return_period = rp, as.data.frame(profile)
    )
  )
  
  return(out)
  
})



quantile.names <- colnames(proflik)
proflik <- proflik[order(proflik[,2]),][,quantile.names[c(2,length(quantile.names))]]
plot(proflik, type = "p", xlab = sprintf("%s (d = %g h)", quantile.names[2], d), main = gg)
segments(x0 = q0, x1 = q0, y0 = max.ll, y1 = max.ll - qchisq(conf, 1)/2, lty = "dashed", col = "grey50")
segments(
  x0 = min(proflik[,1], na.rm = TRUE), 
  y0 = max.ll - qchisq(conf, 1)/2,
  x1 = max(proflik[,1], na.rm = TRUE), 
  y1 = max.ll - qchisq(conf, 1)/2,
  col = "red"
)
