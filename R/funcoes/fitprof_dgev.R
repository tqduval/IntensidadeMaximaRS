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
        control = list(maxit = 1e8, factr = 1e9, pgtol = 1e-8)
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
          method = "Nelder-Mead", control = list(maxit = 1e8)
        ),
        error = function(e) NULL
      )
      if(!is.null(fit.prof) && fit.prof$convergence %in% c(0, 10) && is.finite(fit.prof$value)){
        ll <- fit.prof$value
        if(ll < -1e6) ll <- NA_real_
      }
    }

    full <- numeric(npar)
    if(!is.null(fit.prof) && !is.na(ll)){
      full[-i] <- fit.prof$par
    } else{
      full[-i] <- last.par[-i]
    }
    full[i] <- val.i

    val.i <- val.i + step               # update 'val.i' for next iteration
    slope <- (ll - last.ll)/step        # profile surface slope (check for inconsistency)
    last.par <- full                    # update 'last.par' with current results
    current.iter <- length(result) + 1  # update counter

    result[[current.iter]] <- c(full, ll, slope) # save current results

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
ciprof.par <- function(y.id, ds, par, max.ll, conf, step, lltol, gg){

  lower <- c(-Inf, 1e-6, ifelse(is.null(prior.info), -Inf, -0.5 + 1e-6), 1e-6, 0)
  upper <- c(Inf, Inf, ifelse(is.null(prior.info), Inf, 0.5 - 1e-6), 1 - 1e-6, Inf)
  parnames <- c("xit", "alpha0", "kappa", "H", "theta")
  
  
  ls.profile <- lapply(X = seq_along(par), FUN = function(i){
  
    # Search left (lower CI limit)
    if(par[i] > lower[i] + step){
      left <- roots.proflik(y.id, ds, -step, conf, par, i, max.ll, lower, upper)
    } else{
      # message(sprintf("Optimal `%s` is at the edge of its domain`, returning `ci_lower == par[i]`", parnames[i]))
      left <- matrix(c(par, max.ll, 0), nrow = 1)
      left[1, i] <- lower[i]
      colnames(left) <- c(parnames, "loglik", "slope")
    }
    
    # Search right (upper CI limit)
    if(par[i] < upper[i] - step){
      right <- roots.proflik(y.id, ds, step, conf, par, i, max.ll, lower, upper)
    } else{
      # message(sprintf("Optimal `%s` is at the edge of its domain`, returning `ci_upper == par[i]`", parnames[i]))
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

roots.proflik.q <- function(step, q0, y.id, par, ds, d, p, conf, ll, lower, upper){

  npar <- length(par)
  last.par.free <- par[-1]
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
        control = list(maxit = 1e8, factr = 1e9, pgtol = 1e-8)
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
          method = "Nelder-Mead", control = list(maxit = 1e8)
        ),
        error = function(e) NULL
      )
      if(!is.null(fit.prof) && fit.prof$convergence %in% c(0, 10) && is.finite(fit.prof$value)){
        ll <- -fit.prof$value
        if(ll < -1e6){
          ll <- NA_real_
        }
      } else{
        ll <- NA_real_
      }
    }
    
    if(!is.null(fit.prof) && !is.na(ll)) last.par.free <- fit.prof$par
    
    # # Check GEV limits
    alpha0 <- last.par.free[[1]]
    kappa <- last.par.free[[2]]
    H <- last.par.free[[3]]
    theta <- last.par.free[[4]]

    A <- (d + theta)^H
    if(abs(kappa) < 1e-5){
      xi0 <- q0*A + alpha0 * log(-log(p)) 
    } else{
      xi0 <- q0*A - (alpha0 / kappa) * (1 - (-log(p))^kappa)
    }
    xit <- xi0 / alpha0
    
    current.iter <- length(result) + 1
    result[[current.iter]] <- c(d, q0, xit, last.par.free, ll)
    q0 <- q0 + step
  
    # message(sprintf(
    #   "%s | p: %0.3f | d: %2d | q: %0.3f mm/h | ll: %0.3f | conv: %g",
    #   gg, p, d, q0, ll, fit.prof[["convergence"]]
    # ))
  
  }
  
  rp <- 1/(1 - p)
  if(length(result) == 0) return(NULL)
  prof.par <- do.call(rbind, result)
  colnames(prof.par) <- c("d", sprintf("q_%g", rp), parnames, "loglik")
  
  return(prof.par)

}

# Estimate profile likelihood confidence intervals for quantiles
ciprof.q <- function(y.id,  ds, durations, q0, return.period, par, max.ll, conf, step, lltol, gg){

  # Parameter bounds
  lower <- c(-Inf, 1e-6, ifelse(is.null(prior.info), -Inf, -0.5 + 1e-6), 1e-6, 0)
  upper <- c(Inf, Inf, ifelse(is.null(prior.info), Inf, 0.5 - 1e-6), 1 - 1e-6, Inf)
  
  # All combinations of duration and return period
  if(is.null(durations)){
    grid <- expand.grid(d = ds, rp = return.period) # use durations from the data given
  } else{
    grid <- expand.grid(d = durations, rp = return.period) # use user defined durations
  }

  # Loop through duration and return period combinations
  ls.profile <- lapply(seq_len(nrow(grid)), function(i){
  
    d <- grid[i, "d"]
    rp <- grid[i, "rp"]
    p <- 1 - 1/rp # non-excceedence probability
    q0 <- q.dgev.vec(par, d, p) # original quantile estimate
    step <- if(q0 != 0) step*abs(q0)*0.1 else step
    
    out <- tryCatch({

      left <- roots.proflik.q(-step, q0, y.id, par, ds, d, p, conf, max.ll, lower, upper) # search right
      right <- roots.proflik.q(step, q0, y.id, par, ds, d, p, conf, max.ll, lower, upper) # search left
  
      proflik <- rbind(left, right)
      q.name <- sprintf("q_%g", rp)
      ci.lower <- min(proflik[,q.name], na.rm = TRUE)
      ci.upper <- max(proflik[,q.name], na.rm = TRUE)
      censored <- if(any(is.na(proflik[,"loglik"]))) 1 else 0
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
          censored_ci = censored,
          conf = conf,
          step = step,
          error = 0L,
          error_message = NA_character_
        ),
        # prof.par = unlist(proflik[which.max(proflik[,"loglik"]), seq_len(npar + 1)]),
        profile = cbind(
          gauge_code = gg, d = d,
          return_period = rp, as.data.frame(profile),
          error = 0L
        )
      )
      
      return(out)

    }, error = function(e){

      list(
        ci = data.frame(
          gauge_code = gg, 
          d = d,
          return_period = rp, 
          q = q0,
          ci_lower = NA_real_, 
          ci_upper = NA_real_,
          max = max.ll, 
          lim = max.ll - qchisq(conf, 1)/2,
          prior = deparse(prior.info),
          n_profle = NA_integer_,
          flag_shape = if(abs(par[3]) > 0.5) 1 else 0, # flag '1' whenever |kappa| is larger than 0.
          censored_ci = 1L,
          conf = conf,
          step = step,
          error = 1L,
          error_message = conditionMessage(e)
        ),
        profile = data.frame(
          gauge_code = gg,
          d = d,
          return_period = rp,
          q = NA_real_,
          loglik = NA_real_,
          error = 1L
        )
      )

    })

    return(out)
    
  })
  
  results <- bind_rows(lapply(ls.profile, `[[`, "ci"))
  profile <- bind_rows(lapply(ls.profile, `[[`, "profile"))
  
  return(list(ci = results, profile = profile))

}

# MAIN FUNCTION ----------------------------------------------------------


fitprof.dgev <- function(
  data, prior.info = NULL,
  return.period, durations = NULL,
  cols = c("gauge_code", "year", "duration", "imax"),
  conf = 0.95, lltol = 0.005, step = 0.005
){

  #### CHECAGENS ----

  if(!is.list(data) || length(data) == 0){
    stop("`data` must be a non-empty list with data.frames.")
    warning("`cols` must contain column names matching the following information in the same order: 'gauge_code', 'year', 'duration', 'imax")
  }
  if(!is.numeric(return.period) || length(return.period) == 0) stop("`return.period` must have one or more numeric values.")
  if(!length(prior.info) %in% c(0,2)) stop("`prior.info` must be either 'NULL' or a vector of length == 2 with regional mean and standard deviation on shape")
  if(step <= 0) stop("A positive 'step' must be defined.")
  if(!is.numeric(conf) || length(conf) != 1 || conf <= 0 || conf >= 1) stop("Confidence level `conf` must be between (0,1).")
  if(!is.numeric(lltol) || length(lltol) != 1 || lltol <= 0) stop("Profile likelihood change tolerance `lltol` musc be a positive number.")
  data.cols <- colnames(data[[1]])
  missing.cols <- setdiff(cols, data.cols)
  if(length(missing.cols) > 0) stop(sprintf("Missing columns in `data`: %s", paste(missing.cols, collapse = ", ")))
  
  #### CÁLCULOS ----

  message("Fitting d-GEV distribution, estimating quantiles and confidence intervals with 'profile likelihood'.")
  message(">> No. gauges: ", length(data))
  message(">> Durations [h]: ", paste(if(is.null(durations)) "all available" else durations, collapse = ", "))
  message(">> Return periods [years]: ", paste(return.period, collapse = ", "))
  message(">> Search step: ", step)
  message(">> Confidence level for profile likelihood CIs: ", paste0(conf*100, "%"))

  ls.res <- pbapply::pblapply(data, function(data.gg){

    gg <- unique(data.gg$gauge_code)

    res <- tryCatch({

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
  
      # Confidence intervals
      ls.ci.par <- ciprof.par(y.id, ds, par, max.ll, conf, step, lltol, gg)
      ls.ci.q <- ciprof.q(y.id,  ds, durations, q0, return.period, par, max.ll, conf, step, lltol, gg)
  
      return(list(
        gauge_code = gg, 
        par = ls.ci.par, 
        quantiles = ls.ci.q, 
        error = 0L,
        error_message = NA_character_
      ))

    }, error = function(e){

      return(list(
        gauge_code = gg,
        par = NULL,
        quantiles = NULL,
        error = 1L,
        error_message = conditionMessage(e)
      ))

    })

    return(res)

  }) # fim 'ls.res'

  ls.par <- lapply(ls.res, `[[`, "par")
  ls.quantiles <- lapply(ls.res, `[[`, "quantiles")
  ci.par <- bind_rows(lapply(ls.par, `[[`, "ci")); rownames(ci.par) <- NULL
  ci.q <- bind_rows(lapply(ls.quantiles, `[[`, "ci"))
  prof.par <- bind_rows(lapply(ls.par, `[[`, "profile"))
  prof.q <- bind_rows(lapply(ls.quantiles, `[[`, "profile"))

  return(list(
    par = ci.par,
    idf = ci.q,
    prof.par = prof.par,
    prof.q = prof.q
  ))

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
# durations <- c(1, 2, 4, 6, 8, 10, 12, 18, 24)
durations <- c(1, 6, 12, 24)
return.period <- c(1.1, 1.5, 2, 5, 10, 20, 50, 100, 200, 500)
conf <- 0.95
maxit <- 1e8
step <- 0.005
lltol <- 0.005

# Rodar função
fitprof.ls <- fitprof.dgev(
  data = imax.ls,
  prior.info = prior.info,
  return.period = return.period,
  durations = NULL, # roda com todas as durações disponíveis
  cols = cols
)

saveRDS(fitprof.ls, "base/gerados/ls_fit_dgev_proflik.rds")



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
