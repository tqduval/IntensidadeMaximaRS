
# rm(list = ls()); invisible(gc())

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
roots.proflik <- function(y, step, conf, par, i, ll, lower, upper){
  
  val.i <- par[[i]]
  npar <- length(par)
  last.par <- par
  lim <- ll - qchisq(conf, 1)/2
  parnames <- c("xi", "alpha", "kappa")
  
  result <- list()
  
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
    
    # message(
    #   sprintf(
    #     "%s | ci %-5s | %s: %0.3f | ll: %0.3f | slope: %0.3f", 
    #     gg, ifelse(step < 0, "lower", "upper"), 
    #     parnames[i], val.i, ll, slope)
    # )
    
  } # end of 'while' search
  
  prof.par <- do.call(rbind, result)

  colnames(prof.par) <- c(parnames, "loglik", "slope")
  
  return(prof.par)
  
}

# @title Estimate profile likelihood confidence intervals for a parameter
# @details This functions calls `roots.proflik()` and estimates confidence intervals
# for a given parameter. It's called inside an `lapply()` loop and returns a list with:
# 1) a row tbl_df with original parameters from the log-likelihood maximization and 
# estimated confidence intervals;
# 2) a matrix with all evaluated values of the parameters and respective profile values; and
# 3) a vector with optimal profile log-likelihood parameters and its maximum.
ciprof.par <- function(y, par, max.ll, step, conf, lltol, gg, d){

  lower <- c(-Inf, 1e-6, if (is.null(prior.info)) -Inf else -0.5 + 1e-6)
  upper <- c(Inf, Inf, if (is.null(prior.info)) Inf else 0.5 - 1e-6)
  parnames <- c("xi", "alpha", "kappa")

  # Loop through all GEV parameters
  ls.profile <- lapply(X = seq_along(par), FUN = function(i){

    # Find values to the left and to the right of the original 'par[i]' estimate
    left <- roots.proflik(y, -step, conf, par, i, max.ll, lower, upper)
    right <- roots.proflik(y, step, conf, par, i, max.ll, lower, upper)
    proflik <- rbind(left, right)
    proflik <- proflik[order(proflik[,i]),] # order profile
    slope <- proflik[,"slope"]   # extract profile slope
    ci.lower <- min(proflik[,i]) # lower CI limit
    ci.upper <- max(proflik[,i]) # upper CI limit

    # Update optimal parameters if proflik finds a higher optimal point
    # ll.max.proflik <- max(proflik[,"loglik"])
    # ll.dif <- abs((ll.max.proflik - max.ll)/max.ll) <= lltol
    # if(isFALSE(ll.dif)){
    #   message <- 1 # updated parameters
    #   out.par <- unlist(proflik[which.max(proflik[,"loglik"]), seq_len(npar)])
    # } else{
    #   message <- 0 # original fgev() parameters
    #   out.par <- par
    # }

    # Update optimal parameters if proflik finds a higher optimal point
    ll.max.proflik <- max(proflik[,"loglik"])
    par.max.proflik <- proflik[which.max(proflik[,"loglik"]), seq_len(length(par))]
    if(abs(par.max.proflik[i] - par[i])/par[i] > lltol){
      message <- 1               # new max found
      out.ll <- ll.max.proflik   # update max ll
      out.par <- par.max.proflik # update parameters
    } else{
      message <- 0     # keeping same max
      out.ll <- max.ll # same max ll
      out.par <- par   # same parameters
    }

    profile <- proflik[,c(parnames[i], "loglik", "slope")]
    colnames(profile) <- c("par", "loglik", "slope")
    out <- list(
      ci = data.frame(
        gauge_code = gg,
        d = d,
        parname = parnames[i],
        par = unname(par[i]),
        ci_lower = ci.lower,
        ci_upper = ci.upper,
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
        d = d,
        parname = parnames[i], as.data.frame(profile)
      )
    )

    return(out)

  }) # fim 'ls.profile'

  # names(ls.profile) <- c("xi", "alpha", "kappa") # name parameters

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

    # Check GEV limits
    gev.lim <- opt.xi + opt.alpha/opt.kappa  # xi + alpha/kappa
    if(opt.kappa < 0){
      if(any(y < gev.lim)) ll <- NA_real_ # domain violation for Fréchet law
    }
    if(opt.kappa > 0){
      if(any(y > gev.lim)) ll <- NA_real_ # domain violation for reverse Weibull law
    }

    current.iter <- length(result) + 1 # update counter
    result[[current.iter]] <- c(q0, opt.xi, opt.alpha, opt.kappa, ll) # save current iteration
    q0 <- q0 + step                    # update quantile value

    # message(sprintf(
    #   "%s | rp: %g | q: %0.3f mm/h | ll: %0.3f | conv: %g",
    #   gg, rp, q0, ll,
    #   # ifelse(is.na(gradient), "nm    ", "lbfgsb"), 
    #   fit.prof[["convergence"]]
    # ))

  }

  if (length(result) == 0) return(NULL)
  
  prof.par <- do.call(rbind, result)
  colnames(prof.par) <- c(paste0("q_", rp), "xi", "alpha", "kappa", "loglik")
  
  return(prof.par)

}

# @title Estimate profile likelihood confidence intervals for quantiles
# @details This functions calls `roots.proflik.q()` and estimates confidence intervals
# for quantiles It's called inside an `lapply()` loop and returns a list with:
# 1) a row tbl_df with original parameters from the log-likelihood maximization and 
# estimated confidence intervals;
# 2) a matrix with all evaluated values of the parameters and respective profile values
ciprof.q <- function(y, par, return.period, max.ll, step, conf, gg, d){

  # Limites inferior e superior
  lower <- c(-Inf, 1e-6, if (is.null(prior.info)) -Inf else -0.5 + 1e-6)
  upper <- c(Inf, Inf, if (is.null(prior.info)) Inf else 0.5 - 1e-6)

  # Setup an lapply loop to go through return periods
  ls.profile <- lapply(return.period, function(rp){
    
    # Find roots
    p <- 1 - 1/rp                 # later: this is being done twice: here and inside roots.proflik.q()
    q0 <- qgev(unname(par), p)    # original quantile for 'p'
    step.q <- step*abs(q0)        # new step relative to quantile value ?
    left <- roots.proflik.q(-step.q, conf, par, max.ll, rp, y, prior.info, lower, upper)
    right <- roots.proflik.q(step.q, conf, par, max.ll, rp, y, prior.info, lower, upper)
    proflik <- rbind(left, right)
    
    q.name <- paste0("q_", rp)
    ci.lower <- min(proflik[,q.name]) # lower CI limit
    ci.upper <- max(proflik[,q.name]) # upper CI limit

    if(nrow(proflik) <= 20){
      message <- 1 # possible error, few profile values
    } else{
      message <- 0 # might be ok, idk, check 'kappa' values
    }

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
        lim = max.ll - qchisq(conf, 1)/2,
        prior = deparse(prior.info),
        n_profile = nrow(proflik),
        flag_shape = if(abs(par[3]) > 0.5) 1 else 0, # flag '1' whenever |kappa| is larger than 0.5
        # algo = algo,
        conf = conf,
        step = step
      ),
      # prof.par = unlist(proflik[which.max(proflik[,"loglik"]), seq_len(npar + 1)]),
      profile = cbind(
        gauge_code = gg, 
        return_period = rp, as.data.frame(profile)
      )
    )

    return(out)
    
  }) # fim 'ls.profile'

  results <- bind_rows(lapply(ls.profile, `[[`, "ci")) # bind 'ci' rows for each 'rp'
  # prof.par <- bind_rows(lapply(ls.profile, `[[`, "prof.par")) # bind parameter rows
  profile <- bind_rows(lapply(ls.profile, `[[`, "profile")) # bind 'profile' matrices

  list(ci = results, profile = profile)

}

# MAIN FUNCTION ----------------------------------------------------------

fitprof.gev <- function(
  data, prior.info = NULL,
  return.period,
  cols = c("gauge_code", "year", "duration", "imax"),
  conf = 0.95, lltol = 0.005, step = 0.005
){

  #### CHECAGENS ----

  if(!is.list(data) || length(data) == 0){
    stop("`data` must be a non-empty list with data.frames.")
    warning("`cols` must contain column names with 'gauge_code', 'year', 'duration', and 'imax' information in this order.")
  }
  if(!is.numeric(return.period) || length(return.period) == 0) stop("`return.period` must have one or more numeric values.")
  if(!length(prior.info) %in% c(0, 2)) stop("`prior.info` dever ser 'NULL' ou um vetor de comprimento 2.")
  if(!is.numeric(step) || length(step) != 1 | step <= 0) stop("Profile roots search `step` must be a positive number.")
  if(!is.numeric(conf) || length(conf) != 1 || conf <= 0 || conf >= 1) stop("Confidence level `conf` must be between (0,1).")
  if(!is.numeric(lltol) || length(lltol) != 1 || lltol <= 0) stop("Profile likelihood change tolerance `lltol` musc be a positive number.")
  data.cols <- colnames(data[[1]])
  missing.cols <- setdiff(cols, data.cols)
  if(length(missing.cols) > 0) stop(sprintf("Missing columns in `data`: %s", paste(missing.cols, collapse = ", ")))

  #### CÁLCULOS ----

  message("Fitting GEV by durations, estimating quantiles and confidence intervals with 'profile likelihood'.")
  message(">> No. gauges: ", length(data))
  # message(">> Durations [h]: ", paste(durations, collapse = ", "))
  message(">> Return periods [years]: ", paste(return.period, collapse = ", "))
  message(">> Search step: ", step)
  message(">> Confidence level for profile likelihood CIs: ", paste0(conf*100, "%"))

  ls.res <- pbapply::pblapply(data, function(data.gg){

    gg <- unique(data.gg[[cols[1]]])
    y.id <- data.gg %>%
      pivot_wider(names_from = cols[3], values_from = cols[4], names_prefix = "d_") %>%
      arrange(cols[2]) %>%
      select(-c(cols[1], cols[2])) %>%
      as.matrix() %>%
      na.omit()
  
    # Estimate parameters and quantiles with profile likelihood confidence intervals
    # for every duration in 'y.id'
    # ls.prof.d <- apply(X = y.id, 2, function(y){
    ds <- as.integer(gsub("d_", "", colnames(y.id))) # loop through durations instead of y.id columns
    ls.prof.d <- lapply(ds, function(d){
  
      tryCatch({
  
        # Original optimal parameters and maximum
        y <- y.id[,d]
        fit <- fgev(y, prior.info)
        par <- fit$par
        npar <- length(par)
        convergence <- fit$convergence
        max.ll <- -fit$value
      
        # Run profile likelihood
        ls.ci.par <- ciprof.par(
          y = y, 
          par = par, 
          max.ll = max.ll,
          step = step, 
          conf = conf, 
          lltol = lltol,
          gg, d
        )
        opt.par <- ls.ci.par[["ci"]][["par"]]
        opt.ll <- max(ls.ci.par[["ci"]][["found_max"]])
      
        rerun.ci.par <- ciprof.par(
          y = y,
          par = opt.par,
          max.ll = opt.ll,
          step = step,
          conf = conf,
          lltol = lltol,
          gg, d
        )
        
        # Run profile likelihood for quantiles
        opt.par <- rerun.ci.par[["ci"]][["par"]]
        opt.ll <- max(rerun.ci.par[["ci"]][["found_max"]])
        ls.ci.q <- ciprof.q(
          y = y,
          par = opt.par,
          return.period = return.period,
          max.ll = opt.ll, 
          step = step, 
          conf = conf,
          gg, d
        )
      
        return(list(par = rerun.ci.par, quantiles = ls.ci.q))
  
      }, error = function(e){
  
        warning(sprintf("Gauge %s | duration %s failed: %s", gg, d, conditionMessage(e)))
        return(NULL)
  
      })
  
    })
  
    ls.prof.d <- Filter(Negate(is.null), ls.prof.d)
    if(length(ls.prof.d) == 0) stop("No duration converged for this gauge.")
  
    ls.par <- lapply(ls.prof.d, `[[`, "par")             # extract parameters list
    ls.quantiles <- lapply(ls.prof.d, `[[`, "quantiles") # extract quantiles list
    ci.par <- bind_rows(lapply(ls.par, `[[`, "ci"))      # extract paramaters 'ci' data.frame 
    ci.q <- bind_rows(lapply(ls.quantiles, `[[`, "ci"))  # extract quantiles 'ci' data.frame
    prof.par <- bind_rows(lapply(ls.par, `[[`, "profile"))     # extract parameters profile likelihood
    prof.q <- bind_rows(lapply(ls.quantiles, `[[`, "profile")) # extract quantiles profile likelihood
  
    return(list(
      ci.par = ci.par, 
      ci.q = ci.q, 
      prof.par = prof.par, 
      prof.q = prof.q
    ))
  
  })
  
  
  # Bind everything again
  ci.par <- bind_rows(lapply(ls.res, `[[`, "ci.par")); rownames(ci.par) <- NULL
  ci.q <- bind_rows(lapply(ls.res, `[[`, "ci.q"))
  prof.par <- bind_rows(lapply(ls.res, `[[`, "prof.par"))
  prof.q <- bind_rows(lapply(ls.res, `[[`, "prof.q"))

  return(list(
    par = ci.par,
    idf = ci.q,
    prof.par = prof.par,
    prof.q = prof.q
  ))

}
