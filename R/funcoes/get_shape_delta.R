shape.delta <- function(
    data,
    cols = c("gauge_code", "year", "d", "imax"),
    conf = 0.95,
    maxit = 1e6
) {
  if (any(names(data[[1]]) != cols)) {
    stop(
      "Nomes das colunas e 'cols'não corresponde aos nomes das colunas em 'data'\n",
      "names(data[[1]]): ",
      paste(names(data[[1]]), collapse = ","),
      "\n",
      "cols: ",
      paste(cols, collapse = ", ")
    )
  }
  
  get.shape <- function(data.gg) {
    yd.gg <- data.gg %>%
      pivot_wider(names_from = d, values_from = imax, names_prefix = "d_") %>%
      arrange(wateryear) %>%
      select(-c(gauge_code, wateryear)) %>%
      as.matrix() %>%
      na.omit()
    
    n.years <- nrow(yd.gg)
    
    # Funções
    # Calcular parâmetros iniciais com momentos-L
    flmom <- function(ams) {
      sample.lmom <- lmom::samlmu(ams)
      lmom::pelgev(sample.lmom)
    }
    
    # Ajustar GEV convencional sem dependência temporal
    fgev <- function(param0, y, prior.info = NULL, maxit) {
      nll.gev <- function(param0, y, prior.info) {
        xi <- param0[1]
        alpha <- param0[2]
        kappa <- param0[3]
        n <- length(y)
        
        if (!is.null(prior.info)) {
          mu.kappa <- prior.info[1]
          var.kappa <- prior.info[2]^2
          
          a <- 0.5
          b <- mu.kappa + a
          p <- b^2 * (1 - b) / var.kappa - b
          q <- p * (1 / b - 1)
          
          if (a - abs(kappa) <= 1e-6) {
            return(1e6 + abs(kappa) * 1e4)
          } # impedir log(<0)
          # if(a - abs(kappa) <= 1e-6) return(1e6) # impedir log(<0)
          l.prior <- (p - 1) * log(a + kappa) + (q - 1) * log(a - kappa)
        } else {
          l.prior <- 0
        }
        
        zi <- 1 - kappa / alpha * (y - xi)
        
        # Fallbacks
        if (alpha <= 0) {
          return(1e6 + abs(alpha) * 1e4)
        }
        if (any(zi <= 1e-8)) {
          return(1e6 + abs(min(zi)) * 1e4)
        }
        
        # Negative log-likelihood
        term1 <- (1 / kappa - 1) * log(zi)
        term2 <- zi^(1 / kappa)
        nll <- -(-n * log(alpha) + sum(term1 - term2) + l.prior)
        
        return(nll)
      }
      
      # Função gradiente nll.gev()
      grad.ll.gev <- function(param, y, prior.info) {
        xi <- param[1]
        alpha <- param[2]
        kappa <- param[3]
        N <- length(y)
        
        if (alpha <= 1e-8) {
          return(c(0, 1e6, 0)) # corrigir alpha (param[3]) e não os demais
        }
        
        # Avaliar se há distribuição a priori
        if (!is.null(prior.info)) {
          mu.kappa <- prior.info[1]
          var.kappa <- prior.info[2]^2
          
          a <- 0.5
          b <- mu.kappa + a
          p <- b^2 * (1 - b) / var.kappa - b
          q <- p * (1 / b - 1)
          
          if (a - abs(kappa) <= 1e-6) {
            return(rep(1e6, 3))
          }
          
          dl.prior <- (p - 1) / (a + kappa) - (q - 1) / (a - kappa)
        } else {
          dl.prior <- 0
        }
        
        # Termos simplificadores
        zi <- 1 - kappa / alpha * (y - xi)
        if (any(zi <= 1e-8)) {
          return(rep(1e6, 3))
        }
        inv.z <- 1 / zi
        z.inv.k <- zi^(1 / kappa)
        z.inv.km1 <- zi^(1 / kappa - 1)
        ln.z <- log(zi)
        
        # Derivadas em relação a zi
        dz.dxi <- kappa / alpha
        dz.dalpha <- kappa / alpha^2 * (y - xi)
        dz.dkappa <- -(y - xi) / alpha
        
        # Gradientes
        grad.xi <- -sum(
          (1 / kappa - 1) * inv.z * dz.dxi - 1 / kappa * z.inv.km1 * dz.dxi
        )
        grad.alpha <- N /
          alpha -
          sum(
            (1 / kappa - 1) *
              inv.z *
              dz.dalpha -
              1 / kappa * z.inv.km1 * dz.dalpha
          )
        grad.kappa <- -sum(
          (1 / kappa - 1) *
            inv.z *
            dz.dkappa -
            ln.z / kappa^2 -
            z.inv.k * (1 / kappa * inv.z * dz.dkappa - ln.z / kappa^2)
        ) -
          dl.prior
        
        return(c(grad.xi, grad.alpha, grad.kappa))
      } # fim 'grad.ll.gev()'
      
      # Definir limites para L-BFGS-B
      lower <- c(
        -Inf,
        1e-6,
        ifelse(is.null(prior.info), -Inf, -0.5 + 1e-6),
        1e-6,
        0
      )
      upper <- c(
        Inf,
        Inf,
        ifelse(is.null(prior.info), -Inf, 0.5 - 1e-6),
        1 - 1e-6,
        Inf
      )
      
      # Usando Nelder-Mead
      fit <- tryCatch(
        {
          fit.obj <- optim(
            par = param0,
            fn = nll.gev,
            method = "L-BFGS-B",
            gr = grad.ll.gev,
            lower = lower,
            upper = upper,
            hessian = TRUE,
            y = y,
            prior.info = prior.info,
            control = list(
              maxit = maxit,
              factr = 1e10, # relaxa a tolerância da função
              pgtol = 1e-4, # adiciona tolerância no gradiente
              lmm = 20, # aumenta a memória da hessiana aproximada
              # trace = 1     # acompanhar o progresso do ajuste
            )
          )
        },
        error = function(e) {
          warning("Erro ao maximizar a minimizar 'nll.gev(): ", e)
          return(NULL)
        }
      )
      
      return(fit)
    }
    
    # Parâmetro de forma GEV + MLE
    res.gev.mle <- apply(yd.gg, 2, function(y) {
      start <- flmom(y)
      fit <- try(fgev(start, y, maxit = maxit))
      par <- fit$par
      hes <- fit$hessian
      cov.par <- try(solve(hes), silent = TRUE)
      if (inherits(cov.par, "try-error")) {
        return(matrix(NA, length(par), length(par)))
      }
      sd.kappa <- sqrt(diag(cov.par))[[3]]
      kappa <- par[[3]]
      
      out <- c(
        kappa = kappa,
        ci_lower = kappa - qnorm(1 - (1 - conf) / 2) * sd.kappa,
        ci_upper = kappa + qnorm(1 - (1 - conf) / 2) * sd.kappa
      )
      
      return(out)
    })
    
    # Parâmetro de forma GEV + RGMLE
    res.gev.rgmle <- apply(yd.gg, 2, function(y) {
      start <- flmom(y)
      if (abs(start["k"]) >= 0.5) {
        start["k"] <- 0.0178
      }
      fit <- fgev(
        unname(start),
        y,
        prior.info = c(0.0178, 0.0627),
        maxit = maxit
      )
      par <- fit$par
      hes <- fit$hessian
      cov.par <- try(solve(hes), silent = TRUE)
      if (inherits(cov.par, "try-error")) {
        return(matrix(NA, length(par), length(par)))
      }
      sd.kappa <- sqrt(diag(cov.par))[[3]]
      kappa <- par[[3]]
      
      out <- c(
        kappa = kappa,
        ci_lower = kappa - qnorm(1 - (1 - conf) / 2) * sd.kappa,
        ci_upper = kappa + qnorm(1 - (1 - conf) / 2) * sd.kappa
      )
      
      return(out)
    })
    
    # Parâmetro de forma dGEV + MLE
    fit.dgev.mle <- fdgev(
      data.gg,
      cols = names(data.gg),
      method = "L-BFGS-B",
      scale.inv = "wide"
    )
    if (!is.null(fit.dgev.mle) && fit.dgev.mle$convergence == 0) {
      kappa.dgev.mle <- fit.dgev.mle$par[[3]]
      hes.dgev.mle <- fit.dgev.mle$hes
      cov.par.dgev.mle <- try(solve(hes.dgev.mle), silent = TRUE)
      if (inherits(cov.par.dgev.mle, "try-error")) {
        return(matrix(NA, 5, 5))
      }
      sd.kappa.dgev.mle <- sqrt(diag(cov.par.dgev.mle))[[3]]
      
      res.dgev.mle <- c(
        kappa = kappa.dgev.mle,
        ci_lower = kappa.dgev.mle -
          qnorm(1 - (1 - conf) / 2) * sd.kappa.dgev.mle,
        ci_upper = kappa.dgev.mle +
          qnorm(1 - (1 - conf) / 2) * sd.kappa.dgev.mle
      )
    } else {
      NA
    }
    
    # Parâmetro de forma dGEV + RGMLE
    fit.dgev.rgmle <- fdgev(
      data.gg,
      prior.info = c(0.0178, 0.0627),
      cols = names(data.gg),
      "L-BFGS-B",
      "wide"
    )
    if (!is.null(fit.dgev.rgmle) && fit.dgev.rgmle$convergence == 0) {
      kappa.dgev.rgmle <- fit.dgev.rgmle$par[[3]]
      hes.dgev.rgmle <- fit.dgev.rgmle$hes
      cov.par.dgev.rgmle <- try(solve(hes.dgev.rgmle), silent = TRUE)
      if (inherits(cov.par.dgev.rgmle, "try-error")) {
        return(matrix(NA, 5, 5))
      }
      sd.kappa.dgev.rgmle <- sqrt(diag(cov.par.dgev.rgmle))[[3]]
      
      res.dgev.rgmle <- c(
        kappa = kappa.dgev.rgmle,
        ci_lower = kappa.dgev.rgmle -
          qnorm(1 - (1 - conf) / 2) * sd.kappa.dgev.rgmle,
        ci_upper = kappa.dgev.rgmle +
          qnorm(1 - (1 - conf) / 2) * sd.kappa.dgev.rgmle
      )
    } else {
      NA
    }
    
    gg <- unique(data.gg[[cols[1]]])
    durations <- sort(unique(data.gg[[cols[3]]]))
    shape.names <- c(
      "kappa_gev_mle",
      "kappa_gev_rgmle",
      "kappa_dgev_mle",
      "kappa_dgev_rgmle"
    )
    out.names <- c(
      paste0(rep(shape.names[1:2], each = length(durations)), "_d_", durations),
      shape.names[3:4]
    )
    out <- bind_cols(
      id = out.names,
      rbind(t(res.gev.mle), t(res.gev.rgmle), res.dgev.mle, res.dgev.rgmle)
    )
    
    return(out)
  }
  
  ls.res <- pbapply::pblapply(data, function(data.gg) {
    out <- tryCatch(
      {
        gg <- unique(data.gg[[cols[1]]])
        res <- get.shape(data.gg)
        
        out <- bind_cols(gauge_code = gg, res) %>%
          extract(
            col = id,
            into = c("model", "method", "d"),
            regex = "kappa_([a-z]+)_([a-z]+)(?:_d_([0-9]+))?$",
            remove = FALSE,
            convert = TRUE
          )
        
        return(out)
      },
      error = function(e) {
        warning("Erro na estação ", gg, ".\n", e)
      }
    )
  })
  
  return(bind_rows(ls.res))
}