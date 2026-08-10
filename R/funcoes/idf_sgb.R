# Função p/ ajustar o modelo IDF conforme metodologia do SGB
# apresentada no Atlas Pluviométrico do Brasil

# O modelo ajustado é o convencionalmente utilizado no Brasil:
# i = k*T^lambda/(d + theta)^eta

# Este script especificamente baseia-se nos ajustes da 
# distribuição Generalizada de Valores Extremos (GEV)
# e requer que o script `fgev.R` seja carregado
path.fgev <- "R/funcoes/fgev.R"
if(file.exists(path.fgev)){
  message("RScript `fgev.R` localizado!")
  source(path.fgev)
} else{
  warning(sprintf("RScript `fgev.R` é necessário e não foi encontrado em `%s`.\nVersões estão disponíveis no repositório remoto `IntensidadeMaximaRS`.",
  dirname(path.fgev)))
}

#' @title Estimar de `theta` por tempo de retorno
#' @detail Encontrar theta que lineariza ln(I) x ln(t + theta)
theta.fo <- function(par, ln.imax, ds){
  ln.ds <- log(ds + par)
  lm.fit <- lm(formula = ln.imax ~ ln.ds)
  res <- sigma(lm.fit)
  return(res)
}

#' @title Estimar parâmetros iniciais do modelo IDF convencional
start.idf <- function(duration.vec, qgev.matrix, return.period){
  
  # Construindo uma matriz de quantis estimados p/ diferentes tempos de retorno
  ln.q.matrix <- log(qgev.matrix)

  ## Estimativas de `ln(A)` e `eta`
  est.by.rp <- apply(ln.q.matrix, 2, function(ln.q, ds){
  
    theta.fit <- optimize(
      f = theta.fo,
      interval = c(-min(ds) + 0.0001, max(ds)),
      # interval = c(0, 100),
      tol = 1e-4,
      ln.imax = ln.q,
      ds = ds,
      maximum = FALSE
    )
  
    theta <- theta.fit$minimum
    ln.ds <- log(ds + theta) # ln das durações + theta
    n.ds <- length(ds)       # nro. durações
    sum.ds <- sum(ln.ds)
    sum.ds2 <- sum(ln.ds^2)
    sum.q <- sum(ln.q)
    sum.q.ds <- sum(ln.q*ln.ds)
  
    aux.A <- sum.q*sum.ds2 - sum.ds*sum.q.ds
    aux.eta <- sum.ds*sum.q - n.ds*sum.q.ds
    term <- n.ds*sum.ds2 - sum.ds^2
    
    log.A <- aux.A/term
    eta <- aux.eta/term
  
    return(c(lnA = log.A, eta = eta, theta = theta))
  
  }, ds = duration.vec)
  
  ln.A <- est.by.rp["lnA",]
  eta <- median(est.by.rp["eta",])
  theta <- median(est.by.rp["theta",])
  
  ## Obter estimativa de `ln(k)` e `lambda`
  n.rp <- length(return.period)      # nro. tempos de retorno 
  ln.rp <- log(return.period)
  sum.rp <- sum(ln.rp)
  sum.A <- sum(ln.A)
  sum.rp2 <- sum(ln.rp^2)
  sum.A.rp <- sum(ln.A*ln.rp)
  aux.k <- sum.rp2*sum.A - sum.A.rp*sum.rp
  aux.lambda <- n.rp*sum.A.rp - sum.rp*sum.A
  term <- n.rp*sum.rp2 - sum.rp^2
  
  ln.k <- aux.k/term        # ln(k)
  k <- exp(1)^ln.k          # k
  lambda <- aux.lambda/term # lambda
  
  return(c(k = k, lambda = lambda, eta = eta, theta = theta))

}

#' @title Calcular intensidades
idf.sgb <- function(par, return.period, duration){
  k <- par[[1]]; lambda <- par[[2]]; eta <- par[[3]]; theta <- par[[4]]
  k*return.period^lambda/(duration + theta)^eta
}

#' @title Função objetivo: calcular RMSE
idf.fo <- function(par, q.base.vec, return.period, duration, metric = c("rmse", "log-rmse", "mape")){

  metric <- match.arg(metric, c("rmse", "log-rmse", "mape"))

  tryCatch({
    q.est.vec <- as.vector(sapply(return.period, function(rp) idf.sgb(par, rp, duration)))
    fo <- switch(
      metric,
      "rmse" = {
        res <- q.base.vec - q.est.vec
        (sum(res^2)/(length(res) - 1))^0.5
      },
      "log-rmse" = {
        offset <- 1e-6
        res <- log(q.base.vec + offset) - log(q.est.vec + offset)
        mean(res^2)^0.5
      },
      "mape" = {
        res <- q.base.vec - q.est.vec
        mean(abs(res)/q.base.vec)*100
      }
    )
  }, error = function(e){
    conditionMessage(e)
    NULL
  })

}

fsgb.idf <- function(par.matrix, duration.vec, return.period){

  inf.var <- any(abs(par.matrix[,"kappa"]) > 0.5)

  # Estimar parâmetros iniciais
  p.nexceed.vec <- 1 - 1/return.period
  q.matrix <- t(apply(par.matrix, 1, qgev, p.nexceed = p.nexceed.vec, simplify = TRUE))
  dim.names <- list(duration.vec, return.period)
  dimnames(q.matrix) <- dim.names
  start <- start.idf(duration.vec, q.matrix, return.period)

  # Otimizar parâmetros iniciais
  fit.obj <- tryCatch(
    optim(
      par = start,
      fn = idf.fo,
      method = "L-BFGS-B",
      control = list(maxit = 1e8, factr = 1e7, pgtol = 1e-8),
      lower = c(0.0001, 0.0001, 0.0001, -min(duration.vec) + 0.0001),
      upper = c(Inf, Inf, 0.9999, min(duration.vec)),
      q.base.vec = as.vector(q.matrix),
      return.period = return.period,
      duration = duration.vec,
      metric = "mape"
    ),
    error = function(e){
      warning(sprintf("Erro durante otimização dos parâmetros iniciais:\n%s", conditionMessage(e)))
      NULL
    }
  )

  return(fit.obj$par)

}

#### TESTE ####

# REFAZER NO EXCEL, AVALIAR PQ ALGUMAS EQUAÇÕES ESTÃO TÃO CURVADAS

gg <- "432260802A"
par.ls <- par.lmom.ls
par.matrix <- par.ls[[gg]]
duration.vec <- seq_len(24)
return.period <- 2^seq_len(8)

par.gg <- start
par.gg <- fit.obj$par |> print()

idf.gg <- freq.obs.all.df |> 
  filter(
    gauge_code == gg,
    plot_pos == "weibull"
  ) |> 
  select(gauge_code, wateryear, d, imax, p, return_period) |> 
  left_join(
    as_tibble(par.matrix) |> mutate(gauge_code = gg, d = seq_len(24), .before = 1),
    by = c("gauge_code", "d")
  ) |> 
  rowwise() |> 
  mutate(
    q_idf = idf.sgb(par = par.gg, return.period = return_period, duration = d),
    q_gev = qgev(param = c(xi, alpha, kappa), p.nexceed = 1 - p)
  ) |> print()

idf.gg |> 
  pivot_wider(
    id_cols = d,
    names_from = return_period,
    values_from = imax
  ) |> 
  as.matrix()

idf.gg |> 
  ggplot(aes(x = d, group = factor(return_period))) +
  # geom_line(aes(y = imax, color = factor(return_period))) +
  # geom_line(aes(y = imax, color = "Intensidades observadas")) +
  geom_line(aes(y = q_gev, color = "GEV")) +
  geom_line(aes(y = q_idf, color = "IDF-SGB")) +
  annotate(
    geom = "text", x = Inf, y = Inf, 
    label = sprintf("k: %0.0f\nλ: %0.3f\nη: %0.3f\nθ: %0.3f", par.gg[1], par.gg[2], par.gg[3], par.gg[4]),
    family = font, hjust = 1.1, vjust = 1.1
  ) +
  # scale_color_manual(name = NULL, values = c("Intensidades observadas" = "grey70", "IDF-SGB" = "black")) +
  scale_color_manual(name = NULL, values = c("GEV" = "grey70", "IDF-SGB" = "black")) +
  scale_x_log10() +
  scale_y_log10()


#### FIM TESTE ####
  
#' @title Ajustar o modelo final de IDF
fitsgb.idf <- function(par.ls, duration.vec = NULL, return.period){

  idf.ls <- pbapply::pblapply(seq_along(par.ls), function(i){

    if(is.null(duration.vec)) duration.vec <- as.numeric(rownames(par.matrix))
    gg <- names(par.ls[i])
    par.matrix <- par.ls[[i]]
    inf.var <- any(abs(par.matrix[,"kappa"]) > 0.5)

    # Estimar parâmetros iniciais
    p.nexceed.vec <- 1 - 1/return.period # probabilidade de não excedência
    q.matrix <- t(apply(par.matrix, 1, qgev, p.nexceed = p.nexceed.vec, simplify = TRUE)) # quantis da GEV
    dim.names <- list(duration.vec, return.period) # nome linhas e colunas
    dimnames(q.matrix) <- dim.names
    start <- start.idf(duration.vec, q.matrix, return.period) # estimar parâmetros iniciais

    # Otimizar parâmetros iniciais
    fit.obj <- tryCatch({

      fit1 <- optim(
        par = start,
        fn = idf.fo,
        # method = "L-BFGS-B",
        method = "Nelder-Mead",
        control = list(maxit = 1e8, factr = 1e7, pgtol = 1e-8),
        # lower = c(0.0001, 0.0001, 0.0001, -min(duration.vec) + 0.0001),
        # upper = c(Inf, Inf, Inf, min(duration.vec)),
        q.base.vec = as.vector(q.matrix),
        return.period = return.period,
        duration = duration.vec,
        metric = "rmse"
      )

      fit2 <- optim(
        par = fit1$par,
        fn = idf.fo,
        # method = "L-BFGS-B",
        method = "Nelder-Mead",
        control = list(maxit = 1e8, factr = 1e7, pgtol = 1e-8),
        # lower = c(0.0001, 0.0001, 0.0001, -min(duration.vec) + 0.0001),
        # upper = c(Inf, Inf, Inf, min(duration.vec)),
        q.base.vec = as.vector(q.matrix),
        return.period = return.period,
        duration = duration.vec,
        metric = "mape"
      )

    },
      error = function(e){
        warning(sprintf("Posto %s: erro durante otimização dos parâmetros iniciais :\n%s", gg, conditionMessage(e)))
        NULL
      }
    )

    # Extrair parâmetros finais
    if(!is.null(fit.obj)){
      opt.par <- fit.obj$par
      idf.par <- data.frame(
        gauge_code = gg, 
        t(opt.par), 
        optim_fail = 0, 
        inf_var = inf.var
      )
    } else{
      idf.par <- data.frame(
        gauge_code = gg, 
        k = NA_real_,
        lambda = NA_real_,
        eta = NA_real_,
        theta = NA_real_, 
        optim_fail = 1, 
        inf_var = inf.var
      )
    }

    return(par = idf.par)

  })

  return(bind_rows(idf.ls))

}

#' @title Ajustar o modelo final de IDF
buildsgb.idf <- function(par.ls, duration.vec = NULL, return.period){

  idf.ls <- pbapply::pblapply(seq_along(par.ls), function(i){

    gg <- names(par.ls[i])
    par.matrix <- par.ls[[i]]
    inf.var <- any(abs(par.matrix[,"kappa"]) > 0.5)
    durations.gg <- as.numeric(rownames(par.matrix))

    if(is.null(duration.vec)){
      duration.vec <- durations.gg
    }

    if(isFALSE(any(durations.gg %in% duration.vec))){
      stop(sprintf("Posto %s:\nAlgumas durações em `duration.vec` não estão em `par.ls[['%s']]$d`", gg, gg))
    }

    if(any(!duration.vec %in% durations.gg)){
      warning(sprintf("Posto %s:\nAlgumas durações em `duration.vec` não estão em `par.ls[['%s']]$d` e não serão utilizadas.", gg, gg))
      duration.vec <- duration.vec[duration.vec %in% durations.gg]
    }

    # Estimar parâmetros iniciais
    p.nexceed.vec <- 1 - 1/return.period
    q.matrix <- t(apply(par.matrix, 1, qgev, p.nexceed = p.nexceed.vec, simplify = TRUE))
    dim.names <- list(duration.vec, return.period)
    dimnames(q.matrix) <- dim.names
    start <- start.idf(duration.vec, q.matrix, return.period)

    # Otimizar parâmetros iniciais
    fit.obj <- tryCatch(
      optim(
        par = start,
        fn = idf.fo,
        method = "L-BFGS-B",
        control = list(maxit = 1e8, factr = 1e7, pgtol = 1e-8),
        lower = c(0.0001, 0.0001, 0.0001, -min(duration.vec) + 0.0001),
        upper = c(Inf, Inf, Inf, Inf),
        q.base.vec = as.vector(q.matrix),
        return.period = return.period,
        duration = duration.vec
      ),
      error = function(e){
        warning(sprintf("Posto %s: erro durante otimização dos parâmetros iniciais:\n", conditionMessage(e)))
        NULL
      }
    )

    if(!is.null(fit.obj)){

      opt.par <- fit.obj$par
      quantiles <- sapply(return.period, idf.sgb, par = opt.par, duration = duration.vec)
      dimnames(quantiles) <- dim.names

      quantiles <- as.data.frame(quantiles) |> 
        mutate(gauge_code = gg, .before = 1) |> 
        pivot_longer(
          cols = !gauge_code,
          names_to = "return_period",
          names_transform = as.numeric,
          values_to = "q"
        )
  
      idf.par <- tibble(gauge_code = gg, t(opt.par), optim_fail = 0, inf_var = inf.var)
      quantiles <- quantiles

    } else{

      idf.par <- tibble(gauge_code = gg, t(rep(NA_real_, 3)), optim_fail = 1, inf_var = inf.var)
      quantiles <- NULL

    }

    return(list(
      par = idf.par,
      idf = quantiles
    ))

  })

  par <- bind_rows(lapply(idf.ls, `[[`, "par"))
  idf <- bind_rows(lapply(idf.ls, `[[`, "idf"))

  return(list(par = par, idf = idf))

}


# AJUSTE REGRESSÃO -------------------------------------------------------

rm(list = ls())

rp <- c(2, 5, 10, 25, 50, 100, 200)
p.nexceed.vec <- 1 - 1/rp
d <- seq_len(24)

# Estimativas quantis GEV p/ montar matriz Y = log(qgev)
par.lmom.ls <- split(imax.df, imax.df$gauge_code) |> 
  lapply(function(data.gg){
    data.gg |> 
      select(wateryear, d, imax) |> 
      pivot_wider(names_from = d, values_from = imax) |> 
      column_to_rownames(var = "wateryear") |> 
      as.matrix() |> na.omit() |> 
      apply(2, start.gev) |> 
      t()
  })

fit.theta <- function(theta, df){
  df <- df |> mutate(log_d = log(d + theta))
  fit <- lm(log_y ~ log_rp + log_d, df)
  return(sigma(fit))
}

setup.fit.data <- function(q.matrix, duration.vec, theta){
  as.data.frame(q.matrix) |> 
    mutate(d = duration.vec) |>
    pivot_longer(cols = -d, names_to = "rp", names_transform = as.numeric, values_to = "y") |>
    mutate(log_y = log(y), log_rp = log(rp), log_d = log(d + theta))
}

gg <- "430060402A"
data.gg <- imax.df[imax.df$gauge_code == gg,]
par.gg <- par.lmom.ls[[gg]]

Y <- t(apply(par.gg, 1, qgev, p.nexceed = p.nexceed.vec, simplify = TRUE))
dimnames(Y) <- list(d, rp)

# Ajuste inicial
theta.start <- 0
fit.data.start <- setup.fit.data(Y, d, theta.start) |> print()
fit.lm <- lm(log_y ~ log_rp + log_d, data = fit.data.start)
summary.fit <- summary(fit.lm) |> print()
coef <- fit.lm$coefficients
k <- exp(coef[["(Intercept)"]])
gamma <- coef[["log_rp"]]
eta <- -coef[["log_d"]]
message(sprintf("start: k: %6.3f | gamma: %6.3f | eta: %6.3f | theta: %6.3f", k, gamma, eta, theta.start))

# Ajuste final
theta.opt <- optimize(f = fit.theta, interval = c(0, max(d)/3), df = fit.data.start)
theta.opt <- theta.opt$minimum
fit.data.opt <- setup.fit.data(Y, d, theta.opt) |> print()
fit.lm <- lm(log_y ~ log_rp + log_d, data = fit.data.opt)
summary.fit <- summary(fit.lm) |> print()
coef <- fit.lm$coefficients
k <- exp(coef[["(Intercept)"]])
gamma <- coef[["log_rp"]]
eta <- -coef[["log_d"]]
message(sprintf("start: k: %6.3f | gamma: %6.3f | eta: %6.3f | theta: %6.3f", k, gamma, eta, theta.opt))

est.data <-  expand.grid(d = d, rp = rp) |> 
  select(d, rp) |> 
  mutate(y = k*rp^gamma/(d + theta.opt)^eta) |> 
  print()

ggplot() +
  geom_line(aes(x = d, y = y, group = rp, color = "GEV"), fit.data.opt) +
  geom_line(aes(x = d, y = y, group = rp, color = "IDF convencional OLS"), est.data) +
  scale_x_log10("Duração [h]") +
  scale_y_log10("Intensidade [mm/h]") +
  scale_color_manual(name = NULL, values = c("grey", "black")) +
  theme(
    legend.position = c(0.98, 0.98),
    legend.justification = c(1,1)
  )

  #' @title Ajustar o modelo final de IDF
  fitsgb.idf <- function(par.ls, duration.vec = NULL, return.period){

    fit.theta <- function(theta, df){
      df <- df |> mutate(log_d = log(d + theta))
      fit <- lm(log_y ~ log_rp + log_d, df)
      return(sigma(fit))
    }

    setup.fit.data <- function(q.matrix, duration.vec, theta){
      as.data.frame(q.matrix) |> 
        mutate(d = duration.vec) |>
        pivot_longer(cols = -d, names_to = "rp", names_transform = as.numeric, values_to = "y") |>
        mutate(log_y = log(y), log_rp = log(rp), log_d = log(d + theta))
    }
  
    idf.ls <- pbapply::pblapply(seq_along(par.ls), function(i){
  
      if(is.null(duration.vec)) duration.vec <- as.numeric(rownames(par.matrix))
      gg <- names(par.ls[i])
      par.matrix <- par.ls[[i]]
      inf.var <- any(abs(par.matrix[,"kappa"]) > 0.5)
  
      # Estimar parâmetros iniciais
      p.nexceed.vec <- 1 - 1/return.period # probabilidade de não excedência
      q.matrix <- t(apply(par.matrix, 1, qgev, p.nexceed = p.nexceed.vec, simplify = TRUE)) # quantis da GEV
      dim.names <- list(duration.vec, return.period) # nome linhas e colunas
      dimnames(q.matrix) <- dim.names
      
      fit.par <- tryCatch({
      
        # Ajuste inicial
        theta.start <- 0
        fit.data.start <- setup.fit.data(q.matrix, d, theta.start) 
        
        # Ajuste final
        theta.opt <- optimize(f = fit.theta, interval = c(0, max(d)/3), df = fit.data.start)
        theta.opt <- theta.opt$minimum
        fit.data.opt <- setup.fit.data(q.matrix, d, theta.opt)
        fit.lm <- lm(log_y ~ log_rp + log_d, data = fit.data.opt)
        coef <- coef(fit.lm)
        k <- exp(coef[[1]]) # intercepto
        gamma <- coef[[2]]  # log(rp) 
        eta <- -coef[[3]]   # log(d + theta)
  
        par <- setNames(c(k, gamma, eta, theta.opt), c("k", "gamma", "eta", "theta"))
      
      }, error = function(e){
        
        warning(conditionMessage(e))
        NULL

      })
  
      # Extrair parâmetros finais
      if(!is.null(fit.par)){
        idf.par <- data.frame(
          gauge_code = gg, 
          t(fit.par), 
          fit_fail = 0, 
          inf_var = inf.var
        )
      } else{
        idf.par <- data.frame(
          gauge_code = gg, 
          k = NA_real_,
          gamma = NA_real_,
          eta = NA_real_,
          theta = NA_real_, 
          fit_fail = 1, 
          inf_var = inf.var
        )
      }
  
      return(par = idf.par)
  
    })
  
    return(bind_rows(idf.ls))
  
  }
