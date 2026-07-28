# Estimar parâmetro de forma através de diferentes métodos de máxima verossimilhança
# e avaliar os intervalos de confiança com Método Delta

# Métodos
# - GEV padrão: mle
# - GEV padrão: rgmle
# - dGEV: mle
# - dGEV: rgmle

set.seed(1234)
rm(list = ls())
invisible(gc())


# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, arrow, lmom, future, future.apply)

# Set theme
suppressMessages(extrafont::loadfonts())
font <- "Consolas"
fontsize <- 8
textsize <- (fontsize - 2) / .pt
set_theme(
  theme_minimal() +
    theme(
      text = element_text(color = "black", size = fontsize, family = font),
      axis.text = element_text(color = "black"),
      panel.border = element_rect(color = "black")
    )
)


# FUNÇÕES -----------------------------------------------------------------

source("R/funcoes/aux_functions.R") # funções de maximização da GEV convencional
source("R/funcoes/fdgev.R") # funções de maximização da dGEV


# LER DADOS ---------------------------------------------------------------

# Intensidades máximas anuais
df.imax <- arrow::read_parquet("base/gerados/df_imax.pqt")

# Aplicar critérios de qualidade
data <- df.imax %>%
  filter(na_prct_yr <= 0.2, between(d, 1, 24)) %>%
  group_by(gauge_code) %>%
  filter(n_distinct(d) >= 3, n_distinct(wateryear) >= 8) %>%
  ungroup() %>%
  select(gauge_code, wateryear, d, imax)

ls.data <- split(data, data$gauge_code)


# FUNÇÃO -----------------------------------------------------------------

gg <- "A812"
gg <- "430460601A"
data.gg <- data[data$gauge_code == gg, ]

shape.delta <- function(
  data,
  cols = c("gauge_code", "year", "d", "imax"),
  conf = 0.95,
  maxit = 1e6
){
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

  get.shape <- function(data.gg){

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
    source("R/funcoes/fgev.R")

    # Parâmetro de forma GEV + GMLE
    res.gev.gmle <- apply(yd.gg, 2, function(y) {
      start <- flmom(y)
      fit <- try(fgev(start, y, prior.info = c(-0.10, 0.122), maxit = maxit))
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
      fit <- fgev(unname(start), y, prior.info = c(0.0178, 0.0627), maxit = maxit)
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
    fit.dgev.mle <- fdgev(data.gg, cols = names(data.gg), method = "L-BFGS-B", scale.inv = "wide")
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
        ci_lower = kappa.dgev.mle - qnorm(1 - (1 - conf) / 2) * sd.kappa.dgev.mle,
        ci_upper = kappa.dgev.mle + qnorm(1 - (1 - conf) / 2) * sd.kappa.dgev.mle
      )
    } else {
      rep(NA_real_, 3)
    }

    # Parâmetro de forma dGEV + RGMLE
    fit.dgev.rgmle <- fdgev(data.gg, prior.info = c(0.0178, 0.0627), cols = names(data.gg), "L-BFGS-B", "wide")
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
        ci_lower = kappa.dgev.rgmle - qnorm(1 - (1 - conf) / 2) * sd.kappa.dgev.rgmle,
        ci_upper = kappa.dgev.rgmle + qnorm(1 - (1 - conf) / 2) * sd.kappa.dgev.rgmle
      )
    } else {
      rep(NA_real_, 3)
    }

    gg <- unique(data.gg[[cols[1]]])
    durations <- sort(unique(data.gg[[cols[3]]]))
    shape.names <- c(
      "kappa_gev_gmle",
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
      rbind(t(res.gev.gmle), t(res.gev.rgmle), res.dgev.mle, res.dgev.rgmle)
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

df.shape <- shape.delta(ls.data, cols = names(data), conf = 0.95, maxit = 1e8)
arrow::write_parquet(df.shape, "base/gerados/df_shapes_delta.pqt")

df.shape %>%
  rowwise() %>%
  mutate(width = ci_upper - ci_lower) %>%
  filter(width < 1e-4)
