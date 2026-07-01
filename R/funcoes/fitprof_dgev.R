
# INTERNAL FUNCTIONS -----------------------------------------------------

# Fit d-GEV distribution with maximum likelihood or generalized maximum likelihood
source("R/funcoes/fdgev.R")

# d-GEV negative log-likelihood function for all but one parameter
prof.nll.dgev <- function(par.free, which.par, par.value, y.id, durations, prior.info){

  full <- numeric(length(par.free) + 1) # build 'new' parameters vector
  full[which.par] <- par.value          # fixed value for current parameter
  full[-which.par] <- par.free          # remaining parameters to be optimized

  nll.dgev(full, durations, y.id, prior.info)

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
roots.proflik <- function(y.id, step, conf, par, i, ll, lower, upper){

  val.i <- par[[i]]
  npar <- length(par)
  lim <- ll - qchisq(conf, 1)/2

}

# MAIN FUNCTION ----------------------------------------------------------



fitdelta.dgev <- function(
  data, prior.info = NULL,
  return.period,
  cols = c("gauge_code", "year", "duration", "imax"),
  conf = 0.95, lltol = 0.005, step = 0.005
){

}


# TESTING THE FUNCTION ---------------------------------------------------

rm(list = ls()); invisible(gc())
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

## INSIDE THE FUNCTION ---------------------------------------------------

# Everything below this is inside a 'pblapply()' loop
data.gg <- imax.ls[[1]]
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


# Definir limites para L-BFGS-B
lower <- c(-Inf, 1e-6, ifelse(is.null(prior.info), -Inf, -0.5 + 1e-6), 1e-6, 0)
upper <- c(Inf, Inf, ifelse(is.null(prior.info), Inf, 0.5 - 1e-6), 1 - 1e-6, Inf)

# Test 'roots.proflik()' function
i <- 1 # temp
ll <- max.ll # temp
val.i <- par[[i]]
last.par <- par
npar <- length(par)
lim <- ll - qchisq(conf, 1)/2

prof.nll.dgev(last.par[-i], i, val.i, y.id, ds, prior.info) # funcionando (retorna) o mesmo que 'max.ll'

fit.prof <- tryCatch(
  optim(
    par = last.par[-i],
    fn = function(internal.par){
      prof.nll.dgev(
        par.free = internal.par, which.par = i, par.value = val.i,
        y.id = y.id, durations = ds, prior.info = prior.info
      )
    },
    gr = NULL, method = "L-BFGS-B", lower = lower[-i], upper = upper[-i],
    control = list(maxit = maxit, factr = 1e9, pgtol = 1e-8)
  ), error = function(e) NULL
)

if(!is.null(fit.prof) && fit.prof$convergence == 0){

  ll <- -fit$value # save new profile value 'll'

} else{

  fit.prof <- tryCatch(
    optim(
      par = last.par[-i],
      fn = function(internal.par){
        prof.nll.dgev(
          par.free = internal.par, which.par = i, par.value = val.i,
          y.id = y.id, durations = ds, prior.info = prior.info
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
full <- val.i

# Check GEV limits
xit <- full[1]
alpha0 <- full[2]
kappa <- full[3]
H <- full[4]
theta <- full[5]

scale.fact <- (ds + theta)^H
