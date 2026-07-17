rm(list = ls()); gc()
pacman::p_load(pacman, tidyverse)


imax.df <- arrow::read_parquet("base/gerados/df_imax.pqt")
imax.df <- imax.df |> 
  filter(na_prct_yr <= 0.2,
         between(d, 1, 24)) |> 
  group_by(gauge_code) |>
  filter(n_distinct(d) >= 3,
         n_distinct(wateryear) >= 8) |> 
  ungroup() |> 
  select(gauge_code, wateryear, d, imax)

# fitprof.wprior.ls <- readRDS("base/gerados/ls_fit_dgev_prof_wprior.rds")
fitprof.wprior.ls <- readRDS("base/gerados/ls_fit_dgev_prof_wprior2.rds")
# fitprof.wprior.ls <- fitprof.dgev.wprior.ls
par.df <- fitprof.wprior.ls$par
fitprof.nprior.ls <- fitprof.dgev.nprior.ls
par.df <- fitprof.nprior.ls$par

gg <- unique(par.df$gauge_code)
gg <- "A801"

par.gg <- par.df |> 
  filter(gauge_code == gg) |> 
  select(!c(gauge_code, n_profile, conf, step, rel_width, ci_lower, ci_upper, ci_type, prior)) |> 
  print()

profile.gg <- fitprof.nprior.ls$prof.par |> 
  filter(gauge_code == gg)

par.plot <- "kappa"

profile.plot <- profile.gg |> 
  filter(parname == par.plot) |> 
  select(par, loglik, slope)

max(profile.plot$slope)

plot(profile.plot[,-3], main = gg, sub = par.plot)

idf.df <- fitprof.nprior.ls$idf
profile.idf.gg <- fitprof.nprior.ls$prof.q |> 
  filter(d == 1, return_period == 50) |> 
  select(q, loglik)

profile.idf.gg <- profile.idf.gg[order(profile.idf.gg$q),]

plot(profile.idf.gg, main = gg)

data.gg <- imax.df[imax.df$gauge_code == gg,]
cols <- names(data.gg)

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
  prior.info = c(0.0178, 0.0627),
  cols = cols,
  method = "L-BFGS-B",
  scale.inv = "wide"
) |> print()

step <- 0.005
step.adjusted <- step*abs(fit$par[[3]])
step.adjusted
val.i <- seq(fit$par[[3]], -1, -step.adjusted)
idx <- which.min(abs(val.i))
val.i[idx]
val.i[c(idx-1, idx, idx+1)]
