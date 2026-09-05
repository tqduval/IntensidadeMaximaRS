
# AJUSTE DA GEV POR DURAÇÃO ----------------------------------------------

# Estimativas com perfil de verossimilhança por duração (tem que rodar ainda)
# Por duração sem priori
fitprof.gev.nprior.ls <- readRDS("base/gerados/ls_fit_gev_prof_nprior.rds")
fitprof.gev.wprior.ls <- readRDS("base/gerados/ls_fit_gev_prof_wprior.rds")
fitprof.gev.wprior.ninfo.ls <- readRDS("base/gerados/ls_fit_gev_prof_wprior_ninfo.rds")
fitprof.gev.wprior.martins.ls <- readRDS("base/gerados/ls_fit_gev_prof_wprior_martins.rds")

# Parâmetros
par.gev.df <- bind_rows(
  fitprof.gev.nprior.ls$par,
  fitprof.gev.wprior.ls$par,
  fitprof.gev.wprior.ninfo.ls$par,
  fitprof.gev.wprior.martins.ls$par
)

# Quantis
idf.gev.df <- bind_rows(
  fitprof.gev.nprior.ls$idf,
  fitprof.gev.wprior.ls$idf,
  fitprof.gev.wprior.ninfo.ls$idf,
  fitprof.gev.wprior.martins.ls$idf
)

rm(list = ls()[grepl("fitprof.gev", ls())])
invisible(gc())

write_parquet(par.gev.df, "base/gerados/df_par_gev.pqt")
write_parquet(idf.gev.df, "base/gerados/df_idf_gev.pqt")


# AJUSTE d-GEV -----------------------------------------------------------

# Ler lista resultados perfil verossimilhança
# fitprof_dgev.R
fitprof.dgev.nprior.ls <- readRDS("base/gerados/ls_fit_dgev_prof_nprior.rds") # versão sem priori
fitprof.dgev.wprior.ls <- readRDS("base/gerados/ls_fit_dgev_prof_wprior.rds") # versão com priori regional
fitprof.dgev.wprior.ninfo.ls <- readRDS("base/gerados/ls_fit_dgev_prof_wprior_ninfo.rds") # versão com priori pouco informativa
fitprof.dgev.wprior.martins.ls <- readRDS("base/gerados/ls_fit_dgev_prof_wprior_martins.rds") # versão com priori Martins
invisible(gc())

# Checar resultados da rodada 'wprior' com correções
# fitprof.wprior.ls <- readRDS("base/gerados/ls_fit_dgev_prof_wprior2.rds") # temp
# par.df <- fitprof.wprior.ls$par # temp
# idf.df <- fitprof.wprior.ls$idf # temp
# par.df <- fitprof.nprior.ls$par # temp
# idf.df <- fitprof.nprior.ls$idf # temp

# Parâmetros
par.dgev.df <- bind_rows(
  fitprof.dgev.nprior.ls$par,
  fitprof.dgev.wprior.ls$par,
  fitprof.dgev.wprior.ninfo.ls$par,
  fitprof.dgev.wprior.martins.ls$par
)

# Quantis
idf.dgev.df <- bind_rows(
  fitprof.dgev.nprior.ls$idf,
  fitprof.dgev.wprior.ls$idf,
  fitprof.dgev.wprior.ninfo.ls$idf,
  fitprof.dgev.wprior.martins.ls$idf
)

durations <- c(1, 6, 12, 24)

Profiles
prof.par <- bind_rows(
  fitprof.dgev.nprior.ls$prof.par,       # 1
  fitprof.dgev.wprior.ls$prof.par,       # 2
  fitprof.dgev.wprior.ninfo.ls$prof.par, # 3
  .id = "model"
) |> 
  mutate(
    model = ifelse(model == 1, "MLE", "RGMLE")
  ) |> 
  relocate(model, .after = last_col())

prof.q <- bind_rows(
  fitprof.nprior.ls$prof.q, # 1
  fitprof.wprior.ls$prof.q, # 2
  .id = "model"
) |> 
  mutate(
    model = ifelse(model == 1, "MLE", "RGMLE")
  ) |> 
  relocate(model, .after = last_col())

# Remover `fitprof.dgev` do ambiente
rm(list = ls()[grepl("fitprof.dgev", ls())])
invisible(gc())

write_parquet(par.dgev.df, "base/gerados/df_par_dgev.pqt")
write_parquet(idf.dgev.df, "base/gerados/df_idf_dgev.pqt")