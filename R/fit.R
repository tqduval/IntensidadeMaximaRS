
# PACOTES -----------------------------------------------------------------

# Limpar ambiente
rm(list = ls()); invisible(gc())
# set.seed(1234)

# Pacotes
if(!require(pacman)) install.packages("pacman")
pacman::p_load(pacman,    # gerenciador de pacotes
               arrow,     # leitura de dados 'parquet'
               tidyverse)

# Set theme
suppressMessages(extrafont::loadfonts())
font <- "JetBrainsMono NF"
fontsize <- 8
textsize <- (fontsize - 2)/.pt
set_theme(theme_minimal() +
            theme(text = element_text(color = "black", size = fontsize, family = font),
                  axis.text = element_text(color = "black"),
                  panel.border = element_rect(color = "black")))


# FUNÇÕES -----------------------------------------------------------------

# Ajustar IDF d-GEV e GEV, estimar quantis e intervalos de confiança
source("R/funcoes/fitprof_dgev.R")
source("R/funcoes/fitprof_gev.R")


# LER DADOS ---------------------------------------------------------------

# Intensidades máximas anuais de precipitação
imax.df <- arrow::read_parquet("base/gerados/df_imax.pqt")

# Aplicar critérios de qualidade
na.yr <- 0.2
min.durations <- 3
d.interval <- c(1, 24)
min.yr <- 8

data.df <- imax.df %>% 
  filter(na_prct_yr <= na.yr,
         between(d, d.interval[1], d.interval[2])) %>% 
  group_by(gauge_code) %>% 
  filter(n_distinct(d) >= min.durations,
         n_distinct(wateryear) >= min.yr) %>% 
  ungroup() %>% 
  select(gauge_code, wateryear, d, imax)

# n_distinct(data$gauge_code) # must return 54

data.ls <- split(data.df, data.df$gauge_code)


# CÁLCULOS COM PERFIL DE VEROSSIMILHANÇA ---------------------------------

# Define function arguments
which.prior <- list(NULL, c(0.0178, 0.0627))
cols <- names(data.df)
method <- "L-BFGS-B"
scale.inv <- "wide"
# durations <- c(1, 2, 4, 6, 8, 10, 12, 18, 24)
durations <- NULL  # roda com todas as durações disponíveis
return.period <- c(1.1, 1.5, 2, 5, 10, 20, 50, 100, 200, 500)
conf <- 0.95
maxit <- 1e8
step <- 0.005
lltol <- 0.005

# Rodar função sem priori
prior.info <- which.prior[[1]]
fitprof.dgev.nprior.ls <- fitprof.dgev(
  data = data.ls,
  prior.info = prior.info,
  return.period = return.period,
  durations = durations,
  cols = cols
)

saveRDS(fitprof.dgev.nprior.ls, "base/gerados/ls_fit_dgev_prof_nprior.rds")

# Rodar função com priori
prior.info <- which.prior[[2]]
fitprof.dgev.wprior.ls <- fitprof.dgev(
  data = data.ls,
  prior.info = prior.info,
  return.period = return.period,
  cols = cols
)

invisible(gc())

saveRDS(fitprof.dgev.wprior.ls, "base/gerados/ls_fit_dgev_prof_wprior.rds")


# # CÁLCULOS COM DELTA E BOOTSTRAP -----------------------------------------

# # Argumentos
# prior.info <- c(0.0178, 0.0627)
# cols <- names(data)
# durations <- c(1, 2, 4, 6, 8, 10, 12, 18, 24)
# return.period <- c(1.1, 1.5, 2, 5, 10, 20, 50, 100, 200, 500)
# method <- "L-BFGS-B"
# scale.inv <- "wide"
# B <- 10000


# ## ESTIMATIVA COM INVARIÂNCIA DE ESCALA ----------------------------------

# # Possibilidades
# prior <- list(NULL, prior.info)
# ci <- c("delta", "bootstrap")
# # combinacoes <- expand.grid(prior = prior, ci_type = ci)
# combinacoes <- list(
#   list(prior = NULL, ci_type = "delta"),
#   list(prior = NULL, ci_type = "bootstrap"),
#   list(prior = prior.info, ci_type = "delta"),
#   list(prior = prior.info, ci_type = "bootstrap")
# )

# ls.fit.d <- purrr::map(combinacoes, function(comb){
#   fit <- fit.dgev(data = ls.data,
#                   prior.info = comb$prior,
#                   durations = durations,
#                   return.period = return.period,
#                   cols = cols,
#                   method = method,
#                   scale.inv = scale.inv,
#                   which.ci = comb$ci_type,
#                   B = B,
#                   maxit = 1e8)
# })

# df.par <- bind_rows(lapply(ls.fit.d, `[[`, "parameters"))
# df.idf <- bind_rows(lapply(ls.fit.d, `[[`, "quantiles"))
# arrow::write_parquet(df.par, "base/gerados/df_par.pqt")
# arrow::write_parquet(df.idf, "base/gerados/df_idf.pqt")


# ## ESTIMATIVAS POR DURAÇÃO -----------------------------------------------

# # Possíveis alternativas de restrição no parâmetro de forma
# ls.prior <- list(NULL, c(-0.10, 0.122), c(0.0178, 0.0627))
# ls.fit <- lapply(ls.prior, function(prior){
#   fit <- fitboot.gev(
#     data = ls.data,
#     prior.info = prior,
#     durations = durations,
#     return.period = return.period,
#     cols = names(data),
#     B = B,
#     maxit = 1e8
#   )
# })

# gc()

# saveRDS(ls.fit, "base/gerados/resultados_por_duracao/ls_fit_by_d.rds")

# # Read RDS file and bind results together
# ls.fit <- readRDS("base/gerados/resultados_por_duracao/ls_fit_by_d.rds")
# ls.par <- vector("list", length(ls.fit)) # lista vazia p/ parâmetros
# ls.idf <- vector("list", length(ls.fit)) # lista vazia p/ quantis

# for(i in seq_along(ls.fit)){
#   ls.par[[i]] <- bind_rows(lapply(ls.fit[[i]], `[[`, "par")) 
#   ls.idf[[i]] <- bind_rows(lapply(ls.fit[[i]], `[[`, "idf"))
# }

# df.par <- bind_rows(ls.par)
# df.idf <- bind_rows(ls.idf)

# arrow::write_parquet(df.par, "base/gerados/resultados_por_duracao/df_par_by_d.pqt")
# arrow::write_parquet(df.idf, "base/gerados/resultados_por_duracao/df_idf_by_d.pqt")