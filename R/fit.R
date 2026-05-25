
# PACOTES -----------------------------------------------------------------

# Limpar ambiente
rm(list = ls()); invisible(gc())

# Pacotes
if(!require(pacman)) install.packages("pacman")
pacman::p_load(pacman,    # gerenciador de pacotes
               arrow,     # leitura de dados 'parquet'
               tidyverse)

# Set theme
suppressMessages(extrafont::loadfonts())
font <- "Consolas"
fontsize <- 8
textsize <- (fontsize - 2)/.pt
set_theme(theme_minimal() +
            theme(text = element_text(color = "black", size = fontsize, family = font),
                  axis.text = element_text(color = "black"),
                  panel.border = element_rect(color = "black")))


# FUNÇÕES -----------------------------------------------------------------

# Ajustar IDF d-GEV, estimar quantis e intervalos de confiança c/ bootstrap
source("R/funcoes/fitboot_dgev.R") # faz tudo e estima ICs p/ todas as estações
source("R/funcoes/fit_dgev.R")     # ajusta d-GEV p/ uma estação
source("R/funcoes/q_dgev.R")       # estima quantis p/ uma estação


# LER DADOS ---------------------------------------------------------------

# Intensidades máximas anuais de precipitação
df.imax <- arrow::read_parquet("base/gerados/df_imax.pqt")

# Aplicar critérios de qualidade
na.yr <- 0.2
na.mon <- 0.1
min.durations <- 3
d.interval <- c(1, 24)
min.yr <- 8

data <- df.imax %>% 
  filter(na_prct_yr <= na.yr,
         na_prct_mon <= na.mon,
         between(d, d.interval[1], d.interval[2])) %>% 
  group_by(gauge_code) %>% 
  filter(n_distinct(d) >= min.durations,
         n_distinct(wateryear) >= min.yr) %>% 
  ungroup() %>% 
  select(gauge_code, wateryear, d, imax)

ls.data <- split(data, data$gauge_code)

# CÁLCULO -----------------------------------------------------------------

# Argumentos
prior.info <- c(0.0178, 0.0627)
durations <- c(1, 2, 4, 6, 8, 10, 12, 18, 24)
return.period <- c(1.1, 1.5, 2, 5, 10, 20, 50, 100, 200, 500)
method <- "L-BFGS-B"
scale.inv <- "wide"
B <- 3000

# Ajustar d-GEV e estimar quantis e intervalos de confiança: com priori
ls.fit <- fitboot.dgev(ls.data, prior.info, durations, return.period, cols = names(data),
                       method, scale.inv, B, conf = 0.95, maxit = 1e8)

df.par <- bind_rows(lapply(ls.fit, `[[`, "par"))
df.idf <- bind_rows(lapply(ls.fit, `[[`, "idf"))

saveRDS(ls.fit, "base/gerados/ls_fit.rds")
arrow::write_parquet(df.par, "base/gerados/resultados_com_priori/df_par.pqt")
arrow::write_parquet(df.idf, "base/gerados/resultados_com_priori/df_idf.pqt")

# Ajustar d-GEV e estimar quantis e intervalos de confiança: sem priori
ls.fit <- fitboot.dgev(ls.data, prior.info = NULL, durations, return.period, cols = names(data),
                       method, scale.inv, B, conf = 0.95, maxit = 1e8)

df.par <- bind_rows(lapply(ls.fit, `[[`, "par"))
df.idf <- bind_rows(lapply(ls.fit, `[[`, "idf"))

saveRDS(ls.fit, "base/gerados/resultados_sem_priori/ls_fit.rds")
arrow::write_parquet(df.par, "base/gerados/resultados_sem_priori/df_par.pqt")
arrow::write_parquet(df.idf, "base/gerados/resultados_sem_priori/df_idf.pqt")