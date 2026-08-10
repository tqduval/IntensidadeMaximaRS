
# PACOTES -----------------------------------------------------------------

# Limpar ambiente
rm(list = ls()); invisible(gc())

# Pacotes
if(!require(pacman)) install.packages("pacman")
pacman::p_load(pacman,    # gerenciador de pacotes
               arrow,     # leitura de dados 'parquet'
               tidyverse,
               parallel,
               beepr)


# LER DADOS ---------------------------------------------------------------

# Importar
info.df <- arrow::read_parquet("base/fonte/subdaily/df_subdaily_info.pqt")
data.df <- arrow::read_parquet("base/fonte/subdaily/df_subdaily_data.pqt")

# Critérios p/ filtrar séries
na.accept <- 0 # 0.2
min.years <- 0 # 8

# Funções
source("R/funcoes/fun_filter_set.R") # preencher datas e filtrar
source("R/funcoes/fun_group_ts.R")   # agrupar por time_step
source("R/funcoes/fun_imax_agg_wateryear.R") # agregar por durações e calcular intensidade máxima anual por ano hidrológico

# Preencher datas e filtrar conforme 'na_accept' e 'min_years'
data.df <- fun_filter_set(
  data = data.df,
  daily = FALSE,
  col_names = names(data.df),
  filter = FALSE
); invisible(gc())


# CALCULAR TOTAIS ANUAIS -------------------------------------------------

yr.data <- data.df |> 
  mutate(
    year = year(datetime)
  ) |> 
  group_by(gauge_code, year) |> 
  reframe(
    rain_yr = sum(rain_mm, na.rm = TRUE)
  )

arrow::write_parquet(yr.data, "base/gerados/df_subdaily_yrsum.pqt")

# AGREGAR DURAÇÕES --------------------------------------------------------

# Agrupar séries por 'time_steps'
data.ls <- split(df.data, f = df.data$gauge_code)            # lista de data.frames de duração
data.by.time <- fun_group_ts(data.ls, ts_name = "time_step") # agrupar em durações diferentes numa lista

# Durações
d.subhour <- c(10, 15, 20, 30, 40, 45, 50) # durações sub-horárias
d.subdaily <- seq(1, 23)                   # durações subdiárias
d.daily <- 1:10                            # durações diárias
durations <- c(d.subhour/60, d.subdaily, d.daily*24)  # durações
# which.mon <- c(8, 9, 10, 11, 12, 1, 2, 3, 4, 5, 6, 7) # ano hidrológico
which.mon <- 1:12                                     # ano civil
time.steps <- names(data.by.time)          # lista de time_steps

# Aplicar função 'fun_imax_agg' p/ cada grupo de durações
imax.ls <- lapply(X = time.steps, FUN = function(ts){
  
  current.ts <- data.by.time[[ts]] # extrair listas de estações c/ resolução 'ts'
  ds <- as.numeric(ts)/3600        # segundos p/ horas
  valid.durations <- durations[(durations/ds) %% 1 < 1e-8] # selecionar somente durações compatíveis c/ a resolução 'ts'
  
  message(paste0("\nProcessando conjunto de estações com resolução de ", ds, " minutos..."))
  
  # A função fun_imax_agg() foi substituída pela função imax.wateryear
  # que calcula imax agora considerando o ano hidrológico estabelecido
  # conforme `which.mon` (pegando a primeira observação `which.mon[1]`)
  # como o início do ano hidrológico
  out <- imax.wateryear(
    data = current.ts,
    durations = valid.durations,
    which.mon = which.mon,
    names = c("datetime", "rain_mm")
  )
  
}) # fim 'imax.ls'

# Agosto foi definido como início do ano hidrológico baseado nas avaliações feitas
# em `dry_months_frequency.R`, que possibilitou identificar meses com baixa frequência
# de intensidades máximas anuais

rm(data.df, data.ls, data.by.time); invisible(gc())

# Nomear lista
names(imax.ls) <- time.steps

# Salvar resultados
# saveRDS(object = imax.ls, file = "base/gerados/imax_ls.rds")
imax.df <- bind_rows(imax.ls)
arrow::write_parquet(x = bind_rows(imax.ls), sink = "base/gerados/df_imax.pqt")

# Ler arquivo c/ intensidades máximas
imax.df <- arrow::read_parquet(file = "base/gerados/df_imax.pqt")

quality.df <- imax.df |> 
  mutate(responsible = info.df$responsible[match(gauge_code, info.df$gauge_code)]) |> 
  select(gauge_code, responsible, wateryear, na_prct_yr) |> 
  group_by(gauge_code, responsible, wateryear) |> 
  summarise(na_prct = first(na_prct_yr)*100)

thresholds <- seq(10, 40, 10)
min.years <- seq_len(10)

grid <- expand.grid(threshold = thresholds, min_years = min.years) |> # combinações 
  mutate(n_gauges = purrr::map2_int(.x = threshold, .y = min_years, .f = ~{ # analisar qtas estações p/ cada par
    quality.df |>
      filter(na_prct <= .x) |> # filtrar estações c/ no máximo `threshold` falhas anuais
      count(gauge_code) |>     # contar qtas estações
      filter(n >= .y) |>       # filtrar estações c/ no mínimo `min_years` anos
      nrow()                   # contar qtas estações
  }))

resp.names <- paste(unique(info.df$responsible), collapse = ", ")

# Visualização

# Set theme
suppressMessages(extrafont::loadfonts())
font <- "JetBrainsMono NF"
# font <- "Courier New"
fontsize <- 8
textsize <- (fontsize - 2) / .pt
ggplot2::set_theme(
  theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box.margin = margin(-5.5, 2, 2, 2),
      legend.key.size = unit(fontsize, "pt"),
      text = element_text(color = "black", size = fontsize, family = font),
      axis.text = element_text(color = "black"),
      strip.text = element_text(hjust = 0, size = fontsize - 2, margin = margin(2, 2, 2, 2), color = "black"),
      panel.border = element_rect(color = "black", linewidth = 0.5),
      panel.grid = element_blank(),
      plot.margin = margin(0.5, 0.5, 0.5, 0.5)
    )
)

color <- c("green4", "lightgreen", "orange", "#f10000")
anos.estacao.plot <- grid |> 
  mutate(threshold = factor(threshold, levels = thresholds)) |> 
  ggplot(aes(x = min_years, y = n_gauges, fill = threshold)) +
  facet_wrap(~threshold, nrow = 4, labeller = labeller(threshold = function(x) paste0(x, "% falhas"))) +
  geom_col(show.legend = FALSE) +
  geom_text(aes(label = n_gauges), vjust = -0.3, size = textsize, family = font) +
  scale_x_continuous(breaks = seq(min(grid$min_years), max(grid$min_years), by = 1)) +
  scale_y_continuous(limits = c(0, max(grid$n_gauges) + 50)) +
  scale_fill_manual(values = color) +
  labs(x = "Número de anos por estação",y = "Número de estações"); anos.estacao.plot

ggsave("figuras/preliminares/plot_anos_estacao.png", anos.estacao.plot, width = 15, height = 15, units = "cm", dpi = 300)
