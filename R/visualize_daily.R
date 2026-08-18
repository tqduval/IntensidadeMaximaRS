
# Avaliar se o viés existente em estações longas do conjunto diário está presente
# quando reamostramos subamostras dessas estações curtas

rm(list = ls()); gc()

# PACOTES ----------------------------------------------------------------

pacman::p_load(
  pacman, tidyverse, arrow, readxl, lmom, sf, moments,
  ggspatial, patchwork, xtable, ggtext, ggrepel, geomtextpath
)

# Set theme
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

# LER DADOS --------------------------------------------------------------

# Dados espaciais
macrobh.sf <- read_sf("base/sig/sf_macro_bh_rs.gpkg", layer = "dissolvido")
rs.sf <- read_sf("base/sig/sf_rs.gpkg")
isoietas.pma.sf <- read_sf("base/sig/isoietas_pma_1977_2006/Isoietas_PMA_1977_2006.shp") # precipitação média anual
# isoietas.imax.sf <- 
daily.sf <- read_sf("base/sig/sf_daily.gpkg") |> 
  slice_tail(n = 1, by = gauge_code) |> 
  mutate(gauge_code = sub("^0", "", gauge_code)) |> 
  filter(!gauge_code %in% c("SBPA", "SBUG", "SBSM"))

subdaily.sf <- arrow::read_parquet("base/fonte/subdaily/df_subdaily_info.pqt") |> 
  st_as_sf(coords = c("long", "lat"), crs = st_crs(rs.sf)) |> 
  select(gauge_code, geometry)

# Dados de precipitação máxima diária Rio Grande do Sul
pmax.df <- arrow::read_parquet("base/fonte/daily/df_daily_pmax.pqt")
pmax.df <- pmax.df |> 
  filter(!gauge_code %in% c("SBPA", "SBUG", "SBSM")) |> 
  group_by(gauge_code) |> 
  mutate(n_years = n_distinct(year))

# Calcular `rain_yr_mean`
pmax.original <- readxl::read_excel("base/fonte/daily/Novos_Dados_Diarios_RS.xlsx") |> 
  mutate(gauge_code = sub("^0", "", gauge_code)) |> 
  filter(
    !gauge_code %in% c("SBPA", "SBUG", "SBSM"),
    gauge_code %in% unique(pmax.df$gauge_code)
  )

n_distinct(pmax.original$gauge_code)

info.df <- pmax.df |> 
  mutate(
    rain_yr_mean = pmax.original$yearly_rain[match(gauge_code, pmax.original$gauge_code)]
  ) |> 
  group_by(gauge_code) |> 
  reframe(
    min_yr = min(year, na.rm = TRUE),
    max_yr = max(year, na.rm = TRUE),
    state = first(state),
    elevation = first(elevation),
    long = first(long),
    lat = first(lat),
    responsible = first(responsible),
    n_years = first(n_years),
    rain_yr_mean = first(rain_yr_mean)
  ) |> 
  arrange(desc(n_years)) |> 
    mutate(
      id = row_number(),
      long_lat = sprintf("[%0.3f,%0.3f]", long, lat)
    )

unique(pmax.df$responsible)
n_distinct(pmax.df$gauge_code)
  
# Classes de tamanho das séries
n.years <- info.df$n_years
n.breaks <- 5
probs <- seq(0, 1, length.out = n.breaks + 1)
breaks <- quantile(n.years, probs, type = 1, names = FALSE) # type = 1 garante valores inteiros
# breaks <- c(31, 40, 60, 70, 80, max(n.years))
yr.interval <- cut(n.years, breaks, right = FALSE, include.lowest = TRUE)
info.df$yr_interval <- yr.interval
pmax.df <- pmax.df |> mutate(yr_interval = cut(n_years, breaks, right = FALSE, include.lowest = TRUE))
levels(pmax.df$yr_interval)


## MAPA ESTAÇÕES ---------------------------------------------------------

daily.sf <- daily.sf |> 
  mutate(id = info.df$id[match(gauge_code, info.df$gauge_code)])

bbox.rs <- st_bbox(rs.sf)

map.daily.id.plot <- ggplot() +
# geom_sf(data = rs.sf, fill = NA, color = "black", linewidth = 0.3) +
  geom_sf(data = macrobh.sf, color = "black", fill = NA, linewidth = 0.3) +
  geom_sf(data = subdaily.sf, color = "#1100ff", alpha = 0.2, size = 0.5) +
  geom_sf(data = daily.sf, color = "#ff0000", alpha = 0.8, size = 0.5) +
  geom_text_repel(
    aes(label = id, geometry = geom), daily.sf,
    stat = "sf_coordinates", family = font, size = textsize,
    color = "black", bg.color = "white", bg.r = 0.15, max.overlaps = Inf,
    min.segment.length = 0, segment.color = "grey50", segment.size = 0.2,
    box.padding = 0.2, point.padding = 0.1
  ) +
  coord_sf(crs = st_crs(4674), xlim = bbox.rs[c(1,3)], ylim = bbox.rs[c(2,4)], expand = TRUE, clip = "off") +
  annotation_scale(
    location = "br", width_hint = 0.2, height = unit(1.5, "mm"), text_family = font, text_cex = 0.6,
    bar_cols = c("black", "white"), line_width = 0.4
  ) +
  annotation_north_arrow(
    location = "tr" , width = unit(1, "cm"), height = unit(1, "cm"), which_north = "true",
    style = north_arrow_nautical(text_family = font, text_size = fontsize, fill = c("black", "white"))
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    panel.border = element_rect(linewidth = 0.3)
  ); map.daily.id.plot

ggsave("figuras/preliminares/mapa_daily_id.png", map.daily.id.plot, width = 12, height = 12, units = "cm", dpi = 300)



## TABELA INFORMAÇÕES ESTAÇÕES -------------------------------------------

daily.info.table <- info.df |> select(gauge_code, id, responsible, n_years, rain_yr_mean, elevation, long_lat)
table.names <- c("Código da estação", "ID", "Operador", "N", "PAA$^{*}$", "Elevação", "Coordenadas$^{**}$")
names(daily.info.table) <- table.names
n.col <- ncol(daily.info.table)

# Preparar tabela
daily.info.xtable <- daily.info.table |> 
  xtable(
    caption = "Informações básicas das estações pluviométricas do Rio Grande do Sul utilizadas.",
    label = "tab_daily_info",
    align = rep("c", n.col + 1)
  )

# Adicionar comandos adicionais no latex
extra.rows <- list()
extra.rows$pos <- list(0, nrow(daily.info.table)) # posições das novas linhas
ltx.cmd.units <- sprintf(                         # linha c/ unidades
  " - & - & - & anos & mm & m & º \\\\
  \\midrule
  \\endfirsthead
  \\caption[]{Informações básicas das estações pluviométricas do Rio Grande do Sul utilizadas. (continuação)} \\\\
  \\toprule
  Código da estação & ID & Operador & N & PAA$^{*}$ & Elevação & Coordenadas$^{**}$ \\\\
  - & - & - & anos & mm & m & º \\\\
  \\midrule
  \\endhead
  \\midrule
  \\multicolumn{%d}{r}{{Continua na próxima página...}} \\\\
  \\endfoot
  \\endlastfoot \n", 
  n.col
)
ltx.cmd.footnote <- sprintf( # duas linhas novas c/ notas de rodapé
  "\\bottomrule
  \\multicolumn{%d}{l}{$^{*}$Precipitação acumulada anual média entre anos com percentual de falhas inferior a 20\\%%.} \\\\
  \\multicolumn{%d}{l}{$^*$EPSG 4674.} \\\\ \n",
  n.col, n.col
)
extra.rows$command <- c(ltx.cmd.units, ltx.cmd.footnote)

daily.info.string <- print(
  daily.info.xtable, type = "latex", include.rownames = FALSE, caption.placement = "top",
  sanitize.colnames.function = identity, 
  sanitize.text.function = identity,
  tabular.environment = "longtable", # changes to longtable
  floating = FALSE,                  # required for longtable to work
  add.to.row = extra.rows,
  hline.after = c(-1)
)

cat(
  "\\begingroup\n\\footnotesize\n",
  daily.info.string,
  "\\endgroup\n",
  file = "tabelas/5_daily_info_table.txt",
  sep = ""
)



## TABELA ESTATÍSTICAS BÁSICAS -------------------------------------------

# Função que calcula a tabela resumo igual Papalexiou e Koutsoyiannis (2013)
# Originalmente summary.stats
get.stats <- function(x){
  gg.stats <- c(
    n_years = length(x),
    mean = mean(x, na.rm = TRUE),
    p0 = sum(x == 0)/length(x),
    median = median(x, na.rm = TRUE),
    sd = sd(x, na.rm = TRUE),
    t3 = lmom::samlmu(x)[[3]],
    t4 = lmom::samlmu(x)[[4]]
  )
}

# Originalmente fun.estat.bas
summary.stats <- function(x){
  c(
    min = min(x, na.rm = TRUE),
    q05 = quantile(x, 0.05, na.rm = TRUE),
    q25 = quantile(x, 0.25, na.rm = TRUE),
    mediana = median(x, na.rm = TRUE),
    q75 = quantile(x, 0.75, na.rm = TRUE),
    q95 = quantile(x, 0.95, na.rm = TRUE),
    max = max(x, na.rm = TRUE),
    media = mean(x, na.rm = TRUE),
    assimetria = moments::skewness(x, na.rm = TRUE)
  )
}


# Estatísticas básicas
# Aplicar resumo por estação
stats.df <- pmax.df |> 
  group_by(gauge_code) |> 
  group_modify(~as.data.frame(t(get.stats(.x$peak_24h)))) |>
  ungroup() |> 
  select(-p0)

# Aplicar resumo por estatística
stats.table <- apply(as.matrix(stats.df[,-1]), 2, summary.stats)
rownames(stats.table) <- c("Min", "$Q_5$", "$Q_{25}$", "$Q_{50}$", "$Q_{75}$", "$Q_{95}$", "Max", "Média", "Assimetria")
colnames(stats.table) <- c("$N$", "Média", "Mediana", "$SD$", "$\\hat\\tau_3$", "$\\hat\\tau_4$")
n.col <- ncol(stats.table)

stats.xtable <- xtable::xtable(
  x = stats.table,
  caption = sprintf("Estatísticas básicas do conjunto de %g séries de máximos anuais de precipitação diária", nrow(stats.df)),
  label = "tab_estatisticas_diarias",
  align = c("l", rep("c", ncol(stats.table))),
  digits = 3
)

stats.string <- print(
  stats.xtable, type = "latex", include.rownames = TRUE, caption.placement = "top",
  sanitize.colnames.function = identity,
  sanitize.text.function = identity,
  hline.after = c(-1)
)

cat(
  "\\begingroup\n\\footnotesize\n",
  stats.string,
  "\\endgroup\n",
  file = "tabelas/5_estatisticas_basicas_diarios.txt",
  sep = ""
)

# Localização da estação com registro mais longo
longest <- info.df$gauge_code[which.max(info.df$n_years)]
pmax.original[pmax.original$gauge_code == longest, c("gauge_code", "city")]
info.df$id[info.df$gauge_code == longest]

# Localização estação com maior peak_24h registrado
max.pmax <- pmax.df$gauge_code[which.max(pmax.df$peak_24h)] |> print()
pmax.original[pmax.original$gauge_code == max.pmax, c("gauge_code", "city")]
info.df$id[info.df$gauge_code == max.pmax]

# Estações com os valores mais altos e mais baixos de L-curtose
high.kt <- stats.df$gauge_code[which.min(stats.df$t4)]
pmax.original[pmax.original$gauge_code == high.kt, c("gauge_code", "city")]
info.df$id[info.df$gauge_code == high.kt]

# Razões negativas
stats.df[stats.df$t3 < 0,]
stats.df[stats.df$t4 < 0,]


# L-MOMENTOS -------------------------------------------------------------

flmom <- function(x){
  x <- na.omit(x)
  sample.lmom <- lmom::samlmu(x)
  fit <- lmom::pelgev(sample.lmom)
  return(c(sample.lmom[c(3,4)], fit))
}

pmax.wide <- pmax.df |> 
  select(!c(state, elevation, long, lat, responsible, network, n_years, yr_interval)) |> 
  pivot_wider(
    id_cols = year,
    names_from = gauge_code,
    values_from = peak_24h
  ) |> 
  arrange(year) |> 
  column_to_rownames("year") |> 
  as.matrix()

fit <- apply(pmax.wide, 2, flmom); invisible(gc())

fit.df <- t(fit) |>
  as.data.frame() |> 
  rownames_to_column(var = "gauge_code")


# VISUALIZAR -------------------------------------------------------------

fit.plot <- fit.df |> 
  left_join(
    info.df |> 
      select(gauge_code, n_years, yr_interval),
    by = "gauge_code"
  )

yr.col <- c("#e2fae2", "#B4EEB4", "#698B69", "#00688B", "#191970")
yr.col <- setNames(yr.col, levels(yr.interval)) |> print()
# levels(daily.sf$yr_interval)

# Boxplot kappa por faixa de duração
kappa.boxplot.nyr <- fit.plot |> 
  ggplot(aes(x = yr_interval, y = k, fill = yr_interval)) +
  geom_jitter(color = "grey90", alpha = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey80") +
  geom_boxplot(
    color = "black", median.linewidth = 0.2, staplewidth = 0.5, linewidth = 0.2,
    whisker.linetype = "dashed", outlier.size = 0.5, outlier.alpha = 0.5
  ) +
  scale_fill_manual(
    name = "", values = yr.col
  ) +
  labs(x = "N [anos]", y = "κ") + 
  theme(
    legend.position = "none"
  ); kappa.boxplot.nyr

ggsave("figuras/preliminares/diarios/kappa_boxplot_nyr.png", kappa.boxplot.nyr, width = 14, height = 8, units = "cm", dpi = 300)

# Mapa localização estações por classes de tamanho das séries
daily.sf <- daily.sf |> 
  mutate(yr_interval = info.df$yr_interval[match(gauge_code, daily.sf$gauge_code)])

map.rs.plot <- ggplot() +
  geom_sf(data = rs.sf, fill = "NA", color = "black", linewidth = 0.3) +
  geom_sf(data = daily.sf, aes(color = yr_interval), size = 0.5) +
  coord_sf(crs = st_crs(4674), expand = TRUE) +
  scale_color_manual(name = "N [anos]", values = yr.col) +
  # scale_color_discrete(name = "", palette = yr.col, guide = guide_colorsteps()) +
  annotation_scale(
    location = "br", width_hint = 0.2, height = unit(1.5, "mm"), text_family = font, text_cex = 0.6,
    bar_cols = c("black", "aliceblue"), line_width = 0.4
  ) +
  annotation_north_arrow(
    location = "tr" , width = unit(1, "cm"), height = unit(1, "cm"), which_north = "true",
    style = north_arrow_nautical(text_family = font, text_size = fontsize, fill = c("black", "aliceblue"))
  ) +
  labs(x = "Longitude", y = "Latitude") +
  guides(
    color = guide_legend(
      override.aes = list(shape = 15, size = 4),
    )
  ) +
  theme(
    # legend.title = element_blank(),
    legend.key.spacing.y = unit(-3, "pt"),
    # legend.key.height = unit(fontsize*0.8, "pt"),
    # legend.key.width = unit(fontsize*0.6, "pt"),
    legend.position = c(0, -0.02),
    legend.justification = c(0, 0)
  ); map.rs.plot

ggsave("figuras/preliminares/diarios/map_rs_plot_nyr.png", map.rs.plot, width = 9, height = 9, units = "cm", dpi = 300)

map.boxplot.nyr.kappa <- 
  map.rs.plot + labs(subtitle = "(a)") +
  kappa.boxplot.nyr + labs(subtitle = "(b)") +
  plot_layout(widths = c(1, 0.8), guides = "keep"); map.boxplot.nyr.kappa

ggsave("figuras/preliminares/diarios/map_boxplot_nyr.png", map.boxplot.nyr.kappa, width = 14, height = 7.5, units = "cm", dpi = 300)


## ISOIETAS --------------------------------------------------------------

isoietas.pma.rs.sf <- st_intersection(isoietas.pma.sf, rs.sf)
isoietas.xyz <- st_coordinates(isoietas.pma.rs.sf)

# Isoietas de precipitação anual média
ggplot() +
  geom_sf(data = rs.sf, fill = "NA", color = "black", linewidth = 0.3) +
  # geom_sf(data = isoietas.pma.rs.sf, linewidth = 0.3) +
  geomtextpath::geom_textsf(
    data = isoietas.pma.rs.sf, aes(label = CONTOUR),
    family = font, size = textsize, color = "black", linewidth = 0.3,
    vjust = 0.5, hjust = 0.2, halo = TRUE, upright = TRUE, straight = FALSE
  ) +
  coord_sf(crs = st_crs(4674), expand = TRUE) +
  annotation_scale(
    location = "br", width_hint = 0.2, height = unit(1.5, "mm"), text_family = font, text_cex = 0.6,
    bar_cols = c("black", "aliceblue"), line_width = 0.4
  ) +
  annotation_north_arrow(
    location = "tr" , width = unit(1, "cm"), height = unit(1, "cm"), which_north = "true",
    style = north_arrow_nautical(text_family = font, text_size = fontsize, fill = c("black", "aliceblue"))
  ) +
  labs(x = "Longitude", y = "Latitude")

ggsave("figuras/preliminares/mapa_isoietas_pma.png", width = 12, height = 12, units = "cm", dpi = 300)


## MONTE CARLO -----------------------------------------------------------

set.seed(1234)

N <- 80 # tamanho mínimo das séries
subsample.lengths <- seq(10, max(n.years), 5) # diferentes comprimentos de subamostras
R <- 10000 # número de repetições

gauges.N <- info.df |> 
  filter(n_years >= N) |> 
  arrange(n_years) |>
  pull(gauge_code) |> 
  print()

flmom.boot <- function(indices, x){
  x <- x[indices]
  flmom(x)["k"]
}

m.of.n.boot <- function(x, m, R, replace = FALSE){
  n <- length(x)
  idx <- replicate(R, sample(1:n, size = m, replace))
  x <- x[idx]
  est <- apply(idx, 2, flmom.boot, x = x)
  return(est)  
}

pmax.N <- pmax.wide[,gauges.N]
# pmax.ls <- asplit(pmax.N, 2)
pmax.ls <- apply(pmax.N, 2, na.omit, simplify = FALSE)

# Combinações
grid <- expand.grid(
  gauge_code = gauges.N,
  m = subsample.lengths
)

shape.ls <- pbapply::pblapply(seq_len(nrow(grid)), function(i){

  gg <- grid$gauge_code[i]
  m <- grid$m[i]
  x <- pmax.ls[[gg]]

  est <- tryCatch(
    m.of.n.boot(x, m, R, replace = TRUE),
    error = function(e) rep(NA_real_, R)
  )

  out <- data.frame(
    gauge_code = gg,
    m = m,
    replicate = 1:R,
    kappa = est
  )

  return(out)

})

shape.df <- bind_rows(shape.ls)

# Localizar as estações selecionadas
ggplot() +
  geom_sf(data = rs.sf, fill = "NA", color = "black", linewidth = 0.3) +
  geom_sf(data = daily.sf[daily.sf$gauge_code %in% gauges.N,], aes(color = yr_interval), size = 0.5) +
  coord_sf(crs = st_crs(4674), expand = TRUE) +
  scale_color_manual(name = "", values = yr.col) +
  # scale_color_steps(name = "", values = unname(yr.col)) +
  annotation_scale(
    location = "br", width_hint = 0.2, height = unit(1.5, "mm"), text_family = font, text_cex = 0.6,
    bar_cols = c("black", "aliceblue"), line_width = 0.4
  ) +
  annotation_north_arrow(
    location = "tr" , width = unit(1, "cm"), height = unit(1, "cm"), which_north = "true",
    style = north_arrow_nautical(text_family = font, text_size = fontsize, fill = c("black", "aliceblue"))
  ) +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    legend.title = element_blank(),
    margins = margin(l = 5.5),
    # legend.position = c(1, 0),
    # legend.justification = c(1, 0),
    # legend.direction = "vertical",
    legend.key.height = unit(fontsize*0.8, "pt")
    # legend.spacing.x = unit(fontsize*2, "pt")
    # legend.spacing.y = unit(0, "pt")
  )

# Estão relativamente bem espalhadas pelo estado

# Estimativas originais
shape.original <- t(apply(pmax.N, 2, flmom)) |> 
  as.data.frame() |> 
  rownames_to_column(var = "gauge_code") |> 
  mutate(
    n_years = info.df$n_years[match(gauge_code, info.df$gauge_code)],
    label = sprintf("%s | %g anos", gauge_code, n_years)
  )


# Avaliar o comportamento das séries ao longo dos valores de `m'`
kappa.boot.boxplot <- shape.df |> 
  mutate(
    n_years = info.df$n_years[match(gauge_code, info.df$gauge_code)],
    label = sprintf("%s | %g anos", gauge_code, n_years)
  ) |> 
  filter(m %in% seq(10, max(n.years), 10)) |> 
  ggplot(aes(x = factor(m), y = kappa)) +
  facet_wrap(~fct_reorder(label, n_years, .desc = TRUE), ncol = 1, scales = "free_y") +
  # geom_jitter(color = "grey90", alpha = 0.5, size = 0.8) +
  geom_hline(aes(yintercept = k), shape.original, color = "red", linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey80") +
  geom_boxplot(
    color = "black", fill = "#191970", median.linewidth = 0.2, staplewidth = 0.5, linewidth = 0.2,
    whisker.linetype = "dashed", outlier.size = 0.5, outlier.alpha = 0.5
  ) +
  labs(x = "N [anos]", y = "κ") + 
  theme(
    legend.position = "none"
  ); kappa.boot.boxplot

ggsave("figuras/preliminares/diarios/kappa_boot_boxplot.png", kappa.boot.boxplot, width = 14, height = 20, units = "cm", dpi = 300)

# Avaliar o comportamento do conjunto completo de estações
median(shape.original$k)
kappa.boot.boxplot.grouped <- shape.df |> 
  mutate(
    n_years = info.df$n_years[match(gauge_code, info.df$gauge_code)],
    label = sprintf("%s | %g anos", gauge_code, n_years)
  ) |> 
  filter(m %in% seq(10, max(n.years), 10)) |> 
  ggplot(aes(x = factor(m), y = kappa)) +
  geom_jitter(color = "grey90", alpha = 0.5, size = 0.8) +
  geom_hline(aes(yintercept = median(k)), shape.original, color = "red", linetype = "dashed", linewidth = 0.3) +
  geom_hline(yintercept = c(-0.5, 0.5), linewidth = 0.3, color = "grey80", linetype = "dashed") +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey80") +
  geom_boxplot(
    color = "black", fill = "#191970", median.linewidth = 0.2, staplewidth = 0.5, linewidth = 0.2,
    whisker.linetype = "dashed", outlier.size = 0.5, outlier.alpha = 0.5
  ) +
  labs(x = "N [anos]", y = "κ") + 
  theme(
    legend.position = "none"
  ); kappa.boot.boxplot.grouped

ggsave("figuras/preliminares/diarios/kappa_boot_boxplot_grouped.png", kappa.boot.boxplot.grouped, width = 14, height = 8, units = "cm", dpi = 300)


# GERAR SÉRIES SINTÉTICAS ------------------------------------------------

set.seed(1234)

# Gerar séries sintéticas com diferentes valores de parâmetro de forma e simular
# subamostras de tamanhos diferentes e avaliar o tamanho da amostra necessário
# para se aproximar do valor original do parâmetro de forma
source("R/funcoes/fgev.R")

# Avaliar valores plausíveis 
gg <- "83914" # estação mais longa do conjunto
gg <- info.df$gauge_code[info.df$id == 154]

ams <- pmax.df |> 
  filter(gauge_code == gg) |> 
  pull(peak_24h)

flmom(ams)

fit.df |> reframe(across(everything(), mean))

# Definir parâmetros originais
loc.scale <- c(90, 20)
shape <- c(-0.2, -0.1, -0.05, 0, 0.1)
psi <- matrix(loc.scale, nrow = length(shape), ncol = 2, byrow = TRUE)
psi <- cbind(psi, unname(shape))
colnames(psi) <- c("loc", "scale", "shape")

# Probabilidades p/ gerar a amostra original
N <- 1e4 # anos
p <- runif(N, 0, 1)

# Gerar séries "originais" p/ diferentes parâmetros de forma
samples <- apply(psi, 1, function(par){
  lmom::quagev(p, par)
}) # retornar uma matriz
colnames(samples) <- shape

nrow(samples)

samples.df <- data.frame(
  shape_id = rep(seq_along(shape), each = N),
  shape = rep(shape, each = N),
  value = as.vector(samples)
)

# Histogramas das amostras geradas
ggplot(samples.df, aes(x = value)) +
  facet_wrap(~shape, ncol = 1, labeller = as_labeller(function(x) paste("κ =", x))) +
  geom_histogram(aes(y = ..density..), fill = "white", color = "black")

smin <- min(samples)
smax <- max(samples)
density.grid <- expand.grid(x = seq(smin, smax, length.out = N), shape = psi[,3])
dgev.df <- merge(density.grid, as.data.frame(psi), by = "shape") |> 
  rowwise() |> 
  mutate(
    density = evd::dgev(x, loc, scale, -shape),
    label = fct_reorder(paste("κ =", shape), shape)
  ) |> 
  ungroup()

colors <- c("black", "grey30", "grey70")
ltypes <- c("solid", "dashed", "solid")
dgev.df |> 
  filter(shape %in% c(-0.2, -0.05, 0.1)) |> 
  ggplot(aes(x = x, y = density, color = label, linetype = label)) +
  geom_line(linewidth = 0.3) +
  labs(y = "f<sub>X</sub>(x)") +
  scale_color_manual(values = colors) +
  scale_linetype_manual(values = ltypes) +
  xlim(smin, 300) +
  theme(
    legend.title = element_blank(),
    legend.position = c(1, 0.98),
    legend.justification = c(1, 1),
    axis.title.y = ggtext::element_markdown()
  )

ggsave("figuras/preliminares/diarios/densidades_gev.png", width = 12, height = 6, units = "cm", dpi = 300)

# Tamanhos das subamostras
sample.sizes <- c(300, 200, seq(100, 10, -5))
grid <- expand.grid(idx = seq_len(ncol(samples)), N = sample.sizes)
R <- 1e4

shape.ls <- lapply(seq_len(nrow(grid)), function(i){

  idx <- grid[i,"idx"]
  m <- grid[i, "N"]
  x <- samples[,idx]
  original.shape <- colnames(samples)[idx]
  message(sprintf("κ = %-5s | N = %-3d anos", original.shape, m))

  est <- tryCatch(
    m.of.n.boot(x, m, R, replace = TRUE),
    error = function(e) rep(NA_real_, R)
  )

  out <- data.frame(
    kappa_original = original.shape,
    m = m,
    replicate = 1:R,
    kappa = est
  )

  return(out)

}); invisible(gc())

shape.df <- bind_rows(shape.ls)

shape.df |> 
  mutate(
    kappa_original = as.numeric(kappa_original),
    label = fct_reorder(paste("κ =", kappa_original), kappa_original)
  ) |> 
  filter(m %in% seq(10, max(sample.sizes), 10)) |> 
  ggplot(aes(x = factor(m), y = kappa)) +
  facet_wrap(~label, ncol = 1, scales = "free_y") +
  geom_jitter(color = "grey90", alpha = 0.2, size = 0.5) +
  geom_hline(yintercept = 0, linewidth = 0.3, color = "grey80") +
  geom_boxplot(
    color = "black", fill = "grey80", median.linewidth = 0.2, staplewidth = 0.5, linewidth = 0.2,
    whisker.linetype = "dashed", outlier.size = 0.3, outlier.alpha = 0.5
  ) +
  labs(x = "N [anos]", y = "κ") + 
  theme(
    legend.position = "none"
  )

ggsave("figuras/preliminares/diarios/sim_kappa_boxplot.png", width = 14, height = 20, units = "cm", dpi = 300)
