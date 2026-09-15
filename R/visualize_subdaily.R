
rm(list = ls()); invisible(gc())

pacman::p_load(
  pacman, tidyverse, arrow, 
  ggtext, ggrepel, ggspatial, patchwork, extrafont, ggpattern, ggpp, ggnewscale,
  sf, rnaturalearth, tidyterra, spdep,
  lmom, moments,
  xtable
)


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


# FUNÇÕES ----------------------------------------------------------------


# Ajuste modelos
source("R/funcoes/fdgev.R") # d-GEV
source("R/funcoes/fgev.R")
# source("R/funcoes/idf_sgb.R")

# Quantil d-GEV
qdgev <- function(param, duration, p.nexceed){
  
  param <- unname(param)
  xit <- param[1]
  alpha0 <- param[2]
  kappa <- param[3]
  H <- param[4]
  theta <- param[5]
  
  p <- p.nexceed; d <- duration
  
  # Calcular parâmetros convencionais da GEV
  alpha <- alpha0/(d + theta)^H
  xi <- xit*alpha
  
  yp <- xi + alpha/kappa*(1 - (-log(p))^kappa)
  
  return(yp)
  
}

# Função p/ dissolver bordas
# dissolve_boundaries <- function(x) {
#   if (!requireNamespace("spdep")) {
#     stop("spdep is needed to use this function")
#   }

#   crs <- sf::st_crs(x)
#   nbs <- spdep::poly2nb(x)
#   nb_ids <- attr(nbs, "ncomp")[["comp.id"]]
#   res <- lapply(
#     split(x, nb_ids),
#     sf::st_union,
#     is_coverage = TRUE
#   )

#   sf::st_sfc(
#     unlist(res, recursive = FALSE),
#     crs = crs
#   ) |> sf::st_combine()

# }

# LER DADOS --------------------------------------------------------------

na.yr <- 0.2
min.durations <- 3
d.interval <- c(1, 24)
min.yr <- 8

# Dados observados
imax.df.prefilter <- arrow::read_parquet("base/gerados/df_imax.pqt")
message(sprintf("Número de estações antes de aplicar os filtros: %g", n_distinct(imax.df.prefilter$gauge_code)))

yrsum.df <- arrow::read_parquet("base/gerados/df_subdaily_yrsum.pqt")

# Aplicar critérios de qualidade
imax.df <- imax.df.prefilter |> 
  filter(na_prct_yr <= na.yr,
         between(d, d.interval[1], d.interval[2])) |> 
  group_by(gauge_code) |> 
  mutate(
    n_years = n_distinct(wateryear),
    gauge_label = sprintf("%s | %2g anos", gauge_code, n_years)
  ) |> 
  filter(n_distinct(d) >= min.durations,
         n_years >= min.yr) |> 
  ungroup() |> 
  left_join(yrsum.df, by = join_by(gauge_code, wateryear == year)) |> 
  select(gauge_code, wateryear, d, imax, rain_yr, n_years, gauge_label)

message(sprintf("Número de estações após aplicar os filtros: %g", n_distinct(imax.df$gauge_code)))

gauges <- unique(imax.df$gauge_code)

info.df <- arrow::read_parquet("base/fonte/subdaily/df_subdaily_info.pqt")
info.df <- info.df |> 
  filter(gauge_code %in% gauges)


## GEOESPACIAIS ----------------------------------------------------------

br.sf <- read_sf("base/sig/sf_uf.gpkg")
world.sf <- rnaturalearth::ne_countries(continent = "south america", scale = "large")
rs.sf <- read_sf("base/sig/sf_rs.gpkg")
bhrs.sf <- read_sf("base/sig/sf_hidro_rs.gpkg")
macrobh.sf <- read_sf("base/sig/sf_macro_bh_rs.gpkg", layer = "dissolvido")
macrobh.sf <- macrobh.sf |> mutate(macro_name = ifelse(macro_name == "Coastal", "Litorênea", macro_name))

daily.sf <- read_sf("base/sig/sf_daily.gpkg")
daily.sf <- daily.sf |> slice_head(by = "gauge_code")
subdaily.sf <- info.df |> 
  st_as_sf(coords = c("long", "lat"), crs = st_crs(rs.sf)) |> 
  select(gauge_code, geometry)

raster.files <- list.files(path = "base/sig/anadem_rs/", pattern = ".tif$", full.names = TRUE)
raster.collection <- lapply(raster.files, terra::rast)
raster.anadem <- terra::sprc(raster.collection)
raster.anadem <- terra::mosaic(raster.anadem)

anadem <- terra::crop(raster.anadem, rs.sf, mask = TRUE)
# hist(anadem)
# terra::plot(anadem)


## d-GEV -----------------------------------------------------------------

par.dgev.df <- read_parquet("base/gerados/df_par_dgev.pqt")
idf.dgev.df <- read_parquet("base/gerados/df_idf_dgev.pqt")

# Atribuir nomes aos dados em `prior`
which.prior <- unique(par.dgev.df$prior)
# prior.names <- setNames(c("MLE", "RGMLE", "GMLE"), which.prior) # GMLE aqui é a não informativa
prior.names <- paste(c("MLE", "RGMLE", "GMLE", "GMLE"), "–", gsub("c", "B", which.prior))
prior.names <- gsub(" – NULL", "", prior.names)
prior.names <- setNames(prior.names, which.prior)
prior.names <- prior.names[c(1, 4, 3, 2)]
par.dgev.df <- par.dgev.df |> mutate(model = prior.names[prior])
idf.dgev.df <- idf.dgev.df |> mutate(model = prior.names[prior])


## GEV POR DURAÇÃO -------------------------------------------------------

# Estimativas do parâmetro de forma com bootstrap e L-momentos
# shape_ci_boot_only_lmom.R
shape.boot.df <- arrow::read_parquet("base/gerados/df_shapes_boot_lmom.pqt")
par.gev.df <- read_parquet("base/gerados/df_par_gev.pqt")
idf.gev.df <- read_parquet("base/gerados/df_idf_gev.pqt")
par.gev.df <- par.gev.df |> mutate(model = prior.names[prior])
idf.gev.df <- idf.gev.df |> mutate(model = prior.names[prior])



# VISUALIZAÇÃO -----------------------------------------------------------

# Ordem de plotagem
order.n <- imax.df |> 
  slice_head(n = 1, by = "gauge_code") |> 
  select(gauge_code, n_years, gauge_label) |> 
  arrange(desc(n_years)) |> 
  pull(gauge_label, name = n_years)


# Rótulos p/ facets
gauge.labels <- imax.df |> 
  slice_head(n = 1, by = "gauge_code") |> 
  select(gauge_code, n_years, gauge_label) |> 
  arrange(desc(n_years)) |> 
  pull(gauge_label, name = gauge_code) |> 
  factor(levels = order.n)

# Identificadores p/ estações
gauge.ids <- imax.df |> 
  slice_head(n = 1, by = "gauge_code") |> 
  select(gauge_code, gauge_label, n_years) |> 
  arrange(desc(n_years)) |>
  mutate(id = row_number()) |> print()

# Cores
which.model <- unique(par.dgev.df$model)
col <- c("#648FFF", "#dc267f", "#FE6100")



## MAPA ESTAÇÕES ---------------------------------------------------------

bbox.rs <- st_bbox(rs.sf)
bbox.br <- st_bbox(br.sf)
bbox.br[3] <- -35
macrobh <- macrobh.sf$macro_name

# Brasil
plot.map.br <- ggplot() +
  geom_sf(data = world.sf, color = "grey30", fill = "grey90") +
  geom_sf(data = br.sf[br.sf$SIGLA_UF != "RS",], color = "grey10", fill = "lightgreen", linewidth = 0.2) +
  # geom_sf(data = sf.uf, color = "grey30", fill = NA, linetype = "dashed", linewidth = 0.15, alpha = 0.7) +
  geom_sf(data = rs.sf, color = "grey10", fill = "green3", linewidth = 0.4) +
  geom_sf(data = daily.sf, aes(color = "Pluviômetros"), size = 0.05, alpha = 0.8) +
  geom_sf(data = subdaily.sf, aes(color = "Pluviógrafos"), size = 0.05, alpha = 0.8) +
  guides(color = guide_legend(override.aes = list(size = 1.5))) +
  annotation_scale(
    location = "br", width_hint = 0.2, height = unit(1.5, "mm"), text_family = font, text_cex = 0.6,
    bar_cols = c("black", "aliceblue"), line_width = 0.4
  ) +
  annotation_north_arrow(
    location = "tr" , width = unit(1, "cm"), height = unit(1, "cm"), which_north = "true",
    style = north_arrow_nautical(text_family = font, text_size = fontsize, fill = c("black", "aliceblue"))
  ) +
  scale_color_manual(
    name = "", values = c("Pluviômetros" = "#fc0000", "Pluviógrafos" = "#1100ff"),
    guide = guide_legend(order = 1)
  ) +
  labs(x = "Longitude", y = "Latitude") +
  coord_sf(crs = st_crs(4674), xlim = bbox.br[c(1,3)], ylim = bbox.br[c(2,4)]) +
  theme(plot.margin = margin(t = 10),
        panel.background = element_rect(fill = "aliceblue"),
        legend.key.spacing.x = unit(fontsize*2, "pt"),
        legend.margin = margin(t = -5, l = -2),
        legend.spacing.x = unit(-4, "pt"),
        legend.position = "bottom"); plot.map.br

# anadem.values <- values(anadem)
# anadem.range <- range(anadem.values, na.rm = TRUE)
# # anadem.breaks <- pretty(anadem.range, n = 3)
# anadem.breaks <- round(sort(c(anadem.range, 0, 500, 1000)),0)

plot.map.rs <- ggplot() +
  tidyterra::geom_spatraster(data = anadem, alpha = 0.7) +
  geom_sf(data = rs.sf, fill = "NA", color = "black", linewidth = 0.3) +
  # geom_sf(data = macrobh.sf, color = "black", fill = NA, linewidth = 0.3) +
  geom_sf(data = daily.sf, aes(color = "Pluviômetros"), size = 0.5, show.legend = FALSE) +
  geom_sf(data = subdaily.sf, aes(color = "Pluviógrafos"), size = 0.5, show.legend = FALSE) +
  annotate(geom = "text", x = bbox.rs["xmin"], y = bbox.rs["ymax"], label = "(a)", family = font, size = textsize) +
  coord_sf(crs = st_crs(4674), expand = TRUE) +
  scale_color_manual(
    name = "", values = c("Pluviômetros" = "#fc0000", "Pluviógrafos" = "#1100ff"),
    guide = guide_legend(order = 1)
  ) +
  scale_fill_viridis_c(
    name = "", na.value = NA, option = "cividis", # option = "turbo") +
    guide = guide_colorsteps(show.limits = TRUE, order = 2)
  ) +
  theme(
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0, 0),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    legend.key.height = unit(fontsize*0.9, "pt"),
    plot.margin = margin(l = -10)
    # legend.spacing.x = unit(fontsize*2, "pt")
    # legend.spacing.y = unit(0, "pt")
  ); plot.map.rs

plot.map.bh <- ggplot(macrobh.sf) +
  # geom_sf(color = "black", fill = NA, linewidth = 0.3) +
  ggpattern::geom_sf_pattern(
    aes(pattern_color = macro_name, pattern_angle = macro_name),
    pattern = "line", fill = NA, color = "black", linewidth = 0.3, pattern_size = 0.2, pattern_spacing = 0.01
  ) +
  annotate(geom = "text", x = bbox.rs["xmin"], y = bbox.rs["ymax"], label = "(b)", family = font, size = textsize) +
  scale_pattern_color_manual(name = "", values = setNames(c("#0013a0", "#098f86", "#0a9fc0"), macrobh)) +
  scale_pattern_angle_manual(name = "", values = setNames(c(0, 0, 90), macrobh)) +
  coord_sf(crs = st_crs(4674), expand = TRUE) +
  theme(
    panel.border = element_blank(),
    axis.text = element_blank(),
    axis.title = element_blank(),
    legend.title = element_blank(),
    legend.position = c(0, 0),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    legend.key.height = unit(fontsize*0.9, "pt"),
    legend.key.spacing.y = unit(fontsize*0.2, "pt"),
    plot.margin = margin(l = -10),
    # legend.spacing.x = unit(fontsize*2, "pt")
    # legend.spacing.y = unit(fontsize*0.5, "pt")
  ); plot.map.bh

# plot.map <- plot.map.br + wrap_elements(plot.map.rs) +
#   plot_layout(ncol = 2, guides = "collect") &
#   theme(legend.position = "bottom",
#         legend.margin = margin(t = -5, "mm"),
#         legend.spacing.x = unit(-1, "mm")); plot.map

plot.map <- 
  wrap_elements(plot.map.br) | (wrap_elements(plot.map.rs) / wrap_elements(plot.map.bh)) +
  plot_layout(guides = "collect") &
  theme(
  ); plot.map

  
ggsave("figuras/preliminares/mapa_estacoes.png", plot.map, width = 15.5, height = 12, units = "cm", dpi = 300)
invisible(gc())


# Gráfico separado apresentando os IDs das estações subdiárias
# TALVEZ SEJA INTERESSANTE ADICIONAR CAMADAS DAS REGIÕES HIDROGEOMOEFOLÓGICAS DO SGB
subdaily.sf <- subdaily.sf |> 
  mutate(id = gauge.ids$id[match(gauge_code, gauge.ids$gauge_code)]) |> 
  arrange(id)

map.subdaily.id.plot <- ggplot() +
  # geom_sf(data = rs.sf, fill = NA, color = "black", linewidth = 0.3) +
  geom_sf(data = macrobh.sf, color = "black", fill = NA, linewidth = 0.3) +
  geom_sf(data = daily.sf, color = "#ff0000", alpha = 0.2, size = 0.5) +
  geom_sf(data = subdaily.sf, color = "#1100ff", alpha = 0.8, size = 0.5) +
  geom_text_repel(
    aes(label = id, geometry = geometry), subdaily.sf,
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
  ); map.subdaily.id.plot

ggsave("figuras/preliminares/mapa_estacoes_id.png", map.subdaily.id.plot, width = 12, height = 12, units = "cm", dpi = 300)

# map.gauge.id.main <- ggplot() +
#   geom_sf(data = rs.sf, fill = NA, color = "black", linewidth = 0.3) +
#   geom_sf(data = subdaily.sf, color = "#1100ff", alpha = 0.8, size = 0.5) +
#   geom_text_repel(
#     aes(label = id, geometry = geometry), subdaily.sf,
#     stat = "sf_coordinates", family = font, size = textsize,
#     color = "black", bg.color = "white", bg.r = 0.15, max.overlaps = Inf,
#     min.segment.length = 0, segment.color = "grey50", segment.size = 0.2,
#     box.padding = 0.2, point.padding = 0.1
#   ) +
#   coord_sf(crs = st_crs(4674), xlim = bbox.rs[c(1,3)], ylim = bbox.rs[c(2,4)], expand = TRUE); map.gauge.id.main

# urban.bbox <- subdaily.sf[subdaily.sf$gauge_code == "SBPA",] |> st_bbox() |> 
#   st_as_sfc() |> 
#   st_transform(crs = 31982) |> 
#   st_buffer(dist = 1e4) |> 
#   st_transform(crs = st_crs(subdaily.sf)) |> 
#   st_bbox() |> print()


# map.gauge.id.inset <- ggplot() +
#   geom_sf(data = rs.sf, fill = NA, color = "black", linewidth = 0.3) +
#   geom_sf(data = subdaily.sf, color = "#1100ff", alpha = 0.8, size = 0.5) +
#   geom_text(
#     aes(label = id, geometry = geometry), subdaily.sf,
#     stat = "sf_coordinates", family = font, size = textsize,
#     hjust = -1, vjust = -1
#   ) +
#   coord_sf(crs = st_crs(4674), xlim = urban.bbox[c(1,3)], ylim = urban.bbox[c(2,4)], expand = TRUE); map.gauge.id.inset

# # map.gauge.id.plot <- 
#   map.gauge.id.main + 
#   inset_element(
#     map.gauge.id.inset, left = 0.05, bottom = 0.05, right = 0.4, top = 0.4
#   )



## TABELA INFORMAÇÕES ESTAÇÕES -------------------------------------------

subdaily.info.table <- info.df |> 
  mutate(
    id = gauge.ids$id[match(gauge_code, gauge.ids$gauge_code)],
  ) |> 
  left_join(
    imax.df |> 
      group_by(gauge_code) |> 
      reframe(
        last_year = max(wateryear),
        n_years = n_distinct(wateryear),
        rain_yr_mean = mean(rain_yr, na.rm = TRUE)
      ),
      by = "gauge_code"
  ) |> 
  arrange(id) |> 
  mutate(
    time_step = as.integer(time_step),
    long_lat = sprintf("[%0.3f, %0.3f]", long, lat)
  ) |> 
  select(gauge_code, id, responsible, time_step, n_years,  rain_yr_mean, elevation, long_lat)

table.units <- c("-", "-", "-", "min", "anos", "mm", "m", "º")
table.names <- c("Código estação", "ID", "Operador", "$\\Delta$", "N", "PAA$^{*}$", "Elevação", "Coordenadas$^{**}$")
names(subdaily.info.table) <- table.names

n.col <- ncol(subdaily.info.table)
subdaily.info.xtable <- subdaily.info.table |> 
  xtable(
    caption = "Informações básicas das estações pluviográficas do Rio Grande do Sul utilizadas.",
    label = "tab_subdaily_info",
    align = rep("c", n.col + 1)
  )

extra.rows <- list()
extra.rows$pos <- list(0, nrow(subdaily.info.table)) # posições das novas linhas
ltx.cmd.units <- sprintf(
    "%s \\\\
    \\midrule
    \\endfirsthead
    \\caption[]{Informações básicas das estações pluviográficas do Rio Grande do Sul utilizadas. (continuação)} \\\\
    \\toprule
    %s \\\\
    %s \\\\
     \\midrule
    \\endhead
    \\midrule
    \\multicolumn{%d}{r}{{Continua na próxima página...}} \\\\
    \\endfoot
    \\endlastfoot \n",
    paste(table.units, collapse = " & "),
    paste(table.names, collapse = " & "),
    paste(table.units, collapse = " & "),
    n.col
  )
ltx.cmd.foootnote <- sprintf(
  "\\bottomrule
  \\multicolumn{%d}{l}{$^*$EPSG 4674.} \\\\
  \\multicolumn{%d}{l}{$^{**}$Precipitação acumulada anual média entre anos com percentual de falhas inferior a 20\\%%.} \\\\ \n",
  n.col, n.col
)

extra.rows$command <- c(ltx.cmd.units, ltx.cmd.foootnote)

subdaily.info.string <- print(
  subdaily.info.xtable, type = "latex", include.rownames = FALSE, caption.placement = "top",
  sanitize.colnames.function = identity, 
  sanitize.text.function = identity,
  tabular.environment = "longtable", # changes to longtable
  floating = FALSE,                  # required for longtable to work
  add.to.row = extra.rows,
  hline.after = c(-1)
)

cat(
  "\\begingroup\n\\footnotesize\n",
  subdaily.info.string,
  "\\endgroup\n",
  file = "tabelas/5_subdaily_info_table.txt",
  sep = ""
)


## ANÁLISE EXPLORATÓRIA --------------------------------------------------

# Estimativas por duração com L-momentos
# Razões de L-momentos
ratios.df <- imax.df |> 
  group_by(gauge_code, d) |> 
  arrange(desc(imax)) |> 
  reframe(
    slmom = list(lmom::samlmu(imax)),
    lcv = unlist(slmom)[[2]]/unlist(slmom)[[1]],
    lsk = unlist(slmom)[[3]],
    lkt = unlist(slmom)[[4]] # calcular a 4ª razão precisa de pelo menos 4 observações
  ) |> 
  select(!slmom) |> 
  pivot_longer(
    cols = starts_with("l"),
    names_to = "ratios",
    values_to = "value"
  ) |> 
  mutate(
    label = case_when(
      ratios == "lcv" ~ "L-CV",
      ratios == "lsk" ~ "L-assimetria",
      ratios == "lkt" ~ "L-curtose"
    ),
    label = factor(label, levels = c("L-CV", "L-assimetria", "L-curtose")),
    n_years = imax.df$n_years[match(gauge_code, imax.df$gauge_code)],
  ) |> 
  arrange(gauge_code, d, ratios) 


# Identificar outliers
is.outlier <- function(x) x < quantile(x, 0.25) - 1.5*IQR(x) | x > quantile(x, 0.75) + 1.5*IQR(x)
outlier.gauges <- ratios.df |> 
  group_by(d, ratios) |> 
  mutate(outlier = is.outlier(value)) |> 
  filter(outlier) |> 
  pull(gauge_code) |> 
  unique() |> print()

paste(outlier.gauges, collapse = ", ")

# Criar uma tabela que so inclui os grupos que contiverem 0
#  p/ evitar plotar geom_hline desnecessariamente
zero.lines <- ratios.df |> 
  group_by(label) |> 
  reframe(has_zero = min(value) <= 0 & max(value) >= 0) |> 
  filter(has_zero) |> 
  mutate(y = 0) |> print()

# Boxplot
ratios.d.plot <- ratios.df |> 
  ggplot(aes(x = factor(d), y = value)) +
  facet_wrap(~label, ncol = 1, scale = "free_y", strip.position = "left") +
  geom_hline(aes(yintercept = y), zero.lines, linetype = "dashed", linewidth = 0.2, color = "grey") +
  geom_boxplot(
    fill = "grey90", median.linewidth = 0, staplewidth = 0.5, linewidth = 0.2,
    whisker.linetype = "dashed", outlier.size = 0.5, outlier.alpha = 0.5
  ) +
  stat_summary(
    geom = "crossbar", width = 0.4, linewidth = 0.1, fun.data = function(x){
      c(y = median(x), ymin = median(x), ymax = median(x))
    }
  ) +
  labs(x = "Duração [h]", y = "Razões de L-momentos") +
  theme(
    axis.title.y = element_blank(),
    strip.placement = "outside",
    strip.text.y = element_text(hjust = 0.5, size = fontsize)
  ); ratios.d.plot

# Importante saber quem são as estações "outliers" em cada uma das razões
# Não precisa incluir na figura, só adicionar no texto

ggsave("figuras/preliminares/ratios_d_boxplot.png", ratios.d.plot, height = 10, width = 14, units = "cm", dpi = 300)


# Investigação L-curtose vs. tamanho das séries
# yr.breaks <- quantile(
#   as.integer(ratios.df$n_years), 
#   probs = seq(0, 1, length.out = 3), 
#   type = 1
# ) |> unique() |> print()

# yr.breaks <- c(8, 11, 14, 19)

# # lkt.nyear.plot <- 
# ratios.df |> 
#   filter(
#     d %in% durations,
#     ratios == "lkt"
#   ) |> 
#   mutate(
#     yr_interval = cut(n_years, breaks = yr.breaks, include.lowest = TRUE)
#   ) |> 
#   ggplot(aes(x = yr_interval, y = value)) +
#   facet_wrap(~factor(d)) +
#   geom_boxplot()


# TABELA ESTATÍSTICAS BÁSICAS --------------------------------------------

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
    "Mínimo" = min(x, na.rm = TRUE),
    "$Q_5$" = quantile(x, 0.05, na.rm = TRUE),
    "$Q_{25}$" = quantile(x, 0.25, na.rm = TRUE),
    "$Q_{50}$" = median(x, na.rm = TRUE),
    "$Q_{75}$" = quantile(x, 0.75, na.rm = TRUE),
    "$Q_{95}$" = quantile(x, 0.95, na.rm = TRUE),
    "Máximo" = max(x, na.rm = TRUE),
    "Média" = mean(x, na.rm = TRUE),
    "Assimetria" = moments::skewness(x, na.rm = TRUE)
  )
}

table.duration <- c(1, 6, 12, 24)

# Estatísticas básicas
# Aplicar resumo por estação e duração
stats.df <- imax.df |> 
  filter(d %in% table.duration) |> 
  group_by(gauge_code, d) |> 
  group_modify(~as.data.frame(t(get.stats(.x$imax)))) |> # group_modify requer que o retorno seja um data.frame
  ungroup() |> 
  select(-p0)

# Aplicar resumo por estatísticas
stats.table <- stats.df |> 
  pivot_longer(cols = -c(gauge_code, d), names_to = "statistic", values_to = "value") |> 
  group_by(d, statistic) |> 
  summarise(across(value, ~list(summary.stats(.x))), .groups = "drop") |> # dá pra usar o group_modify tambem
  unnest_wider(value) |> 
  pivot_longer(cols = -c(d, statistic), names_to = "summary_stat", values_to = "value") |> 
  pivot_wider(id_cols = c(d, summary_stat), names_from = statistic, values_from = value) |> 
  relocate(d, .after = last_col()) |> 
  mutate(summary_stat = gsub("\\$[^\\$]*$", "$", summary_stat)) |> 
  select(summary_stat, n_years, mean, median, sd, t3, t4, d) |> 
  rename(
    " " = summary_stat,
    "$N$" = n_years,
    "Média" = mean,
    "Mediana" = median,
    "$SD$" = sd,
    "$\\hat\\tau_3$" = t3,
    "$\\hat\\tau_4$" = t4
  )
  
ncol.d <- n_distinct(stats.table$` `)
nrow.d <- n_distinct(stats.table$` `)

table.ls <- as.list(seq(0,  nrow(stats.table) - 1, by = ncol.d))
ltx.cmd <- sapply(unique(stats.table$d), function(dur){
  sprintf("\\hline \\multicolumn{%d}{c}{\\textbf{Duração: %g h}} \\\\ \\hline\n", ncol(stats.table) - 1, dur)
})

stats.xtable <- stats.table |> 
  # relocate(d, .after = 1) |> 
  select(-d) |> 
  xtable::xtable(
    caption = "Estatísticas básicas do conjunto de séries de intensidades máximas anuais para durações 1, 6, 12 e 24 h.",
    label = "tab_estatistias_subdiarias",
    align = c("l", "l", rep("c", ncol(stats.table) - 2)),
    digits = 3
  )

stats.string <- print(
  stats.xtable,
  type = "latex",
  include.rownames = FALSE,
  caption.placement = "top",
  sanitize.colnames.function = identity,
  sanitize.text.function = identity,
  hline.after = c(-1, nrow(stats.table)),
  add.to.row = list(pos = table.ls, command = ltx.cmd)
)

cat(
  "\\begingroup\n\\footnotesize\n",
  stats.string,
  "\\endgroup\n",
  file = "tabelas/5_estatisticas_basicas_subdiarios.txt",
  sep = ""
)


### Informações complementares -------------------------------------------

# Número de estações c/ N <= 8 e 9 anos
stats.df |> slice_min(d) |> count(n_years == 8) |> mutate(prop = n/(sum(n)))
stats.df |> slice_min(d) |> count(n_years <= 9) |> mutate(prop = n/(sum(n)))
stats.df |> slice_min(d) |> count(n_years == 13) |> mutate(prop = n/(sum(n)))

# Número de estações com theta == 0
par.dgev.df |> filter(model == "MLE", parname == "theta", par == 0) |> nrow() # 20 estações
par.dgev.df |> filter(model == "RGMLE – B(0.0178, 0.0627)", parname == "theta", par == 0) |> nrow() # 18 estações

# Número de estações com theta < 0.1
par.dgev.df |> filter(model == "MLE", parname == "theta", par < 0.1) |> nrow() # 24 estações
par.dgev.df |> filter(model == "RGMLE – B(0.0178, 0.0627)", parname == "theta", par < 0.1) |> nrow() # 23 estações



## DIAGRAMA DE L-MOMENTOS POR TAMANHO DAS SÉRIES -------------------------

table.duration <- c(1, 6, 12, 24)
n.breaks <- 3
n.years <- info.df |>
  mutate(n_years = imax.df$n_years[match(gauge_code, imax.df$gauge_code)]) |> 
  pull(n_years) |> 
  sort(decreasing = TRUE)

probs <- seq(0, 1, length.out = n.breaks + 1)
# breaks <- quantile(n.years, probs, type = 1, names = FALSE) # type = 1 garante valores inteiros
breaks <- c(8,9,12,19)
yr.interval <- cut(n.years, breaks, right = FALSE, include.lowest = TRUE)

# Razões amostrais por duração
ratios.sam <- imax.df |> 
  filter(d %in% table.duration) |> 
  group_by(gauge_code, d) |> 
  summarise(
    origem = "Amostrais",
    t3 = lmom::samlmu(imax)[[3]],
    t4 = lmom::samlmu(imax)[[4]],
    n_years = first(n_years),
    record = cut(
      n_years, breaks = breaks, 
      right = FALSE, include.lowest = TRUE
    ),
    .groups = "drop"
  ) |> 
  arrange(n_years)

# Média amostral por duração
ratios.av <- tibble(
  origem = "Média",
  t3 = mean(ratios.sam$t3, na.rm = TRUE),
  t4 = mean(ratios.sam$t4, na.rm = TRUE)
)

# Razões teóricas --- Gumbel
ratios.gum <- tibble(
  origem = "Gumbel",
  t3 = 0.1699,
  t4 = 0.1504
)

# Razões teóricas --- GEV
ratios.gev <- tibble(kappa = seq(-0.6, 0.6, 0.001)) |> 
  mutate(
    origem = "GEV",
    t3 = 2*(1 - 3^(-kappa))/(1 - 2^(-kappa)) - 3,
    t4 = (1 - 6*2^(-kappa) + 10*3^(-kappa) - 5*4^(-kappa))/(1 - 2^(-kappa))
  )

# Cores
yr.col <- c("#B4EEB4", "#698B69", "#191970")
yr.col <- setNames(yr.col, levels(yr.interval)) |> print()
range.kt <- range(ratios.sam$t4, na.rm = TRUE)

d.label <- tibble(
  d = table.duration,
  d_label = sprintf("Duração: %2g h", d)
)

# Plotar diagrama --- cores por tamanho das séries
gu.col <- "orange2"
av.col <- "red"
# ratio.diagram.plot <- 
ratio.diagram.plot <- ggplot(ratios.sam, aes(x = t3, y = t4)) +
  facet_wrap(~d, nrow = 2) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey80") +
  geom_hline(yintercept = ratios.gum$t4, linewidth = 0.2, color = "orange2", linetype = "dashed") +
  geom_line(aes(linetype = "GEV"), ratios.gev, linewidth = 0.5, color = "grey10") +
  geom_point(aes(color = record), alpha = 0.8, size = 0.8, key_glyph = "rect") +
  geom_point(aes(shape = "Média"), ratios.av, color = av.col, size = 2) +
  geom_point(aes(shape = "Gumbel"), ratios.gum, color = gu.col, size = 2) +
  geom_text(aes(x = Inf, y = Inf, label = d_label), d.label, family = font, size = textsize*1.1, hjust = 1.1, vjust = 2) +
  scale_color_manual(name = "N [anos]", values = yr.col) +
  scale_shape_manual(name = NULL, values = c(17, 16)) +
  scale_linetype_manual(name = NULL, values = "solid") +
  coord_cartesian(ylim = range.kt) +
  labs(x = "L-assimetria", y = "L-curtose") +
  theme(
    strip.text = element_blank(),
    legend.position = c(0, 0.99), legend.justification = c(0,1),
    legend.title = element_text(size = fontsize*0.8),
    legend.spacing.y = unit(-2, "pt")
  ); ratio.diagram.plot

# Boxplot L-curtose por faixa de tamanho de série e duração
ratios.boxplot.nyr <- ratios.sam |> 
  ggplot(aes(x = factor(d), y = t4, fill = record)) +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey80") +
  geom_hline(yintercept = ratios.gum$t4, linewidth = 0.2, color = "orange2", linetype = "dashed") +
  geom_boxplot(
    color = "black", median.linewidth = 0, staplewidth = 0.5, linewidth = 0.2,
    whisker.linetype = "dashed", outlier.size = 0.5, outlier.alpha = 0.5,
    show.legend = FALSE
  ) +
  stat_summary(
    geom = "crossbar", width = 0.4, linewidth = 0.1, fun.data = function(x){
      c(y = median(x), ymin = median(x), ymax = median(x))
    }
  ) +
  scale_fill_manual(name = NULL, values = yr.col) +
  coord_cartesian(ylim = range.kt) +
  labs(y = "L-curtose") +
  theme(
    axis.text.x = element_blank(),
    axis.title.x = element_blank(),
    legend.position = "none",
    plot.margin = margin(l = 5.5)
  ); ratios.boxplot.nyr

# Juntar todo mundo
ratios.plot <- 
  wrap_plots(ratio.diagram.plot + labs(subtitle = "(a)")) /
  wrap_plots(ratios.boxplot.nyr + labs(subtitle = "(b)")) +
  plot_layout(heights = c(1, 0.5)); ratios.plot

ggsave(
  "figuras/preliminares/subdiarios/diagrama_lmomentos.png", ratios.plot, 
  width = 15, height = 15, units = "cm", dpi = 300
)


### Informações complementares -------------------------------------------

# Número de estações c/ t4 < 0.1504 (Gumbel)
ratios.sam |> 
  group_by(d) |> 
  count(t4 < ratios.gum$t4) |> 
  mutate(prop = n/sum(n))


## DIAGRAMA DE L-MOMENTOS POR DURAÇÃO ------------------------------------

# Plotar o diagrama, mas somente as cores das "durações de referência" em `durations`
durations

# Razões amostrais por grupo de duração
ratios.sample <- ratios.df |> 
  pivot_wider(
    id_cols = c(gauge_code, d, n_years),
    names_from = ratios,
    values_from = value
  ) |> 
  mutate(
    d_interval = cut(d, breaks = durations, right = FALSE, include.lowest = TRUE),
    which_ratio = "Amostral"
  )

d.intervals <- unique(ratios.sample$d_interval)

# Quantas estações possuem lkt > 0.4?
ratios.sample |> 
  # group_by(d) |> 
  count(lkt > 0.4) |> 
  mutate(prop = n/sum(n)) |> 
  print(n = 48)

# Comparar razões médias por razões médias por duração
# ratio.average <- ratios.sample |> 
#   group_by(d_interval) |> 
#   reframe(
#     lsk = mean(lsk, na.rm = TRUE),
#     lkt = mean(lkt, na.rm = TRUE),
#     which_ratio = "Média"
#   )

ratio.average <- ratios.sample |> 
  reframe(
    lsk = mean(lsk, na.rm = TRUE),
    lkt = mean(lkt, na.rm = TRUE),
    which_ratio = "Mèdia"
  )

ratio.gumbel <- data.frame(which_ratio = "Gumbel", lsk = 0.1699, lkt = 0.1504)

ratio.gev <- data.frame(kappa = seq(-0.65, 0.8, 0.001)) %>% 
  mutate(
    which_ratio = "GEV",
    lsk =  2*(1 - 3^(-kappa))/(1 - 2^(-kappa)) - 3,
    lkt = (1 - 6*2^(-kappa) + 10*3^(-kappa) - 5*4^(-kappa))/(1 - 2^(-kappa))
  )

# Plot sample l-ratios by intervals of AMS record length
# d.colors <- setNames(c("darkseagreen2", "darkseagreen4", "#5480c8", "midnightblue")[-1], d.intervals)

gu.color <- "grey10"
av.color <- "red"

ratios.diagram.plot <- ggplot(mapping = aes(x = lsk, y = lkt)) +
  geom_point(data = ratios.sample, aes(color = d_interval), alpha = 0.6, size = 0.8) +
  geom_line(data = ratio.gev, aes(linetype = "GEV"), linewidth = 0.8, color = "grey10") +
  geom_point(data = ratio.gumbel, aes(shape = "Gumbel"), color = gu.color, size = 2) +
  geom_point(data = ratio.average, aes(shape = "Average"), color = av.color, size = 2) +
  scale_color_manual(
    name = "Durações [h]", values = setNames(d.colors[-4], d.intervals),
    guide = guide_colorsteps(
      direction = "horizontal", order = 3, label.hjust = 0.5, label.vjust = 0,
      theme = theme(title = element_text(vjust = 0.95, size = fontsize*0.8)), show.limits = TRUE
    )
  ) +
  scale_linetype_manual(name = "", values = c("GEV" = 1), guide = guide_legend(order = 2)) +
  scale_shape_manual(name = "", values = c("Gumbel" = 17, "Average" = 15), guide = guide_legend(order = 1, nrow = 1)) +
  labs(x = "L-assimetria", y = "L-curtose") +
  theme(
    legend.position = c(0, 0.985),
    legend.justification = c(0, 1),
    legend.title.position = "right",
    legend.spacing.y = unit(-fontsize, "pt"),
    legend.text = element_text(vjust = 0.5),
    legend.key.height = unit(fontsize, "pt")
    # legend.title = element_blank()
  ); ratios.diagram.plot

ggsave(
  "figuras/preliminares/subdiarios/ratios_diagram.png", ratios.diagram.plot, 
  width = 15, height = 8, units = "cm", dpi = 300
)

## ESTATÍSTICAS ----------------------------------------------------------

# Calcular estatísticas básicas p/ durações de referência (1, 6, 12, 24)
  

## INCERTEZA EM KAPPA ----------------------------------------------------

select.gg <- c("A801", "A831", "432300204A")
pos.width <- position_dodge(0.6)

col <- c("#648FFF", "#dc267f", "#e1a3aa", "#FE6100")
prior.col <- setNames(col, prior.names)

# Visualizar intervalos de confiança de kappa POR DURAÇÃO com estimativas
# da GEV com L-momentos --- figura c/ todas as estações
shape.boot.plot <- shape.boot.df |> 
  mutate(gauge_label = gauge.labels[gauge_code]) |> 
  select(!c(model, method, id)) |> 
  ggplot(aes(x = factor(d), y = kappa)) +
  facet_wrap(~factor(gauge_label, levels = order.n), ncol = 4, scales = "free_y") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, color = "grey60") +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.2) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper), 
    width = 0, linewidth = 0.3, alpha = 0.7, color = "#2baa68", position = pos.width) +
  geom_point(size = 0.15, position = pos.width, color = "#2baa68") +
  scale_x_discrete(breaks = c(1, 6, 12, 18, 24)) +
  labs(y = "κ") +
  theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.margin = margin(5, 5, 1, 5),
    legend.key.size = unit(3, "mm"),
    text = element_text(size = fontsize*0.8),
    axis.title.x = element_blank(),
    # axis.text.x = element_text(angle = 90, hjust = 1),
    panel.border = element_rect(linewidth = 0.3, color = "black"),
    panel.spacing.y = unit(2, "pt"),
    strip.text = element_text(hjust = 0, size = fontsize*0.7, margin = margin(t = 3, b = 2))
  ); shape.boot.plot

ggsave("figuras/preliminares/subdiarios/ci_kappa_gev_lmom_boot.png", shape.boot.plot, height = 22, width = 16, units = "cm", dpi = 300)


# Visualizar intervalos de confiança de kappa POR DURAÇÃO com estimativas
# da GEV com L-momentos --- figura c/ estações em `select.gg`
shape.boot.plot <- shape.boot.df |> 
  filter(gauge_code %in% select.gg) |> 
  mutate(gauge_label = gauge.labels[gauge_code]) |> 
  select(!c(model, method, id)) |> 
  ggplot(aes(x = factor(d), y = kappa)) +
  facet_wrap(~factor(gauge_label, levels = order.n), nrow = 1) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, color = "grey60") +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.2) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper), 
    width = 0, linewidth = 0.3, alpha = 0.7, color = "#2baa68", position = pos.width) +
  geom_point(size = 0.15, position = pos.width, color = "#2baa68") +
  scale_x_discrete(breaks = c(1, 6, 12, 18, 24)) +
  labs(x = "Duração [h]", y = "κ") +
  theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    panel.border = element_rect(linewidth = 0.3, color = "black"),
    strip.text = element_text(hjust = 0)
  ); shape.boot.plot

ggsave(
  "figuras/preliminares/subdiarios/ci_kappa_gev_lmom_boot_selected.png", shape.boot.plot, 
  height = 7, width = 15, units = "cm", dpi = 300
)


# Visualizar intervalos de confiança p/ kappa com POR DURAÇÃO
# com estimativas com perfil de verossimilhança + priori GMLE (media 0 e -0.1) e RGMLE
# p/ estações selecionadas --- `select.gg`
shape.prof.plot <- par.gev.df |> 
  filter(gauge_code == select.gg[1]) |> 
  mutate(
    gauge_label = gauge.labels[gauge_code],
    prior_label = factor(prior.names[prior], prior.names)
  ) |> 
  filter(
    parname == "kappa",
    prior_label %in% prior.names[2:4] # selecionar GMLE - B(0,0.122); GMLE - B(-0.1,0.122); RGMLE - B(0.0178,0.0627)
  ) |> 
  ggplot(aes(x = factor(d), y = par, color = prior_label)) +
  facet_wrap(~factor(gauge_label, order.n), ncol = 4) +
  # geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, color = "grey60") +
  geom_hline(yintercept = 0, color = "grey80", linewidth = 0.2) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper), 
    width = 0, linewidth = 0.4, alpha = 0.7, position = position_dodge(0.55)
  ) +
  geom_point(size = 0.3, position = position_dodge(0.55)) +
  scale_x_discrete(breaks = c(1, 6, 12, 18, 24)) +
  scale_color_manual(name = NULL, values = prior.col) +
  labs(x = "Durações [h]", y = "κ") +
  theme(
    legend.position = "bottom",
    legend.margin = margin(0, 5, 1, 5),
    legend.key.size = unit(3, "mm"),
    # text = element_text(size = fontsize*0.8),
    panel.border = element_rect(linewidth = 0.3, color = "black"),
    strip.text = element_text(hjust = 0)
  ); shape.prof.plot

ggsave(
  "figuras/preliminares/subdiarios/ci_kappa_gev_gmle_rgmle_selected.png", shape.prof.plot,
  width = 15, height = 8, units = "cm", dpi = 300
)


# Visualizar intervalos de confiança p/ kappa com POR DURAÇÃO
# com estimativas com perfil de verossimilhança + priori GMLE (media 0 e -0.1) e RGMLE
mod <- c("gmle_martins", "gmle", "rgmle")
walk(1:3, function(i){

  par.gev.df |> 
    mutate(gauge_label = gauge.labels[gauge_code], prior_label = factor(prior.names[prior], prior.names)) |> 
    filter(parname == "kappa", prior_label == prior.names[i + 1]) |> 
    ggplot(aes(x = factor(d), y = par, color = prior_label)) +
    facet_wrap(~factor(gauge_label, order.n), ncol = 4) +
    geom_hline(yintercept = 0, color = "grey80", linewidth = 0.2) +
    geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, linewidth = 0.2, alpha = 0.7, position = pos.width) +
    geom_point(size = 0.15, position = position_dodge(0.55)) +
    scale_x_discrete(breaks = c(1, 6, 12, 18, 24)) +
    scale_color_manual(name = NULL, values = prior.col) +
    labs(x = "Durações [h]", y = "κ") +
    theme(
      legend.position = "none",
      axis.title.x = element_blank(),
      text = element_text(size = fontsize*0.8),
      panel.border = element_rect(linewidth = 0.3, color = "black"),
      panel.spacing.y = unit(2, "pt"),
      strip.text = element_text(hjust = 0, size = fontsize*0.7, margin = margin(t = 3, b = 2))
    )

  message(mod[i])
  ggsave(
    sprintf("figuras/preliminares/subdiarios/ci_kappa_gev_prof_%s.png", mod[i]), 
    width = 16, height = 22, units = "cm", dpi = 300
  )

})



# Em que momento na estação A831 `select.gg[2]` kappa deixa de ser negativo
par.gev.df |> 
  filter(
    gauge_code == select.gg[2],
    parname == "kappa",
    prior == "c(0, 0.122)"
  ) |> 
  select(d, par)

# Visualizar intervalos de confiança para kappa com estimativas
# do modelo de invariância de escala com perfil de verossimilhança

pos.width <- position_dodge(0.8)
par.dgev.df |> 
  filter(parname == "kappa") |> 
  left_join(
    imax.df |> 
      slice_head(n = 1, by = "gauge_code") |> 
      select(gauge_code, n_years, gauge_label),
    by = "gauge_code"
  ) |> 
  mutate(
    id = gauge.ids$id[match(gauge_code, gauge.ids$gauge_code)],
    prior_label = factor(prior.names[prior], prior.names)
  ) |>
  filter(prior_label != prior.names[3]) |> # remover GMLE - B(0,0.122)
  # ggplot(aes(x = gauge_label, y = par)) +
  # ggplot(aes(x = fct_reorder(gauge_code, n_years), y = par, color = factor(model, which.model[c(1,3,2)]))) +
  ggplot(aes(x = factor(id), y = par, color = prior_label)) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, color = "grey60") +
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.2, color = "grey60") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, linewidth = 0.3, alpha = 0.7, position = pos.width) +
  geom_point(size = 0.15, position = pos.width) +
  scale_color_manual(name = "", values = prior.col) +
  labs(x = "ID", y = "κ") +
  theme(
    legend.position = c(0, 0),
    legend.justification = c(0, 0),
    legend.margin = margin(5, 5, 1, 5),
    legend.key.size = unit(3, "mm"),
    axis.title.x = element_text(margin = margin(t = 5)),
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.border = element_rect(linewidth = 0.3, color = "black"),
    strip.text = element_text(hjust = 0)
  )

# ggsave("figuras/preliminares/ci_H_dgev_proflik2.png", height = 8, width = 14, units = "cm", dpi = 300)
ggsave("figuras/preliminares/subdiarios/ci_kappa_dgev_proflik.png", height = 8, width = 15, units = "cm", dpi = 300)

# Estações p/ verificar
par.dgev.df |> 
  filter(parname == "kappa") |> 
  select(gauge_code, par, ci_lower, ci_upper, prior) |> 
  filter(abs(par) < 1e-5 | abs(par) > 0.4) |> 
  unique() |> 
  print()


# Informações adicionais
# Quantas estações com RGMLE são positivas e quantas são negativas
par.dgev.df |> 
  filter(
    parname == "kappa",
    prior == which.prior[2] # manter estimativas c/ RGMLE
  ) |> 
  select(gauge_code, par, ci_lower, ci_upper) |> 
  mutate(
    sign = sign(par),                                 # sinal de kappa
    ci_includes_zero = 0 <= ci_upper & 0 >= ci_lower, # IC inclui 0?
  ) |> 
  group_by(sign) |> 
  summarise(
    min = min(par), max = max(par),
    n = n()
  )

# Quantas estações c/ RGMLE ajustaram kappas entre -0,05 e 0,05
par.dgev.df |> 
  filter(
    parname == "kappa",
    prior == which.prior[2], # manter estimativas c/ RGMLE
    par <= 0.05 & par >= -0.05
  ) |> 
  select(gauge_code, par, ci_lower, ci_upper) |>
  left_join(
    imax.df |> select(gauge_code, n_years) |> slice(1, .by = "gauge_code"),
    by = "gauge_code"
  )


# Análises c/ intervalos de confiança
ci.width <- par.dgev.df |> 
  filter(
    parname == "kappa",
    model %in% prior.names[c(1, 4)] # manter somente MLE e RGMLE
  ) |> 
  select(gauge_code, par, ci_lower, ci_upper, model) |> 
  mutate(
    model = case_when(
      model == prior.names[1] ~ "mle",
      model == prior.names[4] ~ "rgmle",
    ),
    ci_width = ci_upper - ci_lower
  ) |> 
  pivot_wider(
    id_cols = gauge_code,
    names_from = model,
    values_from = c(par, ci_width)
  ) |> 
  mutate(
    smaller_ci_w_rgmle = ci_width_mle > ci_width_rgmle
  )

# Quantas estações tiveram ICs menores c/ RGMLE
ci.width |> count(smaller_ci_w_rgmle)

# Quais estações tiveram ICs maiores c/ RGMLE
ci.width |> 
  filter(!smaller_ci_w_rgmle) |> 
  mutate(id = gauge.ids$id[match(gauge_code, gauge.ids$gauge_code)], .after = gauge_code)

# Em quantos porcento o intervalo c/ RGMLE é menor
ci.width |> 
  mutate(width_ratio = (ci_width_rgmle - ci_width_mle)/ci_width_mle*100) |> 
  select(gauge_code, width_ratio) |> 
  arrange(width_ratio) |> 
  print(n = 54)

## PLOTAR ALGUNS PERFIS DE VEROSSIMILHANÇA -------------------------------

# `parname == kappa`
gg <- "A832"
which.par <- "kappa"
# gg <- "430060402A"
proflik.kappa <- fitprof.dgev.wprior.ls$prof.par |> 
# proflik.kappa <- fitprof.nprior.ls$prof.par |> 
  filter(
    parname == which.par, 
    gauge_code == gg
  ) |> 
  select(par, loglik) |> 
  rename(kappa = par) |> 
  filter(loglik < -180)

lim <- fitprof.dgev.wprior.ls$par |> 
  filter(
    gauge_code == gg, 
    parname == which.par
  ) |> 
  pull(lim)

conf <- fitprof.dgev.wprior.ls$par |> 
  filter(
    gauge_code == gg, 
    parname == which.par
  ) |> 
  select(ci_lower, ci_upper) |> 
  unlist()
  
ggplot(proflik.kappa, aes(x = kappa, y = loglik)) +
  geom_line() +
  geom_vline(xintercept = conf, linewidth = 0.3, linetype = "solid") +
  geom_point(shape = 21, fill = "white") +
  geom_hline(yintercept = lim, linewidth = 0.3, linetype = "dashed") +
  annotate(geom = "text", label = round(lim,1), x = -Inf, y = lim, family = font, size = textsize*1.1, hjust = 1.19) +
  labs(x = "κ", y = "l<sub>p</sub>(κ)") +
  coord_cartesian(clip = "off") +
  theme(
    axis.title.y = ggtext::element_markdown(),
    panel.border = element_rect(linewidth = 0.3)
  )

ggsave("figuras/preliminares/subdiarios/perfil_kappa_proflik.png", width = 12, height = 6, units = "cm", dpi = 300)

par.df |> filter(gauge_code == gg, parname == "kappa", prior != "NULL")

# quantil
# proflik.q <- fitprof.wprior.ls$prof.q |> 
gg <- "A801"
proflik.q <- fitprof.dgev.wprior.ls$prof.q |> 
  filter(
    gauge_code == gg,
    d == durations[1], 
    return_period == 100
  ) |> 
  select(q, loglik)

lim <- fitprof.dgev.wprior.ls$idf |> 
  filter(
    gauge_code == gg, 
    d == durations[1], 
    return_period == 100
  ) |> 
  pull(lim)

conf <- fitprof.dgev.wprior.ls$idf |> 
  filter(
    gauge_code == gg, 
    d == durations[1], 
    return_period == 100
  ) |> 
  select(ci_lower, ci_upper) |> 
  unlist()

nrow(proflik.q)

proflik.q <- proflik.q[seq_len(nrow(proflik.q)) %% 2 != 0,]
proflik.q <- proflik.q[seq_len(nrow(proflik.q)) %% 2 != 0,]

nrow(proflik.q)

ggplot(proflik.q, aes(x = q, y = loglik)) +
  geom_line() +
  geom_vline(xintercept = conf, linewidth = 0.3, linetype = "solid") +
  geom_point(shape = 21, fill = "white") +
  geom_hline(yintercept = lim, linewidth = 0.3, linetype = "dashed") +
  annotate(geom = "text", label = round(lim,1), x = -Inf, y = lim, family = font, size = textsize*1.1, hjust = 1.19) +
  coord_cartesian(clip = "off") +
  labs(x = "i<sub>d,p</sub> [mm/h]", y = "l<sub>p</sub>(i<sub>d,p</sub>)") +
    theme(
      axis.title.x = ggtext::element_markdown(),
      axis.title.y = ggtext::element_markdown(),
      panel.border = element_rect(linewidth = 0.3)
    )

ggsave("figuras/preliminares/subdiarios/perfil_q1_100_proflik.png", width = 12, height = 6, units = "cm", dpi = 300)




# IDF --------------------------------------------------------------------

## AVALIAÇÃO ESPACIAL DOS PARÂMETROS DA d-GEV ----------------------------

### Interpolação IDW -----------------------------------------------------


pacman::p_load(sf, terra, tidyterra, gstat)

par.dgev.sf <- par.dgev.df |> 
  left_join(subdaily.sf, by = "gauge_code") |> 
  select(gauge_code, parname, par, ci_lower, ci_upper, model, geometry) |> 
  st_as_sf()

# Plotar c/ xit
par.sf <- par.dgev.sf |> 
  filter(model == prior.names[4])

# ou plotar c/ xi0
par.sf <- par.sf |> 
  # st_drop_geometry() |> 
  pivot_wider(
    id_cols = c(gauge_code, geometry),
    names_from = parname,
    values_from = par
  ) |> 
  mutate(xi0 = xit*alpha0, .after = xit) |> 
  pivot_longer(
    cols = -c(gauge_code, geometry),
    names_to = "parname",
    values_to = "par"
  ) |> 
  filter(parname != "xit")


# Criar uma camada `SpatRaster`
bbox.rs
buffer <- 5
grid <- st_make_grid(rs.sf, cellsize = 0.1, what = "polygons")
grid.sf <- st_sf(geometry = grid) |> st_intersection(rs.sf)

run.idw <- function(points, grid, parname, idp = 2.0){
  form <- as.formula(paste(parname, "~ 1"))
  gstat::idw(formula = form, locations = points, newdata = grid, idp = idp)
}

get.breaks <- function(x, n){
  probs <- seq(0, 1, length.out = n + 1)
  quantile(x, probs, type = 1, names = FALSE)
}

# Definir intervalos prévios
par.breaks <- split(par.sf$par, par.sf$parname) |> 
  lapply(function(x) get.breaks(x, 4)) |> 
  print()

# Redefinir quebras de alguns
par.breaks[["kappa"]] <- c(-0.21, -0.07, 0, 0.07, 0.22)
par.breaks[["theta"]] <- c(0, 0.1, 0.6, 1.560)

# Interpolação IDW
par.r <- par.sf |> 
  group_split(parname) |> 
  map_dfr(function(df){
    par.name <- unique(df$parname)
    run.idw(df, grid.sf, "par", idp = 2) |> 
      mutate(
        parname = par.name,
        interval = cut(var1.pred, par.breaks[[par.name]], right = FALSE, include.lowest = TRUE)
      )
  }); gc()

rownames(par.r) <- NULL
colors.ls <- list(
  # xit = RColorBrewer::brewer.pal(length(par.breaks[["xit"]]), "Greens"),
  xi0 = RColorBrewer::brewer.pal(length(par.breaks[["xi0"]]), "Greens"),
  alpha0 = RColorBrewer::brewer.pal(length(par.breaks[["alpha0"]]), "Blues"),
  kappa = c("#0571b0", "#92c5de", "#f7b5b5", "#db5454"),
  H = RColorBrewer::brewer.pal(length(par.breaks[["H"]]), "Purples"),
  theta = RColorBrewer::brewer.pal(length(par.breaks[["theta"]]), "Greys")
)

subtitle <- c("a", "b", "c", "d", "e")
# par.order <- c("xit", "alpha0", "kappa", "H", "theta")
# par.labels <- c("ξ&#771;", "α<sub>0</sub>", "κ", "H", "Θ")
par.order <- c("xi0", "alpha0", "kappa", "H", "theta")
par.labels <- c("ξ<sub>0</sub>", "α<sub>0</sub>", "κ", "H", "Θ")
par.ls <- split(par.r, par.r$parname)[par.order]
plot.ls <- lapply(seq_along(par.ls), function(i){
  df <- par.ls[[i]]
  par.name <- unique(df$parname)
  ggplot() +
    geom_sf(aes(fill = interval), df, color = NA) +
    geom_sf(data = subdaily.sf, color = "grey10", shape = 3, size = 1, inherit.aes = FALSE) +
    geom_sf(data = rs.sf, color = "black", fill = NA, linewidth = 0.3, inherit.aes = FALSE) +
    scale_fill_manual(name = NULL, values = colors.ls[[par.name]]) +
    coord_sf(crs = st_crs(4674), xlim = bbox.rs[c(1,3)], ylim = bbox.rs[c(2,4)], expand = TRUE, clip = "off") +
    annotation_scale(
      location = "br", width_hint = 0.2, height = unit(1.5, "mm"), text_family = font, text_cex = 0.6,
      bar_cols = c("black", "white"), line_width = 0.4
    ) +
    annotation_north_arrow(
      location = "tr" , width = unit(1, "cm"), height = unit(1, "cm"), which_north = "true",
      pad_x = unit(0.1, "pt"), pad_y = unit(4, "pt"),
      style = north_arrow_nautical(text_family = font, text_size = fontsize, fill = c("black", "white"))
    ) +
    annotate(
      geom = "richtext", x = -Inf, y = Inf, label = sprintf("(%s) %s", subtitle[i], par.labels[i]),
      family = font, size = textsize*1.2, hjust = -0.01, vjust = 1.1, fill = NA, label.color = NA
    ) +
    theme(
      legend.position = c(0, -0.02),
      legend.justification = c(0, 0),
      legend.key.spacing.y = unit(-1, "pt"),
      panel.border = element_rect(linewidth = 0.3),
      plot.margin = margin(3, 3, 3, 3)
    )
}) # fim 

# plot.dgev.par <- wrap_plots(plot.ls) + plot_layout(ncol = 2, axis = "collect")

layout <- "
  1111
  2222
  #33#
"

plot.dgev.par <- wrap_plots(plot.ls[1:2]) / wrap_plots(plot.ls[3:4]) / wrap_plots(plot.ls[5]) + plot_layout(design = layout)

plot.dgev.par

ggsave(
  "figuras/preliminares/subdiarios/mapa_dgev_par.png", plot.dgev.par,
  width = 15, height = 20, units = "cm", dpi = 300
)

# Informações complementares
par.dgev.df |> 
  filter(
    model == prior.names[4],
    parname == "theta",
    par == 0
  ) |> 
  nrow()


### Tabela com parâmetros estimados --------------------------------------

table.dgev <- par.dgev.df |> 
  mutate(id = gauge.ids$id[match(gauge_code, gauge.ids$gauge_code)], .after = gauge_code) |> 
  filter(model == prior.names[4]) |> # manter somente RGMLE
  select(gauge_code, parname, par, ci_lower, ci_upper) |> 
  mutate(
    table_res = sprintf("%4.3f [%3.2f;%3.2f]", par, ci_lower, ci_upper),
    table_res = str_replace(table_res, pattern = "\\.", replacement = ",")
  ) |>  
  pivot_wider(
    id_cols = gauge_code,
    names_from = parname,
    values_from = table_res
  )

table.names <- c("Codigo estacao", "$\\tilde{\\xi}$", "$\\alpha_0$", "$\\kappa$", "$H$", "$\\theta$")

# table.dgev <- par.dgev.df |> 
#   filter(model == prior.names[4]) |> # manter somente RGMLE
#   select(gauge_code, parname, par, ci_lower, ci_upper) |> 
#   mutate(
#     parname = case_when(
#       parname == "xit" ~ "$\\tilde{\\xi}$",
#       parname == "alpha0" ~ "$\\alpha_0$",
#       parname == "kappa" ~ "$\\kappa$",
#       parname == "H" ~ "H",
#       parname == "theta" ~ "$\\theta$"
#     ),
#     across(where(is.numeric), ~ gsub("\\.", ",", sprintf("%0.4f", .x)))
#   )

# table.names <- c("{Código estação}", "Parâmetro", "Estimativa", "Limite inferior", "Limite superior")

names(table.dgev) <- table.names
n.col <- ncol(table.dgev)
cap <- "Parametros ajustados do modelo d-GEV com maxima verossimilhanca generalizada utilizando informacao regional sobre o parametro de forma."
dgev.xtable <- table.dgev |> 
  xtable(
    caption = cap,
    label = "tab_dgev_rgmle",
    align = rep("c", n.col + 1)
  )

extra.rows <- list()
extra.rows$pos <- list(0)
ltx.cmd <- sprintf(
  "\\endfirsthead
  \\caption[]{%s (continuacao)} \\\\
  \\toprule
  %s \\\\
   \\midrule
  \\endhead
  \\midrule
  \\multicolumn{%d}{r}{{Continua na proxima pagina...}} \\\\
  \\endfoot
  \\endlastfoot \n",
  cap,
  paste(table.names, collapse = " & "),
  n.col
)

extra.rows$command <- ltx.cmd

dgev.table.string <- print(
  dgev.xtable, type = "latex", include.rownames = FALSE, caption.placement = "top",
  sanitize.colnames.function = identity, 
  sanitize.text.function = identity,
  tabular.environment = "longtable", # changes to longtable
  floating = FALSE,                  # required for longtable to work
  add.to.row = extra.rows,
  hline.after = c(-1)
)

cat(
  "\\begingroup\n\\scriptsize\n\\renewcommand{\\arraystretch}{0.8}",
  dgev.table.string,
  "\\endgroup\n",
  file = "tabelas/5_dgev_par_rgmle.txt",
  sep = ""
)

## CURVAS DE FREQUÊNCIA --------------------------------------------------

durations <- seq_len(24)

# Posições de plotagem
freq.obs.df <- imax.df |> 
  filter(d %in% durations) |> 
  group_by(gauge_code, d) |> 
  arrange(gauge_code, d, desc(imax)) |> 
  mutate(
    n = n(),
    i = row_number(),
    p = i/(n + 1),
    return_period = 1/p
  ) |> 
  ungroup()
    
freq.est.df <- idf.dgev.df |> 
  filter(
    d %in% durations,
    return_period <= 50
  ) |> 
  mutate(
    n_years = imax.df$n_years[match(gauge_code, imax.df$gauge_code)],
    gauge_label = sprintf("%s | %2g anos", gauge_code, n_years)
  ) |> 
  left_join(
    par.dgev.df |> 
      filter(parname == "kappa") |> 
      rename(kappa = par) |> 
      select(gauge_code, model, kappa),
    by = c("gauge_code", "model")
  )

prior.caption <- c("sem informação a priori", "priori regional informativa", "priori pouco informativa")
file.name <- c("nprior", "wprior", "wprior_ninfo")

model.idx <- 2

# Plot all gauges with facetted plots
ggplot() +
  facet_wrap(~factor(gauge_label, levels = order.n), ncol = 4, scale = "free_y") +
  geom_ribbon(
    aes(x = return_period, ymin = ci_lower, ymax = ci_upper, fill = factor(d)), 
    freq.est.df[freq.est.df$model == which.model[model.idx],],
    alpha = 0.1, linewidth = 0
  ) +
  geom_line(
    aes(x = return_period, y = q, color = factor(d)), 
    freq.est.df[freq.est.df$model == which.model[model.idx],],
    linewidth = 0.2
  ) +
  geom_point(
    aes(x = return_period, y = imax, color = factor(d), shape = factor(d)), freq.obs.df,
    size = 0.6
  ) +
  scale_x_log10() +
  scale_shape_manual(name = "Durações [h]", values = c(15, 16, 17, 18)) +
  scale_color_manual(name = "Durações [h]", values = c("red", "orange", "green", "steelblue")) +
  scale_fill_manual(name = "Durações [h]", values = c("red", "orange", "green", "steelblue")) +
  # labs(x = "T [anos]", y = "Intensidade [mm/h]", caption = prior.caption) +
  labs(x = "T [anos]", y = "Intensidade [mm/h]") +
  theme(
    text = element_text(size = fontsize*0.8),
  )

ggsave(
  sprintf("figuras/preliminares/idf_dgev_%s.png", file.name[model.idx]), 
  height = 22, width = 16, units = "cm", dpi = 300
)

# Plot a single chosen gauge
gg <- "A801"
ggplot() +
  facet_wrap(~factor(gauge_label, levels = order.n), ncol = 4, scale = "free_y") +
  geom_ribbon(
    aes(x = return_period, ymin = ci_lower, ymax = ci_upper, fill = factor(d)), 
    freq.est.df[freq.est.df$gauge_code == gg & freq.est.df$model == which.model[model.idx],],
    alpha = 0.1, linewidth = 0
  ) +
  geom_line(
    aes(x = return_period, y = q, color = factor(d)),
    freq.est.df[freq.est.df$gauge_code == gg & freq.est.df$model == which.model[model.idx],],
    linewidth = 0.2
  ) +
  geom_point(
    aes(x = return_period, y = imax, color = factor(d), shape = factor(d)), freq.obs.df[freq.obs.df$gauge_code == gg,],
    size = 0.6
  ) +
  scale_x_log10() +
  scale_shape_manual(name = "Durações [h]", values = c(15, 16, 17, 18)) +
  scale_color_manual(name = "Durações [h]", values = c("red", "orange", "green", "steelblue")) +
  scale_fill_manual(name = "Durações [h]", values = c("red", "orange", "green", "steelblue")) +
  labs(x = "T [anos]", y = "Intensidade [mm/h]", caption = prior.caption) +
  theme(
    plot.margin = margin(),
    # plot.caption = element_text(margin = margin())
  )

# Plot one figure per duration in `durations`
freq.plot.df <- bind_rows(
  freq.obs.df |> 
    select(gauge_code, gauge_label, d, imax, return_period) |> 
    rename(q = imax) |> 
    mutate(which = "Observado"),
  freq.est.df |> 
    select(gauge_code, gauge_label, d, q, ci_lower, ci_upper, return_period, model) |> 
    mutate(which = "Estimado")
)


# model.idx <- 2
# walk(seq_along(durations), function(i){

#   d <- durations[i]
#   data <- freq.plot.df[freq.plot.df$d == d,]
#   ggplot(data, aes(x = return_period)) + 
#     facet_wrap(~factor(gauge_label, levels = order.n), ncol = 4, scale = "free_y") +
#     geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), ~filter(., which == "Estimado", model == which.model[model.idx]), fill = d.colors[i], linewidth = 0, alpha = 0.1) +
#     geom_line(aes(y = q), ~filter(., which == "Estimado", model == which.model[model.idx]), linewidth = 0.3, color = d.colors[i]) +
#     geom_point(aes(y = q), ~filter(., which == "Observado"), size = 0.6, shape = d.shapes[i], color = d.colors[i]) +
#     scale_x_log10() +
#     labs(x = "T [anos]", y = sprintf("Intensidades para duração de %g h (%s) [mm/h]", d, prior.caption)) +
#     theme(
#       plot.margin = margin(),
#       text = element_text(size = fontsize*0.8),
#       plot.caption = element_text(margin = margin(t = -fontsize*2, unit = "pt"))
#     )
  
#   ggsave(sprintf("figuras/preliminares/idf_dgev_%s_d%g.png", file.name[model.idx], d), height = 22, width = 16, units = "cm", dpi = 300)

# })

walk(seq_along(durations), function(i){

  d <- durations[i]
  data <- freq.plot.df[freq.plot.df$d == d,]
  ggplot(data, aes(x = return_period)) + 
    facet_wrap(~factor(gauge_label, levels = order.n), ncol = 4, scale = "free_y") +
    geom_ribbon(
      aes(ymin = ci_lower, ymax = ci_upper, fill = model, color = model), 
      ~filter(., which == "Estimado"), linewidth = 0.2, linetype = "dashed", alpha = 0.1
    ) +
    geom_line(aes(y = q, color = model), ~filter(., which == "Estimado"), linewidth = 0.3) +
    geom_point(aes(y = q), ~filter(., which == "Observado"), size = 0.6, shape = d.shapes[i], color = "grey20") +
    scale_color_manual(name = "", values = model.col, aesthetics = c("color", "fill")) +
    scale_x_log10() +
    labs(x = "T [anos]", y = sprintf("Intensidades para duração de %g h [mm/h]", d)) +
    theme(
      text = element_text(size = fontsize*0.8),
      plot.caption = element_text(margin = margin(t = -fontsize*2, unit = "pt"))
    )
  
    ggsave(sprintf("figuras/preliminares/idf_dgev_d%g.png", d), height = 22, width = 16, units = "cm", dpi = 300)
  
})



## CURVA `Y` -------------------------------------------------------------


# Derivar 'growth curves'
par.wide.df <- par.dgev.df |> 
  # filter(parname %in% c("H", "theta")) |> 
  select(gauge_code, parname, par, model, prior) |> 
  pivot_wider(
    names_from = parname,
    values_from = par
  )

Y <- idf.dgev.df |> 
  select(gauge_code, d, return_period, q, model) |> 
  left_join(
    par.wide.df |> select(gauge_code, H, theta, model, prior),
    by = c("gauge_code", "model")
  ) |> 
  rowwise() |> 
  mutate(
    b = (d + theta)^H,
    y = q*b
  ) |> 
  ungroup() |> 
  slice_head(n = 1, by = c("gauge_code", "return_period", "model")) |> 
  select(gauge_code, return_period, y, model)

gauge.theta <- par.dgev.df |> 
  filter(
    parname == "theta",
    par > 0
  ) |> 
  arrange(desc(par)) |> 
  pull(gauge_code, name = "par")

gg <- gauge.theta[[1]]
gg <- "A801"

ggplot() +
  facet_wrap(~model) +
  geom_line(
    aes(x = return_period, y = q, color = "i<sub>d,p</sub>", group = factor(d)), 
    idf.dgev.df |> filter(gauge_code == gg)
  ) +
  geom_line(aes(x = return_period, y = y, color  = "a(T)"), Y |> filter(gauge_code == gg), linewidth = 1) +
  scale_color_manual(
    name = "", values = c("a(T)" = "black", "i<sub>d,p</sub>" = "grey"),
  ) +
  scale_x_log10() +
  labs(
    # x = "T [anos]", y = "y<sub>p</sub> = i<sub>d,p</sub>(d + θ)<sup>H</sup> [mm/h<sup>ω</sup>]", 
    x = "T [anos]", y = "a(T) [mm/h<sup>ω</sup>]", 
    subtitle = sprintf("Estação %s", gg)
  ) +
  theme(
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    # legend.position = c(0, 1.06),
    # legend.justification = c(0, 1),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "grey", linetype = "dashed")
  )
  
ggsave(sprintf("figuras/preliminares/growth_curve_%s.png", gg), height = 8, width = 12, units = "cm", dpi = 300)

Y.obs <- freq.obs.df |> 
  filter(gauge_code == gg) |> 
    left_join(
      par.wide.df |> select(gauge_code, H, theta, model, prior),
      by = c("gauge_code")
    ) |> 
    rowwise() |> 
    mutate(
      b = (d + theta)^H,
      y = imax*b
    ) |> 
    ungroup()


ggplot() +
  facet_wrap(~model) +
    # geom_line(aes(x = return_period, y = q, color = "i<sub>d,p</sub>", group = factor(d)), idf.df |> filter(gauge_code == gg)) +
  geom_line(
    aes(x = return_period, y = y, color  = "Y<sub>T</sub>"), 
    Y |> filter(gauge_code == gg, return_period <= 100), 
    linewidth = 0.5
  ) +
  geom_point(
    aes(x = return_period, y = y, color = "y<sub>T</sub><sup>obs</sup>", group = factor(d)), Y.obs, 
    size = 0.8, alpha = 0.8
  ) +
  scale_color_manual(
    name = "", values = c("Y<sub>T</sub>" = "black", "y<sub>T</sub><sup>obs</sup>" = "grey"),
  ) +
    scale_x_log10() +
      labs(
    # x = "T [anos]", y = "y<sub>p</sub> = i<sub>d,p</sub>(d + θ)<sup>H</sup> [mm/h<sup>ω</sup>]", 
    x = "T [anos]", y = "Y<sub>T</sub> [mm/h<sup>ω</sup>]", 
    subtitle = sprintf("Estação %s", gg)
  ) +
  theme(
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    # legend.position = c(0, 1.06),
    # legend.justification = c(0, 1),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "grey", linetype = "dashed")
  )
  
ggsave(sprintf("figuras/preliminares/growth_curve_pts_%s.png", gg), height = 8, width = 12, units = "cm", dpi = 300)

text.plot <- par.wide.df |> 
  filter(gauge_code ==  gg) |> 
  mutate(
    w = 1 - H,
    text = sprintf(
      "κ = %0.3f\nH = %0.3f\nθ = %0.3f\nω = %0.3f", 
      kappa, H, theta, w
    )
  ) |> print()

ggplot() +
  facet_wrap(~model) +
    # geom_line(aes(x = return_period, y = q, color = "i<sub>d,p</sub>", group = factor(d)), idf.df |> filter(gauge_code == gg)) +
  geom_line(
    aes(x = return_period, y = y, color  = "Y<sub>T</sub>"), 
    Y |> filter(gauge_code == gg, return_period <= 100), 
    linewidth = 0.5
  ) +
  geom_point(
    aes(x = return_period, y = y, color = factor(d), group = factor(d)), Y.obs, 
    size = 0.8, alpha = 0.8
  ) +
  geom_text(
    aes(x = 1.05, y = Inf, label = text), text.plot,
    color = "black", family = font, size = textsize, vjust = 1.2, hjust = 0
  ) +
  scale_color_manual(
    name = "", values = c("Y<sub>T</sub>" = "black", setNames(d.colors, durations)),
  ) +
  scale_x_log10() +
  labs(
    # x = "T [anos]", y = "y<sub>p</sub> = i<sub>d,p</sub>(d + θ)<sup>H</sup> [mm/h<sup>ω</sup>]", 
    x = "T [anos]", y = "Y<sub>T</sub> [mm/h<sup>ω</sup>]", 
    subtitle = sprintf("Estação %s", gg)
  ) +
  theme(
    axis.title.y = ggtext::element_markdown(),
    legend.text = ggtext::element_markdown(),
    # legend.position = c(0, 1.06),
    # legend.justification = c(0, 1),
    panel.grid.major.y = element_line(linewidth = 0.2, color = "grey", linetype = "dashed")
  )

ggsave(sprintf("figuras/preliminares/growth_curve_pts_%s.png", gg), height = 8, width = 12, units = "cm", dpi = 300)



## IDF -------------------------------------------------------------------



## MÉTRICAS DE DESEMPENHO ------------------------------------------------


# Frequências observadas para cada estação estão em `freq.obs.df`
# Aqui estamos recalculando as frequências, aproveitando o índice já calculado,
# mas incluindo as posições de weibull, cunnane e gringorten
freq.obs.all.df <- imax.df |> 
  group_by(gauge_code, d) |> 
  arrange(gauge_code, d, desc(imax)) |> 
  mutate(
    n = n(),
    i = row_number(),
    p_weibull = i/(n + 1),
    p_cunnane = (i - 0.4)/(n + 1 - 2*0.4),
    p_gringorten = (i - 0.44)/(n + 1 - 2*0.44)
  ) |> 
  pivot_longer(
    cols = starts_with("p_"),
    names_to = "plot_pos",
    names_prefix = "p_",
    values_to = "p"
  ) |> 
  mutate(
    return_period = 1/p
  ) |> 
  ungroup()


# Organizar os dados de `par.gev.df` e `par.dgev.df`
# Combinar e simplificar a estrutura dos dois
par.df <- bind_rows(
  par.gev.df |> mutate(model = "GEV"),
  par.dgev.df |> mutate(model = "d-GEV")
) |> 
  select(gauge_code, model, d, parname, par, ci_lower, ci_upper, prior) |> 
  group_by(gauge_code, model, d, prior) |> 
  reframe(par = list(setNames(par, parname))) # listar parâmetros


# Falta ajustar a IDF propriamente dita do modelo sem invariância de escala
# os quantis estimados são só as curvas, a equação em si ainda precisa ser calculada.
duration.vec <- unique(imax.df$d)


# Ajustar IDF-SGB c/ código Wilson
# O script retorna o data.frame `idf.sgb.sf`
source("R/fit_sgb_idf.R")
idf.sgb.df <- idf.sgb.df |> select(gauge_code, k, gamma, b, c, dpma)


# Estimativas de intensidade p/ frequências observadas
freq.obs.sgb <- freq.obs.all.df |> 
  left_join(idf.sgb.df, by = "gauge_code") |> 
  select(!c(gauge_label, n, rain_yr))

performance.sgb.df <- freq.obs.sgb |> 
  mutate(
    d_min = d*60, # IDF-SGB usa duração em minutos
    q = purrr::pmap_dbl(
      .l = list(k, gamma, b, c, return_period, d_min), # modelo IDF-SGB usa `b` como parâmetro de deslocamento (`theta`)
      .f = function(k, gamma, theta, eta, rp, ds){     # e `c` como o expoente do denominador (`eta` ou `H`)
        if(any(is.na(c(k, gamma, eta, theta)))){       # a função `idf.sgb()` está escrita com `eta` e `theta`
          return(NA_real_)
        } else{
          idf.sgb(par = c(k, gamma, eta, theta), return.period = rp, duration = ds)
        }
      }
    ),
    residual = imax - q, # calcular resíduos
    model = "IDF-SGB"    # atribuir nome ao modelo p/ combinar com outros depois
  )


# Juntar frequências
# Cada posição de plotagem deve ter 4 linhas, uma p/ cada priori
performance.gev.df <- freq.obs.all.df |> 
  inner_join(
    par.df |> filter(model == "GEV"),
    by = c("gauge_code", "d"),
    relationship = "many-to-many"
  )

performance.dgev.df <- freq.obs.all.df |> 
  inner_join(
    par.df |> filter(model == "d-GEV") |> select(!d),
    by = "gauge_code",
    relationship = "many-to-many"
  )


# Juntar tudo de novo
performance.df <- bind_rows(performance.gev.df, performance.dgev.df) |> 
  mutate(
    q = purrr::pmap_dbl(
      .l = list(model, par, d, p), # lista de parâmetros
      .f = function(mod, param, ds, prob){
        if(mod == "GEV"){
          qgev(param = param, p.nexceed = 1 - prob) # converter em probabilidade de não excedência
        } else{
          qdgev(param = param, duration = ds, p.nexceed = 1 - prob) # converter em probabilidade de não excedência
        }
      }
    ),
    residual = imax - q
  ) |> 
  bind_rows(performance.sgb.df) |> # juntar c/ resultados c/ IDF-SGB
  select(gauge_code, n, d, plot_pos, model, prior, p, return_period, imax, q, residual) |> print()


# Calcular métricas, os somatórios ocorrem nas probabilidades
performance.summary.df <- performance.df |> 
  group_by(gauge_code, d, plot_pos, prior, model) |> 
  summarise(
    return_period = first(return_period),
    n = n(),
    sum_imax = sum(imax, na.rm = TRUE),
    r_rmse = sqrt(n*sum(residual/sum_imax)^2),
    r_bias = n*sum(residual/sum_imax),
    .groups = "drop"
  ) |> 
  mutate(
    model_label = paste0(model, "_", gsub(" .*", "", prior.names[prior]))
  )

# Gráfico de barras das métricas p/ todas as durações
# ordenado pelo valor da métrica de desempenho
col <- c("#648FFF", "#dc267f", "#eb78b0", "#FE6100")
prior.col <- setNames(col, prior.names)

# Plotar gráficos diferentes p/ cada duração em vez de um
# facet por duração
# rRMSE
width <- 0.6
# pos <- position_dodge2(width + 0.2, preserve = "single", padding = 0.05)
pos <- position_dodge(width + 0.2)

# Gráfico barras GEV e d-GEV separados
performance.summary.df |> 
  filter(
    plot_pos == "weibull",
    # # plot_pos == "cunnane",
    # plot_pos == "gringorten",
    model == "d-GEV", # alterar modelo aqui
    d %in% c(1, 6, 12, 24)
  ) |> 
  mutate(prior = factor(prior, which.prior[c(1, 4, 3, 2)])) |> 
  ggplot(aes(x = fct_reorder(gauge_code, r_rmse, .desc = TRUE), y = r_rmse, fill = prior)) +
  facet_wrap(~factor(d), ncol = 1, scales = "free_y", labeller = as_labeller(function(x) sprintf("Duração: %s h", x))) +
  geom_col(position = pos, width = width, linewidth = 0.2) +
  coord_cartesian(ylim = c(0, 0.6)) +
  scale_fill_manual(name = "Priori", values = setNames(col[c(1, 4, 3, 2)], which.prior)) +
  labs(y = "rRMSE") +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = fontsize*0.8),
    axis.title.x = element_blank(),
    legend.box = "vertical",
    axis.ticks.x = element_line(linewidth = 0.2, size = 0.2),
    legend.spacing.y = unit(0, "pt")
  ) +
  guides(
    alpha = guide_legend(order = 1),
    fill = guide_legend(order = 2),
  ); invisible(gc())


# Boxplots

models <- unique(performance.summary.df$model_label)
sgb.col <- "grey50"
which.position <- "weibull"

# Comparar modelos com e sem distribuição a priori
durations.metrics <- c(1, 3, 6, 9, 12, 18, 24)

# Quantas observações acima de um limite
lim.rmse <- 0.7
outliers.rmse <- performance.summary.df |> 
  filter(
    plot_pos == which.position,
    r_rmse > lim.rmse,                   # selecionar 
    model_label %in% models[c(1, 2, 3, 4)], # escolher modelos
    d %in% durations.metrics
  ) |> 
  mutate(
    id = gauge.ids$id[match(gauge_code, gauge.ids$gauge_code)],
    d = fct_reorder(paste(d, "h"), d)
  ) |> 
  select(id, d, prior, model, r_rmse) |> 
  arrange(id) |> 
  group_by(d, prior, model) |> 
  summarise(
    min_r_rmse = min(r_rmse),
    string = sprintf("[%s]", paste(id, collapse = ",")),
    .groups = "drop"
  )

# Comparar com boxplots as estimativas com a GEV ajustada por duração e a d-GEV
# com diferentes níveis de informação adicional
performance.boxplot <- performance.summary.df |> 
  filter(
    plot_pos == which.position,
    model_label %in% models[c(1, 2, 5, 6)],
    d %in% durations.metrics
  ) |> 
  mutate(
    d = fct_reorder(paste(d, "h"), d),
    model_label = factor(model_label, models)
  ) |> 
  select(gauge_code, d, model_label, plot_pos, r_rmse, r_bias) |> 
  ggplot(aes(x = d, y = r_rmse, fill = model_label, alpha = model_label)) +
  geom_boxplot(
    linewidth = 0.3, median.linewidth = 0.3, staplewidth = 0.5, whisker.linetype = "dashed", whisker.linewidth = 0.3,
    outlier.size = 0.5, outlier.alpha = 0.5, key_glyph = "rect"
  ) +
  geom_text(aes(x = d, y = lim.rmse, label = string), outliers.rmse, inherit.aes = FALSE, family = font, size = textsize*0.8) +
  coord_cartesian(ylim = c(0, lim.rmse)) +
  scale_fill_manual(
    name = "GEV",
    values = setNames(rep(col[c(1, 4)], each = 2), models[c(1, 2, 5, 6)]),
    breaks = models[c(1,5)],
    labels = gsub(".*_", "", models[c(1,5)]),
    guide = guide_legend(
      override.aes = list(
        alpha = 0.4, 
        color = "black",
        linewidth = 0.2
      ))
    ) +
  scale_alpha_manual(
    name = "d-GEV",
    values = setNames(c(0.4, 0.95, 0.4, 0.95), models[c(1, 2, 5, 6)]),
    breaks = models[c(2,6)],
    labels = gsub(".*_", "", models[c(2,6)]),
    guide = guide_legend(
      override.aes = list(
        fill = setNames(col[c(1, 4)], models[c(2,6)]), 
        color = "black",
        linewidth = 0.2
      ))
  ) +
  labs(x = NULL, y = "rRMSE", caption = paste("posicao_plotagem:", which.position)) +
  theme(
    legend.position = c(0.99, 0.92),
    legend.justification = c(1, 1),
    legend.spacing.y = unit(-2, "pt"),
    legend.title.position = "top",
    panel.grid.major.y = element_line(linewidth = 0.2, linetype = "solid", color = "grey90"),
    legend.title = element_text(margin = margin(b = 1), hjust = 0)
  ); performance.boxplot

ggsave("figuras/preliminares/subdiarios/boxplot_rmse_gev.png", performance.boxplot, width = 15, height = 10, units = "cm", dpi = 300)
# ggsave(
#   sprintf("figuras/preliminares/subdiarios/analise_pos_plot/boxplot_rmse_gev_%s.png", which.position), performance.boxplot, 
#   width = 15, height = 10, units = "cm", dpi = 300
# )


# Comparar com boxplots as estimativas com a GEV ajustada por duração e a d-GEV
# com diferentes níveis de informação adicional e também com IDF-SGB

# Quantas observações acima de um limite
lim.rmse <- 1
outliers.rmse <- performance.summary.df |> # se não for mostrar a GEV, não precisa mostrar esses outliers
  filter(
    plot_pos == which.position,
    r_rmse > lim.rmse,
    model_label %in% models[c(2, 6, 7)],
    d %in% durations.metrics
  ) |> 
  mutate(
    id = gauge.ids$id[match(gauge_code, gauge.ids$gauge_code)],
    d = fct_reorder(paste(d, "h"), d)
  ) |> 
  select(id, d, prior, model, r_rmse) |> 
  arrange(id) |> 
  group_by(d, prior, model) |> 
  summarise(
    min_r_rmse = min(r_rmse),
    string = sprintf("[%s]", paste(id, collapse = ",")),
    .groups = "drop"
  )

width <- 0.8

# RMSE relativo
rmse.boxplot <- performance.summary.df |> 
  filter(
    plot_pos == which.position,
    model_label %in% models[c(2, 6, 7)],
    d %in% durations.metrics
  ) |> 
  mutate(
    d = fct_reorder(paste(d, "h"), d),
    model_label = factor(model_label, levels = models[c(7,2,6)])
  ) |> 
  select(gauge_code, d, model_label, plot_pos, r_rmse) |> 
  ggplot(aes(x = d, y = r_rmse, fill = model_label)) +
  geom_boxplot(
    linewidth = 0.2, median.linewidth = 0.2, staplewidth = 0.5, whisker.linetype = "dashed", whisker.linewidth = 0.2,
    outlier.size = 0.5, outlier.alpha = 0.5, key_glyph = "rect", width = width
  ) +
  annotate(
    geom = "text_npc", label = "Modelo de referência:\nIntensidades observadas (Weibull)", npcx = 0.02, npcy = 0.85, vjust = 1.2,
    family = font, size = textsize
  ) +
  scale_fill_manual(
    name = NULL, 
    values = setNames(c(sgb.col, prior.col[c(1, 4)]), models[c(7, 2, 6)]),
    labels = c("IDF-SGB", "d-GEV (MLE)", "d-GEV (RGMLE)")
  ) +
  labs(x = NULL, y = "rRMSE") +
  # labs(x = NULL, y = "rRMSE", caption = paste("posicao_plotagem:", which.position)) +
  theme(
    legend.position = c(0.01, 0.975),
    legend.justification = c(0, 1),
    legend.spacing.y = unit(-5, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.title.position = "top",
    panel.grid.major.y = element_line(linewidth = 0.2, linetype = "solid", color = "grey90"),
    legend.title = element_text(margin = margin(b = 1), hjust = 0)
  ); rmse.boxplot

ggsave("figuras/preliminares/subdiarios/boxplot_rmse.png", rmse.boxplot, width = 15, height = 10, units = "cm", dpi = 300)


### Informações complementares -------------------------------------------

# Mediana do rRMSE p/ diferentes durações
performance.summary.df |> 
  # filter(plot_pos == which.position, model == "IDF-SGB") |> 
  filter(
    plot_pos == which.position, model == "d-GEV", prior %in% c("NULL", "c(0.0178, 0.0627)"),
  ) |> 
  group_by(model, prior, d = factor(d)) |> 
  summarise(
    q25 = quantile(r_rmse, probs = 0.25),
    q50 = quantile(r_rmse, probs = 0.5),
    q75 = quantile(r_rmse, probs = 0.75),
    iqr = q75 - q25,
    max = max(r_rmse, na.rm = TRUE),
    max_gauge = gauge_code[which.max(r_rmse)],
    max_id = gauge.ids$id[match(max_gauge, gauge.ids$gauge_code)],
    .groups = "drop"
  ) |> arrange(prior) |> print(n = Inf)

# Mediana do rBIAS p/ diferentes durações
performance.summary.df |> 
  filter(
    plot_pos == which.position, model == "d-GEV", prior %in% c("NULL", "c(0.0178, 0.0627)"),
  ) |> # como são próximos, acho razoável 
  group_by(model, prior, d = factor(d)) |>                       # agrupar os dois modelos da d-GEV 
  summarise(
    q25 = quantile(r_bias, probs = 0.25),
    q50 = quantile(r_bias, probs = 0.5),
    q75 = quantile(r_bias, probs = 0.75),
    iqr = q75 - q25,
    max = max(r_bias, na.rm = TRUE),
    max_gauge = gauge_code[which.max(r_bias)],
    max_id = gauge.ids$id[match(max_gauge, gauge.ids$gauge_code)],
    .groups = "drop"
  ) |> arrange(prior) |> print(n = Inf)

idf.sgb.df[idf.sgb.df$gauge_code == "432260802A",]

# Quantas estações tiverem viés negativo p/ 1 h no modelo convencional
performance.summary.df |> 
  filter(model == "IDF-SGB", d == 24, plot_pos == which.position, r_bias < 0) |> 
  select(gauge_code, d, model, return_period, n, r_rmse, r_bias) |> 
  slice(1, .by = gauge_code) |> 
  print(n = Inf)

# Em quantas estações os modelos d-GEV subestimam (r_bias < 0)
performance.summary.df |> 
  filter(
    plot_pos == which.position, model == "d-GEV", prior %in% c("NULL", "c(0.0178, 0.0627)"),
  ) |> 
  select(gauge_code, d, model, prior, return_period, n, r_rmse, r_bias) |> 
  slice(1, .by = gauge_code)


# ------------------------------------------------------------------------

# BIAS relativo
lim.bias <- 5.3
outliers.bias <- performance.summary.df |>
  filter(
    plot_pos == which.position,
    r_bias > lim.bias | r_bias < -lim.bias,
    model_label %in% models[c(7, 2, 6)],
    d %in% durations.metrics
  ) |> 
  mutate(
    id = gauge.ids$id[match(gauge_code, gauge.ids$gauge_code)],
    d = fct_reorder(paste(d, "h"), d)
  ) |> 
  select(id, d, prior, model, r_bias) |> 
  arrange(id) |> 
  group_by(d, prior, model) |> 
  summarise(
    min_r_bias = min(r_bias),
    max_r_bias = max(r_bias),
    string_pos = sprintf("[%s]", paste(id[r_bias > 0], collapse = ",")),
    string_neg = sprintf("[%s]", paste(id[r_bias < 0], collapse = ",")),
    .groups = "drop"
  )

bias.boxplot <- performance.summary.df |> 
  filter(
    plot_pos == which.position,
    model_label %in% models[c(7, 2, 6)],
    d %in% c(1, 3, 6, 9, 12, 18, 24)
  ) |> 
  mutate(
    d = fct_reorder(paste(d, "h"), d),
    model_label = factor(model_label, levels = models[c(7,2,6)])
  ) |> 
  select(gauge_code, d, model, model_label, plot_pos, r_bias, r_bias) |> 
  ggplot(aes(x = d, y = r_bias, fill = model_label)) +
  geom_boxplot(
    linewidth = 0.2, median.linewidth = 0.2, staplewidth = 0.5, whisker.linetype = "dashed", whisker.linewidth = 0.2,
    outlier.size = 0.5, outlier.alpha = 0.5, key_glyph = "rect", width = width
  ) +
  annotate(
    geom = "text_npc", label = "Modelo de referência:\nIntensidades observadas (Weibull)", npcx = 0.02, npcy = 0.85, vjust = 1.2,
    family = font, size = textsize
  ) +
  scale_fill_manual(
    name = NULL, 
    values = setNames(c(sgb.col, prior.col[c(1, 4)]), models[c(7, 2, 6)]),
    labels = c("IDF-SGB", "d-GEV (MLE)", "d-GEV (RGMLE)")
  ) +
  labs(x = NULL, y = "rBIAS") +
  theme(
    legend.position = c(0.01, 0.975),
    legend.justification = c(0, 1),
    legend.spacing.y = unit(-5, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.title.position = "top",
    panel.grid.major.y = element_line(linewidth = 0.2, linetype = "solid", color = "grey90"),
    legend.title = element_text(margin = margin(b = 1), hjust = 0)
  ); bias.boxplot

ggsave("figuras/preliminares/subdiarios/boxplot_bias.png", bias.boxplot, width = 15, height = 10, units = "cm", dpi = 300)


# Juntar rRMSE e rBIAS na mesma figura
plot.metricas.obs <- 
  wrap_plots(
    rmse.boxplot + 
      labs(subtitle = "(a) rRMSE", y = NULL) + 
      theme(plot.margin = margin(0.5, 0.5, 3, 0.5))
  ) +
  wrap_plots(
    bias.boxplot + 
      labs(subtitle = "(b) rBIAS", y = NULL) +
      theme(plot.margin = margin(3, 0.5, 0.5, 0.5))
  ) +
  plot_layout(ncol = 1); plot.metricas.obs

ggsave(
  "figuras/preliminares/subdiarios/boxplot_metricas_obs.png", plot.metricas.obs,
  width = 15, height = 20, units = "cm", dpi = 300
)

### Comparação entre quantis da GEV (RGMLE) e estimativas da d-GEV -------

# Recalcular métricas rRMSE e rBIAS, mas agora para tempos de retorno maiores
# usando os quantis da GEV (RGMLE) por duração como os "observados"

priors <- unique(par.df$prior)
high.rp <- unique(idf.gev.df$return_period)


# Estimativa dos parâmetros da GEV com L-momentos
# par.gev.lmom <- imax.df |> 
#   mutate(duration_min = d*60) |> 
#   group_by(gauge_code, d, duration_min) |> 
#   summarise(par = map(list(imax), start.gev), .groups = "drop") |> 
#   unnest_wider(par)

# Estimativa dos quantis da GEV com L-mmomenntos
# idf.lmom <- imax.df |> 
#   mutate(duration_min = d*60) |> 
#   group_by(gauge_code, d, duration_min) |> 
#   group_modify(~fgev.lmom(.x$imax, rp)) |> # group_modify usa tbl_df como entrada
#   ungroup()

idf.lmom <- idf.lmom |> # esse arquivo vêm do fit_idf_sgb.R
  rename(q = intensity)

rp <- c(10, 20, 50, 100)

# Juntar `idf.gev.df` e `idf.dgev.df` c/ informações
# dos ajustes c/ incerteza estimada por perfil de verossimilhança
performance.high.rp <- idf.gev.df |>                   # estimativas c/ GEV por duração
  filter(return_period %in% rp, prior == priors[4]) |> # priors[4] é RGMLE
  select(gauge_code, d, return_period, q) |>           # não selecionar coluna `prior`
  rename(q_gev = q) |> 
  left_join(                                           # estimativas c/ d-GEV
    idf.dgev.df |> select(gauge_code, d, return_period, prior, q) |> rename(q_dgev = q),
    by = c("gauge_code", "d", "return_period")
  ) |> 
  left_join(
    idf.sgb.df |> select(-dpma),                       # estimativas IDF-SGB
    by = c("gauge_code")
  ) |> 
  rowwise() |> 
  mutate(
    q_sgb = idf.sgb(par = c(k, gamma, c, b), return.period = return_period, duration = d*60)
  ) |> 
  ungroup() |> 
  select(-c(k, gamma, c, b)) |> 
  pivot_longer(
    cols = -c(gauge_code, d, return_period, prior, q_gev), # estimativas c/ GEV por duração não entram no pivot
    names_to = "model",
    names_prefix = "q_",
    values_to = "q_model"
  ) |> 
  filter(prior %in% priors[c(1, 4)]) |> 
  mutate(
    residual = q_gev - q_model,
    model_label = tolower(paste0(model, "_", gsub(" .*", "", prior.names[prior])))
  ) |> 
  group_by(gauge_code, d, model_label) |> 
  summarise(
    r_rmse = sqrt(n()*sum(residual/sum(q_gev, na.rm = TRUE))^2),
    r_bias = n()*sum(residual/sum(q_gev, na.rm = TRUE)),
    .groups = "drop"
  ) |> 
  filter(model_label != "sgb_mle") |> # só precisa manter 1 modelo sgb, pq são iguais
  mutate(model_label = if_else(model_label == "sgb_rgmle", "sgb", model_label))


# Cores
high.rp.models <- unique(performance.high.rp$model_label)[c(3,1,2)]
high.rp.models <- factor(high.rp.models, high.rp.models)

# Outliers --- rRMSE
lim.rmse <- 1.5 # corta alguns outliers do modelo d-GEV (MLE)
outliers.rmse <- performance.high.rp |> 
  filter(
    r_rmse > lim.rmse,
    d %in% durations.metrics
  ) |> 
  mutate(
    id = gauge.ids$id[match(gauge_code, gauge.ids$gauge_code)],
    d = fct_reorder(paste(d, "h"), d)
  ) |> 
  arrange(id) |> 
  group_by(d, model_label) |> 
  summarise(
    min_r_rmse = min(r_rmse),
    string = sprintf("[%s]", paste(id, collapse = ",")),
    .groups = "drop"
  )

width <- 0.8

# rRMSE
boxplot.high.rp.rmse <- performance.high.rp |> 
  filter(d %in% durations.metrics) |> 
  mutate(
    d = fct_reorder(paste(d, "h"), d)
  ) |> 
  ggplot(aes(x = d, y = r_rmse, fill = factor(model_label, high.rp.models))) +
  geom_boxplot(
    linewidth = 0.2, median.linewidth = 0.2, staplewidth = 0.5, whisker.linetype = "dashed", whisker.linewidth = 0.2,
    outlier.size = 0.5, outlier.alpha = 0.5, key_glyph = "rect", width = width
  ) +
  geom_text(
    aes(x = d, y = lim.rmse, label = string), outliers.rmse, inherit.aes = FALSE,
    family = font, size = textsize*0.8
  ) +
  annotate(
    geom = "text_npc", label = "Modelo de referência: GEV (RGMLE)", npcx = 0.98, npcy = 0.93, vjust = 1.2,
    family = font, size = textsize
  ) +
  coord_cartesian(ylim = c(0, lim.rmse)) +
  scale_fill_manual(
    name = NULL, values = setNames(c(sgb.col, prior.col[c(1, 4)]), high.rp.models),
    labels = c("IDF-SGB", "d-GEV (MLE)", "d-GEV (RGMLE)")
  ) +
  labs(x = NULL, y = "rRMSE") +
  theme(
    legend.position = c(0.01, 0.93),
    legend.justification = c(0, 1),
    legend.spacing.y = unit(-5, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.title.position = "top",
    panel.grid.major.y = element_line(linewidth = 0.2, linetype = "solid", color = "grey90"),
    legend.title = element_text(margin = margin(b = 1), hjust = 0)
  ); boxplot.high.rp.rmse

ggsave(
  "figuras/preliminares/subdiarios/boxplot_rmse_maiores_trs_lmom.png", boxplot.high.rp.rmse,
  width = 15, height = 10, units = "cm", dpi = 300
  )


# Outliers --- rBIAS
lim.bias <- 3 # só d-GEV (MLE)
outliers.bias <- performance.high.rp |> 
  filter(
    r_bias > lim.bias | r_bias < -lim.bias,
    d %in% durations.metrics
  ) |> 
  mutate(
    id = gauge.ids$id[match(gauge_code, gauge.ids$gauge_code)],
    d = fct_reorder(paste(d, "h"), d)
  ) |> 
  arrange(id) |> 
  group_by(d) |>  # não precisa agrupar por modelo, só quando quiser saber qual que é, mas nesse caso são todos d-GEV (MLE)
  summarise(
    min_r_bias = min(r_bias), max_r_bias = max(r_bias),
    string_pos = sprintf("[%s]", paste(id[r_bias > 0], collapse = ",")),
    string_neg = sprintf("[%s]", paste(id[r_bias < 0], collapse = ",")),
    .groups = "drop"
  )

# rBIAS
boxplot.high.rp.bias <- performance.high.rp |> 
  filter(d %in% durations.metrics) |> 
  mutate(
    d = fct_reorder(paste(d, "h"), d)
  ) |> 
  ggplot(aes(x = d, y = r_bias, fill = factor(model_label, high.rp.models))) +
  geom_boxplot(
    linewidth = 0.2, median.linewidth = 0.2, staplewidth = 0.5, whisker.linetype = "dashed", whisker.linewidth = 0.2,
    outlier.size = 0.5, outlier.alpha = 0.5, key_glyph = "rect", width = width
  ) +
  geom_text(
    aes(x = d, y = -lim.bias, label = string_neg), outliers.bias, inherit.aes = FALSE,
    family = font, size = textsize*0.8
  ) +
  annotate(
    geom = "text_npc", label = "Modelo de referência: GEV (RGMLE)", npcx = 0.98, npcy = 0.975, vjust = 1.2,
    family = font, size = textsize
  ) +
  coord_cartesian(ylim = c(-lim.bias, lim.bias)) +
  scale_fill_manual(
    name = NULL, values = setNames(c(sgb.col, prior.col[c(1, 4)]), high.rp.models),
    labels = c("IDF-SGB", "d-GEV (MLE)", "d-GEV (RGMLE)")
  ) +
  labs(x = NULL, y = "rBIAS") +
  theme(
    legend.position = c(0.01, 0.975),
    legend.justification = c(0, 1),
    legend.spacing.y = unit(-5, "pt"),
    legend.key.spacing.y = unit(1, "pt"),
    legend.title.position = "top",
    panel.grid.major.y = element_line(linewidth = 0.2, linetype = "solid", color = "grey90"),
    legend.title = element_text(margin = margin(b = 1), hjust = 0)
  ); boxplot.high.rp.bias

ggsave(
  "figuras/preliminares/subdiarios/boxplot_bias_maiores_trs.png", boxplot.high.rp.bias,
  width = 15, height = 10, units = "cm", dpi = 300
  )

# Juntar todo mundo
plot.boxplot.high.rp <- 
  wrap_plots(
    boxplot.high.rp.rmse + 
      labs(subtitle = "(a) rRMSE", y = NULL) + 
      theme(plot.margin = margin(0.5, 0.5, 3, 0.5))
  ) +
  wrap_plots(
    boxplot.high.rp.bias + 
      labs(subtitle = "(b) rBIAS", y = NULL) +
      theme(plot.margin = margin(3, 0.5, 0.5, 0.5))
  ) +
  plot_layout(ncol = 1)

ggsave(
  "figuras/preliminares/subdiarios/boxplot_metricas_maiores_trs.png", plot.boxplot.high.rp,
  width = 15, height = 20, units = "cm", dpi = 300
)

# Informações complementares --- rRMSE
boxplot.summary.rp <- performance.high.rp |> 
  group_by(model_label, d) |> 
  summarise(
    q25 = quantile(r_rmse, probs = 0.25),
    q50 = quantile(r_rmse, probs = 0.5),
    q75 = quantile(r_rmse, probs = 0.75),
    iqr = q75 - q25,
    max = max(r_rmse, na.rm = TRUE),
    max_gauge = gauge_code[which.max(r_rmse)],
    max_id = gauge.ids$id[match(max_gauge, gauge.ids$gauge_code)],
    .groups = "drop"
  ) |> print(n = Inf)

boxplot.summary.rp |> filter(model_label == "dgev_mle") |> print(n = Inf)   # d-GEV (MLE)
boxplot.summary.rp |> filter(model_label == "dgev_rgmle") |> print(n = Inf) # d-GEV (RGMLE)
boxplot.summary.rp |> filter(model_label == "sgb") |> print(n = Inf)        # IDF-SGB

# Informações complementares --- rBIAs
boxplot.summary.rp <- performance.high.rp |> 
  group_by(model_label, d) |> 
  summarise(
    q25 = quantile(r_bias, probs = 0.25),
    q50 = quantile(r_bias, probs = 0.5),
    q75 = quantile(r_bias, probs = 0.75),
    iqr = q75 - q25,
    max = max(r_bias, na.rm = TRUE),
    max_gauge = gauge_code[which.max(r_bias)],
    max_id = gauge.ids$id[match(max_gauge, gauge.ids$gauge_code)],
    .groups = "drop"
  ) |> print(n = Inf)

boxplot.summary.rp |> filter(model_label == "dgev_mle") |> print(n = Inf)   # d-GEV (MLE)
boxplot.summary.rp |> filter(model_label == "dgev_rgmle") |> print(n = Inf) # d-GEV (RGMLE)
boxplot.summary.rp |> filter(model_label == "sgb") |> print(n = Inf)        # IDF-SGB

# Avaliando estimativas de theta
gg.max.bias <- unique(boxplot.summary.rp$max_gauge)
par.dgev.wide |> 
  filter(gauge_code %in% gg.max.bias, prior %in% prior.names[c(1,4)]) |> 
  print(n = Inf)

# Conferir se há estações com r_bias > 0 p/ todas as durações
performance.high.rp |> 
  group_by(gauge_code, model_label) |> 
  summarise(positive_bias = sum(r_bias > 0)/24 > 0.5, .groups = "drop") |> 
  pivot_wider(
    id_cols = gauge_code,
    names_from = model_label,
    values_from = positive_bias
  ) |> print(n = Inf)

# Correlação entre os parâmetros e rBIAS
par.order <- c("xit", "alpha0", "kappa", "H", "theta")
par.labels <- c("ξ&#771;", "α<sub>0</sub>", "κ", "H", "Θ")
error.par.correlation <- performance.high.rp |> 
  filter(model_label == "dgev_rgmle") |> 
  left_join(
    par.dgev.df |> 
      filter(model == prior.names[4]) |> 
      select(gauge_code, parname, par),
    by = "gauge_code", relationship = "many-to-many"
  ) |> 
  mutate(parname = factor(parname, par.order, par.labels))

p_load(ggpmisc)

error.par.correlation |> 
  ggplot(aes(x = par, y = r_bias)) +
  facet_wrap(~parname, ncol = 1, scales = "free_x") +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey80") +
  geom_point(alpha = 0.5, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4) +
  ggpmisc::stat_poly_eq(
    formula = y ~ x, 
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
    parse = TRUE, family = font, size = textsize, label.x = "right", label.y = "top"
  ) +
  labs(y = "rBIAS", x = NULL) +
  theme(strip.text = element_markdown())

ggsave(
  "figuras/preliminares/subdiarios/correlacao_vies_par_dgev_rgmle.png",
  width = 15, height = 20, units = "cm", dpi = 600
)

error.par.correlation |> 
  filter(parname == par.labels[3]) |> 
  ggplot(aes(x = par, y = r_bias)) +
  facet_wrap(~parname, ncol = 1, scales = "free_x") +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey80") +
  # geom_point(shape = 21, aes(color = d), fill = NA, size = 1) +
  geom_point(shape = 21, color = "black", fill = NA, size = 1) +
  geom_smooth(method = "lm", se = FALSE, linewidth = 0.4, color = "red") +
  ggpmisc::stat_correlation(mapping = use_label("R"), method = "pearson", family = font, size = textsize) +
  ggpmisc::stat_poly_eq(
    formula = y ~ x, 
    aes(label = paste(after_stat(eq.label), after_stat(rr.label), sep = "~~~")), 
    parse = TRUE, family = font, size = textsize, label.x = "right", label.y = "top"
  ) +
  labs(y = "rBIAS", x = par.labels[3]) +
  theme(strip.text = element_blank())

ggsave(
  "figuras/preliminares/subdiarios/correlacao_vies_kappa_dgev_rgmle.png",
  width = 15, height = 10, units = "cm", dpi = 600
)


# ------------------------------------------------------------------------


# Comparar frequências observadas e modeladas plotando as IDFs
# versus as intensidades observadas

# Visualizar IDFs SGB
par.labels <- performance.sgb.df |> 
  slice_head(by = gauge_code) |> 
  mutate(
    gauge_label = gauge.labels[gauge_code],
    eta = c, theta = b
  ) |> 
  select(gauge_label, k, gamma, eta, theta)

# Visualização da IDF de todas as estações e comparar IDF-SGB com quantis da GEV por duração
idf.sgb.plot <- performance.sgb.df |> 
  mutate(
    return_period = factor(return_period),
    gauge_label = gauge.labels[gauge_code]
  ) |> 
  filter(plot_pos == "weibull") |> 
  ggplot(aes(x = d, group = return_period)) +
  facet_wrap(~gauge_label, ncol = 4) +
  geom_line(aes(y = imax, color = "Intensidades observadas"), linewidth = 0.3) +
  geom_line(aes(y = q, color = "IDF-SGB"), linewidth = 0.1) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("k: %7.3f", k)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 1.8, size = textsize*0.55
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("γ: %7.3f", gamma)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 3.3, size = textsize*0.55
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("η: %7.3f", eta)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 4.8, size = textsize*0.55
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("Θ: %7.3f", theta)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 6.3, size = textsize*0.55
  ) +
  scale_x_log10(name = "Duração [h]", breaks = c(1, 3, 9, 24)) +
  scale_y_log10(name = "Intensidade [mm/h]") +
  scale_color_manual(
    name = NULL, 
    values = c("Intensidades observadas" = "grey90", "IDF-SGB" = "grey20")
  ) +
  theme(
    text = element_text(size = fontsize*0.8),
    strip.text = element_text(size = fontsize*0.6, margin = margin(b = 2)),
    panel.border = element_rect(linewidth = 0.3, color = "black"),
    legend.margin = margin()
  ) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 0.4))
  ); idf.sgb.plot

ggsave(
  "figuras/preliminares/subdiarios/idf_sgb.png", idf.sgb.plot, 
  width = 15.5, height = 23, units = "cm", dpi = 600
)


# Visualizar d-GEV
par.labels <- performance.dgev.df |> 
  # filter(model == "d-GEV", prior == names(prior.names)[[1]]) |> # MLE
  filter(model == "d-GEV", prior == names(prior.names)[[4]]) |> # RGMLE
  slice_head(by = gauge_code) |> 
  mutate(gauge_label = gauge.labels[gauge_code]) |> 
  select(gauge_label, model, prior, par) |> 
  unnest_wider(par)

# Visualizar da IDF de todas as estações e comparar IDF-SGB com quantis da GEV por duração
idf.dgev.plot <- performance.df |> 
  filter(
    plot_pos == "weibull",
    model == "d-GEV"
  ) |> 
  mutate(
    gauge_label = gauge.labels[gauge_code]
  ) |> 
  ggplot(aes(x = d, group = factor(return_period))) +
  facet_wrap(~gauge_label, ncol = 4) +
  geom_line(aes(y = imax, color = "Intensidades observadas"), linewidth = 0.3) +
  # geom_line(aes(y = q, color = "d-GEV (MLE)"), ~filter(., prior == names(prior.names)[1]), linewidth = 0.1) +
  geom_line(aes(y = q, color = "d-GEV (RGMLE)"), ~filter(., prior == names(prior.names)[4]), linewidth = 0.1) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("κ: %7.3f", kappa)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 1.8, size = textsize*0.55
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("H: %7.3f", H)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 3.3, size = textsize*0.55
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("θ: %7.3f", theta)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 4.8, size = textsize*0.55
  ) +
  scale_x_log10(name = "Duração [h]") +
  scale_y_log10(name = "Intensidade [mm/h]") +
  scale_color_manual(
    name = NULL,
    # values = c("d-GEV (MLE)" = prior.col[[1]], "Intensidades observadas" = "grey90") # MLE
    values = c("d-GEV (RGMLE)" = prior.col[[4]], "Intensidades observadas" = "grey90") # RGMLE
  ) +
  theme(
    text = element_text(size = fontsize*0.8),
    strip.text = element_text(size = fontsize*0.6, margin = margin(b = 2)),
    panel.border = element_rect(linewidth = 0.3, color = "black"),
    legend.margin = margin()
  ) +
  guides(
    color = guide_legend(override.aes = list(linewidth = 0.4))
  ); idf.dgev.plot

# ggsave(
#   "figuras/preliminares/subdiarios/idf_dgev_mle.png", idf.dgev.plot, 
#   width = 15.5, height = 23, units = "cm", dpi = 600
# )
ggsave(
  "figuras/preliminares/subdiarios/idf_dgev_rgmle.png", idf.dgev.plot, 
  width = 15.5, height = 23, units = "cm", dpi = 600
)


# ------------------------------------------------------------------------

# Avaliar o viés espacialmente



## CURVAS DE FREQUÊNCIA --------------------------------------------------

# Plotar curvas de frequências p/ comparar quantis estimados c/ 

# Verificar estações c/ erros nos ICs
idf.dgev.df |> pull(error_message) |> unique() # NA

# Tempos de retorno p/ quais os quantis c/ incerteza foram calculados
available.rp <- idf.dgev.df |> pull(return_period) |> unique()

par.dgev.wide <- par.dgev.df |> 
  pivot_wider(
    id_cols = c(gauge_code, prior, model),
    names_from = parname,
    values_from = par
  ) |> 
  mutate(
    gauge_label = gauge.labels[gauge_code],
    prior = model, model = "d-GEV", 
    .after = gauge_code
  )

# Juntar resultados
idf.df <- idf.dgev.df |> # quantis estimados c/d-GEV
  select(gauge_code, model, d, return_period, q, ci_lower, ci_upper) |> 
  mutate(
    gauge_label = gauge.labels[gauge_code],
    prior = model, model = "d-GEV", .after = gauge_code
  ) |> 
  left_join( 
    par.dgev.wide, # juntar parâmetros estimados
    by = c("gauge_code", "gauge_label", "prior", "model")
  )

# Função p/ fazer figura
idf.plot <- function(df, par, plot.gauges, plot.model, plot.prior, plot.rp, colors,
  tsize, par.label.x = 1.1, par.label.y = 2.3, prior.label, ncol = 2
){

  # Filtrar parametros
  par <- par |> filter(model == plot.model, prior == plot.prior)

  # Filtrar quantis
  df <- df |> 
    filter(model == plot.model, prior == plot.prior, return_period %in% plot.rp) |> 
    mutate(return_period = factor(return_period))

  # Filtrar estações
  if(!is.null(plot.gauges)) df <- df[df$gauge_code %in% plot.gauges,]
  if(!is.null(plot.gauges)) par <- par[par$gauge_code %in% plot.gauges,]

  model.label <- sprintf("%s (%s)", plot.model, prior.label)

  ggplot(df, aes(x = d, y = q, color = return_period, fill = return_period)) +
    facet_wrap(~gauge_label, ncol = ncol) +
    geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper), alpha = 0.2, linewidth = 0) +
    geom_line(linewidth = 0.2) +
    annotate(
      geom = "text", x = Inf, y = Inf, label = model.label, 
      hjust = par.label.x - nchar(model.label)/400, vjust = par.label.y, family = font, size = tsize
    ) +
    geom_text(
      aes(x = Inf, y = Inf, label = sprintf("κ = %7.3f", kappa)), par, inherit.aes = FALSE,
      hjust = par.label.x, vjust = par.label.y + tsize, family = font, size = tsize
    ) +
    geom_text(
      aes(x = Inf, y = Inf, label = sprintf("H = %7.3f", H)), par, inherit.aes = FALSE,
      hjust = par.label.x, vjust = par.label.y + 2*tsize, family = font, size = tsize
    ) +
    geom_text(
      aes(x = Inf, y = Inf, label = sprintf("Θ = %7.3f", theta)), par, inherit.aes = FALSE,
      hjust = par.label.x, vjust = par.label.y + 3*tsize, family = font, size = tsize
    ) +
    scale_color_manual(
      name = NULL, values = colors, aesthetics = list("color", "fill"),
      labels = paste(plot.rp, "anos")
    ) +
    scale_x_log10(name = "Duração [h]", breaks = c(1, 3, 9, 24)) +
    scale_y_log10(name = "Intensidade [mm/h]") +
    theme(legend.position = c(0, -0.01), legend.justification = c(0, 0))

} # fim `idf.plot`

# Argumentos
plot.rp <- c(5, 20, 100)
plot.gauges <- c("A801", "A803", "A813", "A844")
# colors <- c("black", "grey40", "grey70")
colors <- c("#a53f00", "#FE6100", "#ffa56d")

# Plotar p/ 4 estações --- RGMLE
idf.rgmle.plot <- idf.plot(
  df = idf.df, par = par.dgev.wide,
  plot.gauges = plot.gauges,
  plot.model = "d-GEV",
  plot.prior = prior.names[4], # RGMLE
  prior.label = "RGMLE",
  plot.rp = plot.rp,
  tsize = textsize*0.8,
  colors = colors
); idf.rgmle.plot

ggsave(
  "figuras/preliminares/subdiarios/idf_rgmle_rp.png", idf.rgmle.plot,
  width = 15, height = 12, units = "cm", dpi = 600
)

# Plotar p/ todas as estações --- RGMLE
idf.rgmle.todas.plot <- idf.plot(
  df = idf.df, par = par.dgev.wide,
  plot.gauges = NULL,
  plot.model = "d-GEV",
  plot.prior = prior.names[4], # RGMLE
  prior.label = "RGMLE",
  plot.rp = plot.rp,
  tsize = textsize*0.6, par.label.y = 2,
  colors = colors,
  ncol = 4
) +
  theme(
    text = element_text(size = fontsize*0.8),
    strip.text = element_text(size = fontsize*0.6, margin = margin(b = 2)),
    panel.border = element_rect(linewidth = 0.3, color = "black"),
    legend.margin = margin(),
    legend.position = "bottom", legend.justification = "center",
    axis.title = element_text(margin = margin())
  ); idf.rgmle.todas.plot

ggsave(
  "figuras/preliminares/subdiarios/idf_rgmle_rp_todas.png", idf.rgmle.todas.plot,
  width = 15.5, height = 23, units = "cm", dpi = 600
)


colors <- c("#183889", "#648FFF", "#a2bdff")

# Plotar p/ 4 estações --- MLE
idf.mle.plot <- idf.plot(
  df = idf.df, par = par.dgev.wide,
  plot.gauges = plot.gauges,
  plot.model = "d-GEV",
  plot.prior = prior.names[1], # MLE
  prior.label = "MLE",
  plot.rp = plot.rp,
  tsize = textsize*0.8,
  colors = colors
); idf.mle.plot

ggsave(
  "figuras/preliminares/subdiarios/idf_mle_rp.png", idf.mle.plot,
  width = 15, height = 12, units = "cm", dpi = 600
)

# Plotar p/ todas estações --- MLE
idf.mle.plot.todas <- idf.plot(
  df = idf.df, par = par.dgev.wide,
  plot.gauges = NULL,
  plot.model = "d-GEV",
  plot.prior = prior.names[1], # MLE
  prior.label = "MLE",
  plot.rp = plot.rp,
  tsize = textsize*0.6, par.label.y = 2,
  colors = colors, 
  ncol = 4
) +
  theme(
    text = element_text(size = fontsize*0.8),
    strip.text = element_text(size = fontsize*0.6, margin = margin(b = 2)),
    panel.border = element_rect(linewidth = 0.3, color = "black"),
    legend.margin = margin(),
    legend.position = "bottom", legend.justification = "center",
    axis.title = element_text(margin = margin())
  ); idf.mle.plot.todas

ggsave(
  "figuras/preliminares/subdiarios/idf_mle_rp_todas.png", idf.mle.plot.todas,
  width = 15.5, height = 23, units = "cm", dpi = 600
)


# Informações complementares
rgmle.influence <- idf.dgev.df |> 
  filter(model %in% prior.names[c(1,4)]) |> 
  select(gauge_code, d, return_period, q, ci_lower, ci_upper, model) |> 
  mutate(
    ci_width = ci_upper - ci_lower,
    rel_ci_width = ci_width/q
  ) |> 
  left_join(
    par.dgev.wide |> filter(prior %in% prior.names[c(1,4)]) |> select(gauge_code, prior, kappa),
    by = join_by(gauge_code, model == prior)
  ) |> 
  mutate(model = ifelse(model == "MLE", "mle", "rgmle")) |> 
  pivot_wider(
    id_cols = c(gauge_code, d, return_period),
    names_from = model,
    values_from = c(rel_ci_width, kappa)
  )

# Avaliar estimativas de kappa parecidas
rgmle.influence |> 
  mutate(diff_kappa = abs(kappa_rgmle - kappa_mle), .before = kappa_mle) |> 
  arrange(diff_kappa) |> View()

# Avaliar relação entre kappa e diminuição de IC
changes <- rgmle.influence |> 
  mutate(
    n_years = imax.df$n_years[match(gauge_code, imax.df$gauge_code)],
    diff_kappa = kappa_rgmle < kappa_mle,
    smaller_ci_with_rgmle = rel_ci_width_rgmle < rel_ci_width_mle,
    .after = gauge_code
  ) |> 
  filter(
    n_years %in% c(8,9)
  )

table(changes$diff_kappa, changes$smaller_ci_with_rgmle)


# ------------------------------------------------------------------------



# Exemplo curva IDF Porto Alegre SGB Atlas Pluviométrico
# Esta figura vai p/ exemplificar uma curva IDF no referencial teórico


rps <- c(2, 5, 10, 25, 50, 100)
ds <- c(seq(5, 60, 5), seq(2, 24, 2)*60) # minutos
k <- 812.9
gamma <- 0.1571
eta <- 0.7486
theta <- 6.6

idf.matrix <- sapply(rps, function(rp) k*rp^gamma/(ds + theta)^eta)
dimnames(idf.matrix) <- list(ds, rps)

idf.df <- idf.matrix |> 
  as.data.frame() |> 
  rownames_to_column(var = "d") |>
  pivot_longer(
    cols = !d,
    names_to = "return_period",
    names_transform = as.numeric,
    values_to = "q"
  ) |> 
  mutate(
    d = as.numeric(d),
    rp_label = fct_reorder(sprintf("T = %g anos", return_period), return_period)
  ) |> print()

idf.col <- c("")

idf.df |> 
  ggplot(aes(x = d/60, y = q, color = rp_label)) +
  geom_line(key_glyph = "rect") +
  scale_x_log10(name = "Duração [h]") +
  scale_y_log10(name = "Intensidade [mm/h]") +
  annotate(geom = "text", label = sprintf("k: %-6.1f", k), x = Inf, y = 295, vjust = 0, hjust = 1.2, family = font, size = textsize) +
  annotate(geom = "text", label = sprintf("γ: %-6.4f", gamma), x = Inf, y = 258, vjust = 0, hjust = 1.2, family = font, size = textsize) +
  annotate(geom = "text", label = sprintf("η: %-6.4f", eta), x = Inf, y = 224, vjust = 0, hjust = 1.2, family = font, size = textsize) +
  annotate(geom = "text", label = sprintf("θ: %-6.1f", theta), x = Inf, y = 198, vjust = 0, hjust = 1.2, family = font, size = textsize) +
  scale_color_manual(
    name = NULL, values = rev(c("black", "grey20", "grey40", "grey60", "grey70", "grey85")),
    guide = guide_legend(reverse = TRUE)
  ) +
  theme(
    legend.position = c(0.02, 0.01),
    legend.justification = c(0,0),
    axis.title.x = element_text(margin = margin(t = 5))
  )

ggsave("figuras/preliminares/idf_sgb_porto_alegre_atlas.png", width = 12, height = 8, units = "cm", dpi = 300)
