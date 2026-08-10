
rm(list = ls()); invisible(gc())


# PACOTES ----------------------------------------------------------------

pacman::p_load(
    pacman, tidyverse, 
    geomtextpath, ggspatial, patchwork,
    stars, ncmeta, FNN, cubelyr, # esses depois do stars são pacotes exigidos
    terra
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

# LER DADOS --------------------------------------------------------------

# Arquivos do Xavier de precipitação em NetCDF
nc.files <- list.files("base/sig/precipitacao_xavier/", pattern = "\\.nc$" ,full.names = TRUE)
nc <- read_stars(nc.files, along = "time", var = "pr", proxy = TRUE)
# star.ls <- lapply(nc.files, read_ncdf, proxy = TRUE) # read_netcdf só lê um arquivo por vez

# Carregando como proxy, o R não carrega as informações na memória, apenas aponta para o arquivo
# no diretório informado. Tanto é que se rodar `print(nc)`, não haverá informação, somente o caminho.

# nc <- do.call(c, c(star.ls, along = "time")) # combinando dessa forma faz com que 

min.date <- "1961-01-01"
raw.dates <- st_get_dimension_values(nc, which = "time")
dates <- as.Date(raw.dates, origin = min.date)
nc <- st_set_dimensions(nc, which = "time", values = dates)

# Contorno Rio Grande do Sul
rs.sf <- read_sf("base/sig/sf_rs.gpkg")


# ORGANIZAR NETCDF -------------------------------------------------------

# Atribuir sistema de coordenadas
nc.crs <- 4326 # GEE: print(xavierImageCollection.first().projection().getInfo());
nc <- nc |> st_set_crs(st_crs(nc.crs))

# Transformar rs.sf em WGS 84 (o mais leve)
rs.wgs84.sf <- st_transform(rs.sf, st_crs(nc))

# Recortar extensão RS
rs.wgs84.buffer.sf <- st_buffer(rs.wgs84.sf, dist = 1e4) # 2x a mair que o tamanho do grid do Xavier

ggplot() +
  geom_sf(data = rs.wgs84.sf, fill = NA) +
  geom_sf(data = rs.wgs84.buffer.sf, fill = NA)

nc.rs <- st_crop(nc, rs.wgs84.buffer.sf) # apenas um proxy
# `st_crop()` pode não funcionar quando usa `read_ncdf()`, pois requer que xy sejam retilineos e
# NetCDFs podem apresentar pequenas variações de xy ao longo do tempo. Ler como `read_stars()`
# automaticamente "suaviza" essas diferenças
nc.rs <- st_as_stars(nc.rs); gc()
names(nc.rs) <- "pr"

# RESUMIR EM TOTAIS E MÁXIMOS ANUAIS MÉDIOS ------------------------------

# `aggregate` precisa de uma grade retangular p/ processar, como fizemos um crop
# há muitos pixeis com NA, logo é comum que max(NA, na.rm = TRUE) retorne -Inf
# max(NA, na.rm = TRUE)

# Funções personalizadas p/ tratar NAs
# Máximo
safe.max <- function(x){
  if(all(is.na(x))){
    return(NA_real_)
  } else{
    return(max(x, na.rm = TRUE))
  }
}

# Soma
safe.sum <- function(x){
  if(all(is.na(x))){
    return(NA_real_)
  } else{
    return(sum(x, na.rm = TRUE))
  }
}

# Máximos anuais
nc.annual.max <- nc.rs |> 
  aggregate(by = "years", FUN = safe.max) # como `time` já está em Date, podemos usar by = "years" direto

invisible(gc())

# nc.annual.max <- st_set_dimensions(nc.annual.max, which = "time", names = "year")
# nc.annual.max[["pr"]][is.infinite(nc.annual.max[["pr"]])] <- NA_real_

# Totais anuais
nc.annual.sum <- nc.rs |> 
  aggregate(by = "years", FUN = safe.sum)

invisible(gc())

# Totais mensais
nc.monthly.sum <- nc.rs |> 
  aggregate(by = "3 months", FUN = safe.sum)

invisible(gc())

# Médias
nc.mean.max <- st_apply(nc.annual.max, MARGIN = c("x", "y", "time"), FUN = mean) # para todo par {x,y}, aplique a média ao longo de "time"
nc.mean.sum <- st_apply(nc.annual.sum, MARGIN = c("x", "y"), FUN = mean)
nc.summary <- c(nc.mean.max, nc.mean.sum, nc.mean.mon.sum)
names(nc.summary) <- c("mean_annual_max", "mean_annual_sum", "mean_monthly_sum")

# nc.summary <- nc.summary |>
#   st_warp(crs = st_crs(rs.sf)) |> 
#   st_crop(rs.buffer.sf)

# Resumos anuais
mean.max.df <- as.data.frame(nc.summary["mean_annual_max"])
mean.sum.df <- as.data.frame(nc.summary["mean_annual_sum"])

# Resumos mensais
mon.groups <- unique(month(st_get_dimension_values(nc.monthly.sum, which = "time")))
mon.labels <- setNames(c("JFM", "AMJ", "JAS", "OND"), mon.groups)
mon.sum.df <- as.data.frame(nc.monthly.sum) |> 
  mutate(
    month = month(time),
    period = fct_reorder(mon.labels[factor(month)], month)
  ) |> 
  group_by(period, x, y) |> 
  summarise(mean = mean(pr), .groups = "drop")

# Criar rasters
max.rast <- rast(nc.mean.max)
sum.rast <- rast(nc.mean.sum)

# Interpolar em uma grade mais fina p/ plotar
disagg.fact <- 10

mean.max.smooth.df <- max.rast |> 
  disagg(fact = disagg.fact, method = "bilinear") |> 
  as.data.frame(xy = TRUE, na.rm = TRUE)

mean.sum.smooth.df <- sum.rast |> 
  disagg(fact = disagg.fact, method = "bilinear") |> 
  as.data.frame(xy = TRUE, na.rm = TRUE)


# VISUALIZAR -------------------------------------------------------------

# Primeiro ano da série
plot(slice(nc.annual.max, "time", 1), main = "Máximos anuais em 1961", breaks = "quantile")
plot(slice(nc.annual.sum, "time", 1), main = "Máximos anuais em 1961", breaks = "quantile")

# Médias
plot(nc.summary["mean_annual_max"], main = "Totais médios anuais")
plot(nc.summary["mean_annual_sum"], main = "Média máxima anual")

# Precipitação máxima diária anual
n.bins <- 7
pretty.breaks <- pretty(mean.max.smooth.df$mean, n = n.bins)
q.breaks <- round(quantile(mean.max.smooth.df$mean, probs = seq(0, 1, length.out = n.bins)))
plot.max <- ggplot(data = mean.max.smooth.df, aes(x = x, y = y, z = mean)) +
  geom_contour_filled(
    bins = n.bin, breaks = q.breaks
  ) +
  geom_textcontour(
    bins = n.bin, breaks = q.breaks, family = font, size = textsize*0.65, color = "black", linewidth = 0.3,
    vjust = 0.5, hjust = 0.2, upright = TRUE, straight = TRUE
  ) +
  geom_sf(data = rs.sf, fill = NA, linewidth = 0.4, color = "black", inherit.aes = FALSE) +
  coord_sf(crs = st_crs(rs.sf)) +
  scale_fill_brewer(
    name = "Precipitação\nmáxima diária\nanual [mm]", palette = 1, direction = 1
  ) +
  annotation_scale(
      location = "br", width_hint = 0.2, height = unit(1.5, "mm"), text_family = font, text_cex = 0.6,
      bar_cols = c("black", "aliceblue"), line_width = 0.4
  ) +
  annotation_north_arrow(
      location = "tr" , width = unit(1, "cm"), height = unit(1, "cm"), which_north = "true",
      style = north_arrow_nautical(text_family = font, text_size = fontsize, fill = c("black", "aliceblue"))
  ) +
  # labs(x = "Longitude", y = "Latitude", subtitle = "(a)") +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    legend.position = c(0.01, 0),
    legend.justification = c(0, 0),
    legend.key.spacing.y = unit(-1.4, "pt")
  ); plot.max

ggsave("figuras/preliminares/diarios/mapa_isoietas_max.png", plot.max, height = 12, width = 12, units = "cm", dpi = 300)

# Precipitação acumulada anual
n.bins <- 7
pretty.breaks <- pretty(mean.sum.smooth.df$mean, n = n.bins)
q.breaks <- round(quantile(mean.sum.smooth.df$mean, probs = seq(0, 1, length.out = n.bins)))
manual.breaks <- c(1200, 1300, 1500, 1700, 1800, 1900)
plot.sum <- ggplot(data = mean.sum.smooth.df, aes(x = x, y = y, z = mean)) +
  geom_contour_filled(
    bins = n.bin, breaks = q.breaks
  ) +
  geom_textcontour(
    bins = n.bin, breaks = q.breaks, family = font, size = textsize*0.65, color = "black", linewidth = 0.3,
    vjust = 0.5, hjust = 0.2, upright = TRUE, straight = TRUE
  ) +
  geom_sf(data = rs.sf, fill = NA, linewidth = 0.4, color = "black", inherit.aes = FALSE) +
  coord_sf(crs = st_crs(rs.sf)) +
  scale_fill_brewer(
    name = "Precipitação\nacumulada\nanual [mm]", palette = 1, direction = 1
  ) +
  annotation_scale(
      location = "br", width_hint = 0.2, height = unit(1.5, "mm"), text_family = font, text_cex = 0.6,
      bar_cols = c("black", "aliceblue"), line_width = 0.4
  ) +
  annotation_north_arrow(
      location = "tr" , width = unit(1, "cm"), height = unit(1, "cm"), which_north = "true",
      style = north_arrow_nautical(text_family = font, text_size = fontsize, fill = c("black", "aliceblue"))
  ) +
  # labs(x = "Longitude", y = "Latitude", subtitle = "(b)") +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    legend.position = c(0.01, 0),
    legend.justification = c(0, 0),
    legend.key.spacing.y = unit(-1.4, "pt")
  ); plot.sum

ggsave("figuras/preliminares/diarios/mapa_isoietas_sum.png", plot.sum, height = 12, width = 12, units = "cm", dpi = 300)

plot.isoietas <-  wrap_plots(plot.max) + wrap_plots(plot.sum) + plot_layout(ncol = 1, guides = "keep"); plot.isoietas

# Precipitação total mensal
n.bins <- 7
q.breaks <- round(quantile(mon.sum.df$mean, probs = seq(0, 1, length.out = n.bins), na.rm = TRUE), 0)
manual.breaks <- c(240, 340, 370, 400, 430, 450, 570)
plot.mon.sum <- ggplot(mon.sum.df, aes(x, y, z = mean)) +
  facet_wrap(~period, ncol = 2) +
  geom_contour_filled(
    bins = n.bin, breaks = q.breaks
  ) +
  geom_textcontour(
    bins = n.bin, breaks = q.breaks, family = font, size = textsize*0.65, color = "black", linewidth = 0.3,
    vjust = 0.5, hjust = 0.2, upright = TRUE, straight = TRUE
  ) +
  geom_sf(data = rs.sf, fill = NA, linewidth = 0.4, color = "black", inherit.aes = FALSE) +
  coord_sf(crs = st_crs(rs.sf)) +
  scale_fill_brewer(
    name = "Precipitação\nmensal [mm]", palette = 1, direction = 1
  ) +
  annotation_scale(
      location = "br", width_hint = 0.2, height = unit(1.5, "mm"), text_family = font, text_cex = 0.6,
      bar_cols = c("black", "aliceblue"), line_width = 0.4
  ) +
  annotation_north_arrow(
      location = "tr" , width = unit(1, "cm"), height = unit(1, "cm"), which_north = "true",
      style = north_arrow_nautical(text_family = font, text_size = fontsize, fill = c("black", "aliceblue"))
  ) +
  # labs(x = "Longitude", y = "Latitude", subtitle = "(b)") +
  labs(x = "Longitude", y = "Latitude") +
  theme(
    legend.position = c(0.01, 0),
    legend.justification = c(0, 0),
    legend.key.spacing.y = unit(-1.4, "pt")
  ); plot.mon.sum

ggsave("figuras/preliminares/diarios/mapa_isoietas_meses_sum.png", plot.mon.sum, height = 15.5, width = 15.5, units = "cm", dpi = 300)
