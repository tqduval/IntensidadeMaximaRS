rm(list = ls()); invisible(gc())

pacman::p_load(pacman, tidyverse, arrow, terra, sf, rnaturalearth)

font <- "JetBrainsMono NF"
fontsize <- 8
textsize <- (fontsize - 2)/.pt
set_theme(
    theme_minimal() +
    theme(
      panel.background = element_rect(color = "black", fill = "white"),
      panel.border = element_rect(color = "black", linewidth = 0.3),
      panel.grid = element_blank(),
      axis.text = element_text(color = "black"),
      text = element_text(color = "black", family = font, size = fontsize)
    )
)


# LER DADOS --------------------------------------------------------------

rs.sf <- st_read("base/sig/sf_rs.gpkg")

imax.df <- arrow::read_parquet("base/gerados/df_imax.pqt")

raster.files <- list.files(path = "base/sig/anadem_rs/", pattern = ".tif$", full.names = TRUE)
raster.collection <- lapply(raster.files, terra::rast)
raster.anadem <- terra::sprc(raster.collection)
raster.anadem <- terra::mosaic(raster.anadem)

anadem <- terra::crop(raster.anadem, rs.sf, mask = TRUE)

# Filter dada
na.yr <- 0.2
min.durations <- 3
d.interval <- c(1, 24)
min.yr <- 8

imax.df <- imax.df |> 
  filter(
    na_prct_yr <= na.yr,
    between(d, d.interval[1], d.interval[2])
  ) |> 
  group_by(gauge_code) |> 
  filter(n_distinct(d) >= min.durations,
  n_distinct(wateryear) >= min.yr) |> 
  ungroup() |> 
  select(gauge_code, wateryear, d, imax)
    
imax.ls <- split(imax.df, imax.df$gauge_code)

gauges <- unique(imax.df$gauge_code)

info.df <- arrow::read_parquet("base/fonte/subdaily/df_subdaily_info.pqt")
info.df <- info.df |> 
  filter(gauge_code %in% gauges) |> 
  select(gauge_code, long, lat, elevation)

gauges.sf <- info.df |> st_as_sf(coords = c("long", "lat"), crs = st_crs(rs.sf))

cols <- names(imax.df)
prior.info <- c(0.0178, 0.0627)
method <- "L-BFGS-B"
scale.inv <- "wide"

# CHAMAR FUNÇÕES ---------------------------------------------------------

source("R/funcoes/fdgev.R")


# AJUSTE -----------------------------------------------------------------

fit.ls <- lapply(imax.ls, function(data.gg){

  # data.gg <- imax.ls[[1]]
  # gg <- data.gg[[cols[1]]][1]
  
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
    cols = cols,
    method = method,
    scale.inv = scale.inv
  )
  
  par <- fit$par
  xit <- par[[1]]; alpha0 <- par[[2]]; kappa <- par[[3]]; H <- par[[4]]; theta <- par[[5]]
  xi0 <- xit*alpha0

  out <- c(
    xit = xit,
    xi0 = xi0,
    alpha0 = alpha0,
    kappa = kappa,
    H = H,
    theta = theta
  )
  
  return(out)

})

fit.df <- bind_rows(fit.ls, .id = "gauge_code")


# COVARIÁVEIS ------------------------------------------------------------


## ELEVAÇÃO --------------------------------------------------------------

# Cota das estações estão contidas no próprio conjunto de dados

## DISTÂNCIA DA COSTA ----------------------------------------------------

# Extrair costa
coast.sf <- rnaturalearth::ne_coastline(scale = 110) |> st_transform(crs = st_crs(rs.sf))
bbox.rs <- st_bbox(rs.sf) |> st_as_sfc() |> st_buffer(dist = 0.01)
coast.sf <- st_crop(coast.sf, bbox.rs) |> st_transform(crs = st_crs(gauges.sf))
crs(coast.sf)
crs(gauges.sf)
coast.dist <- st_distance(gauges.sf, coast.sf)
gauges.sf$coast_dist = apply(coast.dist, 1, min)/1000
segments.sf <- st_nearest_points(gauges.sf, coast.sf) |> 
  st_sf(
    gauge_code = gauges.sf$gauge_code,
    coast_dist = gauges.sf$coast_dist,
    geometry = _ # `_` é o placeholder `.` do pipe padrão do R
  )

ggplot() +
  tidyterra::geom_spatraster(data = anadem, alpha = 0.7) +
  geom_sf(data = rs.sf, fill = NA, color = "black") +
  geom_sf(data = segments.sf, color = "#6dbcfc", linetype = "dashed", alpha = 0.5, linewidth = 0.3) +
  geom_sf(data = coast.sf, color = "#6dbcfc", linewidth = 0.8) +
  geom_sf(data = gauges.sf, color = "#6dbcfc", size = 0.8) +
  labs(subtitle = "Distância até a costa") +
  scale_fill_viridis_c(
    name = "", na.value = NA, option = "cividis", # option = "turbo") +
    guide = guide_colorsteps(show.limits = TRUE, order = 2)
  ) +
  theme(
    panel.border = element_blank(),
    legend.position = c(0, 0),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    legend.key.height = unit(fontsize*0.9, "pt")
    # legend.spacing.x = unit(fontsize*2, "pt")
    # legend.spacing.y = unit(0, "pt")
  )

ggsave("figuras/dist_coast_subdaily.png", height = 8, width = 8, units = "cm", dpi = 300)

gg <- "A801"
ggplot() +
  geom_sf(data = rs.sf, fill = "grey90", color = "black") +
  geom_sf(data = segments.sf[segments.sf$gauge_code == gg,], color = "firebrick", linetype = "dashed", alpha = 0.5) +
  geom_sf(data = coast.sf, color = "dodgerblue3", linewidth = 0.8) +
  geom_sf(data = gauges.sf[gauges.sf$gauge_code == gg,], color = "firebrick") +
  geom_sf_text(
    data = gauges.sf[gauges.sf$gauge_code == gg,], aes(label = sprintf("%0.0f km", coast_dist)), 
    hjust = 0.5, vjust = -0.5, family = font, size = textsize
  ) +
  labs(subtitle = "Distância até a costa")


# COMBINAR COVARIÁVEIS E PARÂMETROS --------------------------------------

cor.sf <- fit.df |> 
  left_join(
    gauges.sf |> select(gauge_code, elevation, coast_dist),
    by = "gauge_code"
  ) |> 
  select(gauge_code, xi0, alpha0, H, elevation, coast_dist) |> 
  print()

cor.long.sf <- cor.sf |> 
  pivot_longer(
    cols = !c(gauge_code, elevation, coast_dist),
    names_to = "parname",
    values_to = "par"
  ) |>
  print()

cor.summary.df <- cor.long.sf |> 
  group_by(parname) |> 
  reframe(
    cor_elevation = cor(par, elevation, use = "complete.obs", method = "pearson"),
    cor_coast_dist = cor(par, coast_dist, use = "complete.obs", method = "pearson")
  ) |> 
  print()

# xi0 ~ elevation
cor.sf |> 
  ggplot(aes(x = elevation, y = xi0)) +
  geom_point(shape = 21) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.3, alpha = 0.2) +
  labs(subtitle = sprintf("corr = %0.3f", cor.summary.df[cor.summary.df$parname == "xi0", "cor_elevation"]))

ggsave("figuras/cor_xi0_elevation.png", height = 6, width = 12, units = "cm", dpi = 300)

# alpha0 ~ elevation
cor.sf |> 
  ggplot(aes(x = elevation, y = alpha0)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.3, alpha = 0.2) +
  geom_point(shape = 21) + 
  labs(subtitle = sprintf("corr = %0.3f", cor.summary.df[cor.summary.df$parname == "alpha0", "cor_elevation"]))

ggsave("figuras/cor_alpha0_elevation.png", height = 6, width = 12, units = "cm", dpi = 300)

# H ~ elevation
cor.sf |> 
  ggplot(aes(x = elevation, y = H)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.3, alpha = 0.2) +
  geom_point(shape = 21) +
  labs(subtitle = sprintf("corr = %0.3f", cor.summary.df[cor.summary.df$parname == "H", "cor_elevation"]))

ggsave("figuras/cor_H_elevation.png", height = 6, width = 12, units = "cm", dpi = 300)

# xit ~ dist_coast
cor.sf |> 
  ggplot(aes(x = coast_dist, y = xi0)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.3, alpha = 0.2) +
  geom_point(shape = 21) +
  labs(subtitle = sprintf("corr = %0.3f", cor.summary.df[cor.summary.df$parname == "xi0", "cor_coast_dist"]))

ggsave("figuras/cor_xi0_coast_dist.png", height = 6, width = 12, units = "cm", dpi = 300)

# alpha0 ~ dist_coast
cor.sf |> 
  ggplot(aes(x = coast_dist, y = alpha0)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.3, alpha = 0.2) +
  geom_point(shape = 21) +
  labs(subtitle = sprintf("corr = %0.3f", cor.summary.df[cor.summary.df$parname == "alpha0", "cor_coast_dist"]))

ggsave("figuras/cor_alpha0_coast_dist.png", height = 6, width = 12, units = "cm", dpi = 300)
# H ~ dist_coast

cor.sf |> 
  ggplot(aes(x = coast_dist, y = H)) +
  geom_smooth(method = "lm", color = "black", linewidth = 0.3, alpha = 0.2) +
  geom_point(shape = 21) +
  labs(subtitle = sprintf("corr = %0.3f", cor.summary.df[cor.summary.df$parname == "H", "cor_coast_dist"]))

ggsave("figuras/cor_H_coast_dist.png", height = 6, width = 12, units = "cm", dpi = 300)


## MAPA PARÂMETROS -------------------------------------------------------

fit.sf <- gauges.sf |> 
  left_join(fit.df, by = "gauge_code")

# xi0
ggplot() +
  tidyterra::geom_spatraster(data = anadem, alpha = 0.3) +
  geom_sf(data = rs.sf, fill = NA, color = "black") +
  geom_sf(data = fit.sf, aes(color = xi0), size = 0.8) +
  scale_color_binned(palette = "magma", ) +
  scale_fill_viridis_c(
    name = "", na.value = NA, option = "cividis", # option = "turbo") +
    guide = guide_colorsteps(show.limits = TRUE, order = 2)
  ) +
  theme(
    panel.border = element_blank(),
    legend.position = c(0, 0),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    legend.key.height = unit(fontsize*0.9, "pt"),
    # legend.spacing.x = unit(fontsize*2, "pt")
    legend.spacing.y = unit(0, "pt")
  )

ggsave("figuras/map_xi0.png", height = 10, width = 10, units = "cm", dpi = 300)

# alpha0
ggplot() +
  tidyterra::geom_spatraster(data = anadem, alpha = 0.3) +
  geom_sf(data = rs.sf, fill = NA, color = "black") +
  geom_sf(data = fit.sf, aes(color = alpha0), size = 0.8) +
  scale_color_binned(palette = "magma", ) +
  scale_fill_viridis_c(
    name = "", na.value = NA, option = "cividis", # option = "turbo") +
    guide = guide_colorsteps(show.limits = TRUE, order = 2)
  ) +
  theme(
    panel.border = element_blank(),
    legend.position = c(0, 0),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    legend.key.height = unit(fontsize*0.9, "pt"),
    # legend.spacing.x = unit(fontsize*2, "pt")
    legend.spacing.y = unit(0, "pt")
  )

ggsave("figuras/map_alpha0.png", height = 10, width = 10, units = "cm", dpi = 300)

# kappa
ggplot() +
  tidyterra::geom_spatraster(data = anadem, alpha = 0.3) +
  geom_sf(data = rs.sf, fill = NA, color = "black") +
  geom_sf(data = fit.sf, aes(color = kappa), size = 0.8) +
  scale_color_binned(palette = "magma", ) +
  scale_fill_viridis_c(
    name = "", na.value = NA, option = "cividis", # option = "turbo") +
    guide = guide_colorsteps(show.limits = TRUE, order = 2)
  ) +
  theme(
    panel.border = element_blank(),
    legend.position = c(0, 0),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    legend.key.height = unit(fontsize*0.9, "pt"),
    # legend.spacing.x = unit(fontsize*2, "pt")
    legend.spacing.y = unit(0, "pt")
  )

ggsave("figuras/map_kappa.png", height = 10, width = 10, units = "cm", dpi = 300)

# H
ggplot() +
  tidyterra::geom_spatraster(data = anadem, alpha = 0.3) +
  geom_sf(data = rs.sf, fill = NA, color = "black") +
  geom_sf(data = fit.sf, aes(color = H), size = 0.8) +
  scale_color_binned(palette = "magma", ) +
  scale_fill_viridis_c(
    name = "", na.value = NA, option = "cividis", # option = "turbo") +
    guide = guide_colorsteps(show.limits = TRUE, order = 2)
  ) +
  theme(
    panel.border = element_blank(),
    legend.position = c(0, 0),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    legend.key.height = unit(fontsize*0.9, "pt"),
    # legend.spacing.x = unit(fontsize*2, "pt")
    legend.spacing.y = unit(0, "pt")
  )

ggsave("figuras/map_H.png", height = 10, width = 10, units = "cm", dpi = 300)

# theta
ggplot() +
  tidyterra::geom_spatraster(data = anadem, alpha = 0.3) +
  geom_sf(data = rs.sf, fill = NA, color = "black") +
  geom_sf(data = fit.sf, aes(color = theta), size = 0.8) +
  scale_color_binned(palette = "magma", ) +
  scale_fill_viridis_c(
    name = "", na.value = NA, option = "cividis", # option = "turbo") +
    guide = guide_colorsteps(show.limits = TRUE, order = 2)
  ) +
  theme(
    panel.border = element_blank(),
    legend.position = c(0, 0),
    legend.justification = c(0, 0),
    legend.direction = "vertical",
    legend.key.height = unit(fontsize*0.9, "pt"),
    # legend.spacing.x = unit(fontsize*2, "pt")
    legend.spacing.y = unit(0, "pt")
  )

ggsave("figuras/map_theta.png", height = 10, width = 10, units = "cm", dpi = 300)

