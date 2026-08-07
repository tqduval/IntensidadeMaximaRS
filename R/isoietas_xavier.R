
rm(list = ls()); invisible(gc())


# PACOTES ----------------------------------------------------------------

pacman::p_load(
    pacman, tidyverse, 
    geomtextpath, ggspatial,
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
nc.rs <- st_as_stars(nc.rs)
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

# Médias máximas anuais
nc.mean.max <- st_apply(nc.annual.max, MARGIN = c("x", "y"), FUN = mean) # para todo par {x,y}, aplique a média ao longo de "time"
nc.mean.sum <- st_apply(nc.annual.sum, MARGIN = c("x", "y"), FUN = mean)
nc.summary <- c(nc.mean.max, nc.mean.sum)
names(nc.summary) <- c("mean_annual_max", "mean_annual_sum")

nc.summary <- nc.summary |>
  st_warp(crs = st_crs(rs.sf)) |> 
  st_crop(rs.buffer.sf)

mean.max.df <- as.data.frame(nc.summary["mean_annual_max"])
mean.sum.df <- as.data.frame(nc.summary["mean_annual_sum"])

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

ggplot(data = mean.max.smooth.df, aes(x = x, y = y, z = mean)) +
  geom_contour_filled(bins = 8) +
  geom_sf(data = rs.sf, fill = NA, linewidth = 0.4, color = "black", inherit.aes = FALSE) +
  geom_textcontour(
    bins = 8, family = font, size = textsize, color = "black", linewidth = 0.3,
    vjust = 0.5, hjust = 0.2, halo = TRUE, upright = TRUE, straight = FALSE
  ) +
  coord_sf(crs = st_crs(rs.sf)) +
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
    labs(x = "Longitude", y = "Latitude")

gc()
