
# Plotar figura dos apêndices c/ IDFs em mais de um página
# usando `ggforce::facet_wrap_paginate()`

library(ggplot2)
library(ggforce)
library(tidyr)
library(dplyr)
library(purrr)

# Definições p/ figuras
ncol <- 3
nrow <- 9
textsize <- 5/.pt


# IDF SGB ----------------------------------------------------------------

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
  ggforce::facet_wrap_paginate(~gauge_label, ncol = ncol, nrow = nrow) +
  geom_line(aes(y = imax, color = "Intensidades observadas"), linewidth = 0.3) +
  geom_line(aes(y = q, color = "IDF-SGB"), linewidth = 0.1) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("k: %7.3f", k)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 1.8, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("γ: %7.3f", gamma)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 3.3, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("η: %7.3f", eta)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 4.8, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("Θ: %7.3f", theta)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 6.3, size = textsize
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

plot.seq <- seq_len(ggforce::n_pages(idf.sgb.plot))

path.apendices <- "figuras/preliminares/subdiarios/apendices_idf/"
if(!dir.exists(path.apendices)) dir.create(path.apendices)

walk(plot.seq, function(pg){
  idf.sgb.plot + ggforce::facet_wrap_paginate(~gauge_label, ncol = ncol, nrow = nrow, page = pg)
  f <- sprintf("%sidf_sgb_p%d.png", path.apendices, pg)
  ggsave(f, width = 15, height = 23, units = "cm", dpi = 600)
})


# IDF MLE ----------------------------------------------------------------

# Visualizar d-GEV
par.labels <- performance.dgev.df |> 
  filter(model == "d-GEV", prior == names(prior.names)[[1]]) |> # MLE
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
  geom_line(aes(y = q, color = "d-GEV (MLE)"), ~filter(., prior == names(prior.names)[1]), linewidth = 0.1) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("κ: %7.3f", kappa)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 1.8, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("H: %7.3f", H)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 3.3, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("θ: %7.3f", theta)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 4.8, size = textsize
  ) +
  scale_x_log10(name = "Duração [h]") +
  scale_y_log10(name = "Intensidade [mm/h]") +
  scale_color_manual(
    name = NULL,
    values = c("d-GEV (MLE)" = prior.col[[1]], "Intensidades observadas" = "grey90") # MLE
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

walk(plot.seq, function(pg){
  idf.dgev.plot + ggforce::facet_wrap_paginate(~gauge_label, ncol = ncol, nrow = nrow, page = pg)
  f <- sprintf("%sidf_mle_p%d.png", path.apendices, pg)
  ggsave(f, width = 15, height = 23, units = "cm", dpi = 600)
})


# IDF RGMLE --------------------------------------------------------------

# Visualizar d-GEV
par.labels <- performance.dgev.df |> 
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
  geom_line(aes(y = q, color = "d-GEV (RGMLE)"), ~filter(., prior == names(prior.names)[4]), linewidth = 0.1) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("κ: %7.3f", kappa)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 1.8, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("H: %7.3f", H)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 3.3, size = textsize
  ) +
  geom_text(
    aes(x = Inf, y = Inf, label = sprintf("θ: %7.3f", theta)), par.labels, inherit.aes = FALSE,
    family = font, hjust = 1.1, vjust = 4.8, size = textsize
  ) +
  scale_x_log10(name = "Duração [h]") +
  scale_y_log10(name = "Intensidade [mm/h]") +
  scale_color_manual(
    name = NULL,
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

walk(plot.seq, function(pg){
  idf.dgev.plot + ggforce::facet_wrap_paginate(~gauge_label, ncol = ncol, nrow = nrow, page = pg)
  f <- sprintf("%sidf_rgmle_p%d.png", path.apendices, pg)
  ggsave(f, width = 15, height = 23, units = "cm", dpi = 600)
})



# FUNÇÃO FIGURA IDF ------------------------------------------------------


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


# IDF QUANTIS ESTIMADOS RGMLE --------------------------------------------

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

plot.rp <- c(5, 20, 100)
colors <- c("#a53f00", "#FE6100", "#ffa56d")

# Plotar p/ todas as estações --- RGMLE
idf.rgmle.todas.plot <- idf.plot(
  df = idf.df, par = par.dgev.wide,
  plot.gauges = NULL,
  plot.model = "d-GEV",
  plot.prior = prior.names[4], # RGMLE
  prior.label = "RGMLE",
  plot.rp = plot.rp,
  tsize = textsize, par.label.y = 2,
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

walk(plot.seq, function(pg){
  idf.rgmle.todas.plot + ggforce::facet_wrap_paginate(~gauge_label, ncol = ncol, nrow = nrow, page = pg)
  f <- sprintf("%sidf_rgmle_rp_todas_p%d.png", path.apendices, pg)
  ggsave(f, width = 15, height = 23, units = "cm", dpi = 600)
})


# IDF QUANTIS ESTIMADOS MLE ----------------------------------------------

colors <- c("#183889", "#648FFF", "#a2bdff")

# Plotar p/ todas estações --- MLE
idf.mle.plot.todas <- idf.plot(
  df = idf.df, par = par.dgev.wide,
  plot.gauges = NULL,
  plot.model = "d-GEV",
  plot.prior = prior.names[1], # MLE
  prior.label = "MLE",
  plot.rp = plot.rp,
  tsize = textsize, par.label.y = 2,
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

walk(plot.seq, function(pg){
  idf.mle.plot.todas + ggforce::facet_wrap_paginate(~gauge_label, ncol = ncol, nrow = nrow, page = pg)
  f <- sprintf("%sidf_mle_rp_todas_p%d.png", path.apendices, pg)
  ggsave(f, width = 15, height = 23, units = "cm", dpi = 600)
})
