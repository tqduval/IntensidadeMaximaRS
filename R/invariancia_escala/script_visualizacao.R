# Este script gera figuras para o artigo baseado nos resultados gerados
# em 'script_simple_scaling', que gera dados mais simples, sem análise
# de incerteza com intervalos de confiança.

rm(list = ls()); invisible(gc())

# PACOTES -----------------------------------------------------------------

if(!require(pacman)) install.packages("pacman")
pacman::p_load(pacman, tidyverse, arrow, patchwork, lmom, ggtext)

font <- "JetBrainsMono NF"
fontsize <- 8
textsize <- (fontsize - 2)/.pt

ggplot2::set_theme(
  theme_bw() +
    theme(
      legend.position = "bottom",
      legend.background = element_rect(fill = NA),
      legend.box.margin = margin(-10, 2, 2, 2),
      # legend.box.margin = margin(-5, 2, 2, 2),
      # legend.key.size = unit(fontsize, "pt"),
      legend.key.height = unit(fontsize, "pt"),
      text = element_text(color = "black", family = font, size = fontsize),
      axis.text = element_text(color = "black", family = font),
      strip.text = element_text(color = "black", margin = margin(2, 2, 2, 2), size = fontsize - 2),
      strip.background = element_blank(),
      panel.border = element_rect(color = "black"),
      axis.ticks = element_line(linewidth = 0.3),
      panel.grid = element_blank(),
      plot.margin = margin(),
      # panel.grid = element_line(linewidth = 0.3)
    )
)


# LER DADOS ---------------------------------------------------------------

# Informações
info.df <- arrow::read_parquet("base/df_subdaily_info.pqt") |> 
  filter(responsible != "SGB")

# Máximos anuais
imax.df <- arrow::read_parquet("base/df_imax.pqt")

# Resultados gerados
mom.df <- arrow::read_parquet("base/df_scale_moments.pqt")       # momentos não centrais
scale.df <- arrow::read_parquet("base/df_scale_coefficient.pqt") # expoentes de escala e ICs bootstrap
coefs.df <- arrow::read_parquet("base/df_scale_block_boot.pqt")  # contém coeficientes de desagregação (offset == 0)

# Resultados gerados p/ estação de validação separadamente
mom.val.df <- arrow::read_parquet("base/df_scale_moments_val.pqt")
scale.val.df <- arrow::read_parquet("base/df_scale_coefficient_val.pqt")
coefs.val.df <- arrow::read_parquet("base/df_scale_block_boot_val.pqt")

n_distinct(mom.df$gauge_code) # there should be 25 gauges with at least 10 years

# PREPARAR DADOS ----------------------------------------------------------

# Avaliar estações com mais anos e durações e selecionar uma
gauges <- mom.df |> 
  group_by(gauge_code, which_moment) |> 
  reframe(min_d = min(d),
          max_d = max(d),
          n_durations = n_distinct(d),
          n_year = first(n_year),
          duration_years = n_durations*n_year) |> 
  arrange(desc(n_durations), desc(duration_years)) |> 
  pull(gauge_code) |> 
  unique()

# df.mom %>% 
#   slice_head(by = gauge_code) %>% 
#   arrange(desc(n_year))

# gauge <- "430350901A"
# gauge <- gauges[1]
gauge <- "A801"

# Momentos e expoente de escala local
# Extrair os limites dos intervalos
groups <- unique(scale.df$d)
# group.names <- c("All durations", "Sub-hourly", "Hourly", "Daily")
group.names <- "Hourly"
groups <- setNames(group.names, groups)

# Mesclar df.mom e df.scale e filtrar p/ uma estação
mom.gauge.df <- mom.val.df |> 
  # filter(which_moment == 1, d >= 1, d <= 24) |> # filtrar estação de interesse
  filter(which_moment == 1, d >= 1/24, d <= 1) |> # filtrar estação de interesse
  left_join(scale.val.df |> 
              filter(gauge_code == gauge) %>% 
              mutate(lower = as.numeric(sub("-.*", "", d)),
                     upper = as.numeric(sub(".*-", "", d)),
                     group = factor(recode(d, !!!groups), levels = group.names)) |> 
              select(gauge_code, which_moment, group, scale, rsquared, lower, upper),
            by = c("gauge_code", "which_moment"),
            relationship = "many-to-many") |> 
  filter(group == "Hourly") # garante que as durações estejam compatíveis c/ grupos

unique(mom.df$gauge_code)

mom.reg.df <- mom.df |> 
  filter(d >= 1/24, d <= 1, which_moment == 1) |> 
  group_by(d, which_moment) |> 
  reframe(
    mom = mean(mom)
  ) |> 
  mutate(gauge_code = "Regional") |> 
  left_join(
    y = scale.df |> 
      # filter(d == groups, which_moment == 1, gauge_code == "Regional") |> 
      filter(which_moment == 1, gauge_code == "Regional") |> 
      select(!d),
    by = c("gauge_code", "which_moment")
  )

mom.plot.df <- bind_rows(mom.gauge.df, mom.reg.df) |> 
  select(!group) |> 
  mutate(
    gauge_code = ifelse(gauge_code != "Regional", sprintf("Local (%s)", gauge), gauge_code),
    label = sprintf("%s: H = %0.4f (R² = %0.4f)", gauge_code, scale, rsquared)
  )


# FIGURAS INVARIÂNCIA DE ESCALA ------------------------------------------

labels <- unique(mom.plot.df$label)
colors <- setNames(c("black", "grey50"), labels)
linetypes <- setNames(c("solid", "dashed"), labels)

# Figure a) invariância de escala
plot.a <- mom.plot.df |> 
  ggplot(aes(x = d, y = mom)) +
  geom_smooth(
    aes(color = label, linetype = label),
    se = FALSE, linewidth = 0.3, method = "lm", formula = y ~ x
  ) +
  geom_smooth(
    aes(color = label, linetype = label), 
    se = FALSE, linewidth = 0.3, method = "lm", formula = y ~ x
  ) +
  geom_point(aes(color = label), ~filter(.,gauge_code != "Regional"), shape = 16, size = 0.8) +
  # geom_vline(xintercept = 6) +
  scale_x_log10() +
  scale_y_log10() +
  scale_color_manual(name = "", values = colors) +
  scale_linetype_manual(name = "", values = linetypes) +
  labs(x = "d/D", y = "E[I<sub>d</sub><sup>q</sup>]", subtitle = "(a)") +
  theme(
    axis.title.y = ggtext::element_markdown(),
    legend.position = c(1, 1.06),
    legend.justification = c(1,1),
    legend.key.spacing.y = unit(1, "pt")
  ); plot.a

# Figure b) histograma de R²
plot.b <- scale.df |> 
  filter(
    gauge_code != "Regional",
    which_moment == 1,
    # d == names(groups)[3]
  ) |> 
  ggplot(aes(x = rsquared)) +
  # facet_wrap(~which_moment, nrow = n_distinct(df.scale$which_moment), labeller = as_labeller(function(x) paste0("q = ", x))) +
  geom_histogram(
    aes(y = after_stat(count/sum(count))), bins = 20, 
    fill = "grey90", color = "black", linewidth = 0.3
  ) +
  scale_y_continuous(labels = scales::percent) +
  scale_x_continuous(labels = scales:::label_number(drop0trailing = TRUE)) +
  labs(x = "R²", y = "", subtitle = "(b)") +
  theme(
    legend.position =  "none",
    axis.title.y = ggtext::element_markdown()
  ); plot.b

# Juntar duas figuras
plot.scale.inv <- 
  plot.a + plot.b +
  plot_layout(
    widths = c(0.65, 0.35)
  ); plot.scale.inv

ggsave("artigo/figuras/plot_scale_invariance_hourly.png", plot.scale.inv, width = 14, height = 6, units = "cm", dpi = 300)


# FIGURAS COEFICIENTES DE DESAGREGAÇÃO -----------------------------------

# Coeficientes de desagregação da CETESB
df.cetesb <- data.frame(relacao = c("5 min/30 min", "10 min/30 min", "15 min/30 min", "20 min/30 min", "25 min/30 min", 
                                    "30 min/1 h", "1 h/24 h", "6 h/24 h", "8 h/24 h", "10 h/24 h", "12 h/24 h"),
                        d1 = c(seq(5, 30, 5)/60, 1, 6, 8, 10, 12),
                        d2 = c(rep(0.5, 5), 1, rep(24, 5)),
                        d_target = factor(c(seq(5, 30, 5)/60, 1, 6, 8, 10, 12)),
                        coef_pdmax = c(0.34, 0.54, 0.7, 0.81, 0.91, 0.74, 0.42, 0.72, 0.78, 0.82, 0.85))

df.cetesb$relacao_mod <- as.character(NA)
df.cetesb$coef_cetesb <- as.numeric(NA)

for(i in 1:nrow(df.cetesb)){
  
  if(df.cetesb$d2[i] <= 1){
    
    df.cetesb$relacao_mod[i] <- paste0(df.cetesb$d1[i]*60, " min/1 h")
    df.cetesb$coef_cetesb[i] <- df.cetesb$coef_pdmax[i] * df.cetesb$coef_pdmax[6]
    
  } else{
    
    df.cetesb$relacao_mod[i] <- paste0(df.cetesb$d1[i], " h/24 h")
    df.cetesb$coef_cetesb[i] <- df.cetesb$coef_pdmax[i]
    
  }
}

# Visualizar coeficientes de desagregação
# Definir ordem de plotagem --- ordernar por expoente de escala H
order <- coefs.df |> 
  filter(
    gauge_code %in% gauges,    # df.coefs was run with gauges with at least 8 years
    # d == names(groups)[3],     # filter hourly run
    d == "1-24",     # filter hourly run
    # which_variable == "scale", # order based on scale invariance H
    which_variable == "coef", # order based on disagg coef
    which_moment == 1
  ) |> 
  arrange(statistic) |> 
  pull(gauge_code) |> 
  unique()

filter.target <- c(1, 6, 12)

local.df <- coefs.df |> 
  filter(
    gauge_code != "Regional",
    gauge_code %in% gauges,
    # d == names(groups)[3],
    d == "1-24",
    which_variable == "coef",
    which_moment == 1,
    d_target %in% filter.target
  ) |> 
  mutate(gauge_code = factor(gauge_code, levels = order))

regional.df <- coefs.df |> 
  filter(
    gauge_code == "Regional",
    which_variable == "coef",
    # d == names(groups)[3],
    d == "1-24",
    which_moment == 1,
    d_target %in% filter.target
  ) |> 
  select(!gauge_code)

colors.coef <- setNames(c("black", "grey50", "grey20"), c("Local", "Regional CI<sub>(95%)</sub>", "CETESB"))

plot.coef.disag <- ggplot() +
  facet_wrap(
    ~d_target, 
    ncol = 3,
    labeller = as_labeller(function(x) sprintf("d/D = %s/24 h", x))
  ) +
  geom_rect(
    aes(xmin = 0, xmax = Inf, ymin = ci_lower, ymax = ci_upper, fill = names(colors.coef)[2]), 
    regional.df,
    color = "transparent",
    alpha = 0.4
  ) +
  geom_hline(
    aes(yintercept = coef_cetesb, color = names(colors.coef)[3]),
    df.cetesb |> filter(d_target %in% filter.target),
    linetype = "dashed", linewidth = 0.3
  ) +
  geom_errorbar(
    aes(x = gauge_code, ymin = ci_lower, ymax = ci_upper, color = names(colors.coef)[1]),
    local.df,
    width = 0, linewidth = 0.4
  ) +
  geom_point(
    aes(x = gauge_code, y = statistic, color = "Local"),
    local.df,
    size = 0.6
  ) +
  scale_color_manual(
    name = "",
    values = colors.coef,
    aesthetics = c("color", "fill")
  ) +
  labs(x = "", y = "C<sub>P</sub>") +
  theme(
    # legend.position = "right",
    # legend.margin = margin(l = -5),
    # legend.key.size = unit(fontsize + 2, units = "pt"),
    legend.text = ggtext::element_markdown(),
    axis.title.y = ggtext::element_markdown(),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5, size = fontsize*0.7)
  ); plot.coef.disag

ggsave("artigo/figuras/plot_coef_disag_hourly.png", plot.coef.disag, width = 14, height = 6, units = "cm", dpi = 300)



# DESAGREGAR SÉRIE ORIGINAL ----------------------------------------------

# Coeficientes de desagregação regionais (s/ estação de validação)
coef.hour.reg <- coefs.df |> 
  filter(
    gauge_code == "Regional",
    which_variable == "coef",
    # d == names(groups)[3],
    d == "1-24",
    which_moment == 1,
    d_target %in% filter.target
  ) |> 
  select(d_target, statistic, ci_lower, ci_upper)

# Extrair série de `imax` de d = (1, 6, 12, 24) h da estação de validação
val.data <- imax.df |> 
  filter(
    gauge_code == gauge,
    na_prct_yr <= 0.2,
    d %in% c(filter.target, 24)
    # d == 24
  )

# Extrair série de `pmax` que será desagregada pelos coeficientes
val.ams <- val.data |> 
  filter(d == 24) |> # extrair série de 24 h
  pull(imax)*24      # converter intensidade em precipitação



## DESAREGAR  C/ COEF. REGIONAL ------------------------------------------

# Desagregar série de validação com coeficientes regionais
coef.reg.ls <- split(coef.hour.reg, coef.hour.reg$d_target) # separar os coeficientes tabelas diferentes p/ cada duração alvo
disag.ls <- lapply(coef.reg.ls, function(coef){             # multiplicar cada série pelo seu respectivo coeficiente
  coef <- coef[-1]             # primeira coluna é `d_target`, sobram somente `statistic`, `ci_lower` e `ci_upper`
  disag <- t(coef) %*% val.ams # multiplicar os coeficientes (1x3) pela série (nx1): tabela nx3
  return(as.data.frame(t(disag)))
})

disag.df <- bind_rows(disag.ls, .id = "d") |> 
  mutate(
    wateryear = rep(unique(val.data$wateryear), length(filter.target)),
    d = as.numeric(d),
    disag = statistic
  ) |> 
  left_join(
    val.data |> 
      filter(d %in% filter.target) |> 
      mutate(pmax = imax*d) |> 
      select(d, wateryear, pmax),
    by = c("d", "wateryear")
  ) |> 
  select(d, wateryear, pmax, disag, ci_lower, ci_upper)

disag.fill.df <- disag.df |> 
  complete(d, wateryear = full_seq(wateryear, 1))

colors.disag <- setNames(c("black", "grey20"), c("Observed", "Disaggregated"))
fills.disag <- setNames(c(NA, "grey40"), c("Observed", "Disaggregated"))
linetypes.disag <- setNames(c("solid", "solid"), c("Observed", "Disaggregated"))

plot.disag <- disag.fill.df |> 
  ggplot(aes(x = wateryear)) +
  facet_wrap(~factor(d), scale = "free_y", labeller = as_labeller(function(x) paste0(x, "/24 h"))) +
  geom_ribbon(aes(ymin = ci_lower, ymax = ci_upper, fill = "Disaggregated"), alpha = 0.5) +
  geom_line(aes(y = disag, color = "Disaggregated", linetype = "Disaggregated"), alpha = 0.7, linewidth = 0.3) +
  geom_line(aes(y = pmax, color = "Observed", linetype = "Observed"), linewidth = 0.4) +
  geom_point(aes(y = pmax, color = "Observed"), shape = 21, fill = "white", size = 0.8) +
  scale_color_manual(name = "", values = colors.disag) +
  scale_fill_manual(name = "", values = fills.disag) +
  scale_linetype_manual(name = "", values = linetypes.disag) +
  labs(x = "", y = "Rainfall depths [mm]"); plot.disag


ggsave("artigo/figuras/plot_disag.png", plot.disag, width = 14, height = 6, units = "cm", dpi = 300)

# Estimate l-moments and l-moment ratios for disaggregated series
disag.df[disag.df$d == 1, "pmax"] |> lmom::samlmu()

disag.df |> 
  pivot_longer(
    cols = c(pmax, disag),
    names_to = "which_series",
    values_to = "value"
  ) |> 
  group_by(d, which_series) |> 
  reframe(as_tibble_row(samlmu(value))) |> 
  print()



# DESAGREGAR QUANTIS ESTIMADOS -------------------------------------------

# Fit GEV with L-moments to validation 24 h `pmax` series
val.lmom <- lmom::samlmu(val.ams)
val.gevfit <- lmom::pelgev(val.lmom)

# Empirical quantiles --- Weibull plotting positions
val.ams <- sort(val.ams, decreasing = TRUE)
idx <- seq_along(val.ams)
n <- length(val.ams)
plot.pos <- idx/(n + 1)
rp.obs <- 1/plot.pos
rp.extra <- seq(2, 100, 2)
return.period <- sort(c(rp.obs, rp.extra))

# Fit IDFs to quantiles estimated from the original target duration series
val.ls <- split(val.data, val.data$d)
val.idf.ls <- lapply(val.ls, function(val.df){

  d.ams <- unique(val.df$d)
  val.ams <- sort(val.df$imax*d.ams, decreasing = TRUE)
  idx <- seq_along(val.ams)
  n <- length(val.ams)
  plot.pos <- idx/(n + 1)
  rp.obs <- 1/plot.pos
  return.period <- sort(c(rp.obs, rp.extra))

  val.lmom <- lmom::samlmu(val.ams)
  val.gevfit <- lmom::pelgev(val.lmom)
  val.qgev <- lmom::quagev(f =  1 - 1/return.period, para = val.gevfit)

  res <- data.frame(
    which = "Val.",
    d = d.ams,
    return_period = return.period,
    q = val.qgev
  )

})

val.idf <- bind_rows(val.idf.ls)

# Desagregar com coeficientes de invariância de escala
disag.inv.idf.ls <- lapply(coef.reg.ls, function(coef){
  coef <- coef[-1]
  disag <- t(coef) %*% val.idf$q[val.idf$d == 24]
  return(as.data.frame(t(disag)))
})

disag.inv.idf <- bind_rows(disag.inv.idf.ls, .id = "d") |> 
  mutate(
    d = as.numeric(d),
    which = "Disagg. (R,SS)",
    return_period = val.idf[val.idf$d != 24, "return_period"]
  ) |> 
  rename(q = statistic) |>
  select(which, d, return_period, q,, ci_lower, ci_upper)

# Desagregar com coeficientes da CETESB
coef.cetesb.idf <- df.cetesb |> 
  select(d_target, coef_cetesb) |> 
  filter(d_target %in% filter.target) |> 
  rename(d = d_target, coef = coef_cetesb)

coef.cetesb.ls <- split(coef.cetesb.idf, coef.cetesb.idf$d)
disag.cetesb.idf.ls <- lapply(coef.cetesb.ls, function(coef){
  coef <- coef[-1]
  disag <- t(coef) %*% val.idf$q[val.idf$d == 24]
  return(as.data.frame(t(disag)))
})

disag.cetesb.idf <- bind_rows(disag.cetesb.idf.ls, .id = "d") |> 
  mutate(
    which = "Disagg. (CETESB)",
    d = as.numeric(d),
    return_period = val.idf[val.idf$d != 24, "return_period"]
  ) |> 
  rename(q = coef) |> 
  select(which, d, return_period, q)

# Desagregar com os coeficientes da propria estação de validação
coef.val.idf <- coefs.val.df |> 
  filter(
    d == "1-24",
    d_target %in% filter.target,
    which_variable == "coef",
    which_moment == 1
  ) |> 
  select(d_target, statistic, ci_lower, ci_upper)

coef.val.ls <- split(coef.val.idf, coef.val.idf$d_target)
disag.val.idf.ls <- lapply(coef.val.ls, function(coef){
  coef <- coef[-1]
  disag <- t(coef) %*% val.idf$q[val.idf$d == 24]
  return(as.data.frame(t(disag)))
})

disag.val.idf <- bind_rows(disag.val.idf.ls, .id = "d") |> 
  mutate(
    which = "Disagg. (VAL,SS)",
    d = as.numeric(d),
    return_period = val.idf[val.idf$d != 24, "return_period"]
  ) |> 
  rename(q = statistic) |> 
  select(which, d, return_period, q, ci_lower, ci_upper)

# Juntar todas as desagregações
idf <- bind_rows(
  val.idf[val.idf$d != 24,], # quantis estimados com as séries originais da estação de validação
  disag.inv.idf,             # quantis desagregados com coeficientes regionais
  disag.val.idf,             # quantis estimados com coeficientes da estação de validação
  disag.cetesb.idf           # quantis desagregados com coeficientes da CETESB
) |> 
  mutate( # transformar de precipitação p/ intensidade
    q = q/d,
    ci_lower = ci_lower/d,
    ci_upper = ci_upper/d
  )

# Plotting positions for empirical validation gauge
plot.pos <- val.data |> 
  select(gauge_code, d, imax) |> 
  filter(d %in% filter.target) |> 
  mutate(d = factor(d)) |> 
  group_by(d) |> 
  arrange(desc(imax)) |> 
  mutate(
    i = row_number(),
    p = i/(n() + 1),
    return_period = 1/p
  )

which <- unique(idf$which) |> print()
# which <- which[c(1, 4, 2, 3)]

plot.disag.coef <- 
  bind_rows(
    coef.val.idf |> 
      rename(coef = statistic) |> 
      select(d_target, coef) |> 
      mutate(
        d_target = factor(d_target),
        which = "VAL,SS"
      ),
    coef.hour.reg |> 
      rename(coef = statistic) |> 
      select(d_target, coef) |> 
      mutate(
        d_target = factor(d_target),
        which = "R,SS"
      ),
    coef.cetesb.idf |> 
      rename(d_target = d) |> 
      mutate(
        d_target = factor(d_target),
        which = "CETESB"
      ),
  )

abbrev <- c("VAL,SS" = "VAL", "R,SS" = "R", "CETESB" = "CETESB")
label.order <- c("VAL,SS", "R,SS", "CETESB")

plot.coef.labels <- plot.disag.coef |> 
  mutate(
    which = factor(which, levels = label.order),
    text = sprintf("C<sub>P,%s</sub>: %0.3f", abbrev[as.character(which)], coef)
  ) |> 
  arrange(d_target, which) |>
  group_by(d_target) |> 
  summarise(label = paste(text, collapse = "<br>")) |>
  mutate(d = as.numeric(as.character(d_target)))


colors.idf <- setNames(c("black", "grey30", "grey60", "black"), which)[-1]
fills.idf <- setNames(c(NA, "grey30", "grey70", NA), which)[-1]
linetypes.idf <- setNames(c("solid", "dashed", "dashed", "dotted"), which)[-1]

plot.disag.idf <- idf |> 
  filter(which != "Val.") |> 
  ggplot(aes(x = return_period)) +
  facet_wrap(~factor(d), scale = "free_y", labeller = as_labeller(function(x) sprintf("d/D = %s/24 h", x))) +
  geom_ribbon(
    aes(ymin = ci_lower, ymax = ci_upper, fill = which), 
    ~filter(., which %in% c("Disagg. (R,SS)", "Disagg. (VAL,SS)")),
    color = "transparent", alpha = 0.15
  ) +
  geom_line(aes(y = q, color = which, linetype = which)) +
  geom_point(aes(y = imax, shape = "Observed"), plot.pos, size = 0.6) +
  geom_richtext(
    aes(x = 1.05, y = Inf, label = label), plot.coef.labels,
    inherit.aes = FALSE, hjust = 0, vjust = 1.1, fill = NA, label.color = NA,
    label.padding = unit(c(2,2,2,2), "pt"), size = textsize, family = font
  ) +
  # geom_richtext(
  #   aes(x = Inf, y = 20, label = label), plot.coef.labels,
  #   inherit.aes = FALSE, hjust = 1.1, vjust = 0, fill = NA, label.color = NA,
  #   label.padding = unit(c(2,2,2,2), "pt"), size = textsize, family = font
  # ) +
  scale_x_log10() +
  scale_color_manual(name = "", values = colors.idf, breaks = which[-1], limits = which[-1]) +
  scale_fill_manual(name = "", values = fills.idf, breaks = which[-1], limits = which[-1]) +
  scale_linetype_manual(name = "", values = linetypes.idf, breaks = which[-1], limits = which[-1]) +
  scale_shape_manual(name = "", values = c("Observed" = 9)) +
  labs(x = "T [years]", y = "Intensity [mm/h]") +
  theme(
    legend.spacing.x = unit(-7.5, "pt"),
  ); plot.disag.idf

ggsave("artigo/figuras/plot_disag_idf.png", plot.disag.idf, width = 16, height = 6, units = "cm", dpi = 300)
