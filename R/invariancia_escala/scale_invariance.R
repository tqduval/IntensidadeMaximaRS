
rm(list = ls()); invisible(gc())

# PACOTES -----------------------------------------------------------------

if(!require(pacman)) install.packages("pacman")
pacman::p_load(pacman, tidyverse, arrow, lmom, ggpmisc)

# suppressMessages(extrafont::loadfonts())
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

# LER DADOS ---------------------------------------------------------------

imax.df.prefilter <- arrow::read_parquet(file = "base/gerados/df_imax.pqt")

na.accept <- 0.2    # percentual de falhas (0.2 -> 20%)
min.years <- 8      # mínimo de anos na série

imax.df <- imax.df.prefilter |> 
  filter(na_prct_yr <= na.accept) |> 
  mutate(
    n_years = n_distinct(wateryear),
    gauge_label = sprintf("%s | %2g anos", gauge_code, n_years),
    .by = "gauge_code"
  ) |> 
  filter(
    n_years >= min.years
  )

# Identificadores p/ estações
gauge.ids <- imax.df |> 
  slice_head(n = 1, by = "gauge_code") |> 
  select(gauge_code, n_years) |> 
  arrange(desc(n_years)) |>
  mutate(id = row_number()) |> print()

gauges <- unique(imax.df$gauge_code)

info.df <- arrow::read_parquet("base/fonte/subdaily/df_subdaily_info.pqt")
info.df <- info.df |> 
  filter(gauge_code %in% gauges)


# FUNÇÕES -----------------------------------------------------------------

source("R/funcoes/fun_scale_invariance.R")


# RODAR FUNÇÃO ------------------------------------------------------------

# Argumentos
which.moment <- 1:3 # qual momento usar de base p/ gerar figuras individuais
min.duration <- 3   # mínimo de durações diferentes necessário   

imax.df.prefilter

# Avaliar diferentes modelos
duration.intervals <- list(
  c(10/60, 240), # todas as durações
  c(10, 60)/60,  # sub-horários
  c(1, 24),      # horários
  c(1, 10)*24    # diários
)     

# Regressão + coeficiente de escala + plots
scale.invariance <- lapply(X = duration.intervals, FUN = function(interval){
  
  res <- fit.simple.scaling(
    df.imax = imax.df.prefilter,
    na.accept = na.accept,
    min.years = min.years,
    which.moment = which.moment,
    which.duration = interval,
    min.duration = 3,
    font.family = font,
    font.size = 12,
    plot.dim = c(5,4)
  )
  
}); names(scale.invariance) <- c("all", "subhourly", "hourly", "daily")

# Extrair os momentos
scale.moments <- bind_rows(lapply(scale.invariance, `[[`, "scale.moments")) |> 
  slice(1, .by = c(gauge_code, d, which_moment)) # só precisa de 1 linha de cada duração por estação

# Extrair estimativas de expoente de escala
scale.coefficient <- bind_rows(lapply(scale.invariance, `[[`, "scale.coefficient"))


# MOMENTOS P/ UMA ESTAÇÃO -------------------------------------------------

# Avaliar estações com mais anos e durações e selecionar uma
gauges <- scale.moments %>%
  group_by(gauge_code, which_moment) %>%
  reframe(min_d = min(d),
          max_d = max(d),
          n_durations = n_distinct(d),
          n_year = first(n_year),
          duration_years = n_durations*n_year) %>% 
  arrange(desc(n_durations), desc(duration_years)) %>% 
  pull(gauge_code) %>% 
  unique()

gauge <- gauges[1]
# gauge <- "A826"

duration.intervals <- unique(scale.coefficient$d)
unique(scale.moments$gauge_code)
unique(scale.moments$which_moment)
unique(scale.coefficient$which_moment)

sort(unique(scale.moments$d))

# Momentos e expoente de escala local
df.mom.gauge <- merge(
    x = scale.moments |> filter(gauge_code == gauge),
    y = scale.coefficient |> 
      filter(gauge_code == gauge, d == duration.intervals[1]) |> 
      select(gauge_code, which_moment, scale, rsquared),
      by = c("gauge_code", "which_moment"),
      all.x = FALSE
  ) |> 
  mutate(
    label = sprintf("q = %g | H = %0.4f | R² = %0.4f", which_moment, scale, rsquared),
    duration_ratio = d/max(d)
  )

vlines <- tibble(
  x = c(1, 24)/240,
  label = sprintf("d = %s h", c(1, 24))
)

# Figura momentos p/ uma estação
plot.scale.gauge <- 
  ggplot(data = df.mom.gauge, aes(x = duration_ratio, y = mom)) +
  geom_vline(aes(xintercept = x), vlines, linewidth = 0.2, linetype = "dashed", color = "grey80") +
  geom_text(
    aes(label = label, x = x, y = 0), vlines, family = font, size = textsize, inherit.aes = FALSE,
     vjust = -1, hjust = 1.15
  ) +
  geom_smooth(
    aes(group = which_moment), color = "grey20", linetype = "solid", alpha = 0.8,
    method = "lm", formula = y~x, se = FALSE, linewidth = 0.3
  ) +
  geom_point(aes(shape = label), color = "black", fill = "black", size = 1) +
  scale_x_log10() +
  scale_y_log10(labels = scales::label_number(big.mark = ".")) +
  scale_shape_manual(name = NULL, values = c(21, 24, 22)) +
  labs(
    y = "E[I<sub>d</sub><sup>q</sup>]", x = "d/D",
    # caption = paste0("Estação: ", gauge, " | ", df.mom.gauge$n_year[1], " anos")
  ) +
  theme(
    legend.position =  c(1, 0.99),
    legend.justification = c(1, 1),
    axis.title.y = ggtext::element_markdown()
  ) +
  guides(
    shape = guide_legend(override.aes = list(alpha = 1, color = "black", fill = "black", size = 1.2))
  ); plot.scale.gauge

# ggsave(
#   "figuras/preliminares/subdiarios/invariancia_escala/plot_scale_gauge.png", plot.scale.gauge,
#   width = 15, height = 12, units = "cm", dpi = 300
# )

moment.label <- tibble(
  which_moment = factor(which.moment, which.moment),
  label = paste("q =", which_moment)
)

# Figura histograma R²
plot.r2.hist <- scale.coefficient |> 
  filter(gauge_code != "Regional", d == duration.intervals[1]) |> 
  ggplot(aes(x = rsquared)) +
  facet_wrap(~which_moment, nrow = length(which.moment)) +
  geom_histogram(
    aes(y = after_stat(count/sum(count))),
    bins = 35, fill = "grey90", color = "black", linewidth = 0.2
  ) +
  geom_text(
    aes(label = label, x = Inf, y = Inf), moment.label, family = font, size = textsize,
    inherit.aes = FALSE, hjust = 1.2, vjust = 1.6
  ) +
  scale_y_continuous(labels = scales::label_percent()) +
  labs(x = "R²", y = NULL) +
  theme(
    legend.position =  "none",
    strip.text = element_blank(),
    axis.title.y = ggtext::element_markdown(),
    plot.margin = margin(0.5, 0.5, 0.5, 5.5)
  ); plot.r2.hist

# ggsave(
#   "figuras/preliminares/subdiarios/invariancia_escala/plot_r2_hist.png", plot.r2.hist,
#   width = 8, height = 12, units = "cm", dpi = 300
# )

plot.scale.r2 <- 
  wrap_plots(plot.scale.gauge + labs(subtitle = "(a)")) +
  wrap_plots(plot.r2.hist + labs(subtitle = "(b)")) +
  plot_layout(widths = c(1, 0.5), guides = "keep"); plot.scale.r2

ggsave(
  "figuras/preliminares/subdiarios/invariancia_escala/plot_scale_r2_all.png", plot.scale.r2,
  width = 15, height = 10, units = "cm", dpi = 300
)

# Informações complementares
scale.coefficient |> 
  filter(gauge_code != "Regional", d == duration.intervals[1]) |> 
  arrange(rsquared)

df.mom.gauge$scale |> 
  unique() |> 
  mean()

# INTERVALOS DE CONFIANÇA LOCAIS ------------------------------------------

df.scale.boot <- arrow::read_parquet("base/gerados/scale_invariance/df_scale_block_boot.pqt")
unique(df.scale.boot$d_target)

groups <- unique(df.scale.boot$d)
group.names <- c("Todas", "Sub-horários", "Horários", "Diários")
groups <- setNames(object = group.names, nm = groups)

df.scale.boot <- df.scale.boot |> 
  mutate(d = groups[d])

# Filtrar expoentes de escala 'scale'
df.scale <- df.scale.boot |> 
  filter(
    which_variable == "scale",
    which_moment == 1,
    d != "Todas"
  ) |> 
  select(gauge_code, n_year, d, statistic, ci_lower, ci_upper) |> 
  left_join(
    gauge.ids |> select(gauge_code, id),
    by = "gauge_code"
  ) |> 
  arrange(gauge_code)

colors <- c("Todas" = "black", "Sub-horários" = "orange", "Horários" = "steelblue2", "Diários" =  "green3")

dodge.width <- 0.5

order.gg <- df.scale %>% 
  filter(d == "Horários") %>% 
  arrange(statistic) %>% 
  pull(gauge_code) %>% 
  unique()

order.id <- df.scale %>% 
  filter(d == "Horários") %>% 
  arrange(statistic) %>% 
  pull(id) %>% 
  unique()

# Figura intervalos de confiança locais c/ bootstrap por blocos
plot.scale.ci.block.boot <- df.scale |> 
  # filter(d != "Sub-horários") |> 
  filter(gauge_code != "Regional") |> 
  mutate(
    gauge_code = factor(gauge_code, levels = order.gg),
    id = factor(id, levels = order.id),
    d = factor(d, groups)
  ) |> 
  # mutate(gauge_code = fct_reorder(gauge_code, statistic)) |> 
  ggplot(mapping = aes(y = statistic, ymin = ci_lower, ymax = ci_upper, color = d)) +
  geom_errorbar(aes(x = id), width = 0, linewidth = 0.4, alpha = 0.5, position = position_dodge(width = dodge.width)) +
  ylim(min(df.scale$ci_lower), max(df.scale$ci_upper)) +
  geom_point(aes(x = id), shape = 16, size = 0.6, alpha = 0.9, position = position_dodge(dodge.width)) +
  labs(x = "ID", y = "H", color = "Faixa de duração") +
  scale_color_manual(
    name = NULL, values = colors, aesthetics = c("color", "fill")
  ) + 
  theme(
    legend.position = c(0, -0.02),
    legend.justification = c(0, 0),
    legend.key.width = unit(12, "pt"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  ); plot.scale.ci.block.boot

# Figura intervalos de confiança locais c/ bootstrap por blocos + faixas regional
plot.scale.boot.regional <- plot.scale.ci.block.boot +
  geom_hline(
    aes(yintercept = statistic, color = d), df.scale |> filter(gauge_code == "Regional"),
    linewidth = 0.3, linetype = "dashed"
  ) +
  geom_rect(
    aes(xmin = -Inf, xmax = Inf, ymin = ci_lower, ymax = ci_upper, fill = d),
    df.scale |> filter(gauge_code == "Regional"), 
    alpha = 0.1, color = "transparent", linewidth = 0,
    show.legend = FALSE
  ); plot.scale.boot.regional

ggsave(
  "figuras/preliminares/subdiarios/invariancia_escala/plot_scale_ci_regional.png", plot.scale.boot.regional,
  width = 15, height = 10, units = "cm", dpi = 300, bg = "white"
)

# Figuras por faixa:
df.scale |> 
  filter(gauge_code != "Regional", d %in% c("Sub-horários", "Horários", "Diários")) |> 
  mutate(
    gauge_code = factor(gauge_code, levels = order.gg),
    id = factor(id, levels = order.id),
    d = factor(d, groups)
  ) |> 
  # mutate(gauge_code = fct_reorder(gauge_code, statistic)) |> 
  ggplot(mapping = aes(y = statistic, ymin = ci_lower, ymax = ci_upper, color = d)) +
  # geom_errorbar(aes(x = id), width = 0, linewidth = 0.4, alpha = 0.5, position = position_dodge(width = dodge.width)) +
  ylim(min(df.scale$ci_lower), max(df.scale$ci_upper)) +
  geom_point(aes(x = id), shape = 16, size = 0.6, alpha = 0.9, position = position_dodge(dodge.width)) +
  labs(x = "ID", y = "H", color = "Faixa de duração") +
  scale_color_manual(
    name = NULL, values = colors, aesthetics = c("color", "fill")
  ) + 
  theme(
    legend.position = c(0, -0.02),
    legend.justification = c(0, 0),
    legend.key.width = unit(12, "pt"),
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
  )

ggsave(
  "figuras/preliminares/subdiarios/defesa/plot_scale_ci_todos.png",
  width = 15, height = 10, units = "cm", dpi = 300, bg = "white"
)

# Informações adicionais
# Valores regionais
df.scale |> 
  filter(
    gauge_code == "Regional"
  )

# Quais estações possuem duração base = 30 min
info.df |> 
  filter(time_step == 30)

df.scale |>
  group_by(d) |> 
  summarise(
    min = min(statistic),
    max = max(statistic),
    sd = sd(statistic)
  )

# MODELO REGIONAL HORÁRIO -------------------------------------------------

# Momentos e expoente de escala regionais
df.mom.regional <- 
  merge(x = scale.moments %>%
          mutate(gauge_code = "Regional") %>% 
          filter(d >= 1 & d <= 24),
        y = scale.coefficient %>%
          filter(gauge_code == "Regional", d == duration.intervals[3]) %>% 
          select(gauge_code, which_moment, scale, rsquared),
        by = c("gauge_code", "which_moment"),
        all.x = FALSE)

# Figura modelo regional horário
plot.scale.regional <- ggplot(data = df.mom.regional, aes(x = d, y = mom)) +
  geom_point(color = "grey", alpha = 0.5,
             aes(shape = factor(which_moment, labels = paste0("q = ", unique(which_moment),
                                                              ": H = ", round(unique(scale), 4),
                                                              " | R² = ", round(unique(rsquared), 3))))) +
  geom_smooth(aes(group = which_moment), color = "grey20", method = "lm", formula = y~x, se = FALSE, linewidth = 0.4) +
  scale_x_log10() +
  scale_y_log10(labels = scales::label_comma()) +
  labs(y = "E[I<sub>d</sub><sup>q</sup>]", x = "Durações [h]", shape = "",
       caption = paste0("Modelo regional | Durações: ", duration.intervals[3], " h")) +
  theme(legend.position =  c(0.98, 0.98),
        legend.justification = c(1,0.7),
        axis.title.y = ggtext::element_markdown()) +
  guides(shape = guide_legend(override.aes = list(alpha = 1, color = "grey70"))); plot.scale.regional

ggsave(plot = plot.scale.regional, filename = "figuras/scale_invariance/workshop/plot_scale_regional.png",
       width = 16, height = 12, units = "cm", dpi = 300)

# COEFICIENTES DE DESAGREGAÇÃO --------------------------------------------

source("scripts/script_coef_desagregacao.R")

ggsave(plot = plot.coef.disag, filename = "figuras/scale_invariance/workshop/plot_coef_desagregacao_block.png", bg = "white",
       width = 32, height = 16, units = "cm", dpi = 300)


# AVALIAÇÃO SENTIDO RESTRITO ---------------------------------------------

# Estimativas L-momentos
flmom <- function(x){
  x <- na.omit(x)
  lmom::pelgev(lmom::samlmu(x))
}

strict.sense <- imax.df |> 
  select(gauge_code, gauge_label, d, imax, wateryear) |> 
  group_by(gauge_code, gauge_label, d) |> 
  group_modify(
    ~as_tibble_row(flmom(.x$imax))
  )

strict.sense |> 
  pivot_longer(
    cols = -c(gauge_code, gauge_label, d),
    names_to = "par",
    values_to = "value"
  ) |> 
  filter(
    d >= 1 & d <= 24, gauge_code %in% c("A801", "A803", "A813"),
    par != "k"
  ) |> 
  mutate(
    par = case_when(
      par == "xi" ~ "Parâmetro de posição",
      par == "alpha" ~ "Parâmetro de escala"
    ),
    par = factor(par, levels = c("Parâmetro de posição", "Parâmetro de escala"))
  ) |> 
  ggplot(aes(x = d, y = value)) +
  facet_grid(par~gauge_label, scales = "free_y", switch = "y") +
  geom_smooth(
    aes(group = gauge_code), color = "grey20", linetype = "solid", alpha = 0.8,
    method = "lm", formula = y~x, se = FALSE, linewidth = 0.3
  ) +
  ggpmisc::stat_poly_eq(
    aes(label = after_stat(rr.label)), formula = y~x, rr.digits = 4, family = font, size = textsize, 
    label.x = Inf, hjust = 1.2, fontface = "plain"
  ) +
  geom_point(color = "black", fill = "black", size = 1) +
  scale_x_log10(name = "Durações [h]", breaks = c(1, 3, 10, 24)) +
  scale_y_log10(name = NULL) +
  theme(
    strip.placement = "outside",
    strip.text.y.left = element_text(hjust = 0.5, size = fontsize)
  )

ggsave(
  "figuras/preliminares/subdiarios/invariancia_escala/plot_strict_sense.png",
  width = 15, height = 10, units = "cm", dpi = 300
)
