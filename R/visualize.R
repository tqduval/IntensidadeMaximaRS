# PACOTES -----------------------------------------------------------------

# Limpar ambiente
rm(list = ls())
invisible(gc())

# Pacotes
if (!require(pacman)) {
  install.packages("pacman")
}
pacman::p_load(
  pacman, # gerenciador de pacotes
  arrow, # leitura de dados 'parquet'
  tidyverse
)

# Calcular quantis
source("R/funcoes/q_dgev.R")
get.q <- function(param, duration, p.nexceed) {
  param <- unname(param)
  xit <- param[1]
  alpha0 <- param[2]
  kappa <- param[3]
  H <- param[4]
  theta <- param[5]

  p <- p.nexceed
  d <- duration

  # Calcular parâmetros convencionais da GEV
  alpha <- alpha0 / (d + theta)^H
  xi <- xit * alpha

  yp <- xi + alpha / kappa * (1 - (-log(p))^kappa)

  return(yp)
}

# Set theme
suppressMessages(extrafont::loadfonts())
font <- "Consolas"
fontsize <- 8
textsize <- (fontsize - 2) / .pt
set_theme(
  theme_minimal() +
    theme(
      text = element_text(color = "black", size = fontsize, family = font),
      axis.text = element_text(color = "black"),
      # plot.margin = margin(),
      strip.text = element_text(
        hjust = 0,
        size = fontsize - 2,
        margin = margin(2.5, 2.5, 2.5, 2.5),
        color = "black"
      ),
      panel.border = element_rect(color = "black")
    )
)


# LER DADOS ---------------------------------------------------------------

# Intensidades máximas anuais de precipitação
df.imax <- arrow::read_parquet("base/gerados/df_imax.pqt")

# Ajuste
# ls.fit <- readRDS("base/gerados/resultados_com_priori/ls_fit.rds")
# df.par <- arrow::read_parquet("base/gerados/resultados_com_priori/df_par.pqt")
# df.idf <- arrow::read_parquet("base/gerados/resultados_com_priori/df_idf.pqt")

# Ler resultados com todas as simulações
df.par <- arrow::read_parquet("base/gerados/df_par.pqt")
df.idf <- arrow::read_parquet("base/gerados/df_idf.pqt")
df.shape.boot <- arrow::read_parquet("base/gerados/df_shapes_boot.pqt")
df.shape.delta <- arrow::read_parquet("base/gerados/df_shapes_delta.pqt")

combinacoes <- unique(df.par[, c("prior", "ci_type")])
# combinacao <- combinacoes[1,]
# combinacao <- combinacoes[2,]
# combinacao <- combinacoes[3,]
combinacao <- combinacoes[4, ]

# ORGANIZAÇÃO DOS DADOS ---------------------------------------------------

# durations <- unique(df.idf$d)
# durations <- c(1, 2, 6, 12, 24)
durations <- c(1, 6, 24)

data <- df.imax %>%
  filter(na_prct_yr <= 0.2, between(d, 1, 24)) %>%
  group_by(gauge_code) %>%
  filter(n_distinct(d) >= 3, n_distinct(wateryear) >= 8) %>%
  ungroup() %>%
  select(gauge_code, wateryear, d, imax)


## RAZÕES DE L-MOMENTOS ---------------------------------------------------

# Plotar gráficos de barras p/ razões de l-momentos
df.ratios <- data %>%
  group_by(gauge_code, d) %>%
  arrange(desc(imax)) %>%
  reframe(
    slmom = list(lmom::samlmu(imax)),
    lcv = unlist(slmom)[[2]] / unlist(slmom)[[1]],
    lsk = unlist(slmom)[[3]],
    lkt = unlist(slmom)[[4]]
  ) %>%
  select(-slmom) %>%
  pivot_longer(
    cols = starts_with("l"),
    names_to = "ratios",
    values_to = "value"
  ) %>%
  mutate(
    label = case_when(
      ratios == "lcv" ~ "L-CV",
      ratios == "lsk" ~ "L-assimetria",
      ratios == "lkt" ~ "L-curtose"
    ),
    label = factor(label, levels = c("L-CV", "L-assimetria", "L-curtose"))
  ) %>%
  arrange(gauge_code, d, ratios) %>%
  print()

# Boxplot
plot.ratios.d <- df.ratios %>%
  # filter(d %in% durations) %>% C
  ggplot(aes(x = factor(d), y = value)) +
  facet_wrap(~label, ncol = 1, scale = "free_y") +
  geom_boxplot(
    fill = "grey90",
    median.linewidth = 0.2,
    staplewidth = 0.5,
    whisker.linetype = "dashed",
    outlier.size = 0.5,
    outlier.alpha = 0.5,
    linewidth = 0.2
  ) +
  labs(x = "Duração [h]", y = "Razões de L-momentos")
plot.ratios.d

ggsave(
  "figuras/ratios_d.png",
  plot.ratios.d,
  width = 12,
  height = 14,
  units = "cm",
  dpi = 300
)


## TABELA DE ESTATÍSTICAS -------------------------------------------------

# Calcular estatística básicas p/ durações de referência (1, 6, 12, 24)

## CURVAS DE FREQUÊNCIA ---------------------------------------------------

### Várias durações -------------------------------------------------------

# Quantis empíricos (intensidades observadas)
df.fq.obs <- data %>%
  group_by(gauge_code, d) %>%
  arrange(desc(imax)) %>%
  mutate(
    n = length(imax),
    i = row_number(),
    p_weibull = i / (n + 1),
    tr_weibull = 1 / p_weibull,
    p_cunanne = (i - 0.4) / (n + 0.2),
    tr_cunanne = 1 / p_cunanne,
    p_gring = (i - 0.44) / (n + 0.12),
    tr_gring = 1 / p_gring,
    label = paste(gauge_code, "|", n, "anos")
  ) %>%
  ungroup() %>%
  arrange(gauge_code, d, i)

# Título dos facets
gauge_labels <- df.fq.obs %>%
  slice_head(n = 1, by = "gauge_code") %>%
  select(gauge_code, label, n) %>%
  arrange(desc(n)) %>%
  pull(label, name = gauge_code)

# Ordem dos facets
order <- df.fq.obs %>%
  slice_head(n = 1, by = "gauge_code") %>%
  select(gauge_code, label, n) %>%
  arrange(desc(n)) %>%
  pull(label, name = n)

# Quantis estimados
df.fq.est <- df.idf %>%
  filter(prior == combinacao[[1]], ci_type == combinacao[[2]]) %>%
  mutate(label = gauge_labels[gauge_code]) %>%
  left_join(
    df.par %>%
      filter(parname == "kappa") %>%
      rename(kappa = par) %>%
      select(gauge_code, kappa, prior, ci_type),
    by = c("gauge_code", "prior", "ci_type")
  )

priori.cap <- ifelse(
  isTRUE(combinacao[[1]]),
  "Estimativas com informação a priori",
  "Estimativas sem informação a priori"
)

ggplot() +
  facet_wrap(~ factor(label, levels = order), ncol = 4, scale = "free_y") +
  geom_line(
    aes(x = return_period, y = q, color = factor(d)),
    df.fq.est %>% filter(d %in% durations),
    linewidth = 0.2
  ) +
  geom_ribbon(
    aes(
      x = return_period,
      ymin = ci_lower,
      ymax = ci_upper,
      color = factor(d),
      fill = factor(d)
    ),
    df.fq.est %>% filter(d %in% durations),
    alpha = 0.1,
    linetype = "dashed",
    linewidth = 0.2
  ) +
  geom_point(
    aes(x = tr_gring, y = imax, color = factor(d), shape = factor(d)),
    df.fq.obs %>% filter(d %in% durations),
    size = 0.6
  ) +
  geom_text(
    aes(x = 1, y = Inf, label = sprintf("κ = %0.3f", kappa)),
    df.fq.est,
    stat = "unique",
    family = font,
    size = textsize,
    hjust = 0,
    vjust = 1.8
  ) +
  scale_x_log10() +
  labs(
    x = "T [anos]",
    y = "Intensidade [mm/h]",
    color = "Duração [h]",
    shape = "Duração [h]",
    fill = "Duração [h]"
  ) +
  # labs(x = "T [anos]", y = "Intensidade [mm/h]", color = "Duração [h]", shape = "Duração [h]", fill = "Duração [h]",
  #      caption = paste(priori.cap, "| Método IC: ", toupper(combinacao[[2]]))) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(),
    legend.key.size = unit(fontsize + 1, "pt"),
    strip.text = element_text(color = "black", size = fontsize - 2, hjust = 0)
  )
# theme(legend.position = c(0,1),
#       legend.justification = c(0,1),
#       legend.key.height = unit(fontsize, "pt"),
#       legend.title = element_text(size = fontsize - 2))

# Diretórios de saída
path.freq <- paste0(
  "figuras/resultados/curvas_frequencia/",
  combinacao[[2]],
  "/",
  ifelse(
    isTRUE(combinacao[[1]]),
    "curvas_freq_wpriori.png",
    "curvas_freq_npriori.png"
  )
) %>%
  print()

ggsave(path.freq, width = 15, height = 22, units = "cm", dpi = 300)

# ggsave("figuras/curvas_frequencia_com_priori.png", width = 20, height = 28, units = "cm", dpi = 300)
# ggsave(sprintf("figuras/curvas_frequencia_%s.png", gg), width = 12, height = 8, units = "cm", dpi = 300)

# Curvas de frequência para uma duração só comparando com e sem priori

# IDFs
ggplot() +
  facet_wrap(~ factor(label, levels = order), ncol = 4, scale = "free_y") +
  geom_ribbon(
    aes(
      x = d,
      ymin = ci_lower,
      ymax = ci_upper,
      fill = factor(return_period),
      color = factor(return_period)
    ),
    # df.fq.est %>% filter(return_period %in% c(100)),
    df.fq.est %>% filter(return_period %in% c(10, 50, 100, 500)),
    alpha = 0.2,
    linetype = "dashed",
    linetwidth = 0.1
  ) +
  geom_line(
    aes(x = d, y = q, color = factor(return_period)),
    df.fq.est %>% filter(return_period %in% c(10, 50, 100, 500)),
    # geom_line(aes(x = d, y = q, color = factor(return_period)), df.fq.est %>% filter(return_period %in% c(100)),
    linewidth = 0.3
  ) +
  scale_x_log10() +
  scale_y_log10() +
  labs(
    x = "Duration [h]",
    y = "Intensidade [mm/h]",
    color = "T [anos]",
    fill = "T [anos]",
    caption = paste(priori.cap, "| Método IC: ", toupper(combinacao[[2]]))
  ) +
  theme(
    legend.position = "bottom",
    legend.margin = margin(),
    legend.key.size = unit(fontsize + 1, "pt"),
    strip.text = element_text(color = "black", size = fontsize, hjust = 0)
  )


## RESÍDUOS ---------------------------------------------------------------

# Como o valor das estimativas de kappa em si não mudam com a combinação de método
# tanto faz apresentar com um específico os resultados da análise de resíduo
# Calcular quantis para os tempos de recorrência observados
df.par.comp <- df.par %>%
  filter(prior == combinacao[[1]], ci_type == combinacao[[2]]) %>%
  select(gauge_code, par) %>%
  group_by(gauge_code) %>%
  reframe(par = list(par))

# Calcular média e variância dos resíduos
plot.position <- "gring"
df.residuals <- df.fq.obs %>%
  pivot_longer(
    cols = starts_with(c("p_", "tr_")),
    names_to = c("what", "plot_pos"),
    names_sep = "_",
    values_to = "value"
  ) %>%
  pivot_wider(names_from = what) %>%
  filter(plot_pos == plot.position) %>%
  left_join(df.par.comp, by = "gauge_code") %>%
  rowwise() %>%
  mutate(
    q = get.q(unlist(par), d, 1 - p),
    eqm = (q - imax)^2,
    rel_var = (imax - q)^2 / imax
  ) %>%
  ungroup() %>%
  group_by(gauge_code, label, d, plot_pos) %>%
  reframe(
    eqm = sum(eqm) / n(), # erro quadrático médio (Isim - Iobs)^2
    rel_bias = (mean(q) - mean(imax))^2 / mean(q), # viés relativo
    rel_var = sum(rel_var) / (n() - 1) # nem ideia se iso ta certo ou não
  ) %>%
  print

### Histogramas de média e variância dos resíduos -------------------------

# Média dos resíduos - Viés (BIAS)
plot.rel.bias <- df.residuals %>%
  pivot_longer(
    # cols = starts_with("rel_"),
    cols = c(eqm, rel_bias, rel_var),
    names_to = "statistic"
  ) %>%
  # filter(statistic == "rel_bias") %>%
  filter(statistic == "eqm") %>%
  ggplot(aes(x = d, y = value, fill = statistic)) +
  facet_wrap(~ factor(label, levels = order), ncol = 4, scale = "free_y") +
  geom_col(
    width = 0.8,
    position = position_dodge(),
    color = "black",
    linewidth = 0.2,
    fill = "grey20"
  ) +
  scale_x_continuous(breaks = c(1, 6, 12, 18, 24)) +
  # labs(x = "Durações [h]", y = "Viés relativo médio [%]") +
  labs(x = "Durações [h]", y = "EQM") +
  theme(
    legend.position = c(1, 1.1),
    legend.justification = c(1, 1),
    legend.key.size = unit(3, "mm"),
    legend.key.spacing.y = unit(1, "mm"),

    # strip.text = element_text(hjust = 0, size = fontsize - 2, margin = margin(2, 2, 2, 2))
  )
plot.rel.bias

ggsave(
  # "figuras/resultados/metricas/residuos_rel_bias.png",
  "figuras/resultados/metricas/eqm.png",
  plot.rel.bias,
  width = 15,
  height = 22,
  units = "cm",
  dpi = 300
)

# Variância dos resíduos
plot.rel.bias <- df.residuals %>%
  pivot_longer(
    cols = starts_with("rel_"),
    names_to = "statistic"
  ) %>%
  # filter(gauge_code == "A801") %>%
  filter(statistic == "rel_var") %>%
  ggplot(aes(x = d, y = value, fill = statistic)) +
  facet_wrap(~ factor(label, levels = order), ncol = 4, scale = "free_y") +
  geom_col(
    width = 0.8,
    position = position_dodge(),
    color = "black",
    linewidth = 0.2,
    fill = "grey80"
  ) +
  scale_fill_manual(
    name = "",
    values = c("mean" = "grey20", "var" = "grey80"),
    labels = c("mean" = "Média", "var" = "Variância")
  ) +

  labs(x = "Durações [h]", y = "Variância relativa média dos resíduos") +
  theme(
    legend.position = c(1, 1.1),
    legend.justification = c(1, 1),
    legend.key.size = unit(3, "mm"),
    legend.key.spacing.y = unit(1, "mm")
  )
plot.rel.bias

ggsave(
  "figuras/resultados/metricas/residuos_rel_var.png",
  width = 15,
  height = 22,
  units = "cm",
  dpi = 300
)

# ggsave(paste0("figuras/media_variancia_residuosA801.png"),
#        width = 8, height = 6, units = "cm", dpi = 300)

labels <- c(
  "rel_bias" = "Viés relativo médio [%]",
  "rel_var" = "Variância relativa do resíduo"
)

# Boxplots por duração
df.residuals %>%
  pivot_longer(
    cols = starts_with("rel_"),
    names_to = "statistic"
  ) %>%
  ggplot(aes(x = d, y = value, group = d)) +
  facet_wrap(
    ~statistic,
    scales = "free_y",
    labeller = as_labeller(labels),
    ncol = 1,
    strip.position = "left"
  ) +
  geom_boxplot(
    staplewidth = 0.5,
    whisker.linetype = "dashed",
    fill = "grey90",
    outlier.size = 0.3,
    outlier.alpha = 0.6,
    median.linewidth = 0.3,
    linewidth = 0.3
  ) +
  scale_x_continuous(breaks = c(1, 6, 12, 18, 24)) +
  # coord_cartesian(clip = "off") +
  labs(x = "Durações [h]", y = "") +
  # theme(strip.text = element_text(hjust = 0, size = fontsize))
  theme(
    strip.placement = "outside",
    strip.text = element_text(size = fontsize, hjust = 0.5, vjust = 1.1),
    strip.text.y.left = element_text(margin = margin())
  )

ggsave(
  "figuras/resultados/metricas/residuos_boxplot.png",
  width = 12,
  height = 14,
  units = "cm",
  dpi = 300
)

## INCERTEZA EM KAPPA -----------------------------------------------------

# Visualizar intervalos de confiança para kappa com diferentes métodos
colors <- c("#2baa68", "#648fff", "#fe6100")
color.names <- setNames(colors, c("lmom", "mle", "rgmle"))
fill.names <- setNames(colors[2:3], c("mle", "rgmle"))

pos <- position_dodge(1.2)

df.shape.boot %>%
  filter(model == "gev") %>%
  bind_rows(
    df.shape.boot %>%
      filter(model == "dgev", method != "lmom") %>%
      slice_head(n = 1, by = c("gauge_code", "method")) %>%
      mutate(d = NA_real_)
  ) %>%
  mutate(
    label = gauge_labels[gauge_code],
    label = factor(label, levels = order),
  ) %>%
  # filter(ifelse(method == "mle", model != "gev", TRUE)) %>% # remover GEV-MLE
  arrange(desc(label)) %>%
  # filter(gauge_code == "A844") %>%
  ggplot(aes(x = factor(d), y = kappa, color = method)) +
  facet_wrap(~ factor(label, levels = order), ncol = 4) +
  geom_hline(
    yintercept = c(-0.5, 0.5),
    linetype = "dashed",
    linewidth = 0.2,
    color = "grey60"
  ) +
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = ci_lower,
      ymax = ci_upper,
      fill = method
    ),
    . %>%
      filter(
        model == "dgev",
        # ifelse(gauge_code == "431410002A", method != "mle", TRUE)
      ),
    color = "transparent",
    alpha = 0.3
  ) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    . %>% filter(model == "gev"),
    width = 0,
    position = pos,
    linewidth = 0.3,
    alpha = 0.7
  ) +
  scale_color_manual(
    name = "GEV",
    values = color.names,
    labels = c("lmom" = "LMOM", "mle" = "MLE", "rgmle" = "RGMLE")
  ) +
  scale_fill_manual(
    name = "d-GEV",
    values = fill.names,
    labels = c("mle" = "MLE", "rgmle" = "RGMLE")
  ) +
  geom_point(data = . %>% filter(model == "gev"), position = pos, size = 0.15) +
  labs(x = "Duração [h]", y = "κ") +
  scale_x_discrete(breaks = c(1, 6, 12, 18, 24)) +
  # scale_y_continuous(limits = c(-1.5, 1.5)) +
  coord_cartesian(ylim = c(-1.2, 1.2)) + # usar ylim dentro de coord_cartesian p/ evitar remoção de elementos
  theme(
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),
    legend.margin = margin(),
    panel.grid = element_line(linewidth = 0.2),
    panel.border = element_rect(
      linewidth = 0.3,
      color = "black",
      margin = margin(2, 2, 2, 2)
    ),
  )

ggsave(
  # "figuras/kappa_ci_boot.png",
  "figuras/resultados/kappa_ci_boot.png",
  width = 15,
  height = 22,
  units = "cm",
  dpi = 300
)

df.shape.delta %>%
  mutate(
    label = gauge_labels[gauge_code],
    label = factor(label, levels = order),
  ) %>%
  # filter(ifelse(method == "mle", model != "gev", TRUE)) %>% # remover GEV-MLE
  arrange(desc(label)) %>%
  # filter(gauge_code == "A844") %>%
  ggplot(aes(x = factor(d), y = kappa, color = method)) +
  facet_wrap(~ factor(label, levels = order), ncol = 4) +
  geom_hline(
    yintercept = c(-0.5, 0.5),
    linetype = "dashed",
    linewidth = 0.2,
    color = "grey60"
  ) +
  geom_rect(
    aes(
      xmin = -Inf,
      xmax = Inf,
      ymin = ci_lower,
      ymax = ci_upper,
      fill = method
    ),
    . %>%
      filter(
        model == "dgev",
        # ifelse(gauge_code == "431410002A", method != "mle", TRUE)
      ),
    color = "transparent",
    alpha = 0.3
  ) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper),
    . %>% filter(model == "gev"),
    width = 0,
    position = pos,
    linewidth = 0.3,
    alpha = 0.7
  ) +
  scale_color_manual(
    name = "GEV",
    values = color.names,
    labels = c("lmom" = "LMOM", "mle" = "MLE", "rgmle" = "RGMLE")
  ) +
  scale_fill_manual(
    name = "d-GEV",
    values = fill.names,
    labels = c("mle" = "MLE", "rgmle" = "RGMLE")
  ) +
  geom_point(data = . %>% filter(model == "gev"), position = pos, size = 0.15) +
  labs(x = "Duração [h]", y = "κ") +
  scale_x_discrete(breaks = c(1, 6, 12, 18, 24)) +
  # scale_y_continuous(limits = c(-1.5, 1.5)) +
  coord_cartesian(ylim = c(-1.2, 1.2)) + # usar ylim dentro de coord_cartesian p/ evitar remoção de elementos
  theme(
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),
    legend.margin = margin(),
    panel.grid = element_line(linewidth = 0.2),
    panel.border = element_rect(
      linewidth = 0.3,
      color = "black",
      margin = margin(2, 2, 2, 2)
    ),
  )

ggsave(
  "figuras/resultados/kappa_ci_delta.png",
  width = 15,
  height = 22,
  units = "cm",
  dpi = 300
)
