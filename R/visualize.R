# PACOTES -----------------------------------------------------------------

# Limpar ambiente
rm(list = ls()); invisible(gc())

# Pacotes
if (!require(pacman)) install.packages("pacman")
pacman::p_load(pacman, arrow, tidyverse, patchwork)

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
ggplot2::set_theme(
  theme_minimal() +
    theme(
      legend.position = "bottom",
      legend.box.margin = margin(-5, 2, 2, 2),
      legend.key.size = unit(fontsize, "pt"),
      text = element_text(color = "black", size = fontsize, family = font),
      axis.text = element_text(color = "black"),
      strip.text = element_text(hjust = 0, size = fontsize - 2, margin = margin(2, 2, 2, 2), color = "black"),
      panel.border = element_rect(color = "black"),
      panel.grid = element_blank()
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
df.shape.boot.lmom <- arrow::read_parquet("base/gerados/df_shapes_boot_lmom.pqt")
df.shape.delta <- arrow::read_parquet("base/gerados/df_shapes_delta.pqt")

df.par.byd <- arrow::read_parquet("base/gerados/resultados_por_duracao/df_par_by_d.pqt")
df.idf.byd <- arrow::read_parquet("base/gerados/resultados_por_duracao/df_idf_by_d.pqt")

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

# Facet labels
gauge.labels <- data %>% 
  mutate(
    n_years = n_distinct(wateryear), 
    label = paste(gauge_code, "|", n_years, "anos"),
    .by = c(gauge_code, d)
  ) %>% 
  slice_head(n = 1, by = "gauge_code") %>% 
  select(gauge_code, label, n_years) %>%
  arrange(desc(n_years)) %>%
  pull(label, name = gauge_code)

# Facet order
order <- data %>%
  mutate(
    n_years = n_distinct(wateryear), 
    label = gauge.labels[gauge_code],
    .by = c(gauge_code, d)
  ) %>% 
  slice_head(n = 1, by = "gauge_code") %>%
  select(gauge_code, label, n_years) %>%
  arrange(desc(n_years)) %>%
  pull(label, name = n_years)

# Order labels
gauge.labels <- factor(gauge.labels, levels = order)

# Combine parameter estimates by duration-wise and scale invariance models
df.par.all <- df.par %>% 
  mutate(
    model = "dgev", 
    method = case_when(
      prior == FALSE ~ "mle",
      prior == TRUE ~ "rgmle"
    ), 
    d = NA_real_
  ) %>% 
  select(-c(prior, which_prior)) %>% 
  bind_rows(
    df.par.byd %>% 
      mutate(
        model = "gev", 
        method = case_when(
          prior == "NULL" ~ "mle",
          prior == "c(-0.1, 0.122)" ~ "gmle",
          prior == "c(0.0178, 0.0627)" ~ "rgmle"
        )) %>% 
      select(-prior)
  ) %>% 
  mutate(label = gauge.labels[gauge_code]) %>% 
  arrange(label)

# Combine quantile estimates by duration-wise and scale invariance models

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

ggsave("figuras/ratios_d.png", plot.ratios.d, width = 12, height = 14, units = "cm", dpi = 300)


## TABELA DE ESTATÍSTICAS -------------------------------------------------

# Calcular estatística básicas p/ durações de referência (1, 6, 12, 24)


# VISUALIZAÇÃO -----------------------------------------------------------

## INCERTEZA EM KAPPA -----------------------------------------------------

# Visualizar intervalos de confiança para kappa com diferentes métodos
colors <- c("#2baa68", "#648fff", "#dc267f", "#fe6100")
color.names <- setNames(colors, c("lmom", "mle", "gmle", "rgmle"))
fill.names <- setNames(colors[c(2,4)], c("mle", "rgmle"))
model.labels <- c("gev_mle" = "GEV (MLE)", "gev_gmle" = "GEV (GMLE)", "gev_rgmle" = "GEV (RGMLE)", "dgev_mle" = "d-GEV (MLE)", "dgev_rgmle" = "d-GEV (RGMLE)")
# model.labels <- c("GEV (MLE)" = "gev_mle",  "GEV (GMLE)" = "gev_gmle", "GEV (RGMLE)" = "gev_rgmle", "d-GEV (MLE)" = "dgev_mle", "d-GEV (RGMLE)" = "dgev_rgmle")
model.colors <- c(setNames(colors[2:4], model.labels[1:3]), setNames(fill.names, model.labels[4:5]))

pos <- position_dodge(1.2)

# Plotar intervalos IDF invariância de escala com bootstrap
plot.kappa.ci.boot <- df.par.all %>% 
  mutate(
    model_label = paste0(model, "_", method),
    model_label = model.labels[model_label]
  ) %>%
  filter(parname == "kappa", ci_type == "bootstrap") %>% 
  rename(kappa = par) %>% 
  filter(ifelse(method == "mle", model != "gev", TRUE)) %>% # GEV-MLE estimates are very unstable
  ggplot(aes(x = factor(d), y = kappa, color = model_label)) +
  facet_wrap(~ factor(label, levels = order), ncol = 4) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, color = "grey60") +
  geom_rect(
    aes(xmin = -Inf, xmax = Inf, ymin = ci_lower, ymax = ci_upper, fill = model_label), . %>% filter(model == "dgev", method == "rgmle"),
    color = "transparent", alpha = 0.3
  ) +
  geom_hline(aes(yintercept = kappa, color = model_label), . %>% filter(model == "dgev", method == "rgmle"), linewidth = 0.3) +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper), . %>% filter(model == "gev"),
    width = 0, position = pos, linewidth = 0.3, alpha = 0.7
  ) +
  geom_point(data = . %>% filter(model == "gev"), position = pos, size = 0.15) +
  scale_x_discrete(breaks = c(1, 6, 12, 18, 24)) +
  scale_color_manual(name = "", values = model.colors, breaks = names(model.colors), aesthetics = c("color", "fill")) +
  labs(x = "Durações [h]", y = "κ"); plot.kappa.ci.boot

ggsave("figuras/resultados/kappa_ci_boot.png", plot.kappa.ci.boot, width = 15, height = 22, units = "cm", dpi = 300)

# Plotar apenas estimativas com L-momentos
# Resultados com 5000 réplicas bootstrap
plot.shape.boot.lmom <- df.shape.boot.lmom %>%
  mutate(
    label = gauge_labels[gauge_code],
    label = factor(label, levels = order)
  ) %>%
  arrange(desc(label)) %>%
  ggplot(aes(x = factor(d), y = kappa, color = method)) +
  facet_wrap(~ factor(label, levels = order), ncol = 4) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, color = "grey60") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), . %>% filter(model == "gev"), width = 0, linewidth = 0.3, alpha = 0.7) +
  geom_point(data = . %>% filter(model == "gev"), size = 0.15) +
  scale_color_manual(name = "GEV", values = color.names, labels = c("lmom" = "LMOM")) +
  labs(x = "Duração [h]", y = "κ", caption = ) +
  scale_x_discrete(breaks = c(1, 6, 12, 18, 24)) +
  coord_cartesian(ylim = c(-1.5, 1.5)) + # usar ylim dentro de coord_cartesian p/ evitar remoção de elementos
  theme(
    legend.position = "bottom", legend.location = "panel",
    legend.key.size = unit(3, "mm"),
    panel.grid = element_line(linewidth = 0.2),
    panel.border = element_rect(linewidth = 0.3, color = "black"),
    strip.text = element_text(hjust = 0)
  ); plot.shape.boot.lmom

ggsave("figuras/resultados/kappa_ci_boot_lmom.png", width = 15, height = 22, units = "cm", dpi = 300)

# Plotar estações selecionadas p/ corpo da dissertação
select.gg <- c("A801", "432300204A", "430460601A")

 plot.shape.boot.lmom <- df.shape.boot.lmom %>% 
  filter(gauge_code %in% select.gg) %>% 
  mutate(
    label = gauge_labels[gauge_code],
    label = factor(label, levels = order)
  ) %>%
  arrange(desc(label)) %>%
  ggplot(aes(x = factor(d), y = kappa, color = method)) +
  facet_wrap(~ factor(label, levels = order), ncol = 4) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, color = "grey60") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, linewidth = 0.3, alpha = 0.7) +
  geom_point(size = 0.15) +
  scale_color_manual(name = "GEV", values = color.names, labels = c("lmom" = "LMOM")) +
  labs(x = "Duração [h]", y = "κ", caption = ) +
  scale_x_discrete(breaks = c(1, 6, 12, 18, 24)) +
  # coord_cartesian(ylim = c(-1.5, 1.5)) + # usar ylim dentro de coord_cartesian p/ evitar remoção de elementos
  theme(
    # legend.position = "bottom", legend.location = "panel",
    legend.position = "none",
    legend.key.size = unit(3, "mm"),
    panel.grid = element_line(linewidth = 0.2),
    panel.border = element_rect(linewidth = 0.3, color = "black"),
    strip.text = element_text(hjust = 0)
  ); plot.shape.boot.lmom

ggsave("figuras/resultados/kappa_3_gg_horizontal.png", plot.shape.boot.lmom, width = 15, height = 8, units = "cm", dpi = 300)

# df.shape.delta %>%
#   mutate(
#     label = gauge_labels[gauge_code],
#     label = factor(label, levels = order),
#   ) %>%
#   # filter(ifelse(method == "mle", model != "gev", TRUE)) %>% # remover GEV-MLE
#   arrange(desc(label)) %>%
#   # filter(gauge_code == "A844") %>%
#   ggplot(aes(x = factor(d), y = kappa, color = method)) +
#   facet_wrap(~ factor(label, levels = order), ncol = 4) +
#   geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, color = "grey60") +
#   geom_rect(
#     aes(xmin = -Inf, xmax = Inf, ymin = ci_lower, ymax = ci_upper, fill = method), . %>% filter(model == "dgev"),
#     color = "transparent", alpha = 0.3
#   ) +
#   geom_errorbar(
#     aes(ymin = ci_lower, ymax = ci_upper), . %>% filter(model == "gev"),
#     width = 0, position = pos, linewidth = 0.3, alpha = 0.7
#   ) +
#   geom_point(data = . %>% filter(model == "gev"), position = pos, size = 0.15) +
#   scale_color_manual(name = "GEV", values = color.names, labels = c("lmom" = "LMOM", "mle" = "MLE", "rgmle" = "RGMLE")) +
#   scale_fill_manual(name = "d-GEV", values = fill.names, labels = c("mle" = "MLE", "rgmle" = "RGMLE")) +
#   labs(x = "Duração [h]", y = "κ") +
#   scale_x_discrete(breaks = c(1, 6, 12, 18, 24)) +
#   # scale_y_continuous(limits = c(-1.5, 1.5)) +
#   # coord_cartesian(ylim = c(-1.2, 1.2)) + # usar ylim dentro de coord_cartesian p/ evitar remoção de elementos
#   theme(
#     legend.position = "bottom",
#     legend.key.size = unit(3, "mm"),
#     legend.margin = margin(),
#     panel.grid = element_line(linewidth = 0.2),
#     panel.border = element_rect(linewidth = 0.3, color = "black", margin = margin(2, 2, 2, 2)),
#   )

# ggsave("figuras/resultados/kappa_ci_delta.png", width = 15, height = 22, units = "cm", dpi = 300)

# Resultados dos ajustes por estação
priors <- unique(df.par.byd$prior)
prior.colors <- setNames(color.names[c("mle", "gmle", "rgmle")], priors)

# Juntar todas as estimativas locais em um único data.frame
df.par.byd.all <- df.par.byd %>% 
  filter(parname == "kappa") %>% 
  rename(kappa = par, method = prior) %>% 
  select(gauge_code, method, d, kappa, ci_lower,  ci_upper) %>% 
  bind_rows(
    df.shape.boot.lmom %>% 
      select(gauge_code, method, d, kappa, ci_lower, ci_upper)
  ) %>% 
  mutate(
    label = gauge_labels[gauge_code],
    label = factor(label, levels = order),
    method = case_when(
      method == "NULL" ~ "mle",
      method == "c(-0.1, 0.122)" ~ "gmle",
      method == "c(0.0178, 0.0627)" ~ "rgmle",
      method == "lmom" ~ "lmom"
    )
  ) %>% 
  arrange(desc(label))

# Plotar todas as estações (Apêndice)
df.par.byd.all %>% 
  filter(method != "mle") %>% 
  ggplot(aes(x = factor(d), y = kappa, color = method)) +
  facet_wrap(~factor(label, levels = order), ncol = 4) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, position = pos, linewidth = 0.3, alpha = 0.7) +
  geom_point(position = pos, size = 0.15) +
  # scale_color_manual(name = "GEV", values = prior.colors) +
  scale_color_manual(name = "GEV", values = color.names, labels = c("lmom" = "LMOM", "gmle" = "GMLE", "rgmle" = "RGMLE")) +
  labs(x = "Duração [h]", y = "κ") +
  scale_x_discrete(breaks = c(1, 6, 12, 18, 24)) +
  coord_cartesian(ylim = c(-0.3, 0.3)) + # usar ylim dentro de coord_cartesian p/ evitar remoção de elementos
  # coord_cartesian(ylim = c(-0.5, 0.5)) + # usar ylim dentro de coord_cartesian p/ evitar remoção de elementos
  # coord_cartesian(ylim = c(-1.2, 1.2)) + # usar ylim dentro de coord_cartesian p/ evitar remoção de elementos
  theme(
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),
    panel.grid = element_line(linewidth = 0.2),
    panel.border = element_rect(linewidth = 0.3, color = "black"),
  )

ggsave("figuras/resultados/kappa_ci_boot_por_duracao.png", width = 15, height = 22, units = "cm", dpi = 300)

# Plotar apenas três estações
plot.par.byd <- df.par.byd.all %>% 
  filter(
    method != "mle",
    gauge_code %in% select.gg
  ) %>% 
  ggplot(aes(x = factor(d), y = kappa, color = method)) +
  facet_wrap(~factor(label, levels = order), ncol = 4) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, position = pos, linewidth = 0.3, alpha = 0.7) +
  geom_point(position = pos, size = 0.15) +
  # scale_color_manual(name = "GEV", values = prior.colors) +
  scale_color_manual(name = "GEV", values = color.names, labels = c("lmom" = "LMOM", "gmle" = "GMLE", "rgmle" = "RGMLE")) +
  labs(x = "Duração [h]", y = "κ") +
  scale_x_discrete(breaks = c(1, 6, 12, 18, 24)) +
  coord_cartesian(ylim = c(-0.3, 0.3)) + # usar ylim dentro de coord_cartesian p/ evitar remoção de elementos
  # coord_cartesian(ylim = c(-0.5, 0.5)) + # usar ylim dentro de coord_cartesian p/ evitar remoção de elementos
  # coord_cartesian(ylim = c(-1.2, 1.2)) + # usar ylim dentro de coord_cartesian p/ evitar remoção de elementos
  theme(
    legend.position = "bottom",
    legend.key.size = unit(3, "mm"),
    panel.grid = element_line(linewidth = 0.2),
    panel.border = element_rect(linewidth = 0.3, color = "black"),
  ); plot.par.byd

plot.shape.3.gauges <- 
  (plot.shape.boot.lmom + labs(subtitle = "(a)")) / 
  (plot.par.byd + labs(subtitle = "(b)")) + 
  theme(strip.text = element_blank()) +
  plot_layout(axis_titles = "collect"); plot.shape.3.gauges

ggsave("figuras/resultados/kappa_3_gg.png", plot.shape.3.gauges, width = 15, height = 12, units = "cm", dpi = 300)

# Plotar três estações com faixas da invariância de escala


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

# FAZER: pode ser interessante comparar os resultados com e sem priori regional
# usando o EQM. Selecionar durações de referência 1, 6, 12, 24 h e calcular o EQM
# Fazer um plot p/ cada duração e plotar o EQM (eixo-y) de toda as estações (eixo-x)
# ordenadas por tamanho das séries.

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
    residuo = q - imax, # calcular resíduo p/ ilustrar positivos e negativos
    eqm = residuo^2,
  ) %>%
  ungroup() %>%
  group_by(gauge_code, label, d, plot_pos) %>%
  reframe(
    eqm = sum(eqm) / n(), # erro quadrático médio (Isim - Iobs)^2, não pode ser comparado entre durações
    rel_bias = (mean(q) - mean(imax))^2 / mean(q), # viés relativo
    rel_res = mean(residuo) / (mean(q)) * 100, # resíduo relativo
    rel_var = sum(residuo^2) / (n() - 1) # nem ideia se iso ta certo ou não
  ) %>%
  print

### Histogramas de média e variância dos resíduos -------------------------

# Média dos resíduos relativo
plot.rel.bias <- df.residuals %>%
  pivot_longer(
    # cols = starts_with("rel_"),
    cols = c(eqm, rel_bias, rel_res, rel_var),
    names_to = "statistic"
  ) %>%
  filter(statistic == "rel_res") %>%
  # filter(statistic == "eqm") %>%
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
  labs(x = "Durações [h]", y = "Resíduo relativo [%]") +
  theme(
    legend.position = c(1, 1.1),
    legend.justification = c(1, 1),
    legend.key.size = unit(3, "mm"),
    legend.key.spacing.y = unit(1, "mm"),

    # strip.text = element_text(hjust = 0, size = fontsize - 2, margin = margin(2, 2, 2, 2))
  )
plot.rel.bias

ggsave(
  paste0(
    "figuras/resultados/metricas/residuo_relativo_",
    ifelse(isTRUE(combinacao[[1]]), "_wpriori.png", "_npriori.png")
  ),
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

labels <- c("rel_bias" = "Viés relativo médio [%]", "rel_var" = "Variância relativa do resíduo")

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

ggsave("figuras/resultados/metricas/residuos_boxplot.png", width = 12, height = 14, units = "cm", dpi = 300)




## INCERTEZA NOS QUANTIS -------------------------------------------------

# Calcular  
df.idf
