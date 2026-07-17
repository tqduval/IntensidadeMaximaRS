
rm(list = ls()); invisible(gc())
pacman::p_load(pacman, tidyverse, arrow, ggtext)

# Set theme
suppressMessages(extrafont::loadfonts())
font <- "JetBrainsMono NF"
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
      panel.border = element_rect(color = "black"),
      panel.grid = element_blank(),
      plot.margin = margin()
    )
)


# LER DADOS --------------------------------------------------------------

na.yr <- 0.2
min.durations <- 3
d.interval <- c(1, 24)
min.yr <- 8

# Dados observados
imax.df <- arrow::read_parquet("base/gerados/df_imax.pqt")
imax.df <- imax.df |> 
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
  select(gauge_code, wateryear, d, imax, n_years, gauge_label)


## d-GEV -----------------------------------------------------------------

# Ler lista resultados perfil verossimilhança
# fitprof_dgev.R
fitprof.dgev.nprior.ls <- readRDS("base/gerados/ls_fit_dgev_prof_nprior.rds") # versão sem priori
fitprof.dgev.wprior.ls <- readRDS("base/gerados/ls_fit_dgev_prof_wprior.rds") # versão com priori

# Checar resultados da rodada 'wprior' com correções
# fitprof.wprior.ls <- readRDS("base/gerados/ls_fit_dgev_prof_wprior2.rds") # temp
# par.df <- fitprof.wprior.ls$par # temp
# idf.df <- fitprof.wprior.ls$idf # temp
# par.df <- fitprof.nprior.ls$par # temp
# idf.df <- fitprof.nprior.ls$idf # temp

# Parâmetros
par.dgev.df <- bind_rows(
  fitprof.dgev.nprior.ls$par,
  fitprof.dgev.wprior.ls$par
)

# Quantis
idf.dgev.df <- bind_rows(
  fitprof.dgev.nprior.ls$idf,
  fitprof.dgev.wprior.ls$idf
)

# Atribuir nomes aos dados em `prior`
par.dgev.df <- par.dgev.df |> mutate(model = ifelse(prior == "NULL", "MLE", "RGMLE"))
idf.dgev.df <- idf.dgev.df |> mutate(model = ifelse(prior == "NULL", "MLE", "RGMLE"))

durations <- c(1, 6, 12, 24)

# Profiles
prof.par <- bind_rows(
  fitprof.nprior.ls$prof.par, # 1
  fitprof.wprior.ls$prof.par, # 2
  .id = "model"
) |> 
  mutate(
    model = ifelse(model == 1, "MLE", "RGMLE")
  ) |> 
  relocate(model, .after = last_col())

prof.q <- bind_rows(
  fitprof.nprior.ls$prof.q, # 1
  fitprof.wprior.ls$prof.q, # 2
  .id = "model"
) |> 
  mutate(
    model = ifelse(model == 1, "MLE", "RGMLE")
  ) |> 
  relocate(model, .after = last_col())


## GEV POR DURAÇÃO -------------------------------------------------------

# Estimativas do parâmetro de forma com bootstrap e L-momentos
# shape_ci_boot_only_lmom.R
shape.boot.df <- arrow::read_parquet("base/gerados/df_shapes_boot_lmom.pqt")

# Estimativas com perfil de verossimilhança por duração (tem que rodar ainda)
# Por duração sem priori
fitprof.gev.nprior.ls <- readRDS("base/gerados/ls_fit_gev_prof_nprior.rds")
fitprof.gev.wprior.ls <- readRDS("base/gerados/ls_fit_gev_prof_wprior.rds")
fitprof.gev.wprior.ninfo.ls <- readRDS("base/gerados/ls_fit_gev_prof_wprior_ninfo.rds")

# Parâmetros
par.gev.df <- bind_rows(
  fitprof.gev.nprior.ls$par,
  fitprof.gev.wprior.ls$par,
  fitprof.gev.wprior.ninfo.ls$par
)

which.prior <- unique(par.gev.df$prior)
prior.names <- setNames(c("MLE", "RGMLE", "GMLE"), which.prior)

# Quantis
idf.gev.df <- bind_rows(
  fitprof.gev.nprior.ls$idf,
  fitprof.gev.wprior.ls$idf,
  fitprof.gev.wprior.ninfo.ls$idf
)

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


# Cores
which.model <- unique(par.df$model)
col <- c("#648FFF", "#FE6100")
colors <- c("#2baa68", "#648fff", "#dc267f", "#fe6100")
model.col <- setNames(col, which.model)


# Durações de referência
durations <- c(1, 6, 12, 24)
d.shapes <- setNames(c(15, 16, 17, 18), durations)
d.colors <- setNames(c("red", "orange", "green", "steelblue"), durations)


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
    label = factor(label, levels = c("L-CV", "L-assimetria", "L-curtose"))
  ) |> 
  arrange(gauge_code, d, ratios)


# Identificar outliers
is.outlier <- function(x) x < quantile(x, 0.25) - 1.5*IQR(x) | x > quantile(x, 0.75) + 1.5*IQR(x)
ratios.df |> 
  group_by(d, ratios) |> 
  mutate(outlier = is.outlier(value)) |> 
  filter(outlier)


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
    fill = "grey90", median.linewidth = 0.2, staplewidth = 0.5, linewidth = 0.2,
    whisker.linetype = "dashed", outlier.size = 0.5, outlier.alpha = 0.5
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


## DIAGRAMA DE L-MOMENTOS ------------------------------------------------

# Plotar o diagrama, mas somente as cores das "durações de referência" em `durations`
durations

## ESTATÍSTICAS ----------------------------------------------------------

# Calcular estatísticas básicas p/ durações de referência (1, 6, 12, 24)
  

## INCERTEZA EM KAPPA ----------------------------------------------------

select.gg <- c("A801", "A831", "432300204A")
pos.width <- position_dodge(0.6)


# Visualizar intervalos de confiança de kappa POR DURAÇÃO com estimativas
# da GEV com L-momentos --- figura c/ todas as estações
shape.boot.plot <- shape.boot.df |> 
  mutate(gauge_label = gauge.labels[gauge_code]) |> 
  select(!c(model, method, id)) |> 
  ggplot(aes(x = factor(d), y = kappa)) +
  facet_wrap(~factor(gauge_label, levels = order.n), ncol = 4, scales = "free_y") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, color = "grey60") +
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
    strip.text = element_text(hjust = 0)
  ); shape.boot.plot

ggsave("figuras/preliminares/ci_kappa_gev_lmom_boot.png", shape.boot.plot, height = 22, width = 16, units = "cm", dpi = 300)


# Visualizar intervalos de confiança de kappa POR DURAÇÃO com estimativas
# da GEV com L-momentos --- figura c/ estações em `select.gg`
shape.boot.plot <- shape.boot.df |> 
  filter(gauge_code %in% select.gg) |> 
  mutate(gauge_label = gauge.labels[gauge_code]) |> 
  select(!c(model, method, id)) |> 
  ggplot(aes(x = factor(d), y = kappa)) +
  facet_wrap(~factor(gauge_label, levels = order.n), ncol = 4) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, color = "grey60") +
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
    strip.text = element_text(hjust = 0)
  ); shape.boot.plot

ggsave(
  "figuras/preliminares/ci_kappa_gev_lmom_boot_selected.png", shape.boot.plot, 
  height = 5, width = 14, units = "cm", dpi = 300
)


# Visualizar intervalos de confiança para kappa POR DURAÇÃO
# com estimativas com perfil de verossimilhança + priori Beta(0, 0.122^2)
shape.prof.plot <- par.gev.df |> 
  mutate(
    gauge_label = gauge.labels[gauge_code],
    prior_label = prior.names[prior]
  ) |> 
  filter(
    parname == "kappa",
    prior_label == "GMLE"
  ) |> 
  ggplot(aes(x = factor(d), y = par)) +
  facet_wrap(~factor(gauge_label, order.n), ncol = 4, scales = "free_y") +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, color = "grey60") +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper), 
    width = 0, linewidth = 0.3, alpha = 0.7, color = colors[3], position = pos.width) +
  geom_point(size = 0.15, position = pos.width, color = colors[3]) +
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
    strip.text = element_text(hjust = 0)
  ); shape.prof.plot
  
ggsave(
  "figuras/preliminares/ci_kappa_gev_gmle_prof.png", shape.prof.plot, 
  height = 22, width = 16, units = "cm", dpi = 300
 )
  



# Visualizar intervalos de confiança p/ kappa com POR DURAÇÃO
# com estimativas com perfil de verossimilhança + priori Beta(0, 0.122^2)
# p/ estações selecionadas --- `select.gg`

shape.prof.plot <- par.gev.df |> 
  filter(gauge_code %in% select.gg) |> 
  mutate(
    gauge_label = gauge.labels[gauge_code],
    prior_label = prior.names[prior]
  ) |> 
  filter(
    parname == "kappa",
    prior_label == "GMLE"
  ) |> 
  ggplot(aes(x = factor(d), y = par)) +
  facet_wrap(~factor(gauge_label, order.n), ncol = 4) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, color = "grey60") +
  geom_errorbar(
    aes(ymin = ci_lower, ymax = ci_upper), 
    width = 0, linewidth = 0.3, alpha = 0.7, color = colors[3], position = pos.width) +
  geom_point(size = 0.15, position = pos.width, color = colors[3]) +
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
    strip.text = element_text(hjust = 0)
  ); shape.prof.plot

ggsave(
  "figuras/preliminares/ci_kappa_gev_gmle_prof_selected.png", shape.prof.plot, 
  height = 5, width = 14, units = "cm", dpi = 300
)



# Visualizar intervalos de confiança para kappa com estimativas
# do modelo de invariância de escala com perfil de verossimilhança
 par.dgev.df |> 
  filter(parname == "kappa") |> 
  left_join(
    imax.df |> 
      slice_head(n = 1, by = "gauge_code") |> 
      select(gauge_code, n_years, gauge_label),
    by = "gauge_code"
  ) |> 
  # ggplot(aes(x = gauge_label, y = par)) +
  ggplot(aes(x = fct_reorder(gauge_code, n_years), y = par, color = model)) +
  geom_hline(yintercept = c(-0.5, 0.5), linetype = "dashed", linewidth = 0.2, color = "grey60") +
  geom_hline(yintercept = 0, linetype = "solid", linewidth = 0.2, color = "grey60") +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), width = 0, linewidth = 0.3, alpha = 0.7, position = pos.width) +
  geom_point(size = 0.15, position = pos.width) +
  scale_color_manual(name = "", values = model.col) +
  # labs(y = "H") +
  labs(y = "κ") +
  theme(
    legend.position = c(1, 0),
    legend.justification = c(1, 0),
    legend.margin = margin(5, 5, 1, 5),
    legend.key.size = unit(3, "mm"),
    axis.title.x = element_blank(),
    axis.text.x = element_text(angle = 90, hjust = 1),
    panel.border = element_rect(linewidth = 0.3, color = "black"),
    strip.text = element_text(hjust = 0)
  )

# ggsave("figuras/preliminares/ci_H_dgev_proflik2.png", height = 8, width = 14, units = "cm", dpi = 300)
ggsave("figuras/preliminares/ci_kappa_dgev_proflik.png", height = 8, width = 14, units = "cm", dpi = 300)

# Estações p/ verificar
par.dgev.df |> 
  filter(parname == "kappa") |> 
  select(gauge_code, par, ci_lower, ci_upper, prior) |> 
  filter(abs(par) < 1e-5 | abs(par) > 0.4) |> 
  unique() |> 
  print()


## PLOTAR ALGUNS PERFIS DE VEROSSIMILHANÇA -------------------------------

# `parname == kappa`
gg <- "A801"
# gg <- "430060402A"
proflik.kappa <- fitprof.wprior.ls$prof.par |> 
# proflik.kappa <- fitprof.nprior.ls$prof.par |> 
  filter(
    parname == "kappa", 
    gauge_code == gg
  ) |> 
  select(par, loglik) |> 
  rename(kappa = par) |> 
  filter(loglik < -180)

  
ggplot(proflik.kappa, aes(x = kappa, y = loglik)) +
  geom_line() +
  geom_point(shape = 21, fill = "white") +
  labs(subtitle = gg)

par.df |> filter(gauge_code == gg, parname == "kappa", prior != "NULL")

# quantil
# proflik.q <- fitprof.wprior.ls$prof.q |> 
proflik.q <- fitprof.dgev.nprior.ls$prof.q |> 
  filter(
    gauge_code == gg,
    d == durations[1], 
    return_period == 10
  ) |> 
  select(q, loglik)

ggplot(proflik.q, aes(x = q, y = loglik)) +
  geom_point() +
  labs(
    subtitle = gg,
    caption = sprintf("d = %g h | rp = %g anos", durations[1], 10)
  )

# IDF --------------------------------------------------------------------

## CURVAS DE FREQUÊNCIA --------------------------------------------------

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
  )
    
freq.est.df <- idf.dgev.df |> 
  filter(
    d %in% durations,
    return_period <= 100
  ) |> 
  mutate(
    n_years = imax.df$n_years[match(gauge_code, imax.df$gauge_code)],
    gauge_label = sprintf("%s | %2g anos", gauge_code, n_years)
  ) |> 
  left_join(
    par.df |> 
      filter(parname == "kappa") |> 
      rename(kappa = par) |> 
      select(gauge_code, kappa),
    by = c("gauge_code")
  )

prior.caption <- c("sem informação a priori", "priori regional informativa")
file.name <- c("nprior", "wprior")

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

