
rm(list = ls()); invisible(gc())
set.seed(1234)


# PACOTES -----------------------------------------------------------------

pacman::p_load(pacman, tidyverse, arrow, boot)


# LER DADOS ---------------------------------------------------------------

df.info <- arrow::read_parquet(file = "base/fonte/subdaily/df_subdaily_info.pqt")
val.gauge <- "A801"
imax.df <- arrow::read_parquet(file = "base/gerados/df_imax.pqt")
sim.df <- imax.df |> 
  mutate(responsible = df.info$responsible[match(gauge_code, df.info$gauge_code)]) %>% 
  filter(
    responsible != "SGB",
    # gauge_code != val.gauge # não filtrar, p/ gerar dados p/ dissertação
  )
val.df <- imax.df |> 
  filter(gauge_code == val.gauge)


# FUNÇÃO ------------------------------------------------------------------

source("R/funcoes/fun_scale_block_boot.R")


# RODAR FUNÇÃO ------------------------------------------------------------

# Argumentos
na.accept <- 0.2
min.years <- 8
which.moment <- 1:3
min.duration <- 3

# Intervalos de duração
duration.intervals <- list(c(10/60, 10*24), # todas as durações
                           c(10, 60)/60,    # sub-horários
                           c(1, 24),        # horários
                           c(1, 10)*24)     # diários

block.length <- 1        # tamanho dos blocos: 1 ano
R.boot <- 1e4            # nro. réplicas bootstrap
ci.boot.method <- "perc" # intervalo percentil
return.boot <- TRUE      # retornar objeto boot como atributo (attr)

# Lista de durações alvo p/ desagregação
disag.list <- list(NULL,                        # todas as durações
                   c(5, 10, 15, 20, 25, 30)/60, # sub-horários
                   c(1, 6, 8, 10, 12),          # horários
                   NULL)                        # diários

# Rodar bootstrap por bloco
ls.scale.boot <- lapply(X = seq_along(duration.intervals), FUN = function(i){
  
  simple.scaling.boot(df.imax = sim.df,
                      na.accept = na.accept,
                      min.years = min.years,
                      offset = FALSE, # não estimar parâmetro de offset das durações (theta)
                      which.duration = duration.intervals[[i]],
                      which.moment = which.moment,
                      d.target = disag.list[[i]],
                      ci.boot.method = ci.boot.method,
                      R.boot = R.boot,
                      return.boot = TRUE)
  
})


df.scale.boot <- bind_rows(ls.scale.boot)
dir.create("base/gerados/scale_invariance/")
arrow::write_parquet(x = df.scale.boot, sink = "base/gerados/scale_invariance/df_scale_block_boot.pqt")


# VISUALIZAÇÃO INTERVALOS DE CONFIANÇA ------------------------------------

# Ler e organizar dados de expoente de escala e intervalos de confiança
df.scale.boot <- arrow::read_parquet("base/gerados/scale_invariance/df_scale_block_boot.pqt")

groups <- unique(df.scale.boot$d)
group.names <- c("Todas", "Sub-horários", "Horários", "Diários")
groups <- setNames(object = group.names, nm = groups)

df.scale.boot <- df.scale.boot %>% 
  mutate(d = recode(d, !!!groups),
         d = factor(d, levels = group.names))

# Filtrar expoentes de escala 'scale'
df.scale <- df.scale.boot %>% 
  filter(which_variable == "scale",
         which_moment == 1,
         d != "Todas") %>% 
  select(gauge_code, n_year, d, statistic, ci_lower, ci_upper)

colors <- c("Todas" = "black", "Sub-horários" = "orange", "Horários" = "steelblue2", "Diários" =  "green3")

# Ordem baseado nos scale's horários
order <- df.scale.boot %>% 
  filter(d == "Horários") %>% 
  arrange(statistic) %>% 
  pull(gauge_code) %>% 
  unique()

# Dados locais
df.local <- df.scale %>% 
  mutate(gauge_code = factor(gauge_code, levels = order)) %>% 
  filter(gauge_code != "Regional")

# Dados regionais
df.regional <- df.scale %>% 
  filter(gauge_code == "Regional") %>% 
  select(-gauge_code)

dodge.width <- 1

# Comparações locais
ggplot(mapping = aes(y = statistic, ymin = ci_lower, ymax = ci_upper, color = d)) +
  geom_errorbar(data = df.local, aes(x = gauge_code), width = 0, linewidth = 0.8, alpha = 0.5, position = position_dodge(width = dodge.width)) +
  geom_point(data = df.local, aes(x = gauge_code), shape = 16, size = 1.5, alpha = 0.9, position = position_dodge(dodge.width)) +
  labs(x = "", y = "H", color = "Faixa de duração") +
  scale_color_manual(values = colors, aesthetics = c("color", "fill")) + 
  theme_minimal() +
  theme(legend.position = c(0.02,0.02),
        legend.justification = c(0,0),
        legend.key.spacing.y = unit(-2, "mm"),
        legend.margin = margin(l = 2, r = 2, unit = "mm"),
        legend.background = element_rect(color = "black", linewidth = 0.25),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "serif", color = "black", size = 12))

ggsave(filename = "figuras/scale_invariance/scale_confidence_interval_block_bootstrap_c2024.png",
       width = 16, height = 12, units = "cm", dpi = 300, bg = "white")

# Comparação c/ intervalos regionais
ggplot(mapping = aes(y = statistic, ymin = ci_lower, ymax = ci_upper, color = d)) +
  # facet_wrap(~d, ncol = 1, scales = "free_y") +
  geom_rect(data = df.regional, color = "transparent", aes(xmin = -Inf, xmax = Inf, fill = d), alpha = 0.2, show.legend = FALSE) +
  geom_errorbar(data = df.local, aes(x = gauge_code), width = 0, linewidth = 0.8, alpha = 0.5, position = position_dodge(width = dodge.width)) +
  geom_point(data = df.local, aes(x = gauge_code), shape = 16, size = 1.5, alpha = 0.9, position = position_dodge(width = dodge.width)) +
  labs(x = "", y = "H", color = "Faixa de duração") +
  scale_color_manual(values = colors, aesthetics = c("color", "fill")) + 
  theme_minimal() +
  theme(
    # legend.position = "bottom",
        legend.position = c(0.02,0.02),
        legend.justification = c(0,0),
        legend.key.spacing.y = unit(-2, "mm"),
        legend.margin = margin(l = 2, r = 2, unit = "mm"),
        legend.background = element_rect(color = "black", linewidth = 0.25),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "serif", color = "black", size = 12))

# ggsave(filename = "figuras/scale_invariance/scale_confidence_interval_block_bootstrap_regional.png",
#        width = 16, height = 20, units = "cm", dpi = 300, bg = "white")

ggsave(filename = "figuras/scale_invariance/scale_confidence_interval_block_bootstrap_regional.png",
       width = 16, height = 12, units = "cm", dpi = 300, bg = "white")


# VISUALIZAÇÃO C/ OFFSET --------------------------------------------------

rm(list = ls()); invisible(gc())

# Ler e organizar dados de expoente de escala e intervalos de confiança
df.scale.boot.offset <- arrow::read_parquet("base/gerados/scale_invariance/df_scale_block_boot_offset.pqt")

groups <- unique(df.scale.boot.offset$d)
group.names <- c("Todas", "Sub-horários", "Horários", "Diários")
groups <- setNames(object = group.names, nm = groups)

df.scale.boot.offset <- df.scale.boot.offset %>% 
  mutate(d = recode(d, !!!groups),
         d = factor(d, levels = group.names))

# Filtrar expoentes de escala 'scale'
df.scale <- df.scale.boot.offset %>% 
  filter(which_variable == "scale",
         which_moment == 1,
         d != "Todas") %>% 
  select(gauge_code, n_year, d, statistic, ci_lower, ci_upper)

# Filtrar parâmetros de posiçãõ 'offset'
df.offset <- df.scale.boot.offset %>% 
  filter(which_variable == "offset",
         which_moment == 1,
         d != "Todas") %>% 
  select(gauge_code, n_year, d, statistic, ci_lower, ci_upper)

colors <- c("Todas" = "black", "Sub-horários" = "orange", "Horários" = "steelblue2", "Diários" =  "green3")

# Ordem baseada nos scale's horários
order <- df.scale %>% 
  filter(d == "Horários") %>% 
  arrange(statistic) %>% 
  pull(gauge_code) %>% 
  unique()

# Dados locais
df.scale.local <- df.scale %>% 
  mutate(gauge_code = factor(gauge_code, levels = order)) %>% 
  filter(gauge_code != "Regional")

df.offset.local <- df.offset %>% 
  mutate(gauge_code = factor(gauge_code, levels = order)) %>% 
  filter(gauge_code != "Regional")

# Dados regionais
df.scale.regional <- df.scale %>% 
  filter(gauge_code == "Regional") %>% 
  select(-gauge_code)

df.offset.regional <- df.offset %>% 
  filter(gauge_code == "Regional") %>% 
  select(-gauge_code)

dodge.width <- 1

# Comparações locais
ggplot(mapping = aes(y = statistic, ymin = ci_lower, ymax = ci_upper, color = d)) +
  geom_errorbar(data = df.scale.local, aes(x = gauge_code), width = 0, linewidth = 0.8, alpha = 0.5, position = position_dodge(width = dodge.width)) +
  geom_point(data = df.scale.local, aes(x = gauge_code), shape = 16, size = 1.5, alpha = 0.9, position = position_dodge(dodge.width)) +
  labs(x = "", y = "H", color = "Faixa de duração") +
  scale_color_manual(values = colors, aesthetics = c("color", "fill")) + 
  theme_minimal() +
  theme(legend.position = c(0.02,0.02),
        legend.justification = c(0,0),
        legend.key.spacing.y = unit(-2, "mm"),
        legend.margin = margin(l = 2, r = 2, unit = "mm"),
        legend.background = element_rect(color = "black", linewidth = 0.25),
        legend.title = element_blank(),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        plot.background = element_rect(color = "white"),
        panel.border = element_rect(color = "black", fill = NA),
        text = element_text(family = "serif", color = "black", size = 12))
