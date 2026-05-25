
# PACOTES -----------------------------------------------------------------

# Limpar ambiente
rm(list = ls()); invisible(gc())

# Pacotes
if(!require(pacman)) install.packages("pacman")
pacman::p_load(pacman,    # gerenciador de pacotes
               arrow,     # leitura de dados 'parquet'
               tidyverse)

# Set theme
suppressMessages(extrafont::loadfonts())
font <- "Consolas"
fontsize <- 8
textsize <- (fontsize - 2)/.pt
set_theme(theme_minimal() +
            theme(text = element_text(color = "black", size = fontsize, family = font),
                  axis.text = element_text(color = "black"),
                  panel.border = element_rect(color = "black")))


# LER DADOS ---------------------------------------------------------------

# Intensidades máximas anuais de precipitação
df.imax <- arrow::read_parquet("base/gerados/df_imax.pqt")

# Ajuste
ls.fit <- readRDS("base/gerados/resultados_com_priori/ls_fit.rds")
df.par <- arrow::read_parquet("base/gerados/resultados_com_priori/df_par.pqt")
df.idf <- arrow::read_parquet("base/gerados/resultados_com_priori/df_idf.pqt")


# ORGANIZAÇÃO DOS DADOS ---------------------------------------------------

# durations <- unique(df.idf$d)
durations <- c(1, 2, 6, 12, 24)
durations <- c(1, 6, 24)

data <- df.imax %>% 
  filter(na_prct_yr <= 0.2, na_prct_mon <= 0.1, between(d, 1, 24)) %>% 
  group_by(gauge_code) %>% 
  filter(n_distinct(d) >= 3, n_distinct(wateryear) >= 8) %>% 
  ungroup() %>% 
  select(gauge_code, wateryear, d, imax)


## RAZÕES DE L-MOMENTOS ---------------------------------------------------

# Plotar gráficos de barras p/ razões de l-momentos
df.ratios <- data %>% 
  group_by(gauge_code, d) %>% 
  arrange(desc(imax)) %>% 
  reframe(slmom = list(lmom::samlmu(imax)),
          lcv = unlist(slmom)[[2]]/unlist(slmom)[[1]],
          lsk = unlist(slmom)[[3]],
          lkt = unlist(slmom)[[4]]) %>% 
  select(-slmom) %>% 
  pivot_longer(cols = starts_with("l"),
               names_to = "ratios",
               values_to = "value") %>% 
  mutate(label = case_when(
           ratios == "lcv" ~ "L-CV",
           ratios == "lsk" ~ "L-assimetria",
           ratios == "lkt" ~ "L-curtose"
           ),
         label = factor(label, levels = c("L-CV", "L-assimetria", "L-curtose"))) %>% arrange(gauge_code, d, ratios) %>% print()

# Gráfico de barras
plot.ratios.d <- df.ratios %>% 
  # filter(d %in% durations) %>% C
  ggplot(aes(x = factor(d), y = value)) +
  facet_wrap(~label, ncol = 1, scale = "free_y") +
  geom_boxplot(fill = "grey90", median.linewidth = 0.2, staplewidth = 0.5, outlier.size = 0.5, outlier.alpha = 0.5, linewidth = 0.2) +
  labs(x = "Duração [h]", y = "Razões de L-momentos"); plot.ratios.d

ggsave("figuras/ratios_d.png", plot.ratios.d, width = 12, height = 14, units = "cm", dpi = 300)


## CURVAS DE FREQUÊNCIA ---------------------------------------------------



# Quantis empíricos (intensidades observadas)
gg <- "A801"
data.gg <- data[data$gauge_code == gg,]
y.id.gg <- data.gg %>% 
  pivot_wider(names_from = d,
              values_from = imax,
              names_prefix = "d_") %>% 
  arrange(wateryear) %>% 
  select(-c(gauge_code, wateryear)) %>% 
  na.omit() %>% as.matrix()

df.fq.obs <- data %>% 
  filter(d %in% durations) %>% 
  group_by(gauge_code, d) %>% 
  arrange(desc(imax)) %>% 
  mutate(n = length(imax), i = row_number(),
         p_weibull = i/(n + 1), tr_weibull = 1/p_weibull,
         p_cunanne = (i - 0.4)/(n + 0.2), tr_cunanne = 1/p_cunanne,
         p_gring = (i - 0.44)/(n + 0.12), tr_gring = 1/p_gring) %>% 
  ungroup() %>% arrange(gauge_code, d, i)

n.years <- df.fq.obs %>% 
  slice_head(n = 1, by = "gauge_code") %>% 
  select(gauge_code, n)

labels <- df.fq.obs %>% 
  slice_head(n = 1, by = "gauge_code") %>% 
  select(gauge_code, n) %>% 
  mutate(label = paste(gauge_code, "|", n, "anos")) %>% 
  pull(label, name = gauge_code)

df.fq.est <- df.idf %>% 
  filter(d %in% durations)

ggplot() +
  facet_wrap(~gauge_code, ncol = 4, scale = "free_y", labeller = as_labeller(labels)) +
  geom_line(aes(x = return_period, y = q, color = factor(d)), df.fq.est, linewidth = 0.2) +
  geom_ribbon(aes(x = return_period, ymin = ci_lower, ymax = ci_upper, color = factor(d), fill = factor(d)), df.fq.est,
              alpha = 0.1, linetype = "dashed", linewidth = 0.2) +
  geom_point(aes(x = tr_gring, y = imax, color = factor(d), shape = factor(d)), df.fq.obs, size = 0.6) +
  scale_x_log10() +
  labs(x = "T [anos]", y = "Intensidade [mm/h]", color = "Duração [h]", shape = "Duração [h]", fill = "Duração [h]") +
  theme(legend.position = "bottom",
        legend.key.size = unit(fontsize + 1, "pt"),
        strip.text = element_text(color = "black", size = fontsize, hjust = 0))
  # theme(legend.position = c(0,1),
  #       legend.justification = c(0,1),
  #       legend.key.height = unit(fontsize, "pt"),
  #       legend.title = element_text(size = fontsize - 2))

ggsave("figuras/curvas_frequencia.png", width = 20, height = 28, units = "cm", dpi = 300)
# ggsave(sprintf("figuras/curvas_frequencia_%s.png", gg), width = 12, height = 8, units = "cm", dpi = 300)


# IDFs
