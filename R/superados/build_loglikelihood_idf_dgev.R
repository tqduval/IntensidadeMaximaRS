# Construção da função de log-verossimilhança do modelo de equação de IDF baseado
# na distribuição dGEV (duration-dependent Generalized Extreme Value distribution)

# Provavelmente, vai precisar fazer uma organização dos dados igual na função de
# invariância de escala, organizando em cada coluna a AMS de uma duração diferente
# pq a função vai precisar iterar em um loop das durações e das intensidades
# então, entra na primeira coluna, extrai o primeiro 'd', e os 'y_i' pro somatório


# FUNÇÕES -----------------------------------------------------------------

source("R/funcoes/fit_dgev.R")
source("R/funcoes/q_dgev.R")

rm(list = ls()); invisible(gc())
pacman::p_load(pacman, tidyverse, arrow, numDeriv, evd)

# Set theme
suppressMessages(extrafont::loadfonts())
font <- "Consolas"
fontsize <- 8
textsize <- (fontsize - 2)/.pt
set_theme(theme_minimal() +
            theme(text = element_text(color = "black", size = fontsize, family = font),
                  axis.text = element_text(color = "black"),
                  panel.border = element_rect(color = "black")))


df.imax <- arrow::read_parquet("base/gerados/df_imax.pqt")

data <- df.imax %>% 
  filter(na_prct_yr <= 0.2,
         # na_prct_mon <= 0.1,
         between(d, 1, 24)) %>% 
  group_by(gauge_code) %>% 
  filter(n_distinct(d) >= 3, n_distinct(wateryear) >= 8) %>% 
  ungroup() %>% 
  select(gauge_code, wateryear, d, imax)

gauges <- unique(data$gauge_code)
gg <- "A801"
gg <- "430350901A"
gg <- "430510802A"
gg <- "SBPA"
gg <- "A884" # falha na otimização
gg <- "431020703A"
gg <- "A803"

data.gg <- data %>% 
  filter(gauge_code == gg)

cols <- c("wateryear", "d", "imax")

teste <- data.gg %>% 
  pivot_wider(names_from = cols[1],
              values_from = cols[2],
              names_prefix = "d_")

y.id <- data.gg %>% 
  pivot_wider(names_from = d,
              values_from = imax,
              names_prefix = "d_") %>% 
  arrange(wateryear) %>% 
  select(-c(gauge_code, wateryear)) %>% 
  as.matrix() %>% 
  na.omit()

durations <- as.numeric(stringr::str_remove_all(colnames(y.id), "d_"))

ls.data <- split(data, data$gauge_code)
teste.nprior <- pbapply::pblapply(ls.data, fit.dgev, cols = names(data), method = "Nelder-Mead", scale.inv = "wide")
teste.nprior <- pbapply::pblapply(ls.data, fit.dgev, cols = names(data), method = "L-BFGS-B", scale.inv = "wide")
teste.wprior <- pbapply::pblapply(ls.data, fit.dgev, prior.info = c(0.0178, 0.0627), cols = names(data), method = "Nelder-Mead", scale.inv = "wide")
teste.wprior <- pbapply::pblapply(ls.data, fit.dgev, prior.info = c(0.0178, 0.0627), cols = names(data), method = "L-BFGS-B", scale.inv = "wide")

par.nprior.gg <- teste.nprior[[gg]]$par %>% print()
par.wprior.gg <- teste.wprior[[gg]]$par %>% print()

df.teste <- bind_rows(
  lapply(seq_along(teste), function(i){
    fit <- teste[[i]]
    gg <- names(teste[i])
    row <- tibble(gauge_code = gg, par_name = names(fit$par), par = fit$par, convergence = fit$convergence, message = fit$message,
                  value = fit$value, hessian = list(as_tibble(fit$hessian)))
  })
)


# QUANTIS -----------------------------------------------------------------

# fit <- teste[[gg]]
fit <- teste.nprior[[gg]]
fit <- teste.wprior[[gg]]


data.gg <- data[data$gauge_code == gg,]
cols <- names(data.gg)
y.id.gg <- data.gg %>% 
  pivot_wider(names_from = cols[3],
              values_from = cols[4],
              names_prefix = "d_") %>% 
  arrange(wateryear) %>% 
  select(-cols[1:2]) %>% 
  na.omit() %>% 
  as.matrix()

# Momentos
lratios <- apply(y.id.gg, 2, function(ams){
  s.lmom <- lmom::samlmu(ams)
  t2 <- s.lmom[[2]]/s.lmom[[1]] # L-CV
  t3 <- s.lmom[[3]]             # L-SK
  t4 <- s.lmom[[4]]             # L-KT
  return(c(lcv = t2, lsk = t3, lkt = t4))
})


# Posiçao de plotagem de Gringorten
freq.obs <- apply(y.id.gg, 2, function(ams){
  sorted <- sort(ams, decreasing = TRUE)
  n <- length(ams)
  i <- 1:length(ams)
  p <- (i - 0.44)/(n + 0.12)
})

rp <- c(5, 10, 50, 100, 200)
p <- 1 - 1/rp
rp <- 1/freq.obs[,1]

quantiles <- q.dgev(fit.obj = teste.wprior[[gg]], durations, return.period)

idf <- bind_cols(gauge_code = gg, quantiles) %>% 
  mutate(return_period = rp) %>% 
  pivot_longer(cols = -c(gauge_code, return_period),
               names_to = "d",
               names_transform = list(d = as.numeric),
               values_to = "q") %>% print()

# Frequencias observadas e valores
df.fq.obs <- y.id.gg %>% 
  as_tibble() %>% 
  pivot_longer(cols = where(is.numeric),
               names_prefix = "d_",
               names_to = "d",
               names_transform = list(d = as.numeric),
               values_to = "imax") %>% 
  group_by(d) %>% 
  arrange(desc(imax)) %>% 
  mutate(n = length(imax), i = row_number(),
         p_weibull = i/(n + 1), tr_weibull = 1/p_weibull,
         p_cunanne = (i - 0.4)/(n + 0.2), tr_cunanne = 1/p_cunanne,
         p_gring = (i - 0.44)/(n + 0.12), tr_gring = 1/p_gring) %>% 
  ungroup() %>% print()


# VISUALIZATION -----------------------------------------------------------

# Plotar curva de frequência
ggplot() +
  geom_line(aes(x = return_period, y = q, color = factor(d)), idf %>% mutate(return_period = as.numeric(return_period))) +
  geom_point(aes(x = tr_gring, y = imax, color = factor(d)), df.fq.obs %>% filter(d %in% durations)) +
  scale_x_log10() +
  labs(x = "T [years]", y = "Intensity [mm/h]", color = "Duration [h]", caption = sprintf("Subdaily gauge: %s", gg)) +
  theme(legend.position = c(0,1),
        legend.justification = c(0,1),
        legend.key.height = unit(fontsize, "pt"),
        legend.title = element_text(size = fontsize - 2))

# Plot IDF
ggplot() +
  geom_line(aes(x = d, y = q, color = factor(return_period)), idf, linewidth = 0.5) +
  geom_point(aes(x = d, y = imax), data.gg[data.gg$d %in% durations,], alpha = 0.4, color = "grey10", size = 0.8) +
  scale_x_log10(breaks = durations) +
  # scale_y_log10() +
  labs(x = "Duration [h]", y = "Intensity [mm/h]", color = "T [years]", caption = sprintf("Subdaily gauge: %s", gg)) +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.key.height = unit(fontsize - 1, "pt"),
        legend.title = element_text(size = fontsize - 2))

q.nprior <- lapply(durations, q.dgev, p.nexceed = p, param = par.nprior.gg)
q.wprior <- lapply(durations, q.dgev, p.nexceed = p, param = par.wprior.gg)
idf.nprior.aux <- matrix(unlist(q.nprior), length(p), length(durations), dimnames = list(rp, durations))
idf.wprior.aux <- matrix(unlist(q.wprior), length(p), length(durations), dimnames = list(rp, durations))
idf.nprior <- bind_cols(gauge_code = gg, idf.nprior.aux) %>% 
  mutate(return_period = rp) %>% 
  pivot_longer(cols = -c(gauge_code, return_period),
               names_to = "d",
               names_transform = list(d = as.numeric),
               values_to = "q") %>% print()
idf.wprior <- bind_cols(gauge_code = gg, idf.wprior.aux) %>% 
  mutate(return_period = rp) %>% 
  pivot_longer(cols = -c(gauge_code, return_period),
               names_to = "d",
               names_transform = list(d = as.numeric),
               values_to = "q") %>% print()

ggplot() +
  geom_line(aes(x = d, y = q, color = factor(return_period)), idf.nprior, linewidth = 0.5, linetype = "dashed") +
  geom_line(aes(x = d, y = q, color = factor(return_period)), idf.wprior, linewidth = 0.5, linetype = "solid") +
  # geom_point(aes(x = d, y = imax), data.gg[data.gg$d %in% durations,], alpha = 0.4, color = "grey10", size = 0.8) +
  scale_x_log10(breaks = durations) +
  scale_y_log10() +
  labs(x = "Duration [h]", y = "Intensity [mm/h]", color = "T [years]", caption = sprintf("Subdaily gauge: %s", gg)) +
  theme(legend.position = c(1,1),
        legend.justification = c(1,1),
        legend.key.height = unit(fontsize - 1, "pt"),
        legend.title = element_text(size = fontsize - 2))

ggsave("figuras/idf_compare.png", width = 12, height = 8, units = "cm", dpi = 300)

data.gg %>% 
  filter(d %in% c(1, 2, 6, 12, 24)) %>% 
  ggplot(aes(x = wateryear, y = imax, group = factor(d))) +
  geom_line() +
  labs(caption = sprintf("Gauge: %s", unique(data.gg$gauge_code))) +
  theme(legend.position = "bottom")
  


# VALIDAÇÃO ---------------------------------------------------------------

library(IDF)
IDF.agg(list(as.data.frame(df.data[df.data$gauge_code == gg,])), ds = durations, names = c("datetime", "rain_mm"), na.accept = 0.2)
xdat <- data.gg %>% rename(xdat = imax, ds = d)
fit.ulrich <- IDF::gev.d.fit(xdat = xdat$xdat, ds = xdat$ds)
par.ulrich <- IDF::gev.d.params(fit.ulrich) %>% print()
par.eu <- fit$par %>% print()
