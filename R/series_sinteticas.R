
pacman::p_load(pacman, tidyverse, arrow, evd, lmom, extRemes)
rm(list = ls()); invisible(gc())
set.seed(1234)

# Intensidades máximas anuais
df.imax <- arrow::read_parquet("base/gerados/df_imax.pqt")

# Aplicar critérios de qualidade
data <- df.imax %>%
  filter(na_prct_yr <= 0.2, between(d, 1, 24)) %>%
  group_by(gauge_code) %>%
  filter(n_distinct(d) >= 3, n_distinct(wateryear) >= 8) %>%
  ungroup() %>%
  select(gauge_code, wateryear, d, imax)

# Selecionar estação p/ análise
gg <- "A812" # uma das estações problemáticas
d <- 20      # duração problemática
data.gg <- data[data$gauge_code == gg,]

y.id <- data.gg %>%
  pivot_wider(names_from = d, values_from = imax, names_prefix = "d_") %>%
  arrange(wateryear) %>%
  select(-c(gauge_code, wateryear)) %>%
  as.matrix() %>%
  na.omit()

N <- nrow(y.id)
y <- y.id[,d]

# Estimar L-momentos
sample.mom <- lmom::samlmu(y)
lcv <- sample.mom[[2]]/sample.mom[[1]]
lsk <- sample.mom[[3]]
par0 <- lmom::pelgev(sample.mom)
par0["k"] <- -par0["k"] # cauda leve

n.reps <- 10000 # Martins e Stedinger (2000)
gev.samples <- replicate(n.reps, expr = gev.sam <- evd::rgev(N, loc = par0[1], scale = par0[2], shape = par0[3]), simplify = FALSE)
densidades <- lapply(gev.samples, density)
xlim <- range(sapply(densidades, function(dens) dens$x))
ylim <- range(sapply(densidades, function(dens) dens$y))

plot(density(y), # densidade original
     xlim = xlim,
     ylim = ylim,
     main = "",
     ylab = "Densidade",
     xlab = paste("Intensidades de d =", d, "h"),
     col = adjustcolor("darkblue"),
     lwd = 2)

polygon(density(y), col = rgb(0, 0, 1, 0.3), border = "darkblue")

title("Densidades geradas com séries curtas", adj = 0)

for(i in 1:length(densidades)){
  lines(densidades[[i]],
        col = adjustcolor("darkblue",  alpha.f = 0.1),
        lwd = 1.5)
}

# Investigar ajustes com método delta
simulations <- pbapply::pblapply(gev.samples, function(y){
  
  tryCatch({
    sample.lmom <- lmom::samlmu(y)
    lcv <- sample.lmom[[2]]/sample.mom[[1]]
    lsk <- sample.lmom[[3]]
    fit <- extRemes::fevd(y, type = "GEV", method = "MLE")
    ci.delta <- ci(fit, alpha = 0.05, method = "normal", type = "parameter", which.par = 3)
    
    out <- c(lcv = lcv, lsk = lsk, shape = ci.delta[[2]], ci_lower = ci.delta[[1]], ci_upper = ci.delta[[3]], width = ci.delta[[3]] - ci.delta[[1]])
  }, error = function(e) NULL)
  
})

mc.sims <- do.call(rbind, simulations)
