get.shape <- function(ams){
  sample.lmom <- lmom::samlmu(ams)
  k <- lmom::pelgev(sample.lmom)["k"]
}

rep.shape <- function(y.id.gg){
  idx <- sample(1:nrow(y.id.gg), size = nrow(y.id.gg), replace = TRUE)
  resample <- y.id.gg[idx,]
  kappa.d <- tryCatch(apply(resample, 2, get.shape),
                      error = function(e) return(NULL))
  return(kappa.d)
}

ls.data <- split(data, data$gauge_code)
B <- 2000
boot.shape <- pbapply::pblapply(ls.data, function(data.gg){
  
  gg <- unique(data.gg$gauge_code)
  
  y.id.gg <- data.gg %>% 
    pivot_wider(names_from = d,
                values_from = imax,
                names_prefix = "d_") %>% 
    arrange(wateryear) %>% 
    select(-c(gauge_code, wateryear)) %>% 
    as.matrix() %>% 
    na.omit()
  
  n.years <- nrow(y.id.gg)
  
  # rep.shape(y.id.gg)
  original <- apply(y.id.gg, 2, get.shape)
  boot.shape <- replicate(B, expr = rep.shape(y.id.gg), simplify = FALSE)
  boot.matrix <- do.call(rbind, boot.shape)
  ci <- apply(boot.matrix, 2, quantile, probs = c(0.025, 0.975), na.rm = TRUE)
  names <- colnames(ci)
  ci <- as_tibble(t(ci)) %>%
    mutate(d = as.numeric(gsub("d_", "", names))) %>% 
    rename(ci_lower = `2.5%`, ci_upper = `97.5%`)
  
  out <- bind_cols(gauge_code = gg, kappa = original, ci, n_years = n.years)
  
  return(out)
  
})

boot.shape <- bind_rows(boot.shape)

order <- boot.shape %>%
  slice_head(n = 1, by = "gauge_code") %>% 
  arrange(desc(n_years)) %>% 
  mutate(label = paste(gauge_code, "|", n_years, "anos")) %>% 
  pull(label, name = n_years)

boot.shape %>% 
  mutate(label = paste(gauge_code, "|", n_years, "anos"),
         label = factor(label, levels = order)) %>% 
  ggplot(aes(x = factor(d), y = kappa)) +
  facet_wrap(~label, ncol = 4) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), linewidth = 0.2) +
  geom_point(size = 0.6) +
  scale_x_discrete(breaks = c(1, 6, 12, 18, 24)) +
  labs(x = "Duração [h]", y = "κ", caption = "Estimativas com bootstrap e L-momentos") +
  theme(strip.text = element_text(hjust = 0))

ggsave("figuras/kappa_lmom_bootstrap.png",  width = 20, height = 28, units = "cm", dpi = 300)

ggplot(data[data$gauge_code == "SBPA",], aes(x = wateryear, y = imax, group = d)) +
  geom_line() +
  geom_point()
