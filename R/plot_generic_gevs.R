library(tidyverse)
library(evd)

# Set theme
suppressMessages(extrafont::loadfonts())
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
      plot.margin = margin(2, 0.5, 0.5, 0.5)
    )
)

shape <- c(-0.2, 0, 0.3)
loc <- 55
scale <- 20
x <- seq(3, 220, 0.5)

limit <- loc + scale/0.2

df <- map_dfr(shape, function(k){
  tibble(
    kappa = factor(k),
    x = x,
    density = evd::dgev(x, loc, scale, k)
  )
}) %>% 
  mutate(density = if_else(x > limit & kappa == -0.2, NA_real_, density))

ggplot(df, aes(x = x, y = density, color = kappa, linetype = kappa)) +
  geom_vline(xintercept = limit, linewidth = 0.2, linetype = "dashed", color = "grey80") +
  geom_hline(yintercept = 0, linewidth = 0.2, color = "grey80") +
  geom_line() +
  scale_color_manual(name = NULL, values = c("black", "grey30", "grey70"), labels = paste("κ =", shape*(-1))) +
  scale_linetype_manual(name = NULL, values = c("solid", "dashed", "solid"), labels = paste("κ =", shape*(-1))) +
  labs(y = "Densidade") +
  theme(
    legend.position = c(0.9, 0.2),
    legend.justification = c(1, 0),
    legend.key.width = unit(16, "pt")
  )

ggsave(
  "figuras/preliminares/gevs_introducao.png", width = 10, height = 7, units = "cm", dpi = 600
)

