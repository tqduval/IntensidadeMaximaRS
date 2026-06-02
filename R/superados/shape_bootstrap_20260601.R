# Estimar parâmetro de forma através de todos os métodos e avaliar
# intervalos de confiança com bootstrap percentil.

# Métodos:
# - GEV padrão: l-momentos
# - GEV padrão: mle
# - GEV padrão: rgmle
# - dGEV: mle
# - dGEV: rgmle

set.seed(1234)

# PACOTES -----------------------------------------------------------------

rm(list = ls()); invisible(gc())
pacman::p_load(pacman, tidyverse, arrow, lmom)

# Set theme
suppressMessages(extrafont::loadfonts())
font <- "Consolas"
fontsize <- 8
textsize <- (fontsize - 2)/.pt
set_theme(theme_minimal() +
            theme(text = element_text(color = "black", size = fontsize, family = font),
                  axis.text = element_text(color = "black"),
                  panel.border = element_rect(color = "black")))

# FUNÇÕES -----------------------------------------------------------------

source("R/funcoes/aux_functions.R") # funções de maximização da GEV convencional
source("R/funcoes/fdgev.R")         # funções de maximização da dGEV


# LER DADOS ---------------------------------------------------------------

# Intensidades máximas anuais
df.imax <- arrow::read_parquet("base/gerados/df_imax.pqt")

# Aplicar critérios de qualidade
data <- df.imax %>% 
  filter(na_prct_yr <= 0.2, between(d, 1, 24)) %>% 
  group_by(gauge_code) %>% 
  filter(n_distinct(d) >= 3, n_distinct(wateryear) >= 8) %>% 
  ungroup() %>% 
  select(gauge_code, wateryear, d, imax)

ls.data <- split(data, data$gauge_code)


# CONSTRUÇÃO DA FUNÇÃO ----------------------------------------------------

# gg <- "A801"
# data.gg <- data[data$gauge_code == gg,]

shape.boot <- function(data, cols = names(df.imax), B, conf = 0.95){
  
  get.shape <- function(data.gg, cols = names(data.gg)){
    
    yd.gg <- data.gg %>% 
      pivot_wider(names_from = d,
                  values_from = imax,
                  names_prefix = "d_") %>% 
      arrange(wateryear) %>% 
      select(-c(gauge_code, wateryear)) %>% 
      as.matrix() %>% 
      na.omit()
    
    n.years <- nrow(yd.gg)
    
    # Funções
    flmom <- function(ams){
      sample.lmom <- lmom::samlmu(ams)
      lmom::pelgev(sample.lmom)
    }
    
    fgev <- function(param0, y, prior.info = NULL){
      
      nll.gev <- function(param0, y, prior.info){
        
        xi <- param0[1]
        alpha <- param0[2]
        kappa <- param0[3]
        n <- length(y)
        
        if(!is.null(prior.info)){
          
          mu.kappa <- prior.info[1]
          var.kappa <- prior.info[2]^2
          
          a <- 0.5
          b <- mu.kappa+a
          p <- b^2*(1-b)/var.kappa-b
          q <- p*(1/b-1)
          
          if(a - abs(kappa) <= 1e-6) return(1e6 + abs(kappa)*1e4) # impedir log(<0)
          # if(a - abs(kappa) <= 1e-6) return(1e6) # impedir log(<0)
          l.prior <- (p - 1)*log(a + kappa) + (q - 1)*log(a - kappa)
          
        } else{
          
          l.prior <- 0
          
        }
        
        
        zi <- 1 - kappa/alpha*(y - xi)
        
        # Fallbacks
        if(alpha <= 0) return(1e6)
        if(any(zi <= 1e-8)) return(1e6 + abs(min(zi))*1e4)
        
        # Negative log-likelihood
        term1 <- (1/kappa - 1)*log(zi)
        term2 <- zi^(1/kappa)
        nll <- -(-n*log(alpha) + sum(term1 - term2) + l.prior)
        
        return(nll)
        
      }
      
      # Usando Nelder-Mead
      fit <- tryCatch({
        fit.obj <- optim(par = param0, fn = nll.gev, method = "Nelder-Mead", hessian = TRUE, y = y, prior.info = prior.info) 
      },
      error = function(e){
        warning("Erro ao maximizar a minimizar 'nll.gev(): ", e)
        return(NULL)
      })
      
      return(fit)
      
    }
    
    # Estimar parâmetros de forma
    shape.gev.lmom <- apply(yd.gg, 2, function(y){
      
      par <- flmom(y)
      kappa <- unname(par["k"])
      
    })
    
    shape.gev.mle <- apply(yd.gg, 2, function(y){
      
      start <- flmom(y)
      fit <- fgev(start, y)
      par <- fit$par
      kappa <- par[[3]]
      
    })
    
    shape.gev.rgmle <- apply(yd.gg, 2, function(y){
      
      start <- flmom(y)
      if(abs(start["k"]) >= 0.5) start["k"] <- 0.0178
      fit <- fgev(unname(start), y, prior.info = c(0.0178, 0.0627))
      par <- fit$par
      kappa <- par[[3]]
      
    })
    
    shape.dgev.mle <- fdgev(data.gg, cols = names(data.gg), method = "L-BFGS-B", scale.inv = "wide")[["par"]][["kappa"]]
    shape.dgev.rgmle <- fdgev(data.gg, prior.info = c(0.0178, 0.0627), cols = names(data.gg), "L-BFGS-B", "wide")[["par"]][["kappa"]]
    
    gg <- unique(data.gg$gauge_code)
    durations <- sort(unique(data.gg$d))
    shape.names <- c("kappa_gev_lmom", "kappa_gev_mle", "kappa_gev_rgmle", "kappa_dgev_mle", "kappa_dgev_rgmle")
    out.names <- paste0(rep(shape.names, each = length(durations)), "_d_", durations)
    out <- as.numeric(c(kappa_gev_lmom = shape.gev.lmom,
                        kappa_gev_mle = shape.gev.mle,
                        kappa_gev_rgmle = shape.gev.rgmle,
                        kappa_dgev_mle = rep(shape.dgev.mle, length(durations)),
                        kappa_dgev_rgmle = rep(shape.dgev.rgmle, length(durations))))
    
    return(setNames(out, out.names))
    
  }
  
  rep.shapes <- function(data.gg){
    
    # Reamostrar anos
    years <- unique(data.gg[[cols[2]]])
    n.years <- length(years)
    
    years.resample <- sample(years, size = n.years, replace = TRUE) # reamostrar anos
    data.gg.resample <- bind_rows(lapply(seq_along(years.resample), function(i){
      yr <- years.resample[[i]]
      resample <- subset(data.gg, data.gg[[cols[2]]] == yr)
      resample$id <- i
      return(resample)
    }))
    
    gg <- unique(data.gg$gauge_code)
    data.gg.resample[[cols[2]]] <- data.gg.resample$id
    data.gg.resample$id <- NULL
    shapes <- tryCatch(
      get.shape(data.gg = data.gg.resample),
      error = function(e){
        warning("\nErro na estação ", gg, ": ", e)
        NULL
      }
    )
    
    if(is.null(shapes)) return(NULL)
    
    return(shapes)
    
  }
  
  ls.res <- pbapply::pblapply(data, function(data.gg){
    
    gg <- unique(data.gg$gauge_code)
    original <- get.shape(data.gg)
    boot.shapes <- replicate(expr = rep.shapes(data.gg), simplify = FALSE, n = B)
    boot.shapes <- boot.shapes[!sapply(boot.shapes, is.null)]
    boot.matrix <- do.call(rbind, boot.shapes)
    ci.shapes <- apply(boot.matrix, 2, quantile, probs = c((1 - conf)/2, 1 - (1 - conf)/2), na.rm = TRUE) # parâmetros
    names <- colnames(ci.shapes)
    out <- bind_cols(gauge_code = gg, id = names, kappa = original, as_tibble(t(ci.shapes))) %>% 
      extract(col = id, 
              into = c("model", "method", "d"),
              regex = "kappa_([a-z]+)_([a-z]+)_d_([0-9]+)",
              remove = FALSE, convert = TRUE) %>% 
      rename(ci_lower = `2.5%`, ci_upper = `97.5%`)
    
    return(out)
    
  })
  
  return(bind_rows(ls.res))
  
}

df.shape <- shape.boot(ls.data, B = 3000)
arrow::write_parquet(df.shape, "base/gerados/df_shape.pqt")


# VISUALIZATION -----------------------------------------------------------

n.years <- data %>% 
  group_by(gauge_code, d) %>% 
  reframe(n_years = n_distinct(wateryear)) %>% 
  slice_head(n = 1, by = "gauge_code") %>% 
  select(gauge_code, n_years)

order <- df.shape %>% 
  mutate(n_years = n.years$n_years[match(gauge_code, n.years$gauge_code)]) %>% 
  arrange(desc(n_years)) %>% 
  mutate(label = paste(gauge_code, "|", n_years, "anos")) %>% 
  pull(label, name = n_years)

pos <- position_dodge(0.5)

df.shape %>% 
  filter(model == "dgev") %>% 
  mutate(n_years = n.years$n_years[match(gauge_code, n.years$gauge_code)],
         label = paste(gauge_code, "|", n_years, "anos"),
         label = factor(label)) %>% 
  arrange(desc(n_years), d) %>% 
  ggplot(aes(x = factor(d), y = kappa, color = method)) +
  facet_wrap(~label, ncol = 4) +
  geom_errorbar(aes(ymin = ci_lower, ymax = ci_upper), linewidth = 0.2, position = pos) +
  geom_point(size = 0.6, position = pos) +
  scale_x_discrete(breaks = c(1, 6, 12, 18, 24)) +
  labs(x = "Duração [h]", y = "κ") +
  theme(strip.text = element_text(hjust = 0))
