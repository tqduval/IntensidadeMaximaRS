
pacman::p_load(pacman, tidyverse, arrow)

# Intensidades máximas anuais
imax.df <- arrow::read_parquet("base/gerados/df_imax.pqt")

# Aplicar critérios de qualidade
imax.df <- imax.df %>% 
  filter(na_prct_yr <= 0.2, between(d, 1, 24)) %>% 
  group_by(gauge_code) %>% 
  filter(n_distinct(d) >= 3, n_distinct(wateryear) >= 8) %>% 
  ungroup() %>% 
  select(gauge_code, wateryear, d, imax)

imax.ls <- split(imax.df, imax.df$gauge_code)


# CONSTRUÇÃO DA FUNÇÃO ----------------------------------------------------

# gg <- "A801"
# # gg <- "A844"
# data.gg <- data[data$gauge_code == gg,]

shape.boot <- function(data, cols = c("gauge_code", "year", "d", "imax"), B, conf = 0.95){
  
  if(any(names(data[[1]]) != cols)) stop("Nomes das colunas e 'cols' não corresponde aos nomes das colunas em 'data'")
  
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
    
    # L-momentos
    flmom <- function(ams){
      sample.lmom <- lmom::samlmu(ams)
      lmom::pelgev(sample.lmom)
    }

    # Estimar parâmetro de forma
    shape.gev.lmom <- apply(yd.gg, 2, function(y){
      par <- flmom(y)
      kappa <- unname(par["k"])
    })
    
    gg <- unique(data.gg$gauge_code)
    durations <- sort(unique(data.gg$d))
    shape.names <- c("kappa_gev_lmom")
    out.names <- paste0(rep(shape.names, each = length(durations)), "_d_", durations)
    out <- as.numeric(shape.gev.lmom)
    
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
    shapes <- tryCatch(get.shape(data.gg = data.gg.resample), error = function(e) NULL)
    if(is.null(shapes) || any(is.na(shapes))) return(NULL)
    
    return(shapes)
    
  }
  
  ls.res <- pbapply::pblapply(data, function(data.gg){
  # ls.res <- future.apply::future_lapply(data, function(data.gg){
    
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
    
    # return(out)
    return(list(ci = out, replicates = boot.matrix))
    
  }) # fim 'ls.res'
  # }, future.seed = TRUE) # fim 'ls.res'
  
  # return(bind_rows(ls.res))
  return(ls.res)
  
}

shape.df <- shape.boot(imax.ls, B = 10000, cols = names(imax.df))
arrow::write_parquet(shape.df, "base/gerados/df_shapes_boot_lmom.pqt")

# Run for only one gauge
gg <- "A852"
shape.ls <- shape.boot(imax.ls[gg], B = 10000, cols = names(imax.df))

replicates <- do.call(rbind, lapply(shape.ls, `[[`, "replicates"))

hist(replicates[,24])
d <- 24
qqnorm(replicates[,d]); qqline(replicates[,d])

