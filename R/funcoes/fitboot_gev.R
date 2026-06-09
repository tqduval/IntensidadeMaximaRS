
# ALTERAÇÕES: reestruturar esta função como fitboot.dgev()
# Ajustar GEV por durações e estimar quantis e intervalos de confiança com bootstrap
fitboot.gev <- function(
  data, prior.info = NULL,
  durations, return.period,
  cols = c("gauge_code", "year", "duration", "imax"),
  B = 3000,
  conf = 0.95,
  maxit = 1e6
){
  
  #### CHECAGENS ----
  
  if(!inherits(data, "list")) stop("'data' deve ser uma lista com data.frames.")
  if(!length(prior.info) %in% c(0, 2)) stop("'prior.info' dever ser 'NULL' ou um vetor de comprimento 2.")
  
  #### FUNÇÕES ----

  # Ajustar d-GEV e estimar quantis
  source("R/funcoes/fgev.R")
  
  # Criar amostras, ajustar e estimar quantis
  rep.gev <- function(resample = TRUE, data.gg, prior.info, durations, return.period, cols, maxit){
    
    gg <- unique(data.gg[[cols[1]]])
    y.id <- data.gg %>% 
      pivot_wider(names_from = cols[3],
                  values_from = cols[4],
                  names_prefix = "d_") %>% 
      arrange(cols[2]) %>% 
      select(-c(cols[1], cols[2])) %>% 
      as.matrix() %>% 
      na.omit()
    
    # Reamostrar anos
    if(isTRUE(resample)){
      years <- unique(data.gg[[cols[2]]])
      n.years <- length(years)
      idx <- sample(1:n.years, size = n.years, replace = TRUE)
      y.id <- y.id[idx,]
    }
    
    # Ajustar GEV por duração e estimar quantis
    res <- apply(y.id, 2, function(y){

      fit <- fgev(unname(start), y, prior.info)

      if(is.null(fit) || is.null(fit$par)) {
        return(rep(NA_real_, 3 + length(return.period)))
      }

      par <- fit$par
      p <- 1 - 1/return.period
      idf <- qgev(par, p)
      out <- c(par, idf)
    })

    ds <- gsub("d_", "", colnames(res))
    out.names <- c("xi", "alpha", "kappa", return.period)
    n.out <- length(out.names)
    res.names <- paste0(out.names, "_", rep(ds, each = n.out))

    return(setNames(as.vector(res), res.names))
    
  } # fim 'rep.gev()'
  
  #### CÁLCULOS ----
  
  # Rodar para cada estação
  message("Ajustando GEV por durações, estimando quantis e realizando simulações bootstrap.")
  message(">> Número de estações: ", length(data))
  message(">> Durações [h]: ", paste(durations, collapse = ", "))
  message(">> Tempos de retorno [anos]: ", paste(return.period, collapse = ", "))
  message(">> Réplicas bootstrap: ", B)
  message(">> Nível de confiança para intervalo bootstrap percentil: ", paste0(conf*100, "%"))
  
  ls.res <- pbapply::pblapply(data, function(data.gg){
    
    # Parâmetros e quantis da série original
    gg <- unique(data.gg[[cols[1]]])

    original <- rep.gev(resample = FALSE, data.gg, prior.info, durations, return.period, cols, maxit)
    df.original <- enframe(original, name = "name", value = "value") %>% 
      separate(name, into = c("var", "d"), sep = "_") %>% 
      mutate(d = as.integer(d))

    # Organizar estimativas originais
    par.names <- c("xi", "alpha", "kappa")

    df.par <- bind_cols(gauge_code = gg, df.original) %>% 
      filter(var %in% par.names) %>% 
      rename(parname = var, par = value) %>% 
      select(gauge_code, d, parname, par)

    df.q <- bind_cols(gauge_code = gg, df.original) %>% 
      filter(!var %in% par.names) %>% 
      rename(return_period = var, q = value) %>% 
      select(gauge_code, d, return_period, q) %>%
      mutate(return_period = as.double(return_period))
    
    # Estimar intervalos de confiança bootstrap
    # Réplicas bootstrap
    boot.gev <- replicate(
      n = B,
      expr = tryCatch({
        rep.gev(resample = TRUE, data.gg, prior.info, durations, return.period, cols = cols, maxit)
      }, error = function(e){
        return(NULL)
      }),
      simplify = FALSE
    )
    
    boot.gev <- boot.gev[!sapply(boot.gev, is.null)] # remover iterações que retornaram NULL
    B.success <- length(boot.gev)                      # nro. réplicas boas
    boot.matrix <- do.call(rbind, boot.gev)            # transformar em matriz
    names <- colnames(boot.matrix)
    pattern <- paste0("^(", paste(par.names, collapse = "|"), ")_")
    par.idx <- grepl(pattern, names)
    boot.par <- boot.matrix[,par.idx, drop = FALSE]
    boot.q <- boot.matrix[,!par.idx, drop = FALSE]
    
    # Intervalos de confiança com percentil
    ci.par <- apply(boot.par, 2, quantile, probs = c((1 - conf)/2, 1 - (1 - conf)/2), na.rm = TRUE) # parâmetros
    par.names <- colnames(ci.par)
    df.ci.par <- t(ci.par) %>% 
      as_tibble() %>% 
      mutate(parname = par.names) %>% 
      separate(parname, into = c("parname", "d"), sep = "_", convert = TRUE) %>% 
      rename(ci_lower = `2.5%`, ci_upper = `97.5%`)

    
    ci.q <- apply(boot.q, 2, quantile, probs = c((1 - conf)/2, 1 - (1 - conf)/2), na.rm = TRUE)     # quantis
    q.names <- colnames(ci.q)
    df.ci.q <- t(ci.q) %>% 
      as_tibble() %>% 
      mutate(return_period = q.names) %>% 
      separate(return_period, into = c("return_period", "d"), sep = "_", convert = TRUE) %>% 
      rename(ci_lower = `2.5%`, ci_upper = `97.5%`)

    parameters <- df.par %>% 
      left_join(df.ci.par, by = c("d", "parname")) %>% 
      mutate(prior = deparse(prior.info), ci_type = "bootstrap")
    
    quantiles <- df.q %>% 
      left_join(df.ci.q, by = c("d", "return_period")) %>% 
      mutate(prior = deparse(prior.info), ci_type = "bootstrap")
    
    out <- list(par = parameters,
                idf = quantiles,
                replicates = boot.matrix,
                B.in.out = c(B, B.success))
    
    return(out)
    
  })
  
  return(ls.res)
  
}
