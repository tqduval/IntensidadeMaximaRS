
#' @title Ajustar d-GEV e estimar quantis com intervalos de confiança bootstrap
#' @details
#' Função para ajustar o modelo d-GEV para séries de intensidades máximas anuais de diferentes durações,
#' estimar quantis de chuva intensa e intervalos de confiança baseados em reamostragens bootstrap e percentis.
#' @param data lista com data.frames de intensidades máximas anuais de diferentes estações
#' @param prior.info vetor com média e desvio padrão regionais para o parâmetro de forma da GEV (ou NULL)
#' @param durations vetor com durações [h] para estimar os quantis de chuvas intensas
#' @param return.perid vetor com períodos de retorno [anos] para estimar quantis de chuvas intensas
#' @param cols nome das colunas dos data.frames contidos em `data`
#' @param method otimizador passado em `optim()` (ou 'Nelder-Mead' ou 'L-BFGS-B')
#' @param scale.inv método para gerar estimativa inicial do expoente da relação de invariância de escala (ou 'strict' ou 'wide)
#' @param B número de réplicas bootstrap (recomenda-se B > 2000)
#' @param conf nível de confiaça para intervalo de confiança com percentis
#' @param maxit número máximo de iterações, passado ao `optim()`

fitboot.dgev <- function(data, prior.info = NULL,
                         durations, return.period,
                         cols = c("gauge_code", "year", "duration", "imax"),
                         method = c("Nelder-Mead", "L-BFGS-B"),
                         scale.inv = c("wide", "strict"),
                         B = 3000,
                         conf = 0.95,
                         maxit = 1e6){
  
  #### CHECAGENS ----
  
  if(!inherits(data, "list")) stop("'data' deve ser uma lista com data.frames.")
  if(!length(prior.info) %in% c(0, 2)) stop("'prior.info' dever ser 'NULL' ou um vetor de comprimento 2.")
  
  #### FUNÇÕES ----
  
  # Ajustar d-GEV e estimar quantis
  source("R/funcoes/fdgev.R")
  source("R/funcoes/q_dgev.R")
  
  # Criar amostras, ajustar e estimar quantis
  rep.dgev <- function(data.gg, prior.info, durations, return.period, cols, method, scale.inv, maxit){
    
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
    
    gg <- unique(data.gg[[cols[1]]])                   # nome da estação
    data.gg.resample[[cols[2]]] <- data.gg.resample$id # mudar nome wateryear p/ identificador único
    
    # Ajustar d-GEV
    fit <- tryCatch(
      fdgev(data.gg.resample, prior.info, names(data.gg), method, scale.inv, maxit),
      error = function(e){
        warning("Erro na estação ", gg, ":\n", e)
        NULL
      }
    )
    
    if(is.null(fit)) return(NULL)
    
    # Quantis
    idf <- q.dgev(fit.obj = fit, durations, return.period) # estimar indensidades p/ diferentes durações e tempos de retorno
    
    q.names <- paste0(rep(return.period, length(durations)), "_", rep(durations, each = length(return.period)))
    return(c(fit$par, setNames(as.vector(idf), q.names)))
    
  }
  
  #### CÁLCULOS ----
  
  # Rodar para cada estação
  message(paste0(rep("_", getOption("width")*0.5), collapse = ""), "\n")
  message("Ajustando d-GEV, estimando quantis e realizando simulações bootstrap.")
  message(">> Número de estações: ", length(data))
  message(">> Otimizador: ", method)
  message(">> Durações [h]: ", paste(durations, collapse = ", "))
  message(">> Tempos de retorno [anos]: ", paste(return.period, collapse = ", "))
  message(">> Réplicas bootstrap: ", B)
  message(">> Nível de confiança para intervalo bootstrap percentil: ", paste0(conf*100, "%"))
  
  ls.res <- pbapply::pblapply(data, function(data.gg){
    
    # Parâmetros e quantis da série original
    gg <- unique(data.gg[[cols[1]]])
    fit <- fdgev(data.gg, prior.info, names(data.gg), method, scale.inv, maxit) # ajustar d-GEV
    # fit <- fit.dgev(data.gg, prior.info, names(data.gg), method, scale.inv, maxit) # ajustar d-GEV
    par <- fit$par                                         # parâmetros originais
    idf.aux <- q.dgev(fit, durations, return.period)       # matrix com quantis
    idf <- bind_cols(gauge_code = gg, idf.aux) %>%         # data.frame com quantis, d e tempo de retorno
      mutate(return_period = return.period) %>%
      pivot_longer(cols = -c(gauge_code, return_period),
                   names_to = "d",
                   names_transform = list(d = as.numeric),
                   values_to = "q") %>%
      mutate(id = paste0(return_period, "_", d))
    
    # Estimar intervalos de confiança bootstrap
    # Réplicas bootstrap
    boot.dgev <- replicate(expr = rep.dgev(data.gg, prior.info, durations, return.period, cols = names(data.gg), method, scale.inv, maxit),
                           simplify = FALSE, n = B)
    
    boot.dgev <- boot.dgev[!sapply(boot.dgev, is.null)] # remover iterações que retornaram NULL
    B.success <- length(boot.dgev)
    boot.matrix <- do.call(rbind, boot.dgev)            # transformar em matriz
    npar <- length(par)                                 # número de parâmetros
    nq <- length(durations)*length(return.period)       # número de quantis
    boot.par <- boot.matrix[, 1:npar, drop = FALSE]                # extrair parâmetros
    boot.q <- boot.matrix[, (npar + 1):(npar  + nq), drop = FALSE] # extrair quanti
    
    # Intervalos de confiança com percentil
    ci.par <- apply(boot.par, 2, quantile, probs = c((1 - conf)/2, 1 - (1 - conf)/2), na.rm = TRUE) # parâmetros
    par.names <- colnames(ci.par)
    ci.par <- t(ci.par) %>% 
      as_tibble() %>% 
      mutate(parname = par.names) %>% 
      rename(ci_lower = `2.5%`, ci_upper = `97.5%`)
    
    ci.q <- apply(boot.q, 2, quantile, probs = c((1 - conf)/2, 1 - (1 - conf)/2), na.rm = TRUE)     # quantis
    q.names <- colnames(ci.q)
    ci.q <- t(ci.q) %>% 
      as_tibble() %>% 
      mutate(id = q.names) %>%
      # mutate(aux = q.names) %>%
      # separate_wider_delim(aux, "_", names = c("return_period", "d")) %>% 
      # mutate(across(c(return_period, d), as.numeric)) %>% 
      rename(ci_lower = `2.5%`, ci_upper = `97.5%`)
    
    # Organizar saídas
    parameters <- bind_cols(
      gauge_code = gg, par = unname(par), ci.par, ci_type = "bootstrap",
      prior = ifelse(is.null(prior.info), FALSE, TRUE)) %>% 
      select(gauge_code, parname, par, ci_lower, ci_upper, ci_type, prior)
    
    quantiles <- idf %>% 
      mutate(ci_type = "bootstrap", 
             prior = ifelse(is.null(prior.info), FALSE, TRUE)) %>% 
      left_join(ci.q, by = c("id")) %>% select(-id)
    
    out <- list(par = parameters,
                idf = quantiles,
                replicates = boot.matrix,
                B.in.out = c(B, B.success))
    
    return(out)
    
  })
  
  message(paste0(rep("_", getOption("width")*0.5), collapse = ""), "\n")
  
  return(ls.res)
  
}
