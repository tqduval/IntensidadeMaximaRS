
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
fitboot.dgev <- function(
  data, prior.info = NULL,
  durations, return.period,
  cols = c("gauge_code", "year", "duration", "imax"),
  method = c("Nelder-Mead", "L-BFGS-B"),
  scale.inv = c("wide", "strict"),
  B = 3000,
  conf = 0.95,
  maxit = 1e6
){

  #### CHECKINGS ----

  if(!inherits(data, "list")) stop("'data' deve ser uma lista com data.frames.")
  if(!length(prior.info) %in% c(0, 2)) stop("'prior.info' dever ser 'NULL' ou um vetor de comprimento 2.")

  #### INTERNAL FUNCTIONS ----

  # Fit d-GEV and estimate quantiles
  source("R/funcoes/fdgev.R")

  # Estimate confidence intervals with Bootstrap percentile method
  boot.dgev <- function(data.gg){

    # Resampling framework
    y.id <- data.gg %>%                 # create a matrix year x d, filled with imax
      pivot_wider(names_from = cols[3], # transformar durações em colunas
                  values_from = cols[4],
                  names_prefix = "d_") %>%
      arrange(cols[2]) %>%
      select(-cols[1:2]) %>%
      na.omit() %>% as.matrix()

    N <- nrow(y.id)

    # Generate a matrix of B x n.years
    idx.samples <- matrix(sample(seq_len(N), size = N*B, replace = TRUE), nrow = B, ncol = N) # resampled indices
    idx.original <- seq_len(N)                   # original indices
    indices <- rbind(unname(idx.original), idx.samples)  # bind original indices

    # Fit and estimate quantiles
    boot.dgev <- lapply(seq_len(B + 1), function(i){
      idx <- indices[i,] # extract row with indices, where i = 1 is the original data
      fit <- tryCatch( 
        fdgev(y.id[idx, ], prior.info, cols, method, scale.inv, maxit),
        error = function(e) NULL
      )
      if(is.null(fit)) return(NULL)
      idf <- q.dgev(fit, durations, return.period) # matrix with quantiles
      q.names <- paste0(rep(return.period, length(durations)), "_", rep(durations, each = length(return.period)))
      return(c(fit$par, setNames(as.vector(idf), q.names)))
    }) # end of 'bood.dgev'

    # Organize results
    original <- boot.dgev[[1]]    # results from original data
    boot.samples <- boot.dgev[-1] # resulta from sampled data
    boot.samples <- boot.samples[!sapply(boot.samples, is.null)] # remove NULL results
    B.success <- length(boot.samples)           # number of 'good' replicates
    boot.matrix <- do.call(rbind, boot.samples) # build matrix
    npar <- 5                                   # number of d-GEV parameters
    boot.par <- boot.matrix[,1:npar, drop = FALSE] # parameter matrix
    boot.q <- boot.matrix[,(npar + 1):(npar + length(durations)*length(return.period))] # quantile  matrix

    # Process original estimates
    gg <- unique(data.gg[[cols[1]]])
    par <- original[1:npar]
    idf <- original[(npar + 1):(npar + length(durations)*length(return.period))]

    # Determine percentile confidence invervals
    ci.par <- apply(boot.par, 2, quantile, probs = c((1-conf)/2, 1 - (1-conf)/2), na.rm = TRUE) # parameter CI
    ci.q <- apply(boot.q, 2, quantile, probs = c((1-conf)/2, 1 - (1-conf)/2), na.rm = TRUE)     # quantile CI

    df.ci.par <- t(ci.par) %>% # organize as tbl_df
      as_tibble(.name_repair = "minimal") %>% 
      setNames(c("ci_lower", "ci_upper")) %>% 
      mutate(parname = colnames(ci.par))

    df.ci.q <- t(ci.q) %>% # organize as tbl_df
      as_tibble(.name_repair = "minimal") %>% 
      setNames(c("ci_lower", "ci_upper")) %>% 
      mutate(id = colnames(ci.q)) %>% 
      separate(id, into = c("return_period", "d"), convert = TRUE)

    # Setup outputs
    parameters <- bind_cols(
      gauge_code = gg,
      par = unname(par),
      df.ci.par,
      ci_type = "bootstrap",
      prior = ifelse(is.null(prior.info), FALSE, TRUE),
      which_prior = deparse(prior.info)
    ) %>% 
      select(gauge_code, parname, par, ci_lower, ci_upper, ci_type, prior, which_prior)

    quantiles <- bind_cols(
      gauge_code = gg,
      q = idf,
      df.ci.q,
      ci_type = "bootstrap",
      prior = ifelse(is.null(prior.info), FALSE, TRUE),
      which_prior = deparse(prior.info)
    ) %>% 
      select(gauge_code, d, return_period, q, ci_lower, ci_upper, ci_type, prior, which_prior)

    out <- list(
      par = parameters,
      idf = quantiles,
      replicates = boot.matrix,
      success = sprintf("%g/%g", B.success, B)
    )

    return(out)

  }

  #### RUN ----

  # Messages
  message("Fitting and estimating quantiles. Confidence intervals method: Boostrap.")
  message(">> Number of gauges: ", length(data))
  message(">> Optimizer method: ", method)
  message(">> Durations [h]: ", paste(durations, collapse = ", "))
  message(">> Return period [years]: ", paste(return.period, collapse = ", "))
  message(">> Bootstrap replicates: ", B)
  message(">> Confidence level: ", paste0(conf*100, "%"))

  ls.res <- pbapply::pblapply(data, boot.dgev)

  return(ls.res)

}