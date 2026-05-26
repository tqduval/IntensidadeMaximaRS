fit.dgev <- function(data, prior.info = NULL,
                     durations, return.period,
                     cols = c("gauge_code", "year", "duration", "imax"),
                     method = c("Nelder-Mead", "L-BFGS-B"),
                     scale.inv = c("wide", "strict"),
                     conf = 0.95,
                     which.ci = c("delta", "bootstrap"),
                     B = NULL,
                     maxit = 1e6){
  
  which.ci <- match.arg(which.ci, c("delta", "bootstrap"))
  
  source("R/funcoes/fdgev.R")
  source("R/funcoes/q_dgev.R")
  
  ls.res <- switch(
    which.ci,
    "delta" = {
     
      source("R/funcoes/fitdelta_dgev.R")
      ls.fit <- fitdelta.dgev(data, prior.info, durations, return.period, cols, method, scale.inv, conf, maxit)
      
    },
    "bootstrap" = {
      
      source("R/funcoes/fitboot_dgev.R")
      if(is.null(B)) stop("Para intervalos de confiança com bootstrap, informe o número de réplicas 'B'.")
      ls.fit <- fitboot.dgev(data, prior.info, durations, return.period, cols, method, scale.inv, B, conf, maxit)
      
    })
  
  return(list(parameters = bind_rows(lapply(ls.fit, `[[`, "par")),
              quantiles = bind_rows(lapply(ls.fit, `[[`, "idf"))))
  
}

prior.info <- c(0.0178, 0.0627)
durations <- c(1, 2, 4, 6, 8, 10, 12, 18, 24)
return.period <- c(1.1, 1.5, 2, 5, 10, 20, 50, 100, 200, 500)
method <- "L-BFGS-B"
scale.inv <- "wide"
B <- 3000

teste <- fit.dgev(data = ls.data[1:5],
                  prior.info = prior.info, 
                  durations = durations, 
                  return.period = return.period,
                  cols = names(ls.data[[1]]), 
                  method = method, 
                  scale.inv = scale.inv,
                  which.ci = "bootstrap",
                  B = 3000)