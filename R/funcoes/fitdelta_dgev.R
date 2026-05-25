source("R/funcoes/fit_dgev.R")
source("R/funcoes/q_dgev.R")


gg <- "A801"
data.gg <- data[data$gauge_code == gg,]
durations <- unique(data.gg$d)
prior.info <- c(0.0178, 0.0627)
method <-  "L-BFGS-B"
conf <- 0.95

fit <- fit.dgev(data.gg, prior.info, names(data.gg), method, "wide")
par <- fit$par
hes <- fit$hessian

# Intervalos de confiança parâmetros
# Extrair elementos das diagonais
cov.par <- try(solve(hes), silent = TRUE)
if(inherits(cov.par, "try-error")) return(matrix(NA, length(par), length(par)))
var.par <- diag(cov.par)
sd.par <- sqrt(var.par)

ci.par.l <- par - qnorm(1 - (1 - conf)/2)*sd.par
ci.par.u <- par + qnorm(1 - (1 - conf)/2)*sd.par

# Intervallos de confiança quantis
y.id.gg <- data.gg %>% 
  pivot_wider(names_from = d,
            values_from = imax,
            names_prefix = "d_") %>% 
  arrange(wateryear) %>% 
  select(-c(gauge_code, wateryear)) %>% 
  as.matrix() %>% 
  na.omit()

# Versão vetorizada de q.all() em q.dgev()
q.dgev.vec <- function(param, d, p){
  
  xit <- param[1]; alpha0 <- param[2]; kappa <- param[3]; H <- param[4]; theta <- param[5]
  A <- alpha0/(d + theta)^H
  B <- xit - (1/kappa)*(1 - (-log(p))^kappa)
  
  yp <- A*B
  
  return(unname(yp))
  
}

grad.q.dgev <- function(param, d, p){

  xit <- param[1]; alpha0 <- param[2]; kappa <- param[3]; H <- param[4]; theta <- param[5]
  logp <- -log(p)
  A <- alpha0/(d + theta)^H
  B <- xit - (1/kappa)*(1 - logp^kappa)
  yp <- A*B
  
  dy.dxit <- A
  dy.dalpha0 <- B/(d + theta)^H
  dy.dkappa <- A/kappa^2*(1 + logp^kappa*(kappa*log(logp) - 1))
  dy.dH <- -yp*log(d + theta)
  dy.dtheta <- -yp*H/(d + theta)
  
  rp <- 1/(1 - p)
  grad.names <- paste0(names(param), "_d", rep(d, each = length(param)), "_rp", rp)
  return(setNames(c(dy.dxit, dy.dalpha0, dy.dkappa, dy.dH, dy.dtheta), grad.names))
  
}

return.period <- c(1.05, 1.1, 1.5, 2, 5, 10, 20, 30, 50, 100, 200, 500)
probs <- 1 - 1/return.period

grad.q.dgev(par, d = 1, p = 0.9)

grid <- expand.grid(d = durations, p = probs)

d <- durations[1]
p <- probs[1]

q <- q.dgev.vec(par, d, p)
grad <- grad.q.dgev(par, d, p)

var.q <- t(grad) %*% cov.par %*% grad
sd.q <- sqrt(var.q)

ci.q.l <- q - qnorm(1 - (1 - conf)/2)*sd.q
ci.q.u <- q + qnorm(1 - (1 - conf)/2)*sd.q

out <- c(d = d, return_period = 1/(1 - p), q = q, sd = sd.q, ci_lower = ci.q.l, ci_upper = ci.q.u)

ls.ci.q <- lapply(1:nrow(grid), function(i){
  
  d <- grid[i, "d"]
  p <- grid[i, "p"]
  
  q <- q.dgev.vec(par, d, p)
  grad <- grad.q.dgev(par, d, p)
  
  var.q <- t(grad) %*% cov.par %*% grad
  sd.q <- sqrt(var.q)
  
  ci.q.l <- q - qnorm(1 - (1 - conf)/2)*sd.q
  ci.q.u <- q + qnorm(1 - (1 - conf)/2)*sd.q
  
  out <- c(d = d, return_period = 1/(1 - p), q = q, sd = sd.q, ci_lower = ci.q.l, ci_upper = ci.q.u)
  
})

