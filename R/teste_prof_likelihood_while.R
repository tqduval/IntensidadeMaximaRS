
fit <- fgev(y, prior.info)
par <- fit$par
ll <- -fit$value
del <- 0.005
npar <- length(par)
prof.par <- c(par, -nll.min)
names(prof.par) <- c(names(par), "loglik")

i <- 3
last.par <- par
val.i <- par[[i]]


# Search left
while(ll > lim){

  fit.prof <- tryCatch(
    optim(
      last.par[-i], fn = prof.nll.gev, gr = prof.grad.gev, method = "L-BFGS-B", lower = lower, upper = upper,
      which.par = i, par.value = val.i, prior.info = prior.info, y = y,
      control = list(maxit = maxit, factr = 1e9, pgtol = 1e-8)
    ),
    error = function(e) NULL
  )

  if(!is.null(fit.prof) && fit.prof$convergence == 0){
    ll <- -fit.prof$value
  } else{
    fit.prof <- tryCatch(
      optim(
        last.par[-i], fn = prof.nll.gev, gr = NULL, method = "Nelder-Mead",
        which.par = i, par.value = val.i, prior.info = prior.info, y = y,
        control = list(maxit = maxit)
      ),
      error = function(e) NULL
    )
    if(!is.null(fit.prof) && fit.prof$convergence == 0){
    ll <- -fit.prof$value
      if(ll < -1e6){
        ll <- NA_real_
      }
    } else{
      ll <- NA_real_
    }
  } # end of 'search left'

  full <- numeric(npar)
  if(!is.null(fit.prof) && !is.na(ll)){
    full[-i] <- fit.prof$par # se não falhar, atualizar parâmetros
  } else{
    full[-i] <- last.par[-i]  # se falhar, usar rodada anterior
  }
  full[i] <- val.i
  prof.aux <- c(full, ll)
  prof.par <- rbind(prof.par, prof.aux)

  val.i <- val.i - del
  message(sprintf("Valor da variável: %0.3f | Valor da verossimmilhança: %0.3f", val.i, ll))

}

# Search right
while(ll > lim){

  fit.prof <- tryCatch(
    optim(
      last.par[-i], fn = prof.nll.gev, gr = prof.grad.gev, method = "L-BFGS-B", lower = lower, upper = upper,
      which.par = i, par.value = val.i, prior.info = prior.info, y = y,
      control = list(maxit = maxit, factr = 1e9, pgtol = 1e-8)
    ),
    error = function(e) NULL
  )

  if(!is.null(fit.prof) && fit.prof$convergence == 0){
    ll <- -fit.prof$value
  } else{
    fit.prof <- tryCatch(
      optim(
        last.par[-i], fn = prof.nll.gev, gr = NULL, method = "Nelder-Mead",
        which.par = i, par.value = val.i, prior.info = prior.info, y = y,
        control = list(maxit = maxit)
      ),
      error = function(e) NULL
    )
    if(!is.null(fit.prof) && fit.prof$convergence == 0){
    ll <- -fit.prof$value
      if(ll < -1e6){
        ll <- NA_real_
      }
    } else{
      ll <- NA_real_
    }
  }

  full <- numeric(npar)
  if(!is.null(fit.prof) && !is.na(ll)){
    full[-i] <- fit.prof$par # se não falhar, atualizar parâmetros
  } else{
    full[-i] <- last.par[-i] # se falhar, usar rodada anterior
  }
  full[i] <- val.i
  prof.aux <- c(full, ll)
  prof.par <- rbind(prof.par, prof.aux)

  val.i <- val.i + del
  message(sprintf("Valor da variável: %0.3f | Valor da verossimmilhança: %0.3f", val.i, ll))

} # end of search right

# Construção da função
fit <- fgev(y, prior.info)
par <- fit$par
ll <- -fit$value
lim <- ll - qchisq(0.95, 1)
step <- +0.0025

# Function
val.i <- par[[i]] + step # iniciar contador
last.par <- par          # salvar último conjunto de parâmetros
result <- list()         # iniciar lista com resultados

# Start of 'while search'
while(ll > lim){

  fit.prof <- tryCatch(
    optim(
      last.par[-i], fn = prof.nll.gev, gr = prof.grad.gev, method = "L-BFGS-B", lower = lower, upper = upper,
      which.par = i, par.value = val.i, prior.info = prior.info, y = y,
      control = list(maxit = maxit, factr = 1e9, pgtol = 1e-8)
    ),
    error = function(e) NULL
  )

  if(!is.null(fit.prof) && fit.prof$convergence == 0){
    ll <- -fit.prof$value
  } else{
    fit.prof <- tryCatch(
      optim(
        last.par[-i], fn = prof.nll.gev, gr = NULL, method = "Nelder-Mead",
        which.par = i, par.value = val.i, prior.info = prior.info, y = y,
        control = list(maxit = maxit)
      ),
      error = function(e) NULL
    )
    if(!is.null(fit.prof) && fit.prof$convergence == 0){
    ll <- -fit.prof$value
      if(ll < -1e6){
        ll <- NA_real_
      }
    } else{
      ll <- NA_real_
    }
  }

  full <- numeric(npar)
  if(!is.null(fit.prof) && !is.na(ll)){
    full[-i] <- fit.prof$par # se não falhar, atualizar parâmetros
  } else{
    full[-i] <- last.par[-i]  # se falhar, usar rodada anterior
  }

  full[i] <- val.i
  val.i <- val.i + step
  last.par <- full
  current.iter <- length(result) + 1
  result[[current.iter]] <- c(full, ll)

  cat(sprintf(
    "Intervalo: %-8s | Parâmetro: %s = %0.3f | Loglik: %0.3f\n", 
    ifelse(step < 0, "inferior", "superior"), names(par)[i], val.i, ll))

} # end of 'while search'

prof.par <- do.call(rbind, result) %>% print()
colnames(prof.par) <- c(names(par), "loglik")
prof.par <- as_tibble(prof.par)

# if(length(result) > 0) return(do.call(rbind, result)) else return NULL

ggplot() +
  geom_vline(xintercept = par[i], linetype = "dashed", color = "grey50") +
  geom_line(aes(x = prof.par[[i]], y = prof.par[["loglik"]]), linewidth = 0.8) +
  geom_hline(yintercept = lim, color = "red") +
  geom_text(
    aes(x = par[i], y = lim, label = sprintf("κ = %0.3f", par[i])), 
    family = font, size = textsize, angle = 90, vjust = -1.1, hjust = 1.2
  ) +
  labs(
    x = names(par)[i], y = "Perfil de verossimilhança", 
    subtitle = sprintf("Duração: %g h | Priori: %s | Delta: %0.03f", d, gsub("c", "", deparse(prior.info)), del)
  ) +
  theme(panel.grid = element_blank())
