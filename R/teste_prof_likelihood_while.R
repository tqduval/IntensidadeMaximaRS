
fit <- fgev(y, prior.info)
par <- fit$par
ll <- -fit$value
lim <- ll - qchisq(0.95, 1)
del <- 0.005
npar <- length(par)
prof.par <- c(par, ll)
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



# FUNÇÃO ------------------------------------------------------------------

# Function
roots.proflik <- function(step){
  
  # val.i <- par[[i]] + step      # iniciar contador
  val.i <- par[[i]]             # iniciar contador
  last.par <- par               # salvar último conjunto de parâmetros
  lim <- ll - qchisq(conf, 1)/2 # limite de definição dos ICs
  result <- list()              # iniciar lista com resultados
  
  # Start of 'while search'
  # while(ll > lim){
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
    
    full[i] <- val.i      # 
    val.i <- val.i + step # atualizar contador
    last.par <- full      # atualizar conjunto de parâmetros
    current.iter <- length(result) + 1
    result[[current.iter]] <- c(full, ll) # incluir resultados na lista
    
    # message(sprintf("Intervalo: %-8s | Parâmetro: %s = %0.3f | Loglik: %0.3f\n", 
    #                 ifelse(step < 0, "inferior", "superior"), names(par)[i], val.i, ll))
    
  } # end of 'while search'
  
  prof.par <- do.call(rbind, result)
  colnames(prof.par) <- c(names(par), "loglik")
  
  return(prof.par)
  
}

# prof.par <- do.call(rbind, result) %>% print()
# colnames(prof.par) <- c(names(par), "loglik")
# prof.par <- as_tibble(prof.par)

# Teste função 'roots.proflik()'
# Construção da função
# prior.info <- NULL
# prior.info <- c(0.0178, 0.0627)
d <- 20
i <- 3 # parâmetro
y <- y.id[,d]
fit <- fgev(y, prior.info)
par <- fit$par
npar <- length(par)
ll <- -fit$value
conf <- 0.95
step <- 0.0025
lim <- ll - qchisq(conf, 1)/2
ci.left <- roots.proflik(-step)
ci.right <- roots.proflik(step)

prof.par <- as_tibble(rbind(ci.left, ci.right))
prof.max <- max(prof.par$loglik)
par.prof.max <- prof.par[[i]][which.max(prof.par$loglik)]
ci.lower <- min(prof.par[[i]])
ci.upper <- max(prof.par[[i]])

# if(length(result) > 0) return(do.call(rbind, result)) else return NULL

ggplot() +
  # geom_vline(xintercept = par[i], linetype = "dashed", color = "grey50") +
  geom_segment(aes(x = par[i], y = -Inf, yend = ll), linetype = "dashed", color = "grey50") +
  geom_segment(aes(x = par.prof.max, y = -Inf, yend = prof.max), linetype = "dashed", color = "grey50") +
  geom_line(aes(x = prof.par[[i]], y = prof.par[["loglik"]]), linewidth = 0.6) +
  geom_hline(yintercept = lim, color = "red") +
  geom_point(aes(x = par[i], y = ll), fill = "white", color = "black", shape = 21) +
  geom_point(aes(x = par.prof.max, y = prof.max), fill = "white", color = "black", shape = 21) +
  geom_point(aes(x = c(ci.lower, ci.upper), y = lim), color = "red", shape = 4, size = 2.5) +
  # geom_text(
  #   aes(x = par[i], y = lim, label = sprintf("κ = %0.3f", par[i])), 
  #   family = font, size = textsize, angle = 90, vjust = -1.1, hjust = 1.2
  # ) +
  labs(
    x = names(par)[i], y = "Perfil de verossimilhança", 
    subtitle = sprintf("Duração: %g h | Priori: %s | Delta: %0.4f", d, gsub("c", "", deparse(prior.info)), step)
  ) +
  theme(panel.grid = element_blank())

# Conferir se os máximos batem
prof.max
new.par <- unlist(prof.par[which.max(prof.par$loglik), seq_len(npar)])

-nll.gev(new.par, y, prior.info)
start.gev(y)

prof.par <- prof.par %>%
  mutate(
    lag_kappa = kappa - lag(kappa),
    jumped = if_else(abs(lag_kappa) > step, TRUE, FALSE)
  )

prof.par$kappa[which(abs(prof.par$lag_kappa) > step)]
prof.par$loglik[which(abs(prof.par$lag_kappa) > step)]