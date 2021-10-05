
create_pars <- function(fish_dat, mean_covs, mean_covs_log, mean_covs_logit, catch_covs){
  
  K = nlevels(fish_dat$COMMON_NAME)
  levs = levels(fish_dat$COMMON_NAME)
  lake_index = (fish_dat %>% filter(COMMON_NAME == levs[1]))$DOW # lake_index == lake_id[2]
  lake_id = levels(fish_dat$DOW)
  n_lakes = length(levels(fish_dat$DOW))
  
  pars = list()
  
  # data
  pars$Y = list()
  pars$X = list()
  pars$Z = list()
  pars$effort = list()
  
  for(k in 1:K){
    pars$Y[[k]] = (fish_dat %>% filter(COMMON_NAME == levs[k]) %>% select(TOTAL_CATCH))$TOTAL_CATCH
  }
  
  for(k in 1:K){
    X = fish_dat %>% 
      filter(COMMON_NAME == levs[k]) %>% 
      select(all_of(mean_covs)) %>% 
      mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
      mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
      mutate_at(vars(all_of(mean_covs_log)), ~ log(.))
    # %>% 
    #   mutate(Int = 1)
    
    Z = fish_dat %>% 
      filter(COMMON_NAME == levs[k]) %>%
      select(all_of(catch_covs), GN) %>% 
      mutate_at(vars(all_of(catch_covs)), .funs = list(GN = ~.*GN))
    
    pars$X[[k]] = as.matrix(X)
    pars$Z[[k]] = as.matrix(Z)
  }
  
  for(k in 1:K){
    pars$effort[[k]] = (fish_dat %>% filter(COMMON_NAME == levs[k]) %>% select(EFFORT))$EFFORT
  }
  
  
  # parameters
  pars$n = unlist(lapply(pars$Y, length))
  pars$p_beta = ncol(pars$X[[1]])
  pars$p_phi = ncol(pars$Z[[1]])
  pars$K = K
  
  pars$beta_0 = 0
  pars$beta_0_accept = 0
  pars$beta_0_prior_var = 100
  pars$beta = array(0, dim = c(K, pars$p_beta))
  pars$beta_accept =  array(0, dim = c(K, pars$p_beta))
  pars$beta_prior_var = 100
  pars$phi = array(0, dim = c(K, pars$p_phi))
  pars$phi_accept =  array(0, dim = c(K, pars$p_phi))
  pars$phi_prior_var = 100
  
  # hyperpriors
  pars$Sigma_species = diag(K)
  pars$nu_species = K + 10
  pars$Psi_species = diag(K)
  
  pars$omega = matrix(rep(0, K), ncol = K)
  
  # Proposal variances
  pars$sig_prop_beta = array(2, dim = c(K, pars$p_beta))
  pars$sig_prop_phi = array(2, dim = c(K, pars$p_phi))
  pars$sig_prop_beta_0 = 2
  
  # indexing
  pars$lake_index = lake_index
  pars$lake_id = lake_id
  pars$n_lakes = n_lakes
  pars$fish_names = levels(fish_dat$COMMON_NAME)
  
  return(pars)
  
}

update_beta <- function(pars){
  
  # data
  Y = pars$Y
  X = pars$X
  Z = pars$Z
  effort = pars$effort
  omega = pars$omega
  phi = pars$phi
  
  # parameters
  n = pars$n
  p_beta = pars$p_beta
  p_phi = pars$p_phi
  K = pars$K
  
  # beta monitor values
  beta_accept = array(0, dim = c(K, p_beta))
  beta_curr = pars$beta
  sig_prop_beta = pars$sig_prop_beta
  beta_prior_var = pars$beta_prior_var
  
  # beta_0 monitor
  beta_0 = pars$beta_0
  beta_0_accept = 0
  beta_0_prior_var = pars$beta_0_prior_var
  sig_prop_beta_0 = pars$sig_prop_beta_0
  
  # beta_0
  
  b_prop = rnorm(1, beta_0, sig_prop_beta_0)
  
  
  like_curr = like_prop = 0
  for(k in 1:K){
    like_curr = like_curr + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta_curr[k,] + Z[[k]] %*% phi[k,] + omega[k]), log = T)) + dnorm(beta_0, 0, sig_prop_beta_0, log = T)
    like_prop = like_prop + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(b_prop + X[[k]] %*% beta_curr[k,] + Z[[k]] %*% phi[k,] + omega[k]), log = T)) + dnorm(b_prop, 0, sig_prop_beta_0, log = T)
  }
  # if(is.na(like_prop) | is.na(like_curr)){print(b_prop, beta_0, beta_curr)}
  if((like_prop - like_curr) > log(runif(1))){
    beta_0 = b_prop
    beta_0_accept = 1
  }
  
  # all other betas
  for(i in 1:p_beta){
    for(k in 1:K){
      
      beta_prop = beta_curr
      
      b_prop = beta_prop[k,i] = rnorm(1, beta_curr[k,i], sig_prop_beta[k,i])
      b_curr = beta_curr[k,i]
      
      like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta_curr[k,] + Z[[k]] %*% phi[k,] + omega[k]), log = T)) + dnorm(b_curr, 0, beta_prior_var, log = T)
      like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta_prop[k,] + Z[[k]] %*% phi[k,] + omega[k]), log = T)) + dnorm(b_prop, 0, beta_prior_var, log = T)
      
      # if(is.na(like_prop) | is.na(like_curr)){print(b_prop, beta_0, beta_curr)}
      if((like_prop - like_curr) > log(runif(1))){
        beta_curr[k,i] = b_prop
        beta_accept[k,i] = 1
      }
      
    }
    
  }
  
  pars$beta_0 = beta_0
  pars$beta_0_accept = beta_0_accept
  pars$beta = beta_curr
  pars$beta_accept = beta_accept
  
  return(pars)
  
  
}

update_phi <- function(pars){
  
  # data
  Y = pars$Y
  X = pars$X
  Z = pars$Z
  effort = pars$effort
  omega = pars$omega
  beta = pars$beta
  beta_0 = pars$beta_0
  
  # parameters
  n = pars$n
  p_beta = pars$p_beta
  p_phi = pars$p_phi
  K = pars$K
  
  # beta monitor values
  phi_accept = array(0, dim = c(K, p_phi))
  phi_curr = pars$phi
  sig_prop_phi = pars$sig_prop_phi
  phi_prior_var = pars$phi_prior_var
  
  
  # update_phi
  for(i in 1:p_phi){
    for(k in 1:K){
      
      phi_prop = phi_curr
      
      p_prop = phi_prop[k,i] = rnorm(1, phi_curr[k,i], sig_prop_phi[k,i])
      p_curr = phi_curr[k,i]
      
      like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta[k,] + Z[[k]] %*% phi_curr[k,] + omega[k]), log = T)) + dnorm(p_curr, 0, phi_prior_var, log = T)
      like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta[k,] + Z[[k]] %*% phi_prop[k,] + omega[k]), log = T)) + dnorm(p_prop, 0, phi_prior_var, log = T)
      
      # if(is.na(like_prop) | is.na(like_curr)){print(b_prop, beta_0, beta_curr)}
      if((like_prop - like_curr) > log(runif(1))){
        phi_curr[k,i] = p_prop
        phi_accept[k,i] = 1
      }
      
    }
    
  }
  
  pars$phi = phi_curr
  pars$phi_accept = phi_accept
  
  return(pars)
  
  
}

ll_calc <- function(Y, effort, beta_0, X, beta, Z, phi, omega, K){
  ll = 0
  for(k in 1:K){
    ll = ll + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0 + X[[k]] %*% beta[k,] + Z[[k]] %*% phi[k,] + omega[k]), log = T))
  }
  return(ll)
}

update_omega <- function(pars){
  
  # data
  Y = pars$Y
  X = pars$X
  Z = pars$Z
  effort = pars$effort
  beta_0 = pars$beta_0
  beta = pars$beta
  phi = pars$phi
  omega = pars$omega
  Sigma_species = pars$Sigma_species
  
  # parameters
  n = pars$n
  p = pars$p
  K = pars$K
  
  # choose ellipse
  nu = rmvnorm(1, rep(0, K), Sigma_species)
  nu = nu - mean(nu)
  
  # ll threshold
  logy = ll_calc(Y, effort, beta_0, X, beta, Z, phi, omega, K) + log(runif(1))
  
  # draw initial proposal
  theta = runif(1, 0, 2*pi)
  theta_min = theta - 2*pi
  theta_max = theta
  
  f_proposal = omega * cos(theta) + nu * sin(theta)
  logf = ll_calc(Y, effort, beta_0, X, beta, Z, phi, f_proposal, K)
  
  if(logf > logy){
    keeper = f_proposal
  }else{
    
    while(logf < logy){
      if(theta < 0){
        theta_min = theta
      }else{
        theta_max = theta
      }
      theta = runif(1, theta_min, theta_max)
      f_proposal = omega * cos(theta) + nu * sin(theta)
      logf = ll_calc(Y, effort, beta_0, X, beta, Z, phi, f_proposal, K)
    }
    
    keeper = f_proposal
    
  }
  
  
  pars$omega = keeper
  
  return(pars)
  
  
}

update_sigma_species <- function(pars){
  
  # data
  Sigma_species = pars$Sigma_species
  omega = pars$omega
  
  nu_species = pars$nu_species
  Psi_species = pars$Psi_species
  
  # parameters
  K = pars$K
  
  nu_hat = nu_species + K
  psi_hat = Psi_species + t(omega) %*% omega
  
  pars$Sigma_species = MCMCpack::riwish(nu_hat, psi_hat)
  
  return(pars)
  
}

update_proposal_var_beta <- function(pars, beta_accept_post, i, check_num){
  
  sig_prop = pars$sig_prop_beta
  
  bp = beta_accept_post[,,(i-check_num+1):i]
  accept_rate = apply(bp, c(1,2), mean)
  
  sig_prop = ifelse(accept_rate < 0.2, sig_prop*0.9, sig_prop)
  sig_prop = ifelse(accept_rate > 0.45, sig_prop/0.9, sig_prop)
  
  pars$sig_prop_beta = sig_prop
  
  return(pars)
  
  
}

update_proposal_var_phi <- function(pars, phi_accept_post, i, check_num){
  
  sig_prop = pars$sig_prop_phi
  
  bp = phi_accept_post[,,(i-check_num+1):i]
  accept_rate = apply(bp, c(1,2), mean)
  
  sig_prop = ifelse(accept_rate < 0.2, sig_prop*0.9, sig_prop)
  sig_prop = ifelse(accept_rate > 0.45, sig_prop/0.9, sig_prop)
  
  pars$sig_prop_phi = sig_prop
  
  return(pars)
  
  
}

update_proposal_var_beta_0 <- function(pars, beta_0_accept_post, i, check_num){
  
  sig_prop = pars$sig_prop_beta_0
  
  accept_rate = mean(beta_0_accept_post[(i-check_num+1):i])
  
  if(accept_rate < 0.2){
    sig_prop = sig_prop*0.9
  }else if(accept_rate > 0.45){
    sig_prop = sig_prop/0.9
  }else{
    sig_prop = sig_prop
  }
  
  pars$sig_prop_beta_0 = sig_prop
  
  return(pars)
  
  
}

sampler <- function(fish_dat, mean_covs, mean_covs_log, mean_covs_logit, catch_covs, nits, burnin, thin, check_num = 200, pars = NULL){
  
  if(is.null(pars)){
    pars <- create_pars(fish_dat, mean_covs, mean_covs_log, mean_covs_logit, catch_covs)
  }
  
  p_beta = pars$p_beta
  p_phi = pars$p_phi
  K = pars$K
  
  keep_samp = nits-burnin
  keep_num = length(seq(burnin, nits, by = thin))
  
  beta_post = array(NA, dim = c(K, p_beta, keep_num))
  beta_accept_post = array(NA, dim = c(K, p_beta, nits))
  phi_post = array(NA, dim = c(K, p_phi, keep_num))
  phi_accept_post = array(NA, dim = c(K, p_phi, nits))
  beta_0_post = array(NA, dim = c(keep_num))
  beta_0_accept_post = array(NA, dim = c(nits))
  omega_post = array(NA, dim = c(length(pars$omega), keep_num))
  sigma_species_post = array(NA, dim = c(dim(pars$Sigma_species), keep_num))
  
  pb <- progress_bar$new(
    format = "Burnin Running [:bar] :percent eta: :eta",
    total = (burnin - 1), clear = FALSE, width = 60)
  
  for(i in seq(1, (burnin - 1))){
    
    pars <- update_beta(pars)
    pars <- update_phi(pars)
    pars <- update_omega(pars)
    pars <- update_sigma_species(pars)
    
    beta_accept_post[,,i] = pars$beta_accept
    phi_accept_post[,,i] = pars$phi_accept
    beta_0_accept_post[i] = pars$beta_0_accept
    
    if(i %in% seq(0, burnin-1, by = check_num)){
      pars <- update_proposal_var_beta(pars, beta_accept_post, i, check_num)
      pars <- update_proposal_var_phi(pars, phi_accept_post, i, check_num)
      pars <- update_proposal_var_beta_0(pars, beta_0_accept_post, i, check_num)
    }
    
    pb$tick()
    
  }
  
  pb <- progress_bar$new(
    format = "  Running [:bar] :percent eta: :eta",
    total = (keep_samp + 1), clear = FALSE, width = 60)
  
  j=1
  for(i in seq(burnin, nits)){
    
    pars <- update_beta(pars)
    pars <- update_phi(pars)
    pars <- update_omega(pars)
    pars <- update_sigma_species(pars)
    
    beta_accept_post[,,i] = pars$beta_accept
    phi_accept_post[,,i] = pars$phi_accept
    beta_0_accept_post[i] = pars$beta_0_accept
    
    if(i %in% seq(burnin, nits, by = thin)){
      
      beta_post[,,j] = pars$beta
      phi_post[,,j] = pars$phi
      beta_0_post[j] = pars$beta_0
      omega_post[,j] = pars$omega
      sigma_species_post[,,j] = pars$Sigma_species
      
      j = j+1
      
    }
    
    if(i %in% seq(burnin, nits, by = check_num)){
      pars <- update_proposal_var_beta(pars, beta_accept_post, i, check_num)
      pars <- update_proposal_var_phi(pars, phi_accept_post, i, check_num)
      pars <- update_proposal_var_beta_0(pars, beta_0_accept_post, i, check_num)
    }
    
    
    pb$tick()
    
  }
  
  
  return(list(beta = beta_post,
              beta_accept = beta_accept_post,
              phi = phi_post,
              phi_accept = phi_accept_post,
              sig_prop_beta = pars$sig_prop_beta,
              sig_prop_phi = pars$sig_prop_phi,
              beta_0 = beta_0_post,
              beta_0_accept = beta_0_accept_post,
              sig_prop_beta_0 = pars$sig_prop_beta_0,
              omega = omega_post,
              sig_prop_omega = pars$sig_prop_omega,
              sigma_species = sigma_species_post,
              pars = pars))
  
}

