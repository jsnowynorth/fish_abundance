
create_pars <- function(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs){
  
  K = nlevels(fish_dat$COMMON_NAME)
  levs = levels(fish_dat$COMMON_NAME)
  lake_index = (fish_dat %>% filter(COMMON_NAME == levs[1]))$DOW # lake_index == lake_id[2]
  lake_id = levels(fish_dat$DOW)
  n_lakes = length(levels(fish_dat$DOW))
  n_obs = (fish_dat %>% filter(COMMON_NAME == levs[1]) %>% select(TOTAL_CATCH))$TOTAL_CATCH
  
  pars = list()
  
  # data
  pars$Y = list()
  pars$X = list()
  pars$Z = list()
  pars$effort = list()
  
  
  # projection matrix
  
  # only fixed covs
  P = fish_dat %>% 
    distinct(DOW, .keep_all = T) %>%
    select(all_of(mean_covs), DOW) %>% 
    select(-all_of(temporal_covs), DOW) %>% 
    left_join(fish_dat %>% 
                select(DOW, all_of(temporal_covs)) %>% 
                group_by(DOW) %>% 
                summarise_all(mean), by = 'DOW') %>% 
    select(-DOW) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
    mutate(secchi = secchi - mean(secchi)) %>% 
    as.matrix()
  
  pars$P = diag(nrow(P)) - P %*% solve(t(P) %*% P) %*% t(P)
  
  # ones = matrix(1, nrow = 6, ncol = 6)
  # kronecker(pars$P[1:12,1:12], ones)
  # 
  for(k in 1:K){
    pars$Y[[k]] = (fish_dat %>% filter(COMMON_NAME == levs[k]) %>% select(TOTAL_CATCH))$TOTAL_CATCH
  }
  
  for(k in 1:K){
    
    # X = fish_dat %>% 
    #   filter(COMMON_NAME == levs[k]) %>% 
    #   select(all_of(mean_covs)) %>% 
    #   mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
    #   mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
    #   mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
    #   mutate(Int = 1)
    
    X = fish_dat %>% 
      filter(COMMON_NAME == levs[k]) %>% 
      select(all_of(mean_covs)) %>% 
      mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
      mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
      mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
      mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
      mutate(secchi = secchi - mean(secchi))
    
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
  
  # mle starts
  X = fish_dat %>% 
    filter(COMMON_NAME == 'bluegill') %>% 
    select(all_of(mean_covs)) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ ./100) %>% 
    mutate_at(vars(all_of(mean_covs_logit)), ~ car::logit(., adjust = 0.001)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~ log(.)) %>% 
    mutate_at(vars(all_of(mean_covs_log)), ~. - mean(.)) %>% 
    mutate(secchi = secchi - mean(secchi)) %>% 
    as.matrix()
  
  Y = fish_dat %>%
    select(CPUE, COMMON_NAME, DOW, SURVEYDATE, TN) %>%
    pivot_wider(names_from = 'COMMON_NAME', values_from = 'CPUE') %>% 
    select(-c(DOW, SURVEYDATE, TN)) %>% 
    as.matrix()
  
  pars$beta = t(solve(t(X) %*% X) %*% t(X) %*% log(Y + 0.00001))
  pars$beta_0 = log(apply(Y, 2, mean))
  
  # pars$beta = array(0, dim = c(K, pars$p_beta))
  pars$beta_accept =  array(0, dim = c(K, pars$p_beta))
  # pars$beta_0 = rep(0, K)
  pars$beta_0_accept =  rep(0, K)
  pars$beta_prior_var = 100
  pars$phi = array(0, dim = c(K, pars$p_phi))
  pars$phi_accept =  array(0, dim = c(K, pars$p_phi))
  pars$phi_prior_var = 100
  
  pars$gamma = array(0, dim = c(length(lake_index),K))
  pars$theta = array(0, dim = c(length(lake_index),K))
  
  # hyperpriors
  pars$Sigma_species = diag(K)
  
  # half t priors
  pars$a = rep(1, K)
  pars$A = 1e5
  pars$nu_species = 2
  
  # pars$omega = matrix(rep(0, K), ncol = K)
  pars$omega = matrix(0, nrow = n_lakes, ncol = K)
  pars$omega_star = matrix(0, nrow = n_lakes, ncol = K)
  
  # spatial parameters
  spat_dat = fish_dat %>% 
    distinct(DOW, .keep_all = T) %>% 
    select(DOW, LAKE_CENTER_UTM_EASTING, LAKE_CENTER_UTM_NORTHING)
  
  d = rdist(cbind(spat_dat$LAKE_CENTER_UTM_EASTING, spat_dat$LAKE_CENTER_UTM_NORTHING))/1000
  phi = 10
  pars$Sigma_spatial = Matrix(exp(-d/phi))
  pars$Sigma_spatial_inv = solve(pars$Sigma_spatial)
  # pars$up_chol_spatial = spam::chol(as.spam(pars$Sigma_spatial))
  # pars$up_chol_spatial = t(Matrix::chol(pars$Sigma_spatial))
  pars$up_chol_spatial = t(chol(exp(-d/phi)))
  
  # spat_dat %>%
  #   mutate(lake = pars$Sigma_spatial[1000,]) %>%
  #   ggplot(., aes(x = LAKE_CENTER_UTM_EASTING, y = LAKE_CENTER_UTM_NORTHING, color = lake)) +
  #   geom_point() +
  #   scale_color_gradient(low = 'yellow', high = 'red')
  
  
  # Proposal variances
  pars$sig_prop_beta = array(0.1, dim = c(K, pars$p_beta))
  pars$sig_prop_phi = array(0.1, dim = c(K, pars$p_phi))
  pars$sig_prop_beta_0 = rep(0.1, K)
  
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
  # Z = pars$Z
  effort = pars$effort
  omega = pars$omega_star
  # phi = pars$phi
  theta = pars$theta
  gamma = pars$gamma
  
  # parameters
  n = pars$n
  p_beta = pars$p_beta
  p_phi = pars$p_phi
  K = pars$K
  
  # beta monitor values
  beta_curr = pars$beta
  beta_accept = array(0, dim = c(K, p_beta))
  sig_prop_beta = pars$sig_prop_beta
  
  beta_0_curr = pars$beta_0
  beta_0_accept = rep(0, K)
  sig_prop_beta_0 = pars$sig_prop_beta_0
  
  beta_prior_var = pars$beta_prior_var
  
  
  # construct omega
  ind_array = tibble(id = pars$lake_id, as_tibble(omega, .name_repair = ~LETTERS[1:K]))
  lake_array = tibble(id = pars$lake_index)
  OMEGA = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))
  
  for(k in 1:K){
    
    beta_prop = beta_0_curr
    
    b_prop = beta_prop[k] = rnorm(1, beta_0_curr[k], sig_prop_beta_0[k])
    b_curr = beta_0_curr[k]
    
    like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(b_curr + X[[k]] %*% beta_curr[k,] + theta[,k] + OMEGA[,k]), log = T)) + dnorm(b_curr, 0, beta_prior_var, log = T)
    like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(b_prop + X[[k]] %*% beta_curr[k,] + theta[,k] + OMEGA[,k]), log = T)) + dnorm(b_prop, 0, beta_prior_var, log = T)
    
    if((like_prop - like_curr) > log(runif(1))){
      beta_0_curr[k] = b_prop
      beta_0_accept[k] = 1
    }
    
  }
  
  beta_0 = beta_0_curr
  pars$beta_0 = beta_0
  pars$beta_0_accept = beta_0_accept
  
  # all other betas
  for(k in 1:K){
    
    beta_prop = beta_curr
    
    b_prop = beta_prop[k,] = c(rmvnorm(1, beta_curr[k,], diag(sig_prop_beta[k,])))
    b_curr = beta_curr[k,]
    
    like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0[k] + X[[k]] %*% b_curr + theta[,k] + OMEGA[,k]), log = T)) + dmvnorm(b_curr, rep(0, p_beta), diag(rep(beta_prior_var, p_beta)), log = T)
    like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0[k] + X[[k]] %*% b_prop + theta[,k] + OMEGA[,k]), log = T)) + dmvnorm(b_prop, rep(0, p_beta), diag(rep(beta_prior_var, p_beta)), log = T)
    
    if((like_prop - like_curr) > log(runif(1))){
      beta_curr[k,] = b_prop
      beta_accept[k,] = rep(1, p_beta)
    }
    
  }

  
  
  # for(i in 1:p_beta){
  #   for(k in 1:K){
  #     
  #     beta_prop = beta_curr
  #     
  #     b_prop = beta_prop[k,i] = rnorm(1, beta_curr[k,i], sig_prop_beta[k,i])
  #     b_curr = beta_curr[k,i]
  #     
  #     like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta_curr[k,] + theta[,k] + OMEGA[,k]), log = T)) + dnorm(b_curr, 0, beta_prior_var, log = T)
  #     like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(X[[k]] %*% beta_prop[k,] + theta[,k] + OMEGA[,k]), log = T)) + dnorm(b_prop, 0, beta_prior_var, log = T)
  #     
  #     if((like_prop - like_curr) > log(runif(1))){
  #       beta_curr[k,i] = b_prop
  #       beta_accept[k,i] = 1
  #     }
  #     
  #   }
  #   
  # }
  
  pars$beta = beta_curr
  pars$beta_accept = beta_accept
  
  for(k in 1:K){
    pars$gamma[,k] = pars$beta_0[k] + X[[k]] %*% pars$beta[k,] + OMEGA[,k]
  }
  
  return(pars)
  
  
}

update_phi <- function(pars){
  
  # data
  Y = pars$Y
  # X = pars$X
  Z = pars$Z
  effort = pars$effort
  # omega = pars$omega_star
  # beta = pars$beta
  theta = pars$theta
  gamma = pars$gamma
  
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
      
      like_curr = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(gamma[,k] + Z[[k]] %*% phi_curr[k,]), log = T)) + dnorm(p_curr, 0, phi_prior_var, log = T)
      like_prop = sum(dpois(Y[[k]], lambda = effort[[k]]*exp(gamma[,k] + Z[[k]] %*% phi_prop[k,]), log = T)) + dnorm(p_prop, 0, phi_prior_var, log = T)
    
      if((like_prop - like_curr) > log(runif(1))){
        phi_curr[k,i] = p_prop
        phi_accept[k,i] = 1
      }
      
    }
    
  }
  
  pars$phi = phi_curr
  pars$phi_accept = phi_accept
  
  for(k in 1:K){
    theta[,k] = Z[[k]] %*% pars$phi[k,]
  }
  pars$theta = theta - mean(theta)
  
  return(pars)
  
  
}

ll_calc <- function(Y, effort, beta_0, X, beta, theta, omega, mu, K, lake_id, lake_index){
  
  ind_array = tibble(id = lake_id, as_tibble(omega, .name_repair = ~LETTERS[1:K]))
  lake_array = tibble(id = lake_index)
  OMEGA = as.matrix(lake_array %>% right_join(ind_array, by = 'id') %>% select(-id))
  
  ll = 0
  for(k in 1:K){
    
    ll = ll + sum(dpois(Y[[k]], lambda = effort[[k]]*exp(beta_0[k] + X[[k]] %*% beta[k,] + theta[,k] + OMEGA[,k]), log = T))
  }
  return(ll)
}

update_omega <- function(pars){
  
  # data
  Y = pars$Y
  X = pars$X
  Z = pars$Z
  # theta = pars$theta
  effort = pars$effort
  mu = pars$mu
  beta = pars$beta
  beta_0 = pars$beta_0
  phi = pars$phi
  omega = pars$omega
  Sigma_species = pars$Sigma_species
  up_chol_spatial = pars$up_chol_spatial
  lake_id = pars$lake_id
  lake_index = pars$lake_index
  n_lakes = pars$n_lakes
  P = pars$P
  
  # parameters
  n = pars$n
  p = pars$p
  K = pars$K
  
  # choose ellipse
  # U = Matrix::kronecker(Matrix(t(chol(Sigma_species))), up_chol_spatial)
  U = fastmatrix::kronecker.prod(t(chol(Sigma_species)), as.matrix(up_chol_spatial))
  b = rnorm(n_lakes*K)
  v = as.matrix(U%*%b)
  nu = v - mean(v)
  
  # ll threshold
  logy = ll_calc(Y, effort, beta_0, X, beta, pars$theta, omega, mu, K, lake_id, lake_index) + log(runif(1))
  
  # draw initial proposal
  theta = runif(1, 0, 2*pi)
  theta_min = theta - 2*pi
  theta_max = theta
  
  f_proposal = matrix(c(omega) * cos(theta) + nu * sin(theta), ncol = K)
  logf = ll_calc(Y, effort, beta_0, X, beta, pars$theta, f_proposal, mu, K, lake_id, lake_index)
  
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
      f_proposal = matrix(c(omega) * cos(theta) + nu * sin(theta), ncol = K)
      logf = ll_calc(Y, effort, beta_0, X, beta, pars$theta, f_proposal, mu, K, lake_id, lake_index)
    }
    
    keeper = f_proposal
    
  }
  
  pars$omega = keeper
  pars$omega_star = P %*% pars$omega
  
  return(pars)
  
  
}

update_sigma_species <- function(pars){
  
  # data
  Sigma_species = pars$Sigma_species
  omega = pars$omega
  mu = pars$mu
  Sigma_spatial_inv = pars$Sigma_spatial_inv
  
  nu_species = pars$nu_species
  a = pars$a
  A = pars$A
  
  # parameters
  K = pars$K
  n_lakes = pars$n_lakes
  
  # update sigma species
  nu_hat = n_lakes + nu_species + K - 1
  psi_hat = as.matrix(t(omega) %*% Sigma_spatial_inv %*% (omega) + 2*nu_species*diag(a))
  
  pars$Sigma_species = MCMCpack::riwish(nu_hat, psi_hat)
  
  
  # update a
  I_sig = solve(pars$Sigma_species)
  for(k in 1:K){
    
    a_hat = (nu_species + n_lakes)/2
    b_hat = nu_species * I_sig[k,k] + 1/A^2
    
    pars$a[k] = 1/rgamma(1, a_hat, b_hat)
    
  }
  
  return(pars)
  
}

# https://m-clark.github.io/docs/ld_mcmc/index_onepage.html
update_proposal_var_beta_0 <- function(pars, beta_0_accept_post, i, check_num){
  
  sig_prop = pars$sig_prop_beta_0
  
  bp = beta_0_accept_post[,(i-check_num+1):i]
  accept_rate = apply(bp, 1, mean)
  
  sig_prop = ifelse(accept_rate < 0.25, sig_prop*0.9, sig_prop)
  sig_prop = ifelse(accept_rate > 0.40, sig_prop/0.9, sig_prop)
  
  pars$sig_prop_beta_0 = sig_prop
  
  return(pars)
  
  
}

update_proposal_var_beta <- function(pars, beta_accept_post, i, check_num){
  
  sig_prop = pars$sig_prop_beta
  
  bp = beta_accept_post[,,(i-check_num+1):i]
  accept_rate = apply(bp, c(1,2), mean)
  
  sig_prop = ifelse(accept_rate < 0.15, sig_prop*0.9, sig_prop)
  sig_prop = ifelse(accept_rate > 0.35, sig_prop/0.9, sig_prop)
  
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

sampler <- function(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs, nits, burnin, thin, check_num = 200, pars = NULL){
  
  if(is.null(pars)){
    pars <- create_pars(fish_dat, mean_covs, temporal_covs, mean_covs_log, mean_covs_logit, catch_covs)
  }
  
  
  p_beta = pars$p_beta
  p_phi = pars$p_phi
  K = pars$K
  
  keep_samp = nits-burnin
  keep_num = length(seq(burnin, nits, by = thin))
  
  beta_0_post = array(NA, dim = c(K, keep_num))
  beta_0_accept_post = array(NA, dim = c(K, nits))
  beta_post = array(NA, dim = c(K, p_beta, keep_num))
  beta_accept_post = array(NA, dim = c(K, p_beta, nits))
  phi_post = array(NA, dim = c(K, p_phi, keep_num))
  phi_accept_post = array(NA, dim = c(K, p_phi, nits))
  omega_post = array(NA, dim = c(dim(pars$omega), keep_num))
  sigma_species_post = array(NA, dim = c(dim(pars$Sigma_species), keep_num))
  
  
  pb <- progress_bar$new(
    format = "Burnin Running [:bar] :percent eta: :eta",
    total = (burnin - 1), clear = FALSE, width = 60)
  
  for(i in seq(1, (burnin - 1))){
    
    pars <- update_beta(pars)
    pars <- update_phi(pars)
    pars <- update_omega(pars)
    pars <- update_sigma_species(pars)
    
    beta_0_accept_post[,i] = pars$beta_0_accept
    beta_accept_post[,,i] = pars$beta_accept
    phi_accept_post[,,i] = pars$phi_accept
    
    if(i %in% seq(0, burnin-1, by = check_num)){
      pars <- update_proposal_var_beta_0(pars, beta_0_accept_post, i, check_num)
      pars <- update_proposal_var_beta(pars, beta_accept_post, i, check_num)
      pars <- update_proposal_var_phi(pars, phi_accept_post, i, check_num)
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
    
    beta_0_accept_post[,i] = pars$beta_0_accept
    beta_accept_post[,,i] = pars$beta_accept
    phi_accept_post[,,i] = pars$phi_accept
    
    if(i %in% seq(burnin, nits, by = thin)){
      
      beta_0_post[,j] = pars$beta_0
      beta_post[,,j] = pars$beta
      phi_post[,,j] = pars$phi
      omega_post[,,j] = pars$omega
      sigma_species_post[,,j] = pars$Sigma_species
      
      j = j+1
      
    }
    
    
    if(i %in% seq(burnin, nits, by = check_num)){
      pars <- update_proposal_var_beta_0(pars, beta_0_accept_post, i, check_num)
      pars <- update_proposal_var_beta(pars, beta_accept_post, i, check_num)
      pars <- update_proposal_var_phi(pars, phi_accept_post, i, check_num)
    }
    
    
    pb$tick()
    
  }
  
  return(list(beta_0 = beta_0_post,
              beta_0_accept = beta_0_accept_post,
              beta = beta_post,
              beta_accept = beta_accept_post,
              phi = phi_post,
              phi_accept = phi_accept_post,
              sig_prop_beta = pars$sig_prop_beta,
              sig_prop_phi = pars$sig_prop_phi,
              omega = omega_post,
              sig_prop_omega = pars$sig_prop_omega,
              sigma_species = sigma_species_post,
              pars = pars))
  
}
