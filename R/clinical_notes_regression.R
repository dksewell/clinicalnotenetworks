#' Clinical notes estimation
#' 
#' Estimate the effect of patient features, note importance, and reader 
#' importance on the outcome of unplanned healthcare utilization.  Use 
#' future::plan() to parallelize.
#' 
#' @param data List of lists.  Each element of the outer list corresponds to 
#' running feature_extraction() for a patient and specific time window.  It 
#' should also have appended as an element of the list the outcomes named
#' ed and hospitalization
#' @param outcome_variable character. Which outcome should be modeled?  "ed" or "hospitalization"
#' @param method character. Either "normal_approx" or "adaptive_mcmc".
#' @param mcmc_arguments optional list.  If running adaptive mcmc, supply a list
#' of arguments that correspond to \code{\link[adaptMCMC]{MCMC}}.
#' @param epsilon relative tolerance for convergence assessment
#' @param max_iterations integer.  Maximum number of iterations to run
#' @param phi_boundary_epsilon The optimization of phi will run from 
#' (eps, 1 - eps).
#' 
#' @return Object of class "cnregfit".  Contains a list with the following:
#' \itemize{
#'  \item \code{alpha}: Posterior mode of coefficients for purely patient features
#'  \item \code{beta}: Posterior mode of coefficients for note importance features
#'  \item \code{gamma}: Posterior mode of coefficients for MTS group reading importance features. 
#'  These are bounded between 0 and 1, with higher values indicating more importance.
#'  \item \code{phi}: Posterior mode of penalty for not reading notes.  Bounded between 
#'  0 and 1, with higher values indicating a higher penalty.
#'  \item \code{method}: Which method (normal approximation or MCMC) was used.
#'  \item \code{covariance}: (If \code{method = normal_approx}) Posterior covariance estimate.  Note that for the 
#'  bounded variables, this covariance is on the logit scale.
#'  \item \code{posterior_draws}: (if \code{method = "adaptive_mcmc"}) Object returned from 
#'  \code{\link[adaptMCMC]{MCMC}}.  Retrieve samples via \code{object$posterior_draws$samples}. 
#'  Also, you can feed \code{posterior_draws} into \code{\link[adaptMCMC]{MCMC.add.samples}} 
#'  to obtain more MCMC samples by continuning the MCMC sampler.
#'  \item \code{error}
#' }
#' 
#' @import numDeriv
#' @import future
#' @import future.apply
#' @import holiglm
#' @import adaptMCMC
#' 
#' @export
#' @exportClass cnregfit



clinical_notes_regression = function(data,
                                     outcome_variable = c("ed","hospitalization")[1],
                                     method = c("normal_approx",
                                                "adaptive_mcmc")[1],
                                     mcmc_arguments,
                                     epsilon = 1e-5,
                                     max_iteration = 100,
                                     phi_boundary_epsilon = 1e-2){
  
  warning("Every data frame with factor variables must have levels explicitly specified!")
  
  #--------------------------------------
  # Get baseline features
  #--------------------------------------
  
  N = length(data)
  W = 
    data[[1]]$patient_features %>% 
    mutate(time_in_days = as.integer(end_time - dx_date),
           .keep = "unused")
  for(i in 2:N){
    W %<>%
      bind_rows(
        data[[i]]$patient_features %>% 
          mutate(time_in_days = as.integer(end_time - dx_date),
                 .keep = "unused")
      )
  }
  
  W = model.matrix(~ .,
                   data = W)
  p = ncol(W)
  
  
  #--------------------------------------
  # Get connectivity measures
  #--------------------------------------
  
  # NOTE: We want the covariate matrices for 
  #   \sum_e Z_{ie}(R_{ie} - \phi \ones)'
  #   = Z_i'R_i -\phi \ones\ones'
  
  # Create data matrix and parameter vector for beta | gamma
  ## Create helper for the Z_i'R_i and Z_i'\ones\ones'
  #   (both are qxr matrices).
  helper = function(X){
    Z = 
      model.matrix(~ .,
                   data = X$note_importance)
    
    list(Z_R = crossprod(Z,
                         X$who_read_it), # who_read_it is already multiplied by 2
         Z_ones = crossprod(Z,
                            matrix(1.0,
                                   ncol(Z),
                                   ncol(X$who_read_it)) )
    )
  }
  helper_term1 = function(X){
    Z = 
      model.matrix(~ .,
                   data = X$note_importance)
    crossprod(Z,
              X$who_read_it) # who_read_it is already multiplied by 2
  }
  helper_term2 = function(X){
    Z = 
      model.matrix(~ .,
                   data = X$note_importance)
    
    crossprod(Z,
              matrix(1.0,
                     ncol(Z),
                     ncol(X$who_read_it)) )
  }
  
  Z_R = 
    future_lapply(data, 
                  helper_term1)
  q = nrow(Z_R[[1]])
  r = ncol(Z_R[[1]])
  Z_R %<>%
    array(c(q,r,N))
  
  Z_unscaled =
    future_lapply(data, 
                  helper_term2) %>% 
    unlist() %>% 
    array(c(q,r,N))
  
  
  #--------------------------------------
  # Extract outcome vector
  #--------------------------------------
  y = 
    future_sapply(data,
                  function(X) X[[outcome_variable]])
  
  
  #--------------------------------------
  # Set up objects for optimization
  #--------------------------------------
  
  ## Parameters
  alpha = numeric(p)
  beta = numeric(q)
  gamma = numeric(r) + 1/r
  phi = 0.1
  
  # Tracking for convergence
  old_coefs = new_coefs = numeric(p + q + r + 1) # To assess convergence
  
  ## Create data.frame for hglm
  hglm_df = 
    cbind(W[,-1],
          matrix(0.0,N,r),
          dimnames = list(NULL,
                          paste("constrained",1:r,sep="_"))) %>% 
    as.data.frame() %>% 
    mutate(y = y)
  
  ## Create named constraint vectors
  lower_bounds = numeric(r)
  upper_bounds = rep(1.0, r)
  names(lower_bounds) = 
    names(upper_bounds) =
    paste("constrained",1:r,sep="_")
  
  
  
  #--------------------------------------
  # Perform optimization
  #--------------------------------------
  
  for(iter in 1:max_iteration){
    # Update alpha and beta
    
    ## Compute (Z_R - phi Z\ones \ones')\gamma
    mm = 
      future_apply(Z_R - phi * Z_unscaled,
                   3,
                   function(x) x %*% gamma) %>% 
      t() # Should give Nxq matrix
    
    ## Perform optimization
    fit_given_gamma_phi = 
      glm.fit(x = cbind(W,mm),
              y = y,
              family = "poisson")
    alpha = coef(fit_given_gamma_phi)[1:p]
    beta = coef(fit_given_gamma_phi)[p + 1:q]
    
    
    # Update gamma (and alpha)
    
    ## Compute \beta'(Z_r - \phi Z\ones \ones')
    mm = 
      future_apply(Z_R - phi * Z_unscaled,
                   3,
                   function(x) beta %*% x) %>% 
      t() # Should give Nxr matrix.  And yes, you still need to transpose.
    colnames(mm) = paste("constrained",1:r,sep="_")
    
    ## Modify data.frame for hglm
    for(j in 1:r){
      hglm_df[,p - 1 + j] = # Subtract off one since no intercept
        mm[,j]
    }
    
    ## Perform optimization
    fit_given_beta_phi = 
      hglm(y ~ .,
           data = hglm_df,
           family = "poisson",
           constraints =
             list(lower(lower_bounds),
                  upper(upper_bounds)))
    alpha = coef(fit_give_beta_phi)[1:p]
    gamma = coef(fit_give_beta_phi)[p + 1:r]
    
    
    # Update phi
    llik_given_alpha_beta_gamma = function(x){
      log_lambda =
        W %*% alpha + 
        future_apply(Z_R - x * Z_unscaled,
                     3,
                     function(z) beta %*% z %*% gamma )
      dpois(y,
            exp(log_lambda),
            log = TRUE) %>% 
        sum()
    }
    
    phi = 
      optimize(llik_given_alpha_beta_gamma,
               c(phi_boundary_epsilon, 1.0 - phi_boundary_epsilon),
               maximum = TRUE)$maximum
    
    
    # Assess convergence
    if( max( abs( (old_coefs - new_coefs) / old_coefs) ) < epsilon){
      break
    }else{
      old_coefs = new_coefs
    }
    
  }
  
  
  
  #--------------------------------------
  # Create posterior function 
  #--------------------------------------
  
  # Create helper functions
  logit = function(x) log(x / (1.0 - x))
  expit = function(x) 1.0 / (1.0 + exp(-x))
  
  # Create log-likelihood function
  lpost = function(x){
    alpha = x[1:p]
    beta = x[p + 1:q]
    transformed_gamma = x[p + q + 1:r]
    gamma = expit(transformed_gamma)
    transformed_phi = x[p + q + r + 1]
    phi = expit(transformed_phi)
    
    log_lambda =
      W %*% alpha + 
      future_apply(Z_R - phi * Z_unscaled,
                   3,
                   function(z) beta %*% z %*% gamma )
    
    # Don't forget Jacobian for the prior 
    #   since we are looking at the Hessian wrt logit(theta) and not theta
    # \pi(\phi) = 1
    # Rightarrow \pi(logit(\phi)) = | d \phi / d logit(\phi) |
    #                             = e^\phi / (1 + e^\phi)^2
    sum(dpois(y,
              exp(log_lambda),
              log = TRUE)) +
      sum(transformed_gamma) - 
      2.0 * sum(log(1.0 + exp(transformed_gamma))) +
      transformed_phi - 
      2.0 * log(1.0 + exp(transformed_phi))
  }
  
  
  
  
  
  if(method == "normal_approx"){
    
    #--------------------------------------
    # If normal approximation, find covariance matrix
    #--------------------------------------
    
    cov_estimate = NULL
    try({
      H = 
        numDeriv::hessian(lpost,
                          c(alpha,
                            beta,
                            logit(gamma),
                            logit(phi)))
    },silent = T)
    try({
      cov_estimate = 
        chol2inv(chol(-H))
    },silent = T)
    if(is.null(cov_estimate)){
      try({
        cov_estimate = 
          qr.solve(-H)
      },silent = T)
    }
    if(is.null(cov_estimate)){
      try({
        cov_estimate = 
          solve(-H)
      },silent = T)
    }
    
    
    
    #--------------------------------------
    # Create object to return
    #--------------------------------------
    
    if(is.null(cov_estimate)){
      warning("Could not invert negative Hessian to get covariance matrix. Returning posterior mode only.")
      
      object =
        list(alpha = alpha,
             beta = beta,
             gamma = gamma,
             phi = phi,
             method = "normal_approx",
             error = "Error: singular_hessian")
      
    }else{
      
      object = 
        list(alpha = alpha,
             beta = beta,
             gamma = gamma,
             phi = phi,
             covariance = cov_estimate,
             method = "normal_approx")
      
    }
  
  }
  
  
  if(method == "adaptive_mcmc"){
    
    #--------------------------------------
    # Set arguments for adaptMCMC::MCMC if missing
    #--------------------------------------
    
    if(missing(mcmc_arguments)) mcmc_arguments = list()
    mcmc_arguments$p = lpost
    mcmc_arguments$init = 
      c(alpha,
        beta,
        logit(gamma),
        logit(phi))
    
    if(is.null(mcmc_arguments$n)) mcmc_arguments$n = 5e3
    if(is.null(mcmc_arguments$adapt)) mcmc_arguments$adapt = TRUE
    if(is.null(mcmc_arguments$acc.rate)) mcmc_arguments$acc.rate = 0.234
    
    #--------------------------------------
    # Perform MCMC sampling
    #--------------------------------------
    
    post_samples = NULL
    try({
      post_samples = 
        do.call(adaptMCMC::MCMC,
                mcmc_arguments)
    }, silent = TRUE)
    
    
    #--------------------------------------
    # Create object to return
    #--------------------------------------
    
    if(is.null(post_samples)){
      
      object = 
        list(alpha = alpha,
             beta = beta,
             gamma = gamma,
             phi = phi,
             method = "adaptive_mcmc",
             error = "Error: Error in call to adaptMCMC::MCMC()")
      
    }else{
      
      object = 
        list(alpha = alpha,
             beta = beta,
             gamma = gamma,
             phi = phi,
             posterior_draws = post_samples,
             method = "adaptive_mcmc")
      
    }
    
  }
  
  
  
  class(object) = "cnregfit"
  return(object)
}