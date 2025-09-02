#' Clinical notes estimation
#' 
#' Estimate the effect of patient features, note importance, and reader 
#' importance on the outcome of unplanned healthcare utilization.
#' 
#' @param data List of lists.  Each element of the outer list corresponds to 
#' running feature_extraction() for a patient and specific time window.  It 
#' should also have appended as an element of the list the outcomes named
#' ed and hospitalization
#' @param epsilon relative tolerance for convergence assessment
#' @param max_iterations integer.  Maximum number of iterations to run
#' @param phi_boundary_epsilon The optimization of phi will run from 
#' (eps, 1 - eps).
#' 
#' 
#' @import holiglm
#' @export



clinical_notes_regression = function(data,
                                     outcome_variable = c("ed","hospitalization")[1],
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
    lapply(data, 
           helper_term1)
  q = nrow(Z_R[[1]])
  r = ncol(Z_R[[1]])
  Z_R %<>%
    array(c(q,r,N))
  
  Z_unscaled =
    lapply(data, 
           helper_term2) %>% 
    unlist() %>% 
    array(c(q,r,N))
  
  
  #--------------------------------------
  # Extract outcome vector
  #--------------------------------------
  y = 
    sapply(data,
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
      apply(Z_R - phi * Z_unscaled,
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
      apply(Z_R - phi * Z_unscaled,
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
        apply(Z_R - x * Z_unscaled,
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
  
  
  
  
  
  
}