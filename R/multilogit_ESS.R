multilogit_ESS <- function(Y, X, n_sample = 1000, n_burn = 200, prior_mean = NULL,
                           prior_var = NULL, reference_cat = NULL, probs = TRUE, progress = TRUE){
  # Here, we consider Y to be a matrix of counts with number of rows equal to the number of subjects, and the number of columns equal to the number of categories. 
  # In other words, Y's dimensions are n_sub x n_cat
  Y = as.matrix(Y)
  X = as.matrix(X)
  
  n_sub = nrow(Y)
  n_cat = ncol(Y)
  # X is a n_sub x n_pred matrix of covariates, aka the design matrix. 
  n_pred = ncol(X)
  
  if(nrow(X) != n_sub){
    stop("unequal number of rows in X and Y")
  }
  
  if(!(is.null(prior_mean))){
    
    
    
    if(length(prior_mean) != n_pred){
      stop("prior_mean should have length equal to the number of categories.")
    }
    
  }  
  
  if(!is.null(prior_var)){
    
    
    
    if(!all.equal(dim(prior_var),c(n_pred, n_pred))){
      stop("prior_var should be a square matrix with each dimension equal to the number of categories.")
    }
  }
  # if(!(prior == "flat" | prior == "normal")){
  #   stop("prior must be either flat or normal.")
  #   
  # }
  
  # if(prior == "flat" & ((!is.null(prior_mean)) | (!is.null(prior_var))   )){
  #   stop("If you wish to use a normal prior, use prior = 'normal'")
  # }
  
  if(!is.null(reference_cat)){
    
    if(reference_cat < 1 | reference_cat > n_cat){
      stop("Categories are indexed starting at 1 to the total number of categories.")
    }
    
  }
  # 
  # if(step_size <= 0){
  #   stop("step_size is the standard deviation of the normal random walk proposal distribution. Therefore, it must be a positive constant.")
  #   
  # }
  
  
  out <- multilogit_C_ESS(Y, X, n_sample = n_sample, n_burn = n_burn, 
                      prior = 'normal', prior_mean = prior_mean,
                      prior_var = prior_var, reference_cat = reference_cat, probs = probs, progress = progress)
  
  return(out)
}