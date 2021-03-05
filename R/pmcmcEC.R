#' Title
#'
#' @param theta 
#' @param prior_params 
#'
#' @return
#' @export
#'
#' @examples
prior1 <- function(theta, prior_params){
  p=0
  for (i in 1:length(theta)) {
    p=p+dnorm(theta[i],prior_params[i,1], prior_params[i,2], log = TRUE)
  }
  return(p)
}

#' Title
#'
#' @param theta 
#' @param prior_params 
#' @param i 
#'
#' @return
#' @export
#'
#' @examples
prior2 <- function(theta, prior_params, i){
  p <- dnorm(theta[i],prior_params[i,1], prior_params[i,2], log = TRUE)
  return(p)
}



#' Title
#'
#' @param theta 
#' @param prior_params 
#' @param i 
#'
#' @return
#' @export
#'
#' @examples
gprior <- function(theta,prior_params,i,rate){
  if(i %in% c(1,6,7)){
    p <- dnorm(theta[i],prior_params[i,1], prior_params[i,2], log = TRUE)
  }
  else{p <- dgamma(theta[i], shape = rate*prior_params[i,1], rate, log = TRUE)}
  return(p)
}

#' Title
#'
#' @param theta 
#' @param prop_params 
#'
#' @return
#' @export
#'
#' @examples
proposal1 <- function(theta,prop_params){
  for (i in 1:length(theta)) {
    if(i %in% c(4,5)){theta[i] <- max(0,rnorm(1,mean = theta[i],prop_params[i]))}
    else{theta[i] <- rnorm(1,mean = theta[i],prop_params[i])}
  }
  return(theta)
}

#' Title
#'
#' @param theta 
#' @param prop_params 
#' @param i 
#'
#' @return
#' @export
#'
#' @examples
proposal2 <- function(theta, prop_params, i){
  if(i %in% c(4,5)){theta[i] <- max(1e-10,rnorm(1,mean = theta[i],prop_params[i]))}
  else{theta[i] <- rnorm(1,mean = theta[i],prop_params[i])}
  return(theta)
}



#' Title
#'
#' @param p 
#' @param pinit 
#' @param m 
#' @param theta 
#' @param gamma 
#'
#' @return
#' @export
#'
#' @examples
particle_est <- function(p,pinit, m,theta, gamma){
  par <- c(length(p), theta[1],  theta[2],  theta[3],theta[4],  theta[5], theta[6],1/7, gamma)
  particleFilterCpp(p,pinit,m,par, 1000, theta[7])
}



#' Title
#'
#' @param p 
#' @param pinit 
#' @param m 
#' @param init 
#' @param n_iter 
#' @param prop_params 
#' @param prior_params 
#' @param gamma 
#'
#' @return
#' @export
#'
#' @examples
pmcmc <- function(p,pinit, m ,init, n_iter, prop_params,prior_params, rate, gamma = 43){
  
  theta <- matrix(NA, nrow = length(init), ncol = n_iter*7+1)
  theta[,1] <- init
  ll <- particle_est(p,pinit,m,theta[,1], gamma)$Loglik
  ac <- rep(0,7)
  for (i in 2:n_iter) {
    print(i)
    # take a step in parameter space
    for (j in 1:7) {
      theta_new <- proposal2(theta[,7*(i-2)+j],prop_params, j)
      u <- runif(1)
      ll_new <- particle_est(p,pinit,m,theta_new, gamma)$Loglik
      alpha <- exp( gprior(theta_new,prior_params,j, rate) - gprior(theta[,6*(i-2)+j],prior_params,j, rate) + ll_new - ll )
      if (u < alpha) {
        theta[,7*(i-2)+j+1] <- theta_new
        ll <- ll_new
        print('accepted')
        print(theta_new)
        ac[j] <- ac[j]+1
      }
      else{theta[,7*(i-2)+j+1] <- theta[,7*(i-2)+j]
      print('rejected')}
    }
  }
  
  out <- list(theta = theta[,-1], acr = ac/n_iter)
  
  return(out)
}
