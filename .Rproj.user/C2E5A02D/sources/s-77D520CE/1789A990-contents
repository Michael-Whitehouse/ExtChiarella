
#' Title
#'
#' @param N length of trajectory
#' @param kappa Fundamentalist Coefficient
#' @param beta Chartist Coefficient
#' @param sigma_n Noise trader volitility
#' @param sigma_v Value signal volitility
#' @param g Growth
#' @param v0 Initial fundamental value
#' @param p0 Initial Price
#' @param m0 Initial momentum
#' @param alpha Smoothing parameter
#' @param gamma Trend signal saturation parameter
#'
#' @return
#' @export
#'
#' @examples 
ExtendedChiarella <- function(N, kappa, beta, sigma_n,sigma_v, g, v0, p0,m0, alpha, gamma ){

  v <- c()
  v[1] <- v0
  for (t in 1:(N-1)) {
    v[t+1] <- v[t]+g+rnorm(1,0,sigma_v)
  }

  p <- c()
  m <- c()

  m[1] <- m0
  p[1] <- p0
  p[2] <- p[1] + kappa*(v[1]-p[1])+beta*tanh(gamma*m[1]) + rnorm(1,0,sigma_n)
  m[2] <- (1-(alpha))*m[1] + alpha*(p[2]-p[1])


  for (t in 2:(N-1)) {
    p[t+1] <- p[t] + kappa*(v[t]-p[t]) + beta*tanh(gamma*m[t]) + rnorm(1,0,sigma_n)
    m[t+1] <- (1-alpha)*m[t]+ alpha*(p[t]-p[t-1])
  }

  t <-list(price = p,value = v ,momentum = m, returns = diff(p))

  return(t)

}


#' Title
#'
#' @param kappa Linear coeff
#' @param kappa_3 Cubic term coeff
#' @param x Argument
#'
#' @return
#' @export
#'
#' @examples f
NL <- function(kappa, kappa_3,x){
  kappa*x + kappa_3*x^3
}

#' Title
#'
#' @param N Length of trajectory
#' @param kappa Linear coeff for fundamentalist demand function
#' @param kappa_3 cubic term coefficient for fundamentalist demand function
#' @param beta Chartist Coefficient
#' @param sigma_n Noise trader volitility
#' @param sigma_v Value signal volitility
#' @param g Growth
#' @param v0 Initial fundamental value
#' @param p0 Initial Price
#' @param m0 Initial momentum
#' @param alpha Smoothing parameter
#' @param gamma Trend signal saturation parameter
#'
#' @return
#' @export
#'
#' @examples f
ExtendedNLChiarella <- function(N, kappa, kappa_3, beta, sigma_n,sigma_v, g, v0, p0,m0, alpha, gamma ){

  v <- c()
  v[1] <- v0
  for (t in 1:(N-1)) {
    v[t+1] <- v[t]+g+rnorm(1,0,sigma_v)
  }

  p <- c()
  m <- c()

  m[1] <- m0
  p[1] <- p0
  p[2] <- p[1] + NL(kappa, kappa_3,v[1]-p[1])+beta*tanh(gamma*m[1]) + rnorm(1,0,sigma_n)
  m[2] <- (1-(alpha))*m[1] + alpha*(p[2]-p[1])
  #+ alpha*(p[1])


  for (t in 2:(N-1)) {
    p[t+1] <- p[t] + NL(kappa, kappa_3,v[t]-p[t]) + beta*tanh(gamma*m[t]) + rnorm(1,0,sigma_n)
    m[t+1] <- (1-alpha)*m[t]+ alpha*(p[t]-p[t-1])
  }

  t <-list(price = p,value = v ,momentum = m)

  return(t)

}

#' Title
#'
#' @param alpha Smoothing parameter
#' @param p price
#' @param t0 smoothing start time
#'
#' @return
#' @export
#'
#' @examples f
msignal <- function(alpha, p, t0){
  m <- c()
  m[1] <- (p[t0]-p[t0-1])
  for (t in 1:(length(p)-t0)) {

    m[t+1] <- (1-alpha)*m[t] + alpha*(p[t0 + t]-p[t0 + t-1])
  }
  return(m)
}

#' Title
#'
#' @param p log price
#' @param t0 smoothing start time
#' @param alpha_2 smoothing parameter
#'
#' @return
#' @export
#'
#' @examples f
volsig <- function(p, t0, alpha_2){
  S <- c()
  S[1] <- var(p[1:(t0-1)])
  rmean <- mean(diff(p[1:t0]))

  for (t in 1:(length(p)-t0)) {
    S[t+1] <- (1-alpha_2)*S[t] + alpha_2*((p[t0+t]-p[t+t0-1]) - rmean)^2
    rmean <- mean(diff(p[1:(t+t0)]))
  }
  return(S)
}

#' Title
#'
#' @param p price signal
#' @param kappa Fundamentalist demand coefficient
#' @param beta Chartist coefficient
#' @param gamma Saturation parameter
#' @param sigma_n Noise trader volitility
#' @param sigma_v Value signal volitility
#' @param g Growth
#' @param v0 Initial fundamental value
#' @param t0 Burn in time
#' @param sigma_0 initial variance
#' @param alpha smoothing parameter
#'
#' @return
#' @export
#'
#' @examples f
Kalmanfs <- function(p, kappa, beta, gamma, sigma_n,sigma_v, sigma_0, g, alpha, v0,t0){

  # initialise

  u <- msignal(alpha, p, t0)


  p_ <- p
  p <- p[t0:length(p)]
  vtilde <- c()
  vtilde_ <- c()
  V <- c()
  V_ <- c()

  vtilde[1] <- v0
  V[1] <- sigma_0^2
  vtilde_[1] <- v0
  V_[1] <- sigma_0^2

  # Predict and Filter

  for (t in 2:length(p)) {
    vtilde[t] <- vtilde_[t-1]+g
    V[t] <- V_[t-1] + sigma_v^2

    K <- kappa*V[t]/(kappa^2*V[t]+sigma_n^2)
    vtilde_[t] <- vtilde[t] + K*(p[t]-p[t-1]-kappa*(vtilde[t]-p[t-1])-beta*u[t])
    V_[t] = V[t] - kappa*K*V[t]
  }

  # Smoother
  v_ <- vtilde_
  vT <- c()
  VT <- c()

  vT[length(p)] <- vtilde_[length(p)]
  VT[length(p)] <- V_[length(p)]

  J <- c()
  J[length(p)] <- 0
  for (t in (length(p)):2) {
    J[t-1] = V_[t-1]/V[t]
    vT[t-1] = v_[t-1]+J[t-1]*(vT[t] - v_[t-1])
    VT[t-1] = V_[t-1] + (J[t-1]^2)*(VT[t]-V[t-1])
  }


  out <- list(value = vT,variance = VT)

  return(out)


}




#' Title
#'
#' @param p desc
#' @param vT desc
#' @param VT desc
#' @param u desc
#'
#' @return
#' @export
#'
#' @examples f
Mstep <- function(p,vT,VT,u){

  a <- sum(V+v^2 - 2*v*p_ + p_*2)
  bc <- sum((v-p_)*u)

  A <- matrix(c(a,b,c,d), nrow = 2)
  b <- c(sum((v-p_)*(p-p_)), sum(u*(p-p_)))

  out <- solve(A,b)

  ls <- list(kappa = out[1], beta = out[2])

  return(ls)
}



#' Title
#'
#' @param N Length
#' @param kappa desc
#' @param kappa_3 desc
#' @param rho desc
#' @param beta desc
#' @param sigma_n v
#' @param sigma_v desc
#' @param g desc
#' @param v0 desc
#' @param p0 desc
#' @param m0 desc
#' @param S0 desc
#' @param voltar desc
#' @param alpha_1 desc
#' @param alpha_2 desc
#' @param gamma desc
#'
#' @return
#' @export
#'
#' @examples f
ExtendedvolChiarella <- function(N, kappa, kappa_3, rho, beta, sigma_n,sigma_v, g, v0, p0, m0, S0, voltar , alpha_1, alpha_2, gamma ){

  v <- c()
  v[1] <- v0
  for (t in 1:(N-1)) {
    v[t+1] <- v[t]+g+rnorm(1,0,sigma_v)
  }


  p <- c()
  m <- c()
  S <- c()

  S[1] <- S0
  m[1] <- m0
  p[1] <- p0
  p[2] <- p[1] + NL(kappa, kappa_3,v[1]-p[1])+beta*tanh(gamma*m[1]) + rnorm(1,0,sigma_n)
  m[2] <- (1-(alpha_1))*m[1] + alpha_1*(p[2]-p[1])
  S[2] <- (1-alpha_2)*S[1]
  #+ alpha*(p[1])


  for (t in 2:(N-1)) {
    p[t+1] <- p[t] + NL(kappa, kappa_3,v[t]-p[t]) + beta*tanh(gamma*m[t]) + rho*(tanh(voltar - S[t]))  + rnorm(1,0,sigma_n)
    rmean <- mean(diff(p))
    S[t+1] <- (1-alpha_2)*S[t] + alpha_2*((p[t]-p[t-1]) - rmean)^2
    m[t+1] <- (1-alpha_1)*m[t]+ alpha_1*(p[t]-p[t-1])
    # voltar = voltar*(p[t+1] - p[t])
  }
  retruns <- diff(p)
  volsig <- rep(voltar, length(S)) - S
  t <-list(price = p,value = v ,momentum = m, volsig= volsig , returns = retruns, fundsig = v-p, funddemand = kappa*(v-p) + kappa_3*(v-p)^3 , trenddemand = beta*tanh(gamma*m), voldemand = rho*(tanh(0.5*(volsig)))  )

  return(t)

}
