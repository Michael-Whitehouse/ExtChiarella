#' Title
#'
#' @param kappa description
#' @param kappa3 description
#' @param beta description
#' @param sigma_N description
#' @param sigma_v description
#' @param g description
#' @param tau description
#' @param rho description
#' @param p1 description
#' @param v1 description
#' @param v0 description
#' @param t0 description
#' @param S description
#' @param u description
#'
#' @return
#' @export
#'
#' @examples f
grad_volNLEC <- function(kappa, kappa3 ,beta, sigma_N ,sigma_v  , g , tau, rho, p1, v1,v0,t0,S,u){

  vT <- tail(v1,1)
  nn <- length(p1[t0:length(p1)])
  p <- p1[t0:length(p1)]
  p_ <- p1[(t0-1):(length(p1)-1)]
  v <- v1
  v_ <- head(append(v0,v1),-1)
  vp <- append(v1,vT)[-1]

  m <- msignal(1/7, p1, t0)




  grad <- (1/sigma_v)^2*(v-v_-g-(vp-v-g)) - (1/(sigma_N^2))*(rep(kappa,nn) + 3*kappa3*(v-p_)^2)*(p-p_ - kappa*(v-p_) - kappa3*(v-p_)^3 - beta*u - rho*(tanh(tau - S)))
  return(grad)
}







#' Title
#'
#' @param kappa description
#' @param kappa3 description
#' @param beta description
#' @param sigma_N description
#' @param sigma_v description
#' @param g description
#' @param tau description
#' @param rho description
#' @param p description
#' @param epsilon description
#' @param vinit description
#' @param v0 description
#' @param t0 description
#' @param S description
#' @param u description
#'
#' @return
#' @export
#'
#' @examples f
NonLinear_volNewton <- function(kappa, kappa3 ,beta, sigma_N ,sigma_v  , g , tau, rho ,p,epsilon,vinit, v0,t0,S,u){
  p1 <- p

  nn <- length(p[t0:length(p)])
  Hess <- matrix(nrow = nn, ncol = nn, data = 0)
  p <- p1[t0:length(p1)]
  p_ <- p1[(t0-1):(length(p1)-1)]
  v <- vinit
  v_ <- head(append(v0,v),-1)
  T_ <-  length(p)
  aa <- rep(2/sigma_v^2, nn) + (1/(sigma_N^2))*((rep(kappa,nn) + 3*kappa3*(v-p_)^2)^2 - 6*kappa3*(v-p_)*( p - p_ - kappa*(v - p_) - kappa3*(v - p_)^3 - beta*u- rho*(tanh(rep(tau, length(S)) - S))))
  bb <- rep(-1/sigma_v^2, nn-1)
  diag(Hess) <- aa
  Hess[cbind(1:(nn-1), 2:nn)] <- bb
  Hess[cbind(2:nn, 1:(nn-1))] <- bb
  gr <- grad_volNLEC(kappa, kappa3 ,beta, sigma_N ,sigma_v  , g ,tau, rho, p1, v1 = vinit,v0,t0,S,u)
  iter <- 0
  obj <- c()
  obj[1] <-  -(T_/2)*log(2*pi*sigma_v^2) - (1/(2*sigma_v^2))*sum((v-v_-rep(g, length(v)))^2) -(T_/2)*log(2*pi*sigma_N^2) - (1/(2*sigma_N^2))*sum(( p - p_ - kappa*(v - p_) - kappa3*(v - p_)^3 - beta*u- rho*(tanh(rep(tau, length(S)) - S)))^2)
  while(norm(as.matrix(gr), type='F') > epsilon){
    iter <- iter+1
    gr <- grad_volNLEC(kappa, kappa3 ,beta, sigma_N ,sigma_v  , g ,tau, rho,p1, v1 = v,v0,t0,S,u)
    aa <- rep(2/sigma_v^2, nn) + (1/(sigma_N^2))*((rep(kappa,nn) + 3*kappa3*(v-p_)^2)^2 - 6*kappa3*(v-p_)*( p - p_ - kappa*(v - p_) - kappa3*(v - p_)^3 - beta*u- rho*(tanh(rep(tau, length(S)) - S))))
    v <- v - Solve.tridiag(bb,aa,bb,gr)
    v_ <- head(append(v0,v),-1)
    obj[iter+1] <- -(T_/2)*log(2*pi*sigma_v^2) - (1/(2*sigma_v^2))*sum((v-v_-rep(g, length(v)))^2) -(T_/2)*log(2*pi*sigma_N^2) - (1/(2*sigma_N^2))*sum(( p - p_ - kappa*(v - p_) - kappa3*(v - p_)^3 - beta*u- rho*(tanh(rep(tau, length(S)) - S)))^2)
    if(iter>1000) break
  }
  diag(Hess) <- aa
  out <- list(v_star = v, Hess = Hess, iter = iter, gr = gr, obj = obj)
  return(out)
}


#' Title
#'
#' @param kappa description
#' @param kappa3 description
#' @param beta description
#' @param sigma_N description
#' @param sigma_v description
#' @param tau description
#' @param rho description
#' @param g description
#' @param p1 description
#' @param epsilon description
#' @param vinit description
#' @param v0 description
#' @param t0 description
#' @param alpha1 description
#' @param alpha2 description
#'
#' @return
#' @export
#'
#' @examples f
LaplacevolNLApprox <- function(kappa, kappa3 ,beta, sigma_N ,sigma_v , tau, rho, g ,p1,epsilon,vinit,v0, t0, alpha1,alpha2,gam){
  m <- msignal(alpha = alpha1, p1, t0)
  u <- tanh(gam*m)
  S <- volsig(p1,t0,alpha_2 = alpha2)
  p <- p1[t0:length(p1)]
  p_ <- p1[(t0-1):(length(p1)-1)]
  p0 <- p1[t0-1]
  temp <- NonLinear_volNewton(kappa, kappa3 ,beta, sigma_N ,sigma_v , tau, rho, g ,p1,epsilon,vinit, v0,t0,S,u)

  v <- temp$v_star
  v_ <- append(v0,v[-length(v)])

  T_ <-  length(p)
  lik <- -0.5*logdet(temp$Hess) -(T_/2)*log(2*pi*sigma_v^2) - (1/(2*sigma_v^2))*sum((v-v_-rep(g, length(v)))^2) -(T_/2)*log(2*pi*sigma_N^2) - (1/(2*sigma_N^2))*sum(( p - p_ - kappa*(v - p_) - kappa3*(v - p_)^3 - beta*u- rho*(tanh(rep(tau, length(S)) - S)))^2)

  return(lik)
}

#' Title
#'
#' @param kappa description
#' @param beta description
#' @param sigma_N description
#' @param p1 description
#' @param v1 description
#' @param v0 description
#' @param t0 description
#'
#' @return
#' @export
#'
#' @examples f
grad_EC <- function(kappa,beta,sigma_N,p1, v1,v0,t0){

  vT <- tail(v1,1)

  p <- p1[t0:length(p1)]
  p_ <- p1[(t0-1):(length(p1)-1)]
  v <- v1
  v_ <- head(append(v0,v1),-1)
  vp <- append(v1,vT)[-1]

  m <- msignal(1/7, p, 50)
  u <- tanh(50*m)

  grad <- (1/sigma_v)^2*(v-v_-g-(vp-v-g)) - (1/(sigma_N^2))*kappa*(p-p_ - kappa*(v-p_) - beta*u)
  return(grad)
}


#' Title
#'
#' @param kappa description
#' @param beta description
#' @param sigma_N description
#' @param p description
#' @param epsilon description
#' @param vinit description
#' @param v0 description
#' @param t0 description
#'
#' @return
#' @export
#'
#' @examples f
LinearNewton <- function(kappa,beta,sigma_N,p,epsilon,vinit, v0,t0){
  m <- msignal(1/7, t$price, t0)
  u <- tanh(50*m)
  nn <- length(p[t0:length(p)])
  Hess <- matrix(nrow = nn, ncol = nn, data = 0)
  aa <- rep(2/sigma_v^2 + kappa^2/sigma_N^2, nn)
  bb <- rep(-1/sigma_v^2, nn-1)
  diag(Hess) <- aa
  Hess[cbind(1:(nn-1), 2:nn)] <- bb
  Hess[cbind(2:nn, 1:(nn-1))] <- bb
  v <- vinit
  gr <- grad_EC(kappa,beta,sigma_N,p, v, v0, t0)
  iter <- 0
  while(norm(as.matrix(gr), type='F') > epsilon){
    iter <- iter+1
    gr <- grad_EC(kappa,beta,sigma_N = sigma_N,p, v, v0, t0)
    v <- v - Solve.tridiag(bb,aa,bb,gr)
    if(iter>1000) break
  }
  out <- list(v_star = v, Hess = Hess, iter = iter)
  return(out)
}


#' Title
#'
#' @param kappa description
#' @param beta description
#' @param sigma_N description
#' @param p1 description
#' @param epsilon description
#' @param vinit description
#' @param v0 description
#' @param t0 description
#'
#' @return
#' @export
#'
#' @examples f
LaplaceApprox <- function(kappa,beta, sigma_N,p1,epsilon,vinit,v0, t0){
  m <- msignal(1/7, t$price, t0)
  u <- tanh(50*m)
  p <- p1[t0:length(p1)]
  p_ <- p1[(t0-1):(length(p1)-1)]
  p0 <- p1[t0-1]
  temp <- LinearNewton(kappa,beta,sigma_N,p1,epsilon = epsilon, vinit = vinit, v0 = v0, t0=t0)

  v <- temp$v_star
  v_ <- append(v0,v[-length(v)])

  T_ <-  length(p)
  lik <-  lik <- -0.5*logdet(temp$Hess) -(T_/2)*log(2*pi*sigma_v^2) - (1/(2^sigma_v^2))*sum((v-v_-rep(g, length(v)))^2) -(T_/2)*log(2*pi*sigma_N^2) - (1/(2^sigma_N^2))*sum(( p - p_ - kappa*(v - p_) - beta*u)^2)

  return(lik)
}


#' Title
#'
#' @param kappa description
#' @param kappa_3 description
#' @param x description
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
#' @param kappa description
#' @param kappa3 description
#' @param beta description
#' @param sigma_N description
#' @param sigma_v description
#' @param g description
#' @param p1 description
#' @param v1 description
#' @param v0 description
#' @param t0 description
#'
#' @return
#' @export
#'
#' @examples f
grad_NLEC <- function(kappa, kappa3 ,beta, sigma_N ,sigma_v  , g ,p1, v1,v0,t0){

  vT <- tail(v1,1)
  nn <- length(p1[t0:length(p1)])
  p <- p1[t0:length(p1)]
  p_ <- p1[(t0-1):(length(p1)-1)]
  v <- v1
  v_ <- head(append(v0,v1),-1)
  vp <- append(v1,vT)[-1]

  m <- msignal(1/7, p1, t0)
  u <- tanh(50*m)

  grad <- (1/sigma_v)^2*(v-v_-g-(vp-v-g)) - (1/(sigma_N^2))*(rep(kappa,nn) + 3*kappa3*(v-p_)^2)*(p-p_ - kappa*(v-p_) - kappa3*(v-p_)^3 - beta*u)
  return(grad)
}


#grad_NLEC(kappa = 0.01, kappa3 = 0.08 , beta = 0.1, p = tNL$price, v1 =  tNL$value[-c(1:49)],v0 = tNL$value[49], t0 = 50)



#' Title
#'
#' @param kappa description
#' @param kappa3 description
#' @param beta description
#' @param sigma_N description
#' @param sigma_v description
#' @param g description
#' @param p description
#' @param epsilon description
#' @param vinit description
#' @param v0 description
#' @param t0 description
#'
#' @return
#' @export
#'
#' @examples f
NonLinear_Newton <- function(kappa, kappa3 ,beta, sigma_N ,sigma_v  , g ,p,epsilon,vinit, v0,t0){
  p1 <- p
  m <- msignal(1/7, p=p1, t0)

  u <- tanh(50*m)
  nn <- length(p[t0:length(p)])
  Hess <- matrix(nrow = nn, ncol = nn, data = 0)
  p <- p1[t0:length(p1)]
  T_ <-  length(p)
  p_ <- p1[(t0-1):(length(p1)-1)]
  v <- vinit
  v_ <- head(append(v0,v),-1)
  aa <- rep(2/sigma_v^2, nn) + (1/(sigma_N^2))*((rep(kappa,nn) + 3*kappa3*(v-p_)^2)^2 - 6*kappa3*(v-p_)*( p - p_ - kappa*(v - p_) - kappa3*(v - p_)^3 - beta*u))
  bb <- rep(-1/sigma_v^2, nn-1)
  diag(Hess) <- aa
  Hess[cbind(1:(nn-1), 2:nn)] <- bb
  Hess[cbind(2:nn, 1:(nn-1))] <- bb
  gr <- grad_NLEC(kappa, kappa3 ,beta, sigma_N ,sigma_v  , g ,p1, v1 = vinit,v0,t0)
  obj <- c()
  obj[1] <- -(T_/2)*log(2*pi*sigma_v^2) - (1/(2*sigma_v^2))*sum((v-v_-rep(g, length(v)))^2) -(T_/2)*log(2*pi*sigma_N^2) - (1/(2*sigma_N^2))*sum(( p - p_ - kappa*(v - p_) - kappa3*(v - p_)^3 - beta*u)^2)
  iter <- 0
  while(norm(as.matrix(gr), type='F') > epsilon){
    iter <- iter+1
    gr <- grad_NLEC(kappa, kappa3 ,beta, sigma_N ,sigma_v  , g ,p1, v1 = v,v0,t0)
    aa <- rep(2/sigma_v^2, nn) + (1/(sigma_N^2))*((rep(kappa,nn) + 3*kappa3*(v-p_)^2)^2 - 6*kappa3*(v-p_)*( p - p_ - kappa*(v - p_) - kappa3*(v - p_)^3 - beta*u))
    v <- v - Solve.tridiag(bb,aa,bb,gr)
    v_ <- head(append(v0,v),-1)
    obj[iter+1] <- -(T_/2)*log(2*pi*sigma_v^2) - (1/(2*sigma_v^2))*sum((v-v_-rep(g, length(v)))^2) -(T_/2)*log(2*pi*sigma_N^2) - (1/(2*sigma_N^2))*sum(( p - p_ - kappa*(v - p_) - kappa3*(v - p_)^3 - beta*u)^2)
    if(iter>1000) break
  }
  diag(Hess) <- aa
  out <- list(v_star = v, Hess = Hess, iter = iter, gr = gr, obj=obj)
  return(out)
}


#' Title
#'
#' @param kappa description
#' @param kappa3 description
#' @param beta description
#' @param sigma_N description
#' @param sigma_v description
#' @param g description
#' @param p1 description
#' @param epsilon description
#' @param vinit description
#' @param v0 description
#' @param t0 description
#'
#' @return
#' @export
#'
#' @examples f
LaplaceNLApprox <- function(kappa, kappa3 ,beta, sigma_N ,sigma_v  , g ,p1,epsilon,vinit,v0, t0){
  m <- msignal(1/7, p1, t0)
  u <- tanh(50*m)
  p <- p1[t0:length(p1)]
  p_ <- p1[(t0-1):(length(p1)-1)]
  p0 <- p1[t0-1]
  temp <- NonLinear_Newton(kappa, kappa3 ,beta, sigma_N ,sigma_v  , g ,p1,epsilon,vinit, v0,t0)

  v <- temp$v_star
  v_ <- append(v0,v[-length(v)])

  T_ <-  length(p)
  lik <- -0.5*logdet(temp$Hess) -(T_/2)*log(2*pi*sigma_v^2) - (1/(2*sigma_v^2))*sum((v-v_-rep(g, length(v)))^2) -(T_/2)*log(2*pi*sigma_N^2) - (1/(2*sigma_N^2))*sum(( p - p_ - kappa*(v - p_) - kappa3*(v - p_)^3 - beta*u)^2)

  return(lik)
}
