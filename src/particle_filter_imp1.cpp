// [[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include <RcppArmadilloExtensions/sample.h>
#include <math.h>
using namespace arma;
using namespace Rcpp;

double NL(double kappa1, double kappa3, double x){
  double out = kappa1*x + kappa3*(pow(x,3));
  return(out);
}

void resampleParticles(Rcpp::NumericVector& rParticles,
                       Rcpp::NumericVector& particles,
                       Rcpp::NumericVector& weights,
                       int n_particles){
  // resample according to 
  Rcpp::IntegerVector indices = Rcpp::sample(n_particles,n_particles, true, weights, false);
  for(int i = 0; i < n_particles ; i++) {
    rParticles[i] = particles[indices(i)];
  }
}

void initialiseVariables(Rcpp::NumericVector& resampledParticles ,double& initialState, int n_particles){
  for(int col = 0 ; col < n_particles; col++) {
    resampledParticles[col] = initialState;
  }
}

void importance_sampler(double& particle,
                        double& rParticle,
                        Rcpp::NumericVector params,
                        double mt,
                        double pt1,
                        double pt,
                        double& m_t,
                        double& sig_t2)
{
  
  double v = rParticle;
  double disp = v - pt1;
  
  double pow1 = (v + params[6] - pt1,2)*(v + params[6] - pt1,2);
  double pow2 = (params[1] + 3*params[2]*pow1)*(params[1] + 3*params[2]*pow1);
  sig_t2 = 1/(1/params[4] +pow2*(1/params[5]));
  
  m_t = sig_t2*( (1/params[4])*(v + params[6]) + (params[1] + 3*params[2]*pow1)*(1/params[5])*(pt - ( pt1 + NL(params[1], params[2],disp) + params[3]*(tanh(params[8]*mt))) + (params[1] + 3*params[2]*pow1)*(v+params[6])));
  
  double v1 = R::rnorm(m_t,sig_t2);
  particle = v1;
  
}

void propagateParticles(Rcpp::NumericVector& particles,
                        Rcpp::NumericVector& resampledParticles,
                        Rcpp::NumericVector& params,
                        double mt,
                        double pt1,
                        double pt,
                        Rcpp::NumericVector& m_t,
                        Rcpp::NumericVector& sig_t2,
                        int n_particles){
  for(int i = 0; i < n_particles; i++) {
    importance_sampler(particles[i], resampledParticles[i], params, mt, pt1, pt, m_t[i], sig_t2[i]);
  }
}





void weightParticles(double pt,
                     double pt1,
                     double mt,
                     Rcpp::NumericVector m_t,
                     Rcpp::NumericVector sig_t2,
                     Rcpp::NumericVector& weights,
                     Rcpp::NumericVector& particles,
                     Rcpp::NumericVector& rparticles,
                     Rcpp::NumericVector params,
                     int n_particles){
  for(int i = 0; i < n_particles; i++) {
    
    //cout << m_t;
    //cout << "  ";
    //cout << sig_t2;
    //cout << "  ";
    double v = particles[i];
    double v1 = rparticles[i];
    
    double disp = v - pt1;
    weights[i] = (R::dnorm(pt, pt1 + params[1]*disp + params[2]*(disp*disp*disp) + params[3]*(tanh(params[8]*mt)) , params[5], 1) + R::dnorm(v,v1 + params[6], params[4], 1)) - R::dnorm(v ,m_t[i],sig_t2[i],1);
  }
}


void normaliseWeights(Rcpp::NumericVector& tWeights,
                      Rcpp::NumericVector& nWeights) {
  nWeights = exp(tWeights) / sum(exp(tWeights));
}

double computeLikelihood(Rcpp::NumericVector& tWeights, double lWeightsMax, int n_particles){
  double tLikelihood = lWeightsMax + log(sum(exp(tWeights))) - log(n_particles);
  return tLikelihood;
}



// Particle filter for Extended Chiarella model
// Param vectors take form (N, kappa1, kappa3, beta, sigma_v, sigma_N, g, alpha, gamma)
// State takes the form v
// [[Rcpp::export(name=particleFilterCpp)]]
List particleFilter(Rcpp::NumericVector p,
                    double pinit,
                    Rcpp::NumericVector m,
                    Rcpp::NumericVector params,
                    int n_particles,
                    double initialState
) {
  
  // Initialising
  // Marshalling
  double logLikelihood = 0;
  int n_obs = p.size();
  
  // Create Procedure Variables
  Rcpp::NumericVector particles( n_particles);
  Rcpp::NumericVector particles_1( n_particles);
  Rcpp::NumericVector resampledParticles( n_particles);
  Rcpp::NumericVector weights(n_particles);
  Rcpp::NumericVector tWeights(n_particles);
  Rcpp::NumericVector nWeights(n_particles);
  double lWeightsMax;
  Rcpp::NumericVector sig_t2(n_particles);
  Rcpp::NumericVector m_t(n_particles);
  
  
  Rcpp::NumericVector ess(n_obs);
  Rcpp::NumericVector postmeans(n_obs);
  Rcpp::NumericVector stddev(n_obs);
  initialiseVariables(resampledParticles, initialState,n_particles);
  
  
  for(int t = 0; t < (n_obs); t++) {
    
    ess(t) = 0;
    
    if(t == 0){
      propagateParticles(particles, resampledParticles, params,m[t], pinit, p[t],m_t, sig_t2,n_particles);

      weightParticles(p[t], pinit,m[t],m_t, sig_t2, weights, particles,resampledParticles, params,n_particles);
    } else{
      propagateParticles(particles, resampledParticles, params,m[t], p[t-1], p[t],m_t, sig_t2,n_particles);
      weightParticles(p[t], p[t-1],m[t],m_t, sig_t2, weights, particles,resampledParticles, params,n_particles); 
    }
    
    
    lWeightsMax = Rcpp::max(weights);
    tWeights = weights - lWeightsMax;
    normaliseWeights(tWeights, nWeights);
    
    ess(t) = 1/(sum(pow(nWeights,2)));
    resampleParticles(resampledParticles, particles, nWeights,n_particles);
    
    postmeans(t) = mean(resampledParticles);
    stddev(t) = Rcpp::sd(resampledParticles);
    
    
    
    double ll = computeLikelihood(tWeights, lWeightsMax, n_particles);
    logLikelihood += ll;
  }
  
  List L = List::create(Named("mean") = postmeans, Named("Loglik") = logLikelihood, Named("stddev") = stddev, Named("ess") = ess);
  
  return L;
}