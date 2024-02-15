#include <Rcpp.h>

using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

NumericVector pdf_exp(NumericVector t, NumericVector params) {return dweibull(t,1,params(0));}
NumericVector cdf_exp(NumericVector t, NumericVector params) {return pweibull(t,1,params(0));}
NumericVector pdf_gamma(NumericVector t, NumericVector params) {return dgamma(t,params(1),params(2));}
NumericVector cdf_gamma(NumericVector t, NumericVector params) {return pgamma(t,params(1),params(2));}
NumericVector pdf_unif(NumericVector t, NumericVector params) {return dunif(t,params(6),params(7));}
NumericVector cdf_unif(NumericVector t, NumericVector params) {return punif(t,params(6),params(7));}
NumericVector pdf_gomp(NumericVector t, NumericVector params) {
  Rcpp::Environment flexsurv("package:flexsurv");
  Rcpp::Function dgompertz=flexsurv["dgompertz"];
  return dgompertz(t,params(1),params(2));
}
NumericVector cdf_gomp(NumericVector t, NumericVector params) {
  Rcpp::Environment flexsurv("package:flexsurv");
  Rcpp::Function pgompertz=flexsurv["pgompertz"];
  return pgompertz(t,params(1),params(2));
}
NumericVector pdf_llog(NumericVector t, NumericVector params) {
  Rcpp::Environment flexsurv("package:flexsurv");
  Rcpp::Function dllogis=flexsurv["dllogis"];
  return dllogis(t,params(1),params(3));
}
NumericVector cdf_llog(NumericVector t, NumericVector params) {
  Rcpp::Environment flexsurv("package:flexsurv");
  Rcpp::Function pllogis=flexsurv["pllogis"];
  return pllogis(t,params(1),params(3));
}
NumericVector pdf_lnorm(NumericVector t, NumericVector params) {return dlnorm(t,params(4),params(5));}
NumericVector cdf_lnorm(NumericVector t, NumericVector params) {return plnorm(t,params(4),params(5));}
NumericVector pdf_weibull(NumericVector t, NumericVector params) {return dweibull(t,params(1),params(3));}
NumericVector cdf_weibull(NumericVector t, NumericVector params) {return pweibull(t,params(1),params(3));}

NumericVector pdf(NumericVector t, const char * val, NumericVector params) {
  NumericVector pdf;
  if(std::strcmp(val,"Exponential")==0) {
    pdf=pdf_exp(t,params);
  } else if(std::strcmp(val,"Gamma")==0) {
    pdf=pdf_gamma(t,params);
  }
  else if(std::strcmp(val,"Lognormal")==0) {
    pdf=pdf_lnorm(t,params);
  } else if(std::strcmp(val,"Weibull")==0) {
    pdf=pdf_weibull(t,params);
  } else if(std::strcmp(val,"Uniform")==0) {
    pdf=pdf_unif(t,params);
  }
  return pdf;
}

NumericVector cdf(NumericVector t, const char * val, NumericVector params) {
  NumericVector cdf;
  if(std::strcmp(val,"Exponential")==0) {
    cdf=cdf_exp(t,params);
  } else if(std::strcmp(val,"Gamma")==0) {
    cdf=cdf_gamma(t,params);
  }
  else if(std::strcmp(val,"Lognormal")==0) {
    cdf=cdf_lnorm(t,params);
  } else if(std::strcmp(val,"Weibull")==0) {
    cdf=cdf_weibull(t,params);
  } else if(std::strcmp(val,"Uniform")==0) {
    cdf=cdf_unif(t,params);
  }
  return cdf;
}

// [[Rcpp::export]]
NumericVector S(NumericVector t, const char * S_fun, NumericVector S_params) {
  return 1-cdf(t, S_fun, S_params);
}

// [[Rcpp::export]]
NumericVector f(NumericVector t, const char * S_fun, NumericVector S_params) {
  return(pdf(t, S_fun, S_params));
}

// [[Rcpp::export]]
NumericVector H(NumericVector t, const char * H_fun, NumericVector H_params) {
  return(cdf(t, H_fun, H_params));
}

// [[Rcpp::export]]
NumericVector h(NumericVector t, const char * H_fun, NumericVector H_params) {
  return(pdf(t, H_fun, H_params));
}

double inte(Function f, double lower, double upper, const char * S_fun, NumericVector S_params, const char * H_fun, NumericVector H_params) {
  Function integ("integrate");
  List result=integ(Named("f")=f, Named("lower")=lower, Named("upper")=upper, 
                    Named("S_fun")=S_fun, Named("S_params")=S_params,
                    Named("H_fun")=H_fun, Named("H_params")=H_params);
  return result["value"];
}

NumericVector integrand_H_obs(NumericVector r, const char * S_fun, NumericVector S_params, const char * H_fun, NumericVector H_params) {
  return H(r,H_fun,H_params)*f(r,S_fun,S_params);
}

NumericVector H_obs(NumericVector t, const char * S_fun, NumericVector S_params, const char * H_fun, NumericVector H_params) {
  Function fun("integrand_H_obs");
  double Inf = -std::numeric_limits<double>::infinity();
  double denom=inte(fun,0,Inf,S_fun,S_params,H_fun,H_params);
  NumericVector t_prev=t[Range(0,t.size()-2)];
  t_prev.push_front(0);
  NumericVector result=cumsum(H(t,H_fun,H_params)*(S(t_prev,S_fun,S_params)-S(t,S_fun,S_params)));
  return pmin((result+H(t,H_fun,H_params)*S(t,S_fun,S_params))/denom,1);
}

NumericVector h_obs(NumericVector t, const char * S_fun, NumericVector S_params, const char * H_fun, NumericVector H_params) {
  Function fun("integrand_H_obs");
  double Inf = -std::numeric_limits<double>::infinity();
  double denom=inte(fun,0,Inf,S_fun,S_params,H_fun,H_params);
  NumericVector result=h(t,H_fun,H_params)*S(t,S_fun,S_params)/denom;
  return result;
}

// [[Rcpp::export]]
NumericVector Y(NumericVector t, int n, double pi, const char * S_fun, NumericVector S_params, const char * H_fun, NumericVector H_params, double theta) {
  NumericVector Y_I=pmax(n*pi*S(t,S_fun,S_params)*(1-t/theta),0);
  NumericVector t_prev=t[Range(0,t.size()-2)];
  t_prev.push_front(0);
  double result=sum(S(t,S_fun,S_params)*(H(t,H_fun,H_params)-H(t_prev,H_fun,H_params)));
  NumericVector Y_P=pmax(n*(1-pi)*S(t,S_fun,S_params)*(H(t,H_fun,H_params)-H(t-theta,H_fun,H_params))/result,0);
  return Y_P+Y_I;
}

NumericVector integrand_var_fun(NumericVector r, int n, double pi, const char * S_fun, NumericVector S_params, const char * H_fun, NumericVector H_params, double theta) {
  NumericVector Y_out=Y(r,n,pi,S_fun,S_params,H_fun,H_params,theta);
  return f(r,S_fun,S_params)/(Y_out*S(r,S_fun,S_params));
}

// [[Rcpp::export]]
NumericVector Var_Fun(NumericVector t, int n, double pi, const char * S_fun, NumericVector S_params, const char * H_fun, NumericVector H_params, double theta) {
  NumericVector Y_out=Y(t,n,pi,S_fun,S_params,H_fun,H_params,theta);
  NumericVector t_prev=t[Range(0,t.size()-2)];
  t_prev.push_front(0);
  NumericVector result=cumsum((1/(Y_out*S(t,S_fun,S_params)))*(S(t_prev,S_fun,S_params)-S(t,S_fun,S_params)));
  NumericVector S_temp=S(t,S_fun,S_params);
  result=result*S_temp*S_temp;
  return result;
}

NumericVector sequence(double from, double to, double by) {
  double ln=1+(to-from)/by;
  NumericVector result(ln);
  for(int i=0; i < result.size(); i++) {
    result(i)=from+i*by;
  }
  return result;
}

// [[Rcpp::export]]
NumericVector W(NumericVector t, double w_shape1, double w_shape2, double tau) {
  Rcpp::Environment ExtDist("package:ExtDist");
  Rcpp::Function dBeta_ab=ExtDist["dBeta_ab"];
  return dBeta_ab(t,w_shape1,w_shape2,0,tau);
}

NumericVector W_helper(NumericVector t, double w_shape1, double w_shape2, double tau) {
  Rcpp::Environment ExtDist("package:ExtDist");
  Rcpp::Function pBeta_ab=ExtDist["pBeta_ab"];
  return pBeta_ab(t,w_shape1,w_shape2,0,tau);
}

// [[Rcpp::export]]
double Opt_Fun(double tau, int n, double pi, const char * S_fun, NumericVector S_params, const char * H_fun, NumericVector H_params,
               double w_shape1, double w_shape2, double theta) {
  NumericVector t=sequence(0.0,tau,0.01);
  Rcpp::Environment ExtDist("package:ExtDist");
  Rcpp::Function pBeta_ab=ExtDist["pBeta_ab"];
  NumericVector t_prev=t[Range(0,t.size()-2)];
  t_prev.push_front(0);
  double result=sum(Var_Fun(t,n,pi,S_fun,S_params,H_fun,H_params,theta)*(W_helper(t,w_shape1,w_shape2,tau)-W_helper(t_prev,w_shape1,w_shape2,tau)));
  return result;
}

// [[Rcpp::export]]
double Sol(double tau, int n, const char * S_fun, NumericVector S_params, const char * H_fun, NumericVector H_params, double w_shape1, double w_shape2, double theta) {
  NumericVector pi_vec=sequence(0.01,0.99,0.01);
  NumericVector result(pi_vec.size());
  for(int i=0; i < pi_vec.size(); i++) {
    NumericVector pi_i;
    pi_i=pi_vec(i);
    double pi_temp=sum(pi_i);
    result(i)=Opt_Fun(tau,n,pi_temp,S_fun,S_params,H_fun,H_params,w_shape1,w_shape2,theta);
  }
  double pi_opt=pi_vec[which_min(result)];
  if(min(result) == 1.0/0.0) {
    pi_opt=-1;
  }
  return pi_opt;
}





