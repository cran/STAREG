#include <RcppArmadillo.h>
#include <Rcpp.h>

// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

vec my_pava(vec &values, vec &weight, bool decreasing);
vec replfdr(double &xi00, double &xi01, double &xi10, double &xi11, vec &f1, vec &f2);
vec fdr(vec &Lfdr);

//' @title EM algorithm to estimate local false discovery rate
//' @description Estimate the local false discovery rate across two studies and apply a step-up procedure to control the FDR of replicability null.
//' @param pa_in A numeric vector of p-values from study 1.
//' @param pb_in A numeric vector of p-values from study 2.
//' @param pi0a_in An initial estimate of the null probability in study 1.
//' @param pi0b_in An initial estimate of the null probability in study 2.
//'
//' @return
//' \item{Lfdr}{The estimated local false discovery rate for replicability null.}
//' \item{fdr}{The adjusted values based on local false discovery rate for FDR control.}
//' \item{xi00}{An estimate of the prior probability for joint state (0, 0).}
//' \item{xi01}{An estimate of the prior probability for joint state (0, 1).}
//' \item{xi10}{An estimate of the prior probability for joint state (1, 0).}
//' \item{xi11}{An estimate of the prior probability for joint state (1, 1).}
//' \item{f1}{A non-parametric estimate for the non-null probability density function in study 1.}
//' \item{f2}{A non-parametric estimate for the non-null probability density function in study 2.}
//'
//' @export
// [[Rcpp::export]]
SEXP em_lfdr(SEXP pa_in, SEXP pb_in, SEXP pi0a_in, SEXP pi0b_in) {
  try{
    const int maxIter = 200;
    const double tol = 1e-3;
    const double pvalue_cutoff = 1e-15;
    const double f_cutoff = 1e-15;

    vec pa = as<arma::vec>(pa_in), pb = as<arma::vec>(pb_in);
    const double pi0_pa = Rcpp::as<double>(pi0a_in), pi0_pb = Rcpp::as<double>(pi0b_in);

    int J = pa.size();
    double min_a = pa.elem(find(pa>0)).min(), min_b = pb.elem(find(pb>0)).min();
    pa.elem(find(pa==0)).fill(pvalue_cutoff < min_a ? pvalue_cutoff : min_a);
    pb.elem(find(pb==0)).fill(pvalue_cutoff < min_b ? pvalue_cutoff : min_b);

    vec p1 = sort(pa), p2 = sort(pb);
    uvec ix1 = sort_index(pa), ix2 = sort_index(pb);

    vec p1_diff(J), p2_diff(J);
    p1_diff(0) = p1(0);
    p2_diff(0) = p2(0);
    for(int i = 1; i<J; ++i){
      p1_diff(i) = p1(i) - p1(i-1);
      p2_diff(i) = p2(i) - p2(i-1);
    }

    // Initialization
    double xi00=pi0_pa*pi0_pb, xi01=pi0_pa*(1-pi0_pb), xi10=(1-pi0_pa)*pi0_pb, xi11 = (1-pi0_pa)*(1-pi0_pb);

    vec f0 = ones(J,1), f1(J), f2(J);
    f1 = 1 - pa;
    f2 = 1 - pb;

    vec loglik(maxIter);
    loglik(0) = -datum::inf;

    vec f(J), gamma00(J), gamma01(J), gamma10(J), gamma11(J);
    vec Q1(J), Q2(J), y1(J), y2(J), res1(J), res2(J);
    double loglik_delta;

    //std::cout<<"EM begins:"<<std::endl;
    Rcout << "EM begins:" << "\n";

    for (int i = 1; i < maxIter; i++){
      // E-step
      f = xi00 * f0 % f0 + xi01 * f0 % f2 + xi10 * f1 % f0 + xi11 * f1 % f2;
      gamma00 = xi00 * f0 % f0 / f;
      gamma01 = xi01 * f0 % f2 / f;
      gamma10 = xi10 * f1 % f0 / f;
      gamma11 = 1 - gamma00 - gamma01 - gamma10;

      // M-step
      // update f1 and f2
      Q1 = gamma01 + gamma00;
      Q2 = gamma10 + gamma00;
      Q1 = Q1(ix1);
      Q2 = Q2(ix2);

      vec _Q1 = 1 - Q1, _Q2 = 1 - Q2;
      y1 = - p1_diff * sum(_Q1) / _Q1;
      y2 = - p2_diff * sum(_Q2) / _Q2;

      y1.elem(find_nonfinite(y1)).fill(y1.elem(find_finite(y1)).min());
      y2.elem(find_nonfinite(y2)).fill(y2.elem(find_finite(y2)).min());

      res1 = my_pava(y1, _Q1, true);
      res2 = my_pava(y2, _Q2, true);

      f1 = -1 / res1;
      f1 = f1 / sum(f1 % p1_diff);
      f1(ix1) = f1;
      f1.elem(find_nan(f1)).fill(f1.min());

      f2 = -1 / res2;
      f2 = f2 / sum(f2 % p2_diff);
      f2(ix2) = f2;
      f2.elem(find_nan(f2)).fill(f2.min());

      double min_f1 = f1.elem(find(f1>0)).min(), min_f2 = f2.elem(find(f2>0)).min();
      f1.elem(find(f1<=0)).fill(f_cutoff < min_f1 ? f_cutoff : min_f1);
      f2.elem(find(f2<=0)).fill(f_cutoff < min_f2 ? f_cutoff : min_f2);

      // update xi's
      xi00 = mean(gamma00);
      xi01 = mean(gamma01);
      xi10 = mean(gamma10);
      xi11 = mean(gamma11);

      // calculate the updated log-likelihood
      loglik(i) = sum((gamma10+gamma11)%log(f1)+(gamma01+gamma11)%log(f2)) +
        sum(gamma11*log(xi11)+gamma10*log(xi10)+gamma01*log(xi01)+gamma00*log(xi00));
      loglik_delta = abs((loglik(i) - loglik(i-1))/loglik(i-1));

      //std::cout<<i<<". "<< loglik(i) << ", delta = " << loglik_delta << std::endl;
      Rcout << i << ". " << loglik(i) << ", delta = " << loglik_delta << "\n";

      if(loglik_delta <= tol){
        break;
      }
    }

    vec Lfdr = replfdr(xi00, xi01, xi10, xi11, f1, f2);
    vec lfdr_adj = fdr(Lfdr);

    return Rcpp::List::create(Rcpp::Named("Lfdr") = Lfdr,
                              Rcpp::Named("fdr") = lfdr_adj,
                              Rcpp::Named("loglik") = loglik,
                              Rcpp::Named("xi00") = xi00,
                              Rcpp::Named("xi01") = xi01,
                              Rcpp::Named("xi10") = xi10,
                              Rcpp::Named("xi11") = xi11,
                              Rcpp::Named("f1") = f1,
                              Rcpp::Named("f2") = f2);
  } catch( std::exception &ex ) {
    forward_exception_to_r(ex);
    return Rcpp::List::create();
  } catch(...) {
    ::Rf_error( "C++ exception (unknown reason)..." );
    return Rcpp::List::create();
  }
}

//' @title Plug-in estimates of local false discovery rate
//' @description Calculate the plug-in estimates of local false discovery rate across two studies.
//' @param xi00 The prior probability for joint state (0, 0).
//' @param xi01 The prior probability for joint state (0, 1).
//' @param xi10 The prior probability for joint state (1, 0).
//' @param xi11 The prior probability for joint state (1, 1).
//' @param f1 The non-null probability density function in study 1.
//' @param f2 The non-null probability density function in study 2.
//'
//' @return
//' \item{Lfdr}{The plug-in estimates of local false discovery rate for all features.}
//'
vec replfdr(double &xi00, double &xi01, double &xi10, double &xi11, vec &f1, vec &f2){
  int J = f1.size();
  vec f0 = ones(J,1);

  vec f = xi00 * f0 % f0 + xi01 * f0 % f2 + xi10 * f1 % f0 + xi11 * f1 % f2;
  vec gamma00 = xi00 * f0 % f0 / f;
  vec gamma01 = xi01 * f0 % f2 / f;
  vec gamma10 = xi10 * f1 % f0 / f;
  vec gamma11 = 1 - gamma00 - gamma01 - gamma10;


  vec Lfdr = gamma00 + gamma01 + gamma10;

  return Lfdr;
}

//' @title Calculate the adjusted value of local false discovery rate for FDR control
//' @description Apply a step-up procedure based on the local false discovery rate to obtain the adjusted values for FDR control.
//' @param Lfdr The local false discovery rate for all features.
//'
//' @return
//' \item{lfdr_adj}{The adjusted local false discovery rate values for FDR control.}
//'
vec fdr(vec &Lfdr){
  int m = Lfdr.size();

  vec ordered_lfdr = sort(Lfdr), s = linspace(1,m,m);
  uvec ix_lfdr = sort_index(Lfdr);

  vec lfdr_adj = cumsum(ordered_lfdr)/s;
  lfdr_adj(ix_lfdr) = lfdr_adj;

  return lfdr_adj;
}

//' @title Pool-adjacent-violator-algorithm to fit the isotonic regression.
//' @description The pool-adjacent-violator-algorithm is applied to fit the isotonic regression of a set of data.
//' @param values A numeric vector of data whose isotonic regression is to be calculated.
//' @param weight The weight vector to be used for a weighted isotonic regression.
//' @param decreasing A logical scalar that specifies the order of the isotonic regression (true if decreasing; false otherwise).
//'
//' @return
//' \item{xx}{The fitted values of the data.}
//'
vec my_pava(vec &values, vec &weight, bool decreasing)
{
  if(decreasing){
    values = reverse(values);
    weight = reverse(weight);
  }
  vec w(values.size(), fill::zeros);
  vec x(values.size(), fill::zeros);
  x[0] = values[0];
  w[0] = weight[0];
  unsigned j = 0;
  vec s(values.size(), fill::zeros);

  for (unsigned i = 1; i < values.size(); i++) {
    j += 1;
    x[j] = values[i];
    w[j] = weight[i];
    while (j > 0 && x[j - 1]>x[j]) {
      x[j - 1] = (w[j] * x[j] + w[j - 1] * x[j - 1]) / (w[j] + w[j - 1]);
      w[j - 1] += w[j];
      j -= 1;
    }
    s[j + 1] = i + 1;
  }

  vec ww(values.size(), fill::zeros);
  vec xx(values.size(), fill::zeros);
  for (unsigned k = 0; k < j + 1; k++) {
    for (unsigned i = s[k]; i < s[k + 1]; i++) {
      ww[i] = w[k];
      xx[i] = x[k];
    }
  }

  if(decreasing){
    xx = reverse(xx);
  }

  return xx;
}
