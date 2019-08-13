
#include <Rcpp.h>
#include <math.h>
using namespace Rcpp;

// [[Rcpp::plugins(cpp11)]]

// [[Rcpp::export]]
NumericVector CaviarLoop(NumericVector Beta, NumericVector y, double EmpiricalQuantile)
{
  int p = y.size();
  int i;
  NumericVector VaR(p);
  VaR[0] = EmpiricalQuantile;
  for(i = 1; i < p; i++)
  {
    /* CAViaR */
    VaR[i] = Beta[0] + Beta[1] * VaR[i-1] + Beta[2] * (y[i-1]);
  }
  return(VaR);
}
