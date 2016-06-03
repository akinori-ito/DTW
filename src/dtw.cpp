#include <Rcpp.h>
#include "DTW.h"
#include "Mymat.h"

using namespace Rcpp;

inline void copymat(MyMatrix<double>& a, NumericMatrix b) {
  for (int i = 0; i < b.nrow(); i++)
    for (int j = 0; j < b.ncol(); j++)
      a(i,j) = b(i,j);
}


// [[Rcpp::export]]
List dtw_(NumericMatrix x, NumericMatrix y, int window) {
  MyMatrix<double> mx(x.nrow(),x.ncol());
  MyMatrix<double> my(y.nrow(),y.ncol());
  DTWInfo info;

  copymat(mx,x);
  copymat(my,y);
  dtwbase(mx,my,window,info);

  IntegerMatrix opt(info.opt->nrow(),2);
  for (int i = 0; i < info.opt->nrow(); i++) {
    opt(i,0) = info.opt->at(i,0);
    opt(i,1) = info.opt->at(i,1);
  }
  List ret;
  ret["xsize"] = x.nrow();
  ret["ysize"] = y.nrow();
  //ret["d"] = d;
  //ret["g"] = g;
  //ret["bp"] = bp;
  ret["opt"] = opt;
  //ret["jcenter"] = jcenter;
  return ret;
}
