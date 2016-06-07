#include <cmath>
//#include <Rcpp.h>
#include <vector>
#include <cstdio>
#include "MyMat.h"
#include "dtw.h"

//using namespace Rcpp;

#define j_ind(i,j) ((j)-jcenter(i)+window)
//#define j_ind(i,j) (j)

typedef int intpair[2];

double dist2(MyMatrix<double>& x, MyMatrix<double>& y, int i, int j) {
  double r = 0.0;
  for (int k = 0; k < x.ncol(); k++) {
    double d = x(i,k)-y(j,k);
    r += d*d;
  }
  return(r);
}

int which_min(double *g, int n) {
  int m = 0;
  double mv = g[0];
  for (int i = 1; i < n; i++) {
    if (g[i] < mv) {
      m = i;
      mv = g[i];
    }
  }
  return m;
}

bool inlimit(int x, int limit) {
  return 0 <= x && x < limit;
}

void dtwbase(MyMatrix<double>& x, MyMatrix<double>& y,
             int window, DTWInfo& info) {
  //
  // G[i,j] <=> g[i,j-jcenter[i]+window+1]
  //
  const double Large = 1.0e7; //arbitrary
  int x_size = x.nrow();
  int y_size = y.nrow();
  int wlimit = y_size; //2*window+1;
  MyMatrix<double> g(x_size,wlimit); //<- array(Inf,dim=c(x_size,wlimit))
  MyMatrix<double> d(x_size,wlimit); // <- array(Inf,dim=c(x_size,wlimit))
  MyMatrix<int> bp(x_size,wlimit); // <- array(0,dim=c(x_size,wlimit))
  for (int i = 0; i < x_size; i++) {
    for (int j = 0; j < wlimit; j++) {
      g(i,j) = Large;
      d(i,j) = Large;
      bp(i,j) = 100; //meaningless value
    }
  }
  MyVector<int> jcenter(x_size);
  for (int i = 0; i < x_size; i++) {
    jcenter(i) = floor((double)i/(x_size-1)*(y_size-1)) ;
  }
  g(0,j_ind(0,0)) = d(0,j_ind(0,0)) = dist2(x,y,0,0);
  bp(0,j_ind(0,0)) = 1;
  d(1,j_ind(1,1)) = dist2(x,y,1,1);
  g(1,j_ind(1,1)) = g(0,j_ind(0,0))+d(1,j_ind(1,1));
  bp(1,j_ind(1,1)) = 1;

  for (int i = 2; i < x_size; i++) {

    int jmin = jcenter(i)-window;
    if (jmin < 0) jmin = 0;
    int jmax = jcenter(i)+window;
    if (jmax > y_size-1) jmax = y_size-1;
//printf("%d: jmin=%d jmax=%d\n",i,jmin,jmax);
/*
      int jmin = 0;
      int jmax = y_size-1;
*/
      for (int j = jmin; j <= jmax; j++) {
      int w = j_ind(i,j);
      int w1 = j_ind(i-1,j);
      int w2 = j_ind(i-2,j);
//printf("i=%d j=%d w=%d w1=%d w2=%d\n",i,j,w,w1,w2);
      d(i,w) = dist2(x,y,i,j);
      double gm[3];
      for (int k = 0; k < 3; k++)
        gm[k] = Large;
      if (inlimit(w1,wlimit) && inlimit(w2-1,wlimit)) {
        gm[0] = g(i-2,w2-1)+d(i-1,w1);
      }
      if (inlimit(w1-1,wlimit)) {
        gm[1] = g(i-1,w1-1);
      }
      if (inlimit(w1-2,wlimit) && inlimit(w-1,wlimit)) {
        gm[2] = g(i-1,w1-2)+d(i,w-1);
      }
      int m = which_min(gm,3);
      if (gm[m] >= Large) {
        //printf("Warning: no candidate at %d,%d (physically %d)\n",i,j,j_ind(i,j));
      }
      else {
        g(i,w) = gm[m]+d(i,w);
        bp(i,w) = m;
      }
    }
  }
  //printf("backtracing\n");
  int i = x_size-1;
  int j = y_size-1;
  std::vector<int> optx;
  std::vector<int> opty;
  optx.push_back(i);
  opty.push_back(j);
  while (i > 0 && j > 0) {
//printf("%d %d %d\n",i,j,bp(i,j_ind(i,j)));
    int b = bp(i,j_ind(i,j));
    if (b == 0) {
      optx.push_back(i-1);
      optx.push_back(i-2);
      opty.push_back(j-1);
      opty.push_back(j-1);
      i -= 2;
      j--;
    }
    else if (b == 1) {
      optx.push_back(i-1);
      opty.push_back(j-1);
      i--;
      j--;
    }
    else if (b == 2) {
      optx.push_back(i-1);
      optx.push_back(i-1);
      opty.push_back(j-1);
      opty.push_back(j-2);
        i--;
        j -= 2;
    }
    else {
      printf("Invalid backpointer: %d at %d,%d (physically %d)\n",b,i,j,j_ind(i,j));
      break;
    }
  }
  int opt_size = optx.size();
  MyMatrix<int>* opt = new MyMatrix<int>(opt_size,2);
  for (int i = 0; i < opt_size; i++) {
    opt->at(i,0) = optx[opt_size-i-1]+1;
    opt->at(i,1) = opty[opt_size-i-1]+1;
  }
  info.xsize = x_size;
  info.ysize = y_size;
  info.opt = opt;
}
