//
// DTW.h: Dynamic Time Warping
//
// Akinori Ito, 2016/6/3
//
#ifndef DTW_H
#define DTW_H
#include "MyMat.h"

class DTWInfo {
  public:
    int xsize;
    int ysize;
    MyMatrix<int> *opt;
    DTWInfo() {
      opt = 0;
    }
    ~DTWInfo() {
      if (opt != 0)
        delete opt;
    }
};

void dtwbase(MyMatrix<double>& x, MyMatrix<double>& y, int window, DTWInfo &info);

#endif
