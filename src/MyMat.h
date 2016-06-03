// Original matrix manipulation library
// by Akinori Ito
// 2 June, 2016

#ifndef MYMAT_H
#define MYMAT_H

template<class T>
class MyVector {
    T *_data;
    int _size;
public:
    MyVector(int n) {
        _data = new T[n];
        _size = n;
    }
    ~MyVector() {
        delete[] _data;
    }
    T& at(int i);
    T& operator()(int i) { return at(i); }
    int size() {
        return _size;
   }
};

template<class T>
class MyMatrix {
    T* _data;
    int _nrow;
    int _ncol;
    int index(int i, int j) {
        return i*_ncol+j;
    }
public:
    MyMatrix(int r, int c) {
        _data = new T[r*c];
        _nrow = r;
        _ncol = c;
    }
    ~MyMatrix() {
        delete[] _data;
    }
    T& at(int i, int j);
    T& operator() (int i, int j) { return at(i,j); }
    int nrow() { return _nrow; }
    int ncol() { return _ncol; }
};
#endif
