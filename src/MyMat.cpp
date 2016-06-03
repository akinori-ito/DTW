#include "MyMat.h"
#include <stdio.h>
#include <stdlib.h>

template<class T>
T& MyVector<T>::at(int i) {
    if (i < 0 || _size <= i) {
        fprintf(stderr,"Subscript out of range: %d\n",i);
        abort();
    }
    return _data[i]; 
}

template<class T>
T& MyMatrix<T>::at(int i, int j) {
    if (i < 0 || _nrow <= i || j < 0 || _ncol <= j) {
        fprintf(stderr,"Subscript out of range: %d,%d\n",i,j);
        abort();
    }
    return _data[index(i,j)];
}

template class MyVector<int>;
template class MyVector<double>;
template class MyMatrix<int>;
template class MyMatrix<double>;
