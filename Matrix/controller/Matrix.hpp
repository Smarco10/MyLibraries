//
//  Matrix.h
//
//  Created by Marc-Antoine MARTIN on 09/04/2016.
//  Copyright © 2016 Marc-Antoine MARTIN. All rights reserved.
//

#ifndef Matrix_h
#define Matrix_h

#include "AbstractMatrix.hpp"

template <typename T, ushort R, ushort C>
class Matrix: public AbstractMatrix<T, R, C> {
private:
    T values[R][C];
    
public:
    Matrix() = default;
    Matrix(const T values[R][C]);
    Matrix(const std::initializer_list<std::initializer_list<T>>& values);
    
    template <ushort R2, ushort C2, ushort ROFF=0U, ushort COFF=0U>
    Matrix(const AbstractMatrix<T, R2, C2>& mat);
    
    virtual T(&operator[](const ushort r)) [C] override;
    virtual const T(&operator[](const ushort r) const) [C] override;
    
    Matrix<T, C, R> transpose() const;
    static Matrix<T, R, C> zeros();
    
    template <ushort C2>
    Matrix<T, R, C + C2> increase(const AbstractMatrix<T, R, C2>& mat) const;
    
    template <ushort C2>
    Matrix<T, R, C + C2> operator<<(const AbstractMatrix<T, R, C2>& mat);
    
    template <ushort R2, ushort C2>
    Matrix<T, R2, C2> submatrix(ushort row, ushort col) const;
};

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrix<T,R,C>::Matrix(const T values[R][C]):
    AbstractMatrix<T, R, C>()
{
    AbstractMatrix<T, R, C>::operator=(values);
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrix<T,R,C>::Matrix(const std::initializer_list<std::initializer_list<T>>& values):
    AbstractMatrix<T, R, C>()
{
    AbstractMatrix<T, R, C>::operator=(values);
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
template <ushort R2, ushort C2, ushort ROFF, ushort COFF>
Matrix<T,R,C>::Matrix(const AbstractMatrix<T, R2, C2>& mat):
    AbstractMatrix<T, R, C>()
{
    for(ushort r = ROFF; r < (R < (ROFF + R2) ? R : (ROFF + R2)); ++r)
        for(ushort c = COFF; c < (C < (COFF + C2) ? C : (COFF + C2)); ++c)
            operator[](r)[c] = mat[r - ROFF][c - COFF];
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
T (&Matrix<T,R,C>::operator[](const ushort r)) [C] {
    return values[r];
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
const T (&Matrix<T,R,C>::operator[](const ushort r) const) [C]  {
    return values[r];
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrix<T, C, R> Matrix<T,R,C>::transpose() const {
    Matrix<T, C, R> mat;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            mat[c][r] = operator[](r)[c];
    
    return mat;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrix<T, R, C> Matrix<T, R, C>::zeros() {
    
    Matrix<T, R, C> mout;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            mout[r][c] = static_cast<T>(0);
    
    return mout;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
template <ushort C2>
Matrix<T, R, C + C2> Matrix<T, R, C>::increase(const AbstractMatrix<T, R, C2>& mat) const {
    
    Matrix<T, R, C + C2> mout;
    
    //copie la matrice
    for(ushort j = 0U; j < R; ++j){
        for(ushort i = 0U; i < C; ++i)
            mout[j][i] = operator[](j)[i];
        
        //augmente de la seconde matrice
        for(ushort i = C; i < C+C2; ++i)
            mout[j][i] = mat[j][i - C];
    }
    
    return mout;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
template <ushort C2>
Matrix<T, R, C + C2> Matrix<T, R, C>::operator<<(const AbstractMatrix<T, R, C2>& mat){
    return increase<C2>(mat);
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
template <ushort R2, ushort C2>
Matrix<T, R2, C2> Matrix<T, R, C>::submatrix(const ushort row, const ushort col) const {
    
    Matrix<T, R2, C2> mout;
    
    for(ushort r = row; r < row + R2 && r < R; ++r)
        for(ushort c = col; c < col + C2 && c < C; ++c)
            mout[r - row][c - col] = operator[](r)[c];
    
    return mout;
}

#endif /* Matrix_h */
