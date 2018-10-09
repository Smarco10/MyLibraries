//
//  MatrixMath.hpp
//  testCpp2
//
//  Created by Marc-Antoine MARTIN on 07/10/2018.
//  Copyright © 2018 Marc-Antoine MARTIN. All rights reserved.
//

#ifndef MatrixMath_h
#define MatrixMath_h

#include <limits>

template<typename T, ushort R, ushort C>
class AbstractMatrix;

template<typename T, ushort R, ushort C>
class Matrix;

template<typename T, ushort C>
class MatrixSqr;

namespace MatrixMath {
    template <typename T>
    inline T absoluteValue(const T val) {
        return ((std::numeric_limits<T>::is_signed == true) && ((val) < -std::numeric_limits<T>::epsilon()) ? -(val) : (val));
    }
    
    template <typename T>
    inline bool isGreaterThan(const T val1, const T val2) {
        return (val1) > (val2 + std::numeric_limits<T>::epsilon());
    }
    
    template <typename T>
    inline bool isLesserThan(const T val1, const T val2) {
        return (val1 + std::numeric_limits<T>::epsilon()) < val2;
    }
    
    template <typename T>
    inline T areEqual(const T val1, const T val2) {
        return (isLesserThan(val1, val2) == false) && (isGreaterThan(val1, val2) == false);
    }
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrix<T, R, C> operator+(const AbstractMatrix<T, R, C>& mat, const T value) {
    Matrix<T, R, C> res;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            res[r][c] = mat[r][c] + value;
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrix<T, R, C> operator+(const T value, const AbstractMatrix<T, R, C>& mat) {
    return (mat + value);
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrix<T, R, C> operator+(const AbstractMatrix<T, R, C>& mat1, const AbstractMatrix<T, R, C>& mat2) {
    Matrix<T, R, C> res;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            res[r][c] = mat1[r][c] + mat2[r][c];
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrix<T, R, C> operator-(const AbstractMatrix<T, R, C>& mat1, const AbstractMatrix<T, R, C>& mat2) {
    Matrix<T, R, C> res;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            res[r][c] = mat1[r][c] - mat2[r][c];
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrix<T, R, C> operator-(const AbstractMatrix<T, R, C>& mat, const T value) {
    Matrix<T, R, C> res;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            res[r][c] = mat[r][c] - value;
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrix<T, R, C> operator-(const T value, const AbstractMatrix<T, R, C>& mat) {
    Matrix<T, R, C> res;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            res[r][c] = value - mat[r][c];
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrix<T, R, C> operator*(const AbstractMatrix<T, R, C>& mat, const T value) {
    Matrix<T, R, C> res;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            res[r][c] = mat[r][c] * value;
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrix<T, R, C> operator*(const T value, const AbstractMatrix<T, R, C>& mat) {
    return (mat * value);
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C, ushort OC>
Matrix<T, R, OC> operator*(const AbstractMatrix<T, R, C>& mat1, const AbstractMatrix<T, C, OC>& mat2) {
    // [R, OC] = [R, C] x [C, OC]
    Matrix<T, R, OC> res;
    
    for(ushort l = 0U; l < R; ++l){
        for(ushort oc = 0U; oc < OC; ++oc){
            res[l][oc] = static_cast<T>(0);
            for(ushort c = 0U; c < C; ++c)
                res[l][oc] += mat1[l][c] * mat2[c][oc];
        }
    }
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrix<T, R, C> operator/(const AbstractMatrix<T, R, C>& mat1, const MatrixSqr<T, C>& mat2) {
    return mat1 * mat2.invert();
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrix<T, R, C> operator/(const AbstractMatrix<T, R, C>& mat1, const AbstractMatrix<T, C, C>& mat2) {
    return mat1 / static_cast<MatrixSqr<T,C>>(mat2);
}

#endif /* MatrixMath_h */
