//
//  AbstractMatrix.hpp
//
//  Created by Marc-Antoine MARTIN on 07/10/2018.
//  Copyright © 2018 Marc-Antoine MARTIN. All rights reserved.
//

#ifndef AbstractMatrix_h
#define AbstractMatrix_h

#include <initializer_list>
#include <iostream>

template <typename T, ushort R, ushort C>
class AbstractMatrix {
public:
    AbstractMatrix() = default;
    virtual ~AbstractMatrix() = default;
    
    void operator=(const T values[R][C]);
    void operator=(const std::initializer_list<std::initializer_list<T>> & values);
    void operator=(const AbstractMatrix<T, R, C>& mat);
    
    constexpr ushort rowCount() const;
    constexpr ushort columnCount() const;
    
    virtual T (&operator[](const ushort r)) [C] = 0;
    virtual const T (&operator[](const ushort r) const) [C] = 0;
    
    bool operator==(const AbstractMatrix<T, R, C>& mat2) const;
    
    AbstractMatrix<T, R, C>& operator+=(const AbstractMatrix<T, R, C>& mat2);
    AbstractMatrix<T, R, C>& operator+=(const T value);
    
    AbstractMatrix<T, R, C>& operator-=(const AbstractMatrix<T, R, C>& mat2);
    AbstractMatrix<T, R, C>& operator-=(T value);
    
    T gauss_pivot();
};

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
void AbstractMatrix<T,R,C>::operator=(const T values[R][C]) {
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            operator[](r)[c] = values[r][c];
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
void AbstractMatrix<T,R,C>::operator=(const std::initializer_list<std::initializer_list<T>> & values) {
    ushort r = 0U;
    for(const std::initializer_list<T>* row = values.begin(); (row != values.end()) && (r < R); ++row) {
        ushort c = 0U;
        for(const T* vals = row->begin(); (vals != row->end()) && (c < C); ++vals) {
            operator[](r)[c] = *vals;
            c++;
        }
        r++;
    }
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
void AbstractMatrix<T,R,C>::operator=(const AbstractMatrix<T, R, C>& mat) {
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            operator[](r)[c] = mat[r][c];
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
constexpr ushort AbstractMatrix<T,R,C>::rowCount() const {
    return R;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
constexpr ushort AbstractMatrix<T,R,C>::columnCount() const {
    return C;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
bool AbstractMatrix<T, R, C>::operator==(const AbstractMatrix<T, R, C>& mat2) const {
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            if(operator[](r)[c] != mat2[r][c])
                return false;
    
    return true;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
AbstractMatrix<T, R, C>& AbstractMatrix<T, R, C>::operator+=(const AbstractMatrix<T, R, C>& mat) {
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            operator[](r)[c] += mat[r][c];
    
    return *this;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
AbstractMatrix<T, R, C>& AbstractMatrix<T, R, C>::operator+=(const T value) {
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            operator[](r)[c] += value;
    
    return *this;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
AbstractMatrix<T, R, C>& AbstractMatrix<T, R, C>::operator-=(const AbstractMatrix<T, R, C>& mat) {
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            operator[](r)[c] -= mat[r][c];
    
    return *this;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
AbstractMatrix<T, R, C>& AbstractMatrix<T, R, C>::operator-=(T value) {
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            operator[](r)[c] -= value;
    
    return *this;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
T AbstractMatrix<T, R, C>::gauss_pivot(){
    // /!\ modifie la matrice out
    
    T det = static_cast<T>(1); //permet le calcul du determinent
    ushort p = 0U; //nombre de permutation, permet le calcul du determient
    
    int r = -1; //r est l'indice de ligne du dernier pivot trouvé
    static const ushort minSide = R < C ? R : C;
    
    for(ushort j = 0U; j < minSide; ++j){ //(j décrit tous les indices de colonnes)
        //Rechercher max(|mat[i][j]|, r+1 <= i < n)
        ushort k = static_cast<T>(r+1); //Noter k l'indice de ligne du maximum
        for(ushort i = r+1; i < R; ++i)
            if(MatrixMath::isGreaterThan(MatrixMath::absoluteValue(operator[](i)[j]), MatrixMath::absoluteValue(operator[](k)[j])) == true)
                k = i;
        
        //mat[k][j] est le pivot
        if(MatrixMath::areEqual(operator[](k)[j], static_cast<T>(0)) == false){
            ++r;
            det *= operator[](k)[j];
            
            //Diviser la ligne k par mat[k][j] (On normalise la ligne de pivot de façon que le pivot prenne la valeur 1)
            T tmp = operator[](k)[j];
            for(ushort i = 0U; i < C; ++i)
                operator[](k)[i] /= tmp;
            
            //Échanger les lignes k et r (On place la ligne du pivot en position r)
            if(k != r){
                std::swap(operator[](k), operator[](r));
                ++p;
            }
            
            //(On simplifie les autres lignes)
            for(ushort i = 0; i < R; ++i){
                tmp = operator[](i)[j];
                if((i != r) && (MatrixMath::areEqual(tmp, static_cast<T>(0)) == false))
                    for(k = 0U; k < C; ++k)
                        //Soustraire à la ligne i la ligne r multipliée par A[i,j] (de façon à annuler A[i,j])
                        operator[](i)[k] -= operator[](r)[k] * tmp;
            }
        }
    }
    
    return static_cast<T>((p % 2U == 0U) ? 1 : -1) * det;
}

/********************************************************************************************************/
template<typename T, ushort R, ushort C>
std::ostream& operator<<(std::ostream& out, const AbstractMatrix<T, R, C>& mat){
    for(ushort r = 0U; r < R; ++r){
        for(ushort c = 0U; c < C; ++c)
            out << mat[r][c] << "\t";
        out << std::endl;
    }
    return out;
}

#endif /* AbstractMatrix_h */
