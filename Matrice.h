//
//  Matrice.h
//  test
//
//  Created by Marc-Antoine MARTIN on 09/04/2016.
//  Copyright © 2016 Marc-Antoine MARTIN. All rights reserved.
//

#ifndef MATRICE_H
#define MATRICE_H

#include <initializer_list>
#include <iostream>
#include <utility>

#ifndef ushort
typedef unsigned short ushort;
#endif

namespace MatriceMath {
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
/******************************************** AbstractMatrice *******************************************/
/********************************************************************************************************/
template <typename T, ushort R, ushort C>
class AbstractMatrice {
public:
    AbstractMatrice() = default;
    virtual ~AbstractMatrice() = default;
    
    void operator=(const T values[R][C]);
    void operator=(const std::initializer_list<std::initializer_list<T>> & values);
    void operator=(const AbstractMatrice<T, R, C>& mat);
    
    constexpr ushort rowCount() const;
    constexpr ushort columnCount() const;
    
    virtual T (&operator[](const ushort r)) [C] = 0;
    virtual const T (&operator[](const ushort r) const) [C] = 0;
    
    bool operator==(const AbstractMatrice<T, R, C>& mat2) const;
    
    AbstractMatrice<T, R, C>& operator+=(const AbstractMatrice<T, R, C>& mat2);
    AbstractMatrice<T, R, C>& operator+=(const T value);
    
    AbstractMatrice<T, R, C>& operator-=(const AbstractMatrice<T, R, C>& mat2);
    AbstractMatrice<T, R, C>& operator-=(T value);
    
    T gauss_pivot();
};

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
void AbstractMatrice<T,R,C>::operator=(const T values[R][C]) {
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            operator[](r)[c] = values[r][c];
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
void AbstractMatrice<T,R,C>::operator=(const std::initializer_list<std::initializer_list<T>> & values) {
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
void AbstractMatrice<T,R,C>::operator=(const AbstractMatrice<T, R, C>& mat) {
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            operator[](r)[c] = mat[r][c];
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
constexpr ushort AbstractMatrice<T,R,C>::rowCount() const {
    return R;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
constexpr ushort AbstractMatrice<T,R,C>::columnCount() const {
    return C;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
bool AbstractMatrice<T, R, C>::operator==(const AbstractMatrice<T, R, C>& mat2) const {
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            if(operator[](r)[c] != mat2[r][c])
                return false;
    
    return true;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
AbstractMatrice<T, R, C>& AbstractMatrice<T, R, C>::operator+=(const AbstractMatrice<T, R, C>& mat) {
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            operator[](r)[c] += mat[r][c];
    
    return *this;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
AbstractMatrice<T, R, C>& AbstractMatrice<T, R, C>::operator+=(const T value) {
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            operator[](r)[c] += value;
    
    return *this;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
AbstractMatrice<T, R, C>& AbstractMatrice<T, R, C>::operator-=(const AbstractMatrice<T, R, C>& mat) {
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            operator[](r)[c] -= mat[r][c];
    
    return *this;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
AbstractMatrice<T, R, C>& AbstractMatrice<T, R, C>::operator-=(T value) {
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            operator[](r)[c] -= value;
    
    return *this;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
T AbstractMatrice<T, R, C>::gauss_pivot(){
    // /!\ modifie la matrice out
    
    T det = static_cast<T>(1); //permet le calcul du determinent
    ushort p = 0U; //nombre de permutation, permet le calcul du determient
    
    int r = -1; //r est l'indice de ligne du dernier pivot trouvé
    static const ushort minSide = R < C ? R : C;
    
    for(ushort j = 0U; j < minSide; ++j){ //(j décrit tous les indices de colonnes)
        //Rechercher max(|mat[i][j]|, r+1 <= i < n)
        ushort k = static_cast<T>(r+1); //Noter k l'indice de ligne du maximum
        for(ushort i = r+1; i < R; ++i)
            if(MatriceMath::isGreaterThan(MatriceMath::absoluteValue(operator[](i)[j]), MatriceMath::absoluteValue(operator[](k)[j])) == true)
                k = i;
        
        //mat[k][j] est le pivot
        if(MatriceMath::areEqual(operator[](k)[j], static_cast<T>(0)) == false){
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
                if((i != r) && (MatriceMath::areEqual(tmp, static_cast<T>(0)) == false))
                    for(k = 0U; k < C; ++k)
                        //Soustraire à la ligne i la ligne r multipliée par A[i,j] (de façon à annuler A[i,j])
                        operator[](i)[k] -= operator[](r)[k] * tmp;
            }
        }
    }
    
    return static_cast<T>((p % 2U == 0U) ? 1 : -1) * det;
}

/********************************************************************************************************/
/*********************************************** MatriceRef ************************************************/
/********************************************************************************************************/
template <typename T, ushort R, ushort C>
class MatriceRef: public AbstractMatrice<T, R, C> {
private:
    T (&values)[R][C];
    
public:
    MatriceRef(T (&values)[R][C]);
    
    virtual T (&operator[](const ushort r)) [C] override;
    virtual const T (&operator[](const ushort r) const) [C] override;
};

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
MatriceRef<T,R,C>::MatriceRef(T (&values)[R][C]):
    AbstractMatrice<T, R, C>(),
    values(values)
{
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
T (&MatriceRef<T,R,C>::operator[](const ushort r)) [C] {
    return values[r];
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
const T (&MatriceRef<T,R,C>::operator[](const ushort r) const) [C]  {
    return values[r];
}

/********************************************************************************************************/
/*********************************************** Matrice ************************************************/
/********************************************************************************************************/
template <typename T, ushort R, ushort C>
class Matrice: public AbstractMatrice<T, R, C> {
private:
    T values[R][C];
    
public:
    Matrice() = default;
    Matrice(const T values[R][C]);
    Matrice(const std::initializer_list<std::initializer_list<T>>& values);
    
    template <ushort R2, ushort C2, ushort ROFF=0U, ushort COFF=0U>
    Matrice(const AbstractMatrice<T, R2, C2>& mat);
    
    virtual T(&operator[](const ushort r)) [C] override;
    virtual const T(&operator[](const ushort r) const) [C] override;
    
    Matrice<T, C, R> transpose() const;
    static Matrice<T, R, C> zeros();
    
    template <ushort C2>
    Matrice<T, R, C + C2> increase(const AbstractMatrice<T, R, C2>& mat) const;
    
    template <ushort C2>
    Matrice<T, R, C + C2> operator<<(const AbstractMatrice<T, R, C2>& mat);
    
    template <ushort R2, ushort C2>
    Matrice<T, R2, C2> submatrix(ushort row, ushort col) const;
};

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrice<T,R,C>::Matrice(const T values[R][C]):
    AbstractMatrice<T, R, C>()
{
    AbstractMatrice<T, R, C>::operator=(values);
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrice<T,R,C>::Matrice(const std::initializer_list<std::initializer_list<T>>& values):
    AbstractMatrice<T, R, C>()
{
    AbstractMatrice<T, R, C>::operator=(values);
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
template <ushort R2, ushort C2, ushort ROFF, ushort COFF>
Matrice<T,R,C>::Matrice(const AbstractMatrice<T, R2, C2>& mat):
    AbstractMatrice<T, R, C>()
{
    for(ushort r = ROFF; r < (R < (ROFF + R2) ? R : (ROFF + R2)); ++r)
        for(ushort c = COFF; c < (C < (COFF + C2) ? C : (COFF + C2)); ++c)
            operator[](r)[c] = mat[r - ROFF][c - COFF];
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
T (&Matrice<T,R,C>::operator[](const ushort r)) [C] {
    return values[r];
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
const T (&Matrice<T,R,C>::operator[](const ushort r) const) [C]  {
    return values[r];
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrice<T, C, R> Matrice<T,R,C>::transpose() const {
    Matrice<T, C, R> mat;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            mat[c][r] = operator[](r)[c];
    
    return mat;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrice<T, R, C> Matrice<T, R, C>::zeros() {
    
    Matrice<T, R, C> mout;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            mout[r][c] = static_cast<T>(0);
    
    return mout;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
template <ushort C2>
Matrice<T, R, C + C2> Matrice<T, R, C>::increase(const AbstractMatrice<T, R, C2>& mat) const {
    
    Matrice<T, R, C + C2> mout;
    
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
Matrice<T, R, C + C2> Matrice<T, R, C>::operator<<(const AbstractMatrice<T, R, C2>& mat){
    return increase<C2>(mat);
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
template <ushort R2, ushort C2>
Matrice<T, R2, C2> Matrice<T, R, C>::submatrix(const ushort row, const ushort col) const {
    
    Matrice<T, R2, C2> mout;
    
    for(ushort r = row; r < row + R2 && r < R; ++r)
        for(ushort c = col; c < col + C2 && c < C; ++c)
            mout[r - row][c - col] = operator[](r)[c];
    
    return mout;
}

/********************************************************************************************************/
/****************************************** MatriceSqr **************************************************/
/********************************************************************************************************/
template<typename T, ushort R>
class MatriceSqr : public Matrice<T, R, R>{
private:
    template <ushort R2 = R - 1U, typename DUMMY = void> struct Func {
        static T cofactor(const MatriceSqr<T, R>& mat, ushort ri, ushort ci);
    };
    
    template <typename DUMMY> struct Func<0U, DUMMY> {
        static T cofactor(const MatriceSqr<T, R>& mat, ushort ri, ushort ci);
    };
    
    template <ushort R2 = R - 1>
    MatriceSqr<T, R2> cofactorMat(ushort ri, ushort ci) const;
    
    template <ushort R2 = R - 1U>
    T cofactor(ushort ri, ushort ci) const;
    
public:
    MatriceSqr() = default;
    MatriceSqr(const T values[R][R]):Matrice<T, R, R>(values){}
    MatriceSqr(const std::initializer_list<std::initializer_list<T>> values):Matrice<T, R, R>(values){}
    MatriceSqr(const MatriceSqr<T, R>& mat):Matrice<T, R, R>(mat){}
    
    template <ushort R2, ushort C2, ushort ROFF=0U, ushort COFF=0U>
    MatriceSqr(const AbstractMatrice<T, R2, C2>& mat):Matrice<T, R, R>(mat){}
    
    template <ushort R2 = R - 1U>
    T determinant() const;
    
    MatriceSqr<T, R> invert() const;
    
    template <ushort C>
    Matrice<T, R, C> solve(const AbstractMatrice<T, R, C> &res) const;
    
    MatriceSqr<T, R> invert_gauss() const;
    
    void decomposition_lu(MatriceSqr<T, R>& l, MatriceSqr<T, R>& u) const;
    
    template <ushort C>
    MatriceSqr<T,R> solve_lu(const AbstractMatrice<T,R,C> & res) const;
    MatriceSqr<T,R> invert_lu() const;
    
    static MatriceSqr<T,R> diag(const Matrice<T,R,static_cast<ushort>(1U)> & res);
    static MatriceSqr<T, R> identity();
};

/********************************************************************************************************/
template <typename T, ushort R>
template <ushort R2, typename DUMMY>
T MatriceSqr<T,R>::Func<R2, DUMMY>::cofactor(const MatriceSqr<T, R>& mat, ushort ri, ushort ci) {
    return static_cast<T>(((ri + ci) % 2U == 0U) ? 1 : -1) * mat.cofactorMat(ri, ci).determinant();
}

/********************************************************************************************************/
template <typename T, ushort R>
template <typename DUMMY>
T MatriceSqr<T,R>::Func<0U, DUMMY>::cofactor(const MatriceSqr<T, R>& mat, ushort ri, ushort ci) {
    return static_cast<T>(1);
}

/********************************************************************************************************/
template <typename T, ushort R>
template <ushort R2>
MatriceSqr<T, R2> MatriceSqr<T,R>::cofactorMat(ushort ri, ushort ci) const {
    MatriceSqr<T, R2> res;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < R; ++c)
            if(r != ri && c != ci)
                res[(r > ri) ? (r - 1U) : r][(c > ci) ? (c - 1U) : c] = MatriceSqr<T, R>::operator[](r)[c];
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R>
template <ushort R2>
T MatriceSqr<T,R>::cofactor(ushort ri, ushort ci) const {
    return Func<R2>::cofactor(*this, ri, ci);
}

/********************************************************************************************************/
template <typename T, ushort R>
template <ushort R2>
T MatriceSqr<T,R>::determinant() const {
    T det = static_cast<T>(0);
    
    for(ushort c = 0U; c < R; ++c)
        det += MatriceSqr<T, R>::operator[](0U)[c] * cofactor(0U, c);
    
    return det;
}

/********************************************************************************************************/
template <typename T, ushort R>
MatriceSqr<T, R> MatriceSqr<T,R>::invert() const {
    T detRev = static_cast<T>(1) / determinant();
    
    MatriceSqr<T, R> mat;
    
    //calc commatrice transposed
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < R; ++c)
            mat[c][r] = detRev * cofactor(r, c);
    
    return mat;
}

/********************************************************************************************************/
template <typename T, ushort R>
template <ushort C>
Matrice<T, R, C> MatriceSqr<T,R>::solve(const AbstractMatrice<T, R, C> &mat) const {
    //solve only systems of n equations and n unknowns
    
    //augmente d'une matrice identite de taille eqs[0].length x eqs[0].length
    Matrice<T, R, R + C> tmp = this->increase(mat); //creer une copie TODO: utiliser une MatriceRef
    tmp.gauss_pivot(); //modifie tmp
    
    return tmp.template submatrix<R, C>(0U, R);
}

/********************************************************************************************************/
template <typename T, ushort R>
MatriceSqr<T, R> MatriceSqr<T,R>::invert_gauss() const {
    return MatriceSqr<T, R>::solve<R>(identity());
}

/********************************************************************************************************/
template <typename T, ushort R>
void MatriceSqr<T,R>::decomposition_lu(MatriceSqr<T, R>& l, MatriceSqr<T, R>& u) const {
    //Decomposition en matrice L et U: mat = l*u
    
    u[0U][0U] = MatriceSqr<T, R>::operator[](0U)[0U];
    
    for(ushort j = 1U; j < R; ++j){
        u[0U][j] = MatriceSqr<T, R>::operator[](0U)[j]; // la première ligne de U est la première ligne de mat
        l[j][0U] = MatriceSqr<T, R>::operator[](j)[0U] / MatriceSqr<T, R>::operator[](0U)[0U]; // la première colonne de L est formée des pivots
    }
    
    for(ushort i = 1U; i < R-1U; ++i){
        u[i][i] = MatriceSqr<T, R>::operator[](i)[i];
        
        // on fabrique les coefficients de la diagonale de U sauf u[n][n]
        for(ushort k = 0U; k < i; ++k)
            u[i][i] -= l[i][k] * u[k][i];
        
        
        // on fabrique les coefficients de U et L en dehors de la diagonale
        for(ushort j = i+1U; j < R; ++j){
            u[i][j] = MatriceSqr<T, R>::operator[](i)[j];
            l[j][i] = MatriceSqr<T, R>::operator[](j)[i];
            
            for(ushort k = 0U; k < i; ++k){
                u[i][j] -= l[i][k] * u[k][j];
                l[j][i] -= l[j][k] * u[k][i];
            }
            l[j][i] /= u[i][i];
        }
    }
    
    u[R - 1U][R - 1U] = MatriceSqr<T, R>::operator[](R - 1U)[R - 1U];
    for(ushort k = 0U; k < R-1U; ++k)
        u[R - 1U][R - 1U] -= l[R - 1U][k] * u[k][R - 1U];
}

/********************************************************************************************************/
template <typename T, ushort R>
template <ushort C>
MatriceSqr<T,R> MatriceSqr<T,R>::solve_lu(const AbstractMatrice<T,R,C> & res) const {
    //Base sur: http://www.lycee-pothier.com/LYCEE/pcsi1/file/ipt/TP/tp10/Matrice.pdf
    
    //l est une matrice identite et u une matrice zero
    
    //mat*X=res => resoudre: l*z=res puis u*x=z => x
    MatriceSqr<T,R> u = MatriceSqr<T,R>::zeros();
    MatriceSqr<T,R> l = MatriceSqr<T,R>::identity();
    
    decomposition_lu(l, u);
    
    //substitution avant: l*z=res    (mxm * m*n = m*n)
    Matrice<T,R,C> z = Matrice<T,R,C>::zeros();
    for(ushort k = 0U; k < R; ++k){
        z[k][0U] = res[k][0U];
        for(ushort i = 1U; i < C; ++i){
            z[k][i] = res[k][i];
            for(ushort j = 0U; j < i; ++j)
                z[k][i] -= l[i][j] * z[k][j];
        }
    }
    
    //substitution arriere: u*x=z    (mxm * m*n = m*n)
    Matrice<T,R,C> x = Matrice<T,R,C>::zeros();
    for(ushort k = 0U; k < R; ++k){
        for(short i = C - 1U; i > -1U; --i){
            x[k][i] = z[k][i];
            for(ushort j = i+1U; j < C; ++j)
                x[k][i] -= u[i][j] * x[k][j];
            x[k][i] /= u[i][i];
        }
    }
    
    return x.transpose();
}

/********************************************************************************************************/
template <typename T, ushort R>
MatriceSqr<T,R> MatriceSqr<T,R>::invert_lu() const {
    return solve_lu(identity());
}

/********************************************************************************************************/
template <typename T, ushort R>
MatriceSqr<T,R> MatriceSqr<T,R>::diag(const Matrice<T,R,static_cast<ushort>(1U)> & mat) {
    MatriceSqr<T, R> res;
    
    for(ushort i = 0U; i < R; ++i)
        res[i][i] = mat[i][0U];
    
    return res;
}


/********************************************************************************************************/
template <typename T, ushort R>
MatriceSqr<T, R> MatriceSqr<T,R>::identity() {
    constexpr T ONE = static_cast<T>(1);
    return diag({{ONE},{ONE},{ONE}});
}

/********************************************************************************************************/
/****************************************** OTHER *******************************************************/
/********************************************************************************************************/

/********************************************************************************************************/
template<typename T, ushort R, ushort C>
std::ostream& operator<<(std::ostream& out, const AbstractMatrice<T, R, C>& mat){
    for(ushort r = 0U; r < R; ++r){
        for(ushort c = 0U; c < C; ++c)
            out << mat[r][c] << "\t";
        out << std::endl;
    }
    return out;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrice<T, R, C> operator+(const AbstractMatrice<T, R, C>& mat, const T value) {
    Matrice<T, R, C> res;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            res[r][c] = mat[r][c] + value;
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrice<T, R, C> operator+(const T value, const AbstractMatrice<T, R, C>& mat) {
    return (mat + value);
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrice<T, R, C> operator+(const AbstractMatrice<T, R, C>& mat1, const AbstractMatrice<T, R, C>& mat2) {
    Matrice<T, R, C> res;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            res[r][c] = mat1[r][c] + mat2[r][c];
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrice<T, R, C> operator-(const AbstractMatrice<T, R, C>& mat1, const AbstractMatrice<T, R, C>& mat2) {
    Matrice<T, R, C> res;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            res[r][c] = mat1[r][c] - mat2[r][c];
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrice<T, R, C> operator-(const AbstractMatrice<T, R, C>& mat, const T value) {
    Matrice<T, R, C> res;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            res[r][c] = mat[r][c] - value;
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrice<T, R, C> operator-(const T value, const AbstractMatrice<T, R, C>& mat) {
    Matrice<T, R, C> res;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            res[r][c] = value - mat[r][c];
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrice<T, R, C> operator*(const AbstractMatrice<T, R, C>& mat, const T value) {
    Matrice<T, R, C> res;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < C; ++c)
            res[r][c] = mat[r][c] * value;
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrice<T, R, C> operator*(const T value, const AbstractMatrice<T, R, C>& mat) {
    return (mat * value);
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C, ushort OC>
Matrice<T, R, OC> operator*(const AbstractMatrice<T, R, C>& mat1, const AbstractMatrice<T, C, OC>& mat2) {
    // [R, OC] = [R, C] x [C, OC]
    Matrice<T, R, OC> res;
    
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
Matrice<T, R, C> operator/(const AbstractMatrice<T, R, C>& mat1, const MatriceSqr<T, C>& mat2) {
    return mat1 * mat2.invert();
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
Matrice<T, R, C> operator/(const AbstractMatrice<T, R, C>& mat1, const AbstractMatrice<T, C, C>& mat2) {
    return mat1 / static_cast<MatriceSqr<T,C>>(mat2);
}

#endif /* MATRICE_H */
