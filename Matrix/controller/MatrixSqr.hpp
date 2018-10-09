//
//  MatrixSqr.hpp
//
//  Created by Marc-Antoine MARTIN on 07/10/2018.
//  Copyright © 2018 Marc-Antoine MARTIN. All rights reserved.
//

#ifndef MatrixSqr_h
#define MatrixSqr_h

#include "Matrix.hpp"

template<typename T, ushort R>
class MatrixSqr : public Matrix<T, R, R>{
private:
    template <ushort R2 = R - 1U, typename DUMMY = void> struct Func {
        static T cofactor(const MatrixSqr<T, R>& mat, ushort ri, ushort ci);
    };
    
    template <typename DUMMY> struct Func<0U, DUMMY> {
        static T cofactor(const MatrixSqr<T, R>& mat, ushort ri, ushort ci);
    };
    
    template <ushort R2 = R - 1>
    MatrixSqr<T, R2> cofactorMat(ushort ri, ushort ci) const;
    
    template <ushort R2 = R - 1U>
    T cofactor(ushort ri, ushort ci) const;
    
public:
    MatrixSqr() = default;
    MatrixSqr(const T values[R][R]):Matrix<T, R, R>(values){}
    MatrixSqr(const std::initializer_list<std::initializer_list<T>> values):Matrix<T, R, R>(values){}
    MatrixSqr(const MatrixSqr<T, R>& mat):Matrix<T, R, R>(mat){}
    
    template <ushort R2, ushort C2, ushort ROFF=0U, ushort COFF=0U>
    MatrixSqr(const AbstractMatrix<T, R2, C2>& mat):Matrix<T, R, R>(mat){}
    
    template <ushort R2 = R - 1U>
    T determinant() const;
    
    MatrixSqr<T, R> invert() const;
    
    template <ushort C>
    Matrix<T, R, C> solve(const AbstractMatrix<T, R, C> &res) const;
    
    MatrixSqr<T, R> invert_gauss() const;
    
    void decomposition_lu(MatrixSqr<T, R>& l, MatrixSqr<T, R>& u) const;
    
    template <ushort C>
    MatrixSqr<T,R> solve_lu(const AbstractMatrix<T,R,C> & res) const;
    MatrixSqr<T,R> invert_lu() const;
    
    static MatrixSqr<T,R> diag(const Matrix<T,R,static_cast<ushort>(1U)> & res);
    static MatrixSqr<T, R> identity();
};

/********************************************************************************************************/
template <typename T, ushort R>
template <ushort R2, typename DUMMY>
T MatrixSqr<T,R>::Func<R2, DUMMY>::cofactor(const MatrixSqr<T, R>& mat, ushort ri, ushort ci) {
    return static_cast<T>(((ri + ci) % 2U == 0U) ? 1 : -1) * mat.cofactorMat(ri, ci).determinant();
}

/********************************************************************************************************/
template <typename T, ushort R>
template <typename DUMMY>
T MatrixSqr<T,R>::Func<0U, DUMMY>::cofactor(const MatrixSqr<T, R>& mat, ushort ri, ushort ci) {
    return static_cast<T>(1);
}

/********************************************************************************************************/
template <typename T, ushort R>
template <ushort R2>
MatrixSqr<T, R2> MatrixSqr<T,R>::cofactorMat(ushort ri, ushort ci) const {
    MatrixSqr<T, R2> res;
    
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < R; ++c)
            if(r != ri && c != ci)
                res[(r > ri) ? (r - 1U) : r][(c > ci) ? (c - 1U) : c] = MatrixSqr<T, R>::operator[](r)[c];
    
    return res;
}

/********************************************************************************************************/
template <typename T, ushort R>
template <ushort R2>
T MatrixSqr<T,R>::cofactor(ushort ri, ushort ci) const {
    return Func<R2>::cofactor(*this, ri, ci);
}

/********************************************************************************************************/
template <typename T, ushort R>
template <ushort R2>
T MatrixSqr<T,R>::determinant() const {
    T det = static_cast<T>(0);
    
    for(ushort c = 0U; c < R; ++c)
        det += MatrixSqr<T, R>::operator[](0U)[c] * cofactor(0U, c);
    
    return det;
}

/********************************************************************************************************/
template <typename T, ushort R>
MatrixSqr<T, R> MatrixSqr<T,R>::invert() const {
    T detRev = static_cast<T>(1) / determinant();
    
    MatrixSqr<T, R> mat;
    
    //calc commatrice transposed
    for(ushort r = 0U; r < R; ++r)
        for(ushort c = 0U; c < R; ++c)
            mat[c][r] = detRev * cofactor(r, c);
    
    return mat;
}

/********************************************************************************************************/
template <typename T, ushort R>
template <ushort C>
Matrix<T, R, C> MatrixSqr<T,R>::solve(const AbstractMatrix<T, R, C> &mat) const {
    //solve only systems of n equations and n unknowns
    
    //augmente d'une matrice identite de taille eqs[0].length x eqs[0].length
    Matrix<T, R, R + C> tmp = this->increase(mat); //creer une copie TODO: utiliser une MatrixRef
    tmp.gauss_pivot(); //modifie tmp
    
    return tmp.template submatrix<R, C>(0U, R);
}

/********************************************************************************************************/
template <typename T, ushort R>
MatrixSqr<T, R> MatrixSqr<T,R>::invert_gauss() const {
    return MatrixSqr<T, R>::solve<R>(identity());
}

/********************************************************************************************************/
template <typename T, ushort R>
void MatrixSqr<T,R>::decomposition_lu(MatrixSqr<T, R>& l, MatrixSqr<T, R>& u) const {
    //Decomposition en matrice L et U: mat = l*u
    
    u[0U][0U] = MatrixSqr<T, R>::operator[](0U)[0U];
    
    for(ushort j = 1U; j < R; ++j){
        u[0U][j] = MatrixSqr<T, R>::operator[](0U)[j]; // la première ligne de U est la première ligne de mat
        l[j][0U] = MatrixSqr<T, R>::operator[](j)[0U] / MatrixSqr<T, R>::operator[](0U)[0U]; // la première colonne de L est formée des pivots
    }
    
    for(ushort i = 1U; i < R-1U; ++i){
        u[i][i] = MatrixSqr<T, R>::operator[](i)[i];
        
        // on fabrique les coefficients de la diagonale de U sauf u[n][n]
        for(ushort k = 0U; k < i; ++k)
            u[i][i] -= l[i][k] * u[k][i];
        
        
        // on fabrique les coefficients de U et L en dehors de la diagonale
        for(ushort j = i+1U; j < R; ++j){
            u[i][j] = MatrixSqr<T, R>::operator[](i)[j];
            l[j][i] = MatrixSqr<T, R>::operator[](j)[i];
            
            for(ushort k = 0U; k < i; ++k){
                u[i][j] -= l[i][k] * u[k][j];
                l[j][i] -= l[j][k] * u[k][i];
            }
            l[j][i] /= u[i][i];
        }
    }
    
    u[R - 1U][R - 1U] = MatrixSqr<T, R>::operator[](R - 1U)[R - 1U];
    for(ushort k = 0U; k < R-1U; ++k)
        u[R - 1U][R - 1U] -= l[R - 1U][k] * u[k][R - 1U];
}

/********************************************************************************************************/
template <typename T, ushort R>
template <ushort C>
MatrixSqr<T,R> MatrixSqr<T,R>::solve_lu(const AbstractMatrix<T,R,C> & res) const {
    //Base sur: http://www.lycee-pothier.com/LYCEE/pcsi1/file/ipt/TP/tp10/Matrix.pdf
    
    //l est une matrice identite et u une matrice zero
    
    //mat*X=res => resoudre: l*z=res puis u*x=z => x
    MatrixSqr<T,R> u = MatrixSqr<T,R>::zeros();
    MatrixSqr<T,R> l = MatrixSqr<T,R>::identity();
    
    decomposition_lu(l, u);
    
    //substitution avant: l*z=res    (mxm * m*n = m*n)
    Matrix<T,R,C> z = Matrix<T,R,C>::zeros();
    for(ushort k = 0U; k < R; ++k){
        z[k][0U] = res[k][0U];
        for(ushort i = 1U; i < C; ++i){
            z[k][i] = res[k][i];
            for(ushort j = 0U; j < i; ++j)
                z[k][i] -= l[i][j] * z[k][j];
        }
    }
    
    //substitution arriere: u*x=z    (mxm * m*n = m*n)
    Matrix<T,R,C> x = Matrix<T,R,C>::zeros();
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
MatrixSqr<T,R> MatrixSqr<T,R>::invert_lu() const {
    return solve_lu(identity());
}

/********************************************************************************************************/
template <typename T, ushort R>
MatrixSqr<T,R> MatrixSqr<T,R>::diag(const Matrix<T,R,static_cast<ushort>(1U)> & mat) {
    MatrixSqr<T, R> res;
    
    for(ushort i = 0U; i < R; ++i)
        res[i][i] = mat[i][0U];
    
    return res;
}


/********************************************************************************************************/
template <typename T, ushort R>
MatrixSqr<T, R> MatrixSqr<T,R>::identity() {
    constexpr T ONE = static_cast<T>(1);
    return diag({{ONE},{ONE},{ONE}});
}

#endif /* MatrixSqr_h */
