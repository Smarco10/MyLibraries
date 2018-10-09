//
//  KalmanFilterExtended.h
//
//  Created by Marc-Antoine MARTIN on 10/09/2017.
//  Copyright © 2017 Marc-Antoine MARTIN. All rights reserved.
//

#ifndef KalmanFilterExtended_h
#define KalmanFilterExtended_h

#include "matrix"

#ifndef ushort
typedef unsigned short ushort;
#endif

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL=1U>
class KalmanFilterExtended {
protected:
    //System-state:
    //	X(k+1)	= f(X(k))	//equation d'etat
    //	Y		= h(X(k),B)	//equation de mesure
    
    const MatrixSqr<T, NB_STATES>& A; //matrice jacobienne de transition (nb_etats x nb_etats)
    const Matrix<T, NB_MESURES, NB_STATES>& H; //matrice jacobienne d'observation (nb_mesures x nb_etats)
    const MatrixSqr<T, NB_STATES>& Q; //matrice de covariance des erreurs de modelisation d'etat (nb_etats x nb_etats)
    const MatrixSqr<T, NB_MESURES>& R; //matrice de covariance des bruits de mesures (nb_mesures x nb_mesures)
    
    const MatrixSqr<T, NB_STATES> At; //accelere les calculs
    const Matrix<T, NB_STATES, NB_MESURES> Ht; //accelere les calculs
    
    Matrix<T, NB_STATES, NB_ST_MES_COL> X; //matrice d'etats (nb_etats x NB_ST_MES_COL)
    MatrixSqr<T, NB_STATES> P; //matrice de prediction (nb_etats x nb_etats)
    
    virtual Matrix<T, NB_STATES, NB_ST_MES_COL> a() = 0; //fonction de generation de la matrice d'etats  (nb_etats x N) avec la matrice de transition (nb_etats x nb_etats)
    virtual Matrix<T, NB_MESURES, NB_ST_MES_COL> h() = 0; //fonction de generation de la matrice de mesure (nb_mesures x N) avec la matrice d'observation (nb_mesures x nb_etats)
    
    virtual void updateA() = 0;
    virtual void updateH() = 0;
    
public:
    KalmanFilterExtended(Matrix<T, NB_STATES, NB_ST_MES_COL>& X,
                         const MatrixSqr<T, NB_STATES>& A,
                         const Matrix<T, NB_MESURES, NB_STATES>& H,
                         const MatrixSqr<T, NB_STATES>& Q,
                         const MatrixSqr<T, NB_MESURES>& R
                         );
    virtual ~KalmanFilterExtended() = default;
    
    void predict(void);
    void update(const AbstractMatrix<T, NB_MESURES, NB_ST_MES_COL>& Y);
    
    inline void predict_update(const AbstractMatrix<T, NB_MESURES, NB_ST_MES_COL>& Y);
    
    template <ushort ROW, ushort COL=0U>
    inline T getState() const;
};

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL>
KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::KalmanFilterExtended(Matrix<T, NB_STATES, NB_ST_MES_COL>& X,
                                                                                    const MatrixSqr<T, NB_STATES>& A,
                                                                                    const Matrix<T, NB_MESURES, NB_STATES>& H,
                                                                                    const MatrixSqr<T, NB_STATES>& Q,
                                                                                    const MatrixSqr<T, NB_MESURES>& R
                                                                                    ):A(A),H(H),Q(Q),R(R),At(A.transpose()),Ht(H.transpose()),X(X),P(MatrixSqr<T, NB_STATES>::zeros())
{
};

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL>
void KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::predict(void)
{
    updateA();
    
    //X = f(X)
    X = a();
    
    //P = A * P * A' + Q
    P = (A * P) * At;
    P += Q; //optimisation (do not create a new temporary matrice)
}

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL>
void KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::update(const AbstractMatrix<T, NB_MESURES, NB_ST_MES_COL>& Y)
{
    //Y est la matrice des nouvelles mesures
    
    updateH();
    
    //K = P * H' * inv(R + H * P * H')
    MatrixSqr<T, NB_MESURES> tmp(((H * P) * Ht));
    tmp += R;
    const Matrix<T, NB_STATES, NB_MESURES> K((P * Ht) * tmp.invert_gauss());
    //Matrix<T, NB_STATES, NB_MESURES> K((P * Ht) / (R + ((H * P) * Ht)));
    
    //P = P - K * H * P
    P -= (K * H) * P;
    
    //X = X + K * (Y - h(X,B))
    X += K * (Y - h());
}

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL>
void KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::predict_update(const AbstractMatrix<T, NB_MESURES, NB_ST_MES_COL>& Y)
{
    predict();
    update(Y);
}

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL>
template <ushort ROW, ushort COL>
T KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::getState() const
{
    return X[ROW][COL];
}

#endif /* KalmanFilterExtended_h */
