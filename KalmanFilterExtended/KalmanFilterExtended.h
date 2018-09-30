//
//  KalmanFilterExtended.h
//  testCpp
//
//  Created by Marc-Antoine MARTIN on 10/09/2017.
//  Copyright © 2017 Marc-Antoine MARTIN. All rights reserved.
//

#ifndef KalmanFilterExtended_h
#define KalmanFilterExtended_h

#include "Matrice.h"

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL=1U>
class KalmanFilterExtended {
protected:
    //System-state:
    //	X(k+1)	= f(X(k))	//equation d'etat
    //	Y		= h(X(k),B)	//equation de mesure
    
    const MatriceSqr<T, NB_STATES>& A; //matrice jacobienne de transition (nb_etats x nb_etats)
    const Matrice<T, NB_MESURES, NB_STATES>& H; //matrice jacobienne d'observation (nb_mesures x nb_etats)
    const MatriceSqr<T, NB_STATES>& Q; //matrice de covariance des erreurs de modelisation d'etat (nb_etats x nb_etats)
    const MatriceSqr<T, NB_MESURES>& R; //matrice de covariance des bruits de mesures (nb_mesures x nb_mesures)
    
    const MatriceSqr<T, NB_STATES> At; //accelere les calculs
    const Matrice<T, NB_STATES, NB_MESURES> Ht; //accelere les calculs
    
    Matrice<T, NB_STATES, NB_ST_MES_COL> X; //matrice d'etats (nb_etats x NB_ST_MES_COL)
    MatriceSqr<T, NB_STATES> P; //matrice de prediction (nb_etats x nb_etats)
    
    virtual Matrice<T, NB_STATES, NB_ST_MES_COL> a() = 0; //fonction de generation de la matrice d'etats  (nb_etats x N) avec la matrice de transition (nb_etats x nb_etats)
    virtual Matrice<T, NB_MESURES, NB_ST_MES_COL> h() = 0; //fonction de generation de la matrice de mesure (nb_mesures x N) avec la matrice d'observation (nb_mesures x nb_etats)
    
    virtual void updateA() = 0;
    virtual void updateH() = 0;
    
public:
    KalmanFilterExtended(Matrice<T, NB_STATES, NB_ST_MES_COL>& X,
                         const MatriceSqr<T, NB_STATES>& A,
                         const Matrice<T, NB_MESURES, NB_STATES>& H,
                         const MatriceSqr<T, NB_STATES>& Q,
                         const MatriceSqr<T, NB_MESURES>& R
                         );
    virtual ~KalmanFilterExtended() = default;
    
    void predict(void);
    void update(const AbstractMatrice<T, NB_MESURES, NB_ST_MES_COL>& Y);
    
    inline void predict_update(const AbstractMatrice<T, NB_MESURES, NB_ST_MES_COL>& Y);
    
    template <ushort ROW, ushort COL=0U>
    inline T getState() const;
};

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL>
KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::KalmanFilterExtended(Matrice<T, NB_STATES, NB_ST_MES_COL>& X,
                                                                                    const MatriceSqr<T, NB_STATES>& A,
                                                                                    const Matrice<T, NB_MESURES, NB_STATES>& H,
                                                                                    const MatriceSqr<T, NB_STATES>& Q,
                                                                                    const MatriceSqr<T, NB_MESURES>& R
                                                                                    ):A(A),H(H),Q(Q),R(R),At(A.transpose()),Ht(H.transpose()),X(X),P(MatriceSqr<T, NB_STATES>::zeros())
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
void KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::update(const AbstractMatrice<T, NB_MESURES, NB_ST_MES_COL>& Y)
{
    //Y est la matrice des nouvelles mesures
    
    updateH();
    
    //K = P * H' * inv(R + H * P * H')
    MatriceSqr<T, NB_MESURES> tmp(((H * P) * Ht));
    tmp += R;
    const Matrice<T, NB_STATES, NB_MESURES> K((P * Ht) * tmp.invert_gauss());
    //Matrice<T, NB_STATES, NB_MESURES> K((P * Ht) / (R + ((H * P) * Ht)));
    
    //P = P - K * H * P
    P -= (K * H) * P;
    
    //X = X + K * (Y - h(X,B))
    X += K * (Y - h());
}

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL>
void KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::predict_update(const AbstractMatrice<T, NB_MESURES, NB_ST_MES_COL>& Y)
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

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL=1U>
class KalmanFilter: public KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL> {
private:
    //System-state:
    //	X(k+1)	= A*X(k)	//equation d'etat
    //	Y		= H*X(k)+B	//equation de mesure
    
    inline virtual Matrice<T, NB_STATES, NB_ST_MES_COL> a() override;
    inline virtual Matrice<T, NB_MESURES, NB_ST_MES_COL> h() override;
    
    inline virtual void updateA() override {} //do nothing
    inline virtual void updateH() override {}; //do nothing
    
public:
    
    KalmanFilter(Matrice<T, NB_STATES, NB_ST_MES_COL>& X,
                 const MatriceSqr<T, NB_STATES>& A,
                 const Matrice<T, NB_MESURES, NB_STATES>& H,
                 const MatriceSqr<T, NB_STATES>& Q,
                 const MatriceSqr<T, NB_MESURES>& R
                 );
    
    virtual ~KalmanFilter() override = default;
};

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL>
KalmanFilter<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::KalmanFilter(Matrice<T, NB_STATES, NB_ST_MES_COL>& X,
                                                                    const MatriceSqr<T, NB_STATES>& A,
                                                                    const Matrice<T, NB_MESURES, NB_STATES>& H,
                                                                    const MatriceSqr<T, NB_STATES>& Q,
                                                                    const MatriceSqr<T, NB_MESURES>& R
                                                                    ):KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>(X,A,H,Q,R)
{
};

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL>
Matrice<T, NB_STATES, NB_ST_MES_COL> KalmanFilter<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::a()
{
    return KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::A * KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::X;
}

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL>
Matrice<T, NB_MESURES, NB_ST_MES_COL> KalmanFilter<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::h()
{
    return KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::H * KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::X;
}

#endif /* KalmanFilterExtended_h */
