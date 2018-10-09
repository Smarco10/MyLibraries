//
//  KalmanFilter.hpp
//
//  Created by Marc-Antoine MARTIN on 07/10/2018.
//  Copyright © 2018 Marc-Antoine MARTIN. All rights reserved.
//

#ifndef KalmanFilter_h
#define KalmanFilter_h

#include "KalmanFilterExtended.hpp"

template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL=1U>
class KalmanFilter: public KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL> {
private:
    //System-state:
    //    X(k+1)    = A*X(k)    //equation d'etat
    //    Y        = H*X(k)+B    //equation de mesure
    
    inline virtual Matrix<T, NB_STATES, NB_ST_MES_COL> a() override;
    inline virtual Matrix<T, NB_MESURES, NB_ST_MES_COL> h() override;
    
    inline virtual void updateA() override {} //do nothing
    inline virtual void updateH() override {}; //do nothing
    
public:
    
    KalmanFilter(Matrix<T, NB_STATES, NB_ST_MES_COL>& X,
                 const MatrixSqr<T, NB_STATES>& A,
                 const Matrix<T, NB_MESURES, NB_STATES>& H,
                 const MatrixSqr<T, NB_STATES>& Q,
                 const MatrixSqr<T, NB_MESURES>& R
                 );
    
    virtual ~KalmanFilter() override = default;
};

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL>
KalmanFilter<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::KalmanFilter(Matrix<T, NB_STATES, NB_ST_MES_COL>& X,
                                                                    const MatrixSqr<T, NB_STATES>& A,
                                                                    const Matrix<T, NB_MESURES, NB_STATES>& H,
                                                                    const MatrixSqr<T, NB_STATES>& Q,
                                                                    const MatrixSqr<T, NB_MESURES>& R
                                                                    ):KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>(X,A,H,Q,R)
{
};

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL>
Matrix<T, NB_STATES, NB_ST_MES_COL> KalmanFilter<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::a()
{
    return KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::A * KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::X;
}

/********************************************************************************************************/
template <typename T, ushort NB_STATES, ushort NB_MESURES, ushort NB_ST_MES_COL>
Matrix<T, NB_MESURES, NB_ST_MES_COL> KalmanFilter<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::h()
{
    return KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::H * KalmanFilterExtended<T, NB_STATES, NB_MESURES, NB_ST_MES_COL>::X;
}

#endif /* KalmanFilter_h */
