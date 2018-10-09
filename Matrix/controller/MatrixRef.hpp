//
//  MatrixRef.hPP
//
//  Created by Marc-Antoine MARTIN on 09/10/2018.
//  Copyright © 2018 Marc-Antoine MARTIN. All rights reserved.
//

#ifndef MatrixRef_h
#define MatrixRef_h

#include "AbstractMatrix.hpp"

template <typename T, ushort R, ushort C>
class MatrixRef: public AbstractMatrix<T, R, C> {
private:
    T (&values)[R][C];
    
public:
    MatrixRef(T (&values)[R][C]);
    
    virtual T (&operator[](const ushort r)) [C] override;
    virtual const T (&operator[](const ushort r) const) [C] override;
};

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
MatrixRef<T,R,C>::MatrixRef(T (&values)[R][C]):
    AbstractMatrix<T, R, C>(),
    values(values)
{
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
T (&MatrixRef<T,R,C>::operator[](const ushort r)) [C] {
    return values[r];
}

/********************************************************************************************************/
template <typename T, ushort R, ushort C>
const T (&MatrixRef<T,R,C>::operator[](const ushort r) const) [C]  {
    return values[r];
}

#endif /* MatrixRef_h */
