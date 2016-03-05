/* 
 * File:   DEIntegrator.h
 * Author: qfeuille
 *
 * Created on 05 May 2013, 00:33
 */

#ifndef DEINTEGRATOR_H
#define	DEINTEGRATOR_H

#include "DESolversInc.h" 
namespace NumMethod {

 template <typename T>
    class TrapeziumDE {
        Eigen::Matrix <T, N, 1 > k;
    public:

        template<typename FunctPtr>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::Matrix<T, N, 1 > & q, const Eigen::Matrix<T, N, 1 > & initialdq, const T& stepSize) {

            q += stepSize * initialdq;
        };
        

    };
}

#endif	/* DEINTEGRATOR_H */

