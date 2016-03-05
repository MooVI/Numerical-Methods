/* 
 * File:   MonteCarloIntegration.h
 * Author: qfeuille
 *
 * Created on 11 September 2012, 15:36
 */

#include <Eigen/Dense>
#include "mtrand.h" 
#ifndef MONTECARLOINTEGRATION_H
#define	MONTECARLOINTEGRATION_H
namespace NumMethod {
   
#include <Cuba/cuba.h>

template<typename T,int N>
struct SimpleMCIntegration{
    template<typename FunctPtr>
    T integrate (FunctPtr f, const Eigen::Matrix<T,N+1,1>& lower,const  Eigen::Matrix<T,N+1,1>& upper, int nsteps){
        MTRand * drand = new MTRand ( time (0));
        T integral=0;
        Eigen::Matrix<T,N+1,1> values; Eigen::Matrix<T,N+1,1> range = upper -lower;
        for (int i=0;i<nsteps;i++){
            for (int j=0;j<N+1;j++)
                values[j]= ((*drand)())*range [j];
            auto output = f(values.template head<N> ());
            if (values[N]<output)
                integral+=1;
        }
        return (integral/nsteps)*range.prod();
        
    }
};

}

#endif	/* MONTECARLOINTEGRATION_H */

