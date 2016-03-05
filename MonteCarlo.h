/* 
 * File:   MonteCarlo.h
 * Author: qfeuille
 *
 * Created on 11 September 2012, 15:31
 */
#include "mtrand.h"

#ifndef MONTECARLO_H
#define	MONTECARLO_H



namespace NumMethod {

 
    
struct SimpleMC {
    MTRand * drand;
    template <typename FunctPtr>
    void operator () (FunctPtr f, int N){

    for (int i=0;i<N;i++){
        if (!f((*drand)())) break;
    }
    }
    SimpleMC (){
     drand = new MTRand (time (0));
    }
    ~SimpleMC(){
        delete drand;
    }
};



}

#endif	/* MONTECARLO_H */

