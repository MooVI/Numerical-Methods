/* 
 * File:   Differentiation.h
 * Author: qfeuille
 *
 * Created on 24 August 2012, 01:53
 */

#ifndef DIFFERENTIATION_H
#define	DIFFERENTIATION_H
#include <Eigen/Dense>
#include "StdVectorMath.h"
namespace NumMethod{
    
struct DiffByDef { //Really not advised, unless f is really expensive
    template<typename FunctPtr, typename T>
    T operator () (FunctPtr f,const T x,const T h){
        return (f(x+h)-f(x))/h;
    }
    template<typename FunctPtr, typename T>
    T operator () (FunctPtr f,const T x, const T h, const T fx){ //When you already have f(x)
        return (f(x+h)-fx)/h;
    }
    double advisedStep (double x){//Only for double accuracy, assuming error 
         return sqrt(1.11e-16)*x; // in f is ~ machine precision
    }
};

struct DiffSymmetrical { //OK
    template<typename FunctPtr, typename T>
    T operator () (FunctPtr f,const T x,const T h){
        return (f(x+h)-f(x-h))/(2.0*h);
    }
    double advisedStep (double x){ //Only for double accuracy, assuming error 
        return pow(1.11e-16,0.333333)*x; // in f is ~ machine precision
    }
};

template<typename T, int TabSize=10>
class Ridders { //By far the best
    const T shrink,shrinksq;
    const T large=1e100;
    const T safety;
    Eigen::Matrix<T,TabSize,TabSize> tab;
public:
    template<typename FunctPtr>
    T operator () (FunctPtr f,const T x,T h, T& error){
        if (h==0) return 0;
        T ret;
        tab (0,0)=(f(x+h)-f(x-h))/(2.0*h);
        error=large;
        for (int i=1;i<TabSize;i++){
            h /= shrink;
            tab (0,i) = (f(x+h)-f(x-h))/(2.0*h);
            T fac = shrinksq;
            for (int j=1;j<=i;j++){
             tab (j,i) =( tab (j-1,i)*fac - tab(j-1,i-1)) /(fac-1.0);
             fac*=shrinksq;
             T newerror = maxAbs (tab(j,i)-tab(j-1,i),tab(j,i)-tab(j-1,i-1));
             if (newerror<= error)
             {error =newerror;
             ret = tab(j,i);
             }
            }
            if (fabs(tab(i,i)-tab(i-1,i-1))>= safety*error) break;
        }
        return ret;
    }
    Ridders (T ishrink=1.4 ,T isafety=2.0): shrink (ishrink),safety(isafety), shrinksq(ishrink*ishrink){
        
                
    }
};

struct NumerovDifferentiate{
    enum {centre = 2};
    const double outc,a1c,a2c,b1c,b2c; 
    template<typename FunctPtr, typename B, typename T>
    T operator ()(FunctPtr f,const Eigen::ArrayBase<B>& dvalues,const T& x,const T& h){
        double hsqtwelve = h*h/ 12.0;
        double a1 =0.5*( dvalues [centre+1]- dvalues[centre-1] );
        double a2 =0.5*( dvalues [centre+2]- dvalues[centre-2] );
        double b1 =0.5*hsqtwelve*( f(x+h)* dvalues [centre+1]- f(x-h)*dvalues[centre-1] );
        double b2 =0.5*hsqtwelve*( f(x+2*h)*dvalues [centre+2]- f(x-2*h)*dvalues[centre-2]);
        return outc * (a1c*a1+b1c*b1+a2c*a2+b2c*b2)/ h;
    }
    NumerovDifferentiate ():a1c(-1),a2c(37.0/32.0),b1c(-37.0/5.0),b2c(-17.0/40.0), outc(16.0/21.0){}
    
};



};
#endif	/* DIFFERENTIATION_H */

