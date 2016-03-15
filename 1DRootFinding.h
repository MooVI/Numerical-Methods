/* 
 * File:   1DRootFinding.h
 * Author: Jack Kemp
 *
 * Created on 19 October 2012, 12:36
 */

#include <math.h>
#include <limits>
#include "Useful.h"

#include <iostream>
#ifndef ONEDROOTFINDING_H
#define ONEDROOTFINDING_H
namespace NumMethod{
 
    struct Bisection {
        
        
        template <typename FunctPtr, typename T, int maxStep = 40>
         T operator () (FunctPtr& f, const T& lower,const T& upper,const T& desAcc){
            T value = f (lower);
            T midvalue = f(upper);
            T step, ret, mid;
            if (value*midvalue >= 0.0) {
                std::cerr<<"Bisection: Root not bracketed!"<<std::endl;
                return T(0);
            }
            ret = value < 0.0 ? ((step=upper-lower),lower)
                                 : (step=(lower - upper),upper);
            
            for (int i=0;i<maxStep;i++){
                midvalue = f(mid = ret+ (step*=0.5));
                if (midvalue<=0.0) ret = mid;
                if (fabs(step)<desAcc||midvalue==0.0) return ret;
            }
            std::cerr<<"Bisection: Maximum number of iterations reached."<<std::endl;
            return T(0);
        }
                
  };
           //Find a root between upper and lower using Brent's method 
        struct Brent {
            template <typename FunctPtr, typename T, int maxStep = 100>
        T operator () (FunctPtr& f, const T& lower,const T& upper,const T& desAcc){
                const T eps = std::numeric_limits<T>::epsilon();
                T a = lower, b =upper, c= upper, d,e,min1,min2;
                T fa = f(a), fb = f(b), fc,p,q,r,s,tol1,xm;
                
                if ((fa>0.0&&fb>0.0)||(fa<0.0&&fb<0.0)){
                    std::cerr<<"Brent: Root not bracketed!"<<std::endl;
                }
                fc=fb;
                for (int i=0;i<maxStep;i++){
                    if ((fb>0.0 &&fc>0.0)||fb <0.0 &&fc<0.0)
                    {
                        //Root is not between c and b, so move interval
                        c=a;
                        fc=fa;
                        e=d=b-a;
                    }
                    if (fabs(fc)<fabs(fb)){
                        //Relabel (a,b,c) -> (c,b,a)
                        a=b;
                        b=c;
                        c=a;
                        fa=fb;
                        fb=fc;
                        fc=fa;   
                    }
                    //Test convergence, adding a machine precision safety factor
                    tol1 = 2.0*eps*fabs(b)+0.5*desAcc;
                    xm = 0.5*(c-b);
                    if (fabs(xm)<=tol1 || fb==0.0) return b;
                    if (fabs(e)>=tol1 &&fabs(a)>fabs(b)){
                        //Attempt inverse quadratic interpolation
                        s = fb/fa;
                        if (a==c){
                            p = 2.0*xm*s;
                            q =1.0-s; 
                        }
                        else {
                            q = fa/fc;
                            r = fb/fc;
                            p = s*(2.0*xm*q*(q-r)-(b-a)*(r-1.0));
                            q = (q-1.0)*(r-1.0)*(s-1.0);   
                        }
                        if (p>0.0) q = -q;
                        p = fabs(p);
                        min1 = 3.0*xm*q-fabs(tol1*q);
                        min2 = fabs(e*q);
                        if (2.0*p<min(min1,min2)){
                            //Accept interpolation
                            e=d;
                            d=p/q;
                        }
                        else {
                            //Reject interpolation, bisect
                            d=xm;
                            e=d;
                        }     
                    }
                    else {
                         //Interpolation too slow, bisect
                            d=xm;
                            e=d;
                    }
                    a=b;
                    fa=fb;
                    if (fabs(d)>tol1)
                        b+=d;
                    else 
                        b+= sign (tol1,xm);
                    fb = f(b);   
                    
                }
                std::cerr<<"Brent: Maximum iterations reached!"<<std::endl;
                return NAN;     
            }
            
            
        };    

    
        //Safe secant applies secant method but checks to make sure step not too large
    struct SafeSecant {
        template <typename FunctPtr, typename T, int maxStep = 100 >
                T operator () (FunctPtr&& f, const T& guess, const T& initialStep, const T& desAcc, const T& maxRelInc=1000){
            T x1 = guess;
            T x2 = guess + initialStep;
            T y1 = f(x1);
            T y2 = f(x2);
            //Orient so that y2 is closer to root.
            if (fabs(y2) > fabs(y1)) {
                std::swap(y1, y2);
                std::swap(x1, x2);
            }
            int i = 0;
            T dx = initialStep;
            while (i < maxStep) {
                T newstep = (x1 - x2) * y2 / (y2 - y1);
                dx = sign(newstep, dx);
                dx = minAbs(newstep, maxRelInc*dx);
                x1 = x2;
                y1 = y2;
                x2 += dx;
                y2 = f(x2);
                if (fabs(dx)<desAcc) return x2;
                i++;
                
            }
            std::cerr<<"Secant: Maximum number of iterations reached!"<<std::endl;
            return NAN;
        } 
        };
        
        
    //From an initial guess and an idea of appropriate length scales 
    //(initialstep), find a root using the secant method  
    struct Secant {
        template <typename FunctPtr, typename T, int maxStep = 100 >
                T operator () (FunctPtr&& f, const T& guess, const T& initialStep, const T& desAcc){
            T x1 = guess;
            T x2 = guess + initialStep;
            T y1 = f(x1);
            T y2 = f(x2);
            //Orient so that y2 is closer to root.
            if (fabs(y2) > fabs(y1)) {
                std::swap(y1, y2);
                std::swap(x1, x2);
            }
            int i = 0;

            while (i < maxStep) {
                T dx = (x1 - x2) * y2 / (y2 - y1);
                x1 = x2;
                y1 = y2;
                x2 += dx;
                y2 = f(x2);
                if (fabs(dx)<desAcc) return x2;
                i++;
                
            }
            std::cerr<<"Secant: Maximum number of iterations reached!"<<std::endl;
            return NAN;
        } 
        //An overload for when f is expensive. Sets fmin to  y minimum found.
        template <typename FunctPtr, typename T, int maxStep = 100 >
                T operator () (FunctPtr& f, const T& guess, const T& initialStep, const T& desAcc, T& fmin){
            T x1 = guess;
            T x2 = guess + initialStep;
            T y1 = f(x1);
            T y2 = f(x2);
            //Orient so that y2 is closer to root.
            if (fabs(y2) > fabs(y1)) {
                std::swap(y1, y2);
                std::swap(x1, x2);
            }
            int i = 0;

            while (i < maxStep) {
                T dx = (x1 - x2) * y2 / (y2 - y1);
                x1 = x2;
                y1 = y2;
                x2 += dx;
                y2 = f(x2);
                if (fabs(dx)<desAcc){
                    fmin = y2;
                    return x2;
                }
                i++;
                
            }
            std::cerr<<"Secant: Maximum number of iterations reached!"<<std::endl;
            return T(0);
        } 
        };
        
        
         //From an initial guess and an idea of appropriate length scales 
    //(initialstep), find a root using the secant method  
    struct DirectedSecant {
        enum { TOO_HIGH=1,TOO_LOW=-1,OK=0,MAX_IT=2};
        template <typename FunctPtr, typename T, int maxStep = 100 >
                int operator () (FunctPtr& f, const T& guess, const T& initialStep, const T& desAcc, T& xout){
            T x1 = guess;
            T x2 = guess + initialStep;
            int check;
            T y1 = f(x1,check);
            T y2 = f(x2,check);
            //Orient so that y2 is closer to root.
            if (fabs(y2) > fabs(y1)) {
                std::swap(y1, y2);
                std::swap(x1, x2);
            }
            int i = 0;

            while (i < maxStep) {
                T dx = (x1 - x2) * y2 / (y2 - y1);
                x1 = x2;
                y1 = y2;
                x2 += dx;
                y2 = f(x2,check);
                if (check!=OK) {
                    xout = x2;
                    return check;
                }
                if (fabs(dx)<desAcc){
                    xout= x2;
                    return check;
                }
                i++;    
            }
            std::cerr<<"Maximum number of iterations reached!"<<std::endl;
            return MAX_IT;
        }
        
    };  
        
        
    struct SecantBracket {
        template <typename FunctPtr, typename T, int maxStep = 4 >
                void operator () (FunctPtr& f, const T& guess, const T& initialStep, T& lower, T & upper){
            T x1 = guess;
            T x2 = guess + initialStep;
            T y1 = f(x1);
            T y2 = f(x2);
            //Orient so that y2 is closer to root.
            if (fabs(y2) > fabs(y1)) {
                std::swap(y1, y2);
                std::swap(x1, x2);
            }
            int i = 0;

            while (y1 * y2 >= 0.0 && i < maxStep) {
                T dx = (x1 - x2) * y2 / (y2 - y1);
                x1 = x2;
                y1 = y2;
                //2* as we are endeavouring to bracket, not solve
                x2 += 2*dx;
                y2 = f(x2);
                i++;
            }
            if (i == maxStep) {
                std::cerr << "SecantBracket: Failed to bracket root.\n";
                lower = T(0);
                upper = T(0);
            } else if (x1 > x2) {
                upper = x1;
                lower = x2;
            } else {
                upper = x2;
                lower = x1;
            }
        }    
        };
        
        
        
    struct SecantStep {
        template <typename T>
                void operator () (const T& f1,const T& f2, const T& guess,const T& initialStep){
            T x1 = guess;
            T x2 = guess + initialStep;
            T y1 = f1;
            T y2 = f2;
            //Orient so that y2 is closer to root.
            if (fabs(y2) > fabs(y1)) {
                std::swap(y1, y2);
                std::swap(x1, x2);
            }
                T dx = (x1 - x2) * y2 / (y2 - y1);
                x1 = x2;
                y1 = y2;
                x2 += dx;
                return x2;
        }
    };    
    
    
    // A basically mindless function that walks intialStep by intialStep until it finds a root    

    struct ConservativeBracket {

        template <typename FunctPtr, typename T, int maxStep = 1000000 >
                void operator () (FunctPtr& f, const T& guess, 
                const T& initialStep, T& lower, T & upper) {
            T x1 = guess;
            T x2;
            T y1 = f(x1);
            T y2;
            int i;
            for (i = 0; i < maxStep; i++) {
                x2 = x1 + initialStep;
                y2 = f(x2);
                if (y2 * y1 <= 0)
                    break;
                x1 = x2;
                y1 = y2;
            }
            if (i == maxStep) {
                std::cerr << "ConservativeBracket: Failed to bracket root.\n";
                lower = T(0);
                upper = T(0);
            } else if (x1 > x2) {
                upper = x1;
                lower = x2;
            } else {
                upper = x2;
                lower = x1;
            }
        }
    };
    
  
    
}



#endif	/* 1DROOTFINDING_H */

