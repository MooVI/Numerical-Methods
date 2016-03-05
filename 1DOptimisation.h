/* 
 * File:   1DOptimisation.h
 * Author: qfeuille
 *
 * Created on 02 May 2013, 15:54
 */

#include "Cyclic.h"
#include <math.h>
#ifndef ONEDOPTIMISATION_H
#define	ONEDOPTIMISATION_H


namespace NumMethod {
   
    class GoldenSection {

        void shiftleft(double & a, double& b, double&c, const double& in) {
            a = b;
            b = c;
            c = in;
        }
    
        const double R = 0.61803398875;
        const double C = 1.0 - R;
        
    public:
        template <typename FunctPtr, typename T, int maxStep = 100 >
                T operator () (FunctPtr& f, const T& lower, const T& upper, const T& mid, const T& desAcc, T &xmin){
            double x0, x1, x2, x3;

            x0 = lower;
            x3 = upper;
            if (fabs(upper - mid) > fabs(mid - lower)) {
                x1 = mid;
                x2 = mid + C * (upper - mid);
            } else {
                x2 = mid;
                x1 = mid - C * (mid - lower);
            }
            double f1 = f(x1), f2 = f(x2);

            while (fabs(x3 - x0) > desAcc * (fabs(x1) + fabs(x2))) {
                if (f2 < f1) {
                    shiftleft(x0, x1, x2, R * x2 + C * x3);
                    f1 = f2;
                    f2 = f(x2);
                } else {
                    shiftleft(x3, x2, x1, R * x1 + C * x0);
                    f2 = f1;
                    f1 = f(x1);
                }
            }
            if (f1 < f2) {
                xmin = x1;
                return f1;
            } else {
                xmin = x2;
                return f2;
            }
        }

GoldenSection(){};

    };
 
    
    
    
}


#endif	/* ONEDOPTIMISATION_H */

