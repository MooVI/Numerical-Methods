/* 
 * File:   PolynomialInterpolation.h
 * Author: qfeuille
 *
 * Created on 24 August 2012, 00:20
 */
#include<Eigen/Dense>
#ifndef POLYNOMIALINTERPOLATION_H
#define	POLYNOMIALINTERPOLATION_H
namespace NumMethod {

    //Optimised Version for when you know the order you'll be using
    template<typename T, int Order>
    struct FixedPolynomialInterpolation {
        typedef Eigen::Matrix<T, Order, 1 > Vector;

        inline T operator ()(const Vector& x,
                const Vector& y,
                const T& in, T & error) {
            Vector c(y), d(y);
            T diff = 0;
            int index = 0;
            T olddiff = fabs(in - x[0]);
            for (int i = 0; i < Order; i++) {
                if ((diff = fabs(in - x[i])) < olddiff) {
                    index = i;
                    olddiff = diff;
                }
            }
            T out = y [index--];
            for (int j = 1; j < Order; j++) {
                for (int i = 0; i < Order - j; i++) {
                    T step = x[i] - in;
                    T stepj = x[i + j] - in;
                    T temp = (c[i + 1] - d[i]) / (step - stepj);
                    d[i] = stepj*temp;
                    c[i] = step*temp;
                }
                out += (error = (2 * (index + 1)< (Order - j) ? c[index + 1] : d[index--]));
            }
            return out;
        };
    };
    
    
    //For when the order is unknown. Also has the advantage it can take most
    //containers as input.
     template<typename T>
    struct PolynomialInterpolation {
         template<typename X,typename Y>
        inline T operator ()(const X& x,
                const Y& y,
                const T& in, const int& Order, T & error) {
             std::vector<double> c (Order),d(Order);
             for (int i=0;i<Order;i++)
                 c[i]=d[i]=y[i];
            T diff = 0;
            int index = 0;
            T olddiff = fabs(in - x[0]);
            for (int i = 0; i < Order; i++) {
                if ((diff = fabs(in - x[i])) < olddiff) {
                    index = i;
                    olddiff = diff;
                }
            }
            T out = y [index--];
            for (int j = 1; j < Order; j++) {
                for (int i = 0; i < Order - j; i++) {
                    T step = x[i] - in;
                    T stepj = x[i + j] - in;
                    T temp = (c[i + 1] - d[i]) / (step - stepj);
                    d[i] = stepj*temp;
                    c[i] = step*temp;
                }
                out += (error = (2 * (index + 1)< (Order - j) ? c[index + 1] : d[index--]));
            }
            return out;
        };
    };
    
    
}
#endif	/* POLYNOMIALINTERPOLATION_H */

