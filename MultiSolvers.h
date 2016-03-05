/* 
 * File:   MultiSolvers.h
 * Author: qfeuille
 *
 * Created on February 16, 2012, 5:05 PM
 */

#ifndef MULTISOLVERS_H
#define	MULTISOLVERS_H
#include "DESolversInc.h" 
namespace NumMethod {
 template <typename T, int N>
    class MultiVelocityVerlet {
        Eigen::Matrix <T, N, 2 > k;
        NumMethod::Cyclic < 2 > cyclic;

        enum {
            Q = N / 2
        };
    public:

        template<typename FunctPtr>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::Matrix<T, N, 1 > & q, const T& stepSize) {
            q.template tail<Q > () += (stepSize * q.template head<Q > () + 0.5 * stepSize * stepSize * k.col(cyclic()).template head<Q > ());
            dqfunc(t + stepSize, q, k.col(cyclic + 1));
            q.template head<Q > () += 0.5 * stepSize * (k.col(cyclic()).template head<Q > () + k.col(cyclic + 1).template head<Q > ());
            cyclic++;


        };

        template<typename FunctPtr>
        void start(FunctPtr& dqfunc, const T& t, Eigen::Matrix<T, N, 1 > & q, const T& stepSize) {
            dqfunc(t, q, k.col(0));
            cyclic = 0;
        }
    };
    
   


    

    template <typename T, int N>
    class Beeman {
        Eigen::Matrix <T, N, 3 > k;

        enum {
            Q = N / 2
        };
        const double sixth;
        NumMethod::Cyclic < 3 > cyclic;
    public:

        template<typename FunctPtr>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::Matrix<T, N, 1 > & q, const T& stepSize) {
            q.template tail<Q > () += stepSize * q.template head<Q > ()
                    + stepSize * stepSize * sixth *
                    (4 * k.col(cyclic()).template head<Q > ()
                    - k.col(cyclic - 1).template head<Q > ()
                    );

            dqfunc(t + stepSize, q, k.col(cyclic + 1));
            q.template head<Q > () += sixth * stepSize *
                    (2 * k.col(cyclic + 1).template head<Q > ()
                    + 5 * k.col(cyclic()).template head<Q > ()
                    - k.col(cyclic - 1).template head<Q > ()
                    );
            cyclic++;


        };

        template<typename FunctPtr>
        void start(FunctPtr& dqfunc, const T& t, Eigen::Matrix<T, N, 1 > & q, const T& stepSize) {
            dqfunc(t, q, k.col(1));
            dqfunc(t - stepSize, q - stepSize * k.col(1), k.col(0));
            cyclic = 1;
        }

        Beeman() : sixth(1 / (double) 6) {
        };
    };
    
    
}
#endif	/* MULTISOLVERS_H */

