/* 
 * File:   Numerovlike.h
 * Author: qfeuille
 *
 * Created on 17 March 2013, 18:46
 */

#ifndef NUMEROVLIKE_H
#define	NUMEROVLIKE_H
#include "DESolversInc.h"
namespace NumMethod {

    template <typename T, int N>

    
    class Numerov {
        Eigen::Matrix <T, N, 1 > qtemp;
        const T fivesixth, twelth;

    public:
        //Take a single step using the method using the supplied values of f
        // and y_n (q) and y_n-1 (past)
        template<typename FunctPtr>
        void operator() (FunctPtr& ddqfunc, const T& t, 
                Eigen::Matrix<T, N, 1 > & q, const T& stepSize,
                Eigen::Matrix<T, N, 1 > & past, const Eigen::Matrix <T, N, 3 >& f,
                Cyclic<3>& cyclic) {

            qtemp = q;
            T dsq = stepSize*stepSize;
            q = (((2 + dsq * fivesixth * f.col(cyclic).array()) * q.array())
                - (1 - dsq * twelth * f.col(cyclic - 1).array()) * past.array())
                    / (1 - dsq * twelth * f.col(cyclic + 1).array());
            past = qtemp;
            cyclic++;
        };

        //Return the relative error per step for the variable stepper
        T getError(const T& f, const T& t, const T& stepSize) {
            T error = f * twelth * stepSize*stepSize;
            return fabs(error);
        }

        //Estimate the relative accuracy needed for a desired global acc
        //from max value of f in region and the integration range.
        T getRelAcc(const T& maxf, const T& range, const T& desAcc, 
                int * numSteps = NULL) {
            T minStep = pow(fabs(desAcc * 240.0 / (maxf * maxf * maxf * range))
                            , 0.25);
            int num = (int) (range / minStep);
            if (numSteps != NULL)
                *numSteps = num;
            return desAcc / (T) num;
        }

       //Estimate necessary initialstepsize
        template<typename FunctPtr>
        T getInitialStep(FunctPtr& ddqfunc, const T& begin, const T& relAcc) {
            T temp = fabs(12.0 * pow(relAcc / 7.2, 1 / 3.0) / ddqfunc(begin));
            temp = temp < 1.0 ? temp : 1.0;
            return sqrt(temp);
        }

        Numerov() : fivesixth(5.0 / 6.0), twelth(1.0 / 12.0) {
        };
    };
}

#endif	/* NUMEROVLIKE_H */

