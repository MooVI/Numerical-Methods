/* 
 * File:   VarSteppers.h
 * Author: qfeuille
 *
 * Created on February 16, 2012, 4:39 PM
 */

#ifndef VARSTEPPERS_H
#define	VARSTEPPERS_H
#include "DESolversInc.h" 
#include "PolynomialInterpolation.h"
namespace NumMethod {
   template <typename T, int N>
    class StepDoublerODESolver {
        Eigen::Matrix<T, N, 1 > initialdq, initialdq2;
        Eigen::Matrix<T, N, 1 > q1, q2;
       // Eigen::Array<T, N, 1 > qscal;
        double error, trialStep, power, errcon;
    public:
        template<typename FunctPtr, typename Interrupt, typename Method,typename Errorf>
        void solve(Method& step, FunctPtr& dqfunc, Eigen::Matrix<T, N, 1 > & q, T &t, const T& end, Interrupt& interrupt,
                const int& order, double& stepSize1, const double& desAcc,Errorf calcError,
                int& finalNumSteps, const double& pgrow = 5, const double& pshrink = 0.1, const double& safety = 0.9);
    };

    template <typename T, int N>
    template<typename FunctPtr, typename Interrupt, typename Method,typename Errorf>
    inline void StepDoublerODESolver<T, N>::solve(Method& step, FunctPtr& dqfunc, Eigen::Matrix<T, N, 1 > & q, T &t, const T& end, Interrupt& interrupt,
            const int& order, double& stepSize, const double& desAcc,Errorf calcError, int& finalNumSteps,
            const double& pgrow, const double& pshrink, const double& safety) {
        power = -1 / (double) (order + 1);
        errcon = pow((pgrow / safety), 1 / power);
        const int orient = t > end ? -1 : 1;
        if (!interrupt(t, q, stepSize)) {
        while (t*orient < end*orient) {
            dqfunc(t, q, std::move(initialdq));
            while (true) {
                q1 = q;
                q2 = q;
                step(dqfunc, t, q1, initialdq, stepSize);
                step(dqfunc, t, q2, initialdq, stepSize / 2.0);
                dqfunc(t + stepSize / 2.0, q2, std::move(initialdq2));
                step(dqfunc, t + stepSize / 2.0, q2, initialdq2, stepSize / 2.0);
              //  qscal = q.array().abs() + stepSize * initialdq.array().abs() + 1.0e-30;
                error = calcError((q2 - q1).array().abs(),q,initialdq, stepSize,desAcc,t);
             //   error /= desAcc;
                if (error <= 1.0) break;
                trialStep = safety * stepSize * pow(error, power);
                stepSize = NumMethod::signmaxAbs(trialStep, pshrink * stepSize);
                assert(t + stepSize != t);
            }
            finalNumSteps++;
            t += stepSize;
            q = q2;
            if (interrupt(t, q, stepSize)) break;
            stepSize = error > errcon ? safety * stepSize * pow(error, power) : pgrow*stepSize;
            stepSize = NumMethod::signminAbs(stepSize, end - t);
        }
        }
    }

      template <typename T, int N>
    class EmbeddedPairODESolver {
        Eigen::Matrix<T, N, 1 > initialdq;
        Eigen::Matrix<T, N, 1 > q1;
      //  Eigen::Array<T, N, 1 > qscal;
        double error, trialStep, spower,gpower, errcon;
    public:
        std::vector<double> errors;
        template<typename FunctPtr, typename Interrupt, typename Method,typename Errorf>
        void solve(Method& step, FunctPtr& dqfunc, Eigen::Matrix<T, N, 1 > & q, T &t, const T& end, Interrupt& interrupt,
                const int& order, double& stepSize1, const double& desAcc,Errorf calcError,
                int& finalNumSteps, const double& pgrow = 5, const double& pshrink = 0.1, const double& safety = 0.9);
    };

    template <typename T, int N>
    template<typename FunctPtr, typename Interrupt, typename Method,typename Errorf>
    inline void EmbeddedPairODESolver<T, N>::solve(Method& step, FunctPtr& dqfunc, Eigen::Matrix<T, N, 1 > & q, T &t, const T& end, Interrupt& interrupt,
            const int& order, double& stepSize, const double& desAcc,Errorf calcError, int& finalNumSteps,
            const double& pgrow, const double& pshrink, const double& safety) {
        spower = -1 / (double) (order );
        gpower = -1 / (double) (order);
        errcon = pow((pgrow / safety), 1 / gpower);
        const int orient = t > end ? -1 : 1;
        if (!interrupt(t, q, stepSize)) {
        while (orient*t < end*orient) {
            dqfunc(t, q, std::move(initialdq));
            while (true) {
                q1 = q;
                step(dqfunc, t, q1, initialdq, stepSize);
               // qscal = q.array().abs() + stepSize * initialdq.array().abs() + 1.0e-30;
                error = calcError(step.getError().abs(), q,initialdq, stepSize, desAcc,t);
               // error /= desAcc;
                if (error <= 1.0) break;
                trialStep = safety * stepSize * pow(error, spower);
                stepSize = NumMethod::maxAbs(trialStep, pshrink * stepSize);
                assert(t + stepSize != t);
            }
            finalNumSteps++;
            t += stepSize;
            q = q1;
            if (interrupt(t, q, stepSize)) break;
            errors.push_back(step.getError().abs().maxCoeff());
           // stepSize = error > errcon ? safety * stepSize * pow(error, gpower) : pgrow*stepSize;
            stepSize=NumMethod::minAbs ( safety * stepSize * pow(error, gpower) , pgrow*stepSize);
            stepSize = NumMethod::signminAbs(stepSize, end - t);
        }
       }
    }

  template <typename T>
    class NumerovVarODESolver {
      FixedPolynomialInterpolation<T,4> interpolate;
      Eigen::Matrix<T,4,1> storeq;
      Eigen::Matrix<T,4,1> storet;
      Eigen::Matrix<T,1,1> past2;
      Eigen::Matrix <T, 1, 3 > f;
      NumMethod::Cyclic < 3 > cyclic;
      const T twelth = 1.0/12.0;
      template<typename FunctPtr>
      void midpoint (FunctPtr& ddqfunc, Eigen::Matrix<T,1,1>& past, Eigen::Matrix<T,1,1>&q, const T& t, const T& stepSize){
          T hsqtwelve = twelth*stepSize*stepSize;
          T fhalf  = ddqfunc (t);
          past[0] = ((1-hsqtwelve*f[cyclic])*q[0] + (1 - hsqtwelve*f[cyclic-1])*past[0])/ (2 + 10* hsqtwelve*fhalf);
          f[cyclic-1] = fhalf;
      };
      
    public:
        template<typename FunctPtr, typename Interrupt, typename Method>
        void solve(Method& step, FunctPtr& dqfunc,T& stepSize, Eigen::Matrix<T, 1, 1 > & q,T &t, const T& end,const T & relAcc, Eigen::Matrix<T, 1, 1 > & past, Interrupt& interrupt);
        
    };

    template <typename T>
    template<typename FunctPtr, typename Interrupt, typename Method>
    inline void NumerovVarODESolver<T>::solve(Method& step, FunctPtr& ddqfunc,
   T& stepSize, Eigen::Matrix<T, 1, 1 > & q, T &t, const T& end, const T & relAcc,Eigen::Matrix<T, 1,1> & past, Interrupt& interrupt) {
        ddqfunc (t, f.col(0));
        ddqfunc ((t+=stepSize), f.col(1));
        cyclic=1;
        T del = pow (relAcc, 1/3.0)/7.2;
        T error;
        bool stepdouble = false;
        const int orient = t > end ? -1 : 1;
        if (!interrupt(t, q, stepSize)) {
            while (t * orient <= (end - 2 * stepSize) * orient) {
                error = step.getError(f[cyclic], t, stepSize);
                if (5 * error < del and stepdouble) {
                    stepSize *= 2;
                    past = past2;
                    f.col(cyclic - 1) = f.col(cyclic + 1);
                    stepdouble = false;
                } else {
                    stepdouble = true;
                    if (error > del) {
                        do {
                            stepSize /= 2;
                            midpoint(ddqfunc, past, q, t-stepSize, stepSize);                            
                        } while ((error /= 4) > del);
                    }
                }
                past2 = past;
                ddqfunc(t += stepSize, f.col(cyclic + 1));
                step(ddqfunc, t, q, stepSize, past, f,cyclic);
                if (interrupt(t, q, stepSize)) break;
            }
        storet [0] = t;
        storeq [0] = q[0];
        for (int i=1;i<4;i++) {
            ddqfunc(t + stepSize, f.col(cyclic + 1));
            step(ddqfunc, t, q, stepSize,past,f,cyclic);
            t += stepSize;
            storet [i] = t;
            storeq [i] = q[0];
        }
        t = end-stepSize;
        q[0] =interpolate (storet,storeq,t,error);
        if (interrupt (t,q,stepSize)) return;
        t+=stepSize;
        past =q;
        q[0] =interpolate (storet,storeq,t,error);
        interrupt (t,q,stepSize);
        }
    }


struct NRVarStepError {
    
    template<typename T, typename A, int N> 
    inline T operator ()(const Eigen::ArrayBase<A>&& error,
    const Eigen::Matrix <T,N,1>& q,
    const Eigen::Matrix <T,N,1>& dq,
    const T& stepSize,
    const T& desAcc,
    const T& t
    )
    {
    return (error/(q.array().abs() + stepSize * dq.array().abs() + 1.0e-30)).maxCoeff()/desAcc;
    }
    
    
};
}
#endif	/* VARSTEPPERS_H */

