/* 
 * File:   FixedSteppers.h
 * Author: qfeuille
 *
 * Created on February 16, 2012, 4:37 PM
 */

#ifndef FIXEDSTEPPERS_H
#define	FIXEDSTEPPERS_H

#include "DESolversInc.h" 
namespace NumMethod {

    template <typename T, int N>
    class FixedODESolver {
        Eigen::Matrix<T, N, 1 > initialdq;
    public:
        template<typename FunctPtr, typename Interrupt, typename Method>
        void solve(Method& step, FunctPtr& dqfunc, const int& numSteps, Eigen::Matrix<T, N, 1 > & q,
                T &t, const T& end, Interrupt& interrupt);
    };

    template <typename T, int N>
    template<typename FunctPtr, typename Interrupt, typename Method>
    inline void FixedODESolver<T, N>::solve(Method& step, FunctPtr& dqfunc,
            const int& numSteps, Eigen::Matrix<T, N, 1 > & q, T &t, const T& end, Interrupt& interrupt) {
        T stepSize = (end - t) / (T) numSteps;
        if (!interrupt(t, q, stepSize)) {
            for (int i = 0; i < numSteps; ++i) {
                dqfunc(t, q, std::move(initialdq));
                step(dqfunc, t, q, initialdq, stepSize);
                t += stepSize;
                if (interrupt(t, q, stepSize)) break;
            }
        }
    }

    template <typename T, int N>
    class PartFixedODESolver {

        enum {
            Q = N / 2
        };
        Eigen::Matrix<T, Q, 1 > idp;
    public:
        template<typename FunctPtr, typename Interrupt, typename Method>
        void solve(Method& step, FunctPtr& dqfunc, const int& numSteps, Eigen::Matrix<T, N, 1 > & q,
                T &t, const T& end, Interrupt& interrupt);
    };

    template <typename T, int N>
    template<typename FunctPtr, typename Interrupt, typename Method>
    inline void PartFixedODESolver<T, N>::solve(Method& step, FunctPtr& dqfunc,
            const int& numSteps, Eigen::Matrix<T, N, 1 > & q, T &t, const T& end, Interrupt& interrupt) {
        T stepSize = (end - t) / (T) numSteps;
        if (!interrupt(t, q, stepSize)) {
            for (int i = 0; i < numSteps; ++i) {

                dqfunc.dp(t, q.template tail<Q > (), std::move(idp));
                step(dqfunc, t, q.template head<Q > (), q.template tail<Q > (), idp, stepSize);
                t += stepSize;
                if (interrupt(t, q, stepSize)) break;
            }
        }
    }

    template <typename T, int N>
    class FSALPartFixedODESolver {

        enum {
            Q = N / 2
        };
        Eigen::Matrix<T, Q, 1 > idp;
    public:
        template<typename FunctPtr, typename Interrupt, typename Method>
        void solve(Method& step, FunctPtr& dqfunc, const int& numSteps, Eigen::Matrix<T, N, 1 > & q,
                T &t, const T& end, Interrupt& interrupt);
    };

    template <typename T, int N>
    template<typename FunctPtr, typename Interrupt, typename Method>
    inline void FSALPartFixedODESolver<T, N>::solve(Method& step, FunctPtr& dqfunc,
            const int& numSteps, Eigen::Matrix<T, N, 1 > & q, T &t, const T& end, Interrupt& interrupt) {
        T stepSize = (end - t) / (T) numSteps;
        dqfunc.dp(t, q.template tail<Q > (), std::move(idp));
        if (!interrupt(t, q, stepSize)) {
            for (int i = 0; i < numSteps; ++i) {
                step(dqfunc, t, q.template head<Q > (), q.template tail<Q > (), idp, stepSize);
                t += stepSize;
                if (interrupt(t, q, stepSize)) break;
            }
        }
    }

    template <typename T, int N>
    class MultiStepFixedODESolver {
        Eigen::Matrix<T, N, 1 > initialdq;
    public:
        template<typename FunctPtr, typename Interrupt, typename Method>
        void solve(Method& step, FunctPtr& dqfunc, const int& numSteps, Eigen::Matrix<T, N, 1 > & q,
                T &t, const T& end, Interrupt& interrupt);
    };

    template <typename T, int N>
    template<typename FunctPtr, typename Interrupt, typename Method>
    inline void MultiStepFixedODESolver<T, N>::solve(Method& step, FunctPtr& dqfunc,
            const int& numSteps, Eigen::Matrix<T, N, 1 > & q, T &t, const T& end, Interrupt& interrupt) {
        T stepSize = (end - t) / (T) numSteps;
        step.start(dqfunc, t, q, stepSize);
        if (!interrupt(t, q, stepSize)) {
            for (int i = 0; i < numSteps; ++i) {
                step(dqfunc, t, q, stepSize);
                t += stepSize;
                if (interrupt(t, q, stepSize)) break;
            }
        }
    }

    template <typename T, int N>
    class NumerovFixedODESolver {
      Eigen::Matrix <T, 1, 3 > f;
      NumMethod::Cyclic < 3 > cyclic;
    public:
        template<typename FunctPtr, typename Interrupt, typename Method>
        void solve(Method& step, FunctPtr& dqfunc, const int& numSteps, Eigen::Matrix<T, N, 1 > & q,
                T &t, const T& end, Eigen::Matrix<T, N, 1 > & past, Interrupt& interrupt, bool cont = false);

    };

    template <typename T, int N>
    template<typename FunctPtr, typename Interrupt, typename Method>
    inline void NumerovFixedODESolver<T, N>::solve(Method& step, FunctPtr& ddqfunc,
            const int& numSteps, Eigen::Matrix<T, N,1> & q, T &t, const T& end,Eigen::Matrix<T, N, 1 > & past, Interrupt& interrupt,bool cont) {
        T stepSize = (end - t) / (T) numSteps;
        if (!cont){ 
        ddqfunc (t, f.col(0));
        ddqfunc ((t+=stepSize), f.col(1));
        cyclic=1;
        }
        if (!interrupt(t, q, stepSize)) {
            for (int i = 0; i < (numSteps - (!cont)); ++i) {
                ddqfunc(t += stepSize, f.col(cyclic + 1));
                step(ddqfunc, t, q, stepSize,past,f,cyclic);
                if (interrupt(t, q, stepSize)) break;
            }
        }
    }


}

#endif	/* FIXEDSTEPPERS_H */

