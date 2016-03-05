/* 
 * File:   RKMethods.h
 * Author: qfeuille
 *
 * Created on February 16, 2012, 4:41 PM
 */

#ifndef RKMETHODS_H
#define	RKMETHODS_H
#include "DESolversInc.h" 
namespace NumMethod {
 template <typename T, int N>
    class EulerForward {
        Eigen::Matrix <T, N, 1 > k;
    public:

        template<typename FunctPtr>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::Matrix<T, N, 1 > & q, const Eigen::Matrix<T, N, 1 > & initialdq, const T& stepSize) {

            q += stepSize * initialdq;
        };
        EulerForward();

    };

    template<typename T, int N>
    EulerForward<T, N>::EulerForward() {
    };

    template <typename T, int N>
    class UnpSymplecticEuler {
        Eigen::Matrix <T, N, 1 > k;

        enum {
            Q = N / 2
        };
    public:

        template<typename FunctPtr>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::Matrix<T, N, 1 > & q, const Eigen::Matrix<T, N, 1 > & k0, const T& stepSize) {
            q.template head<Q > () += stepSize * k0.template head<Q > ();
            q.template tail<Q > () += stepSize * q.template head<Q > ();
        };
    };
    
    
    template <typename T, int N>
    class UnpOxfordImproved {
        Eigen::Matrix <T, N, 2 > k;

        enum {
            Q = N / 2
        };
    public:

        template<typename FunctPtr>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::Matrix<T, N, 1 > & q, const Eigen::Matrix<T, N, 1 > & k0, const T& stepSize) {
            dqfunc(t + stepSize, q + stepSize*k0, k.col(1));
            q.template head<Q > () += 0.5 * stepSize * (k0.template head<Q > () + k.col(1).template head<Q > ());
            q.template tail<Q > () += 0.5 * stepSize * (k0.template tail<Q > () + q.template head<Q > ());
        };

    };
    
    
    template <typename T, int N>
    class Improved {
        Eigen::Matrix <T, N, 2 > k;
    public:

        template<typename FunctPtr>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::Matrix<T, N, 1 > & q, const Eigen::Matrix<T, N, 1 > & k0, const T& stepSize) {

            dqfunc(t + stepSize, q + stepSize*k0, k.col(1));

            q += 0.5 * stepSize * (k0 + k.col(1));
        };

    };

    template <typename T, int N>
    class RK4 {
        Eigen::Matrix <T, N, 4 > k;
    public:

        template<typename FunctPtr>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::Matrix<T, N, 1 > & q, const Eigen::Matrix<T, N, 1 > & k0, const T& stepSize) {
            dqfunc(t + 0.5 * stepSize, q + 0.5 * stepSize*k0, k.col(1));
            dqfunc(t + 0.5 * stepSize, q + 0.5 * stepSize * k.col(1), k.col(2));
            dqfunc(t + stepSize, q + stepSize * k.col(2), k.col(3));
            q += (1 / (double) 6) * stepSize * (k0 + 2 * k.col(1) + 2 * k.col(2) + k.col(3));
        };
    };
    
    

    template<typename T, int N, int DEGREE, int J, int END>
    struct RKLoopUnroller {

        template<typename FunctPtr>
                inline static void apply(FunctPtr& dqfunc, const T& stepSize, const T& t,
                const Eigen::Matrix<T, DEGREE - 1, DEGREE - 1 > & coeff, const Eigen::Matrix<T, DEGREE - 1, 1 > & tcoeff,
                const Eigen::Matrix<T, N, 1 > & q,
                Eigen::Matrix<T, N, DEGREE>& k) {
            dqfunc(t + tcoeff(J) * stepSize,
                    q + stepSize * (k.template leftCols < J + 1 > () * coeff.col(J).template head < J + 1 > ()),
                    k.col(J + 1));
            RKLoopUnroller<double, N, DEGREE, J + 1, END>::apply(dqfunc, stepSize, t, coeff, tcoeff, q, k);
        }
    };

    template<typename T, int N, int DEGREE, int END>
    struct RKLoopUnroller<T, N, DEGREE, END, END> {

        template<typename FunctPtr>
                inline static void apply(FunctPtr& dqfunc, const T& stepSize, const T& t,
                const Eigen::Matrix<T, DEGREE - 1, DEGREE - 1 > & coeff, const Eigen::Matrix<T, DEGREE - 1, 1 > & tcoeff,
                const Eigen::Matrix<T, N, 1 > & q,
                Eigen::Matrix<T, N, DEGREE>& k) {
        }

    };

    template<typename T, int N, int DEGREE >
    class RungeKuttaGeneral {
        Eigen::Matrix<T, N, DEGREE> k;
        const Eigen::Matrix <T, DEGREE - 1, DEGREE - 1 > coeff;
        const Eigen::Matrix <T, DEGREE - 1, 1 > tcoeff;
        const Eigen::Matrix <T, DEGREE, 1 > fcoeff;

    public:

        template<typename FunctPtr>
        void operator()(FunctPtr& dqfunc, const T& t, Eigen::Matrix<T, N, 1 > & q, const Eigen::Matrix<T, N, 1 > & k0, const T& stepSize) {

            k.col(0) = std::move(k0);
            RKLoopUnroller<double, N, DEGREE, 0, DEGREE - 1 > ::apply(dqfunc, stepSize, t, coeff, tcoeff, q, k);
#if 0     
            for (int j = 0; j < (DEGREE - 1); j++)
                dqfunc(t + tcoeff(j) * stepSize,
                    q + stepSize * (k.leftCols(j + 1) * coeff.col(j).head(j + 1)),
                    k.col(j + 1));
#endif
            q.noalias() += stepSize * (k * fcoeff);
        }


        RungeKuttaGeneral(const Eigen::Matrix<T, DEGREE, DEGREE> icoeff);
    };

    template<typename T, int N, int DEGREE >
    RungeKuttaGeneral<T, N, DEGREE>::RungeKuttaGeneral(const Eigen::Matrix<T, DEGREE, DEGREE> icoeff) :
    coeff(icoeff.template topRightCorner<DEGREE - 1, DEGREE - 1 > ().transpose()),
    tcoeff(icoeff.col(0).template head<DEGREE - 1 > ()),
    fcoeff(icoeff.template bottomRows < 1 > ().transpose()) {
    }


 template <typename T, int N>
    class CashKarp54 {
        Eigen::Matrix <T, N, 6> k;
        Eigen::Matrix<T,N,1> error;
        static const T a21,
                       a31,a32,
                       a41,a42,a43,
                       a51,a52,a53,a54, 
                       a61,a62,a63,a64,a65;
        static const T c2,c3,c4,c5,c6;
        static const T b1,b3,b4,b6;
        static const T e1,e3,e4,e5,e6;
    public:

        template<typename FunctPtr>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::Matrix<T, N, 1 > & q, const Eigen::Matrix<T, N, 1 > & k0, const T& stepSize) {
            
            dqfunc(t + c2 * stepSize, 
                    q + a21 * stepSize*k0,
                    k.col(1));
            dqfunc(t + c3 * stepSize,
                    q + stepSize *(a31* k0+a32*k.col(1)),
                    k.col(2));
            dqfunc(t + c4 * stepSize,
                    q + stepSize *(a41*k0+a42*k.col(1)+a43*k.col(2)),
                    k.col(3) );
            dqfunc(t + c5 * stepSize,
                    q + stepSize *(a51*k0+a52*k.col(1)+a53*k.col(2)+a54*k.col(3)),
                    k.col(4));
             dqfunc(t + c6 * stepSize,
                    q + stepSize *(a61*k0+a62*k.col(1)+a63*k.col(2)+a64*k.col(3)+a65*k.col(4)),
                    k.col(5));
     
            q += stepSize * (b1*k0 + b3 * k.col(2) + b4* k.col(3) + b6*k.col(5));
            error = stepSize*(e1*k0 + e3 * k.col(2) + e4* k.col(3) + e5* k.col(4) + e6*k.col(5));
        };
        
        Eigen::Array <T,N,1> getError (){
            return error;
        }
    };
    
template <typename T, int N> const T CashKarp54<T, N>::a21 = 1.0/5.0;
template <typename T, int N> const T CashKarp54<T, N>::a31 = 3.0/ 40.0;
template <typename T, int N> const T CashKarp54<T, N>::a32 = 9.0 / 40.0;
template <typename T, int N> const T CashKarp54<T, N>::a41 = 3.0 / 10.0;
template <typename T, int N> const T CashKarp54<T, N>::a42 = -9.0 / 10.0;
template <typename T, int N> const T CashKarp54<T, N>::a43 = 6.0 / 5.0;
template <typename T, int N> const T CashKarp54<T, N>::a51 = -11.0 / 54.0;
template <typename T, int N> const T CashKarp54<T, N>::a52 = 5.0 / 2.0;
template <typename T, int N> const T CashKarp54<T, N>::a53 = -70.0 / 27.0;
template <typename T, int N> const T CashKarp54<T, N>::a54 = 35.0 / 27.0;
template <typename T, int N> const T CashKarp54<T, N>::a61 = 1631.0 / 55296.0;
template <typename T, int N> const T CashKarp54<T, N>::a62 = 175.0 / 512.0;
template <typename T, int N> const T CashKarp54<T, N>::a63 = 575.0 / 13824.0;
template <typename T, int N> const T CashKarp54<T, N>::a64 = 44275.0 / 110592.0;
template <typename T, int N> const T CashKarp54<T, N>::a65 = 253.0 / 4096.0;
template <typename T, int N> const T CashKarp54<T, N>::c2 = a21;
template <typename T, int N> const T CashKarp54<T, N>::c3 = a31+a32;
template <typename T, int N> const T CashKarp54<T, N>::c4 = a41+a42+a43;
template <typename T, int N> const T CashKarp54<T, N>::c5 = a51+a52+a53+a54;
template <typename T, int N> const T CashKarp54<T, N>::c6 = a61+a62+a63+a64+a65;
template <typename T, int N> const T CashKarp54<T, N>::b1 = 37.0/ 378.0;
template <typename T, int N> const T CashKarp54<T, N>::b3 = 250.0 / 621.0;
template <typename T, int N> const T CashKarp54<T, N>::b4 = 125.0 / 594.0;
template <typename T, int N> const T CashKarp54<T, N>::b6 = 512.0 / 1771.0;
template <typename T, int N> const T CashKarp54<T, N>::e1 = b1-2825.0 / 27648.0;
template <typename T, int N> const T CashKarp54<T, N>::e3 = b3-18575.0 / 48384.0;
template <typename T, int N> const T CashKarp54<T, N>::e4 = b4-13525.0 / 55296.0;
template <typename T, int N> const T CashKarp54<T, N>::e5 = -277.0 / 14336.0;
template <typename T, int N> const T CashKarp54<T, N>::e6 = b6-0.25;
    
}
#endif	/* RKMETHODS_H */

