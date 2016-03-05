/* 
 * File:   PRKMethods.h
 * Author: qfeuille
 *
 * Created on February 16, 2012, 4:58 PM
 */

#ifndef PRKMETHODS_H
#define	PRKMETHODS_H
#include "DESolversInc.h" 
namespace NumMethod {
 template <typename T, int N>
class SymplecticEuler {
        enum {Q = N / 2};
        Eigen::Matrix <T, Q, 1 > k;

        
    public:

        template<typename FunctPtr,typename A, typename B>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::MatrixBase<A> && p, Eigen::MatrixBase<B> && q, const Eigen::Matrix<T, Q, 1 > & kp, const T& stepSize) {
            p += stepSize * kp;
          
            q += stepSize*dqfunc.dq(p);
        };
    };

   template <typename T, int N>
    class OxfordImproved {
          enum {
            Q = N / 2
        };
        Eigen::Matrix <T, Q, 2 > k;
        
    public:

        template<typename FunctPtr,typename A, typename B>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::MatrixBase<A> && p, Eigen::MatrixBase<B> && q, const Eigen::Matrix<T, Q, 1 > & kp, const T& stepSize)  {
            k.col(0)=dqfunc.dq(p);
            dqfunc.dp(t + stepSize, q + stepSize*p, k.col(1));
            p += 0.5 * stepSize * (kp + k.col(1));
            q+= 0.5 * stepSize * (dqfunc.dq(k.col(0)) + dqfunc.dq(p));
        };

    };

    
     template <typename T, int N>
    class Ruth3 {
          enum {
            Q = N / 2
        };
        Eigen::Matrix<T,Q,1> k;
    static const T  a1,a2,a3,b1,b2,b3,c1,c2;    
        
    public:

        template<typename FunctPtr,typename A, typename B>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::MatrixBase<A> && p, Eigen::MatrixBase<B> && q, const Eigen::Matrix<T, Q, 1 > & kp, const T& stepSize)  {
            p += a1*stepSize * kp;
            q+= b1*stepSize * dqfunc.dq(p);
            dqfunc.dp(t +c1*stepSize, q, std::move(k));
             p += a2*stepSize * k;
            q+= b2*stepSize * dqfunc.dq(p);
            dqfunc.dp(t +c2*stepSize, q, std::move(k));
             p += a3*stepSize * k;
            q+= b3*stepSize * dqfunc.dq(p);
            
        };
        
    };
    
 template <typename T, int N> const T Ruth3<T, N>::a1 = 7.0 / 24.0;
 template <typename T, int N> const T Ruth3<T, N>::a2= 3.0/4.0;
 template <typename T, int N> const T Ruth3<T, N>::a3= -1.0/24.0;
 template <typename T, int N> const T Ruth3<T, N>::b1=  2.0/3.0;
 template <typename T, int N> const T Ruth3<T, N>::b2=  -2.0/3.0;
 template <typename T, int N> const T Ruth3<T, N>::b3= 1;
 template <typename T, int N> const T Ruth3<T, N>::c1= b1;
 template <typename T, int N> const T Ruth3<T, N>::c2= b1+b2;
   
 
   template <typename T, int N>
    class Ruth4 {
          enum {
            Q = N / 2
        };
         Eigen::Matrix<T,Q,1> k;
    static const T  a1,a2,a3,a4,a5,a6,b1,b2,b3,b4,b5,b6,c1,c2,c3,c4,c5;    
        
    public:

        template<typename FunctPtr,typename A, typename B>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::MatrixBase<A> && p, Eigen::MatrixBase<B> && q, const Eigen::Matrix<T, Q, 1 > & kp, const T& stepSize)  {
            
            q+= b1*stepSize * dqfunc.dq(p);
            dqfunc.dp(t +c1*stepSize, q, std::move(k));
            p += a1*stepSize * k;
            q+= b2*stepSize * dqfunc.dq(p);
            dqfunc.dp(t +c2*stepSize, q, std::move(k));
            p += a2*stepSize * k;
              q+= b3*stepSize * dqfunc.dq(p);
            dqfunc.dp(t +c3*stepSize, q, std::move(k));
            p += a3*stepSize * k;
              q+= b4*stepSize * dqfunc.dq(p);
            dqfunc.dp(t +c4*stepSize, q, std::move(k));
            p += a4*stepSize * k;
              q+= b5*stepSize * dqfunc.dq(p);
            dqfunc.dp(t +c5*stepSize, q, std::move(k));
            p += a5*stepSize * k;
              q+= b6*stepSize * dqfunc.dq(p);
            
        };
        
    };
    
 template <typename T, int N> const T Ruth4<T, N>::a1 = -1.0 / 24.0;
 template <typename T, int N> const T Ruth4<T, N>::a2= 3.0/4.0;
 template <typename T, int N> const T Ruth4<T, N>::a3= 7.0/12.0;
  template <typename T, int N> const T Ruth4<T, N>::a4 = 3.0 / 4.0;
 template <typename T, int N> const T Ruth4<T, N>::a5= -1.0/24.0;
 template <typename T, int N> const T Ruth4<T, N>::a6= 0;
  template <typename T, int N> const T Ruth4<T, N>::b1= 1;
 template <typename T, int N> const T Ruth4<T, N>::b2= - 2.0/3.0;
 template <typename T, int N> const T Ruth4<T, N>::b3=  2.0/3.0;
 template <typename T, int N> const T Ruth4<T, N>::b4=  2.0/3.0;
 template <typename T, int N> const T Ruth4<T, N>::b5=  -2.0/3.0;
  template <typename T, int N> const T Ruth4<T, N>::b6=  1;
 template <typename T, int N> const T Ruth4<T, N>::c1= b1;
 template <typename T, int N> const T Ruth4<T, N>::c2= b1+b2;
  template <typename T, int N> const T Ruth4<T, N>::c3= b1+b2+b3;
 template <typename T, int N> const T Ruth4<T, N>::c4= b1+b2+b3+b4;
  template <typename T, int N> const T Ruth4<T, N>::c5= b1+b2+b3+b4+b5;
}

#endif	/* PRKMETHODS_H */

