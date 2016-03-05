/* 
 * File:   RKNMethods.h
 * Author: qfeuille
 *
 * Created on February 16, 2012, 5:00 PM
 */

#ifndef RKNMETHODS_H
#define	RKNMETHODS_H
#include "DESolversInc.h" 
namespace NumMethod {
 
   template <typename T, int N>
    class Calvo4 {
          enum {
            Q = N / 2
        };
       
    static const T  a1,a2,a3,a4,a5,b1,b2,b3,b4,b5,g1,g2,g3,g4,g5,c1,c2,c3,c4,c5;    
        
    public:

        template<typename FunctPtr,typename A, typename B>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::MatrixBase<A> && p, Eigen::MatrixBase<B> && q, Eigen::Matrix<T, Q, 1 > & kp, const T& stepSize)  {
            
            p += b1*stepSize * kp;
            q+= a2*stepSize * dqfunc.dq(p);
            dqfunc.dp(t +c2*stepSize, q, std::move(kp));
            p += b2*stepSize * kp;
              q+= a3*stepSize * dqfunc.dq(p);
            dqfunc.dp(t +c3*stepSize, q, std::move(kp));
            p += b3*stepSize * kp;
              q+= a4*stepSize * p;
            dqfunc.dp(t +c4*stepSize, q, std::move(kp));
            p += b4*stepSize * kp;
              q+= a5*stepSize * p;
            dqfunc.dp(t +c5*stepSize, q, std::move(kp));
            p += b5*stepSize * kp;
            dqfunc.dp(t +c5*stepSize, q, std::move(kp));
            
        };
        
    };
   
    

 template <typename T, int N> const T Calvo4<T, N>::g2= 0.2051776615422863869;
 template <typename T, int N> const T Calvo4<T, N>::g3= 0.6081989431465009739;
  template <typename T, int N> const T Calvo4<T, N>::g4 = 0.4872780668075869657;
 template <typename T, int N> const T Calvo4<T, N>::g5= 1.0;  
 template <typename T, int N> const T Calvo4<T, N>::a2= g2;
 template <typename T, int N> const T Calvo4<T, N>::a3= g3-g2;
  template <typename T, int N> const T Calvo4<T, N>::a4 = g4-g3;
 template <typename T, int N> const T Calvo4<T, N>::a5= g5-g4;
  template <typename T, int N> const T Calvo4<T, N>::b1= 0.061758858135626325;
 template <typename T, int N> const T Calvo4<T, N>::b2= 0.3389780265536433551;
 template <typename T, int N> const T Calvo4<T, N>::b3=  0.614791307175577662;
 template <typename T, int N> const T Calvo4<T, N>::b4=  -0.1405480146593733802;
 template <typename T, int N> const T Calvo4<T, N>::b5=  0.1250198227945261338;
 template <typename T, int N> const T Calvo4<T, N>::c1= b1;
 template <typename T, int N> const T Calvo4<T, N>::c2= b1+b2;
  template <typename T, int N> const T Calvo4<T, N>::c3= b1+b2+b3;
 template <typename T, int N> const T Calvo4<T, N>::c4= b1+b2+b3+b4;
  template <typename T, int N> const T Calvo4<T, N>::c5= b1+b2+b3+b4+b5;
  
    template <typename T, int N>
    class Calvo8 {
          enum {
            Q = N / 2
        };
       
    static const T  a2,a3,a4,a5,a6,a7,a8,a9,a10,a11,a12,a13;
     static const T  b1,b2,b3,b4,b5,b6,b7,b8,b9,b10,b11,b12,b13;
      static const T  c1,c2,c3,c4,c5,c6,c7,c8,c9,c10,c11,c12;
      static const T g2,g3,g4,g5,g6,g7,g8,g9,g10,g11,g12,g13;
        
    public:

        template<typename FunctPtr,typename A, typename B>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::MatrixBase<A> && p, Eigen::MatrixBase<B> && q, Eigen::Matrix<T, Q, 1 > & kp, const T& stepSize)  {
            
            p += b1*stepSize * kp;
            q+= a2*stepSize * dqfunc.dq(p);
            dqfunc.dp(t +c2*stepSize, q, std::move(kp));
            p += b2*stepSize * kp;
              q+= a3*stepSize * dqfunc.dq(p);
            dqfunc.dp(t +c3*stepSize, q, std::move(kp));
            p += b3*stepSize * kp;
              q+= a4*stepSize * p;
            dqfunc.dp(t +c4*stepSize, q, std::move(kp));
            p += b4*stepSize * kp;
              q+= a5*stepSize * p;
            dqfunc.dp(t +c5*stepSize, q, std::move(kp));
            p += b5*stepSize * kp;
             q+= a6*stepSize * p;
            dqfunc.dp(t +c6*stepSize, q, std::move(kp));
            p += b6*stepSize * kp;
             q+= a7*stepSize * p;
            dqfunc.dp(t +c7*stepSize, q, std::move(kp));
            p += b7*stepSize * kp;
             q+= a8*stepSize * p;
            dqfunc.dp(t +c8*stepSize, q, std::move(kp));
            p += b8*stepSize * kp;
             q+= a9*stepSize * p;
            dqfunc.dp(t +c9*stepSize, q, std::move(kp));
            p += b9*stepSize * kp;
             q+= a10*stepSize * p;
            dqfunc.dp(t +c10*stepSize, q, std::move(kp));
            p += b10*stepSize * kp;
             q+= a11*stepSize * p;
            dqfunc.dp(t +c11*stepSize, q, std::move(kp));
            p += b11*stepSize * kp;
             q+= a12*stepSize * p;
            dqfunc.dp(t +c12*stepSize, q, std::move(kp));
            p += b12*stepSize * kp;
               q+= a13*stepSize * p;
            dqfunc.dp(t +c12*stepSize, q, std::move(kp));
            p += b13*stepSize * kp;
             dqfunc.dp(t +c12*stepSize, q, std::move(kp));
        };
        
    };
   
    

 template <typename T, int N> const T Calvo8<T, N>::g2= 0.60715821186110352503;
 template <typename T, int N> const T Calvo8<T, N>::g3= 0.96907291059136392378;
  template <typename T, int N> const T Calvo8<T, N>::g4 =- 0.10958316365513620399;
 template <typename T, int N> const T Calvo8<T, N>::g5= 0.05604981994113413605;  
 template <typename T, int N> const T Calvo8<T, N>::g6= 1.3088652991863123401;
 template <typename T, int N> const T Calvo8<T, N>::g7= -0.11642101198009154794;
  template <typename T, int N> const T Calvo8<T, N>::g8 = -0.29931245499473964831;
 template <typename T, int N> const T Calvo8<T, N>::g9= -0.16586962790248628655;  
 template <typename T, int N> const T Calvo8<T, N>::g10= 1.22007054181677755238;
 template <typename T, int N> const T Calvo8<T, N>::g11= 0.20549254689579093228;
template <typename T, int N> const T Calvo8<T, N>::g12= 0.86890893813102759275;
template <typename T, int N> const T Calvo8<T, N>::g13= 1;
 template <typename T, int N> const T Calvo8<T, N>::a2= g2;
 template <typename T, int N> const T Calvo8<T, N>::a3= g3-g2;
  template <typename T, int N> const T Calvo8<T, N>::a4 = g4-g3;
 template <typename T, int N> const T Calvo8<T, N>::a5= g5-g4;
 template <typename T, int N> const T Calvo8<T, N>::a6= g6-g5;
  template <typename T, int N> const T Calvo8<T, N>::a7 = g7-g6;
 template <typename T, int N> const T Calvo8<T, N>::a8= g8-g7;
 template <typename T, int N> const T Calvo8<T, N>::a9= g9-g8;
  template <typename T, int N> const T Calvo8<T, N>::a10 = g10-g9;
 template <typename T, int N> const T Calvo8<T, N>::a11= g11-g10;
 template <typename T, int N> const T Calvo8<T, N>::a12= g12-g11;
  template <typename T, int N> const T Calvo8<T, N>::a13= g13-g12;
  template <typename T, int N> const T Calvo8<T, N>::b1= g2/2.0;
 template <typename T, int N> const T Calvo8<T, N>::b2= (g3)/2.0;
 template <typename T, int N> const T Calvo8<T, N>::b3= (g4-g2)/2.0;
 template <typename T, int N> const T Calvo8<T, N>::b4= (g5-g3)/2.0;
 template <typename T, int N> const T Calvo8<T, N>::b5= (g6-g4)/2.0;
 template <typename T, int N> const T Calvo8<T, N>::b6= (g7-g5)/2.0;
 template <typename T, int N> const T Calvo8<T, N>::b7= (g8-g6)/2.0;
 template <typename T, int N> const T Calvo8<T, N>::b8= (g9-g7)/2.0;
 template <typename T, int N> const T Calvo8<T, N>::b9= (g10-g8)/2.0;
 template <typename T, int N> const T Calvo8<T, N>::b10= (g11-g9)/2.0;
 template <typename T, int N> const T Calvo8<T, N>::b11= (g12-g10)/2.0;
 template <typename T, int N> const T Calvo8<T, N>::b12= (g13-g11)/2.0;
 template <typename T, int N> const T Calvo8<T, N>::b13= (1-g12)/2.0;
 template <typename T, int N> const T Calvo8<T, N>::c1= b1;
 template <typename T, int N> const T Calvo8<T, N>::c2= c1+b2;
  template <typename T, int N> const T Calvo8<T, N>::c3= c2+b3;
 template <typename T, int N> const T Calvo8<T, N>::c4= c3+b4;
  template <typename T, int N> const T Calvo8<T, N>::c5= c4+b5;
  template <typename T, int N> const T Calvo8<T, N>::c6= c5+b6;
  template <typename T, int N> const T Calvo8<T, N>::c7= c6+b7;
  template <typename T, int N> const T Calvo8<T, N>::c8= c7+b8;
  template <typename T, int N> const T Calvo8<T, N>::c9= c8+b9;
  template <typename T, int N> const T Calvo8<T, N>::c10= c9+b10;
  template <typename T, int N> const T Calvo8<T, N>::c11= c10+b11;
   template <typename T, int N> const T Calvo8<T, N>::c12= c11+b12;
  
    template <typename T, int N>
    class VelocityVerlet {
         enum {
            Q = N / 2
        };
        Eigen::Matrix <T, Q, 2 > k;
        
    public:

        template<typename FunctPtr,typename A,typename B>
        void operator() (FunctPtr& dqfunc, const T& t, Eigen::MatrixBase<A> && p, Eigen::MatrixBase<B> && q, 
                Eigen::Matrix<T, Q, 1 > & dp, const T& stepSize) {
            q += (stepSize * p + 0.5 * stepSize * stepSize * dp);
            dqfunc.dp(t + stepSize, q, k.col(1));
            p += 0.5 * stepSize * (dp + k.col(1));
            dp=k.col(1);
        };

    };
}
#endif	/* RKNMETHODS_H */

