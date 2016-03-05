/* 
 * File:   HermitePoly.h
 * Author: qfeuille
 *
 * Created on 23 March 2013, 00:44
 */
#include "Polynomial.h"
#include "Cyclic.h"
#ifndef HERMITEPOLY_H
#define	HERMITEPOLY_H
namespace NumMethod {
    
    //Calculates  hDegree th Hermite Polynomial coefficient for the pOrderth
    //power at compile time.
    template <int hDegree, int pOrder> struct hermiteCoeff{
            enum {value = !(pOrder<=hDegree) ? 0 : (2* hermiteCoeff <hDegree-1,pOrder-1>::value - 2*(hDegree-1)*hermiteCoeff<hDegree-2,pOrder>::value)};
        };
        template <>
        struct hermiteCoeff<0,0> {
            enum {value = 1};
        };
       
         template <>
        struct hermiteCoeff<1,0> {
            enum {value = 0};
        };
        
         template <>
        struct hermiteCoeff<1,1> {
            enum {value = 2};
        };
        
          template <int hDegree>
        struct hermiteCoeff<hDegree,-1> {
            enum {value = 0};
        };
         template <int pOrder>
        struct hermiteCoeff<-1,pOrder> {
            enum {value = 0};
        };
        
         template <int pOrder>
        struct hermiteCoeff<-2,pOrder> {
            enum {value = 0};
        };
        
         template<int J, int END, typename T>
    struct Filler {

                inline static void apply(Eigen::Matrix<T,END,1>& coeff) {
                  coeff [J] = (T) hermiteCoeff<END-1,J>::value;
            Filler<J + 1, END,T>::apply(coeff);
        }
    };
    template <int END,typename T>
    struct Filler<END,END,T> {
        inline static void apply(Eigen::Matrix<T,END,1>& coeff) {};
    };
    
    //Fills a Polynomial class with Hermite coefficients calculated at compile time. 
    template <typename T, int Degree>
    struct ComplileTimeHermitePoly{
        
    
       static Polynomial<T> construct (){
            Eigen::Matrix<T,Degree,1> coeff;
            Filler <0,Degree+1,double>::apply(coeff);
            return Polynomial<T> (coeff);
        }
        
    
    };
    
    template<typename T>
    class HermitePoly {
        
        int n;
        
        public:
        
            //Use this for the coefficients, not for evaluation!
        static Polynomial<T> construct (const int& Degree){
            Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic> H 
                    = Eigen::Matrix<T,Eigen::Dynamic,Eigen::Dynamic>::Zero(Degree+1,3);
            H (0,0) =1;
            if (Degree==0) return Polynomial<T> (H.col(0))  ;
            H (1,1) =2;
            if (Degree==1) return Polynomial<T> (H.col(1));
            NumMethod::Cyclic<3> cyclic = 1;
            for (int i = 2; i<=Degree; i++){
                H(0,cyclic+1) = -2*(i-1.0)* H (0,cyclic-1);
                for (int j=1;j<=Degree;j++){
                    H(j,cyclic+1) = 2*(H (j-1,cyclic) - (i-1.0)*H(j,cyclic-1));
                }
                cyclic++;
                }
            
            return Polynomial<T> (H.col(cyclic));
            
            
        }
        
        //Use this for evaluation
        T operator () (T x) {
            T value [3]= {1,2*x,0};
            NumMethod::Cyclic<3> cyclic = 1;
            for (int i=1; i<n; i++){
                value [cyclic+1] =2*x*value [cyclic] - 2*i*value[cyclic-1];
                cyclic++;
            }
            return n > 1 ?  value [cyclic] : value [n];
        }
        
        void setDegree (int in){ n = in;}
        
        HermitePoly (int in): n(in) {}
    }; 
    
    
    
    
    
    
    
    
    
    
    
}


#endif	/* HERMITEPOLY_H */

