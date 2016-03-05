/* 
 * File:   Polynomial.h
 * Author: qfeuille
 *
 * Created on 22 March 2013, 17:07
 */
#include <Eigen/Dense>
#include <ostream>
#ifndef POLYNOMIAL_H
#define	POLYNOMIAL_H
namespace NumMethod {
    
    //This class is to be used for fast, convenient evaluation of polynomials
    //It is not to be used for CAS purposes
    template <typename T>
    class Polynomial {
       
    protected:
        Eigen::Matrix<T, Eigen::Dynamic, 1 > coeff;
    public:

        void set (T incoeff, int index);

        T safeget (int index) const ;

        T get (int index) const;
        
       Eigen::Matrix<T, Eigen::Dynamic, 1 > getCoeff() const;
         
       template < typename A>
        Polynomial(const Eigen::MatrixBase<A> & icoeff);
       
        Polynomial <T> operator * (T scal) const;
        
        

        //Via Horner's Method

        T evaluate(double x) const;
        
        T operator () (double x) const;
        
        

        Polynomial(int n=1) ;
        
        template<typename A>
        friend Polynomial <A> operator * (A scal,Polynomial <A> in);
         
        template<typename Y>
       friend std::ostream& operator<<(std::ostream& os, const Polynomial<Y>& poly);

    };
   
    
     template<typename T>
   std::ostream& operator<<(std::ostream& os, const Polynomial<T>& poly){     
      os << poly.coeff;
      return os; 
   }
    
    template<typename T>
    void Polynomial<T>::set(T incoeff, int index) {
            coeff[index] = incoeff;
        }

        template<typename T>
        T  Polynomial<T>::safeget(int index) const {
            return (coeff.size()-1)>= index ? coeff[index] : (T) 0.0;
        }
        
   
        template<typename T>
        T  Polynomial<T>::get(int index) const {
            return coeff[index];
        }
        
        template<typename T>
       Eigen::Matrix<T, Eigen::Dynamic, 1 >  Polynomial<T>::getCoeff() const {return coeff;}
         
        template<typename T>
       template < typename A>
         Polynomial<T>::Polynomial(const Eigen::MatrixBase<A> & icoeff) : coeff(icoeff) {
        };

       template<typename T>
        Polynomial <T>  Polynomial<T>::operator * (T scal) const {
            
            return Polynomial<T> (scal*this->coeff);
             }
        
        

        //Via Horner's Method
template<typename T>
        T  Polynomial<T>::evaluate(double x) const{
            double ret = 0;
            for (int i = coeff.size()-1; i >= 0; i--)
                ret = coeff[i] + x * ret;
            return ret;
        }
        template<typename T>
        T  Polynomial<T>::operator () (double x) const{
            double ret = 0;
            for (int i = coeff.size()-1; i >= 0; i--)
                ret = coeff[i] + x * ret;
            return ret;
            
        }
        
   
    
    template <typename T>
    Polynomial<T>::Polynomial(int n){
        coeff = Eigen::Matrix<T, Eigen::Dynamic,1>::Zero(n);
        
    }
    
    
    
    template<typename T>
    Polynomial <T> operator * (T scal,Polynomial <T> in) {
            
            return Polynomial<T> (scal*in.coeff);
             }
}





#endif	/* POLYNOMIAL_H */

