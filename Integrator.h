/* 
 * File:   Integrator.h
 * Author: qfeuille
 *
 * Created on 24 May 2012, 01:08
 */


#ifndef INTEGRATOR_H
#define	INTEGRATOR_H
#include "PolynomialInterpolation.h"

namespace NumMethod{
template<typename T>
struct TrapeziumRule {
    
    protected:
    int numpoints=0;
    T range;
    T upper;
    T lower;
    T integral;

    public:
        
    template <typename FunctPtr >
            T operator () (FunctPtr & f) {
        if (numpoints == 0){
            numpoints++;
            return integral = (f(lower) + f(upper))*0.5 * (range);
        }
        else {
            T spacing= range/ (T)numpoints;
            T x = lower + 0.5* spacing;
            T sum = 0;
            for (int i=0;i<numpoints;i++,x+=spacing)
                sum += f(x);
            integral= 0.5*(integral+ range*sum/ (T) numpoints);
            numpoints<<=1;
            return integral;
            
        }
    }
    

    void setUp(const T ilower, const T iupper) {
        upper = iupper;
        lower = ilower;
        range = upper - lower;
        numpoints=0;
    }

    T getStepSize () const {return range/((T) numpoints);}
};



template<typename T>
struct MidPoint {

    protected:
    int numpoints=0;
    T range;
    T upper;
    T lower;
    T integral;

    public:
    template <typename FunctPtr >
            T operator () (FunctPtr & f) {
        if (numpoints == 1){
            numpoints++;
            return integral = f(0.5*(lower+upper))*range;
        }
        else {
            int add = pow (3,numpoints-1);
            T spacing= range/ (3*(T)add);
            T doubleSpacing = spacing+spacing;
            T x = lower + 0.5* spacing;
            T sum = 0;
            for (int i=0;i<add;i++,x+=spacing){
               sum += f(x);
               x+=doubleSpacing;
               sum +=f(x);
            }
            numpoints++;
            return integral=(integral+ range*sum/ (T) add)/3.0;
            }      
        }
    

   void setUp(const T ilower, const T iupper) {
        upper = iupper;
        lower = ilower;
        range = upper - lower;
        numpoints=1;
    }

    T getStepSize () const {return range/((T) numpoints);}
};

template<typename T>
class ChangeVariableMidPoint {
protected:
    int numpoints=0;
    T range;
    T upper;
    T lower;
    T integral;

    public:
        
    template <typename FunctPtr >   
    inline T changevar (FunctPtr& f, T var) {
        return f(var);     
    }
    template <typename FunctPtr >
     T operator () (FunctPtr & f) {
        if (numpoints == 1){
            numpoints++;
            return integral = changevar(f,0.5*(lower+upper))*range;
        }
        else {
            int add = pow (3,numpoints-1);
            T spacing= range/ (3*(T)add);
            T doubleSpacing = spacing+spacing;
            T x = lower + 0.5* spacing;
            T sum = 0;
            for (int i=0;i<add;i++,x+=spacing){
               sum += changevar(f,x);
               x+=doubleSpacing;
               sum +=changevar(f,x);
            }
            numpoints++;
            return integral=(integral+ range*sum/ (T) add)/3.0;
            }      
        }
   

   virtual void setUp(const T ilower, const T iupper) {
        upper = iupper;
        lower = ilower;
        range = upper - lower;
        numpoints=1;
    }

    T getStepSize () const {return range/((T) numpoints);}
};

template<typename T>
class SemiInfiniteMidPoint{
    protected:
    int numpoints=0;
    T range;
    T upper;
    T lower;
    T integral;

    public:
        
     template <typename FunctPtr >   
    inline T changevar (FunctPtr& f, T var) {
        return f(1/var)/ (var*var);    
     }
        
    template <typename FunctPtr >
     T operator () (FunctPtr & f) {
        if (numpoints == 1){
            numpoints++;
            return integral = changevar(f,0.5*(lower+upper))*range;
        }
        else {
            int add = pow (3,numpoints-1);
            T spacing= range/ (3*(T)add);
            T doubleSpacing = spacing+spacing;
            T x = lower + 0.5* spacing;
            T sum = 0;
            for (int i=0;i<add;i++,x+=spacing){
               sum += changevar(f,x);
               x+=doubleSpacing;
               sum +=changevar(f,x);
            }
            numpoints++;
            return integral=(integral+ range*sum/ (T) add)/3.0;
            }      
        }
    
    void setUp(const T ilower, const T iupper) {
            assert(iupper * ilower > 0.0);
            if (iupper < ilower) {
                upper = 1 / iupper;
                lower = 1 / ilower;
            } else {
                upper = 1 / ilower;
                lower = 1 / iupper;
            }
            range = upper - lower;
            numpoints = 1;
        }

    };


template<typename T>
class ExpMidPoint{
    protected:
    int numpoints=0;
    T range;
    T upper;
    T lower;
    T integral;

    public:
        
         template <typename FunctPtr >   
    inline T changevar (FunctPtr& f, T var) {
        return f(-log(var))/ var;   
    }   
        
   template <typename FunctPtr >
     T operator () (FunctPtr & f) {
        if (numpoints == 1){
            numpoints++;
            return integral = changevar(f,0.5*(lower+upper))*range;
        }
        else {
            int add = pow (3,numpoints-1);
            T spacing= range/ (3*(T)add);
            T doubleSpacing = spacing+spacing;
            T x = lower + 0.5* spacing;
            T sum = 0;
            for (int i=0;i<add;i++,x+=spacing){
               sum += changevar(f,x);
               x+=doubleSpacing;
               sum +=changevar(f,x);
            }
            numpoints++;
            return integral=(integral+ range*sum/ (T) add)/3.0;
            }      
        }
   
    
    void setUp(const T ilower, const T iupper){
        upper = exp(-ilower);
        lower = 0;
        range = upper - lower;
        numpoints=1;
    }

};

template<typename T>
struct ExpSqMidPoint {
    int sign;
    protected:
    int numpoints=0;
    T range;
    T upper;
    T lower;
    T integral;

    public:
        template <typename FunctPtr >   
    inline T changevar (FunctPtr& f, T var) {
            double inter = sqrt(-log(var));
        return 0.5*f(sign*inter)/ (var*inter);   
    }   
        
    template <typename FunctPtr >
     T operator () (FunctPtr & f) {
        if (numpoints == 1){
            numpoints++;
            return integral = changevar(f,0.5*(lower+upper))*range;
        }
        else {
            int add = pow (3,numpoints-1);
            T spacing= range/ (3*(T)add);
            T doubleSpacing = spacing+spacing;
            T x = lower + 0.5* spacing;
            T sum = 0;
            for (int i=0;i<add;i++,x+=spacing){
               sum += changevar(f,x);
               x+=doubleSpacing;
               sum +=changevar(f,x);
            }
            numpoints++;
            return integral=(integral+ range*sum/ (T) add)/3.0;
            }      
        }
    
    
   void setUp(const T ilower, const T iupper){
        assert (iupper*ilower>=0.0);
        if (fabs(iupper)>fabs(ilower))
        {
        upper = exp(-ilower*ilower);
        lower = 0;
        sign = 1;
        }
        else {
        upper = exp(-iupper*iupper);
        lower = 0;
        sign = 1;
        }
        range = upper - lower;
        numpoints=1;
    }

};




template<typename T>
class TrapeziumMethod {
    TrapeziumRule<T> trap;
    T prec;
    int maxSteps; //Where actual number of function evaluations is 2^maxsteps
    public:
        template <typename FunctPtr >
    T operator ()(FunctPtr f, T lower, T upper) {
            T value,oldvalue;
            trap.setUp(lower,upper);
            for (int i=0;i<maxSteps;i++){
                value = trap(f);
                if (i>5)
                    if (fabs(value-oldvalue)<prec*fabs(oldvalue)
                            || (value==0 && oldvalue==0))
                        return value;
                oldvalue=value;
            }
                std::cerr<<"Trapezium: Too many steps to reach required precision."<<std::endl;
                return NAN;
                
            };
            
            T getStepSize () const {return trap.getStepSize();}

    TrapeziumMethod(T iprec, int imaxSteps = 20) : prec(iprec), maxSteps(imaxSteps) {
    };
        
    
};


template<typename T>
class SimpsonsMethod {
    TrapeziumRule<T> trap;
    T prec;
    int maxSteps; //Where actual number of function evaluations is 2^maxsteps
    public:
        template <typename FunctPtr >
    T integrate(FunctPtr f, T lower, T upper) {
            T value,oldvalue,temp,oldtemp;
            trap.setUp(lower,upper);
            for (int i=0;i<maxSteps;i++){
                temp = trap(f);
                value = (4.0*temp-oldtemp)/3.0;
                if (i>5)
                    if (fabs(value-oldvalue)<prec*fabs(oldvalue)
                            || (value==0 && oldvalue==0))
                        return value;
                oldvalue=value;
                oldtemp=temp;
            }
                std::cerr<<"Simpson: Too many steps to reach required precision."<<std::endl;
                return NAN;
                
            };

    SimpsonsMethod(T iprec, int imaxSteps = 20) : prec(iprec), maxSteps(imaxSteps) {
    };
        
    
};

    template<typename T, int Order = 5, int maxSteps = 20 >
    class Romberg {
        FixedPolynomialInterpolation<T, Order> interpolate;
        TrapeziumRule<T> trap;
        T prec;
        //Where actual number of function evaluations is 2^maxsteps
    public:

        template <typename FunctPtr >
        T operator()(FunctPtr f, T lower, T upper) {
            Eigen::Matrix<T, maxSteps, 1> values;
            Eigen::Matrix<T, maxSteps + 1, 1 > steps;
            T out, error;
            trap.setUp(lower, upper);
            steps[0] = 1;
            for (int i = 0; i < maxSteps; i++) {
                values[i] = trap(f);
                if (i > Order) {
                    out = interpolate(steps. template segment<Order> (i - Order),
                            values. template segment<Order> (i - Order),
                            0.0, error);
                    if (fabs(error) <= prec * fabs(out))return out;
                }
                steps[i + 1] = 0.25 * steps[i];
            }
            std::cerr << "Romberg: Too many steps to reach required precision." << std::endl;
            return out;
        }

        Romberg(T iprec) : prec(iprec) {
        }
    };

   
    
 template<typename T, int Order = 5, int maxSteps = 20>
    class RombergOpen {
        FixedPolynomialInterpolation<T, Order> interpolate;
        T prec;
        //Where actual number of function evaluations is 2^maxsteps
    public:

        template <typename FunctPtr, typename Method>
        T operator()(FunctPtr& f, Method& method, T lower, T upper) {
            Eigen::Matrix<T, maxSteps, 1> values;
            Eigen::Matrix<T, maxSteps + 1, 1 > steps;
            T out, error;
            method.setUp(lower, upper);
            steps[0] = 1;
            for (int i = 0; i < maxSteps; i++) {
                values[i] = method(f);
                if (i > Order) {
                    out = interpolate(steps. template segment<Order> (i - Order),
                            values. template segment<Order> (i - Order),
                            0.0, error);
                    if (fabs(error) <= prec * fabs(out))return out;
                }
                steps[i + 1] = steps[i] / 9.0;
            }
            std::cerr << "Romberg Open Too many steps to reach required precision." << std::endl;
            return out;

        }

        RombergOpen(T iprec) : prec(iprec) {
        }
    };

}

#endif	/* INTEGRATOR_H */

