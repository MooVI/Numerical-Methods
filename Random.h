/* 
 * File:   Random.h
 * Author: qfeuille
 *
 * Created on 15 September 2012, 23:52
 */

#ifndef RANDOM_H
#define	RANDOM_H



#include <time.h>
#include <math.h>
 #include "mtrand.h"
#include "sobol.hpp"
#include <gsl/gsl_qrng.h>


  

  

  
 
namespace NumMethod {
   

    
    template<int nDim>
    class Sobol {
        long long int lseed;
        int seed;
        gsl_qrng * q; 
        unsigned long getseed (){
      time_t seconds;
      time (&seconds); 
      return (unsigned long) seconds;
    } 
        public:
            void operator ()(double []);
            Sobol (){q = gsl_qrng_alloc (gsl_qrng_sobol, nDim);}
            ~Sobol(){gsl_qrng_free(q);}
    };
    
    template<int nDim> inline
    void Sobol<nDim>::operator ()(double out []){
        gsl_qrng_get(q, out);
    }
    
     
    
            
            
        
        
        
        
        
   
    class RandomBase{
    
protected:
    
    MTRand_int32 * gen;
    unsigned long seed (){
      time_t seconds;
      time (&seconds); 
      return (unsigned long) seconds;
    } 
    
public:
    double randomgen (double lowerbound, double upperbound){
        return (upperbound-lowerbound)*gen->get() + lowerbound;
    }
    
    double operator () (){return gen->get();}
    
    virtual ~RandomBase (){delete gen;}
    virtual void reseed ()=0;
};

class RandomOpen: public RandomBase{
public:
    RandomOpen (){gen= new MTRand_open (seed());}
    void reseed () { delete gen;
                        gen= new MTRand_open (seed());}
};

class RandomClosed: public RandomBase{
public:
    RandomClosed (){gen= new MTRand_closed (seed());}
    void reseed () { delete gen;
                        gen= new MTRand_closed (seed());}
};
    
    
    
//Takes random numbers [0,1) and converts them to various ranges
class OCRandomConversion {
    
    OCRandomConversion ();
public:
    //Closed Interval between min and max [min,max]
    static int Cint (double in, int min, int max);
    //Semi-Open Interval between min and max [min,max)
    static int OCint (double in, int min, int max);
    //Open Interval between min and max (min,max)
    static int Oint (double in, int min, int max);
    
    //Semi-Open Interval between min and max [min,max)
    static double OCdouble (double in, double min, double max);
  
    
    
};


 //Closed Interval between min and max [min,max]
    inline int OCRandomConversion::Cint (double in, int min, int max){
        return floor((max-min+1)*in+min);
    }
    //Semi-Open Interval between min and max [min,max)
   inline int OCRandomConversion::OCint (double in, int min, int max){
        return floor((max-min)*in+min);
    }
    //Open Interval between min and max (min,max)
    inline int OCRandomConversion::Oint (double in, int min, int max){
        return floor((max-min+1)*in+min);
    }
    
    //Semi-Open Interval between min and max [min,max)
    inline double OCRandomConversion::OCdouble (double in, double min, double max){
        return floor((max-min)*in+min);
    }
   
}
#endif	/* RANDOM_H */

