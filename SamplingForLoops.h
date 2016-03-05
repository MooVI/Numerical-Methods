/* 
 * File:   SamplingForLoops.h
 * Author: qfeuille
 *
 * Created on 25 October 2012, 16:02
 */


#include "ForLoop.h"

#ifndef SAMPLINGFORLOOPS_H
#define	SAMPLINGFORLOOPS_H

namespace NumMethod{

template<typename T>
struct ForLoopParams {
    int numPoints;
    T start;
    T end;
    
 };

struct EqualSpaceFor {
    
    template<typename T,typename FunctPtr>
    static void loop(FunctPtr& f,const ForLoopParams<T>& params);
    
    
    
 };

 template<typename T,typename FunctPtr>
 void EqualSpaceFor::loop(FunctPtr& f,const ForLoopParams<T>& params){
     ForLoopProgress<int> prog (0,params.numPoints);
     T step = (params.end-params.start)/((T) (params.numPoints -1));
     T x = params.start;
     for (int i =0;i<params.numPoints;i++){
         prog(i);
          if (f (x,i)) break;   
         x+=step;     
     }
 }
 
 
 
struct LogFor {
    
    template<typename T,typename FunctPtr>
    static void loop(FunctPtr& f,const ForLoopParams<T>& params);
    
    
    
 };

 template<typename T,typename FunctPtr>
  void LogFor::loop(FunctPtr& f,const ForLoopParams<T>& params){
     ForLoopProgress<int> prog (0,params.numPoints);
     T powdiff = log10 (params.end/params.start);
     T powgap = pow(10,powdiff/(params.numPoints-1));
     T x = params.start;
     for (int i =0;i<params.numPoints;i++){
         prog(i);
          if (f (x,i)) break;   
         x*=powgap;     
     }
 }
 
 struct CubicFor {
    
    template<typename T,typename FunctPtr>
    static void loop(FunctPtr& f,const ForLoopParams<T>& params);
    
    
    
 };

 template<typename T,typename FunctPtr>
 void CubicFor::loop(FunctPtr& f,const ForLoopParams<T>& params){
     ForLoopProgress<int> prog (0,params.numPoints);
     T x = params.start;
     for (int i =0;i<params.numPoints;i++){
         x = (pow (( (2)*( (double) i/ (double) (params.numPoints-1) -0.5) ),3)
             +1)
             /2*(params.end-params.start)
             +params.start;
         prog(i);
         if (f (x,i)) break;    
     }
 }

 
 
}
 


#endif	/* SAMPLINGFORLOOPS_H */

