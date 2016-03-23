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

#include <cstdlib>
    
template<typename T>
struct ForLoopParams {
    int numPoints;
    T start;
    T end;
    
 };
 
 int default_argv_for_loop = {1,2,3};
 
 template<typename T>
 ForLoopParams<T> get_for_from_cmd(char** argv, int positions[3] = default_argv_for_loop){
     ForLoopParams<T> ret;
     ret.numPoints = atoi(argv[positions[0]]);
     ret.start = atof(argv[positions[1]]);
     ret.end = atof(argv[positions[2]]);
     return ret;
 }

 template<typename T>
 class GetXFor{
     std::vector<T> x;
     public:
 bool operator ()(T value, int i){
     x.push_back(value);
     return true;
 }
 
 std::vector<T> get_x(){
     return x;
 }
     
 };
 
 struct Range {
    template<typename FunctPtr>
    static void loop(FunctPtr& f, int begin, int exend, int step =1); 
    
    static inline std::vector<int> get_x(int begin, int exend, int step =1); 
 
 };
 
 template<typename FunctPtr>
 void Range::loop(FunctPtr& f, int begin, int exend, int step){
     ForLoopProgress<int> prog (begin, exend);
     for (int i =0, x=begin; x<exend ;i++, x+=step){
         prog(i);
         if (f (x,i)) break;        
     }
 }
 
 inline std::vector<int> Range::get_x(int begin, int exend, int step){
     std::vector<int> xs((exend-begin)/step);
     for (int i =0, x=begin; x<exend ;i++, x+=step){
         xs[i] = x;       
     }
     return xs;
 }
 
 
 
struct EqualSpaceFor {
    
    template<typename T,typename FunctPtr>
    static void loop(FunctPtr& f,const ForLoopParams<T>& params);
    
    template<typename T>
    static std::vector<T> get_x(const ForLoopParams<T>& params);
    
 };

 template<typename T>
 std::vector<T> EqualSpaceFor::get_x(const ForLoopParams<T>& params){
    std::vector<T> xvec (params.numPoints);        
     T step = (params.end-params.start)/((T) (params.numPoints -1));
     T x = params.start;
     for (int i =0;i<params.numPoints;i++){
         xvec[i] = x;
         x+=step;
     }
     return xvec;
 }
 

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

