/* 
 * File:   newfile.h
 * Author: qfeuille
 *
 * Created on 07 October 2012, 15:03
 */

#ifndef FORLOOP_H
#define	FORLOOP_H
#include <math.h>
#include <iostream>
#include <time.h>
 namespace NumMethod {
template<typename T>
struct ForLoopProgress{
    
    T begin,range;
void operator () (T i){
    std::cout <<100*((i-begin)/(double)(range)) <<"%\r"<<std::flush;
}
ForLoopProgress (T ibegin,T iend): begin(ibegin),range(iend-ibegin){};
    
};



 template<typename T>
 struct FastForLoopProgress{
     int oldtime, newtime, count;
     T begin,range;
 void operator () (T i){
	 newtime = time(NULL);
	 if (newtime-oldtime>60){
     std::cout <<count<<": "<<100*((i-begin)/(double)(range)) <<"%\r"<<std::endl;
     count++;
     oldtime=newtime;
	 }
 }
 FastForLoopProgress (T ibegin,T iend): begin(ibegin),range(iend-ibegin){
	 oldtime=time(NULL);
     count=1;};

 };

 
 // I can't see a way of generalising this for arbitrary functions. Just copy
 // and paste and write needed code in apply.
#if 0 
 
  template<int J, int END>
    struct TemplateForLoop {

                inline static void apply() {
            TemplateForLoop<J + 1, END>::apply();
        }
    };
    template <int END>
    struct TemplateForLoop<END,END> {
        inline static void apply() {};
    };
    
#endif
 
 }
#endif	/* NEWFILE_H */

