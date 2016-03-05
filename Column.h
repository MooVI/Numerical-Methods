/* 
 * File:   Column.h
 * Author: qfeuille
 *
 * Created on 07 October 2012, 15:09
 */

#ifndef COLUMN_H
#define	COLUMN_H
#include <iostream>
#include <fstream>
#include <math.h>
#include <assert.h>
namespace NumMethod{
    
    //
    // Copies the contents of two containers, from from to to.
    // The from container *must* implement .size() and they must
    // both overload [] to allow dereferencing.
    inline void cColumnWrite (const int& i, std::ofstream& file){}
    
       template <typename A,typename... Args>
    void cColumnWrite (const int& i, std::ofstream& file, const A& first,const Args&... rest){
          assert(i<first.size());
            file<<first[i]<<" ";
           cColumnWrite (i,file,rest...);
    }
    
    template <typename A,typename... Args>
    void cColumnWrite (std::ofstream& file, const A& first,const Args&... rest){
        for (int i=0;i<first.size();i++){
            file<<first[i]<<" ";
            cColumnWrite (i,file,rest...);
            file<<std::endl;   
        }
    }
    
    
    
    
 
    template <typename... Args>
    void cCount (int& size, const Args&... rest){
        size = sizeof...(Args);
    }
   
     
    template <typename A, typename B>
    void cCopy (const A& from, B& to) {
        for (int i=0;i<from.size();i++)
            to[i]=from[i];
    }
    
    template <typename A>
    auto cMax (const A& cont) -> decltype (cont[0]) {
        double max = cont[0];
        for (int i=1;i<cont.size();i++)
            if (cont[i]>max)
                max = cont[i];
        return max;          
    }
    template <typename A>
    auto cMin (const A& cont) -> decltype (cont[0]) {
        double min = cont[0];
        for (int i=1;i<cont.size();i++)
            if (cont[i]<min)
                min = cont[i];
        return min;          
    }
    
      template <typename A, typename T>
    void cMinMax (const A& cont,T& min,T& max ){
        min = cont[0];
        max = cont[0];
        for (int i=1;i<cont.size();i++){
            if (cont[i]>max)
                max = cont[i];
            else if (cont[i]<min)
                min = cont[i];
            
        }
      }
    
    
};

#endif	/* COLUMN_H */

