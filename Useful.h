/* 
 * File:   Useful.h
 * Author: qfeuille
 *
 * Created on 07 October 2012, 15:07
 */
#include <math.h>
#include <algorithm>
#include <array>

#ifndef USEFUL_H
#define	USEFUL_H
namespace NumMethod{
  
   
    
    template <typename T>
    inline T signmaxAbs (const T& x, const T& y){
       return (x < y ^ x > 0 ? x: y);
    }
    template <typename T>
    inline T signminAbs (const T& x, const T& y){
       return (x < y ^ x > 0 ? y: x);
    }
    
    template<typename T>
    inline T maxAbs (const T& x, const T& y){
        return fabs(x) > fabs(y) ? x: y;
    }
    
     template<typename T>
    inline T minAbs (const T& x, const T& y){
        return fabs(x) < fabs(y) ? x: y;
    }
     
      template<typename T>
    inline T max (const T& x, const T& y){
        return x > y ? x: y;
    }
    
     template<typename T>
    inline T min (const T& x, const T& y){
        return x < y ? x: y;
    }
     
       template<typename T>
    inline const T sign (const T& x, const T& y){
        return y >= 0.0 ? (x>=0.0 ? x: -x) : (x>=0.0 ? -x:x);
    }
     
     
inline  double roundtodp (double n,double decimalplace){
	decimalplace++;
	double multiple=pow(10,decimalplace);
	int x=(int) ceil(n*multiple);
	int remainder=x%10;
	double ret=n*multiple/10;
	if (remainder>4) n=(ceil(ret));
	else n=(floor(ret));
	n=n*10/multiple;
		return (n);}





template <int N>
struct NearestPowerof2
{
    enum { val = NearestPowerof2<N/2>::val * 2 };
};

template <>
struct NearestPowerof2<0>
{
    enum { val = 1 };
    
};


template<typename A, typename B>
std::array<int,2> find_bounding_indices(const B& val, const A& array, const int& begin, const int& end){
  if(val == array[begin]) return {{begin, begin+1}};
  int d = end-begin;  
  if(d==1) return {{begin, begin+1}};
  int centre = int(begin + floor(d/2));
  auto centreval = array[centre];
  if(val < centreval) return find_bounding_indices(val, array, begin, centre);
  return find_bounding_indices(val, array, centre, end);
};  


}

#endif	/* USEFUL_H */

