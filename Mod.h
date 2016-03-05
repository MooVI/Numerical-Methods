/* 
 * File:   mod.h
 * Author: qfeuille
 *
 * Created on 07 October 2012, 15:05
 */

#ifndef MOD_H
#define	MOD_H
#include <math.h>
namespace NumMethod {
inline int posmod (int num,int mod){
    return (num % mod + mod) % mod;
	}

//Look up table for + 1 mod N for when you have more memory than CPU speed
template <int N>
class add1mod {
    int table [N];
    public:
    int operator () (int i){return table [i];}
    const int * getTable (){return table; };
    add1mod(){
        
        for (int i = 0; i<N-1;i++){
        table [i]=i+1;
        }
        table [N-1]=0;    
    }
};


//Look up table for - 1 mod N for when you have more memory than CPU speed
template <int N>
class sub1mod {
    int table [N];
    public:
    int operator () (int i){return table [i];}
    const int * getTable (){return table; };
    sub1mod(){
        table [0]=N-1;
        for (int i = 1; i<N;i++){
        table [i]=i-1;
        }
        
    }
};

};

#endif	/* MOD_H */

