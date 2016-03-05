/* 
 * File:   Cyclic.h
 * Author: qfeuille
 *
 * Created on 07 October 2012, 15:11
 */

#ifndef CYCLIC_H
#define	CYCLIC_H
namespace NumMethod{
 template<int N>
    class Cyclic{
       enum {P=(N-1)}; 
       int j;
    public:
        constexpr operator int () const{return j;}
        constexpr int operator ()() const {return j;};
        constexpr int operator +  (const int & i){
        return (i+j> P) ? (i+j-N):i+j;}
        constexpr int operator - (const int & i){
        return (j-i< 0) ? (N+j-i):j-i;}
        void operator ++ (int){
        j=*this+1;}
        void operator -- (int){
        j=*this-1;}
    
     void operator = (int in){
        j=in;}
     Cyclic (): j(0) {}
     Cyclic (const int& in): j(in){}
    };
}

#endif	/* CYCLIC_H */

