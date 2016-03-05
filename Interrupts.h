/* 
 * File:   Interrupts.h
 * Author: qfeuille
 *
 * Created on December 29, 2011, 11:29 PM
 */
#include "StdVectorMath.h"
#include "Cyclic.h"
#include<vector>
#include <Eigen/Dense>
#include<Eigen/StdVector>



#ifndef INTERRUPTS_H
#define	INTERRUPTS_H
using std::vector;
//typedef  std::function<bool (const T& dq,const Eigen::Matrix <T,N,1>& q,const T& stepSize)> Interruptx;

namespace NumMethod
{


class DoNothing {
public:
    template<typename T, int N> 
    bool operator () (const T& dq,const Eigen::Matrix <T,N,1>& q,const T& stepSize){return false;};
};


template<typename T,int N, typename Interrupt=DoNothing>
class SamplePerStep{
    vector<Eigen::Matrix<T,N,1>,
            Eigen::aligned_allocator<Eigen::Matrix<T,N,1> > > storeq;
    vector<T> storedq;
public:
    Interrupt interrupt;
    Eigen::Array<T,Eigen::Dynamic,1>  getq (int ident);   
    Eigen::Array<T,Eigen::Dynamic,1> getdq ();
    void clear() {
        storeq.resize(0);
        storedq.resize(0);
    };
    bool operator () (const T& dq,const Eigen::Matrix <T,N,1>& q,const T& stepSize){
        storeq.push_back(q);
        storedq.push_back(dq);
        
        return interrupt (dq,q,stepSize);};
        template<typename ...Args>
        SamplePerStep (int numSteps,const Args& ...args);
        SamplePerStep (int numSteps=100000);
        ~SamplePerStep();
};


template<typename T,int N,typename Interrupt>
template<typename ...Args>
SamplePerStep<T,N,Interrupt>::SamplePerStep (int numSteps,const Args& ...args):
interrupt (args ...){
            storeq.reserve(numSteps);
            storedq.reserve(numSteps);
        }

template<typename T,int N,typename Interrupt>
SamplePerStep<T,N,Interrupt>::SamplePerStep (int numSteps):
interrupt (){
            storeq.reserve(numSteps);
            storedq.reserve(numSteps);
        }

template<typename T,int N,typename Interrupt>
Eigen::Array<T,Eigen::Dynamic,1> SamplePerStep<T,N,Interrupt>::getdq(){
    Eigen::Array<T,Eigen::Dynamic,1> ret (storedq.size());
        NumMethod::cCopy(storedq,ret);
        return ret;
}

template<typename T,int N, typename Interrupt>
SamplePerStep<T,N,Interrupt>::~SamplePerStep(){}

template<typename T,int N,typename Interrupt>
Eigen::Array<T,Eigen::Dynamic,1> SamplePerStep<T,N,Interrupt>::getq(int ident){
        Eigen::Array<T,Eigen::Dynamic,1> ret(storeq.size());
        assert (ident<storeq[0].size());
        for (int i=0;i<storeq.size();++i){
            ret [i] = (storeq [i] [ident]);
        }
        return ret;
        }

template<typename T,int N,typename Interrupt=DoNothing>
class FixedSample{
    vector<Eigen::Matrix<T,N,1>,
            Eigen::aligned_allocator<Eigen::Matrix<T,N,1> > > storeq;
    vector<T> storedq;
    T _period;
    T _periodstore;
    double orient;
   
public:
    Interrupt interrupt;
    Eigen::Array<T,Eigen::Dynamic,1>  getq (int ident);
    Eigen::Array<T,Eigen::Dynamic,1> getdq ();
    Eigen::Array<T,Eigen::Dynamic,1> gett ();
    bool operator () (const T& dq,const Eigen::Matrix <T,N,1>& q,const T& stepSize){
        if ((_periodstore-=orient*stepSize) <= 0){
        storeq.push_back(q);
        storedq.push_back(dq);
        _periodstore=_period;
        }
        
        return interrupt(dq,q,stepSize);};
        
        void clear(int iorient=0);
        const int getSize () const  {return storeq.size();}
        void setOrient (int iorient) { orient = iorient;}
        
        
        template<typename ...Args>
        FixedSample (T period, int numSteps, const Args&... args);
        
        ~FixedSample();
};

template<typename T,int N,typename Interrupt>
void FixedSample<T,N,Interrupt>::clear(int iorient){
    storeq.resize(0);
    storedq.resize(0);
    _periodstore=0;
    if (iorient!=0)
        orient = iorient;

}


template<typename T,int N,typename Interrupt>
template<typename ...Args>
FixedSample<T,N,Interrupt>::FixedSample (T period, int numSteps,const Args&... args):
_period(period),_periodstore(0),interrupt (args...){
            if (numSteps>100000) 
                numSteps=100000;
            storeq.reserve(numSteps);
            storedq.reserve(numSteps);
            if (period<0){
                period*=-1;
                orient = -1;
            }
            else orient =1;    
        }


template<typename T,int N,typename Interrupt>
FixedSample<T,N,Interrupt>::~FixedSample(){}

template<typename T,int N,typename Interrupt>
Eigen::Array<T,Eigen::Dynamic,1> FixedSample<T,N,Interrupt>::getdq(){
    Eigen::Array<T,Eigen::Dynamic,1> ret (storedq.size());
        NumMethod::cCopy(storedq,ret);
        return ret;
}

template<typename T,int N,typename Interrupt>
Eigen::Array<T,Eigen::Dynamic,1> FixedSample<T,N,Interrupt>::gett(){
    Eigen::Array<T,Eigen::Dynamic,1> ret (storedq.size());
        NumMethod::cCopy(storedq,ret);
        return ret;
}

template<typename T, int N,typename Interrupt>
Eigen::Array<T,Eigen::Dynamic,1> FixedSample<T,N,Interrupt>::getq(int ident){
        Eigen::Array<T,Eigen::Dynamic,1> ret (storeq.size());
        assert (ident<storeq[0].size());
        for (int i=0;i<storeq.size();++i){
            ret[i]=storeq [i] [ident];
        }
        return ret;
        }
    


struct StopFirst
  {
 bool first;
  public:
      
    template<typename T, int N> 
    bool operator () (const T& dq,const Eigen::Matrix <T,N,1>& q,const T& stepSize){
        if (first){
            first=false;
            return false;
        }
        return true;}
        
       
        StopFirst(){first=true;}
  };

  template<typename T>
  struct SampleStepSize {
      vector<T> storestepSize;
  public:
      template<int N> 
    bool operator () (const T& dq,const Eigen::Matrix <T,N,1>& q,const T& stepSize){
          storestepSize.push_back(stepSize);
          return false;
      }
      Eigen::Array<T,Eigen::Dynamic,1> getStepSize ();
  };
  
 template<typename T>
Eigen::Array<T,Eigen::Dynamic,1> SampleStepSize<T>::getStepSize(){
    Eigen::Array<T,Eigen::Dynamic,1> ret (storestepSize.size());
        NumMethod::cCopy(storestepSize,ret);
        return ret;
}

 /*
template<class F1,class F2>
class JoinIntHelper{
    F1* f1; F2 * f2;
public:
    template<typename T,int N>
    bool operator () (const T& dq,const Eigen::Matrix <T,N,1>& q,const T& stepSize){
        if (*f1(dq,q,stepSize)||*f2(dq,q,stepSize)) return true;
}
       JoinIntHelper (F1& if1,F2& if2){f1=&if1;f2=&if2;}
};

template<typename T,i>
struct JoinInt {
    
   template <class F1,class F2>
   std::function<bool (const T& dq,const Eigen::Matrix <T,N,1>& q,const T& stepSize)> operator () (F1 f1, F2 f2){
       return [&](const T& dq,const Eigen::Matrix <T,N,1>& q,const T& stepSize){
           if (f1(dq,q,stepSize)||f2(dq,q,stepSize)) return true;
       }; 
};
};
*/
 
 //The second interrupt cannot have any constructor arguments!
 template<class Interrupt1,class Interrupt2>
 struct JoinInt {
   Interrupt1 interrupt1; Interrupt2 interrupt2;
   template<int N,typename T> 
    bool operator () (const T& dq,Eigen::Matrix <T,N,1>& q,const T& stepSize){      
       return (interrupt1 (dq,q,stepSize) || interrupt2 (dq,q,stepSize));
 }
   
   template<typename ...Args>
   JoinInt (Args& ... args): interrupt1 (args...) ,interrupt1(){};
     
     
   
 };
 
}

#endif	/* INTERRUPTS_H */

