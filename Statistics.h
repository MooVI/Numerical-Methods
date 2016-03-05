/* 
 * File:   Statistics.h
 * Author: qfeuille
 *
 * Created on 17 September 2012, 11:18
 */

#include<vector>
#include<math.h>
#ifndef STATISTICS_H
#define	STATISTICS_H
namespace NumMethod{
    double dummy;
template<typename T>
struct Average {
template<typename FunctPtr, typename ...Args>
static auto average (FunctPtr& f, int numRuns,const Args&... args)->decltype(f(args...)){
    auto accumulate = f(args...);
    for (int i=1;i<numRuns;i++){       
        accumulate= accumulate + f(args...);
    }
    return  ( accumulate/ (T) numRuns);
}

// f must be of form auto f (T& weight, args)
template<typename FunctPtr, typename ...Args>
static auto weightedaverage (FunctPtr& f, int numRuns,const Args&... args)->decltype(f(dummy,args...)){
    T weight;
    auto accumulate = f(weight, args...)*weight;
    T norm=weight;
    for (int i=1;i<numRuns;i++){
        auto result= f(weight,args...);
        accumulate= accumulate + result*weight;
        norm+=weight;
    }
    return  ( accumulate/ norm);
}
};

template <typename T>
class RunningStats
{
public:
    RunningStats();
    template<typename A>
    RunningStats(A array);
    void Clear();
    void Push(T x);
    long long NumDataValues() const;
    T Mean() const;
    T Variance() const;
    T StandardDeviation() const;
    T Skewness() const;
    T Kurtosis() const;
    template<class S>
    friend RunningStats<S> operator+(const RunningStats<S> a, const RunningStats<S> b);
    RunningStats<T>& operator+=(const RunningStats<T> &rhs);

private:
    long long n;
    T M1, M2, M3, M4, nd;
};

template <typename T>
RunningStats<T>::RunningStats() 
{
    Clear();
}

template <typename T>
template<typename A>
RunningStats<T>::RunningStats(A array) 
{
    Clear();
    for (int i=0; i<array.size();i++) this->Push(array[i]);
}

template <typename T>
void RunningStats<T>::Clear()
{
    n = 0;
    nd = 0;
    M1 = M2 = M3 = M4 = 0.0;
}

template <typename T>
void RunningStats<T>::Push(T x)
{
    T delta, delta_n, delta_n2, term1;
    T n1 = nd;
    n++;
    nd = (T) ((double) n);
    delta = x - M1;
    delta_n = delta / nd;
    delta_n2 = delta_n * delta_n;
    term1 = delta * delta_n * n1;
    M1 += delta_n;
    M4 += term1 * delta_n2 * (nd*nd - 3*nd + 3) + 6 * delta_n2 * M2 - 4 * delta_n * M3;
    M3 += term1 * delta_n * (nd - 2) - 3 * delta_n * M2;
    M2 += term1;
}

template <typename T>
long long RunningStats<T>::NumDataValues() const
{
    return n;
}

template <typename T>
T RunningStats<T>::Mean() const
{
    return M1;
}

template <typename T>
T RunningStats<T>::Variance() const
{
    return M2/(nd-1.0);
}

template <typename T>
T RunningStats<T>::StandardDeviation() const
{
    return sqrt( Variance() );
}

template <typename T>
T RunningStats<T>::Skewness() const
{
    return sqrt(nd) * M3/ pow(M2, 1.5);
}

template <typename T>
T RunningStats<T>::Kurtosis() const
{
    return nd*M4 / (M2*M2) - 3.0;
}

template <typename T>
RunningStats<T> operator+(const RunningStats<T> a, const RunningStats<T> b)
{
    RunningStats<T> combined;
    
    combined.n = a.n + b.n;
    combined.nd = a.nd +b.nd;
    
    T delta = b.M1 - a.M1;
    T delta2 = delta*delta;
    T delta3 = delta*delta2;
    T delta4 = delta2*delta2;
    
    combined.M1 = (a.n*a.M1 + b.n*b.M1) / combined.n;
    
    combined.M2 = a.M2 + b.M2 + 
                  delta2 * a.n * b.n / combined.n;
    
    combined.M3 = a.M3 + b.M3 + 
                  delta3 * a.n * b.n * (a.n - b.n)/(combined.n*combined.n);
    combined.M3 += 3.0*delta * (a.n*b.M2 - b.n*a.M2) / combined.n;
    
    combined.M4 = a.M4 + b.M4 + delta4*a.n*b.n * (a.n*a.n - a.n*b.n + b.n*b.n) / 
                  (combined.n*combined.n*combined.n);
    combined.M4 += 6.0*delta2 * (a.n*a.n*b.M2 + b.n*b.n*a.M2)/(combined.n*combined.n) + 
                  4.0*delta*(a.n*b.M3 - b.n*a.M3) / combined.n;
    
    return combined;
}

template <typename T>
RunningStats<T>& RunningStats<T>::operator+=(const RunningStats<T>& rhs)
{ 
	RunningStats<T> combined = *this + rhs;
	*this = combined;
	return *this;
}


}
#endif	/* STATISTICS_H */
