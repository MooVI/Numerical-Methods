#ifndef STATISTICS_H
#define	STATISTICS_H
namespace NumMethod{


#include <algorithm>
#include <iostream>
#include <numeric>
#include <vector>


template <typename T>
std::vector<T> linear_fit(const std::vector<T>& x, const std::vector<T>& y) {
    const auto n    = x.size();
    const auto s_x  = std::accumulate(x.begin(), x.end(), 0.0);
    const auto s_y  = std::accumulate(y.begin(), y.end(), 0.0);
    const auto s_xx = std::inner_product(x.begin(), x.end(), x.begin(), 0.0);
    const auto s_xy = std::inner_product(x.begin(), x.end(), y.begin(), 0.0);
    const auto slope    = (n * s_xy - s_x * s_y) / (n * s_xx - s_x * s_x);
    const auto intercept = (s_y - slope * s_x) / n;
    return {slope, intercept};
}

}

#endif
