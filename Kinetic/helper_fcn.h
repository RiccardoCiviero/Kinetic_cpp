#pragma once
#include<Eigen>


template <typename T> inline int sgn(T val) {
	return (T(0) < val) - (val < T(0));
}

//Eigen::ArrayXd cumsum(Eigen::ArrayXd a) {
//	double temp = 0;
//	double temp2 = 0;
//	for (auto i = 0; i < a.rows(); i++) {
//		temp2 = a(i);
//		a(i) += temp;
//		temp += temp2;
//	}
//	return a;
//}

template <typename Derived> inline
Eigen::Array<typename Derived::Scalar, Eigen::Dynamic, 1> cumsum(const Eigen::ArrayBase<Derived>& array)
{
    using Scalar = typename Derived::Scalar;
    Eigen::Array<Scalar, Eigen::Dynamic, 1> result(array.size());

    if (array.size() == 0)
        return result;

    result(0) = array(0);
    for (int i = 1; i < array.size(); ++i)
    {
        result(i) = result(i - 1) + array(i);
    }

    return result;
}
