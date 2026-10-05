// vim: fdm=syntax:ts=4:sw=4:syntax=c.doxygen:tw=70:fo+=Mm

#ifndef __Transform_h_
#define __Transform_h_

#include <vector>
#include <cmath>
#include <cstdlib>
#include <limits>
#include <map>

#include <Eigen/Dense>

template <class T>
class Transform {
    private:
        /**
         * 积分点。
         */
        static std::map<unsigned int, std::vector<T> > quad_point;

        /**
         * Hermite 变换矩阵。
         */
        static std::map<unsigned int, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > trans_matrix;

        /**
         * Hermite 逆变换矩阵。
         */
        static std::map<unsigned int, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > inv_trans_matrix;

        /**
         * 生成 @p n 乘 @p n 的 Hermite 变换矩阵。
         */
        static void _GenMatrix(unsigned int n) {
            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix = trans_matrix[n];
            matrix.resize(n, n);
            std::vector<T> quad_point;
            HermiteQuadPoint(n, quad_point);

            for (unsigned int i = 0; i < n; i++)
                for (unsigned int j = 0; j < n; j++) {
                    if (j == 0)
                        matrix(i, j) = 1;
                    else if (j == 1)
                        matrix(i, j) = quad_point[i];
                    else
                        matrix(i, j) = quad_point[i] * matrix(i, j-1) - (j - 1) * matrix(i, j-2);
                }
        }

        /**
         * 生成 @p n 乘 @p n 的 Hermite 逆变换矩阵。
         */
        static void _GenInvMatrix(unsigned int n) {
            if (trans_matrix.find(n) == trans_matrix.end())
                _GenMatrix(n);

            Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>& matrix = inv_trans_matrix[n];
            matrix = trans_matrix[n].inverse();
        }

    public:
        /**
         * 求 @p __deg 次 Hermite 多项式的最大根。
         */
        static T MaxRootOfHermitePolynomial(unsigned int __deg)
        {
            static std::map<unsigned int, T> __roots;
            typename std::map<unsigned int, T>::iterator __it;
            if ((__it = __roots.find(__deg)) != __roots.end()) 
                return __it->second;

            T __res, __r1, __N = 2 * __deg + 1;
            static const T __pipm4 = sqrt(M_2_SQRTPI * 0.5);
            __res = sqrt(__N) - 1.85575 * pow(__N, -1.0 / 6); 
            do {
                unsigned int __j;
                double __p1, __p2 = 0, __p3 = __pipm4;
                for (__j = 0; __j <= __deg - 1; __j++) {
                    __p1 = __p2;
                    __p2 = __p3;
                    __p3 = __p2 * __res * sqrt(2.0 / (__j + 1))
                        - __p1 * sqrt(__j / (__j + 1.0));
                }
                double __dp3 = sqrt(2.0 * __j) * __p2;
                __r1 = __res;
                __res -= __p3 / __dp3;
            } while (fabs(__res - __r1) >= 100 * std::numeric_limits<T>::epsilon() * (1 + fabs(__res)));
            __res *= M_SQRT2;
            __roots[__deg] = __res;
            return __res;
        }

        /**
         * 求 Gauss 积分点。积分点的个数由 @p n 给出，计算出的积分点存在数组
         * @p x 中。
         */
        static void HermiteQuadPoint(unsigned int n, std::vector<T>& x)
        {
            typename std::map<unsigned int, std::vector<T> >::iterator it;
            if ((it = quad_point.find(n)) != quad_point.end()) {
                x = it->second;
                return;
            }

            unsigned int i, j;
            static const T pipm4 = pow(T(M_PI), -0.25);
            T r, r1, p1, p2, p3, dp3;

            x.resize(n);
            for (i = 0; i <= (n+1)/2-1; i++) {
                if (i == 0)
                    r = sqrt(T(2*n+1))-1.85575*pow(T(2*n+1), -T(1)/T(6));
                else {
                    if (i == 1)
                        r = r-1.14*pow(T(n), 0.426)/r;
                    else if (i == 2)
                        r = 1.86*r-0.86*x[0];
                    else if (i == 3)
                        r = 1.91*r-0.91*x[1];
                    else
                        r = 2*r-x[i-2];
                }

                do {
                    p2 = 0;
                    p3 = pipm4;
                    for (j = 0; j <= n-1; j++) {
                        p1 = p2;
                        p2 = p3;
                        p3 = p2*r*sqrt(T(2)/T(j+1))-p1*sqrt(T(j)/T(j+1));
                    }

                    dp3 = sqrt(T(2*j))*p2;
                    r1 = r;
                    r = r-p3/dp3;
                } while (fabs(r-r1) >= std::numeric_limits<T>::epsilon()*(1+fabs(r))*100);
                x[i] = r;
                x[n-1-i] = -x[i];
            }

            for (i = 0; i < n; i++) x[i] *= M_SQRT2;
            quad_point[n] = x;
        }

        /**
         * 做 Hermite 多项式变换。把多项式的系数值 @p c 变为插值点上的
         * 值 @p x。
         */
        static void HermiteTransform(const std::vector<T>& c, std::vector<T>& x)
        {
            typename std::map<unsigned int, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >::iterator it;
            unsigned int n = c.size();
            std::vector<T> result(n, 0);

            if ((it = trans_matrix.find(n)) == trans_matrix.end()) {
                _GenMatrix(n);
                it = trans_matrix.find(n);
            }

            for (unsigned int i = 0; i < n; i++)
                for (unsigned int j = 0; j < n; j++)
                    result[i] += c[j] * (it->second)(i, j);

            std::swap(result, x);
        }

        /**
         * 做 Hermite 多项式逆变换。把插值点上的值 @p x 变为多项式的系
         * 数值 @p c。
         */
        static void InvHermiteTransform(const std::vector<T>& x, std::vector<T>& c)
        {
            typename std::map<unsigned int, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> >::iterator it;
            unsigned int n = x.size();
            std::vector<T> result(n, 0);

            if ((it = inv_trans_matrix.find(n)) == inv_trans_matrix.end()) {
                _GenInvMatrix(n);
                it = inv_trans_matrix.find(n);
            }

            result.assign(n, 0);
            for (unsigned int i = 0; i < n; i++)
                for (unsigned int j = 0; j < n; j++)
                    result[i] += x[j] * (it->second)(i, j);

            std::swap(result, c);
        }

};

template <class T> std::map<unsigned int, std::vector<T> > Transform<T>::quad_point;
template <class T> std::map<unsigned int, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > Transform<T>::trans_matrix;
template <class T> std::map<unsigned int, Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> > Transform<T>::inv_trans_matrix;

template <class T> T MaxRootOfHermitePolynomial(unsigned int __deg) {
    return Transform<T>::MaxRootOfHermitePolynomial(__deg);
}

template <class T> void HermiteQuadPoint(unsigned int n, std::vector<T>& x) {
    Transform<T>::HermiteQuadPoint(n, x);
}

template <class T>
void HermiteTransform(const std::vector<T>& c, std::vector<T>& x) {
    Transform<T>::HermiteTransform(c, x);
}

template <class T>
void InvHermiteTransform(const std::vector<T>& x, std::vector<T>& c) {
    Transform<T>::InvHermiteTransform(x, c);
}


#endif // __Transform_h_
