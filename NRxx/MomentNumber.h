// vim: fdm=syntax:ts=4:sw=4:syntax=c.doxygen:tw=70:fo+=Mm

#ifndef __MomentNumber_h_
#define __MomentNumber_h_

#include "config.h"

NRXX_NAMESPACE_OPEN

/**
 * 这里先声明 @ref CoefTable 和 @ref Distribution 模板类，方便下面使用。
 */
template <int DIM, class _T> class CoefTable;
template <int DIM, class _T, int K> class Distribution;

/**
 * 用于计算某一个矩数的一个类。
 */
class MomentNumber {
    private:
        /**
         * 程序中所涉及的最高维数。
         */
        static const int __dim = MAX_VECTOR_DIMENSION;

		/**
         * 计算 1 到 @ref __dim 维速度空间中 1 至 @p __n 阶以下矩的个
         * 数。设返回值为 @p result，则 <tt>result[i-1+j*__dim]</tt>
         * 表示对于@p i 维的速度空间，@p j 阶以下矩的个数。
		 */
		static std::vector<unsigned int>
		_GetMomentNumber(unsigned int __n) {
			if (!__n) return std::vector<unsigned int>();
			std::vector<std::vector<unsigned int> > __n_mnt(__dim, std::vector<unsigned int>(__n));

			// 一维情形，矩的个数等于阶数加 1
			for (unsigned int __j = 0; __j < __n; __j++)
				__n_mnt[0][__j] = __j + 1;

			// 高维情形，矩的个数等于组合数 C(__i + __j, __i)
			for (int __i = 1; __i < __dim; __i++) {
				__n_mnt[__i][0] = 1;
				for (unsigned int __j = 1; __j < __n; __j++)
					__n_mnt[__i][__j] =
						__n_mnt[__i - 1][__j] + __n_mnt[__i][__j - 1];
			}

			// 整理成一维数组
			std::vector<unsigned int> __result(__dim * __n);
			for (int __i = 0; __i < __dim; __i++)
				for (unsigned int __j = 0; __j < __n; __j++)
					__result[__j * __dim + __i] = __n_mnt[__i][__j];
			return __result;
		}

		/**
		 * 获取 @p __i 维 @p __j 阶以下矩的个数。
		 */
		static int
		_GetMomentNumber(size_t __i, size_t __j) {
			static std::vector<unsigned int> __n = _GetMomentNumber(MAX_ORDER);
			return __n[__i - 1 + __j * __dim];
		}

        /**
         * 使得 @p CoefTable 类能调用该类中的函数。
         */
        template<int, class> friend class CoefTable;

        /**
         * 使得 @p Distribution 类能调用该类中的函数。
         */
        template<int, class, int> friend class Distribution;
};

NRXX_NAMESPACE_CLOSE

#endif // __MomentNumber_h_
