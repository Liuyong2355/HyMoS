// vim: ts=4:sw=4:tw=70:fo+=Mm:filetype=c.doxygen

#ifndef __Vector_h_
#define __Vector_h_

#include <valarray>
#include <cmath>
#include <numeric>
#include <boost/static_assert.hpp>
#include <boost/preprocessor/arithmetic/add.hpp>
#include <boost/preprocessor/repetition.hpp>
#include "config.h"

NRXX_NAMESPACE_OPEN

template <size_t DIM, typename _T>
class Vector : public std::valarray<_T> {
	private:
		typedef std::valarray<_T> _Base;

	public:
		Vector() : _Base(DIM) {}
		template <class _C> Vector(const _C& __c) : _Base(__c) {}
		Vector(const _T& __x) : _Base(__x, DIM) {}
		Vector(const _T* __px) : _Base(__px, DIM) {}
		Vector(_T* __px) : _Base(__px, DIM) {}

		/**
		 * 构造函数。直接对每一维的速度赋值。这组函数由 boost 库的预处
		 * 理器实现，最多支持的向量维数为 @p MAX_VECTOR_DIMENSION。该
		 * 值默认为 3，若要增加支持的维数，则可在包含 headers.h 之前对
		 * 此宏重新进行定义。这个构造函数的用法如下所示：
		 * <pre>
		 *     Vector<3, float> v(3, 4, 5);
		 * </pre>
		 */
#define VECTOR_ASSIGN(z, n, unused) (*this)[n] = __x##n;
#define ADD2(n) BOOST_PP_ADD(n, 2)
#define VECTOR_CONSTRUCTOR(z, n, unused) \
		Vector( \
				BOOST_PP_ENUM_PARAMS(ADD2(n), const _T& __x) \
			  ) : _Base(DIM) { \
			BOOST_STATIC_ASSERT(DIM == n + 2); \
			BOOST_PP_REPEAT(ADD2(n), VECTOR_ASSIGN, ~) \
		}
		BOOST_PP_REPEAT(MAX_VECTOR_DIMENSION, VECTOR_CONSTRUCTOR, ~)
#undef VECTOR_ASSIGN
#undef VECTOR_CONSTRUCTOR

		_T Length() const {
			_T __ip = std::inner_product(&(*this)[0],
					&(*this)[DIM], &(*this)[0], _T(0));
			return sqrt(__ip);
		}

		template <class _C> Vector& operator=(_C __c) { _Base::operator=(__c); return *this; }
		template <class _C> Vector& operator+=(_C __c) { _Base::operator+=(__c); return *this; }
		template <class _C> Vector& operator-=(_C __c) { _Base::operator-=(__c); return *this; }
		template <class _C> Vector& operator*=(_C __c) { _Base::operator*=(__c); return *this; }
		template <class _C> Vector& operator/=(_C __c) { _Base::operator/=(__c); return *this; }
		template <class _C> Vector& operator^=(_C __c) { _Base::operator^=(__c); return *this; }
		template <class _C> Vector& operator&=(_C __c) { _Base::operator&=(__c); return *this; }
		template <class _C> Vector& operator|=(_C __c) { _Base::operator|=(__c); return *this; }
		template <class _C> Vector& operator%=(_C __c) { _Base::operator%=(__c); return *this; }
		template <class _C> Vector& operator<<=(_C __c) { _Base::operator<<=(__c); return *this; }
		template <class _C> Vector& operator>>=(_C __c) { _Base::operator>>=(__c); return *this; }

	private:
		void resize(size_t);
};

NRXX_NAMESPACE_CLOSE

#endif // __Vector_h_
