// vim: fdm=syntax:ts=4:sw=4:syntax=c.doxygen:tw=70:fo+=Mm

#ifndef __CoefTable_h_
#define __CoefTable_h_

#include <cstring>
#include <cstdarg>
#include <iterator>
#include <vector>
#include "config.h"
#include "MomentNumber.h"

NRXX_NAMESPACE_OPEN

/**
 * 用于表示一个 @p DIM 维表下标的类。
 */
template <int DIM>
class CoefTableIndices {
	private:
		/**
		 * 下标各分量的值。
		 */
		int __x[DIM];

	public:
		/**
		 * 默认构造函数。将下标值清零。
		 */
		CoefTableIndices() {
			memset(__x, 0, DIM * sizeof(int));
		}

		/**
		 * 复制构造函数。拷贝各下标。
		 */
		CoefTableIndices(const CoefTableIndices& __c) {
			memcpy(__x, __c.__x, DIM * sizeof(int));
		}

		/**
		 * 构造函数，设置各下标。
		 */
		CoefTableIndices(int __x1, ...) {
			va_list ap;
			va_start(ap, __x1);
			__x[0] = __x1;
			for (unsigned int __i = 1; __i < DIM; __i++) {
				__x1 = va_arg(ap, int);
				__x[__i] = __x1;
			}
		}

		/*@{*/
		/**
		 * 取出下标值。
		 */
		int& operator[](size_t __i) { return __x[__i]; }
		const int& operator[](size_t __i) const { return __x[__i]; }
		/*@}*/

		/**
		 * 赋值运算符。
		 */
		CoefTableIndices& operator=(const CoefTableIndices& __c) {
			memcpy(__x, __c.__x, DIM * sizeof(int));
			return *this;
		}
};

/**
 * 用于存储谱展开系数的模版类。参数 @p DIM 表示速度所在空间的维数，一
 * 般为 3。参数 @p _T 表示系数的类型，一般为 @p double 或 @p float。
 */
template <int DIM, class _T>
class CoefTable {
	public:
		/**
		 * 下标类型。
		 */
		typedef CoefTableIndices<DIM> Indices;

		/**
		 * 系数表中元素的类型。
		 */
		typedef _T value_type;

	private:
		/**
		 * 迭代器。模版参数 @p T1 用于区分 @p non-const 迭代器和 @p
		 * const 迭代器。对于我们的算法，一个单向迭代器已经够用了。
		 */
		template <class T1>
		class _IteratorBase : public std::iterator<std::forward_iterator_tag, T1> {
			public:
				/**
				 * 基类类型。
				 */
				typedef std::iterator<std::forward_iterator_tag, T1> _Base;

				/**
				 * 元素值类型。
				 */
				typedef typename std::iterator_traits<_Base>::value_type value_type;

				/**
				 * 元素指针类型。
				 */
				typedef typename std::iterator_traits<_Base>::pointer pointer;

				/**
				 * 元素引用类型。
				 */
				typedef typename std::iterator_traits<_Base>::reference reference;

			private:
				/**
				 * 元素下标，可通过 @ref GetIndices 取出。
				 */
				Indices __i;

				/**
				 * 指向当前元素的指针。
				 */
				pointer __p;

				/**
				 * 构造函数。直接对内部成员赋值，仅由 @ref CoefTable 的 @ref
				 * begin 和 @ref end 成员函数使用。
				 */
				_IteratorBase(const Indices& __i, pointer __p) :
					__i(__i), __p(__p) {}

			public:
				/**
				 * 默认构造函数。
				 */
				_IteratorBase() : __p(NULL) {}

				/**
				 * 前缀增量操作。
				 */
				_IteratorBase& operator++() {
					__p++;

					int __s = __i[DIM - 1], __n = DIM - 1;
					__i[DIM - 1] = 0;
					for (; __n > 0; __n--) {
						if (__i[__n - 1]) {
							__i[__n] = __s + 1;
							__i[__n - 1] --;
							return *this;
						}
					}
					__i[0] = __s + 1;
					return *this;
				}

				/**
				 * 后缀增量操作。
				 */
				_IteratorBase operator++(int) {
					_IteratorBase __bak(*this);
					operator++();
					return __bak;
				}

				/**
				 * 判断两个迭代器是否相等。为了节省时间，只对指针进行判断。
				 */
				bool operator==(const _IteratorBase& __it) const {
					return __p == __it.__p;
				}

				/**
				 * 判断两个迭代器是否不等。为了节省时间，只对指针进行判断。
				 */
				bool operator!=(const _IteratorBase& __it) const {
					return __p != __it.__p;
				}

				/**
				 * Dereference 操作。
				 */
				reference operator*() {
					return *__p;
				}
				
				/**
				 * 取指针操作。
				 */
				pointer operator->() {
					return __p;
				}

				/**
				 * 取出元素下标。
				 */
				const Indices& GetIndices() const {
					return __i;
				}

				/**
				 * 获取当前矩所在的阶数。
				 */
				int Order() const {
					int __order = 0;
					for (int __n = 0; __n < DIM; __n++)
						__order += __i[__n];
					return __order;
				}

				/**
				 * 允许 @ref CoefTable 类调用私有构造函数。
				 */
				friend class CoefTable<DIM, _T>;
		};

	public:
		/**
		 * @p non-const 迭代器。
		 */
		typedef _IteratorBase<_T> Iterator;

		/**
		 * @p const 迭代器。
		 */
		typedef _IteratorBase<const _T> ConstIterator;

	private:
		/**
		 * 获取 @p __i 维 @p __j 阶以下矩的个数。
		 */
		static int
		_GetMomentNumber(size_t __i, size_t __j) {
			return MomentNumber::_GetMomentNumber(__i, __j);
		}

		/**
		 * 展开的阶数。
		 */

		unsigned int __M;
		/**
		 * 展开的系数。
		 */
		_T* __x;

	private:
		/**
		 * 把用 @ref CoefTableIndices 表示的指标转化为数组 @ref __x
		 * 中的下指标。
		 */
		static int _IndexConvert(const Indices& __i) {
			int __d[DIM], __index = 0;
			__d[DIM - 1] = __i[DIM - 1];
			for (int __m = DIM - 1; __m > 0; __m--)
				__d[__m - 1] = __d[__m] + __i[__m - 1];

			for (int __m = 0; __m < DIM; __m++) {
				if (!__d[__m]) return __index;
				__index += _GetMomentNumber(DIM - __m, __d[__m] - 1);
			}
			return __index;
		}

	public:
		/**
		 * 默认构造函数。置指针为 NULL。
		 */
		CoefTable() : __M((unsigned int)(-1)), __x(NULL) {}

		/**
		 * 复制构造函数。
		 */
		CoefTable(const CoefTable& __c) {
			if (__c.__x == NULL) {
				__x = NULL; __M = (unsigned int)(-1);
			} else {
				int __n = _GetMomentNumber(DIM, __c.__M);
				__M = __c.__M;
				__x = new _T[__n];
				memcpy(__x, __c.__x, sizeof(_T) * __n);
			}
		}

		/**
		 * 构造函数。指定展开的阶数，并将分布函数清零。
		 */
		explicit CoefTable(unsigned int __m) {
			__M = __m;
			__x = new _T[_GetMomentNumber(DIM, __M)];
			memset(__x, 0, sizeof(_T) * _GetMomentNumber(DIM, __M));
		}

		/**
		 * 析构函数。释放空间。
		 */
		~CoefTable() {
			if (__x) delete [] __x;
		}

		/**
		 * 获取展开的阶数。
		 */
		unsigned int GetOrder() const {
			return __M;
		}

		/**
		 * 获取指标 @p __i 对应的元素。
		 */
		const _T& operator()(const Indices& __i) const {
			return __x[_IndexConvert(__i)];
		}

		/**
		 * 获取指标 @p __i 对应的元素。
		 */
		_T& operator()(const Indices& __i) {
			return __x[_IndexConvert(__i)];
		}

		/**
		 * 置整个系数表为一常数。为防止这一操作产生歧义，参数 @p __x
		 * 只允许为零。
		 */
		CoefTable& operator=(const _T& __zero) {
			if (__x)
				memset(__x, 0, sizeof(_T) * _GetMomentNumber(DIM, __M));
			return *this;
		}

		/**
		 * 赋值运算符。
		 */
		CoefTable& operator=(const CoefTable& __c) {
			if (__x != __c.__x) {
				if (__x) delete [] __x;
				if (__c.__x == NULL) {
					__x = NULL; __M = (unsigned int)(-1);
				} else {
					__M = __c.__M;
					int __n = _GetMomentNumber(DIM, __M);
					__x = new _T[__n];
					memcpy(__x, __c.__x, sizeof(_T) * __n);
				}
			}
			return *this;
		}

		/**
		 * 重新设定展开的阶数，清除原有的所有信息，并将分布函数清零。
		 */
		void Reinit(unsigned int __m) {
			__M = __m;
			if (__x) delete [] __x;
			__x = new _T[_GetMomentNumber(DIM, __M)];
			memset(__x, 0, sizeof(_T) * _GetMomentNumber(DIM, __M));
		}

		/**
		 * 重设阶数，保留原来的所有系数。
		 */
		void ResetOrder(unsigned int __m) {
			unsigned int __min = std::min(__M, __m);
			__M = __m;
			_T* __x1 = new _T[_GetMomentNumber(DIM, __M)];
			memset(__x1, 0, sizeof(_T) * _GetMomentNumber(DIM, __M));
			if (__x) {
				memcpy(__x1, __x, _GetMomentNumber(DIM, __min) * sizeof(_T));
				delete [] __x;
			}
			__x = __x1;
		}

		/**
		 * 返回指向首元素的迭代器。
		 */
		Iterator begin() {
			Indices __it;
			return Iterator(__it, __x);
		}

		/**
		 * 返回指向首元素的 @p const 迭代器。
		 */
		ConstIterator begin() const {
			Indices __it;
			return ConstIterator(__it, __x);
		}

		/**
		 * 返回指向末尾元素的迭代器。
		 */
		Iterator end() {
			Indices __it;
			__it[0] = __M + 1;
			return Iterator(__it, __x + _GetMomentNumber(DIM, __M));
		}

		/**
		 * 返回指向末尾元素的 @p const 迭代器。
		 */
		ConstIterator end() const {
			Indices __it;
			__it[0] = __M + 1;
			return ConstIterator(__it, __x + _GetMomentNumber(DIM, __M));
		}

		/**
		 * 返回指向 @p __m 阶矩首个系数的迭代器。
		 */
		Iterator begin(unsigned int __m) {
			Indices __it;
			if (!__m) return Iterator(__it, __x);
			__it[0] = __m;
			return Iterator(__it, __x + _GetMomentNumber(DIM, __m - 1));
		}

		/**
		 * 返回指向 @p __m 阶矩首个系数的 @p const 迭代器。
		 */
		ConstIterator begin(unsigned int __m) const {
			Indices __it;
			if (!__m) return ConstIterator(__it, __x);
			__it[0] = __m;
			return ConstIterator(__it, __x + _GetMomentNumber(DIM, __m - 1));
		}

		/**
		 * 返回系数表中的元素个数。
		 */
		unsigned int size() const {
			return (__x ? _GetMomentNumber(DIM, __M) : 0);
		}

		/**
		 * 取首系数。
		 */
		_T& front() { return __x[0]; }

		/**
		 * 取首系数的 @p const 版本。
		 */
		const _T& front() const { return __x[0]; }
};

NRXX_NAMESPACE_CLOSE

#endif // __CoefTable_h_
