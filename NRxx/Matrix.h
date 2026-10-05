// vim: ts=4:sw=4:tw=70:fo+=Mm:filetype=cpp.doxygen

#ifndef __Matrix_h_
#define __Matrix_h_

#include <cstring>

/**
 * 矩阵类。
 */
template <typename _T> class Matrix
{
	private:
		/**
		 * 矩阵行数。
		 */
		unsigned int __n_row;

		/**
		 * 矩阵列数。
		 */
		unsigned int __n_column;

		/**
		 * 矩阵元素，按行存储。
		 */
		_T *__x;

	public:
		/**
		 * 矩阵元素类型。
		 */
		typedef _T value_type;

		/**
		 * 构造函数。构造一个空矩阵。
		 */
		Matrix() : __n_row(0), __n_column(0), __x(NULL) {}

		/**
		 * 构造函数。设置行数和列数。
		 */
		Matrix(unsigned int __n_row, unsigned int __n_column) :
			__n_row(__n_row), __n_column(__n_column)
		{
			if (__n_row != 0 && __n_column != 0) {
				__x = new _T[__n_row * __n_column];
				memset(__x, 0, sizeof(_T) * __n_row * __n_column);
			} else
				__x = NULL;
		}

		/**
		 * 构造函数。生成一个方阵。
		 */
		Matrix(unsigned int __n_row) : __n_row(__n_row), __n_column(__n_row)
		{
			if (__n_row != 0 && __n_column != 0) {
				__x = new _T[__n_row * __n_column];
				memset(__x, 0, sizeof(_T) * __n_row * __n_column);
			} else
				__x = NULL;
		}

		/**
		 * 析构函数。
		 */
		~Matrix() {
			delete[] __x;
		}

		/**
		 * 重设行列数。
		 */
		void Reinit(unsigned int __n_row, unsigned int __n_column)
		{
			if (__x != NULL) delete[] __x;
			if (__n_row != 0 && __n_column != 0) {
				this->__n_row = __n_row; this->__n_column = __n_column;
				__x = new double[__n_row * __n_column];
				memset(__x, 0, sizeof(_T) * __n_row * __n_column);
			} else
				__x = NULL;
		}

		/**
		 * 重设行列数。
		 */
		void Reinit(unsigned int __n_row)
		{
			if (__x != NULL) delete[] __x;
			if (__n_row != 0 && __n_column != 0) {
				this->__n_row = __n_row; this->__n_column = __n_row;
				__x = new double[__n_row * __n_column];
				memset(__x, 0, sizeof(_T) * __n_row * __n_column);
			} else
				__x = NULL;
		}

		/**
		 * 获取行数。
		 */
		unsigned int m() const { return __n_row; }

		/**
		 * 获取列数。
		 */
		unsigned int n() const { return __n_column; }

		/**
		 * 获取矩阵元素。
		 */
		_T& operator()(unsigned int __i, unsigned int __j)
		{
			return __x[__i * __n_column + __j];
		}

		/**
		 * 获取矩阵元素。
		 */
		const _T& operator()(unsigned int __i, unsigned int __j) const
		{
			return __x[__i * __n_column + __j];
		}
};

/**
 * 生成一个单位矩阵。
 */
template <class _T>
Matrix<_T> IdentityMatrix(unsigned int __size)
{
	Matrix<_T> __matrix(__size);
	for (unsigned int __i = 0; __i < __size; __i++)
		__matrix(__i, __i) = 1;
	return __matrix;
}

#endif // __Matrix_h_
