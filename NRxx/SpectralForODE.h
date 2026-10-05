#ifndef __SPECTRAL_H__
#define __SPECTRAL_H__

#include <vector>
#include <map>
#include <fftw3.h>
#include "config.h"

#ifdef MULTITHREAD_NRXX
#include <boost/thread.hpp>
#endif // MULTITHREAD_NRXX

#define MAX_SM_NODE 64

inline double even(int m)
{
	if(m%2) return -1.;
	return 1.;
}

/** 
 * @brief 用于进行 DCT 变换，只需进行一次初始化即可，这样可以节省很多时间。
 */
class DCT
{
	private:
		std::vector<double> in, out;
		fftw_plan p;
		bool initialized;
	public:
		DCT() : initialized(false) {}

		~DCT() {
			if (initialized) fftw_destroy_plan(p);
		}

		void Reinit(unsigned int n) {
			if (n + 1 == in.size()) return;
			if (initialized) fftw_destroy_plan(p);
			in.resize(n+1); out.resize(n+1);
			p = fftw_plan_r2r_1d(n+1, &in[0], &out[0], FFTW_REDFT00, FFTW_MEASURE);
			initialized = true;
		}

		template <typename T1, typename T2>
			void Transform(const std::vector<T1>& _in, std::vector<T2>& _out) {
				std::vector<double> local_in(in.size());
				std::vector<double> local_out(out.size());

				typename std::vector<T1>::const_iterator it_in = _in.begin();
				std::vector<double>::iterator it = local_in.begin();
				std::vector<double>::iterator it_end = local_in.end();
				for (; it != it_end; ++it, ++it_in) *it = *it_in;

				fftw_execute_r2r(p, &local_in[0], &local_out[0]);

				_out.resize(local_out.size());
				typename std::vector<T2>::iterator it_out = _out.begin();
				it = local_out.begin(), it_end = local_out.end();
				for (; it != it_end; ++it, ++it_out) *it_out = *it;
			}
};

#ifdef MULTITHREAD_NRXX
class DCTMap : private std::map<unsigned int, DCT> {
	private:
		boost::mutex mutex;
	public:
		template <typename T1, typename T2>
			void Transform(unsigned int n, const std::vector<T1>& _in, std::vector<T2>& _out) {
				iterator it;
				boost::unique_lock<boost::mutex> lock(mutex);
				if ((it = find(n)) != end()) {
					lock.unlock();
					it->second.Transform(_in, _out);
				} else {
					DCT& dct = (*this)[n];
					dct.Reinit(n);
					lock.unlock();
					dct.Transform(_in, _out);
				}
			}
};
#else
class DCTMap : private std::map<unsigned int, DCT> {
	public:
		template <typename T1, typename T2>
			void Transform(unsigned int n, const std::vector<T1>& _in, std::vector<T2>& _out) {
				iterator it;
				if ((it = find(n)) != end()) {
					it->second.Transform(_in, _out);
				} else {
					DCT& dct = (*this)[n];
					dct.Reinit(n);
					dct.Transform(_in, _out);
				}
			}
};
#endif // MULTITHREAD_NRXX

inline std::vector<std::vector<double> > COS()
{
	std::vector<std::vector<double> > vec(MAX_SM_NODE);
	for(int __i = 1; __i <= MAX_SM_NODE; ++__i)
	{
		vec[__i-1].resize(__i + 1);
		for(int __j = 0; __j <= __i; ++__j)
		{
			vec[__i-1][__j] = cos(__j * M_PI / __i);
		}
	}
	return vec;
}

/** 
 * @brief 谱方法实现的主要部分, 使用谱方法求解 ODE 方程
 * @p y_0 为初值，@p y 为右端项, 最后解仍存在 @p y 中。
 */
template <class _T>
void SpectralForODE(_T y_0, std::vector<_T>& y)
{
	static DCTMap dct_map;
	static std::vector<std::vector<double> > _cos( COS() );
	int N = y.size() - 1;

	std::vector<_T> C(N+1);
	std::vector<_T> tmp(N+1);
	for(int i = 0; i <= N; ++i)
	{
		tmp[i] = y[i] * 0.5;
	}

	dct_map.Transform(N, tmp, C);

	C[0] /= 2. * N;
	for(int i = 1; i < N; ++i)
	{
		C[i] /= N;
	}
	C[N] /= 2. * N;

	tmp[0] = 0;
	tmp[1] = (2*C[0] - C[2]) * 0.5;
	for(int i = 2; i < N; ++i)
	{
		tmp[i] = ( C[i-1] - C[i+1] ) / i * 0.5;
	}
	tmp[N] = C[N-1] / N;

	dct_map.Transform(N, tmp, y);

	_T connum = y_0;
	_T tmpnum = C[N] / (N+1);
	connum += C[0] * 2. - C[2];
	for(int i = 2; i < N; ++i)
	{
		connum += (C[i-1] - C[i+1]) / i * even(i-1);
	}
	connum += C[N-1] / N * even(N-1);
	connum += tmpnum * even(N);

	for(int i = 0; i <= N; ++i)
	{
		y[i] += connum + tmpnum * even(i) * _cos[N-1][i];
	}
}

#endif //__SPECTRAL_H__
