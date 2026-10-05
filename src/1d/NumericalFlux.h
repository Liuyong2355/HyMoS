/**
 * @file   NumericalFlux.h
 * @author HU Zhicheng <huzhicheng1986@gmail.com>
 * @date   Tue Jan  7 11:32:21 2014
 * 
 * @brief  数值通量计算 (包含某个矩数扰动后增量的快速计算)
 * 
 * 此函数与 NRxx/NumericalFlux.h 文件同名 (想不出其他名字), HLL Flux 在
 * NRxx/NumericalFlux.h 中已经给出, 这里主要是提供矩数扰动后增量的快速
 * 计算公式.
 * 
 */
#ifndef __SSP_NUMERICAL_FLUX_H__
#define __SSP_NUMERICAL_FLUX_H__

#include <NRxx/NumericalFlux.h> // 这个头文件就只放这里好了
#include <NRxx/Transform.h>

/**
 * @brief 数值通量类型枚举
 * 
 * 目前仅实现了 HLL 通量. 每种实现的数值通量都应该提供两个函数, 一个是
 * 解在边界的通量, 一个是扰动解后通量的增量计算函数.
 */
enum FluxType {HLL};

/** 
 * @brief HLL Flux 扰动分布函数时的增量
 *
 * @note 注意的几点:
 *
 * - 函数的实现公式, 参见 essay/fast-residual.tex 文档. 关于计算公式,
 *   有如下几点:
 * 
 *   - 这里仅针对不考虑速度和温度的矩系数扰动 ( 除 \f$ f_{e_j} \f$ 和
 *     \f$ f_{2e_j} \f$. 至于会改变速度温度的扰动变化会比较复杂, 目前考
 *     虑重新计算通量再减去旧通量的方式.
 * 
 *   - 本函数的实现没有考虑正则化部分扰动. 成对非守恒型正则化项通量形式
 *     并不了解, 不过这个问题不大. 事实上, 一方面正则化项仅依赖于最后二
 *     阶 (M-1, M 阶) 的矩系数, 因而正则化项扰动增量的影响有限. 另一方
 *     面, 目前 (local Newton) 我们还只考虑当前单元分布函数 (@p dis0)
 *     扰动对该单元对流项增量的影响, 若正则化项采用简单的中心差分离散,
 *     当前单元分布函数的扰动便不影响该单元对流项整体的增量.
 * 
 * - 与 NRxx 中 HLLFlux 函数一样, 这里要求传入的分布函数在标准展开下.
 */
template <typename DISTRIBUTION>
void DeltaHLLFlux( const DISTRIBUTION& dis0, 
                   const DISTRIBUTION& dis1, 
                   const std::vector<double>& normal,
                   const double df, 
                   const typename DISTRIBUTION::Indices& ind,
                   DISTRIBUTION& dflux_l2r ){
  typedef typename DISTRIBUTION::value_type value_t;

  static const int dim = DISTRIBUTION::dim;

  const unsigned int order = dis0.GetOrder();
  value_t C_M = Transform<value_t>::MaxRootOfHermitePolynomial( order+1 );

  const typename DISTRIBUTION::Velocity& c0 = dis0.Center();
  const typename DISTRIBUTION::Velocity& c1 = dis1.Center();
  const value_t& s0 = dis0.Scaling();
  const value_t& s1 = dis1.Scaling();    
  value_t u0 = std::inner_product(normal.begin(), normal.end(), &c0[0], (value_t)0);
  value_t u1 = std::inner_product(normal.begin(), normal.end(), &c1[0], (value_t)0);

  value_t lambda_l = std::min<value_t>( u0 - C_M * s0, u1 - C_M * s1 );
  value_t lambda_r = std::max<value_t>( u0 + C_M * s0, u1 + C_M * s1 );

  /**
   * 这前面一段算特征值倒也可以外头计算, 只是虽则有些效率, 程序更难看
   * 懂了.
   *
   */
  dflux_l2r.Reinit( dis0 );

  if( lambda_r <= 0 )   return;

  value_t lambda_2 = lambda_r / (lambda_r - lambda_l);
//    value_t lambda_2 = lambda_l / (lambda_r - lambda_l);
  value_t lambda_3 = lambda_l * lambda_r / ( lambda_r - lambda_l );

  typename DISTRIBUTION::Indices _ind = ind;
  unsigned int ind_order=0;
  for( int i=0; i<dim; ++i )   ind_order += _ind[i];

  for( unsigned int i=0; i<normal.size(); ++i ){
    if( lambda_l >= 0 ){
      dflux_l2r ( _ind ) += normal[i] * c0[i] * df;
      _ind[i] += 1;
      if( ind_order != order )   dflux_l2r( _ind ) += normal[i] * df * s0 * s0;
      _ind[i] -= 2;
      if( _ind[i] >= 0 )
        dflux_l2r( _ind ) += normal[i] * df * (_ind[i]+1);
    }
    else{
      /// 注意 lambda_3 这一项是不乘以法向的!
      dflux_l2r( _ind ) += ( normal[i] * lambda_2 * df * c0[i] - df * lambda_3 );   
      _ind[i] += 1;
      if( ind_order != order )
        dflux_l2r( _ind ) += normal[i] * df * s0 * s0 * lambda_2;
      _ind[i] -= 2;
      if( _ind[i] >= 0 )
        dflux_l2r( _ind ) += normal[i] * df * (_ind[i]+1) * lambda_2; 
    }
    _ind = ind;
  }
} // function DeltaHLLFlux()

#endif // __SSP_NUMERICAL_FLUX_H__

/**
 * end of file
 * 
 */
