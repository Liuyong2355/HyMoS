/**
 * @file   OperatorInterpolate.h
 * @author HU Zhicheng <huzhicheng1986@gmail.com>
 * @date   Tue Sep  9 21:02:28 2014
 * 
 * @brief  粗细网格间解的插值算子(限制与延拓)
 * 
 * 
 */

#ifndef __OPERATOR_INTERPOLATE_H__
#define __OPERATOR_INTERPOLATE_H__

#include <iostream>
#include <vector>

#include "Solution.h"

namespace Operator {
  typedef double value_t;

  /** 
   * @brief 两个细网格单元上的分布函数在粗网格单元固定参数系下的投影
   * 
   * @param dis_H 返回投影后粗网格上在固定参数系下的分布函数
   * @param dx_H 
   * @param dis_h_l 返回细网格上在此固定参数系下的分布函数
   * @param dx_h_l 
   * @param dis_h_r 同 @p dis_h_l
   * @param dx_h_r 以上参数同上一个函数
   * @param u_H 粗网格上的固定参数系
   * @param s_H 
   */
  template <typename DISTRIBUTION>
  void Restriction( DISTRIBUTION &dis_H, const value_t dx_H, 
                    DISTRIBUTION &dis_h_l, const value_t dx_h_l, 
                    DISTRIBUTION &dis_h_r, const value_t dx_h_r, 
                    const typename DISTRIBUTION::Velocity& u_H, const value_t s_H ) {
    dis_H.Reinit( u_H, s_H, std::min(dis_h_l.GetOrder(),dis_h_r.GetOrder()) );

    DISTRIBUTION tmp_l( u_H, s_H, dis_h_l.GetOrder() );
    DISTRIBUTION tmp_r( u_H, s_H, dis_h_r.GetOrder() );
    Project( dis_h_l, tmp_l );
    Project( dis_h_r, tmp_r );
    /**
     * 下面是计算粗网格分布的矩系数. 如果 u_H, s_H 是标准展开参数系, 那
     * 么其实只要算剩余的矩系数就可以了 (一阶矩不需算了). 但是将所有矩
     * 系数再按下面的方式计算, 效果还是一样, 因为粗网格的一阶矩系数按此
     * 式还是 0. 而且, 如果 u_H, s_H 不是标准展开参数系, 此时粗网格上未
     * 必是标准展开, 这个时候一阶矩系数等就需要用下面的式子计算了.
     * 
     */
    typename DISTRIBUTION::Iterator it_dis = dis_H.begin();
    typename DISTRIBUTION::Iterator the_end = dis_H.end();
    typename DISTRIBUTION::Iterator it_dis_l = tmp_l.begin();
    typename DISTRIBUTION::Iterator it_dis_r = tmp_r.begin();
    for(; it_dis != the_end; ++it_dis, ++it_dis_l, ++it_dis_r )
    {
      *it_dis = ((*it_dis_l) * dx_h_l + (*it_dis_r) * dx_h_r) / dx_H;
    }
    /// 返回粗细网格间的差异
    std::swap( tmp_l, dis_h_l );
    std::swap( tmp_r, dis_h_r );
  }

  /** 
   * @brief 两个细网格单元与所拼成的粗网格单元之间分布函数的投影
   *
   * 我们要求返回的粗网格单元上的分布函数 @p dis_H 是一个标准展开, 因而
   * 首先需要根据细网格分布函数计算粗网格上的展开参数系. 然后为了节省投
   * 影的次数, 此函数将修改两个细网格上分布函数的展开参数系为 @p dis_H
   * 的参数系, 即 @p dis_h_l 和 @p dis_h_r 返回在 @p dis_H 参数系下展开
   * 的细网格分布函数.
   * 
   * 由于现在每个单元上是分片常分布函数, 所以单元信息只要单元长度就够了.
   *
   * @param dis_H 返回粗网格上的分布函数, 标准展开
   * @param dx_H 粗网格单元长度
   * @param dis_h_l 
   * @param dx_h_l 
   * @param dis_h_r 
   * @param dx_h_r 
   */
  template <typename DISTRIBUTION>
  inline void Restriction( DISTRIBUTION& dis_H, const value_t dx_H, 
                           DISTRIBUTION& dis_h_l, const value_t dx_h_l,
                           DISTRIBUTION& dis_h_r, const value_t dx_h_r ){
    static unsigned int dim = DISTRIBUTION::dim;
    value_t pv_l[dim+2], pv_r[dim+2];

    dis_h_l.Moments( pv_l ); /**
                              * pv_l 里是 \rho, \rho * u, E = 0.5 *
                              * (\rho u^2 + D \rho \theta ) 向量.
                              */
    dis_h_r.Moments( pv_r );
    
    value_t rho_H = ( pv_l[0] * dx_h_l + pv_r[0] * dx_h_r ) / dx_H;
    typename DISTRIBUTION::Velocity u_H;
    value_t tmp = 0.;
    for( unsigned int i=0; i<dim; ++i )
    {
      u_H[i] = ( pv_l[i+1] * dx_h_l + pv_r[i+1] * dx_h_r ) / dx_H / rho_H;
      tmp += u_H[i] * u_H[i];
    }
    value_t theta_H = ( pv_l[dim+1] * dx_h_l + pv_r[dim+1] * dx_h_r ) / dx_H / rho_H - 0.5 * tmp;
    theta_H *= 2. / dim; 
    theta_H = sqrt( theta_H );
    
    Restriction( dis_H, dx_H, dis_h_l, dx_h_l, dis_h_r, dx_h_r, u_H, theta_H );
  }

  /** 
   * @brief 解的粗化
   *
   * 只处理粗网格由细网格两两合并而成的情形, 即细网格单元数是粗网格的两
   * 倍.
   *
   * 此函数有两个用处: 一是 multigrid 中的 Restriction, 此时为了减少投
   * 影次数, 将会更改细网格解的展开参数系; 二是作为普通的网格间的插值函
   * 数, 通常此时原网格的解就无用了, 所以更改细网格解也无妨.
   *
   * @param coarser_sol 粗网格上的解
   * @param coarser_mesh 
   * @param finer_sol 细网格上的解
   * @param finer_mesh 
   */
  template <typename MESH, typename SOLUTION>
  void Restriction( SOLUTION& coarser_sol, const MESH& coarser_mesh,
                    SOLUTION& finer_sol, const MESH& finer_mesh ) {
    for( unsigned int i=0; i<coarser_mesh.n_ele(); ++i )
    {
      Restriction( coarser_sol[i], coarser_mesh.Length(i), 
                   finer_sol[2*i], finer_mesh.Length(2*i), 
                   finer_sol[2*i+1], finer_mesh.Length(2*i+1) );
    }
  }


  /**
   * @brief 数值解的插值函数
   * 
   * 首先要求新旧网格的区域是相同的. 还有考虑两种特殊情形: 
   * - 新网格是旧网格加密一倍得到
   * - 新网格是旧网格稀疏一倍得到
   *
   * 对这两种特殊情形, 仅做一个简单的判断, 即单元数是否满足相关关系 (一
   * 倍或二分), 所以在使用时需要注意, 如果单元数满足相关关系, 而实际上
   * 新网格非由旧网格加密或合并而得, 那么返回的插值可能会是一个不正确的
   * 结果.
   *
   * 对于其他情形, 有一个非常土的插值方式.
   * 
   * @note 注意此函数可能改变旧网格上的解
   *
   * @param old_mesh 旧网格
   * @param old_sol 旧网格上的解
   * @param new_mesh 新网格
   * @param new_sol 新网格上的解, 要求空间已经分配, 与新网格长度一致.
   */
  template <typename MESH, typename SOLUTION, typename _BV>
  //template <typename MESH, typename SOLUTION>
  void Interp( const MESH& old_mesh, SOLUTION& old_sol,
               const MESH& new_mesh, SOLUTION& new_sol ) {
    if( 2*old_sol.size() == new_sol.size() ) 
    {
      /**
       * 此时就认为新网格是旧网格的加密
       *
       * 此时, 旧网格的 i 单元对应于新网格上的 2i 和 2i+1 这两个单元. 
       */
#if 0
      /**
       * 先考虑恒等插值算子, 即细网格两个单元都直接赋值为旧网格相关单元
       * 上的解
       * 
       */
      for( unsigned int i=0; i<old_sol.size(); ++i ) 
      {
        new_sol[2*i] = old_sol[i];
        new_sol[2*i+1] = old_sol[i];
      }
#else
    typedef typename SOLUTION::dis_t dis_t;
      /*
       * 做一个线性重构的插值可能效果更好, 可以使得初始残量较小, 不过似
       * 乎对于mg时无所谓, 一步 mg 就可以迅速将残量减小下来.
       */
      dis_t lbv, rbv;
      _BV::BoundaryValue( lbv, old_sol[0], false );
      _BV::BoundaryValue( rbv, old_sol[old_sol.size()-1], true );

      for( size_t i=0; i<old_sol.size(); ++i )
      {
        /// 单元中心与左单元中心的间距
        const value_t dx_l = (i!=0) ? 0.5*(old_mesh.Length(i-1)+old_mesh.Length(i)) : old_mesh.Length(i); 
        /// 单元中心与右单元中心的间距
        const value_t dx_r = (i!=old_sol.size()-1) ? 0.5*(old_mesh.Length(i)+old_mesh.Length(i+1)) : old_mesh.Length( i );

        const dis_t& dis_l = (i!=0) ? old_sol[i-1] : lbv;
        const dis_t& dis_r = (i!=old_sol.size()-1) ? old_sol[i+1] : rbv;
        const dis_t& dis = old_sol[i];
        dis_t pro_dis_l( dis.Center(), dis.Scaling(), dis.GetOrder() );
        Project( dis_l, pro_dis_l );
        dis_t pro_dis_r( dis.Center(), dis.Scaling(), dis.GetOrder() );
        Project( dis_r, pro_dis_r );

        dis_t recon_dis_l( dis.Center(), dis.Scaling(), dis.GetOrder() );
        dis_t recon_dis_r( dis.Center(), dis.Scaling(), dis.GetOrder() );

        typename dis_t::ConstIterator it_dis = dis.begin();
        const typename dis_t::ConstIterator it_end = dis.end();
        typename dis_t::Iterator it_pro_l = pro_dis_l.begin();
        typename dis_t::Iterator it_pro_r = pro_dis_r.begin();
        typename dis_t::Iterator it_recon_l = recon_dis_l.begin();
        typename dis_t::Iterator it_recon_r = recon_dis_r.begin();
        for( ; it_dis != it_end; ++it_dis, ++it_pro_l, ++it_pro_r, ++it_recon_l, ++it_recon_r )
        {
          value_t slope_l = (*it_dis - *it_pro_l) / dx_l;
          value_t slope_r = (*it_pro_r - *it_dis) / dx_r;
        
          slope_l = minmod( slope_l, slope_r );

          *it_recon_l = *it_dis - slope_l * new_mesh.Length( 2*i ) * 0.5;
          *it_recon_r = *it_dis + slope_l * new_mesh.Length( 2*i+1 ) * 0.5;
        }

        dis_t& new_dis_l = new_sol[ 2*i ];
        dis_t& new_dis_r = new_sol[ 2*i+1 ];
        new_dis_l.Reinit( dis.Center(), dis.Scaling(), dis.GetOrder() );
        new_dis_r.Reinit( dis.Center(), dis.Scaling(), dis.GetOrder() );
        new_dis_l.ProjectToStdSpace( recon_dis_l );
        new_dis_r.ProjectToStdSpace( recon_dis_r );

        /// 标准基下, f_1 = 0.0; 下面这句不过是去掉那些机器精度下的扰动,
        /// 是可以不要的.
        // 这只是对一维的语句
        // *new_dis_l.begin(1) = 0.0;  *new_dis_l.begin(2)=0.0;
        // *new_dis_r.begin(1) = 0.0;  *new_dis_r.begin(2)=0.0;
      }
#endif
      return;
    }
  
    if( old_sol.size() == 2*new_sol.size() ) 
    {
      /**
       * 但凡此种情形, 就认为旧网格是新网格的加密
       * 
       * 此时, 新网格的 i 单元对应于旧网格的 2i 和 2i+1 这两个单元
       * 
       */
      Restriction( new_sol, new_mesh, old_sol, old_mesh );
      return;
    }

    /**
     * 其他情形. 下面是最简单的策略, 就是新网格中点落在旧网格中哪个单元,
     * 就用哪个单元的解.
     * 
     */
    for( size_t i=0, j=0; i<new_sol.size(); ++i )
    {
      double x = new_mesh.Center(i);
      /// 寻找 x 所在旧网格的单元, x_r 是第 j 个单元的有边界.
      double x_r = old_mesh[j+1];
      while( x > x_r )
      {
        ++j;
        x_r = old_mesh[j+1];
      }
      new_sol[i] = old_sol[ j ];
    }
  }

  /** 
   * @brief 分布函数的外推
   * 
   * @param new_dis 
   * @param dis_H 
   * @param dis_h 
   * @param _ratio 
   */
  template <typename DISTRIBUTION>
  void Extrap( DISTRIBUTION& new_dis, const DISTRIBUTION& dis_H, 
               const DISTRIBUTION& dis_h, const double& _ratio ) {
    DISTRIBUTION tmp_H( dis_h.Center(), dis_h.Scaling(), dis_h.GetOrder() );
    DISTRIBUTION new_tmp( dis_h.Center(), dis_h.Scaling(), dis_h.GetOrder() );

    Project( dis_H, tmp_H );
    typename DISTRIBUTION::Iterator the_dis = new_tmp.begin();
    typename DISTRIBUTION::Iterator the_end = new_tmp.end();
    typename DISTRIBUTION::ConstIterator the_dis_H = tmp_H.begin();
    typename DISTRIBUTION::ConstIterator the_dis_h = dis_h.begin();
    for( ; the_dis != the_end; ++the_dis, ++the_dis_H, ++the_dis_h )
    {
      * the_dis = Extrap( *the_dis_H, *the_dis_h, _ratio );
    }
    new_dis.Reinit( dis_h );
    new_dis.ProjectToStdSpace( new_tmp );
  }

  /** 
   * @brief 通过外推法提高解的精度
   * 
   * 完全类似于数值积分中的外推方法和数值微分中的外推方法. 假设在某点有
   * \f$ f - f_h = C h^m + O(h^{m+1}), \f$ 相应的就有 \f$ f - f_{ph} =
   * Cp^m h^m + O(h^{m+1}), \f$, 组合两式可得: \f[ f - \frac{f_{ph} -
   * p^m f_h }{ 1-p^m } = O( h^{m+1} ). \f] 也就是说, 通过两个低精度的
   * 解的组合可以得到高一阶精度的解. 这就是外推的基本递推公式. 通常
   * \f$ p \f$ 取为 0.5.
   * 
   * 再对整个区域的解做外推的时候, 虽然对是否是均匀网格以及 \f$ p \f$的
   * 比例, 没有必要的限制. 但由于其它情形计算同一个位置上两个网格尺度上
   * 的解相对麻烦, 因而使用此函数时最好使用 \f$ p=0.5 \f$ 的两个网格(这
   * 样的话即使非均匀网格, 只要细网格是粗网格的加密, 就容易实现 )
   *
   * @param new_sol 外推的解, 与 @p sol_h 尺度相同
   * @param sol_H 网格尺度为 H 的解
   * @param sol_h 网格尺度为 h 的解
   * @param accuracy_order 即公式中的阶数 \f$ m \f$.
   * @param ratio_of_h 即公式中的 \f$ p \f$, 等于 h/H
   */
  template <typename DISTRIBUTION>
  void Extrap( Solution<DISTRIBUTION>& new_sol, 
               const Solution<DISTRIBUTION>& sol_H, 
               const Solution<DISTRIBUTION>& sol_h, 
               const unsigned int& accuracy_order, 
               const double& ratio_of_h=0.5 ){
    const unsigned int n_h = sol_h.size();
    if( n_h != sol_H.size() * 2 )
      std::cerr << "Extrapolation might have problem since ratio of h is not  0.5..." << std::endl;

    double p_m = pow( ratio_of_h, accuracy_order );
    new_sol.resize( n_h );
    for(unsigned int i=0; i<n_h; ++i )
    {
      Extrap( new_sol[i], sol_H[ i/2 ], sol_h[i], p_m );
    }
  }

} // namespace Operator
 

#endif // __OPERATOR_INTERPOLATE_H__

/**
 * end of file
 * 
 */
