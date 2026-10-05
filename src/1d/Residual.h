/**
 * @file   Residual.h
 * @author HU Zhicheng <huzhicheng1986@gmail.com>
 * @date   Mon Dec  9 00:02:54 2013
 * 
 * @brief  矩方程组残量计算框架
 * 
 * 
 */

#ifndef __RESIDUAL_H__
#define __RESIDUAL_H__

#include <omp.h>

#include "config.h"
#include "NumericalFlux.h"
#include "OpConvection.h"
#include "scat.h"
#include "tools/tools.h"

/**
 * @class Residual
 * @brief 矩方程组残量计算框架
 *
 * 类中不在定义数据成员, 仅仅是一些函数集合.
 *
 * 残量的计算需要一个具体问题原型 PROBLEM. 这个 PROBLEM 要提供网格类型,
 * 数值通量计算, 外力计算, 碰撞项处理, 以及边界条件等. 各部分的计算都应
 * 包含某个矩数扰动后增量的快速计算函数. 由于外力和碰撞项都由单独的类处
 * 理, PROBLEM 类的那层封装似无特别需要. 现在就直接把各部分作为模板参数
 * 传给 Residual 类. 也就是 Residual 类的模板参数有:
 *
 * - MESH: 单层网格类型
 * - SOLUTION: Boltzmann 方程解类型 (现在是 DISTRIBUTION 的一个 vector)
 * - FluxType: 通量枚举类型
 * - FORCE: 外力类, 这里需要 Boltzmann 方程外力项处理部分函数
 * - ColModel: 碰撞模型枚举类
 * - _BV: 提供初边值条件的具体问题类, 还包括特殊问题的参数等信息
 *
 * 由于函数模板不能部分特化, 暂时又不想对通量部分再封装一层类, 所以当前
 * FluxType 这个模板类可以说是无用的. 现在我们要求在 _BV 中提供
 * NumericalFlux() 和 GetDeltaFlux() 这两个函数, 当前可参考
 * cases/1d/CouetteCase.h
 * 
 */
template <typename MESH, typename SOLUTION, FluxType FT, typename FORCE, ColModel CM, typename _BV>
class Residual {
public:
  typedef MESH mesh_t;  // 这里再加一个 typename 编译竟然错了
  typedef SOLUTION sol_t;
  typedef FORCE force_t;
  typedef _BV PROBLEM;

  typedef typename sol_t::dis_t dis_t;
  typedef typename dis_t::value_type value_t;
  typedef typename dis_t::Indices index_t;
  typedef typename dis_t::Velocity velocity_t;

  enum { dim = dis_t::dim  };

public:
  static void GetResidualL2Norm(value_t& l2_norm, const sol_t& res, 
                                const mesh_t& mesh ) {
    l2_norm = 0.0;
    for( unsigned int i=0; i<res.size(); ++i )
    {
      double local_norm = L2NormDis( res[i] );
      l2_norm += local_norm * local_norm * mesh.Length( i );
    }
    // l2_norm = sqrt( l2_norm );
    /**
     * 除非区域做了无量纲化, 使得长度为 1, 否则这样得到的 l2_norm 将比
     * 局部 l2_norm 大 根号 L 倍, 其中 L 是区域长度. 于是在 SLGSNewton
     * 中用局部 l2_norm 判断已达局部稳态, 不再计算, 而全局来看, 仍未达
     * 稳态, 于是就有了一个困扰, SLGSNewton 停在某处不收敛.
     * 
     * 基于此, 对这全局 l2_norm 做个归一化, 或者说是求平均的 l2_norm.
     * 
     */    
    l2_norm = sqrt( l2_norm / mesh.DomainSize() );
  }
  static void GetResidualL2Norm( value_t& l2_norm, const sol_t& sol,
                                  const sol_t& rhs, const mesh_t& mesh, 
                                  const force_t& force_handle ) {
    sol_t res( sol.size() );
    GetResidual( res, sol, rhs, mesh, force_handle );
    GetResidualL2Norm( l2_norm, res, mesh );
  }

  /** 
   * 当 M > 10 (特别是大于 20) 时, 由于机器精度等因素, L2 norm 未必能准
   * 确计算, 这极有可能是此前观察到 L2 残量不下降而 max 残量仍然下降的
   * 原因. 考虑到实际我们观测的都是那些低阶矩信息, 所以就用求和式中前面
   * 几个之和代替 L2 norm 也是可取的. 按 PL2NormDis 中的定义, 这里求和
   * 式中仅计算指标小于 \alpha 小于等于 3 的那些项
   *
   */
  static void GetResidualPL2Norm(value_t& pl2_norm, const sol_t& res, 
                                 const mesh_t& mesh ) {
    pl2_norm = 0.0;
    for( unsigned int i=0; i<res.size(); ++i )
    {
      double local_norm = PL2NormDis( res[i] );
      pl2_norm += local_norm * local_norm * mesh.Length( i );
    }
    /**
     * 对这全局 l2_norm 做个归一化, 或者说是求平均的 l2_norm.
     * 
     */    
    pl2_norm = sqrt( pl2_norm / mesh.DomainSize() );
  }
  static void GetResidualPL2Norm( value_t& pl2_norm, const sol_t& sol,
                                  const sol_t& rhs, const mesh_t& mesh, 
                                  const force_t& force_handle ) {
    sol_t res( sol.size() );
    GetResidual( res, sol, rhs, mesh, force_handle );
    GetResidualPL2Norm( pl2_norm, res, mesh );
  }

  /** 
   * @brief 方程组残量残量最大值
   *
   * 为判断算法收敛与否, 需要检查所有单元上的残量是否足够小. 所有单元上
   * 残量最大值就存储在 @p max_res 中.
   * 
   * 这是一个通用函数, 甚至可以写到通用的 tools.h 中去
   * 
   * @param max_res 返回最大残量
   * @param res 残量
   */
  static void GetResidualMaxNorm( value_t& max_res, 
                                  const sol_t& res ) {
    max_res = 0.0;
    for( size_t i=0; i<res.size(); ++i )
    {
      MaxMomentsValue( res[i], max_res );
    }
  }
  /**
   * @brief 计算残量最大值用以判断是否收敛的函数接口
   *
   * 如此分开程序, 要多线程的话就要超多的 Helper 函数了 
   */
  static void GetResidualMaxNorm( value_t& max_res, const sol_t& sol,
                                  const sol_t& rhs, const mesh_t& mesh, 
                                  const force_t& force_handle ) {
    sol_t res( sol.size() );
    GetResidual( res, sol, rhs, mesh, force_handle );
    GetResidualMaxNorm( max_res, res );
  }
  /**
   * @brief 求解时间发展方程的函数接口
   * 
   */
  static void GetResidual( value_t& max_res, sol_t& res, 
                           const sol_t& sol, const sol_t& rhs, 
                           const mesh_t& mesh, 
                           const force_t& force_handle ) {
    GetResidual( res, sol, rhs, mesh, force_handle );
    //   GetResidualMaxNorm( max_res, res );
    //   GetResidualL2Norm( max_res, res, mesh );
    GetResidualPL2Norm( max_res, res, mesh );
  }
  /**
   * @brief 右端项为0时方程残量计算接口
   *
   * 即此函数计算 -R( @p sol ) 存在 @p res 并返回
   */
  static void GetResidual( sol_t& res, const sol_t& sol, 
                           const mesh_t& mesh, 
                           const force_t& force_handle ) {
    #pragma omp parallel for
    for( int i=0; i<(int)sol.size(); ++i )
    {
      res[i].Reinit( sol[i] );
      GetResidual( res[i], sol, mesh, force_handle, i );
    }
  }
  /** 
   * @brief 方程组残量计算
   *
   * 给定解与右端项, 返回方程残量 @p rhs - R( @p sol )
   * 
   * @param res 返回方程残量
   * @param sol 解
   * @param rhs 方程右端项
   * @param mesh 网格
   * @param force_handle 外力项处理句柄, 提供 @p mesh 上的外力等信息
   */
  static void GetResidual( sol_t& res, 
                           const sol_t& sol, const sol_t& rhs, 
                           const mesh_t& mesh,
                           const force_t& force_handle ) {
    #pragma omp parallel for
    for( int i=0; i<(int)sol.size(); ++i )
    {
      res[i].Reinit( sol[i] );
      Project( rhs[i], res[i] );
      GetResidual( res[i], sol, mesh, force_handle, i ); 
    }
  }
  
  /** 
   * @brief 矩方程组某个单元残量计算接口
   * 
   * 按照统一的方式, 各网格上要解 R(f) = r. 最细网格上右端 r 为 0, 粗网
   * 格 未必是 0. 而左边算子 R 已经包含了对流, 外力, 碰撞算子. 本函数将
   * 传入右端 r, 再程序体内计算 R(f), 然后将 r-R(f) 返回 (最后方程残量).
   *
   * 特别注意这个函数与之前实现的差别! 之前右端为 0, 故而计算 R(f), 现
   * 在实际是 -R(f) 了.
   *
   * R(f) 计算公式为 \f[ R(f) = \frac{1}{\Delta x} \left( \hat{F}(f_l,
   * f) - \hat{F}(f,f_r) \right) - G(f) - S(f). \f] 于是计算就分为了三
   * 部分:
   *
   * - 对流部分 \f$ \frac{\hat{F}(f_l, f) - \hat{F}(f,f_r)}{\Delta x}\f$
   * - 外力部分 \f$ G(f) \f$
   * - 散射部分 \f$ S(f) \f$
   *
   * @param res_dis 传入右端项, 返回方程的残量. 参数系与第 k 单元解相同 
   * @param sol 解
   * @param mesh 网格
   * @param force_handle 外力项处理句柄, 提供 @p mesh 上的外力等信息
   * @param i_ele 当前处理单元指标
   */
  static void GetResidual( dis_t& res_dis, const sol_t& sol, 
                           const mesh_t& mesh, 
                           const force_t& force_handle,
                           const unsigned int i_ele ) {
    /// 对流项
    Convection( res_dis, sol[i_ele], sol, mesh, i_ele );

    /// 外力
    force_handle.GetAcceleration( res_dis, sol[i_ele], i_ele );

    /// 碰撞
    Collision<CM, dis_t, force_t, _BV>( res_dis, sol[i_ele], force_handle, i_ele );
//    PROBLEM::Collision( res_dis, sol[i_ele], i_ele );
  }

  static void GetResidual( dis_t& res_dis, 
                           const dis_t& new_dis,
                           const sol_t& sol, 
                           const mesh_t& mesh, 
                           const force_t& force_handle,
                           const unsigned int i_ele ) {
    /// 对流项
    Convection( res_dis, new_dis, sol, mesh, i_ele );

    /// 外力
    force_handle.GetAcceleration( res_dis, new_dis, i_ele );

    /// 碰撞
    Collision<CM, dis_t, force_t, _BV>( res_dis, new_dis, force_handle, i_ele );
//    PROBLEM::Collision( res_dis, sol[i_ele], i_ele );
  }

  /** 
   * @brief 计算一个单元扰动分布函数后的该单元上残量增量
   *
   * 返回变量 @p delta_rhs 每个分量存储的是: \f$ R(f1) - R(f) \f$, 其中
   * \f$ f1 \f$ 是扰动后的分布函数.
   *
   * 最简单的方式是: 扰动 @p dis 的某个分量, 就计算新的残量. 修改
   * GetResidual() 函数接口为原来传左中右分布函数的形式, 此程序就会简单
   * 很多 (本程序差不多要重新 copy 一份 GetResidual() 函数的执行过程).
   * 但是, 由于每次都仅仅扰动 @p dis 的一个分量, 在使用 GetResidual()
   * 函数计算时自然重复计算了许多信息, 效率也就下来了. 因而最终将要弃用
   * 此方式.
   *
   * 更快速的计算是扰动 @p dis 以后, 直接计算相应的增量. 只是这样
   * 需要具体的增量表达式, 这个推导比较复杂, 而且对不同的碰撞项都需要做
   * 推导! 具体实现见 GetFastDeltaResidual() 函数
   *
   * 对于扰动变量, 我们选择 \f$ w = [f_0, f_1,f_2,\ldots, f_M]^T
   * \f$. 另一种选择 \f$ w = [f_0, u, s, f_3, \ldots, f_M]^T \f$ 如今已
   * 经弃用, 因为两者相差的只有三四次(看具体实现, 还要乘以2)投影.
   *
   * @param dw 返回值, 扰动向量, 第 i 个值对应第 i 次对第 i 个变量的扰动值
   * @param delta_rhs 返回值, 包含 M+1 个分布函数, 对应于 M+1 个扰动的残
   *                  量.
   *
   * @param sol 解
   * @param mesh 网格
   * @param force_handle 外力项处理句柄, 提供 @p mesh 上的外力等信息
   * @param i_ele 当前处理单元指标, 要扰动的分布函数指标
   * @param perturbation_ratio 扰动比例
   */
  static void GetDeltaResidual( std::vector<value_t>& dw, 
                                std::vector<dis_t>& delta_rhs, 
                                const sol_t& sol, const mesh_t& mesh,
                                const force_t& force_handle,
                                const unsigned int i_ele, 
                                const value_t perturbation_ratio ) {
    const dis_t &dis = sol[i_ele];
    const unsigned int n_moments = dis.size();
    dw.resize( n_moments );
    delta_rhs.resize( n_moments );

    /// 离散方程 R(f) = r, 右端项 r 与 f 无关, 故增量与 r 无关.
    dis_t res_dis( dis.Center(), dis.Scaling(), dis.GetOrder() );
    /// 计算未扰动前的残量 -R(f) (因为右端项设为0).
    GetResidual( res_dis, sol, mesh, force_handle, i_ele );
    
    dis_t delta_dis = dis; // dis 是 const, 所以要先复制一份

    /**
     * 扰动分布函数以后的残量计算. 虽然 GetResidual() 现在的传入参数的
     * 形式利于空间高维问题一致的接口, 但此处却不适合调用该函数了. 因此
     * 函数效率毕竟不高, 以后将用更高效率的函数替换, 所以这里就照搬
     * GetResidual() 函数体中的内容, 仅做简单的修改.
     * 
     */
    /// 矩系数扰动
    typename dis_t::Iterator the_dis = delta_dis.begin();
    for( size_t i=0; i<n_moments; ++i, ++the_dis )
    {
      delta_rhs[i] = res_dis; /// -R(f): f 是扰动前的分布函数

      const value_t dis_val_before_ran = *the_dis;
      // 用这个bool变量来判断扰动是否仍在标准展开下, 若否, 使用常规的残
      // 量计算方式
      bool is_standard_expansion = true;
      GetDeltaValue( dw[i], is_standard_expansion, i, 
                     dis_val_before_ran, perturbation_ratio, dis );

      *the_dis += dw[i];    

      dis_t *p_delta_dis = &delta_dis;
      dis_t tmp( dis.Center(), dis.Scaling(), dis.GetOrder() );
      if( !is_standard_expansion )
      {// 可能改变了标准展开参数系
        tmp.ProjectToStdSpace( delta_dis );
        p_delta_dis = &tmp;
      }

      dis_t new_res_dis( p_delta_dis->Center(), p_delta_dis->Scaling(), p_delta_dis->GetOrder() );

      Convection( new_res_dis, *p_delta_dis, sol, mesh, i_ele );
      force_handle.GetAcceleration( new_res_dis, *p_delta_dis, i_ele );
      Collision<CM, dis_t, force_t, _BV>( new_res_dis, *p_delta_dis, force_handle, i_ele );
//      PROBLEM::GetAcceleration( new_res_dis, *p_delta_dis, i_ele );
//      PROBLEM::Collision( new_res_dis, *p_delta_dis, i_ele );
      /// 至此, new_res_dis = - R(f1), 其中 f1 是扰动后的分布函数
      if( !is_standard_expansion )
      {
        tmp.Reinit( res_dis );
        Project( new_res_dis, tmp );
        delta_rhs[i] -= tmp;
      }
      else
        delta_rhs[i] -= new_res_dis;

      *the_dis = dis_val_before_ran;
    }
  }

  /** 
   * @brief 快速计算一个单元扰动分布函数后的该单元上残量增量
   * 
   * @see GetDeltaResidual() 函数
   *
   * @param dw 返回值, 扰动向量, 第 i 个值对应第 i 次对第 i 个变量的扰动值
   * @param delta_rhs 返回值, 包含 M+1 个分布函数, 对应于 M+1 个扰动的残
   *                  量.
   *
   * @param sol 解
   * @param mesh 网格
   * @param force_handle 外力项处理句柄, 提供 @p mesh 上的外力等信息
   * @param i_ele 当前处理单元指标, 要扰动的分布函数指标
   * @param perturbation_ratio 扰动比例
   */
  static void GetFastDeltaResidual( std::vector<value_t>& dw, 
                                    std::vector<dis_t>& delta_rhs, 
                                    const sol_t& sol, const mesh_t& mesh,
                                    const force_t& force_handle,
                                    const unsigned int i_ele, 
                                    const value_t perturbation_ratio ) {
    const dis_t &dis = sol[i_ele];
    const unsigned int n_moments = dis.size();
    dw.resize( n_moments );
    delta_rhs.resize( n_moments );

    /**
     * 对目前的算例, 即使扰动改变了温度和速度, 也用相同的快速计算公式,
     * 对迭代步数几乎没有影响, 耗时却能省不少 (M=3, 可以省 50%). 因而程
     * 序都改用快速计算公式. 但是也要防止以后遇到的极端情形, 需要考虑更
     * 好的近似 (原来未简化的方式). 因而把原些的程序也留在这里. 要返回
     * 用原来的方式计算, 把本函数中两个 #if 0 改成 #if 1 即可.
     * 
     */

#if 0
    /// 离散方程 R(f) = r, 右端项 r 与 f 无关, 故增量与 r 无关.
    dis_t res_dis( dis.Center(), dis.Scaling(), dis.GetOrder() );
    /// 计算未扰动前的量 rhs-A, 其中 A 是l对流项离散 (rhs =0). 此量做再
    /// 减去外力和碰撞部分, 就是残量.
    Convection( res_dis, dis, sol, mesh, i_ele );
    
    dis_t delta_dis = dis; // dis 是 const, 所以要先复制一份

    /// 矩系数扰动, 不改变展开中心
    typename dis_t::Iterator the_dis = delta_dis.begin();
#else
    typename dis_t::ConstIterator the_dis = dis.begin();
#endif
    for( size_t i=0; i<n_moments; ++i, ++the_dis )
    {
      const index_t& ind = the_dis.GetIndices();
      /// 用这个bool变量来判断扰动是否仍在标准展开下, 若否, 使用常规的残
      /// 量计算方式
      bool is_standard_expansion = true;    
      GetDeltaValue( dw[i], is_standard_expansion, i, 
                     *the_dis, perturbation_ratio, dis );

#if 0
      if( !is_standard_expansion )
      {
        const value_t dis_val_before_ran = *the_dis;
        *the_dis += dw[i]; // 此时, delta_dis 即是扰动后的分布函数
        dis_t tmp( dis.Center(), dis.Scaling(), dis.GetOrder() );
        tmp.ProjectToStdSpace( delta_dis );
        dis_t new_res_dis( tmp.Center(), tmp.Scaling(), tmp.GetOrder() );

        Convection( new_res_dis, tmp, sol, mesh, i_ele );
        tmp.Reinit( res_dis );
        Project( new_res_dis, tmp );

        delta_rhs[i] = res_dis;
        delta_rhs[i] -= tmp; /**
                              * 此时得到的是 B-A, 其中 B 是扰动分布函数
                              * 对应的对流项, A 是扰动前的对流项
                              */

        *the_dis = dis_val_before_ran;
      }
      else
#endif
      {// 可以使用快速计算的部分
        /**
         * 先算对流部分的增量. 由于现在问题空间一维(形参那样一写就已经
         * 确定了), 所以这里的增量的计算就直接代以左右通量增量来算了.
         * 
         */
        static std::vector<value_t> unit_normal(1);
        unit_normal[0] = 1.0;
        // delta_rhs[i] 的初始化也放到这个函数里头了
//        GetDeltaFlux( dis, sol[i_ele+1], unit_normal, dw[i], ind, delta_rhs[i] );
        PROBLEM::GetDeltaFlux( dis, sol[i_ele+1], unit_normal, dw[i], ind, delta_rhs[i] );
        unit_normal[0] = -1.0;
        dis_t dflux;
        // GetDeltaFlux( dis, sol[i_ele-1], unit_normal, dw[i], ind, dflux );
        PROBLEM::GetDeltaFlux( dis, sol[i_ele-1], unit_normal, dw[i], ind, dflux );
        /**
         * 这个操作也占了很多计算量, 因为实际, delta_rhs[i] 有很多非 0
         * 项, 然而这里本来常数的计算量本拉到 N_M(未知量个数) 的计算量,
         * 于是结合外层循环, 就有了 N_M^2 的计算量. 程序分析说明这个函
         * 数以及分布函数迭代器 ++ 操作占据大量时间, 使得总耗时随 M 的
         * 增长远大于时间发展格式随 M 的增长! 所以必须要优化这部分程序
         * (不知 小蔡 优化了 ++ 操作以后这里影响还有多大)!
         * 
         */
        delta_rhs[i] += dflux;
        delta_rhs[i] /= mesh.Length( i_ele );
      }
      /**
       * 当扰动改变标准展开速度和温度时, 外力和碰撞部分仍可用下式计算.
       * 
       */
      force_handle.GetDeltaAcceleration( delta_rhs[i], dis, dw[i], ind, i_ele );
      GetDeltaCollision<CM, dis_t, force_t, _BV>( delta_rhs[i], dis, dw[i], ind, force_handle, i_ele );
//      PROBLEM::GetDeltaAcceleration( delta_rhs[i], dis, dw[i], ind );
      // PROBLEM::GetDeltaCollision( delta_rhs[i], dis, dw[i], ind );
    }
  }

//private:
public:
  /** 
   * @brief 还是把对流项整体单独封装出来
   * 
   * 对流部分实际是 \f[ v \frac{\partial f}{\partial x} \f] 的离散形式,
   * 将其值记为 A, 于是由于残量计算的公式, 本函数实际是把 - A 加到右端
   * 项 @p res_dis 上去并返回.
   * 
   * 关于 A 的计算, 准确地讲, 是向外流量之和. 那么这里的计算我们既然都
   * 是取当前单元流出的通量, 则就应该是两部分相加作为对流的近似, 所以这
   * 里会和我们写出来的离散形式有个符号的差别.
   * 
   * 最后, 这个封装于高维可能并不很好.
   *
   * @param res_dis 返回: -= A, 参数系与 @p dis 同
   * @param dis 当前单元的分布函数, 单独传入是因为允许它是一个扰动后的
   *            分布函数, 要求是标准展开
   * @param sol 解, 用来取相邻单元的分布, 暂时还不考虑相邻单元的扰动
   * @param mesh 网格
   * @param i_ele 当前单元指标
   */
  static void Convection( dis_t& res_dis, const dis_t& dis, 
                          const sol_t& sol, const mesh_t& mesh,
                          const unsigned int i_ele ) {
#ifdef RECONSTRUCT
    OpConvection<mesh_t, sol_t, _BV>::Convection(res_dis, dis, sol, mesh, i_ele);
#else
    dis_t flux_l2r, flux_r2l;
    std::vector<value_t> uno( 1, 1. ); /// 这个信息应该改为由网格提供
    value_t dx = mesh.Length(i_ele);

    dis_t tmp( res_dis.Center(), res_dis.Scaling(), res_dis.GetOrder() );
    if( i_ele > 0 && i_ele < sol.size()-1 )
    {
//      NumericalFlux( sol[i_ele-1], dis, uno, flux_l2r, flux_r2l );
      PROBLEM::NumericalFlux( sol[i_ele-1], dis, uno, flux_l2r, flux_r2l );
      Project( flux_r2l, tmp );
      res_dis.Add( -1./dx, tmp );
//      NumericalFlux( dis, sol[i_ele+1], uno, flux_l2r, flux_r2l );
      PROBLEM::NumericalFlux( dis, sol[i_ele+1], uno, flux_l2r, flux_r2l );
      Project( flux_l2r, tmp );
      res_dis.Add( -1./dx, tmp );
    }
    else if( i_ele == 0 )
    {
      dis_t bv;
      _BV::BoundaryValue( bv, dis, false );
//      NumericalFlux( bv, dis, uno, flux_l2r, flux_r2l );
      // PROBLEM::BoundaryValue( bv, dis, false );
      PROBLEM::NumericalFlux( bv, dis, uno, flux_l2r, flux_r2l );
      Project( flux_r2l, tmp );
      res_dis.Add( -1./dx, tmp );
//      NumericalFlux( dis, sol[i_ele+1], uno, flux_l2r, flux_r2l );
      PROBLEM::NumericalFlux( dis, sol[i_ele+1], uno, flux_l2r, flux_r2l );
      Project( flux_l2r, tmp );
      res_dis.Add( -1./dx, tmp );
    }
    else if( i_ele == sol.size()-1 )
    {
//      NumericalFlux( sol[i_ele-1], dis, uno, flux_l2r, flux_r2l );
      PROBLEM::NumericalFlux( sol[i_ele-1], dis, uno, flux_l2r, flux_r2l );
      Project( flux_r2l, tmp );
      res_dis.Add( -1./dx, tmp );
      dis_t bv;
      _BV::BoundaryValue( bv, dis, true );
//      NumericalFlux( dis, bv, uno, flux_l2r, flux_r2l );
      
      PROBLEM::NumericalFlux( dis, bv, uno, flux_l2r, flux_r2l );
      Project( flux_l2r, tmp );
      res_dis.Add( -1./dx, tmp );
    }
#endif
  }
  /** 
   * @brief 获得扰动量
   *
   * 获得扰动量通常的策略是原始值 (@p _moment) 乘以扰动比例 (@p
   * perturbation_ratio). 由于我们现在考虑矩系数的扰动, 而被扰动分布函
   * 数通常在自身的标准展开下. 这样一阶矩系数是 0, 于是按上式计算的扰动
   * 量也就是 0. 为了避免这一常 0 的出现, 并考虑到一阶矩大致与 \f$
   * \rho u \f$ 同阶, 在考虑一阶矩的扰动时便使用 \f$ \rho u \f$ 作为扰
   * 动的基量. 同样对二阶矩, 采用相同的基量 \f$ \rho \theta \f$ 做扰动.
   *
   * 当扰动这些特殊的矩系数时, 最后获得的扰动分布函数可能改变速度和温度
   * (标准展开参数系变了). 此时残量的扰动增量需要另行处理, 因而本函数还
   * 返回这样一个标志.
   * 
   * @param dw 扰动量
   * @param is_standard_expansion 这个扰动是否会改变标准展开参数系的标识
   * @param _i 扰动的矩系数在一维数组中的指标
   * @param _moment 扰动的基量
   * @param perturbation_ratio 扰动比例
   * @param dis 分布函数, 从这里获得密度速度温度. 事实上, 扰动基量也可
   *            以通过 @p _i 在这里获得
   */
  static void GetDeltaValue( value_t& dw, bool& is_standard_expansion, 
                             const unsigned int _i, const value_t& _moment,
                             const value_t perturbation_ratio,
                             const dis_t& dis ) {
    if( _i>=1 && _i<= dim )
    {// 由于一阶矩系数与 rho * u 同阶, 故用它做扰动基量更好
      dw = perturbation_ratio * dis.front() * dis.Center()[_i-1] + 1e-12;
      is_standard_expansion = false; // 扰动改变平均速度
    }
    else if( _i>dim && _i<=EndSecondMomentsNumber<dim>() )
    { // 同样的道理, 用 rho * theta 作为二阶矩系数的扰动基量. 按照标准
      // 投影的公式, 倒并不是所有二阶矩的扰动都需要标准投影, 不过这个就
      // 不来分别了.
      dw = perturbation_ratio * dis.front() * dis.Scaling() * dis.Scaling() + 1e-12;
      is_standard_expansion = IsStillStdExpansion<dim>( _i );
    }
    else
    {
      dw = perturbation_ratio * _moment + 1e-12;
      is_standard_expansion = true;
    }    
    // // 然后发现, 无视是否改变标准展开参数系, 都使用快速计算方式 (只要下
    // // 面设置 is_standard_expansion 为 true 即可), 对于测试的 Couette
    // // 总的迭代步数都无影响, 而效率更高了. 当然, 若细致起来, 无视的话要
    // // 改前面程序.
    // is_standard_expansion = true;
  }
};

#endif // __RESIDUAL_H__

/**
 * end of file
 * 
 */
