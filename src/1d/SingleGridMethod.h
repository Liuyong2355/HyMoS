/**
 * @file   SingleGridMethod.h
 * @author HU Zhicheng <huzhicheng1986@gmail.com>
 * @date   Wed Dec  4 12:14:01 2013
 * 
 * @brief 稳态 Boltzmann 方程矩方程组给定网格上的基本迭代法
 * 
 * 给定网格上的基本迭代法, 目前想到可以实现的有 SGS - local Newton,
 * global Newton - block SGS (1D 时追赶法), block SGS - local Newton,
 * 已经实现的只有 SGS - local Newton. 将这些方法实现在一个类中还是封装
 * 成不同的类, 这点还在考虑中, 但是这各个方法都需要提供这样一个基本解法
 * 的接口, 以供多重网格法调用. 这个接口便是: 给定一个解, 右端项, 网格,
 * 基本迭代法的步数, 返回一个新的解.
 *
 * 各种迭代法的一个基础是方程残量的计算, 这个封装开来, 由一个计算残量的
 * 类提供.
 *
 * 为提供数据的比较, 这里还提供一个显格式推进的求解器.
 *
 * @note 这里提供的基本迭代法并不考虑外力方程的求解. 换句话说, 这里的求
 *       解器求解的是固定外力下的 Boltzmann 方程. 关于外力的方程的求解,
 *       由额外的程序提供. 
 *  
 * 最后, 所有 Boltzmann 方程的这些求解器 (包括这些基本迭代法, 多重网格
 * 法, 显格式时间推进等) 都将统一提供一个函数接口 (实际参数已有变化): 
 * @code 
 * void Solve( sol_t& sol, const sol_t& rhs, const mesh_t& mesh, const unsigned int steps = 1 ) 
 * @endcode 
 * 以便外力耦合问题的程序实现. 还有一个主要的统一接口是后处理: 
 * @code void PostProcessing( sol_t& sol ) @endcode
 */


#ifndef __SINGLE_GRID_METHOD_H__
#define __SINGLE_GRID_METHOD_H__

#include <memory>
#include <vector>
#include <sstream>
#include <omp.h>

#include <Eigen/Dense>

#include <NRxx/Distribution.h>
#include <NRxx/Transform.h>
#include <NRxx/SpectralForODE.h>
#include <NRxx/NumericalFlux.h>

#include "ConfigFile.h"
#include "tools/tools.h"

/**
 * 声明几个本文件中使用的公共函数, 其定义在此文件最后.
 * 
 */

template <typename DISTRIBUTION, typename _VEC, typename PROBLEM>
void _PreservingPositivityStepLength( double &mu_max, const DISTRIBUTION &dis, const _VEC &_df );

/**
 * 函数声明结束
 * 
 */


/**
 * @class SingleGridMethod
 * @brief 单层网格基本迭代法基类
 *
 * 想来还是将各种方法封装开来, 这样就该将相同部分剥离出来. 这个基类不能
 * 直接使用.
 *
 * 暂时, 基类只包含两个数据成员以及两个纯虚的接口函数等. 实现都放到唯一
 * 的继承类中, 待实现其他继承类时, 再将相同部分移到此类.
 *
 * 似乎改用几个指针存储 Solve 传进来的参数 sol, rhs, mesh, 这样后头函数
 * 就可以少写参数了.
 *
 * 各种需要的类型都可以从残量计算类 _RESIDUAL 中获取, 所以这里将模板参
 * 数改成一个:
 * - RESIDUAL: 方程 (有无外力, 各种碰撞项) 残量计算, 提供一个残量计算接
 *        口 (这个接口对高维可能不好, 可能残量还是在此类中组装为是), 还
 *        有各种需要的类型
 */
template <typename _RESIDUAL>
// template <typename DISTRIBUTION, typename MESH, typename _RESIDUAL, typename PROBLEM>
class SingleGridMethod {
public:
  typedef _RESIDUAL RESIDUAL;
  typedef typename RESIDUAL::mesh_t mesh_t;
  typedef typename RESIDUAL::sol_t sol_t;
  typedef typename RESIDUAL::force_t force_t;
  typedef typename RESIDUAL::PROBLEM PROBLEM;

  typedef typename sol_t::dis_t dis_t;
  typedef typename dis_t::value_type value_t;
  typedef typename dis_t::Indices index_t;
  typedef typename dis_t::Velocity velocity_t;
  
  typedef std::shared_ptr<mesh_t> ptr_mesh_t;
  typedef std::shared_ptr<force_t> ptr_force_t;

protected:
  double perturbation_ratio;    /**< 计算 Jacobian 矩阵数值微分扰动比例 */
  double lambda;                /**< 正则化 Jacobian 矩阵的比例系数 */
  double mu;                    /**< 更新解步长, 准确的 Newton 迭代应为 1 */
  double tolerance;             /**< 判断是否达到稳态 (收敛) 的残量忍量 */

public:
  SingleGridMethod() : perturbation_ratio(0.10), lambda(2.0),
                       mu(1.0), tolerance(1e-6) { }

  SingleGridMethod( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) {
    read_config( p_mesh, cf );
  }
  SingleGridMethod( const ptr_mesh_t& p_mesh, std::istream& is ) {
    is.read((char *)&perturbation_ratio, sizeof(double));
    is.read((char *)&lambda, sizeof(double));
    is.read((char *)&mu, sizeof(double));
    is.read((char *)&tolerance, sizeof(double) );
  }

  void Reinit( const ptr_mesh_t& p_mesh ) { };

  void read_config( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) {
    perturbation_ratio = (double)cf.Value( "SingleGridMethod", "Perturbation_ratio" );
    lambda = (double)cf.Value( "SingleGridMethod", "RegularizeJacobiMatrixScaling" );
    mu = (double)cf.Value( "SingleGridMethod", "UpdatingDisStepScaling");

    tolerance = (double)cf.Value( "Error", "Tol" );
  }
  void print_config() const {
    std::cout << "SingleGridMethod: perturbation_ratio = " << perturbation_ratio << std::endl;
    std::cout << "SingleGridMethod: lambda = " << lambda << std::endl;
    std::cout << "SingleGridMethod: mu = " << mu << std::endl;
//    std::cout << "SingleGridMethod: tolerance = " << tolerance << std::endl;
  };
  void dump_data( std::ostream& os ) const {
    os.write((const char *)&perturbation_ratio, sizeof(double));
    os.write((const char *)&lambda, sizeof(double));
    os.write((const char *)&mu, sizeof(double));
    os.write((const char *)&tolerance, sizeof(double) );
  }
  void reset_tol( const double _tol ) { 
    tolerance = _tol;
  }
  double& reset_tol() {
    return tolerance;
  }
  const double& reset_tol() const {
    return tolerance;
  }


  void SetParameters( const double _ratio, const double _lambda, 
                      const double _mu, const double _tolerance ) {
    SetParameters( _ratio, _lambda, _mu );
    reset_tol( _tolerance );
  }

  void SetParameters( const double _ratio, const double _lambda, const double _mu ) {
    perturbation_ratio = _ratio;
    lambda = _lambda;
    mu = _mu;
  }

  /** 
   * @brief 给定网格上基本迭代法求解到稳态.
   *
   * 也即本函数求解直到残量小于 @p tolerance 时为止. 于是为了避免这个函
   * 数求解的耗时过多, 在多重网格方法调用中, 唯有单元数很少的时候才适合
   * 调用此函数.
   *
   * 这个函数放在基类中未必恰当, 因为全局 Newton 迭代最初粗网格上并不需
   * 要此函数. 还有, 这个函数放在基类中, 函数名就不能是 Solve 了, 否则
   * 由于和虚函数同名, 这个重载函数将在派生类中隐藏.
   * 
   * @param sol 
   * @param rhs 
   * @param p_mesh 
   * @param p_force_handle 外力项处理句柄, 提供 @p mesh 上的外力等信息
   * @param max_steps 当作为 MG 的最粗网格调用时, 最好有这个参数, 但为
   *                  同时顾及单独使用情况, 给一个默认参数为 0. 当然目
   *                  前在 SteadyState.h 中组的框架并用不到此函数.
   */
  void ExactSolve( sol_t& sol, const sol_t& rhs, 
                   const ptr_mesh_t& p_mesh, 
                   const ptr_force_t& p_force_handle,
                   const unsigned int max_steps=0 ) {
    double initial_res;
    RESIDUAL::GetResidualPL2Norm( initial_res, sol, rhs, *p_mesh, *p_force_handle );
    // 以下的设置是不是效率又防止最后残量降不下来.
    const double _Tol = (initial_res > 10*tolerance) ? tolerance : 0.1 * tolerance;

    double res;
    unsigned int step = 0;
    do {
      ++step;
      Solve( sol, rhs, p_mesh, p_force_handle, 1 );
//      RESIDUAL::GetResidualMaxNorm( res, sol, rhs, *p_mesh, *p_force_handle );
      RESIDUAL::GetResidualPL2Norm( res, sol, rhs, *p_mesh, *p_force_handle );
/// 根据判断原则, 单独使用而非在 MG 使用时, while 的判断效率略低些.
    } while (res >= _Tol && (step<max_steps || max_steps==0) );
//    std::cout<< "coarsest steps: = " << step << ", res = " << res << std::endl;
  }

  /** 
   * @brief 纯虚函数: 给定网格基本迭代法函数接口
   *
   * @param sol 解, 传入初值, 返回新解. 暂时考虑传入返回的解在各单元上
   *            都是标准展开
   * @param rhs 方程右端项. 最密网格上是 0, 粗网格上则是细网格上残量的
   *            投影, 各单元上展开参数系与 @p sol 相同, 也即未必是标准
   *            展开. 由于迭代一次后, rhs 和 sol 的展开参数系就不同了,
   *            所以这个要求未必需要, 还是在方程求解时多做一步投影罢
   *            (rhs 在 sol 参数系重新展开)
   * @param p_mesh 求解的网格
   * @param p_force_handle 外力项处理句柄, 提供 @p mesh 上的外力等信息
   * @param steps 基本迭代法求解步数. 注意继承类实现时需要设定重新设定
   *              默认值为 1
   */
  virtual void Solve( sol_t& sol, const sol_t& rhs, 
                      const ptr_mesh_t& p_mesh, 
                      const ptr_force_t& p_force_handle,
                      const unsigned int steps = 1 ) = 0;

  /** 
   * @brief 纯虚函数: 基本迭代法求解的后处理
   * 
   * 在 SGS-local Newton 法中, 由于迭代法不保证守恒性, 故需要一个后处理
   * 保证此性质. 其他方法可能也会需要, 所以用一个统一的函数调用.
   *
   * @param sol 解
   * @param p_mesh 网格
   */
  virtual void PostProcessing( sol_t& sol, const ptr_mesh_t& p_mesh ) = 0;

  /** 
   * @brief 初值的准备
   * 
   * 也是统一的接口, 除了 Full multi grid 需要用一个算法获得初值, 其余
   * 方法的初值都直接由问题给定.
   * 
   * @param sol 
   * @param p_mesh 
   * @param p_force_handle 
   */
  void initial_value( sol_t& sol, const ptr_mesh_t& p_mesh, 
                      const ptr_force_t& p_force_handle ) {
    PROBLEM::initial_value( sol, *p_mesh );
  }
};   // class SingleGridMethod



/**
 * @class SLGSNewton
 * @brief 单层网格基本迭代法: symmetric local Gauss-Seidel - Newton
 *
 * 这是目前实现的基本迭代法. 继承自 SingleGridMethod 类
 *
 * @note 考虑多重网格, 所以这里先只把原些程序单线程部分取来. 若考虑多线
 * 程, 多线程参似乎封装到 mesh_t 中为好?
 */
template <typename _RESIDUAL>
class SLGSNewton : public SingleGridMethod<_RESIDUAL> {
public:
  /**
   * 似乎在模板类继承中, c++ 标准不能直接使用基类的 typedef 的类型, 数
   * 据成员和函数, 而是需要使用 this-> 或者 using 基类, 或者就是基类
   * BASE::. 但是如此派生类对象要调用时可怎么办?
   * 
   */
  typedef SingleGridMethod<_RESIDUAL> _BASE;
//  using _BASE;
  typedef _RESIDUAL RESIDUAL;  // 多重网格类又是不能直接取 RESIDUAL
  typedef typename RESIDUAL::mesh_t mesh_t;
  typedef typename RESIDUAL::sol_t sol_t;
  typedef typename RESIDUAL::force_t force_t;
  typedef typename RESIDUAL::PROBLEM PROBLEM;

  typedef typename sol_t::dis_t dis_t; // 或者 _BASE::sol_t, 
  typedef typename dis_t::value_type value_t;
  typedef typename dis_t::Indices index_t;
  typedef typename dis_t::Velocity velocity_t;

  typedef std::shared_ptr<mesh_t> ptr_mesh_t;
  typedef std::shared_ptr<force_t> ptr_force_t;

  SLGSNewton() : _BASE() {  };
  SLGSNewton( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) : _BASE(p_mesh, cf) { };
  SLGSNewton( const ptr_mesh_t& p_mesh, std::istream& is ) : _BASE( p_mesh, is) { }

public:
  /// 函数接口: 继承的虚函数
  void Solve( sol_t& sol, const sol_t& rhs, 
              const ptr_mesh_t& p_mesh,
              const ptr_force_t& p_force_handle,
              const unsigned int steps = 1 ) {
    for( size_t i=0; i<steps; ++i )
    {
      SymmetricGaussSeidel( sol, rhs, *p_mesh, *p_force_handle );
      p_force_handle->CalcExternalForce( sol );
    }
  }
  
  /// 由于问题不同, 保持平均密度为 1, 有些苛刻, 还是让用户提供此函数罢
  void PostProcessing( sol_t& sol, const ptr_mesh_t& p_mesh ) {
    PROBLEM::PostProcessing( sol, *p_mesh );
  }
    
private:
  /**
   * @brief symmetric Gauss-Seidel 迭代法的一步
   * 
   */
  void SymmetricGaussSeidel( sol_t& sol, const sol_t& rhs, 
                             const mesh_t& mesh, 
                             force_t& force_handle ) {
    GaussSeidelOneStep( sol, rhs, mesh, force_handle );
    ReverseGaussSeidelOneStep( sol, rhs, mesh, force_handle );
  }
  /**
   * @brief 按单元增序一步 Gauss-Seidel 迭代
   * 
   */
  void GaussSeidelOneStep( sol_t& sol, const sol_t& rhs,
                           const mesh_t& mesh, 
                           force_t& force_handle ) {
    for( size_t k=0; k<sol.size(); ++k ){
      GaussSeidelOneStep( sol, rhs, mesh, force_handle, k );
      force_handle.CalcExternalForce( sol, k );
    }      
  }
  /**
   * @brief 按单元递减序一步 Gauss-Seidel 迭代
   * 
   */
  void ReverseGaussSeidelOneStep( sol_t& sol, const sol_t& rhs, 
                                  const mesh_t& mesh,
                                  force_t& force_handle ) {
    for( size_t k=sol.size()-1; k != 0; --k ){
      GaussSeidelOneStep( sol, rhs, mesh, force_handle, k );
      force_handle.CalcReverseExternalForce( sol, k );
    }      
  }
  /** 
   * @brief Gauss-Seidel 迭代某个单元上的求解
   * 
   * @param k 单元指标
   */
  void GaussSeidelOneStep( sol_t& sol, const sol_t& rhs,
                           const mesh_t& mesh,
                           const force_t& force_handle,
                           const unsigned int k ) {    
    dis_t &dis = sol[k];
    /**
     * 一维的时候, 未知量就是 order+1:
     *
     *   const size_t order = dis.GetOrder();
     *
     * 考虑维数无关就不应这样了.
     */
    const unsigned int n_moments = dis.size();
    Eigen::MatrixXd jacobi_matrix( n_moments, n_moments );
    
    /**
     * 下面是一个局部单元上的 Newton 迭代. 其迭代步数选取需要综合考量该
     * 局部问题的求解精度和整体的计算效率. 数值实验表明, 这个局部迭代精
     * 确求解 (增加迭代步数) 并不会带来整体 SGS 迭代步数的显著减少, 反
     * 而大大增加计算时间. 并且每次 SGS 仅做一步 Newton 迭代效果便已不
     * 错(大概最省时). 但也不甚放心, 觉得有时候多几步也是有好处的. 于是
     * 采用如下的标准判断 Newton 迭代是否终止: 1. 该单元残量小于
     * tolerance (比全局tolerance 更小为好), 此时认为该单元已是局部稳态,
     * 不需再算; 2. 或者残量已下降到比初始残量的一半或以下, 通常一步就
     * 可以了; 3. 达到预设最大迭代步数.
     * 
     * 使用 BackTracking 选更新步长以后, 下面这几行就不大好放到 for 循
     * 环中了.
     */
    dis_t res_dis( dis.Center(), dis.Scaling(), dis.GetOrder() );
    Project( rhs[k], res_dis );

    /// 切记: res_dis = rhs - R(f)
    RESIDUAL::GetResidual( res_dis, sol, mesh, force_handle, k );
//      res_norm = L1NormMomentsValue( res_dis );
    double res_norm = PL2NormDis( res_dis );
    double initial_res = res_norm;
    // 为了与 ExactSolve() 中一致, 否则将是死循环
    const double _Tol = (initial_res > 10*_BASE::tolerance) ? _BASE::tolerance : 0.1 * _BASE::tolerance;
    
    double _factor = 1; bool is_success = true;
    static const unsigned int max_newton_steps = 5;
    for( size_t _i=0; _i<max_newton_steps; ++_i )
    { /// 就单层网格的迭代法言, 这里可能只算一步效率更好
      if( ((res_norm < initial_res * 0.5) || (res_norm < _Tol)) )
      {/**
        * 这是局部 newton 迭代终止的两个条件, 第一个说明残量已经减半,第
        * 二个说明已经局部稳态. 局部稳态的 tolerance 如果采用与全局相同
        * 的 tolerance 的话, 将有可能影响最后几步的收敛速度. 特别是在多
        * 重网格方法调用中, 前后光滑步是 1 的时候, 算例在最后几步残量甚
        * 至会暂时增大. 当这里的 tolerance 比全局更小的时候 (如全局的
        * 0.1), 该算例在最后几步残量也是稳固减小. 或者就是让此迭代至少
        * 做一次, 这也可以保证残量稳固减小.
        *
        * 不过, 上述两中调整方式, 在多重网格迭代法前后光滑3次和2次时,
        * 总的迭代步都回多了一步, (如光滑3次, 15步时比tolerance 还是大
        * 一点点, 要16步才收敛), 因而耗时多了十来秒. 对于前后光滑1次,
        * 耗时没有基本没有变化, 而由于残量的稳固减小, 总迭代步从 30 降
        * 为 28.
        * 
        * 以上算例使用 V-型多重网格方法, 当使用 W-型时, 不管何种调整方
        * 式, 都将使得耗时较多的增加 (虽然总迭代步可能减少).
        * 
        * 这里做一个调整, 让程序至少迭代一步 (会有点效率损失). 记得当光
        * 滑次数大于 1 时, 可以去掉此调整.
        *
        * 新跑 couette 算例, M=6, N=2048, tol=1e-8, V-cycle, 前后光滑
        * 1 次时, 当残量到 1.5e-8 之后迅速减慢, 而且趋势是会稳定在某个
        * 大于 1e-8 的残量上. 调小局部 tol 或至少迭代一次都行不通. 所以
        * 还是尽量不要用 smoothing = 1 做吧, 虽然这样省时. 这好难理解.
        */
        break;
      }

      GetJacobiMatrix( jacobi_matrix, sol, mesh, force_handle, k );
      /**
       * @brief 用残量修正 Jacobian 矩阵
       *
       * 将散射项计算改对以后, 上面得到的 jacobi_matrix 使用 gauss_jordan求
       * 逆容易出现 nan, 所以需要强行做到主对角元占优, 这就是一个局部松弛技
       * 术. 具体的就是把线性方程组修改为 \f[ \left( \Lambda +
       * \frac{\delta R(f^{(m)})}{\delta f} \right) \Delta f^{(m)} = -
       * R(f^{(m)}), \f] 其中 \f$ \Lambda = \lambda I \f$, \f$ I \f$ 是单位
       * 阵, \f$ \lambda \f$ 是松弛参数, 参考建议是 \f[ \lambda \propto
       * \Vert R(f^{(k,m)}) \Vert. \f] 这有点相当于取时间步长为 \f$
       * \frac{1}{\lambda} \f$ 的概念, 于是残量越小, 松弛越少, 对应就可以说
       * 时间步长可以取得更大.
       *
       * 至于这里范数的取法, 可能需要去测试. 或者可以进一步细化为 \f[
       * \Lambda = \lambda R(f^{(k,m)}), \f] 这样的效果也待检测.
       *
       * 这里的 lambda 取值会影响迭代步数和收敛速度, 需要按当前算例调节.
       */
      jacobi_matrix.diagonal().array() += _factor * _BASE::lambda * res_norm;

	      if ( !is_success )
	      {
	        std::cerr << "Local Newton iteration did not converge at element "
	                  << k << ", position " << mesh.Center(k) << std::endl;
	      }

      /// 求逆
      jacobi_matrix = jacobi_matrix.inverse().eval();

#if 0 // 看矩阵求逆是否出了问题
      for( size_t i=0; i<n_moments; ++i )
        for( size_t j=0; j<n_moments; ++j )
        {
          if( isnan(jacobi_matrix(j,i)) )
          {
            std::cout << "matrix ele nan found after gauss_jordan in elements k =" 
                      << k << ", j = " << j << ", i = " << i << std::endl;
//          getchar();
          }
        }
#endif

      Eigen::VectorXd deal_res( n_moments );
      typename dis_t::Iterator it_res = res_dis.begin();
      const typename dis_t::Iterator it_res_end = res_dis.end();
      for( size_t i=0; it_res != it_res_end; ++it_res, ++i )
        deal_res( i ) = *it_res;
      Eigen::VectorXd deal_df( n_moments );

      /**
       * 现在 GetResidual() 的计算方式, 这里求解的 deal_df 正好是 \f$
       * \Delta f^{(m)} \f$, 而不是以前的 \f$ -\Delta f^{(m)} \f$. 所以
       * 要注意后面程序的相应修改.
       * 
       */
      deal_df = jacobi_matrix * deal_res; 

      /// 更新步长的选取
      double mu_max = _BASE::mu;
      GetDisUpdateStepLength( mu_max, dis, deal_df );
      if( mu_max < dis.front() * 1.0e-14  ) /// 这个时候认为解不该再有变化.
        break;  

      // UpdateDistribution( dis, mu_max, deal_df );

      is_success = BackTracking( res_norm, res_dis, sol, mu_max, deal_df, rhs, mesh, force_handle, k );

      if ( !is_success )
      {
	        std::cerr << "Backtracking failed at element " << k
	                  << ", position " << mesh.Center(k)
	                  << "; applying a positivity-preserving fallback step."
	                  << std::endl;
#if 0
        // _factor *= 10;
        // getchar();
//        continue;
#else 
            value_t C_M = Transform<value_t>::MaxRootOfHermitePolynomial( dis.GetOrder()+1 );
    // 下式位置空间一维时才正确
            double local_dt = 0.5 * mesh.Length(k) / (fabs( dis.Center()[0] ) + C_M * dis.Scaling());
//CFL * _grid_size / lambda;
            _PreservingPositivityStepLength<dis_t,dis_t,PROBLEM>( local_dt, dis, res_dis);
            dis_t new_dis( dis.Center(), dis.Scaling(), dis.GetOrder() );
            dis.Add( local_dt, res_dis );
            new_dis.ProjectToStdSpace( dis );
            std::swap( dis, new_dis );
            
            res_dis.Reinit( dis );
            Project( rhs[k], res_dis );
            RESIDUAL::GetResidual( res_dis, dis, sol, mesh, force_handle, k );
            res_norm = PL2NormDis( res_dis );
#endif
      }
    }
  }
  /** 
   * @brief local jacobian 矩阵
   * 
   * 使用数值微分计算此局部 jacobian 矩阵. 数值实验表明, jacobian 矩阵
   * 的计算不能太粗糙.
   *
   * 按, 数值微分的计算方式必须优化, 也就是扰动残量需要快速计算, 而不是
   * 仍然调用通常的残量计算公式. 具体原因是: 计算结果已经显示, 矩数M=3
   * 时, 简单的 SGS + Newton 总的迭代步数比时间发展少 1/10, 但是耗时更
   * 久. 细细想来也是合理, 因为如果调用常规的残量计算方式 SGS + Newton
   * 方式需要计算残量的总次数比时间发展多三四十倍....
   * 所以呢, 必须优化, 把 jacobi 矩阵的计算时间降下来.
   *
   * @param jacobi_matrix 返回值, 要求矩阵阶数已由外部初始化好.
   *
   * @param sol 解
   * @param mesh 网格
   * @param force_handle 外力项处理句柄, 提供 @p mesh 上的外力等信息
   * @param i_ele 计算 jacobian 矩阵的当前单元指标
   */
  void GetJacobiMatrix( Eigen::Matrix<value_t, Eigen::Dynamic, Eigen::Dynamic>& jacobi_matrix, 
                        const sol_t& sol, const mesh_t& mesh,
                        const force_t& force_handle,
                        const unsigned int i_ele ) {
    std::vector<value_t> dw;
    std::vector<dis_t> delta_rhs;
    if( i_ele >=1 && i_ele < sol.size()-1 )
    {/**
      * 关于使用这个快速计算公式和旧的非快速计算方式的比较, 可见
      * essay/fast-residual.tex 中的相关记录.
      */
      RESIDUAL::GetFastDeltaResidual( dw, delta_rhs, sol, mesh, force_handle, i_ele, _BASE::perturbation_ratio );
    }
    else
    {/**
      * 由于 ghost cell 分布函数依赖于边界单元的分布函数. 对于边界单元
      * jacobi 矩阵的计算, 快速计算公式并不容易推导, 因而仍采用如下方式
      * 计算: 先计算扰动后的残量, 再与扰动前残量做差. 这个实现在下面这
      * 个函数中. 
      *
      * 当然考虑到 ghost cell 的分布函数仅影响残量的通量部分, 故而外力
      * 项和碰撞项带来的扰动增量也可以参考 GetFastDeltaResidual() 函数
      * 那样计算. 也就是说仅对流部分使用先计算扰动前后对流项值再做差.
      * 
      * 毕竟 GetDeltaResidual() 这个函数仅仅边界单元才调用, 而边界单元
      * 毕竟少数, 所以目前未做修改, 还是计算扰动前后整体的残量, 再做差
      * 的方式.
      * 
      */
      RESIDUAL::GetDeltaResidual( dw, delta_rhs, sol, mesh, force_handle, i_ele, _BASE::perturbation_ratio ); 
    }

    for( size_t i=0; i<delta_rhs.size(); ++i ){
      typename dis_t::Iterator it_rhs = delta_rhs[i].begin();
      const typename dis_t::Iterator it_rhs_end = delta_rhs[i].end();
      for( size_t j=0; it_rhs != it_rhs_end; ++ it_rhs, ++j ){
        jacobi_matrix( j, i ) = *it_rhs / dw[i];
      }
    }

  }
  
  /** 
   * @brief 更新解后密度和温度恒正的步长选取
   *
   * @param mu_max 步长, 返回值
   * @param dis 原分布函数
   * @param deal_df \f$ \Delta f \f$
   */
  void GetDisUpdateStepLength( value_t &mu_max, const dis_t &dis, 
                               const Eigen::Matrix<value_t, Eigen::Dynamic, 1> &deal_df ) {
    _PreservingPositivityStepLength<dis_t, Eigen::Matrix<value_t, Eigen::Dynamic, 1>, PROBLEM>( mu_max, dis, deal_df );
  }

  /** 
   * @brief 根据步长和线性方程组解向量更新单元上的分布函数
   * 
   * 本函数操作实际就是 new_dis += mu_max * delta_df 
   *
   * @param new_dis 需要初始化为用来计算 jacobi 矩阵的当前分布, 亦是返回值.
   * @param mu_max 更新步长, @see GetDisUpdateStepLength 函数
   * @param deal_df 解向量的增量, 就是矩系数 \$ [f_0,f_{e_1}, \ldots,
   *                f_{Me_3}] \$
   *
   * 做 newton 迭代的时候, 有两种方式选择分布函数依赖的变量. 其一就是这
   * 里展开的矩系数, 其二则是用展开中心 u 和 scaling 的温度替换矩系数中
   * 的 f_{e_i}, i=1,2,3 和 f_{2e_1}. 总算厘清, 第二种方式仍无法避免本
   * 函数最后的标准投影 (速度一维除外), 使用第二种方式也就是 newton 迭
   * 代求 jacobian 矩阵时省了四次标准投影 (三次一阶矩的扰动, 一次温度的
   * 扰动). 
   * 
   * 本质上来说, 不管以哪组变量求解, 结果没有区别, 所以还是最终确定下来,
   * 就直接使用矩系数做为变量了!
   */
  void UpdateDistribution( dis_t& new_dis, const double mu_max, 
                           const Eigen::VectorXd& deal_df ) {
//    dis_t tmp_dis = new_dis;
    typename dis_t::Iterator the_dis = new_dis.begin();
    const typename dis_t::Iterator end_dis = new_dis.end();
    for( size_t i=0; the_dis != end_dis; ++the_dis, ++i )
    {
      *the_dis += mu_max * deal_df( i );
    }
    /// 投影回标准空间  
//    new_dis.ProjectToStdSpace( tmp_dis );
    dis_t tmp_dis( new_dis.Center(), new_dis.Scaling(), new_dis.GetOrder() );
    tmp_dis.ProjectToStdSpace( new_dis );
    std::swap( new_dis, tmp_dis );
  }

  /** 
   * @brief 松弛 Newton 迭代参数选取
   *
   * 在用多重网格求解 Poisson 耦合的半导体器件问题时, 计算总是失败, 主
   * 要原因应该就是原些的步长选取方式, 使得某些单元局部 Newton 迭代发散
   * 了.
   *
   * 按照理论, Newton 方向总是一个残量下降的方向, 所以通过不断缩短步长,
   * 总能找到一个新的解, 使得残量是下降的. 这个技术称为线性搜索, 或者回
   * 溯 (BackTracking). 这里采用简单的减半步长搜索方式.
   *
   * 这个函数实际是通过测试新的解是否已经使得残量下降来判断我们要不要这
   * 个新的解, 所以程序结束之后, 解已经更新, 新的残量也已经算出来, 于是
   * 把这些信息返回 Newton 迭代函数, 可以节省一些计算量
   * 
   * @param res_norm 残量范数, 传进旧残量范数, 返回新残量的范数
   * @param res_dis 新残量, 返回值
   * @param sol 解, 其中当前单元的解一般会被更新
   * @param mu_max 初始更新步长 (已经做了温度, 密度保正)
   * @param deal_df 更新方向, 或者 newton 方向
   * @param rhs 
   * @param mesh 
   * @param force_handle 
   * @param k 
   * 
   * @return 若是 false, 说明回溯失败, 提示迭代可以放弃这个单元的
   *         Newton 了
   */
  bool BackTracking( double& res_norm, dis_t& res_dis, sol_t& sol, 
                     double& mu_max, const Eigen::VectorXd& deal_df,
                     const sol_t& rhs, const mesh_t& mesh, 
                     const force_t& force_handle, 
                     const unsigned int k ) {
    /**< 为了效率起见, 设置一个最大的回溯步数. 步长减半搜索的方式, 5 步
     * 以后, 步长已经是初始的 2^(-4) = 0.0625. 这个数已经很小了, 所以如
     * 果还不成功, 也确实可以放弃了. */
    static const unsigned int max_steps = 5; 

    bool is_success = false;
    for( unsigned int i=0; i<max_steps; ++i )
    {
      dis_t new_dis = sol[k];
      UpdateDistribution( new_dis, mu_max, deal_df );

      dis_t _res_dis( new_dis.Center(), new_dis.Scaling(), new_dis.GetOrder() );
      Project( rhs[k], _res_dis );
      RESIDUAL::GetResidual( _res_dis, new_dis, sol, mesh, force_handle, k );
      double new_res_norm = PL2NormDis( _res_dis );

      if( new_res_norm > res_norm ){
        mu_max *= 0.5;
        continue;
      }
      else{
        res_norm = new_res_norm;
        std::swap( new_dis, sol[k] );
        std::swap( _res_dis, res_dis );
        is_success = true;
        break;
      }
    }

    return is_success;    
  }
};  // class SLGSNewton


/**
 * @class ExplicitTimeForward
 * @brief 显格式时间推进法求解器
 *
 * 此方法主要用于为稳态求解器提供参考数值解以及计算效率的比较数据. 
 *
 * 这个不适合从 SingleGridMethod 类派生, 但为了调用的方便, 仍然提供两个
 * 相同的函数接口 Solve() 和 PostProcessing().
 *
 * 而为此函数添加一个 ExactSolve() 函数接口以后, 就可用于 MultiGrid 方
 * 法了 (为了清爽, 需要注释掉其中打印的信息). 作为稳态求解器, 这个算法
 * 相当于在全局时间步长下的一个 jacobi 迭代 (文献中可能更多的叫
 * Richardson 迭代). 由于不需要抓住每个时间步的信息, 这个算法可以使用局
 * 部时间步长等技术加速, 为了界面清爽, 还是重新写两个相关的类: 
 * @see
 * - ExplicitRichardson
 * - ExplicitSGS
 *
 * @note 
 *
 * - 从这个类来看, 或许所有的求解器都应该存下这三个指针: 解, 方程
 * 右端项, 解所在的网格. 这样调用函数可以少传递参数, 也便于统一接口函数
 * 的实现.
 *
 * - 注意到检测是否达到稳态需要计算残量, 显格式推进实际就是将此残量加到
 * 解上去. 若每个时间步都检查是否达到稳态, 则如今的封装将使得残量的计算
 * 次数加倍. 也就是效率上会差一些, 不过从测试数据来看, 总的时间并没有加
 * 倍, 说明还有其他部分计算量占了很多.
 *
 * - 其实还应该把类中参数设置的接口也统一起来, 比如所有参数都从一个参数
 *   类 Para 中读取, 这过些时候考虑.
 */
template <typename _RESIDUAL>
// template <typename DISTRIBUTION, typename MESH, typename RESIDUAL, typename PROBLEM>
class ExplicitTimeForward {
public:
  typedef _RESIDUAL RESIDUAL;  // 多重网格类又是不能直接取 RESIDUAL
  typedef typename RESIDUAL::mesh_t mesh_t;
  typedef typename RESIDUAL::sol_t sol_t;
  typedef typename RESIDUAL::force_t force_t;
  typedef typename RESIDUAL::PROBLEM PROBLEM;

  typedef typename sol_t::dis_t dis_t; 
  typedef typename dis_t::value_type value_t;
  typedef typename dis_t::Indices index_t;
  typedef typename dis_t::Velocity velocity_t;

  typedef std::shared_ptr<mesh_t> ptr_mesh_t;
  typedef std::shared_ptr<force_t> ptr_force_t;

private: 
  double CFL;
  double tolerance;             /**< 判断是否达到稳态 (收敛) 的残量忍量 */

public:
  ExplicitTimeForward() : CFL(0.5), tolerance(1e-6) {}
  ExplicitTimeForward( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) {
    read_config( p_mesh, cf );
  }
  ExplicitTimeForward( const ptr_mesh_t& p_mesh, std::istream& is ) {
    is.read((char *)&CFL, sizeof(double));
    is.read((char *)&tolerance, sizeof(double));
  }

  void Reinit( const ptr_mesh_t& p_mesh ) { };

  void read_config( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) {
    CFL = (double)cf.Value( "Time", "CFL" );
    tolerance = (double)cf.Value( "Error", "Tol" );
  }

  void print_config() const {
    std::cout << "ExplicitTimeForward: CFL = " << CFL << std::endl;
    std::cout << "ExplicitTimeForward: Tol = " << tolerance << std::endl;
  };

  void dump_data( std::ostream& os ) const {
    os.write((const char *)&CFL, sizeof(double));
    os.write((const char *)&tolerance, sizeof(double));
  }


  void reset_tol( const double _tol ) { 
    tolerance = _tol;
  }
  double& reset_tol() {
    return tolerance;
  }
  const double& reset_tol() const {
    return tolerance;
  }

  void ExactSolve( sol_t& sol, const sol_t& rhs, 
                   const ptr_mesh_t& p_mesh, 
                   const ptr_force_t& p_force_handle, 
                   const unsigned int max_steps=0 ) {
    double initial_res;
    RESIDUAL::GetResidualPL2Norm( initial_res, sol, rhs, *p_mesh, *p_force_handle );
    // 以下的设置是不是效率又防止最后残量降不下来.
    const double _Tol = (initial_res > 10*tolerance) ? tolerance : 0.1 * tolerance;

    double res;
    unsigned int step = 0;
    do {
      ++step;
      Solve( sol, rhs, p_mesh, p_force_handle, 1 );
//      RESIDUAL::GetResidualMaxNorm( res, sol, rhs, *p_mesh, *p_force_handle );
      RESIDUAL::GetResidualPL2Norm( res, sol, rhs, *p_mesh, *p_force_handle );
    } while (res >= _Tol && (step<max_steps || max_steps==0) );
//    std::cout<< "coarsest steps: = " << step << ", res = " << res << std::endl;
  }


  /** 
   * @brief 时间发展一步的函数接口
   *
   * @param sol 解, 传入初值, 返回新解. 暂时考虑传入返回的解在各单元上
   *            都是标准展开
   * @param rhs 方程右端项. 最密网格上是 0, 粗网格上则是细网格上残量的
   *            投影, 各单元上展开参数系与 @p sol 相同, 也即未必是标准
   *            展开. 由于迭代一次后, rhs 和 sol 的展开参数系就不同了,
   *            所以这个要求未必需要, 还是在方程求解时多做一步投影罢
   *            (rhs 在 sol 参数系重新展开)
   * @param p_mesh 求解的网格
   * @param p_force_handle 外力项处理句柄, 提供 @p mesh 上的外力等信息
   * @param steps 基本迭代法求解步数
   */
  void Solve( sol_t& sol, const sol_t& rhs, 
              const ptr_mesh_t& p_mesh, 
              const ptr_force_t& p_force_handle,
              const unsigned int steps = 1 ) {
    static double t = 0;
	 double T = 0.3;
    for( unsigned int i=0; i<steps; ++i )
    {
      // 时间步长
      double dt = TimeStepLength( sol, *p_mesh );
		//double dt = 0.0002;
      t += dt;
      //std::cout << "t = " << t << ", dt = " << dt << std::endl;
		
      // 残量计算
      sol_t res_sol( sol.size() );
      RESIDUAL::GetResidual( res_sol, sol, rhs, *p_mesh, *p_force_handle );
      // 解的更新
      ForwardStep( sol, res_sol, dt );

      p_force_handle->CalcExternalForce( sol );

		//if (t>=T && T < 0.30004){
		//	 std::stringstream file_dis;
		//    file_dis << "shocktube_T.dat";
		//	 SaveMacroVars(sol, *p_mesh, file_dis.str());
		//}
    }
  };

  /** 
   * @brief 基本迭代法求解的后处理
   * 
   * 显格式时间推进是守恒型的, 因而这里不需处理
   *
   * @param sol 解
   * @param p_mesh 网格
   */
  void PostProcessing( sol_t& sol, const ptr_mesh_t& p_mesh ) { };  

  /** 
   * @brief 初值的准备
   * 
   * 也是统一的接口, 除了 Full multi grid 需要用一个算法获得初值, 其余
   * 方法的初值都直接由问题给定.
   * 
   * @param sol 
   * @param p_mesh 
   * @param p_force_handle 
   */
  void initial_value( sol_t& sol, const ptr_mesh_t& p_mesh, 
                      const ptr_force_t& p_force_handle ) {
    PROBLEM::initial_value( sol, *p_mesh );
  }

private:
  double TimeStepLength( const sol_t& sol, const mesh_t& mesh ) {
    const unsigned int n_sol = sol.size();
    double dt = 1000.;
#pragma omp parallel for reduction(min:dt)
    for( int i=0; i<(int)n_sol; ++i )
    {
      value_t C_M = Transform<value_t>::MaxRootOfHermitePolynomial( sol[i].GetOrder()+1 );
      // 位置空间只是一维, 所以这里才可以只取 Center 的第一个分量
      double lambda = fabs( sol[i].Center()[0] ) + C_M * sol[i].Scaling();
      dt = std::min( dt, CFL * mesh.Length(i) / lambda );
    }
    return dt;
  }
  
  /** 
   * @brief 时间向前发展一步时解的更新
   *
   * 显格式方程是 \f$ f^{n+1} = f^{n} + \Delta t r^{n} \f$, 其中 \f$
   * r^{n} \f$ 是当前解计算出来的方程残量 rhs-R(f).
   * 
   * 这里显格式的残量 @p res_sol 由外部计算传入, 而不是在此函数体中实现,
   * 是基于这样的考虑: 若每个时间步都检验是否到达稳态, 调用此函数可以去
   * 除残量的重复计算.
   *
   * @param sol 
   * @param res_sol 显格式解的残量. 以参数传入而不是在本函数中计算.
   * @param dt 
   */
  void ForwardStep( sol_t& sol, const sol_t& res_sol, const value_t dt ) {
    const unsigned int n_ele = sol.size();
    sol_t new_sol(n_ele);
#pragma omp parallel for
    for( int i=0; i<(int)n_ele; ++i )
    {
      new_sol[i].Reinit( sol[i] );
      sol[i].Add( dt, res_sol[i] );
      new_sol[i].ProjectToStdSpace( sol[i] );
    }
    std::swap( sol, new_sol );
  }
};  // class ExplicitTimeForward


/**
 * @class ExplicitRichardson
 * @brief 基于显格式时间推进的稳态求解器
 *
 * 成和某些文献习惯称之为 Jacobi 迭代, 但似乎文献中更正确的名称是
 * Richardson 迭代. Richardson 迭代是这样的格式
 *
 * \f[ u^{n+1} = u^{n} + w ( r-R(u^{n}) ) \f]
 *
 * 因而, 若取 \f$ w \f$ 为时间步长, 就是一个 ExplicitTimeForward 方法了.
 * 不过既然是稳态求解器, 可以考虑一些加速技术. 暂时想到的有: 
 *
 * - 使用 local time step.
 *
 * 与 ExplicitTimeForward 类一样, 这个不适合从 SingleGridMethod 类派生,
 * 但为了调用的方便, 仍然要提供相同的函数接口: ExactSolve(), Solve() 和
 * PostProcessing(). 不过可能可以考虑 ExplicitTimeForward,
 * ExplicitRichardson 和 ExplicitSGS 这三个基于时间推进的格式弄一个共同
 * 的基类出来.
 *
 * @see
 * - ExplicitTimeForward
 * - ExplicitSGS
 *
 * @note
 * 对于 Si 器件的测试算例 (Vbias=0.5V, N=1920, M=5, 100-400-100nm 器件)
 * 而言, local time step 和 global time step 比较, 此 local time step
 * 方法的迭代总步数只减少了 20% 左右, 总耗时也就少 30% 左右.
 */
template <typename _RESIDUAL>
// template <typename DISTRIBUTION, typename MESH, typename RESIDUAL, typename PROBLEM>
class ExplicitRichardson {
public:
  typedef _RESIDUAL RESIDUAL;  // 多重网格类又是不能直接取 RESIDUAL
  typedef typename RESIDUAL::mesh_t mesh_t;
  typedef typename RESIDUAL::sol_t sol_t;
  typedef typename RESIDUAL::force_t force_t;
  typedef typename RESIDUAL::PROBLEM PROBLEM;

  typedef typename sol_t::dis_t dis_t; 
  typedef typename dis_t::value_type value_t;
  typedef typename dis_t::Indices index_t;
  typedef typename dis_t::Velocity velocity_t;

  typedef std::shared_ptr<mesh_t> ptr_mesh_t;
  typedef std::shared_ptr<force_t> ptr_force_t;

private: 
  double CFL;
  double tolerance;             /**< 判断是否达到稳态 (收敛) 的残量忍量 */

public:
  ExplicitRichardson() : CFL(0.5), tolerance(1e-6) {}
  ExplicitRichardson( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) {
    read_config( p_mesh, cf );
  }
  ExplicitRichardson( const ptr_mesh_t& p_mesh, std::istream& is ) {
    is.read((char *)&CFL, sizeof(double));
    is.read((char *)&tolerance, sizeof(double));
  }

  void Reinit( const ptr_mesh_t& p_mesh ) { };

  void read_config( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) {
    CFL = (double)cf.Value( "Time", "CFL" );
    tolerance = (double)cf.Value( "Error", "Tol" );
  }

  void print_config() const {
    std::cout << "ExplicitRichardson: CFL = " << CFL << std::endl;
    std::cout << "ExplicitRichardson: Tol = " << tolerance << std::endl;
  };

  void dump_data( std::ostream& os ) const {
    os.write((const char *)&CFL, sizeof(double));
    os.write((const char *)&tolerance, sizeof(double));
  }


  void reset_tol( const double _tol ) { 
    tolerance = _tol;
  }
  double& reset_tol() {
    return tolerance;
  }
  const double& reset_tol() const {
    return tolerance;
  }

  /// 这个函数大家都一样, 是否也弄到一起呢?
  void ExactSolve( sol_t& sol, const sol_t& rhs, 
                   const ptr_mesh_t& p_mesh, 
                   const ptr_force_t& p_force_handle, 
                   const unsigned int max_steps=0 ) {
    double initial_res;
    RESIDUAL::GetResidualPL2Norm( initial_res, sol, rhs, *p_mesh, *p_force_handle );
    // 以下的设置是不是效率又防止最后残量降不下来.
    const double _Tol = (initial_res > 10*tolerance) ? tolerance : 0.1 * tolerance;

    double res;
    unsigned int step = 0;
    do {
      ++step;
      Solve( sol, rhs, p_mesh, p_force_handle, 1 );
//      RESIDUAL::GetResidualMaxNorm( res, sol, rhs, *p_mesh, *p_force_handle );
      RESIDUAL::GetResidualPL2Norm( res, sol, rhs, *p_mesh, *p_force_handle );
    } while (res >= _Tol && (step<max_steps || max_steps==0) );
//    std::cout<< "coarsest steps: = " << step << ", res = " << res << std::endl;
  }


  /** 
   * @brief 时间发展一步的函数接口
   *
   * @param sol 解, 传入初值, 返回新解. 暂时考虑传入返回的解在各单元上
   *            都是标准展开
   * @param rhs 方程右端项. 最密网格上是 0, 粗网格上则是细网格上残量的
   *            投影, 各单元上展开参数系与 @p sol 相同, 也即未必是标准
   *            展开. 由于迭代一次后, rhs 和 sol 的展开参数系就不同了,
   *            所以这个要求未必需要, 还是在方程求解时多做一步投影罢
   *            (rhs 在 sol 参数系重新展开)
   * @param p_mesh 求解的网格
   * @param p_force_handle 外力项处理句柄, 提供 @p mesh 上的外力等信息
   * @param steps 基本迭代法求解步数
   */
  void Solve( sol_t& sol, const sol_t& rhs, 
              const ptr_mesh_t& p_mesh, 
              const ptr_force_t& p_force_handle,
              const unsigned int steps = 1 ) {
    for( unsigned int i=0; i<steps; ++i )
    {
      // 残量计算
      sol_t res_sol( sol.size() );
      RESIDUAL::GetResidual( res_sol, sol, rhs, *p_mesh, *p_force_handle );
      // 解的更新
      ForwardStep( sol, res_sol, *p_mesh );

      p_force_handle->CalcExternalForce( sol );
    }
  };

  /** 
   * @brief 基本迭代法求解的后处理
   *
   * 使用局部时间步长加速以后, 算法就未必守恒, 故而这里该加一个后处理的.
   * 先留个标记.
   *
   * @param sol 解
   * @param p_mesh 网格
   */
  void PostProcessing( sol_t& sol, const ptr_mesh_t& p_mesh ) { 
    PROBLEM::PostProcessing( sol, *p_mesh );
  };  

  /** 
   * @brief 初值的准备
   * 
   * 也是统一的接口, 除了 Full multi grid 需要用一个算法获得初值, 其余
   * 方法的初值都直接由问题给定.
   * 
   * @param sol 
   * @param p_mesh 
   * @param p_force_handle 
   */
  void initial_value( sol_t& sol, const ptr_mesh_t& p_mesh, 
                      const ptr_force_t& p_force_handle ) {
    PROBLEM::initial_value( sol, *p_mesh );
  }

private:
  /** 
   * @brief local time step based on CFL
   *
   * 或许这里的局部时间步长还应该考虑相邻单元的特征速度情况.
   * 
   * @param _dis 
   * @param _grid_size 
   * 
   * @return 
   */
  double LocalTimeStepLength( const dis_t& _dis, const double& _grid_size ) {
    value_t C_M = Transform<value_t>::MaxRootOfHermitePolynomial( _dis.GetOrder()+1 );
    // 下式位置空间一维时才正确
    double lambda = fabs( _dis.Center()[0] ) + C_M * _dis.Scaling();
    return CFL * _grid_size / lambda;
  }
  
  /** 
   * @brief 向前发展一步
   *
   * 显格式方程是 \f$ f^{n+1} = f^{n} + \Delta t r^{n} \f$, 其中 \f$
   * r^{n} \f$ 是当前解计算出来的方程残量 rhs-R(f).
   * 
   * 这里显格式的残量 @p res_sol 由外部计算传入, 而不是在此函数体中实现,
   * 是基于这样的考虑: 若每个时间步都检验是否到达稳态, 调用此函数可以去
   * 除残量的重复计算.
   *
   * @param sol 
   * @param res_sol 显格式解的残量. 以参数传入而不是在本函数中计算.
   * @param mesh
   */
  void ForwardStep( sol_t& sol, const sol_t& res_sol, const mesh_t& mesh ) {
    const unsigned int n_ele = sol.size();
    sol_t new_sol(n_ele);
#pragma omp parallel for
    for( int i=0; i<(int)n_ele; ++i )
    {
      double local_dt = LocalTimeStepLength( sol[i], mesh.Length(i) );
      PositivityStepLength( local_dt, sol[i], res_sol[i] );
      
      new_sol[i].Reinit( sol[i] );
      sol[i].Add( local_dt, res_sol[i] );
      new_sol[i].ProjectToStdSpace( sol[i] );
    }
    std::swap( sol, new_sol );
  }

  void PositivityStepLength( value_t& local_dt, const dis_t &dis, const dis_t &res_dis ){
    _PreservingPositivityStepLength<dis_t, dis_t, PROBLEM>( local_dt, dis, res_dis );
  }
};  // class ExplicitRichardson


/**
 * @class ExplicitSGS
 * @brief 基于显格式时间推进的对称 Gauss-Seidel 稳态求解器
 *
 * 这实际就是 ExplicitRichardson 类的一个变体, 每个分量都采用最新的解来
 * 计算. 注意: 由于这是一个基于 Richardson 变体的 Gauss-Seidel 迭代, 因
 * 而没法做 SOR 了 (若在此迭代上做一个线性方程组迭代的 SOR 那种线性组合,
 * 无非就是在时间步长上乘了一个系数)
 *
 * 与 ExplicitTimeForward 类一样, 这个不适合从 SingleGridMethod 类派生,
 * 但为了调用的方便, 仍然要提供相同的函数接口: ExactSolve(), Solve() 和
 * PostProcessing(). 不过可能可以考虑 ExplicitTimeForward,
 * ExplicitRichardson 和 ExplicitSGS 这三个基于时间推进的格式弄一个共同
 * 的基类出来.
 *
 * @see
 * - ExplicitTimeForward
 * - ExplicitRichardson
 *
 * @note
 *
 * 从测试的算例来看, 这个算法 (CFL=0.5) 较 SLGSNewton 方法总迭代步数多
 * 四倍左右, 总耗时多了三分之一左右. 若取 CFL=0.75, 总步数多不到 2.5 倍,
 * 总耗时比 SLGSNewton 方法略少了. 这样的话, 用它来做光滑子, 未必比
 * SLGSNewton 总耗时会多了吧.
 *
 * 与 ExplicitRichardson 方法比较, 在相同 CFL 数下, CFL=0.5 时总迭代步
 * 数大概是 ExplicitRichardson 的三分之一多些; CFL=0.75 时, 总迭代步数
 * 大概是 ExplicitRichardson 的三分之一. 但是注意到此类的一步是对称的,
 * 相当于 ExplicitRichardson 中的两步, 所以实际仅仅是加速了一点点. 在
 * CFL=0.75 时总耗时仅是 ExplicitRichardson 的一半, 而这两个方式计算量
 * 实际应该是几乎一样的, 所以其他部分的程序也省了不少的计算量.
 */
template <typename _RESIDUAL>
// template <typename DISTRIBUTION, typename MESH, typename RESIDUAL, typename PROBLEM>
class ExplicitSGS {
public:
  typedef _RESIDUAL RESIDUAL;  // 多重网格类又是不能直接取 RESIDUAL
  typedef typename RESIDUAL::mesh_t mesh_t;
  typedef typename RESIDUAL::sol_t sol_t;
  typedef typename RESIDUAL::force_t force_t;
  typedef typename RESIDUAL::PROBLEM PROBLEM;

  typedef typename sol_t::dis_t dis_t; 
  typedef typename dis_t::value_type value_t;
  typedef typename dis_t::Indices index_t;
  typedef typename dis_t::Velocity velocity_t;

  typedef std::shared_ptr<mesh_t> ptr_mesh_t;
  typedef std::shared_ptr<force_t> ptr_force_t;

private: 
  double CFL;
  double tolerance;             /**< 判断是否达到稳态 (收敛) 的残量忍量 */

public:
  ExplicitSGS() : CFL(0.5), tolerance(1e-6) {}
  ExplicitSGS( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) {
    read_config( p_mesh, cf );
  }
  ExplicitSGS( const ptr_mesh_t& p_mesh, std::istream& is ) {
    is.read((char *)&CFL, sizeof(double));
    is.read((char *)&tolerance, sizeof(double));
  }

  void Reinit( const ptr_mesh_t& p_mesh ) { };

  void read_config( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) {
    CFL = (double)cf.Value( "Time", "CFL" );
    tolerance = (double)cf.Value( "Error", "Tol" );
  }

  void print_config() const {
    std::cout << "ExplicitSGS: CFL = " << CFL << std::endl;
    std::cout << "ExplicitSGS: Tol = " << tolerance << std::endl;
  };

  void dump_data( std::ostream& os ) const {
    os.write((const char *)&CFL, sizeof(double));
    os.write((const char *)&tolerance, sizeof(double));
  }

  void reset_tol( const double _tol ) { 
    tolerance = _tol;
  }
  double& reset_tol() {
    return tolerance;
  }
  const double& reset_tol() const {
    return tolerance;
  }

  /// 这个函数大家都一样, 是否也弄到一起呢?
  void ExactSolve( sol_t& sol, const sol_t& rhs, 
                   const ptr_mesh_t& p_mesh, 
                   const ptr_force_t& p_force_handle, 
                   const unsigned int max_steps=0 ) {
    double initial_res;
    RESIDUAL::GetResidualPL2Norm( initial_res, sol, rhs, *p_mesh, *p_force_handle );
    // 以下的设置是不是效率又防止最后残量降不下来.
    const double _Tol = (initial_res > 10*tolerance) ? tolerance : 0.1 * tolerance;

    double res;
    unsigned int step = 0;
    do {
      ++step;
      Solve( sol, rhs, p_mesh, p_force_handle, 1 );
//      RESIDUAL::GetResidualMaxNorm( res, sol, rhs, *p_mesh, *p_force_handle );
      RESIDUAL::GetResidualPL2Norm( res, sol, rhs, *p_mesh, *p_force_handle );
    } while (res >= _Tol && (step<max_steps || max_steps==0) );
//    std::cout<< "coarsest steps: = " << step << ", res = " << res << std::endl;
  }

  void Solve( sol_t& sol, const sol_t& rhs, 
              const ptr_mesh_t& p_mesh, 
              const ptr_force_t& p_force_handle,
              const unsigned int steps = 1 ) {
    for( unsigned int i=0; i<steps; ++i ){
      SymmetricGaussSeidel( sol, rhs, *p_mesh, *p_force_handle );
      p_force_handle->CalcExternalForce( sol );
    }
  };

  void PostProcessing( sol_t& sol, const ptr_mesh_t& p_mesh ) { 
    PROBLEM::PostProcessing( sol, *p_mesh );
  };  

  /** 
   * @brief 初值的准备
   * 
   * 也是统一的接口, 除了 Full multi grid 需要用一个算法获得初值, 其余
   * 方法的初值都直接由问题给定.
   * 
   * @param sol 
   * @param p_mesh 
   * @param p_force_handle 
   */
  void initial_value( sol_t& sol, const ptr_mesh_t& p_mesh, 
                      const ptr_force_t& p_force_handle ) {
    PROBLEM::initial_value( sol, *p_mesh );
  }

private:
  void SymmetricGaussSeidel( sol_t& sol, const sol_t& rhs, 
                             const mesh_t& mesh, 
                             force_t& force_handle ) {
    GaussSeidelOneStep( sol, rhs, mesh, force_handle );
    ReverseGaussSeidelOneStep( sol, rhs, mesh, force_handle );
  }
  
  void GaussSeidelOneStep( sol_t& sol, const sol_t& rhs,
                           const mesh_t& mesh, 
                           force_t& force_handle ) {
    for( size_t k=0; k<sol.size(); ++k ){
      GaussSeidelOneStep( sol, rhs, mesh, force_handle, k );
      force_handle.CalcExternalForce( sol, k );
    }      
  }

  void ReverseGaussSeidelOneStep( sol_t& sol, const sol_t& rhs, 
                           const mesh_t& mesh,
                           force_t& force_handle ) {
    for( size_t k=sol.size()-1; k != 0; --k ){
      GaussSeidelOneStep( sol, rhs, mesh, force_handle, k );
      force_handle.CalcReverseExternalForce( sol, k );
    }      
  }

  /** 
   * @brief Gauss-Seidel 迭代某个单元上的求解
   * 
   * @param k 单元指标
   */
  void GaussSeidelOneStep( sol_t& sol, const sol_t& rhs,
                   const mesh_t& mesh,
                   const force_t& force_handle,
                   const unsigned int k ){
    dis_t &dis = sol[k];
    
    dis_t res_dis( dis.Center(), dis.Scaling(), dis.GetOrder() );
    Project( rhs[k], res_dis );

    /// 切记: res_dis = rhs - R(f)
    RESIDUAL::GetResidual( res_dis, sol, mesh, force_handle, k );
    
    double local_dt = LocalTimeStepLength( dis, mesh.Length(k) );
    PositivityStepLength( local_dt, dis, res_dis );

    dis_t new_dis( dis.Center(), dis.Scaling(), dis.GetOrder() );
    dis.Add( local_dt, res_dis );
    new_dis.ProjectToStdSpace( dis );
    std::swap( dis, new_dis );
  }

  double LocalTimeStepLength( const dis_t& _dis, const double& _grid_size ) {
    value_t C_M = Transform<value_t>::MaxRootOfHermitePolynomial( _dis.GetOrder()+1 );
    // 下式位置空间一维时才正确
    double lambda = fabs( _dis.Center()[0] ) + C_M * _dis.Scaling();
    return CFL * _grid_size / lambda;
  }

  void PositivityStepLength( value_t& local_dt, const dis_t &dis, const dis_t &res_dis ){
    _PreservingPositivityStepLength<dis_t, dis_t, PROBLEM>( local_dt, dis, res_dis );
  }
};  // class ExplicitSGS


/**
 * @class HeunTimeForward
 * @brief 时间推进的高级格式: Heun's method
 *
 * @see 
 * http://en.wikipedia.org/wiki/Heun%27s_method
 * 
 * @see
 * - ExplicitTimeForward
 * - ExplicitRichardson
 * - ExplicitSGS
 *
 * 真该给他们写个统一的基类, 重复代码太多了. 
 */
template <typename _RESIDUAL>
// template <typename DISTRIBUTION, typename MESH, typename RESIDUAL, typename PROBLEM>
class HeunTimeForward {
public:
  typedef _RESIDUAL RESIDUAL;  // 多重网格类又是不能直接取 RESIDUAL
  typedef typename RESIDUAL::mesh_t mesh_t;
  typedef typename RESIDUAL::sol_t sol_t;
  typedef typename RESIDUAL::force_t force_t;
  typedef typename RESIDUAL::PROBLEM PROBLEM;

  typedef typename sol_t::dis_t dis_t; 
  typedef typename dis_t::value_type value_t;
  typedef typename dis_t::Indices index_t;
  typedef typename dis_t::Velocity velocity_t;

  typedef std::shared_ptr<mesh_t> ptr_mesh_t;
  typedef std::shared_ptr<force_t> ptr_force_t;

private: 
  double CFL;
  double tolerance;             /**< 判断是否达到稳态 (收敛) 的残量忍量 */

public:
  HeunTimeForward() : CFL(0.5), tolerance(1e-6) {}
  HeunTimeForward( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) {
    read_config( p_mesh, cf );
  }
  HeunTimeForward( const ptr_mesh_t& p_mesh, std::istream& is ) {
    is.read((char *)&CFL, sizeof(double));
    is.read((char *)&tolerance, sizeof(double));
  }

  void Reinit( const ptr_mesh_t& p_mesh ) { };

  void read_config( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) {
    CFL = (double)cf.Value( "Time", "CFL" );
    tolerance = (double)cf.Value( "Error", "Tol" );
  }

  void print_config() const {
    std::cout << "ExplicitTimeForward: CFL = " << CFL << std::endl;
    std::cout << "ExplicitTimeForward: Tol = " << tolerance << std::endl;
  };

  void dump_data( std::ostream& os ) const {
    os.write((const char *)&CFL, sizeof(double));
    os.write((const char *)&tolerance, sizeof(double));
  }


  void reset_tol( const double _tol ) { 
    tolerance = _tol;
  }
  double& reset_tol() {
    return tolerance;
  }
  const double& reset_tol() const {
    return tolerance;
  }

  void ExactSolve( sol_t& sol, const sol_t& rhs, 
                   const ptr_mesh_t& p_mesh, 
                   const ptr_force_t& p_force_handle, 
                   const unsigned int max_steps=0 ) {
    double initial_res;
    RESIDUAL::GetResidualPL2Norm( initial_res, sol, rhs, *p_mesh, *p_force_handle );
    // 以下的设置是不是效率又防止最后残量降不下来.
    const double _Tol = (initial_res > 10*tolerance) ? tolerance : 0.1 * tolerance;

    double res;
    unsigned int step = 0;
    do {
      ++step;
      Solve( sol, rhs, p_mesh, p_force_handle, 1 );
//      RESIDUAL::GetResidualMaxNorm( res, sol, rhs, *p_mesh, *p_force_handle );
      RESIDUAL::GetResidualPL2Norm( res, sol, rhs, *p_mesh, *p_force_handle );
    } while (res >= _Tol && (step<max_steps || max_steps==0) );
//    std::cout<< "coarsest steps: = " << step << ", res = " << res << std::endl;
  }


  /** 
   * @brief 时间发展一步的函数接口
   *
   * 
   *
   * @param sol 解, 传入初值, 返回新解. 暂时考虑传入返回的解在各单元上
   *            都是标准展开
   * @param rhs 方程右端项. 最密网格上是 0, 粗网格上则是细网格上残量的
   *            投影, 各单元上展开参数系与 @p sol 相同, 也即未必是标准
   *            展开. 由于迭代一次后, rhs 和 sol 的展开参数系就不同了,
   *            所以这个要求未必需要, 还是在方程求解时多做一步投影罢
   *            (rhs 在 sol 参数系重新展开)
   * @param p_mesh 求解的网格
   * @param p_force_handle 外力项处理句柄, 提供 @p mesh 上的外力等信息
   * @param steps 基本迭代法求解步数
   */
  void Solve( sol_t& sol, const sol_t& rhs, 
              const ptr_mesh_t& p_mesh, 
              const ptr_force_t& p_force_handle,
              const unsigned int steps = 1 ) {
//    static double t = 0;
    for( unsigned int i=0; i<steps; ++i )
    {
      // 时间步长
      double dt = TimeStepLength( sol, *p_mesh );
      //    t += dt;
//      std::cout << "t = " << t << ", dt = " << dt << std::endl;
      // 残量计算
      sol_t res_sol( sol.size() );
      RESIDUAL::GetResidual( res_sol, sol, rhs, *p_mesh, *p_force_handle );
      // 解的更新
      sol_t updated_sol=sol;
      ForwardStep( updated_sol, res_sol, dt );

      sol_t new_res_sol( sol.size() );
      RESIDUAL::GetResidual( new_res_sol, updated_sol, rhs, *p_mesh, *p_force_handle );

      ForwardStep( sol, res_sol, new_res_sol, dt );

      p_force_handle->CalcExternalForce( sol );
    }
  };

  /** 
   * @brief 基本迭代法求解的后处理
   * 
   * 显格式时间推进是守恒型的, 因而这里不需处理
   *
   * @param sol 解
   * @param p_mesh 网格
   */
  void PostProcessing( sol_t& sol, const ptr_mesh_t& p_mesh ) {
    // 如果只当时间发展来用而不放到 multigrid 中去, 这个要去掉, 不过不
    // 去掉应该也没关系, 本来质量就保守恒么
    PROBLEM::PostProcessing( sol, *p_mesh );
  };  

  /** 
   * @brief 初值的准备
   * 
   * 也是统一的接口, 除了 Full multi grid 需要用一个算法获得初值, 其余
   * 方法的初值都直接由问题给定.
   * 
   * @param sol 
   * @param p_mesh 
   * @param p_force_handle 
   */
  void initial_value( sol_t& sol, const ptr_mesh_t& p_mesh, 
                      const ptr_force_t& p_force_handle ) {
    PROBLEM::initial_value( sol, *p_mesh );
  }

private:
  double TimeStepLength( const sol_t& sol, const mesh_t& mesh ) {
    const unsigned int n_sol = sol.size();
    double dt = 1000.;
    for( unsigned int i=0; i<n_sol; ++i )
    {
      value_t C_M = Transform<value_t>::MaxRootOfHermitePolynomial( sol[i].GetOrder()+1 );
      // 位置空间只是一维, 所以这里才可以只取 Center 的第一个分量
      double lambda = fabs( sol[i].Center()[0] ) + C_M * sol[i].Scaling();
      dt = std::min( dt, CFL * mesh.Length(i) / lambda );
    }
    return dt;
  }
  
  /** 
   * @brief 时间向前发展一步时解的更新
   *
   * 显格式方程是 \f$ f^{n+1} = f^{n} + \Delta t r^{n} \f$, 其中 \f$
   * r^{n} \f$ 是当前解计算出来的方程残量 rhs-R(f).
   * 
   * 这里显格式的残量 @p res_sol 由外部计算传入, 而不是在此函数体中实现,
   * 是基于这样的考虑: 若每个时间步都检验是否到达稳态, 调用此函数可以去
   * 除残量的重复计算.
   *
   * @param sol 
   * @param res_sol 显格式解的残量. 以参数传入而不是在本函数中计算.
   * @param dt 
   */
  void ForwardStep( sol_t& sol, const sol_t& res_sol, const value_t dt ) {
    const unsigned int n_ele = sol.size();
    sol_t new_sol(n_ele);
    for( unsigned int i=0; i<n_ele; ++i )
    {
      new_sol[i].Reinit( sol[i] );
      sol[i].Add( dt, res_sol[i] );
      new_sol[i].ProjectToStdSpace( sol[i] );
    }
    std::swap( sol, new_sol );
  }

  void ForwardStep( sol_t& sol, const sol_t& res_sol, 
                    const sol_t& new_res_sol, const value_t dt ) {
    const unsigned int n_ele = sol.size();
    sol_t new_sol(n_ele);
    for( unsigned int i=0; i<n_ele; ++i )
    {
      sol[i].Add( 0.5*dt, res_sol[i] );
      dis_t tmp( sol[i].Center(), sol[i].Scaling(), sol[i].GetOrder() );
      Project( new_res_sol[i], tmp );
      sol[i].Add( 0.5*dt, tmp );
      new_sol[i].Reinit( sol[i] );
      new_sol[i].ProjectToStdSpace( sol[i] );
    }
    std::swap( sol, new_sol );
  }
};  // class HeunTimeForward

/**
 * 本文件公共函数定义
 * 
 */

/** 
 * @brief 控制更新步长以使更新后分布函数密度和温度保正
 *
 * 若来自 SLGSNewton 方法, @p _df 则是一个 dealii 的向量, 而若来自
 * ExplicitTime 相关方法, 则还是一个分布函数, 所以这里用了两个模板参数.
 * 最后一个则是问题相关参数的模板类 (应该提供一个温度的下限
 * Temperature)
 * 
 * @param mu_max 更新步长, 传进来时候应该已经有一个缺省的步长
 * @param dis 分布函数
 * @param _df 更新的增量, \f$ \Delta f \f$
 *
 * 以下的解释来自于原来 SLGSNewton 方法测试时关于公式和测试记录的说明.
 *
 * 新的解为 \f$ f + \mu \Delta f \f$.
 *
 * 首先, 若 \f$ f \f$ 在当前位置的速度 \f$ u \f$ 和温度 \f$ \theta
 * \f$ 上展开, 必有 \f$ f_1=0,~f_2=0 \f$. 当 \f$ f \f$ 在其它的速度
 * \f$ u' \f$ 和温度 \f$ \theta' \f$ 上展开时, \f$ f_1, f_2 \f$ 不为
 * 0. 它们间的关系为 \f{align*}{ \rho (u-u') = f_1', \\ \rho \vert
 * u-u' \vert^2 + \rho (\theta - \theta') = 2 f_2', \f} 于是有 \f[
 * \theta = \frac{ 2 f_2' - (f_1')^2 / \rho }{\rho} + \theta'. \f] 现
 * 在已知第 m 步的分布 \f$ f^{(m)} \f$ (所以这里的循环是不是也多几步
 * 会好一些?) 是在标准展开的, 速度和温度分别记为\f$ u^{(m)},
 * \theta^{(m)} \f$. 然后根据这里的计算, 再次标准投影前, 有 \f[ \rho
 * = \rho^{(m)} + \mu \Delta f_0^{(m)}, u' = u^{(m)} , \quad \theta'
 * = \theta^{(m)}, \quad f_1' = \mu \Delta f_1^{(m)}, \quad f_2' =
 * \mu \Delta f_2^{(m)}. \f] 可得\f[ \theta^{(m+1)} = \frac{2 \mu
 * \Delta f_2^{(m)} - \mu^2 (\Delta f_1^{(m)})^2 / (\rho^{(m)} + \mu
 * \Delta f_0^{(m)})}{\rho^{(m)} + \mu \Delta f_0^{(m)}} +
 * \theta^{(m)}. \f] 计算过程要求 \f$ \theta > 0 \f$, 即有 \f[ \mu^2
 * (\Delta f_1^{(m)})^2 - 2 \mu (\rho^{(m)} + \mu \Delta f_0^{(m)})
 * \Delta f_2^{(m)} - (\rho^{(m)} + \mu \Delta f_0^{(m)})^2
 * \theta^{(m)} < 0. \f] 化简得 \f[ A \mu^2 + B \mu + C <0, \f] 其中,
 * \f[ A = (\Delta f_1^{(m)})^2- 2 \Delta f_0^{(m)} \Delta f_2^{(m)}
 * - \theta^{(m)} (\Delta f_0^{(m)})^2, \quad B = -2 \rho^{(m)}
 * (\Delta f_2^{(m)} + \Delta f_0^{(m)} \theta^{(m)}), \quad C =
 * -(\rho^{(m)})^2 \theta^{(m)}.\f] 讨论其解要分为三种情形:
 *
 * - \f$ A=0 \f$. 这种情形程序里一般不大可能出现, 此时退化为线性问
 *    题. 因为 \f$ C \f$ 必小于 0, 若 \f$ B \leq 0 \f$, 显然对所有
 *    \f$ \mu >0 \f$ 不等式恒成立. 若 \f$ B>0 \f$, 则有 \f$ 0 <
 *    \mu < - \frac{C}{B} \f$.
 * - \f$ A>0 \f$. 此时, 还需要考虑 \f$ \Delta = B^2 - 4 AC \f$.
 *  - 若 \f$ \Delta >0 \f$, 有 \f[ \frac{-B - \sqrt{B^2-4AC}}{2A}
 *    < \mu < \frac{-B + \sqrt{B^2-4AC}}{2A}, \f] 当然还要结合上
 *    \f$ \mu >0 \f$.
 *  - 若 \f$ \Delta \leq 0 \f$, 此时无解, 应该打印出一个错误信息! 当
 *    然因为 \f$ \theta^{(m)} >0 \f$, 所以这种情形根本就不应该出现.
 * - \f$ A<0 \f$.  
 *  - 若 \f$ \Delta \geq 0 \f$, 有 \f[ \mu > \mu_2 = \frac{-B -
 *    \sqrt{ B^2-4AC}}{2A} \quad \text{or} \quad \mu < \mu_1 =
 *    \frac{-B + \sqrt{B^2-4AC}}{2A}. \f] 因为 \f$ C <0 \f$, 有
 *    \f$ \sqrt{\Delta} < |B| \f$. 于是可得,
 *   - 若 \f$ B>0 \f$, \f$ \mu_2 > \mu_1 >0 \f$, 这种情形我们就取 \f$
 *     0< \mu < \mu_1 \f$.
 *   - 若 \f$ B<0 \f$, \f$ 0 > \mu_2 > \mu_1 \f$, 这种情形可取所有的
 *     \f$ \mu >0 \f$.
 *  - 若 \f$ \Delta < 0 \f$, 对所有 \f$ \mu>0 \f$ 恒成立.
 * 
 * 从上面的温度保证可以得到参数可选的最大值 \f$ \mu_{\max} \f$, 于
 * 是最后拿来更新的参数选为 \f$ \mu = \min{\{\beta \mu_{\max},
 * \tilde{\mu}\}} \f$. 其中, 参数 \f$ \beta \f$ 暂时取为 0.5. 至于
 * 步长 \f$ \tilde{\mu} \f$, 在无松弛的 newton 迭代 (也就是 jacobi
 * 矩阵 \f$ \Lambda \f$) 时, 步长显然应该选 1.0. 而对于现在带松弛的
 * 情形, 选 1.0 的话, 几次迭代以后就会出现负温度 (这个原因后面解释),
 * 而 0.618 却可以比较稳定的算到最后. 但是这个 0.618 显然是很无理的!
 * 另一方面, 随着残量的减小, 我们的松弛逐渐减弱, 因而步长 \f$
 * \tilde{\mu} \f$ 可以逐渐趋近于 1.0. 因而现在使用公式 \f[
 * \tilde{\mu} = \frac{\rho} {\rho+ \lambda} \f] 计算步长. 需要说明
 * 地是, 对这个公式成并不能给出合理的解释! 按照无松弛步长为 1.0 的
 * 原则, 成形式上推出来的公式是这个公式的倒数. 成对大于 1.0 的数没
 * 有太多信心 (用这个数的话, 大概就会上溢出? 从而出现了负温度), 所
 * 以就使用了其倒数.  不过对于正在测试的例子, 这个公式确实加速了不
 * 少! 原本三十多分钟的程序快了近十分钟, 总的迭代步数也少了许多.
 *
 * @note 
 * - 类似地, 对于不分裂的时间发展格式, 也可以通过减小时间步长来达
 * 到避免出现 nan.
 * - 还有, 这个求根公式可能会出现 \f$ 4 AC << B^2 \f$ 而导致大
 * 数吃小数的问题, 至少 octave 中就没掉了, 尚不知该如何解决.
 * 
 * 现在回到步长取为 1.0 时还是会出现负温度的情形. 经过仔细检查和打
 * 出数据比较, 现在的防止出现负温度计算公式本身是没有问题的! 原因就
 * 是取 \f$ \tilde{\mu}=1.0 \f$ 时, 由于迭代之初变化剧烈而触碰到了
 * 机器精度的概念. 对于 double 型, 其机器精度是 2.2e-16, 比较出现
 * NaN 时成的温度计算值和 NRxx 中的计算值, 两者并不相等! 比如前者
 * 2.4e-17 是一个正数, 后者 -7.04731e-19 则是负数, 而当温度远远大于
 * 这个量级时, 两者的计算结果却是一致的! 另外由于晶格温度是 1e-2 量
 * 级的数, 所以这里不一致的原因极有可能是机器精度的扰动所致.
 *
 * 按说, 由于 \f$ \mu \leq \beta \mu_{\max} \f$, 应该是不会出现这么
 * 小的温度值的. 而且到了这里, 似乎也没有好的策略绝对的防止负温度的
 * 出现了...
 *
 * 那么回到要求温度大于 0 这里, 我们加强这个条件, 把它修改为 \f[
 * \theta > \tilde{\theta}, \f] 其中 \f$ \tilde{\theta} >0 \f$ 是一个
 * 事先给定的常数, 根据目前得到的这些稳定解来看, 我们可以选择\f$
 * \tilde{\theta} = 0.1 \theta_0 \f$, 而 \f$ \theta_0 \f$ 是晶格温度.
 * 于是, 选择步长的一元二次方程变为 \f[ \mu^2 (\Delta f_1^{(m)})^2 -
 * 2 \mu (\rho^{(m)} + \mu \Delta f_0^{(m)}) \Delta f_2^{(m)} -
 * (\rho^{(m)} + \mu \Delta f_0^{(m)})^2 (\theta^{(m)} -
 * \tilde{\theta}) < 0. \f] 即有 \f[ A = (\Delta f_1^{(m)})^2- 2
 * \Delta f_0^{(m)} \Delta f_2^{(m)} - (\theta^{(m)}-\tilde{\theta})
 * (\Delta f_0^{(m)})^2, \quad B = - 2 \rho^{(m)} (\Delta f_2^{(m)}
 * + \Delta f_0^{(m)} ( \theta^{(m)}-\tilde{\theta})), \quad C =
 * -(\rho^{(m)})^2 (\theta^{(m)}-\tilde{\theta}).\f] 如此应该是不会碰
 * 到机器精度的问题了.
 *
 * @note 更细致地, 可能还要判断修改后的密度仍然大于 0, 这个一般没有
 * 问题的吧.
 *
 * 即使这个时候, \f$ \tilde{\mu} =1.0 \f$ 仍然不是一个合适的选择,
 * 发现这样的话会到其它数值解上去? 残量基本不变, osci 在 5.67 左右
 * 也变化很小! 将这个图画出来, 大概还是在原先可能出现 NaN 的那几个
 * 单元, 解的性质极差!
 *
 * @note 多维的时候, 只要把 \f$ (f_1')^2 \f$ 和 \f$ f_2' \f$ 分别换成
 * \f$ \sum_d (f_{e_d}')^2 \f$ 和 \f$ \sum_d f_{e_{2d}}' \f$. 还有, 不
 * 要忘了, 原始等式中是 \f$ D \rho (\theta - \theta') \f$, 也就是 \f$
 * \theta \f$ 前面要有个系数 \f$ D \f$.
 *
 * 虽然用扰动变量 \f$ [\rho, u, s, f_3, f_4,\ldots,f_M] \f$ 的方式已经
 * 弃用, 但这个记录仍先留着. 这里 \f$ s=\sqrt{\theta} \f$, 温度作为
 * \f$ s \f$ 的平方, 显然恒正. 但若步长不加选择, 对某些情形, 比如迭代刚
 * 开始的时候, 或是矩数很高的时候, 对步长不加选择会出现 \f$ s < 0 \f$,
 * 甚或 \f$ \rho < 0 \f$, 虽然程序并不会因此马上出错停止计算, 但这往往
 * 最终会导致 nan 的出现. 因而考虑到这些, 这里还是需要做一个小的修正.
 */
template <typename DISTRIBUTION, typename _VEC, typename PROBLEM>
void _PreservingPositivityStepLength( double &mu_max, 
                                      const DISTRIBUTION &dis, 
                                      const _VEC &_df ) {
  const unsigned int dim = DISTRIBUTION::dim;
  // 密度保正: \f$ \rho + \mu \Delta f_0 > 0 \f$
  const double _drho = *(_df.begin());
  if( _drho < 0 )
    mu_max = std::min( mu_max, - 0.8 * dis.front() / _drho );

  // 下面则是温度保正需要

  // 这个参数可能还是外部传进来最好 
  const static double _factor = 0.5;
  // 下式则限制了温度要始终大于等于 0.5 * Temperature, 所以注意类中
  // Temperature 的选择, 要使得 _theta > 0
  const double _theta = dim * ( dis.Scaling() * dis.Scaling() - 0.5 * PROBLEM::Temperature );
  double val_moments_1 = SumSquareOfFirstMoments( dim, _df.begin());
  /// \f$ \sum f_{e_{2d}} \f$
  double val_moments_2 = SumDiagOfSecondMoments<dim, _VEC>( _df );

  // 一元二次方程的系数 A, B, C
  double A = val_moments_1 - 2 * _drho * val_moments_2 - _theta * _drho * _drho;
  double B = - 2 * dis.front() * (val_moments_2 + _drho * _theta);
  double C = - dis.front() * dis.front() * _theta;
  // 一元二次方程判别式, 由于大数吃小数的缘故, 这个数可能就是 B*B.
  double Delta = B*B - 4 * A * C;
  
  if( A > dis.front() * 1.0e-15 )
  {// double 型的机器精度就是 2.22045e-16, 怎么程序还有这么多 -20, -30
    // 次的?
    if( Delta > 0 )
      if( B<=0 )
        mu_max = std::min( mu_max, _factor * (-B + sqrt(Delta)) / A * 0.5 );
      else/// 此时分子加法有可能就出现大数吃小数
        mu_max = std::min( mu_max, _factor * (-2 * C) / (B + sqrt(Delta)) );
    else 
    {// 此时意味着 0.5*Temperature 可能取大了, 实际温度要更小
	      std::cerr << "No admissible positive update step; theta = "
	                << _theta << std::endl;
	      mu_max = 0;
    }
  }
  else if( A < - dis.front() * 1.0e-15 )
  {
    if( Delta >= 0 && B > 0 )
      mu_max = std::min( mu_max, _factor * (-2 * C) / (B + sqrt(Delta)) );
  }
  else
  {/// 此时认为退化为线性, 还是打出一个信息, 看 1.0e-16 是否足够小
    // std::cerr << "A = " << A << ", using linear solution!" << std::endl;
    if( B > 0 )
      mu_max = std::min( mu_max, _factor * (-C) / B );
  }  
}

/**
 * 公共函数定义结束
 * 
 */


#endif // __SINGLE_GRID_METHOD_H__

/**
 * end of file
 * 
 */
