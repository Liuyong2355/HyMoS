/**
 * @file   SteadyState.h
 * @author HU Zhicheng <huzhicheng1986@gmail.com>
 * @date   Mon Nov 12 10:06:22 2012
 * 
 * @brief  稳态 Boltzmann 方程矩方程组迭代解法类
 *
 * 重新设计此类, 试图作为各种稳态问题的基本框架, 包括气体动理学,
 * Wigner-Poisson 方程等. 
 *
 * 关于方法, 目前采用的实际上是称作 Gummel 的迭代方法, 即 Boltzmann 方
 * 程和外力方程 (通常是 Poisson 方程) 解耦求解. 此类总假设外力的存在,
 * 但是具体的实现需要有 typename PROBLEM 类型提供.
 *
 * @note 这里只处理空间一维的情形, 速度高维应该不成问题.
 */

#ifndef __STEADY_STATE_H__
#define __STEADY_STATE_H__

#include <atomic>
#include <memory>
#include <vector>
#include <functional>
#include <sstream>
#include <omp.h>

// time_t 给出的时间精度为 s, 为了获得更好的精度, 需要用以下头文件中定
// 义的 timeval 类以及 gettimeofday(&timeval, NULL) 函数.
#include <sys/time.h> 

#include <boost/thread.hpp>
#include <boost/archive/text_oarchive.hpp>
#include <boost/serialization/vector.hpp>


#include <NRxx/Distribution.h>
#include <NRxx/Transform.h>
#include <NRxx/SpectralForODE.h>
#include <NRxx/NumericalFlux.h>

#include "config.h"
#include "ConfigFile.h"
#include "OperatorInterpolate.h"
#include "tools/tools.h"

/**
 * @class SteadyState
 * @brief 稳态 Boltzmann - 外力耦合问题求解器
 *
 * 把各部分程序拆开之后, 这个类实际仅剩下调用求解器运行程序的最外层循环
 * 以及一些输入输出函数.
 *
 * 以下是重新设定的模板参数类型:
 * - MESH: 单层网格类
 * - SOLUTION: Boltzmann 方程解类型 (现在是 DISTRIBUTION 的一个 vector)
 * - FORCE: 外力类, 这里需要其中外力方程求解器
 * - SOLVER: Boltzmann 方程求解器
 * - _BV: 提供初边值条件的具体问题类, 还包括特殊问题的参数等信息
 * 
 * 事实上, 模板类型 SOLVER 的构造已经包含了其它几个模板类型. 但毕竟
 * SOLVER 只是设定为 Boltzmann 方程的求解器. 而这里封装的是 Boltzmann
 * -- 外力耦合的问题, 从 SOLVER 取外力类型 FORCE 总觉得怪怪的, 因而这里
 * 做为独立模板类型传入. 至于其它几个模板参数, 也是类似的道理.
 * 
 */
template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
class SteadyState{
public:
//@{ 类型别名
  typedef MESH mesh_t;
  typedef SOLUTION sol_t;
  typedef FORCE force_t;
  typedef SOLVER solver_t;
  typedef _BV PROBLEM;

  typedef typename sol_t::dis_t dis_t;
  typedef typename dis_t::value_type value_t;
  typedef typename dis_t::Indices index_t;
  typedef typename dis_t::Velocity velocity_t;
  // 或许这个也放到 SteadyState 的模板中?
  typedef typename SOLVER::RESIDUAL RESIDUAL; 

  typedef std::shared_ptr<mesh_t> ptr_mesh_t;
  typedef std::shared_ptr<force_t> ptr_force_t;

  enum { dim = dis_t::dim };
//@}

  /**
   * 是否结束 run 函数的运行的旗标。当 run 函数在运行的时候，其值为
   * false，如果其值被设置成为 true，则 run 函数会结束运行。
   */
  bool is_stop_run; 
  bool is_exit ;                /**< 与 is_stop_run 都用于交互控制程序,
                                 * 用于退出整个程序 */
  bool is_read; /**< 设定为是否已读参数或初值 */

/// @name 构造函数与析构函数  
//@{ 
  /** 
   * @brief 缺省构造函数
   * 
   * 由于我们把参数的读取放入各类中, 所以各个类的指针要事先构造出来. 当
   * 然 read_config() 函数的实现避免了这个限制. 这要求各个类提供一个使
   * 用 ConfigFile 的构造函数.
   */
  SteadyState( ) : is_stop_run(true), 
                   is_exit(false), is_read(false),
                   n_thread(1),
                   p_solver(new solver_t()), 
                   p_force_handle(new force_t()), 
                   p_mesh( new mesh_t() ){};
//@}

private:
/// @name 参数
//@{
  int ORDER;                    /**< 使用的矩方法的阶数 */
  
  int n_thread;					 /**< 并行计算线程数 ==2025.06.17(星期二)== */
  value_t tol;                  /**< 求解终止条件. */
//@}

/// @name 变量
//@{
  // std::shared_ptr<PROBLEM> p_problem; /**< 或许 PROBLEM 中的信息也用指针来调用 */
  std::unique_ptr<solver_t> p_solver;
  ptr_force_t p_force_handle;
  ptr_mesh_t p_mesh;               /**< 内部网格, 顺序存储网格点坐标. */

  std::unique_ptr<sol_t> sol;    /**< 存储解变量的空间 */
  std::unique_ptr<sol_t> rhs;    /**< 方程右端项, 最密网格上是 0 */

  double _res;                  /**< Boltzmann 方程残量最大模, 当前仅用
                                 * _res 是否小于 tol 判断程序收敛 */
//@}

public:
//@{
  /** 
   * @brief 初值函数
   *
   * 初值总由问题提供, 但是对于 full multigrid 等希望得到一个好的初值的
   * 方法, 除了从问题获取初值, 还有一个算法来准备初值. 这两个方便都在这
   * 个函数接口内完成, 当然为了使用后者, 就需要所有的 solver 提供一个
   * initial_value() 相关的接口.
   */
  void initial_value() {
    if( sol.get()==NULL ) alloc_memory();

    p_solver->initial_value( *sol, p_mesh, p_force_handle );
    rhs->Reinit( *sol );
  }

  /// 运行主要接口
  void run();
  void run_in_memory(unsigned int max_steps, std::vector<double>& residuals,
                     const std::atomic<bool>* stop_requested = NULL,
                     const std::function<void(unsigned int, double)>& progress =
                         std::function<void(unsigned int, double)>()) {
    run_impl(false, max_steps, &residuals, stop_requested, progress);
  }
  const sol_t& solution_view() const { return *sol; }
  const mesh_t& mesh_view() const { return *p_mesh; }

private:
  void run_impl(bool legacy_io, unsigned int max_steps, std::vector<double>* residuals,
                const std::atomic<bool>* stop_requested = NULL,
                const std::function<void(unsigned int, double)>& progress =
                    std::function<void(unsigned int, double)>());
public:

  void set_maxwellian() {
    if( sol.get()==NULL ) alloc_memory();

    std::cerr << "Set the solution as Maxwellian ... " << std::endl;
    /// 目前只考虑 Maxwellian 的初值吧
    PROBLEM::initial_value( *sol, *p_mesh );

    /// 右端项是 0
    rhs->Reinit( *sol );

    std::cerr << "Set initial value OK! " << std::endl;
  }
//@}

/// @name IO 有关的函数
//@{
  /** 
   * @brief 从配置参数文件中读入参数。
   * 
   * @param filename 配置参数文件的文件名
   */
  void read_config(const std::string& filename);
  void read_config(ConfigFile& cf);
  /*
   *为了方便调试, 把当前的可调数据调出来看一下
   */
  void print_config() const ;
  /**
   * 将运行环境参数写入到流中和从流中载入运行环境参数。
   */
  void dump_config(std::ostream& os) const;
  void load_config(std::istream& is);

  /**
   * 将解向量写入到流中。
   * 
   * @param os 输出流
   * @param _sol 解向量
   */
  void dump_solution(std::ostream& os, const sol_t& _sol) const;

  /** 
   * 从流中载入解向量。
   * 
   * @param is 输入流
   * @param _sol 解向量
   */
  void load_solution(std::istream& is, sol_t& _sol) const;

  /** 
   * 将计算环境和解函数以二进制形式写入文件中，以后可将数据读入重新开始
   * 计算。
   * 
   * @param filename 输出数据的文件名。
   */
  void dump_data(const std::string& filename) const;
  /** 
   * 将二进制形式的数据文件读入内存，为开始新的计算设置好计算环境和初
   * 值。
   * 
   * @param filename 载入数据的文件名。
   */
  void load_data(const std::string& filename);

  /** 
   * 输出数据, 便于小蔡的程序读取. 与上面的函数功能有所重叠.
   * 
   * @param filename 
   */
  void save_solution( const std::string& filename ) const {
    std::ofstream os(filename.c_str());
    boost::archive::text_oarchive archive(os);
    archive & (*sol);
  }
//@}

/// @name 设置几个重要参数.
//@{ 
  /** 
   * @brief 重新设置计算的矩数。
   * 
   * @param order 矩数
   */
  void reset_order(int order){
    std::cerr << "Reset order as " << order << std::endl;
    ORDER = order;

    for (size_t i = 0; i < sol->size(); ++ i)
    {
      (*sol)[i].ResetOrder(ORDER);
      (*rhs)[i].ResetOrder(ORDER);
    }
  };
  /** 
   * @brief 重设网格
   *
   * @param _Nx 单元数
   */
  void reset_mesh( int _Nx = -1 );
  
  /** 
   * @brief 加密网格
   *
   * 只考虑对网格单元均分为二的加密
   * 
   * @param _levels 加密次数
   */
  void refine_mesh( int _levels = 1 );

  /**
   * 显然这个函数只是给多重网格方法设置参数的. 要不要把 p_solver 设为公
   * 有成员, 把这个函数交给具体的 cpp 去做? 不然在单网格的 solver 做这
   * 样一个接口总是比较怪异.
   * 
   */
  void reset_mesh_levels( int levels ) {
#ifdef MULTI_GRIG_SOLVER
    std::cerr << "Reset multi mesh as " << levels << " levels" << std::endl;
    p_solver->reset_levels( levels );
#else
    std::cerr << "Could not reset mesh levels: The solver is not a multigrid solver!" << std::endl;
#endif
  }
  void reset_mg_smoothing( int _n_smooth ) {
#ifdef MULTI_GRIG_SOLVER
    std::cerr << "Reset multi grid's smoothing step as " << _n_smooth << std::endl;
    p_solver->reset_smoothing( _n_smooth );
#else
    std::cerr << "Could not reset multi grid's smoothing step: The solver is not a multigrid solver!" << std::endl;
#endif
  }

  /**
   * 设置判断收敛的忍量. 注意 DualGate 类中的 _tol 暂时已经没有意义了
   */
  void reset_tol(const double _tol){
    std::cout << " Reset tol as :" << _tol << std::endl ;
    tol = _tol;
    
    p_solver->reset_tol( _tol );
  }

  void stop_run() { is_stop_run = true; }
//@}

private:
/// @name 判断是否收敛到稳定解的变量和函数
//@{
  /**
   * 返回 (残量 < 忍量)。
   */
  bool is_converged() const { return  _res < tol; }
//@}

/// @name 为计算环境分配和释放内存
//@{ 
  /** 
   * 分配内存。要求在参数配置文件的信息被读入了以后，应该就可以调用此
   * 函数建立起完整的计算环境。
   */
  void alloc_memory();
  /// 释放内存
  void free_memory(){   /* 清理环境需要释放哪些内存？*/
    // 现在使用的都是智能指针, 所以释放空间函数都不用了.
  };  
//@}
};

#include "SteadyState.template.h"

#endif // __STEADY_STATE_H__

/**
 * end of file
 * 
 */
