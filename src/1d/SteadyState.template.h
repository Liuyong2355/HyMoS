/**
 * @file   SteadyState.template.h
 * @author HU Zhicheng <huzhicheng1986@gmail.com>
 * @date   Mon Nov 12 14:37:41 2012
 * 
 * @brief  SteadyState 类的实现
 * 
 * 
 */

#ifndef __STEADY_STATE_TEMPLATE_H__
#define __STEADY_STATE_TEMPLATE_H__

#include "SteadyState.h"
#include <omp.h>
#include <cmath>
#include <stdexcept>
#include <atomic>

template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
void SteadyState<MESH, SOLUTION, FORCE, SOLVER, _BV>::run() {
  run_impl(true, 0, NULL, NULL);
}

template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
void SteadyState<MESH, SOLUTION, FORCE, SOLVER, _BV>::run_impl(
    bool legacy_io, unsigned int max_steps, std::vector<double>* residuals,
    const std::atomic<bool>* stop_requested,
    const std::function<void(unsigned int, double)>& progress)
{
  // Disabled local streams avoid touching the host's global stream buffers.
  std::ostream silent(NULL);
  std::ostream& log_out = legacy_io ? std::cout : silent;
  std::ostream& log_err = legacy_io ? std::cerr : silent;
  if (residuals) residuals->clear();
  if (!legacy_io && !is_read) throw std::runtime_error("solver is not configured");
  if( !is_read ) 
  {
    log_err << "Please read the config file (Default: input.txt) or load data \
(Default: DATA.dat) from previous rtest first!" << std::endl;
    return;
  }

  log_err << "Entering run ... " << std::endl;
  is_stop_run = false;

  // Wall-clock timing includes the actual elapsed runtime of the solver.
  timeval start_time;
  gettimeofday( &start_time, NULL ); // 获取当前时间

  static unsigned int legacy_n_run = 0;
  const unsigned int n_run = legacy_io ? ++legacy_n_run : 1;
  if( legacy_io && n_run == 1 )
  {
    std::ofstream clear_file( "res_step", std::ios::out );
    clear_file.close();
  }
  std::ofstream os_res;
  if (legacy_io) os_res.open("res_step", std::ios::app);
  os_res.precision( 12 );
  os_res << "% the n-th run: n = " << n_run << std::endl;
 
  /**
   * 若已分配内存, 就假设初值已经准备好了. 
   */
  if (sol.get() == NULL) initial_value();
  /// initial_value() 函数仅设置了 Boltzmann 方程解的初值. 为了计算初值
  /// 的残量大小, 需要先根据初值计算外力.
  p_force_handle->initial_value( *sol );

  unsigned int steps = 0;
  //OPENMP parallel
  omp_set_num_threads(n_thread);

  // 计算初值残量大小
  /**
   * @brief 残量计算方式选取
   *
   * M 增大时, 由于机器精度的限制, 标准的 L2 范数计算可能不再准确. 由于
   * 归一化这一步可能就有问题, 所以函数 GetResidualMaxNorm 也一样会不适
   * 用. 因而整个代码中将残量的计算方式改为使用 GetResidualPL2Norm 函数.
   * 示例: 在 Couette 流 M=26 时, 即使用 PL2 范数来判断, 残量的极限也就
   * 到 3e-9 的样子, 不过此前, 多重网格算法中此残量范数还是稳步递减的.
   * 而 MaxNorm 和 L2 范数, 一开始便不稳定, 忽大忽小. 到后期倒是也会慢
   * 慢减小, 但是 PL2 范数是 3e-9 的时候, L2 和 MaxNorm 都早就不能再降
   * 低, 已经稳定在 2e-4 左右了. 所以为了更清楚的观察 多重网格算法的性
   * 质, 使用 PL2 范数还是合理的吧. 
   *
   * 以下是测试比较各个范数行为的代码
   * @code
   *   sol_t res_sol( sol->size() );
   *   RESIDUAL::GetResidual( res_sol, *sol, *rhs, *p_mesh, *p_force_handle );
   *   double max_norm=0;
   *   RESIDUAL::GetResidualMaxNorm( max_norm, res_sol );
   *   double pl2_norm = 0;
   *   RESIDUAL::GetResidualPL2Norm( pl2_norm, res_sol, *p_mesh );
   *   RESIDUAL::GetResidualL2Norm( _res, res_sol, *p_mesh );
   * @endcode
   * 
   */
  RESIDUAL::GetResidualPL2Norm( _res, *sol, *rhs, *p_mesh, *p_force_handle );
    if (residuals) residuals->push_back(_res);
    if (!legacy_io && !std::isfinite(_res))
      throw std::runtime_error("non-finite residual");
  if (progress) progress(steps, _res);
  /// _res 的结果已经在 calculate_residualHelper 函数中计算出来

  timeval prepare_initial_time;
  gettimeofday( &prepare_initial_time, NULL );
  
  os_res << "% Elapsed time for preparing initial value of " << n_run << "-th run: " << getElapsedTime( start_time, prepare_initial_time ) << "s\n\n% " << n_run <<"-th run residual: \n";
  os_res << steps << " " << _res << " " << 0 << std::endl;

  log_out.precision( 12 );
  log_out << "steps = " << steps << ": residual = " << _res << std::endl;

  while( (! (is_stop_run || is_converged() || is_exit ||
              (stop_requested != NULL && stop_requested->load()) ) ) &&
         (legacy_io || steps < max_steps) )
	{
#if 0
    SaveMacroVars( *sol, *p_mesh, "Sol0th.dat" );

    p_force_handle->save_results( 0 );
    _BV::save_current( *sol, *p_mesh, 0 );
    getchar();
#endif
    ++steps;
    // if( steps == 2 ) break;
    /**
     * Boltzmann 方程求解器. 在某些问题中, 外力变化缓慢, 所以可以考虑给
     * 定外力下 Boltzmann 方程多算几步. 但是显然在一个不是最终外力下求
     * 解 Boltzmann 方程尽量精确在效率上也是不划算的.
     * 
     */
    for( unsigned int _i =0; _i<1; ++_i)
    {  
    p_solver->Solve( *sol, *rhs, p_mesh, p_force_handle );
    // 若考虑 full mg, 若这个 PostProcessing 确实在做事情, 可能还是要放
    // 到 Solve 里头去.
    p_solver->PostProcessing( *sol, p_mesh );
    /// 外力计算, 已经放到 p_solver->Solve() 中去, 所以 PostProcessing 是否也如此呢?
    // p_force_handle->CalcExternalForce( *sol );

    }
    /**
     * 此前, 有一个版本的程序在这里计算三次才保存一次数据, 在测试 Si 器
     * 件的时候, 单网格显式时间推进格式竟然能节省三分之一时间 (单网格迭
     * 代内部有一些计算范数的地方, 但也省了近五分之一), 那么说明这个计
     * 算范数的代价还是太大了些. 不过这并非改善效率的本质, 所以先把其他
     * 地方调好, 再来考虑这个问题.
     * 
     */

    // Compute the residual after the solver and post-processing step.
    RESIDUAL::GetResidualPL2Norm( _res, *sol, *rhs, *p_mesh, *p_force_handle );
    if (residuals) residuals->push_back(_res);
    if (!legacy_io && !std::isfinite(_res))
      throw std::runtime_error("non-finite residual");
    if (progress) progress(steps, _res);
	 timeval end_time;
    gettimeofday( &end_time, NULL );
    log_out << "steps = " << steps << ": residual = " << _res << std::endl;
    os_res << steps << " " << _res << " " << getElapsedTime( prepare_initial_time, end_time ) << std::endl;
//    _BV::save_current( *sol, *p_mesh, n_run );
    if (legacy_io) getControl();
	}
//  is_exit=true;
  if( is_exit )
  {
    log_err << "The Pro is stopped!" << std::endl;
  }
  else if (is_stop_run)
  {
    log_err << "Run is stopped." << std::endl;
  }

  timeval end_time;
  gettimeofday( &end_time, NULL );

  os_res << "% Elapsed time for " << n_run << "-th run (exclude the initialization time): " << getElapsedTime( prepare_initial_time, end_time ) << "s!\n% ";
  os_res << "% Total Elapsed time for " << n_run << "-th run: " << getElapsedTime( start_time, end_time ) << "s!\n% ";

  log_err << "The elapsed steps is " << steps << std::endl;

  if (!legacy_io) {
    is_stop_run = true;
    return;
  }

  if( is_converged() )
  {
    std::stringstream file_sol;
    file_sol << "Sol" << n_run << "th.dat";
    SaveMacroVars( *sol, *p_mesh, file_sol.str() );
    
    /// 以下三行存下整个分布函数给 BoundaryLayer.cpp 读取设置狄氏边界条件用.
    std::stringstream file_dis;
    file_dis << "DisSol" << n_run << "th.dat";
    SaveSolution( *sol, *p_mesh, file_dis.str() );

    p_force_handle->save_results( n_run );
    _BV::save_current( *sol, *p_mesh, n_run );
//    _BV::save_one_over_tau( *sol, *p_mesh, *p_force_handle, n_run);

    is_stop_run = true;
    /*
     *这一行用来使得收敛条件达到后，sol 还未释放前，把config以及 解的各
     *个信息dump 出来其实不需要，避免忘记dump_data用
     */
    std::stringstream file_data;
    file_data << "DATA" << n_run << "th.dat";
    dump_data( file_data.str() ) ;
    
    double max_norm =0;
    RESIDUAL::GetResidualMaxNorm( max_norm, *sol, *rhs, *p_mesh, *p_force_handle );
    log_out << "max_norm is " << max_norm << std::endl;
    os_res << "% max_norm is " << max_norm << "\n" << std::endl;

//    log_err << "please rename DATA_TMP.dat!!" <<std::endl ;
    log_err << "The elapsed time is " << getElapsedTime( start_time, end_time ) << "s!" << std::endl;
    log_out << "Leaving run due to convergence." << std::endl;
  }
  os_res << std::endl;
  os_res.close();
}

/**
 * 暂时只考虑重设网格为均匀网格.
 *
 * 有两种策略重设网格: 
 *    - 由外部设定具体的单元数, 不更改矩数
 *    - 根据保持耗散不变的原则, 更改矩数相应增加网格. 此种情况 @p _Nx
 *      用默认值表示, 并且需要先运行 reset_order() 函数. 还有, 此种方式
 *      产生的网格数可能不适于现在多重网格方法, 所以在用多重网格方法时,
 *      应避免使用此方式, 或者选择邻近的网格数, 使得可以粗化网格到只有
 *      几个单元.
 *
 * 若要保持相同耗散, 则要求 \f$ h \lambda_{\max} \equiv const \f$. 其中
 * \f$ \lambda_{\max} \f$ 便是 order + 1 阶 Hermite 多项式的最大根, 而
 * \f$ h \f$ 则是网格步长. 由于现在假设网格可能非均匀的, 此时考虑网格步
 * 长没什么意义,所以下面改用 \f$ \lambda_{\max} / Nx \equiv const \f$
 * 来衡量. 另外, 在矩数趋于无穷时, 大致有 \f$ Nx \f$ 和 \f$ \sqrt{M}
 * \f$ 成比例的关系.
 * 
 */
template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
void SteadyState<MESH, SOLUTION, FORCE, SOLVER, _BV>::reset_mesh( int _Nx )
{
  int Nx = _Nx; 
  if( _Nx <= 0 )
  {
    /**
     * 先设一个同时加密的参考值. 暂时是 \f$ M=4 \f$ 时, Nx=
     * 512, 即在均匀网格下, 参考的步长设为 0.5.
     * 
     */
    const static value_t C_M_ref = Transform<value_t>::MaxRootOfHermitePolynomial( 4+1 );
    const static size_t _Nx_ref = (int) (512);
    
    value_t C_M = Transform<value_t>::MaxRootOfHermitePolynomial( ORDER+1 );
      Nx = C_M * _Nx_ref / C_M_ref;
    std::cout << "C_M_ref = " << C_M_ref << ", C_M = " << C_M 
              << ", new Nx = " << Nx << std::endl;
  }

  ptr_mesh_t p_old_mesh = p_mesh;
  std::unique_ptr<sol_t> _sol = std::move( sol );

  free_memory( );
  /// 然后重新分配解空间
  p_mesh.reset( new mesh_t(p_old_mesh, Nx) );
  alloc_memory( );

  Operator::Interp<mesh_t, sol_t, _BV>( *p_old_mesh, *_sol, *p_mesh, *sol );
  (*rhs).Reinit( *sol );

  PROBLEM::Reinit( p_mesh );

  p_force_handle->Reinit( p_mesh );
  // 现在的程序设计, 在最密层网格重设时, 需要计算与V, rho 无关的右端项
  // 部分, 于是没有下面这句, 可能使程序使用此函数时导致一个 bug.
  p_force_handle->initial_value( *sol );

  p_solver->Reinit( p_mesh );
  
  std::cerr << "Reset the number of elements as " << Nx << std::endl;
}

template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
void SteadyState<MESH, SOLUTION, FORCE, SOLVER, _BV>::refine_mesh( int _levels )
{
  for( unsigned int _level = 0; _level < _levels; ++_level ){
    ptr_mesh_t p_old_mesh = p_mesh;
    std::unique_ptr<sol_t> _sol = std::move( sol );
    /**
     * free_memory() 和 alloc_memory() 做的事情实际很少, 所以这样放在此
     * 循环中, 也不会太没效率
     * 
     */
    free_memory( );
    /// 然后重新分配解空间
    p_mesh.reset( new mesh_t() );
    p_mesh.RefineFrom( *p_old_mesh );

    alloc_memory( );

    Operator::Interp<mesh_t, sol_t, _BV>( *p_old_mesh, *_sol, *p_mesh, *sol );
  }
  (*rhs).Reinit( *sol );

  PROBLEM::Reinit( p_mesh );

  p_force_handle->Reinit( p_mesh );
  // 现在的程序设计, 在最密层网格重设时, 需要计算与V, rho 无关的右端项
  // 部分, 于是没有下面这句, 可能使程序使用此函数时导致一个 bug.
  p_force_handle->initial_value( *sol );

  p_solver->Reinit( p_mesh );
  
  std::cerr << "Refine mesh is OK!" << std::endl;
}


template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
void SteadyState<MESH, SOLUTION, FORCE, SOLVER, _BV>::alloc_memory()
{
   // 在 read_config 阶段, 已经为网格 @p p_mesh 等分配好空间.
   unsigned int n_ele = p_mesh->n_ele();
   sol.reset(new sol_t(n_ele,ORDER));
   rhs.reset(new sol_t(n_ele,ORDER));
}

template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
void SteadyState<MESH, SOLUTION, FORCE, SOLVER, _BV>::read_config(const std::string & filename)
{
  std::cerr << "Reading config file " << filename << " ... " << std::flush;
  ConfigFile cf( filename );
  try
  {
    read_config(cf);
    std::cerr << "OK!" << std::endl;
  }
  catch (std::string s)
  {
    std::cerr << s << std::endl;
    exit(1);
  }
}

template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
void SteadyState<MESH, SOLUTION, FORCE, SOLVER, _BV>::read_config(ConfigFile& cf) {
    ORDER = (int)cf.Value("System","ORDER");
    n_thread = std::max(1, (int)cf.Value("System", "n_thread", 1.0));
    tol = (value_t)cf.Value("Error", "Tol" );

    if( p_mesh.get() == NULL )
      p_mesh.reset( new mesh_t( cf ) );
    else
      p_mesh->read_config( cf );

    // 重新设置 Poisson 方程右端项计算以后, 这一句要提前
    PROBLEM::read_config( p_mesh, cf );
    
    if( p_force_handle.get() == NULL ) 
      p_force_handle.reset( new force_t( p_mesh, cf ) );
    else
      p_force_handle->read_config( p_mesh, cf );

    if( p_solver.get() == NULL )
      p_solver.reset( new solver_t( p_mesh, cf ) );
    else
      p_solver->read_config( p_mesh, cf );

    is_read = true;
}

template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
void SteadyState<MESH, SOLUTION, FORCE, SOLVER, _BV>::print_config() const 
{
  std::cerr << " Order = " << ORDER << std::endl;
  std::cerr << " Nx = " << p_mesh->n_ele() << std::endl;
  std::cerr << " n_thread = " << n_thread << std::endl;
  std::cerr << " tol = " << tol << std::endl;

  _BV::print_config();

  p_solver->print_config();
  p_force_handle->print_config();

  std::cerr << std::flush ;
}

template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
void SteadyState<MESH, SOLUTION, FORCE, SOLVER, _BV>::dump_config(std::ostream& os) const {
  os.write((const char *)&ORDER, sizeof(int));
  /// 这里, dim 这个数应该是不需要存储的.
  os.write((const char *)&tol, sizeof(double));
}

template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
void SteadyState<MESH, SOLUTION, FORCE, SOLVER, _BV>::load_config(std::istream& is) {
  is.read((char *)&ORDER, sizeof(int));
  is.read((char *)&tol, sizeof(double));
}

template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
void SteadyState<MESH, SOLUTION, FORCE, SOLVER, _BV>::dump_solution(std::ostream& os,
                                         const sol_t& _sol) const
{
  size_t n = _sol.size();
  std::cout << " dump n : " << n << std::endl ; 
  os.write((const char *)&n, sizeof(size_t));
  for (u_int i = 0; i < n; ++ i)
  {
    const dis_t& dis = _sol[i];
    for( size_t j=0; j<dim; ++j )
      os.write((const char *)&(dis.Center()[j]), sizeof(double));
    os.write((const char *)&(dis.Scaling()), sizeof(double));
    typename dis_t::ConstIterator the_entry = dis.begin(), end_entry = dis.end();
    for (; the_entry != end_entry; ++ the_entry)
    {
      os.write((const char *)&(*the_entry), sizeof(double));
    }
  }
}

template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
void SteadyState<MESH, SOLUTION, FORCE, SOLVER, _BV>::load_solution(std::istream& is,
                                         sol_t& _sol) const
{
  size_t n;
  is.read((char *)&n, sizeof(size_t));
  std::cout << "Nx : " << p_mesh->n_ele() << " n" << n << std::endl ;
  assert((n == _sol.size()));
  for (u_int i = 0; i < n; ++ i)
  {
    dis_t& dis = _sol[i];
    for( size_t j=0; j<dim; ++j )
      is.read((char *)&(dis.Center()[j]), sizeof(double));
    is.read((char *)&(dis.Scaling()), sizeof(double));
    dis.ResetOrder(ORDER);
    typename dis_t::Iterator the_entry = dis.begin(), end_entry = dis.end();
    for (; the_entry != end_entry; ++ the_entry)
    {
      is.read((char *)&(*the_entry), sizeof(double));
    }
  }
}


template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
void SteadyState<MESH, SOLUTION, FORCE, SOLVER, _BV>::dump_data(const std::string& filename) const
{
  std::cerr << "Dumping data to file " << filename << " ... " << std::flush;
  std::ofstream os(filename.c_str(), std::ios::binary);
  p_mesh->dump_data( os );

  dump_config(os);
  dump_solution(os, *sol);

  _BV::dump_data( os );

  p_force_handle->dump_data( os );
  p_solver->dump_data( os );

  std::cerr << "OK!" << std::endl;
}

template <typename MESH, typename SOLUTION, typename FORCE, typename SOLVER, typename _BV>
void SteadyState<MESH, SOLUTION, FORCE, SOLVER, _BV>::load_data(const std::string& filename)
{
  std::cerr << "Loading data from file " << filename << " ... " << std::flush;
  free_memory();

  std::ifstream is(filename.c_str(), std::ios::binary);

  p_mesh.reset( new mesh_t(is) );

  load_config(is);
  alloc_memory();
  load_solution(is, *sol);
  (*rhs).Reinit( *sol );

  _BV::load_data( p_mesh, is );

  p_force_handle.reset( new force_t( p_mesh, is ) );
  p_solver.reset( new solver_t(p_mesh, is) );

  is_read = true;
  std::cerr << "OK!" << std::endl;
}


#endif /// __STEADY_STATE_TEMPLATE_H__

/**
 * end of file
 * 
 */
