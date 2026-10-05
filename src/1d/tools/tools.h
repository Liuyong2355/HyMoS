/**
 * @file   tools.h
 * @author HU Zhicheng <huzhicheng1986@gmail.com>
 * @date   Fri Nov 29 15:38:43 2013
 * 
 * @brief  不适宜作为类成员函数的部分通用函数
 * 
 */


#ifndef __STEADY_STATE_TOOLS_H__
#define __STEADY_STATE_TOOLS_H__

#include <iostream>
#include <vector>
#include <sys/time.h>

#include <Eigen/Dense>

#include <NRxx/CoefTable.h>
#include <NRxx/Distribution.h>

#include "../Solution.h"


/** 
 * @brief 计算两个时间的差异
 *
 * time_t 精度只有 s, 为了获得更好的运行时间精度, 需要用 sys/time.h 中
 * 的 timeval 类以及 gettimeofday( &timeval, NULL ) 函数. timeval 有两
 * 个成员, tv_sec 和 tv_usec, 需要组合才能算出具体运行时间. 注意
 * tv_usec 是微秒 (1/1000000 s).
 * 
 * @param _start 起始时间
 * @param _end 终止时间
 * 
 * @return 返回时间差 (单位 s)
 */
inline double getElapsedTime( const timeval& _start, const timeval& _end ) { 
  return _end.tv_sec - _start.tv_sec + (_end.tv_usec - _start.tv_usec) / 1000000.;  
}

/** 
 * @brief 计算 i 的阶乘
 * 
 * @note 成笔记本的 64 位机器中, int 型能正确表示的最大阶乘是 12!, 而其
 * 它整型 (long int, long long int 以及它们的 unsigned 版), 也只能正确
 * 表示到 20!. 因而实际计算中, 应该避免使用此函数, 以免舍入误差可能带来
 * 的影响. 避免的方式自然是考虑将阶乘的计算和其它部分的计算穿插起来, 使
 * 得目标值总不会太大.
 *
 * 由于阶乘非常快就超出了整型所能表示的范围, 因而这里设置了两个模版参数: 
 *
 * - INT 用来表示输入的参数, 普通的整型就足够了 (总是要求这个参数是整型);
 *
 * - LINT 用来表示返回值类型, 由于所有整型所能表示的阶乘都非常小, 所以
 *   这个类型更多的应该设为 double 型以上. 虽然此时只有十几位的有效精度
 *   (i>20 的阶乘就不精确成立了), 但是所能表示的阶乘数却大了许多, 即使
 *   是 \f$ i=100\f$, 也能基本正确的表示为 9.3326e+157. 事实上, 对
 *   double 型而言, 能表示的最大阶乘是 170!, 而 long double 型所能表示
 *   的阶乘就能更大了.
 *
 * @note gsl 库也有计算阶乘的函数, 不知效率如何.
 * 
 * @param i 要求 \f$ i\geq 0 \f$
 * 
 * @return 
 */
template <typename INT, typename LINT>
inline LINT FactorialNumber( const INT i ){
  // if ( i>20 ){// 成 64位机器能算的最大阶乘
  //   std::cerr << "FactorialNumber: Out of range " << std::endl;
  //   getchar();
  // }

  static std::map<INT,LINT> factorial;
  typename std::map<INT,LINT>::iterator it=factorial.find(i);
  if ( it != factorial.end() )
    return it->second;
  
  if( i > 1 )
  {
    factorial.insert( typename std::map<INT,LINT>::value_type( i, i*FactorialNumber<INT,LINT>( i-1 )) );    
  }
  else if ( (i==0) || (i==1) )
  {
    factorial.insert( typename std::map<INT,LINT>::value_type(i, 1) );
  }
  else
  {
    std::cerr << "FactorialNumber: No definition for negative number!  i = " << i << std::endl;
    return 0;
  }
  
  return FactorialNumber<INT,LINT>( i );
}

/** 
 * @brief \f$ \sqrt{n !} \f$ 的计算函数
 * 
 * 开根号后往往不再是整数, 所以返回值模板类换个名字, 为 DOUBLE, 不过在
 * 成的系统里, double, long double, long long double, 所能表示的精度都
 * 是一样的. 
 * 
 * 单独写一个此函数, 而不是从 FactorialNumber<INT,LINT>() 函数先得到n!
 * 再开根号是有好处的, 至少, 所能计算的最大 n 可能可以更大 (double 型能
 * 算的最大阶乘是 170!, 但至少能表示到 \f$\sqrt{200 !}\f$). 不过不知道
 * 多次开根号相乘和相乘后再开根号的精度如何.
 * 
 * @param i 即 \f$ n \f$
 * 
 * @return 
 */
template <typename INT, typename DOUBLE>
inline DOUBLE SqrtFactorialNumber( const INT i ){
  static std::map<INT,DOUBLE> sqrt_factorial;
  typename std::map<INT,DOUBLE>::iterator it=sqrt_factorial.find(i);
  if ( it != sqrt_factorial.end() )
    return it->second;
  
  if( i > 1 )
  {
    sqrt_factorial.insert( typename std::map<INT,DOUBLE>::value_type( i, sqrt(i) * SqrtFactorialNumber<INT,DOUBLE>( i-1 )) );    
  }
  else if ( (i==0) || (i==1) )
  {
    sqrt_factorial.insert( typename std::map<INT,DOUBLE>::value_type(i, 1) );
  }
  else
  {
    std::cerr << "SqrtFactorialNumber: No definition for negative number!  i = " << i << std::endl;
    return 0;
  }
  
  return SqrtFactorialNumber<INT,DOUBLE>( i );
}

/** 
 * @brief 分布函数的归一化
 *
 * 在当前的展开表示中, 由于基函数 \f$ \mathcal{H}_{\theta,
 * \alpha} \f$ 吸收了一个 \f$ \theta^{-|\alpha| / 2}\f$ 的因子, 故而,
 * 很明显有这样的量纲关系: \f$ f_{\alpha} \sim \theta^{|\alpha|/2} \f$
 * 或 \f$ f_{\alpha} \sim \rho \theta^{|\alpha|/2} \f$ (\f$ f_0 \f$ 即
 * 是 \f$ \rho \f$). 进一步地, 若考虑采用单位化的标准正交基, 则还会有一
 * 个关于阶乘的因子, 也就是 \f$ f_{\alpha} \sim \frac{\rho
 * \theta^{|\alpha|}}{\sqrt{\alpha !}} \f$. 于是为了更好的比较矩系数,
 * 应该对其做一个归一化, 也即
 * \f[
 * \hat{f}_{\alpha} = f_{\alpha} \theta^{-|\alpha|/2} \sqrt{\alpha !}
 * \f]
 * 至少也应是 \f$\hat{f}_{\alpha} = f_{\alpha}
 * \theta^{-|\alpha|/2}\f$, 这里不除以 \f$ \rho \f$ 是因为传入的 @p
 * dis 未必能取到, 而且上面这样已经做到各矩系数在同一量纲上了.
 * 
 * @param normalized_dis 有时候, 我们仅需要前面几阶矩的信息, 所以
 *        normalized_dis 的阶数放到程序外部定义, 注意啊其阶数必须小于等
 *        于 @p dis 的阶数.
 * @param dis
 */
template <typename DISTRIBUTION>
inline void NormalizedDistribution( DISTRIBUTION& normalized_dis, 
                                    const DISTRIBUTION& dis ) {
  static int dim = DISTRIBUTION::dim;
  double one_over_s = 1. / dis.Scaling(); // theta^{-1/2}
  double coe = 1.;
  unsigned int old_order = 0;

  normalized_dis.Reinit( dis.Center(), dis.Scaling() );
  typename DISTRIBUTION::ConstIterator it_dis = dis.begin();
  //typename DISTRIBUTION::ConstIterator it_end = dis.end();
  typename DISTRIBUTION::Iterator it_normalized_dis = normalized_dis.begin();
  const typename DISTRIBUTION::Iterator the_end = normalized_dis.end();
  for( ; it_normalized_dis != the_end; ++ it_dis, ++ it_normalized_dis )
  {
    unsigned int order = it_dis.Order();
    if( order > old_order )
    {// 根据分布函数的存储方式, order 总是递增, 而且, order 比
     // old_order 至多大 1
      coe *= one_over_s;
      old_order = order;
    }
    double local_value = coe;
    const typename DISTRIBUTION::Indices& ind = it_dis.GetIndices();
    for( int i=0; i<dim; ++i )
      local_value *= SqrtFactorialNumber<int, double>( ind[i] );

    /**
     * 这种先将 @p local_value 算好, 再乘到矩系数上去的方式, 对于高阶矩
     * 可能会损失较多的精度. 之前的二范数在 GetOrder() 为 30 时精度只有
     * 1e-12 可能就是这个原因. 然而, 除非仅考虑 dim=1 的情形, 否则很难
     * 高效的完成计算, 又顾及精度的问题 (避免机器精度损失). 所以还是先
     * 以这种土土的方式来计算.
     *
     * 如果不乘那个阶乘的开根号, 不知道各阶矩大小数量级又如何. 
     */
    *it_normalized_dis = *it_dis * local_value;
  }
}

/** 
 * @brief 比较矩系数绝对值的最大值, 也即矩系数向量的无穷范数.
 *
 * 范数计算的准确与否, 与实质上方法的收敛速度并无关系, 仅仅可能影响的是
 * 算法收敛的判别. 虽然如此, 还是希望有一个比较好的范数来衡量所有变量.
 * 观察数据可以看到, @p dis 中高阶矩的矩系数是低阶矩的高阶小量, 因而直
 * 接的比较实际上仅仅体现了低阶矩的信息. 所以应先对分布函数做归一化后再
 * 行判断.
 *
 * @param dis 分布函数, 往往不在自己的标准基下展开
 * @param _max 存储矩系数的最大值, 需要先做初始化再传进来. 一般初始化
 * 为 0.0, 此时返回值就是矩系数的最大值. 其他情形, 返回的是矩系数和这
 * 个初始化值中绝对值最大者.
 */
template <typename DISTRIBUTION>
void MaxMomentsValue( const DISTRIBUTION& dis, double& _max=0. ) {
  DISTRIBUTION normalized_dis( dis.GetOrder() );
  NormalizedDistribution( normalized_dis, dis );

  typename DISTRIBUTION::Iterator it_dis = normalized_dis.begin();
  typename DISTRIBUTION::Iterator it_end = normalized_dis.end();
  for( ; it_dis != it_end; ++it_dis )
    _max = std::max( _max, fabs(*it_dis) );
}

/** 
 * @brief 矩系数的 L1 范数
 *
 * 简单地讲, 就是归一化后矩系数的绝对值之和.
 *
 * @note 注意这并不是分布函数的 L1 范数. 若根据范数定义和基函数的正交性,
 * 分布函数的 L1 范数就是 0 阶矩系数的某个常数倍.
 * 
 * @param dis 同样, 分布函数往往不在自己的标准基下展开
 * 
 * @return L1 范数
 */
template <typename DISTRIBUTION>
double L1NormMomentsValue( const DISTRIBUTION& dis ) {
  DISTRIBUTION normalized_dis( dis.GetOrder() );
  NormalizedDistribution( normalized_dis, dis );

  double l1_norm = 0.0;
  typename DISTRIBUTION::Iterator it_dis = normalized_dis.begin();
  typename DISTRIBUTION::Iterator it_end = normalized_dis.end();
  for( ; it_dis != it_end; ++it_dis )
    l1_norm += fabs(*it_dis);

  return l1_norm;
}

/** 
 * @brief 分布函数的 L2 范数
 *
 * 主要是重复代码提取, 这个函数并不做归一化, 仅是矩系数平方和开根号再乘
 * 以一个常数. 于是在外部使用并没有什么意义, 目前仅在函数 L2NormDis 和
 * PL2NormDis 中调用. 为区别 L2NormDis 等函数名, 这里起名做点差异, 为
 * L2NormSDis, 其中 S 就表示 Simplified 吧.
 * 
 * @param dis 
 * 
 * @return 
 */
template <typename DISTRIBUTION>
inline double L2NormSDis( const DISTRIBUTION& normalized_dis ) {
  double l2_norm = 0.0;
  typename DISTRIBUTION::ConstIterator it_dis = normalized_dis.begin();
  typename DISTRIBUTION::ConstIterator the_end = normalized_dis.end();
  for( ; it_dis != the_end; ++it_dis )
  {
    l2_norm += (*it_dis) * (*it_dis);
  }
  // 归一化的时候只是做到了量纲统一, 所以要得到真正的二范数, 还要乘一个
  // 常数.
  static int dim = DISTRIBUTION::dim;
  static double pi_c = pow( 2*M_PI, - 0.25*dim );

  return sqrt(l2_norm) * pi_c * pow( normalized_dis.Scaling(), -dim );  
}

/** 
 * @brief 归一化后分布函数的 L2 范数
 *
 * L2 范数应为 sqrt( sum (C_k * f_k * f_k) ), C_k 的公式参见
 * essay/method.tex. 我们知道, 当一个函数在标准正交基下展开时, 其二范数
 * 的平方就是简单的展开系数平方和. 因而, NormalizedDistribution() 以后
 * 分布函数矩系数的平方和与这里的 L2 范数的平方, 实际上就差一个常数.
 *
 * @note
 *
 * 在器件模拟的时候, 无量纲的 theta 是一个小于 1 的正数, 于是按照公式计
 * 算的 C_k 将是一个极大的数 (当 M 较大时, 是一个极大的阶乘再乘以
 * theta^{-|\alpha|}), 难道高阶矩真的就如此小么? 显然由于机器精度的原因
 * 二范数的计算很难保证准确了 (Couette 中 M 大点也就如此), 这个时候应该
 * 考虑使用 PL2NormDis() 函数来替代计算.
 *
 * @param dis 
 * 
 * @return 
 */
template <typename DISTRIBUTION>
double L2NormDis( const DISTRIBUTION& dis ) {
  DISTRIBUTION normalized_dis( dis.GetOrder() );
  NormalizedDistribution( normalized_dis, dis );
  return L2NormSDis( normalized_dis );
}
/** 
 * @brief 参见 L2NormDis 函数
 * 
 * 只是 L2NormDis 求和中的部分和 (按程序中的设置, 就是仅计算指标和小于
 * 等于3的那些系数)
 *
 * @param dis 
 * 
 * @return 
 */
template <typename DISTRIBUTION>
double PL2NormDis( const DISTRIBUTION& dis ) {
  // 下面这个构造函数参数必须是 unsigned int
  DISTRIBUTION normalized_dis( std::min( 3u, dis.GetOrder() ) );
  NormalizedDistribution( normalized_dis, dis );
  return L2NormSDis( normalized_dis );
}

/** 
 * @brief 多维指标的阶数
 *
 * NRxx 中指标类并未提供其和函数, 只在矩系数类中为迭代器提供了 Order()
 * 函数.
 * 
 * @param _vec 
 * 
 * @return 
 */
template <int dim> 
inline unsigned int SumOfIndices( const NRXX::CoefTableIndices<dim>& _vec ) {
  int order = 0;
  for( int i=0; i<dim; ++i )
    order += _vec[i];
  return order;
}

/** 
 * @brief 返回二阶矩最后一个在一维数组中的指标
 */
template <int dim> inline unsigned int EndSecondMomentsNumber();
template <> inline unsigned int EndSecondMomentsNumber<1> () { return 2; }
template <> inline unsigned int EndSecondMomentsNumber<2> () { return 5; }
template <> inline unsigned int EndSecondMomentsNumber<3> () { return 9; }

/** 
 * @brief 判断所扰动的二阶矩是否会改变分布函数标准展开的参数系
 *
 * 二阶矩指标范围是 i >=dim+1 且 i <= EndSecondMomentsNumber<dim>(). 其
 * 中, 只有扰动的是 \$ f_{2e_d} \$ 时才会改变温度, 需要重新标准展开. 在
 * NRxx 中, 二阶矩的排序是 (三维) 是: (0,0,0), (1,0,0), (0,1,0),
 * (0,0,1), (2,0,0), (1,1,0), (1,0,1), (0,2,0), (0,1,1), (0,0,2), 由此
 * 可知改变标准展开的指标 (见程序实现) 
 * 
 * @param i 要求用户保证是二阶矩的指标
 * 
 * @return 若扰动不改变参数系, 返回 true; 否则, 返回 false.
 */
template <int dim> inline bool IsStillStdExpansion( unsigned int i );
template <> inline bool IsStillStdExpansion<1> ( unsigned int i ) {
  // 此时二阶矩只有一个, 故而必改变标准展开参数系
  return false;
}
template <> inline bool IsStillStdExpansion<2> ( unsigned int i ) {
  return !( (i==3) || (i==5) ); // dim+1 = 3
}
template <> inline bool IsStillStdExpansion<3> ( unsigned int i ) {
  return !( (i==4) || (i==7) || (i==9) );
}

/** 
 * @brief 一阶矩系数平方和
 *
 * 由于迭代器可能还来自其他向量, 所以做这个模板函数
 * 
 * @param first 这是所有矩系数的起始迭代器 (0阶矩系数)
 * 
 * @return 
 */
template <typename Iterator>
inline double SumSquareOfFirstMoments( unsigned int dim, Iterator first ){
  ++first;
  double _sum=0.;
  for( unsigned int i=0; i<dim; ++i, ++first )
    _sum += (*first) * (*first);
  return _sum;
}

/** 
 * @brief 计算 \f$ \sum f_{e_{2d}} \f$
 * 
 * @param vec 关于矩系数的一个向量. 目前使用此函数的 vec 可能是一个
 * dealii的 Vector, 也可能是 NRxx 的 Distribution.
 * 
 * @return 求和值
 */
template <unsigned int dim, typename VEC> inline double SumDiagOfSecondMoments( const VEC& vec );
template <> 
inline double SumDiagOfSecondMoments<1, Eigen::VectorXd > ( const Eigen::VectorXd& vec ) { return vec( 2 ); } // dim+1
template <> 
inline double SumDiagOfSecondMoments<2, Eigen::VectorXd > ( const Eigen::VectorXd& vec ) { return vec(3) + vec(5); }
template <> 
inline double SumDiagOfSecondMoments<3, Eigen::VectorXd > ( const Eigen::VectorXd& vec ) { return vec(4) + vec(7) + vec(9); }
// 由于模板函数不能部分特化, 否则的话这两函数是完全可以写到一起的, 只要
// 用迭代器就可以了. 好罢, 小蔡的迭代器没有 + 的运算
template <> 
inline double SumDiagOfSecondMoments<1, NRXX::Distribution<1,double,0> > ( const NRXX::Distribution<1,double,0>& _dis ) { return *(_dis.begin(2)); } // dim+1
template <> 
inline double SumDiagOfSecondMoments<2, NRXX::Distribution<2,double,0> > ( const NRXX::Distribution<2,double,0>& _dis ) { 
  static typename NRXX::Distribution<2,double,0>::Indices ind1(0,2);
  return (*(_dis.begin(2))) + _dis(ind1); 
} 
template <> 
inline double SumDiagOfSecondMoments<3, NRXX::Distribution<3,double,0> > ( const NRXX::Distribution<3,double,0>& _dis ) { 
  static typename NRXX::Distribution<3,double,0>::Indices ind1(0,2,0); 
  static typename NRXX::Distribution<3,double,0>::Indices ind2(0,0,2); 
  return (*(_dis.begin(2))) + _dis(ind1) + _dis(ind2); 
} 

template <typename _T>
inline _T minmod( const _T& _x1, const _T& _x2 ) {
  if( _x1 * _x2 <= 0 )
    return 0.0;

  return ( fabs(_x1) < fabs(_x2) ) ? _x1 : _x2;
}

/** 
 * @brief 基本外推公式
 *
 * @see Extrap( sol_t&, const sol_t&, const sol_t&, const value_t&,
 * const value_t& );
 * 
 * @param val_H 
 * @param val_h 
 * @param _ratio \f$ p^m \f$
 * 
 * @return 
 */
template <typename _T>
inline _T Extrap( const _T& val_H, const _T& val_h, const _T& _ratio ) {
  return ( val_h - _ratio * val_H ) / (1 - _ratio);
}

template <typename Iterator>
void save_vector( std::ostream& os, Iterator first, const Iterator last )
{
  for (; first != last; ++first )
  {
    os << *first << ' ';
  }
}

/**
 * @brief 分布函数输出运算符 << 重载
 *
 * 若将 NRXX::Distribution<dim, double, 0> 整个作为模板类型, 与内部定义
 * 的 << 就会混淆. 为了函数更加通用, 将其中文字输出去掉, 输出格式是按序
 * 输出一行, 先是展开参数系 u, s=sqrt(theta), 然后是展开系数.
 * 
 */
template <int dim>
std::ostream& operator<<( std::ostream& os, 
                          const NRXX::Distribution<dim, double, 0>& _dis )
{
  typedef NRXX::Distribution<dim, double, 0> dis_t;
  // os << " Center is : ";
  const typename dis_t::Velocity& u = _dis.Center();
  for( size_t i=0; i<u.size(); ++i )
    os << u[i] << " ";
  // os << "\n Scaling is : " << _dis.Scaling() << std::endl;
  // os << "The moments is : " ;
  os << _dis.Scaling() << " ";
  save_vector( os, _dis.begin(), _dis.end() );
  // typename dis_t::ConstIterator the_dis = _dis.begin();
  // const typename dis_t::ConstIterator the_end = _dis.end();
  // for( ; the_dis != the_end; ++ the_dis )
  //   os << *the_dis << " ";
  // // os << std::endl;

  return os;
}

template <typename DISTRIBUTION, typename MESH>
void SaveSolution( const Solution<DISTRIBUTION>& _sol, 
                   const MESH& mesh, const std::string& file_sol )
{
  std::ofstream os( file_sol.c_str(), std::ios::out );
  os.precision( 12 );

  os << "% 1D distribution solution generated by HyMoS" << std::endl;
  os << "% column 1 = x, remaining values = full distribution data in internal storage order" << std::endl;

  for ( size_t i=0; i<_sol.size(); ++i )
  {
    os << mesh.Center(i) << ' ' << _sol[i] << std::endl;
  }
  os.close() ;
}

template <typename DISTRIBUTION, typename MESH>
void SaveMacroVars( const Solution<DISTRIBUTION>& _sol, 
                    const MESH& mesh, const std::string& file_sol )
{
  std::ofstream os( file_sol.c_str(), std::ios::out );
  os.precision( 12 );

  double total_mass = 0.;
  double pv[DISTRIBUTION::dim + 2];
  double heat_flux[DISTRIBUTION::dim];

  os << "% 1D macroscopic solution generated by HyMoS" << std::endl;
  os << "% columns: x rho u_1 u_2 u_3 T tauxx tauxy tauyy qx qy" << std::endl;
  os << " variables= \"x\", \"rho\", "
     << " \"u_1\", \"u_2\", \"u_3\", "
     << " \"T\", \"tauxx\", \"tauxy\", \"tauyy\", "
     << " \"qx\", \"qy\" " << std::endl;
  os << "zone i = " << _sol.size() << ", datapacking = \"point\" " << std::endl;

  for( size_t i = 0 ; i < _sol.size() ; ++i ){
    _sol[i].PrimitiveVars( pv ) ;
    total_mass += pv[0] * mesh.Length(i)  ;
    os << mesh.Center(i) << ' ';
    os << pv[0] << ' '
       << pv[1] << ' '
       << pv[2] << ' '
       << pv[3] << ' '
       << pv[4] << ' ';

    _sol[i].HeatFlux( heat_flux );
    typename DISTRIBUTION::Indices ind_sigma_xx(2, 0, 0);
    typename DISTRIBUTION::Indices ind_sigma_xy(1, 1, 0);
    typename DISTRIBUTION::Indices ind_sigma_yy(0, 2, 0);
    os << _sol[i](ind_sigma_xx, 0) << ' '
       << _sol[i](ind_sigma_xy, 0) << ' '
       << _sol[i](ind_sigma_yy, 0) << ' '
       << heat_flux[0] << ' '
       << heat_flux[1];

    os << std::endl;
  }
  /// octave 能认 # 开头的注释, 但 matlab 似乎只认 % 开头的注释
  os << "%# , total_mass = " << total_mass;
  os.close() ;
}

/**
 * @note 非均匀网格这个函数可能有些问题
 * 
 */
template <typename DISTRIBUTION, typename MESH>
void SaveCurrent( const Eigen::VectorXd &current, const MESH& mesh, 
                    const Solution<DISTRIBUTION>& _sol, 
                    const std::string& file_sol )
{
/**
 * 但是, 当考虑速度多维时, 电流也应该是多维的, 所以电流也要修改
 * 
 */

  std::ofstream os( file_sol.c_str(), std::ios::out );
  os.precision( 12 );
  double curr_all = 0., rho_recon = 0.;

  for( size_t i=0; i<current.size(); ++i )
  {
    if( i==0 )
    {
      /**
       * 边界条件还在类中, 作为全局函数, 这里不是很好取, 所以先做一个简
       * 单处理, 所以区域边界, 就用当前单元的密度, 不在重构.
       * 
       */
      // DISTRIBUTION lbv; 
      // LeftBoundaryValue( lbv, _sol );
      // rho_recon = (lbv.front() + _sol[i].front() ) / 2.;
      rho_recon = _sol[i].front();
    }
    else if ( i==_sol.size() )
    {
      // DISTRIBUTION rbv;
      // RightBoundaryValue( rbv, _sol );
      // rho_recon = (_sol[i-1].front() + rbv.front() ) / 2.;
      rho_recon = _sol[i-1].front();
    }
    else 
    {
      rho_recon = ( _sol[i-1].front() + _sol[i].front() ) / 2.;
    }
    curr_all += current(i);
    os << mesh[i] << " " << current(i) << " ";
    /*
     * _sol 存储的速度会出现负值, 这说是不合理的. 李老师说重构面上的密
     * 度根据电流计算速度看看
     * 
     */
    os << current(i) / rho_recon << std::endl;
  }
  /// octave 能认 # 开头的注释, 但 matlab 似乎只认 % 开头的注释
  os << "%# , curr_ave = " << curr_all / current.size() ;
  os.close() ;
}

template <typename DISTRIBUTION>
bool ISNAN( const DISTRIBUTION& dis )
{
  bool flag = false;
  const typename DISTRIBUTION::Velocity& center = dis.Center();
  const typename DISTRIBUTION::value_type& scaling = dis.Scaling();

  for( int j=0; j<DISTRIBUTION::dim; ++j )
    if( isnan( center[j] ) )
      flag = true;
  if( isnan(scaling) )
    flag = true;

  typename DISTRIBUTION::ConstIterator it_dis = dis.begin();
  const typename DISTRIBUTION::ConstIterator it_end = dis.end();
  for( size_t j=0; it_dis != it_end; ++it_dis, ++j )
  {
    if( isnan( *it_dis ) )
      flag = true;
  }

  return flag;
}

template <typename DISTRIBUTION>
bool ISNAN( const Solution<DISTRIBUTION>& sol ) {
  for( unsigned int i=0; i<sol.size(); ++i )
  {
    if( ISNAN( sol[i] ) )
    {
      std::cout << "NAN value in i = " << i << "\n" << sol[i];
      return true;
    }
  }
  return false;
}

#if 0 // 以后可能还会用到的实现
/** 
 * 这个函数实际上求解的是这样一个方程 \f[ -\Delta \tilde{V} =
 * \tilde{f}_i, \quad \tilde{V}(x_l)=\tilde{V}(x_r)=0, \f] 然后再把电
 * 势的增量更新到各单元上. 这里的 \f$ \tilde{f}_i \f$ 就是传入的参数
 * drho. 这个问题的差不多可以写出准确的解, @see essay/poisson.tm.
 *
 * 由于 Poisson 方程的全局依懒性, 一个单元上的扰动可引起整个区域上电势
 * 的变化. 但是按照这个准确的计算, 当所有单元都有密度变化的时候,
 * Poisson 方程的计算量就达到 \f$ O(N^2) \f$, 这应该是不可忍受的. 所以
 * 引入 width 这个参数, 每个单元的密度变化, 仅更新其周围 2*width+1 个单
 * 元上的电势. 这样就把这里最终的复杂度控制为 \f$
 * O(2\times\text{width}\times N) \f$. 当然这样会引入不少误差, 不知道导
 * 致的电势会不会很难看, 所以可能在所有单元密度都变化以后, 需要去解一个
 * 全局的三对角线性方程组, 把这里丢掉的电势部分的补回来. 因为追赶法是
 * \f$ O(N) \f$ 的计算量, 所以这个如果有效果的话, 计算量还是可接受的.
 *
 * 还有就是这里可能可以使用 FFT 计算, 但这样的复杂度还是有 \f$ O(N
 * \log N). \f$
 *
 * @note 注意这个函数是在 V 定义在单元上时的计算, 如果是定义在点上, 可
 * 能就不能用了
 *
 */
void UpdatePoissonPotential( const value_t& drho, const size_t _i, const int width )
{
  const value_t& len_domain = Len;
  const value_t& a = mesh[0];
  const value_t& b = mesh[Nx];
  const value_t& len_i = mesh.Length(_i);
  const value_t& center_i = mesh.Center(_i);

  size_t IL = 0;
  size_t IR = Nx-1; /// 最大的单元指标
  if( width >= 0 )
  {
    IL = std::max( (int)_i-width, 0 );
    IR = std::min( _i+width, IR );
  }
  /// 更新指标在 IL 和 IR 之间的电势, 分为三段
  for(size_t j=IL; j<_i; ++j )
  {
    V[j] += (mesh.Center(j) - a) / len_domain * (b - center_i) * len_i * drho;
  }
  for( size_t j=_i+1; j<=IR; ++j )
  {
    V[j] += (mesh.Center(j) - b) / len_domain * (a - center_i) * len_i * drho;
  }
  V[_i] += drho * len_i * ( ((b+a-center_i)*center_i - a*b) / len_domain - len_i / 6 );
  
  /**
   * @note FFT 的加速实质就是减少重复的计算量. 观察上面这三个式子, i 和
   * j 两个指标确实存在某种对称关系, 于是当两个指标同时在所有单元循环的
   * 时候, 使用 FFT 确实可以加速计算. 但是现在的境况每当一个单元密度发
   * 生变化, 便及时更新下一个单元的电势, 从而影响下一个单元乃至后续单元
   * 的密度. 于是不管是这里的策略还是在计算第 j 个单元的 Boltzmann 部分
   * 再一次把所有单元上的电势的贡献算出来, 似乎都用不上 FFT 了.
   * 
   */
}
#endif

/**
 * @brief 多线程片段
 * 
 */
// boost::thread_group threads;
// for( size_t i=0; i<n_threads-1; ++i )
// {
//   threads.create_thread(boost::bind(&SteadyState<PROBLEM>::calculate_residualHelper, this, i ));
// }

// calculate_residualHelper( n_threads-1 );
// threads.join_all();

//  boost::mutex rw_mutex[100];   /**< boost 互斥锁 */
// /// 将迭代残量写到 SteadyState::_res 中去
// boost::lock_guard<boost::mutex> lock( rw_mutex[0] );
// // mutex[0].acquire();
// _res = std::max( _res, local_res );
// // mutex[0].release();


// DEBUG 这个宏似乎可以从编译选项中传入
#ifdef DEBUG 
template <int dim>
void print_dis( const NRXX::Distribution<dim, double, 0>& _dis ){
  std::cout << "The details of dis: \n  ";
  std::cout.precision( 16 );
  std::cout << _dis << std::endl;
}

void print_jacobi_matrix( Eigen::MatrixXd& _matrix, char* _str ){
  std::ofstream os( _str );
  os.precision(16);
  for( size_t i=0; i<static_cast<size_t>(_matrix.cols()); ++i )
  {
    for(size_t j=0; j<static_cast<size_t>(_matrix.rows()); ++j )
    {
      os << _matrix(i, j) << " ";
    }
    os << std::endl;
  }
  os.close();
  std::cout << "Print Jacobi matrix is OK!" << std::endl;
}
#endif // DEBUG


#endif // __STEADY_STATE_TOOLS_H__


/**
 * end of file
 * 
 */
