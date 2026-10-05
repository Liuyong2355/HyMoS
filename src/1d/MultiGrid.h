/**
  * @file   MultiGrid.h
  * @author HU Zhicheng <huzhicheng1986@gmail.com>
  * @date   Thu Nov  7 13:41:06 2013
  * 
  * @brief  矩方程组的多重网格求解器
  * 
  * 最终想将这个算法实现结构清晰, 漂亮, 通用, 但是目前看来, 我还做不到...
  *
  * 想将网格信息从问题中剥离出来. 然后对于多重网格, 我们构建和存储多层
  * 网格的信息, 因而这个类最好继承自基础的 mesh_t ?
  * 
  * @note 还是空间一维问题
  */

#ifndef __MULTI_GRID_H__
#define __MULTI_GRID_H__

#include <memory>
#include <vector>
#include <sstream>

#include "config.h"
#include "ConfigFile.h"
#include "OperatorInterpolate.h"


/**
 * @class MultiLevelDataHandle 
 * @brief 多重网格方法中多层网格数据处理句柄类型
 *
 * 在多重网格方法中, 我们从最密的网格出发, 通过粗化得到各层网格并存储起
 * 来. 若使用的 Gummel 算法总是在给定外力下使用多重网格方法求解
 * Boltzmann 方程, 对于外力的处理, 也需要通过一个粗化插值得到各层网格上
 * 的外力信息. 这个操作应该在每次改变外力时才操作一次, 因而也要存储起来.
 * 至于碰撞项中涉及的 OneOverTau 等信息, 或许也该类似操作. 在器件模拟的
 * 算例测试表明, 固定外力下用多重网格更新 Boltzmann 方程是不好的, 对这
 * 个算例, 外力需要及时更新, 因而最新考虑的方法是, Boltzmann 方程和外力
 * 方程一起做多重网格. 在这个设计下, 外力也是要有一个粗化程序的.
 *
 * 由于这个事情只是对多重网格方法才做, 在算例中定义方法类型, 又定义多层
 * 数据类型并不恰当. 我们把这个多层数据的处理隐藏到具体的多重网格方法中,
 * 外部只要设定单层网格, 单层网格上的数据类型即可.
 *
 * 然后我们实际发现, 之前实现的多重网格几何信息类 MultiMesh1D 与将要实
 * 现多重网格外力类 MultiLevelForce 有共通之处, 应该可以通过类似的模板
 * 类生成. 
 *
 * 本类的模板参数个数可变. 当前特化了单个模板参数和双个模板参数两种情形:
 * @code
 * template <typename DATAHANDLE> class MultiLevelDataHandle<DATAHANDLE>;
 * template <typename MESH, typename DATAHANDLE> class MultiLevelDataHandle<MESH, DATAHANDLE>
 * @endcode
 *
 * 这里, 模板参数类型 DATAHANDLE 就是要在多重网格方法中应用的那些单层网
 * 格上的数据类型, 如单层网格类, 单层网格上外力类等. 这个模板类需要提供
 * 一个从细网格粗化的接口函数 CoarseFrom( DATAHANDLE ). 而模板参数MESH
 * 自然是单层的网格类型. 如果 DATAHANDLE 类提供的 CoarseFrom() 函数实现
 * 完全 (当 DATAHANDLE 附有网格信息时也包含网格的粗化), 那么用单模板参
 * 数情形封装多层数据即可. 但是外力随着计算推进会发生变化, 故而会频繁调
 * 用 CoarseFrom() 函数, 而多层网格信息相对变化缓慢 (也就是自适应等重设
 * 网格时变化). 为了效率起见, 在外力等 DATAHANDLE 的CoarseFrom() 函数中,
 * 仅粗化除网格外的其他必要成员, 至于粗化的DATAHANDLE 中的网格, 则用一
 * 个已经生成的多重网格几何信息类赋值, 而 DATAHANDLE 中那些只与网格有关
 * 的变量, 在这个网格赋值同时生成. 于是, 外力等 DATAHANDLE 的多层封装
 * 还需要一个网格参数. 并且, 这个双模板参数情形依赖于一个单模板参数情形
 * MultiLevelDataHandle<MESH>.
 * 
 * 考虑到效率问题, 这些数据成员应尽量使用指针.
 *
 * @note 
 * - 需要注意的是, 多重网格算法叙述中, 通常 T_0, T_1, \ldots, T_M是一族
 *   从粗到密的网格. 由于当前考虑算法时总是从最密的网格出发, 采用多重网
 *   格策略求解, 因而本类关于网格的存储恰好相反, 网格在 std::vector中的
 *   指标, 0 是最密的网格 T_M 上的信息, 1 是 T_{M-1}, 依次类推. 当然为
 *   了与叙述一致, 本类将隐藏这点差异, 也就是按叙述的方式返回各层网格,
 *   参见类中函数 GetDataHandle()
 * - 如果 DATAHANDLE 封装得当, 或许这里仍然可以只使用单个模板参数. 但毕
 *   竟这样分成两种情形特化比较清晰. 是否合并以后考虑.
 */
template <typename ... Args> class MultiLevelDataHandle;
template <typename DATAHANDLE> class MultiLevelDataHandle<DATAHANDLE>;
template <typename MESH, typename DATAHANDLE> class MultiLevelDataHandle<MESH, DATAHANDLE>;

/**
 * @class _MultiLevelDataHandleBase
 * @brief 多层网格数据处理句柄封装基类
 *
 * 由于 MultiLevelDataHandle 模板类提供了单个模板参数和双个模板参数两种
 * 特化, 但是获得数据处理句柄的函数是一样的, 因而用这样一个基类先封装获
 * 得数据类型句柄的相关函数.
 */
template <typename DATAHANDLE>
class _MultiLevelDataHandleBase : public std::vector<std::shared_ptr<DATAHANDLE> > {
public:
  typedef std::shared_ptr<DATAHANDLE> ptr_handle_t;
  typedef std::vector<ptr_handle_t> BASE; /**< 基类的 size 就是网格的层数 */

public:
  _MultiLevelDataHandleBase( unsigned int _n_levels = 1 ) : BASE(_n_levels) { };

  /** 
   * @brief 返回第 @p level 层的网格上信息句柄
   *
   * @param level 网格层级, 为了和多重网格方法习惯的叙述一致, level ==
   * this->size()-1 时是最密的网格上信息, level==0 时是最粗的网格, 不检
   * 查参数是否有效.
   * 
   * @return 第 level 层网格上信息句柄, 返回share_ptr引用好还是数据引用
   * 好?
   */
  ptr_handle_t& GetDataHandlePtr( const unsigned int level ) {
    return (*this)[ this->size() - 1 - level ];
  }
  const ptr_handle_t& GetDataHandlePtr( const unsigned int level ) const {
    return (*this)[ this->size() - 1 - level ];
  }

  /** 
   * 为了与习惯一致, 最密网格的层级是 this->size()-1, 这个数似乎没法作
   * 为默认参数, 而最密网格又是常用的, 于是重载下此函数
   *
   * @return 返回最密网格上信息指针
   */
  ptr_handle_t& GetDataHandlePtr() { return (*this)[0]; }
  const ptr_handle_t& GetDataHandlePtr() const { return (*this)[0]; }

  DATAHANDLE& GetDataHandle( const unsigned int level ) {
    return *((*this)[ this->size() - 1 - level ]);
  }
  const DATAHANDLE& GetDataHandle( const unsigned int level ) const {
    return *((*this)[ this->size() - 1 - level ]);
  }
  DATAHANDLE& GetDataHandle() { return *((*this)[0]); }
  const DATAHANDLE& GetDataHandle() const { return *((*this)[0]); }  
}; // class _MultiLevelDataHandleBase

/// 单个模板参数情形特化
template <typename DATAHANDLE>
class MultiLevelDataHandle<DATAHANDLE> : public _MultiLevelDataHandleBase<DATAHANDLE> {
public:
  typedef std::shared_ptr<DATAHANDLE> ptr_handle_t;
  typedef _MultiLevelDataHandleBase<DATAHANDLE> BASE;

public:
  MultiLevelDataHandle( unsigned int _n_levels = 1 ) : BASE(_n_levels) { };
  /**
   *  这个构造函数若是 const 编译要出错. 还有, 经过简单测试 shared_ptr
   *  指针的 reset 等函数, 发现当两个以上的 shared_ptr 指针用同一个
   *  DATAHANDLE& 或者 DATAHANDLE* 的数据初始化时, 会出现各种 segment
   *  fault. 所以把这两个构造函数去掉, 避免该种错误.
   *
   * @code 
   *   MultiLevelDataHandle( DATAHANDLE& _handle, 
                             unsigned int _n_levels = 1 ) : BASE(_n_levels) { 
         Reinit( ptr_handle_t( &_handle ) );
       };
       MultiLevelDataHandle( DATAHANDLE* p_handle, 
                             unsigned int _n_levels = 1 ) : BASE(_n_levels) { 
         Reinit( ptr_handle_t( p_handle ) );
       };
   * @endcode
   */

  MultiLevelDataHandle( const ptr_handle_t& p_handle, 
                        unsigned int _n_levels = 1 ) : BASE(_n_levels) {
    Reinit( p_handle );
  };

  void Reinit( const ptr_handle_t& p_handle ) {
    (*this)[0] = p_handle;
    SetMultiLevelHandle ();
  }

  void Reinit( const ptr_handle_t& p_handle, const unsigned int _n_levels ) {
    this->resize(_n_levels);
    Reinit( p_handle );
  }

  void SetMultiLevelHandle () {
    for( size_t i=1; i<this->size(); ++i )
    {
      (*this)[i].reset( new DATAHANDLE() );
      (*this)[i]->CoarseFrom( *((*this)[i-1]) );
    }
  }
  void SetMultiLevelHandle (const unsigned int _n_level ) {
    unsigned int current_levels = this->size();
    this->resize( _n_level );
    /// 若 current_levels < _n_level, 新增加的粗网格上信息由下面的方式生成填
    /// 充.
    for( size_t i=current_levels; i < this->size(); ++i )
    {
      (*this)[i].reset( new DATAHANDLE() );
      (*this)[i]->CoarseFrom( *((*this)[i-1]) );
    }
  }
}; // class MultiLevelDataHandle<DATAHANDLE>

/**
 * 两个模板参数情形特化
 * 
 * 我们最新考虑的是外力和 Boltzmann 方程一起做多重网格方法 (每层网格上
 * 的光滑步都是一个 Gummel 迭代). 这种方式下, DATAHANDLE 的
 * CoarseFrom() 将在每步多重网格内部去调用, 而在构造
 * MultiLevelDataHandle<MESH, DATAHANDLE> 类时倒是不需要调用
 * CoarseFrom() 函数了. 不过, 此类的构造还是需要多层网格的信息, 以及
 * DATAHANDLE 中那些仅与网格相关的变量. 因为这些信息在多重网格方法中是
 * 不变的, 所以还是可以先计算存起来.
 */
template <typename MESH, typename DATAHANDLE>
class MultiLevelDataHandle<MESH, DATAHANDLE> : public _MultiLevelDataHandleBase<DATAHANDLE> {
public:
  typedef std::shared_ptr<DATAHANDLE> ptr_handle_t;
  typedef _MultiLevelDataHandleBase<DATAHANDLE> BASE;
  typedef MultiLevelDataHandle<MESH> multi_mesh_t;
  typedef std::shared_ptr<multi_mesh_t> ptr_multi_mesh_t;

public:
  MultiLevelDataHandle( ) : BASE(1) { };
  MultiLevelDataHandle( const ptr_multi_mesh_t& _p_multi_mesh ) 
    : BASE(_p_multi_mesh->size()), p_multi_mesh(_p_multi_mesh) {
    SetMultiLevelHandle();
  };
  MultiLevelDataHandle( const ptr_multi_mesh_t& _p_multi_mesh, 
                        const ptr_handle_t& p_handle ) 
    : BASE(_p_multi_mesh->size()), p_multi_mesh(_p_multi_mesh) {
    SetMultiLevelHandle();
    Reinit( p_handle );
  };

  /** 
   * 这个函数也与 MultiLevelDataHandle<DATAHANDLE> 类不同了, 这里仅赋值
   * 给最密网格了. 而粗网格上的信息到具体 MultiGrid 算法中再处理.
   * 
   * @param p_handle 
   */
  void Reinit( const ptr_handle_t& p_handle ) {
    (*this)[0] = p_handle;
//    SetMultiLevelHandle ();
  }

  // void Reinit( const ptr_multi_mesh_t& _p_multi_mesh, 
  //              const ptr_handle_t& p_handle ) {
  //   p_multi_mesh = _p_multi_mesh;
  //   this->resize(p_multi_mesh->size());
  //   Reinit( p_handle );
  // }
  
  void Reinit( const ptr_multi_mesh_t& _p_multi_mesh ) {
    p_multi_mesh = _p_multi_mesh;
    this->resize( p_multi_mesh->size() );
    // 计算 DATAHANDLE 中仅与网格相关的量
    SetMultiLevelHandle();
  }
  /** 
   * 这个函数与 MultiLevelDataHandle<DATAHANDLE> 类不同, 这里只计算
   * DATAHANDLE 中那些仅与网格相关的信息, 而不调用 DATAHANDLE 的
   * CoarseFrom() 函数.
   * 
   */
  void SetMultiLevelHandle () {
    for( size_t i=0; i<this->size(); ++i )
    {
      (*this)[i].reset( new DATAHANDLE( (*p_multi_mesh)[i] ) );
    }
  }
  /** 
   * 要求 @p p_multi_mesh 和 @p _p_multi_mesh 是一致的, 即 @p
   * _p_multi_mesh 和 @p p_multi_mesh 是相互包含关系.
   * 
   * @param _p_multi_mesh 
   */
  void SetMultiLevelHandle ( const ptr_multi_mesh_t& _p_multi_mesh ) {
    p_multi_mesh = _p_multi_mesh;

    unsigned int current_levels = this->size();
    this->resize( p_multi_mesh->size() );
    /// 只处理多出来的几层网格相关信息
    for( size_t i=current_levels; i < this->size(); ++i )
    {
      (*this)[i].reset( new DATAHANDLE( (*p_multi_mesh)[i] ) );
    }
  }

private:
  ptr_multi_mesh_t p_multi_mesh; /**< 使用时需要注意其 size 与 BASE 的一样 */
}; // class MultiLevelDataHandle<MESH, DATAHANDLE>

/**
 * @class NonlinearMultiGrid
 * @brief 非线性多重网格方法类
 *
 * 此方法在粗网格上仍然求解一个非线性问题. 当基本迭代法采用 symmetric
 * local Gauss-Seidel - Newton 方法时, 只能采用此多重网格方法. 
 *
 * 原些的设计, 本类需要一个基本迭代法的模板类, 该模板类需要提供以下几个
 * 函数接口:
 *
 * - 光滑步: void Solve( sol_t& sol, const sol_t& rhs, const mesh_t&
 *                       mesh, const unsigned int steps )
 * - 粗网格上精确求解器: void ExactSolve( sol_t& sol, const sol_t& rhs,
 *                                        const mesh_t& * mesh )
 *
 * 同时, 与基本迭代法一样, 本类将提供外部调用的主要函数接口是 @code
 * void Solve( sol_t& sol, const sol_t& rhs, const mesh_t& mesh, const
 * unsigned int steps = 1 ) @endcode
 *
 * 当一个多重网格类 _MULTIMOMENTSMETHOD 同时提供上述光滑步求解器和粗网
 * 格上精确求解器时, 便可作为此类的基本迭代法而使用. 这样便实现了一种速
 * 度和位置都做 MG 的算法: 空间多重网格的每个基本迭代步都是一个矩方向的
 * 多重网格, 这正好符合最早的想法: SGS-Newton 算法的局部线性方程组用MG
 * 求解 (实际与此想法也略有不同, 局部线性方程组用 MG 是一个局部的概念,
 * 而这里实际还是将整个空间上的问题当做整体来降阶的). 但此算法当前的问
 * 题有二, 一是计算量实在太大, 二是矩数稍微大些算法就会严重退化.
 *
 * 参考 p-MG 的实现, 一个简洁些的算法是空间 MG 的光滑步仍使用单网格的迭
 * 代法, 而最粗网格上用矩方向的 MG. 为此将模板参数拆分为三个. 第一
 * 个_BASICMETHOD 通常就是一个单网格的基本迭代法, 从中提取此类的光滑步
 * 算法. 第二个 _MULTIMOMENTSMETHOD 提供最粗网格的求解器, 为保持原些接
 * 口不变, 默认与 _BASICMETHOD 相同. 考虑到纯空间的 MG 在最粗网格可能可
 * 以调用 ExactSolve, 所以引入一个 FLAG 辅助模板参数, 其值为 0 表示从
 * _BASICMETHOD 提取最粗网格求解器 (此时总是使用 ExactSolve), 值为 1 表
 * 示从_MULTIMOMENTSMETHOD 提取最粗网格求解器 (若
 * _MULTIMOMENTSMETHOD==_BASICMETHOD, 则是最粗网格上仅使用光滑步迭代数
 * 步). 当然, 按照设计原则, 当 _MULTIMOMENTSMETHOD 与 _BASICMETHOD 不同
 * 时, FLAG 一定该设为 1 (设为0 没什么意义了). 但是需要注意的是, 以前的
 * 调用方式或者偷懒的方式(后两个模板参数都省略了), 又不想在最粗网格求解
 * 时精确求解, 那就需要手动改下程序了.
 *
 * 忽然想到, 按 SingleGridMethod 类中 ExactSolve 的实现, 设置
 * input.txt 中的最粗网格迭代步数 (0 则精确求解) 即可实现上述 FLAG 想要
 * 实现的功能. 不过, 当前已经有一个更好的想法, 最终将各种 double MG 实
 * 现到统一的框架中去. 也就是说, 以后这个代码可能也会弃用, 因而就不再调
 * 整 FLAG 等接口了.
 */
template <typename _BASICMETHOD, typename _MULTIMOMENTSMETHOD=_BASICMETHOD, int FLAG=0>
class NonlinearMultiGrid {
public:
  typedef typename _BASICMETHOD::RESIDUAL RESIDUAL;
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

  typedef MultiLevelDataHandle<mesh_t> multi_mesh_t;
  typedef MultiLevelDataHandle<mesh_t, force_t> multi_force_t;

  enum { dim = dis_t::dim  };

protected:
  std::shared_ptr<multi_mesh_t> p_multi_mesh; /**< 多重网格几何信息. 现
                                               * 在 multi_mesh_t 是一个
                                               * 指针的 vector, 但为了
                                               * multi_force_t 中的使用,
                                               * 还是定义其为指针. */

  multi_force_t multi_force;    /**< 多重网格方法中外力项处理句柄 */

  unsigned int n_coarsest_smooth; /**< 最粗网格的最多迭代步数, 其值为
                                   * 0 有特殊含义, 是迭代直至收敛的意思,
                                   * @see
                                   * _BASICMETHOD::ExactSolve() */

  /**
   * 关于前后光滑步, 可能应该考虑不同层的网格光滑步数不同的策略, 这样效
   * 率会更高
   * 
   */
  unsigned int n_pre_smooth;    /**< 前光滑次数 */
  unsigned int n_post_smooth;   /**< 后光滑次数, 一般与前光滑同 */

  // unsigned int n_levels;        /**< 多重网格嵌套的最大重数, 1: 无多重,
  //                                * 2: 双网格法 */

  CycleType vw_cycle;        /**< V 循环为 1, W 循环为 2, F 循环为 3.... */
  int fmg_level_steps;       /**< FMG 中间层 startup 步数; <0 表示沿用旧行为 */
  int fmg_finest_steps;      /**< FMG 最细层 startup 步数; <0 表示不额外启动 */

  std::shared_ptr<_BASICMETHOD> p_basic_method;
  std::shared_ptr<_MULTIMOMENTSMETHOD> p_mm_method;

public:
  NonlinearMultiGrid() : p_multi_mesh(new multi_mesh_t()),
                         fmg_level_steps(-1),
                         fmg_finest_steps(-1) {};
#if 0 /**
       * 以下两个构造函数并不是与其他方法类相通的公共接口, 于是在目前的
       * 框架下基本也用不上了, 所以先不编译了.
       * 
       */
  NonlinearMultiGrid( const ptr_mesh_t& p_mesh, 
                      const unsigned int _n_levels, 
                      const unsigned int _n_pre_smooth, 
                      const unsigned int _n_post_smooth, CycleType _vw_cycle )
    : p_multi_mesh(new multi_mesh_t(p_mesh, _n_levels)),
      multi_force(multi_force_t(p_multi_mesh)),
      n_pre_smooth(_n_pre_smooth), n_post_smooth(_n_post_smooth), 
      vw_cycle(_vw_cycle), p_basic_method(new _BASICMETHOD()) { };

  NonlinearMultiGrid( const ptr_mesh_t& p_mesh, 
                      const unsigned int _n_levels, 
                      const unsigned int _n_pre_smooth, 
                      const unsigned int _n_post_smooth, CycleType _vw_cycle, 
                      const std::shared_ptr<_BASICMETHOD>& _basic_method ) 
    : p_multi_mesh(new multi_mesh_t(p_mesh, _n_levels)),
      multi_force(multi_force_t(p_multi_mesh)),
      n_pre_smooth(_n_pre_smooth), n_post_smooth(_n_post_smooth), 
      vw_cycle(_vw_cycle), p_basic_method(_basic_method) { };
#endif
  NonlinearMultiGrid( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) : p_multi_mesh(new multi_mesh_t()) {
    read_config( p_mesh, cf );
  }
  NonlinearMultiGrid( const ptr_mesh_t& p_mesh, std::istream& is ) {
    unsigned int n_levels;//p_multi_mesh->size();
    is.read((char *)&n_levels, sizeof(unsigned int));
    is.read((char *)&n_coarsest_smooth, sizeof(unsigned int));
    is.read((char *)&n_pre_smooth, sizeof(unsigned int));
    is.read((char *)&n_post_smooth, sizeof(unsigned int));
    is.read((char *)&vw_cycle, sizeof(unsigned int));
    is.read((char *)&fmg_level_steps, sizeof(int));
    is.read((char *)&fmg_finest_steps, sizeof(int));
    p_multi_mesh.reset( new multi_mesh_t(p_mesh, n_levels) );
    multi_force.Reinit(p_multi_mesh);
    p_basic_method.reset( new _BASICMETHOD( p_mesh, is ) );
    p_mm_method.reset( new _MULTIMOMENTSMETHOD( p_mesh, is ) );
  }

  void Reinit( const ptr_mesh_t& p_mesh ) {
    p_multi_mesh->Reinit( p_mesh );
    multi_force.Reinit( p_multi_mesh );
    // // 若以后单网格方法存储网格了, 就要添加下面这句了.
    // p_basic_method.Reinit( p_mesh );
  }

  void read_config( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) {
    n_coarsest_smooth = (int)cf.Value("MultiGridMethod", "CoarsestGridSmoothingSteps");
    n_pre_smooth = (int)cf.Value("MultiGridMethod", "PreSmoothingSteps");
    n_post_smooth = (int)cf.Value("MultiGridMethod", "PostSmoothingSteps");
    try {
      vw_cycle = static_cast<CycleType>((int)cf.Value("MultiGridMethod", "CycleType"));
    } catch (const std::string&) {
      vw_cycle = static_cast<CycleType>((int)cf.Value("MultiGridMethod", "VW_CYCLE"));
    }
    if (vw_cycle < VCYCLE || vw_cycle > FVCYCLE) {
      throw std::string("MultiGridMethod/CycleType must be 1 (V-cycle), 2 (W-cycle), or 3 (F-cycle)");
    }
    try {
      fmg_level_steps = std::max(0, (int)cf.Value("MultiGridMethod", "FMGLevelStartupSteps"));
    } catch (const std::string&) {
      fmg_level_steps = -1;
    }
    try {
      fmg_finest_steps = std::max(0, (int)cf.Value("MultiGridMethod", "FMGFinestStartupSteps"));
    } catch (const std::string&) {
      fmg_finest_steps = -1;
    }

    unsigned int n_levels = 0;
    try {
      const unsigned int coarsest_grid = (int)cf.Value("MultiGridMethod", "CoarsestGridSize");
      unsigned int current_grid = p_mesh->n_ele();
      if (coarsest_grid == 0 || current_grid < coarsest_grid) {
        throw std::string("MultiGridMethod/CoarsestGridSize is invalid");
      }
      n_levels = 1;
      while (current_grid > coarsest_grid) {
        if (current_grid % 2 != 0) {
          throw std::string("Mesh/Nx must stay even while coarsening toward MultiGridMethod/CoarsestGridSize");
        }
        current_grid /= 2;
        ++n_levels;
      }
      if (current_grid != coarsest_grid) {
        throw std::string("Mesh/Nx and MultiGridMethod/CoarsestGridSize are incompatible");
      }
    } catch (const std::string&) {
      n_levels = (int)cf.Value("MultiGridMethod", "MeshLevels");
    }
    p_multi_mesh->Reinit( p_mesh, n_levels );
    multi_force.Reinit(p_multi_mesh);
    if( p_basic_method.get()==NULL )    
      p_basic_method.reset( new _BASICMETHOD( p_mesh, cf ) );
    else 
      p_basic_method->read_config( p_mesh, cf );
    
    if( p_mm_method.get() == NULL )
      p_mm_method.reset( new _MULTIMOMENTSMETHOD( p_mesh, cf ) );
    else
      p_mm_method->read_config( p_mesh, cf );
  }
  void print_config() const {
    std::cout << "NonlinearMultiGrid: n_coarsest_smooth = " << n_coarsest_smooth << std::endl;
    std::cout << "NonlinearMultiGrid: n_pre_smooth = " << n_pre_smooth << std::endl;
    std::cout << "NonlinearMultiGrid: n_post_smooth = " << n_post_smooth << std::endl;
    std::cout << "NonlinearMultiGrid: multi_levels = " << p_multi_mesh->size() << std::endl;
    std::cout << "NonlinearMultiGrid: CycleType (V-, W-, FV-cycle, ...: 1, 2, 3, ...) = " << vw_cycle << std::endl;
    std::cout << "FullNonlinearMultiGrid: LevelStartupSteps = " << fmg_level_steps << std::endl;
    std::cout << "FullNonlinearMultiGrid: FinestStartupSteps = " << fmg_finest_steps << std::endl;
    p_basic_method->print_config();
    std::cout << "IfDoubleMultiGrid or ExactSolveCoarestProblem: FLAG = " << FLAG << std::endl;
    p_mm_method->print_config();
  }

  void dump_data( std::ostream& os ) const { 
    unsigned int n_levels = p_multi_mesh->size();
    os.write((const char *)&n_levels, sizeof(unsigned int));
    os.write((const char *)&n_coarsest_smooth, sizeof(unsigned int));
    os.write((const char *)&n_pre_smooth, sizeof(unsigned int));
    os.write((const char *)&n_post_smooth, sizeof(unsigned int));
    os.write((const char *)&vw_cycle, sizeof(unsigned int) );
    os.write((const char *)&fmg_level_steps, sizeof(int));
    os.write((const char *)&fmg_finest_steps, sizeof(int));

    p_basic_method->dump_data( os );
    p_mm_method->dump_data( os );
  }

  void SetParameters( const unsigned int _n_levels,
                      const unsigned int _n_pre_smooth, 
                      const unsigned int _n_post_smooth, 
                      const CycleType _vw_cycle ) {
    n_pre_smooth = _n_pre_smooth;
    n_post_smooth = _n_post_smooth;
    vw_cycle = _vw_cycle;
    reset_levels( _n_levels );
  }

  void reset_levels( const unsigned int levels ) {
    p_multi_mesh->SetMultiLevelHandle( levels );
    multi_force.SetMultiLevelHandle( p_multi_mesh );
  }

  void reset_tol( const double _tol ) {
    p_basic_method->reset_tol( _tol );
    p_mm_method->reset_tol( _tol );
  }
  void reset_smoothing( const unsigned int _n_smooth ) {
    n_pre_smooth = _n_smooth;
    n_post_smooth = _n_smooth;
  }

  /** 
   * @brief 多重网格法求解的函数接口
   * 
   * 由于基本迭代法没有存解和网格的信息, 因而接口函数需要传入网格信息.
   * 而在这里, 网格的信息显然有些多余, 也就是这个接口实际还不能很好的统
   * 一起来.
   *
   * @param sol 
   * @param rhs 
   * @param p_mesh 无用信息, 或者不考虑接口的统一性而去掉此参数?
   * @param force_handle 外力项处理句柄, 提供 @p mesh 上的外力等信息
   * @param steps 
   */
  void Solve( sol_t& sol, const sol_t& rhs, const ptr_mesh_t& p_mesh,
              const ptr_force_t& p_force_handle,
              const unsigned int steps = 1 ) {
    /**
     * 由于之前考虑了固定外力下, Boltzmann 方程的多重网格方法, 于是代码
     * 的遗留, 这个接口不是很好, 为了最密网格上的外力与 SteadyState 类
     * 中的外力一致, 下面这句话还是要的.
     * 
     */
    multi_force.Reinit( p_force_handle );

    for( size_t i=0; i<steps; ++i )
      MultiGridOneStep( sol, rhs, p_multi_mesh->size()-1 );
  }
  /** 
   * @brief 后处理: 质量守恒校正等
   * 
   * @param sol 
   * @param p_mesh 为了接口和单网格方法统一, 这里先加上这个参数
   */
  void PostProcessing( sol_t& sol, const ptr_mesh_t& p_mesh ) {
    // 后处理只是对最细网格上的解做的, 传入的 sol 正是 mesh 上的解. 也
    // 就是说传入的 p_mesh 应该和 p_multi_mesh->GetDataHandlePtr() 是一致
    // 的. 当然粗网格上也可以做这个后处理, 只是没必要.
    p_basic_method->PostProcessing( sol, p_mesh );
  }

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

protected:
  /** 
   * @brief 成员模板函数的重载版本
   *
   * 模板函数的统一调用接口. 可能最好是把循环类型在编译时给定, 但目前的
   * 代码是从配置文件 (input.txt) 中读取具体的循环类型, 存在 vw_cycle
   * 变量中. 这样就需要程序运行时用 switch 选择循环类型, 所幸收敛的多重
   * 网格迭代步数通常不会太多, 把 switch 放在最外层迭代中判断, 也不算费
   * 时.
   *
   * 还有, 这个函数仅是当前最密网格上的调用 (这不是一个递归函数, 成员模
   * 板函数才是递归调用版本)
   *
   * @param _sol 
   * @param _rhs 
   * @param _level 
   */
  void MultiGridOneStep( sol_t &_sol, const sol_t &_rhs, 
                         const unsigned int _level ){
    switch( vw_cycle ) {
    case VCYCLE:
      MultiGridOneStep<VCYCLE>( _sol, _rhs, _level );
      break;
    case WCYCLE:
      MultiGridOneStep<WCYCLE>( _sol, _rhs, _level );
      break;
    case FVCYCLE:
      MultiGridOneStep<FVCYCLE>( _sol, _rhs, _level );
      break;
    default:
      throw std::string("MultiGridMethod/CycleType must be 1 (V-cycle), 2 (W-cycle), or 3 (F-cycle)");
    }
  }
  
  /** 
   * @brief 一步多重网格迭代
   * 
   * 模板参数用来决定粗网格校正时所使用的循环类型. 由于 FV 型等循环方式
   * 会用到 V 型循环, 因而这个这个模板参数不大适合作为整个类的模板参数.
   * 由于模板参数应该在编译期给定, 而具体循环类型是从配置文件
   * (input.txt) 中读取, 因此重载了一个同名函数作为调用的统一接口.
   * 
   * 此函数可以设为 private, 因为派生的 FullNonlinearMultiGrid 类也不直
   * 接调用它了.
   * 
   * @param _sol 
   * @param _rhs 
   * @param _level 
   */
  template <CycleType _cycle>
  void MultiGridOneStep( sol_t &_sol, const sol_t &_rhs, 
                         const unsigned int _level ){
    if( _level == 0 )
    {/**
      * 最粗网格上也不用前后光滑了, 所以把这递归调用终止条件放到最前面.
      * 
      */
      CoarsestSolver( _sol, _rhs );
      return;
    }
      
    
    /// pre-smooth
    Smoothing( _sol, _rhs, _level, n_pre_smooth );
    
    /**
     * sub grid
     * 
     */
    {
      sol_t finer_res( _sol.size( ) ); // 细网格上的残量
      RESIDUAL::GetResidual( finer_res, _sol, _rhs, p_multi_mesh->GetDataHandle(_level), multi_force.GetDataHandle(_level) );
      /**
       * 可以观察到几步以后在最粗网格上只要迭代一步就能得到精确解, 或者
       * 说可以认为实际这个时候粗网格上都可以不用求解了? 于是这里顺带加
       * 入这个细网格残量大小的判断, 若细网格上的已经是精确解, 那么就不
       * 再生事去粗网格校正了. 
       *
       * 可是若最粗網格仅仅求解到最密网格的 tolerance, 有时候会使得最密
       * 网格的 residual 降不下来. 所以这里要判断的话, 选用的tolerance
       * 需要比最密网格再小些才行, 但是未必就有效率的说. (还有, 残量计
       * 算效率似乎也不算太高.)
       * 
       */

      ptr_mesh_t& coarser_mesh = p_multi_mesh->GetDataHandlePtr(_level-1);
      std::unique_ptr<sol_t> coarser_sol( new sol_t( coarser_mesh->n_ele() ) );
      /**
       * 近似解的投影 (限制), 细网格解变量的展开参数系变为与粗网格对应
       * 单元相同
       * 
       */
      Restriction( *coarser_sol, _sol, _level );
      // 保存一份粗网格初值以做回代用.
      sol_t init_coarser_sol=*coarser_sol; 

      /// 细网格残量的投影, 使用粗网格上解的参数系表示, 注意 finer_res 也被修改
      sol_t coarser_rhs( coarser_mesh->n_ele() );
      Restriction( coarser_rhs, finer_res, _level, *coarser_sol );

      // 外力的限制
      RestrictionForce( _level );
      // 保存粗网格上的初始外力
      force_t old_coarser_force = multi_force.GetDataHandle( _level-1 );
      // 外力方程右端项
      CoarseForceRHS( *coarser_sol, _level );

      /// 粗网格上方程左边的残量 - R(f_H)
      sol_t coarser_res( coarser_mesh->n_ele() );
      RESIDUAL::GetResidual( coarser_res, *coarser_sol, *coarser_mesh, multi_force.GetDataHandle(_level-1) );

      /// 最终粗网格上的右端项 coarser_rhs + R(f_H)
      coarser_rhs -= coarser_res;

      /**
       * 粗网格上求解, 但是 FV 等循环不太适合直接的递归调用, 因而重新封
       * 装一个函数.
       * 
       */
      SubGridIteration<_cycle>( *coarser_sol, coarser_rhs, _level-1 );

      // 计算粗网格的变化量, 并存储在 @p *coarser_sol 中
      CalcCorrection( init_coarser_sol, *coarser_sol );
      /// 回代, _sol += mu * correction, 步长 mu 用来防止这一步出现 nan
      BackSubstitution( *coarser_sol, _sol );

      /**
       *  外力的回代. Poisson 方程附着与 Boltzmann 方程, 所以总是可以直
       *  接计算出来. 若按多重网格的延拓算子来回代, 则当单网格方法使用
       *  SLGSNewton<RESIDUAL> 时, 多重网格就基本会失败. 而若这里外力直
       *  接由 Poisson 方程计算, SLGSNewton<RESIDUAL> 对某些情形, 还是
       *  能够用多重网格方法计算下来的.
       *
       * 由于这个外力并不是通过迭代计算得来, 所以这两种方式仅仅会影响其
       * 中一步 Boltzmann 方程求解的不同. 看来还是
       * SLGSNewton<RESIDUAL> 太敏感了.
       * 
       */
//      BackSubstitutionForce( old_coarser_force, _level );
      multi_force.GetDataHandlePtr( _level )->CalcExternalForce( _sol );
    }

    /// post-smooth
    Smoothing( _sol, _rhs, _level, n_post_smooth );    
  }
  
  /** 
   * @brief 粗网格上的求解
   * 
   * FV 等循环类型的多重网格不太适合直接的递归调用, 因此将粗网格的求解
   * 封装出来, 成为一个独立的模板函数. 这样即可实现各个 type 的循环, 重
   * 复代码也尽可能不增加太多. 注意, 虽说这是一个模板函数, 但是除了 V
   * 型与 W 型, 很难统一写在一起. 也就是说这个函数实际需要进行特化.
   *
   * 问题是 g++ 不允许在类定义内部特化成员模板函数. 另一方面, 即使再类
   * 外部特化其成员模板函数, 其前提要求也是类本身是一个特定的类型. 也就
   * 是说对于类模板而言要先进行具体的特化, 才能进而特化其成员模板函数.
   * 这样显然不大合实际, 我们需要的仅仅是特化这一成员函数, 而不需特化整
   * 个类.
   *
   * 后来找到一个替代模板函数的重载函数方式. 就是定义一个结构体 (见
   * config.h), 可以将 CycleType 的枚举值转换一个个独立的数据类型. 这样
   * 就可以重载函数, 使得该函数参数个数相同, 参数类型仅差一个枚举值. 这
   * 样调用的时候仍然有一个统一的接口. 这个函数暂时名为
   *
   * void SubGridIteration( sol_t&, const sol_t&, const unsigned int,
   * CycleTypeIdentity<CycleType> );
   *
   * 实际上有了这样一个重载函数, 本模板函数倒是可以不需要了. 这里也仅仅
   * 是直接调用重载函数而已.
   *
   * 本函数以及重载的同名函数都更应该属于 private.
   * 
   * @param _sol 
   * @param _rhs 
   * @param _level 
   */
  template <CycleType _cycle>
  void SubGridIteration( sol_t &_sol, const sol_t &_rhs, 
                         const unsigned int _level ){
    SubGridIteration( _sol, _rhs, _level, CycleTypeIdentity<_cycle>() );
  }
  /** 
   * @brief V 型多重网格方法的粗网格求解
   *
   * 直接的递归调用自己而已.
   */
  void SubGridIteration( sol_t &_sol, const sol_t &_rhs, 
                         const unsigned int _level, CycleTypeIdentity<VCYCLE> ){
    MultiGridOneStep<VCYCLE>( _sol, _rhs, _level );
  }
  /** 
   * @brief W 型多重网格方法的粗网格求解
   * 
   * 递归调用自己两次.
   */
  void SubGridIteration( sol_t &_sol, const sol_t &_rhs, 
                         const unsigned int _level, CycleTypeIdentity<WCYCLE> ) {
    MultiGridOneStep<WCYCLE>( _sol, _rhs, _level );
    MultiGridOneStep<WCYCLE>( _sol, _rhs, _level );
  }
  /** 
   * @brief FV 型多重网格方法的粗网格求解
   * 
   * 调用两次: 第一次调用自己 (递归), 第二次调用 V 型迭代.
   */
  void SubGridIteration( sol_t &_sol, const sol_t &_rhs, 
                         const unsigned int _level, CycleTypeIdentity<FVCYCLE> ) {
    MultiGridOneStep<FVCYCLE>( _sol, _rhs, _level );
    MultiGridOneStep<VCYCLE>( _sol, _rhs, _level );
  }
  /** 
   * @brief VF 型多重网格方法的粗网格求解
   * 
   * 调用两次: 第一次调用 V 型迭代, 第二次调用自己 (递归).
   */
  void SubGridIteration( sol_t &_sol, const sol_t &_rhs, 
                         const unsigned int _level, CycleTypeIdentity<VFCYCLE> ) {
    MultiGridOneStep<VCYCLE>( _sol, _rhs, _level );
    MultiGridOneStep<VFCYCLE>( _sol, _rhs, _level );
  }
  /** 
   * @brief VV 型多重网格方法的粗网格求解
   * 
   * 调用两次 V 型迭代.
   */
  void SubGridIteration( sol_t &_sol, const sol_t &_rhs, 
                         const unsigned int _level, CycleTypeIdentity<VVCYCLE> ) {
    MultiGridOneStep<VCYCLE>( _sol, _rhs, _level );
    MultiGridOneStep<VCYCLE>( _sol, _rhs, _level );
  }
  /** 
   * @brief VFFV 型多重网格方法的粗网格求解
   * 
   * 第一次和第三次调用 V 型迭代, 第二次调用自己 (递归). 画一个示意图可
   * 以看到就是 VF 和 FV 拼起来的
   */
  void SubGridIteration( sol_t &_sol, const sol_t &_rhs, 
                         const unsigned int _level, CycleTypeIdentity<VFFVCYCLE> ) {
    MultiGridOneStep<VCYCLE>( _sol, _rhs, _level );
    MultiGridOneStep<VFFVCYCLE>( _sol, _rhs, _level );
    MultiGridOneStep<VCYCLE>( _sol, _rhs, _level );
  }
  /** 
   * @brief VFVV 型多重网格方法的粗网格求解
   * 
   * 第一次和第三次调用 V 型迭代, 第二次调用FV型.
   */
  void SubGridIteration( sol_t &_sol, const sol_t &_rhs, 
                         const unsigned int _level, CycleTypeIdentity<VFVVCYCLE> ) {
    MultiGridOneStep<VCYCLE>( _sol, _rhs, _level );
    MultiGridOneStep<FVCYCLE>( _sol, _rhs, _level );
    MultiGridOneStep<VCYCLE>( _sol, _rhs, _level );
  }
  /** 
   * @brief VVFV 型多重网格方法的粗网格求解
   * 
   * 第一次和第三次调用 V 型迭代, 第二次调用VF型.
   */
  void SubGridIteration( sol_t &_sol, const sol_t &_rhs, 
                         const unsigned int _level, CycleTypeIdentity<VVFVCYCLE> ) {
    MultiGridOneStep<VCYCLE>( _sol, _rhs, _level );
    MultiGridOneStep<VFCYCLE>( _sol, _rhs, _level );
    MultiGridOneStep<VCYCLE>( _sol, _rhs, _level );
  }


  /** 
   * @brief 前后光滑步迭代
   * 
   * @param _sol 给定迭代初值, 计算后返回值
   * @param _rhs 多重网格求解方程右端项, 最细网格上是 0
   * @param _level 当前求解方程的网格层级, 最密的网格是 \f$ M\f$
   * @param n_smooth 光滑的次数, 1次是一个 SGS 迭代步
   */
  void Smoothing( sol_t & _sol, const sol_t & _rhs,
                  const unsigned int _level, const unsigned int n_smooth ){
    /**
     * 怎样做这个接口竟难煞成也
     * 
     */
    p_basic_method->Solve( _sol, _rhs, p_multi_mesh->GetDataHandlePtr(_level), multi_force.GetDataHandlePtr(_level), n_smooth );
  }

  /** 
   * @brief 最粗网格上的求解器
   *
   * @param _sol 
   * @param _rhs 
   */
  void CoarsestSolver( sol_t &_sol, const sol_t &_rhs ) {
    /**
     * 对最粗网格, 也只能通过迭代 (Smoothing) 得到精确解. 当网格单元数
     * 还是比较大的时候, 这样的计算量还是很可观的. 某种程度来讲, 最粗网
     * 格也没必要精确求解. 因而当网格单元数较大时, 这里的求解器不妨换成
     * 普通的 Smoothing, 仅仅是比细网格上多几步 Smoothing 就可以了.
     * 
     */
    CoarsestSolver( _sol, _rhs, ArtificialTypeIdentity<FLAG>() );
  }
  void CoarsestSolver( sol_t &_sol, const sol_t &_rhs, ArtificialTypeIdentity<0> ) {
    /**
     *  按照 SingleGridMethod 中的代码, 如果 n_coarsest_smooth == 0, 则
     *  是精确求解, 若 n_coarsest_smooth > 0, 那么最多还是仅迭代
     *  n_coarsest_smooth 步.
     * 
     */
    p_basic_method->ExactSolve( _sol, _rhs, p_multi_mesh->GetDataHandlePtr(0), multi_force.GetDataHandlePtr(0), n_coarsest_smooth );
  }
  void CoarsestSolver( sol_t &_sol, const sol_t &_rhs, ArtificialTypeIdentity<1> ) {
    p_mm_method->Solve( _sol, _rhs, p_multi_mesh->GetDataHandlePtr(0), multi_force.GetDataHandlePtr(0), n_coarsest_smooth );
  }

  /** 
   * @brief 解或残量从细网格到粗网格的限制
   * 
   * 为了少回代少些投影的操作. 返回的细网格上的解也将变为在粗网格相应单
   * 元参数系下的展开.
   *
   * @see essay/multigrid.tex and doc/Multigrid1D
   *
   * 该文档中的方式实现虽然简洁些, 但没有考虑回代可能出现的 nan, 所以这
   * 里还是按照最原始的方式来, 即: 计算粗网格的展开参数系, 将细网格的分布函
   * 数通过同伦投影在相应粗网格参数系下展开 (返回). 获得粗网格上的解
   * (返回).
   *
   * @param coarser_sol 粗网格上的解
   * @param finer_sol 细网格上传入解, 返回在粗网格参数系下展开的解
   * @param _level 细网格在多重网格中的层级, 所以取值是 \f$ M,M-1,\ldots,1\f$
   */
  void Restriction( sol_t &coarser_sol, sol_t &finer_sol, 
                    const unsigned int _level ) {
    const ptr_mesh_t& coarser_mesh = p_multi_mesh->GetDataHandlePtr( _level-1 );
    const ptr_mesh_t& finer_mesh = p_multi_mesh->GetDataHandlePtr( _level );

    Operator::Restriction( coarser_sol, *coarser_mesh, finer_sol, *finer_mesh );
  }
  
  /** 
   * @brief 使用粗网格固有参数系, 获得细网格上残量在粗网格上的限制
   *
   * 注意这里会改变 @p finer_res 
   * 
   * @param coarser_rhs 粗网格残量, 展开参数系由 @p coarser_sol
   * @param finer_res 细网格残量
   * @param _level 细网格在多重网格中的层级
   * @param coarser_sol 粗网格上的解, 提供粗网格残量展开参数系
   */
  void Restriction( sol_t& coarser_rhs, sol_t& finer_res, 
                    const unsigned int _level, 
                    const sol_t& coarser_sol ) {
    const ptr_mesh_t& coarser_mesh = p_multi_mesh->GetDataHandlePtr( _level-1 );
    const ptr_mesh_t& finer_mesh = p_multi_mesh->GetDataHandlePtr( _level );

    for( size_t i=0; i<coarser_mesh->n_ele(); ++i ) 
    {
      Operator::Restriction( coarser_rhs[i], coarser_mesh->Length(i),
                             finer_res[2*i], finer_mesh->Length(2*i),
                             finer_res[2*i+1], finer_mesh->Length(2*i+1),
                             coarser_sol[i].Center(), coarser_sol[i].Scaling() );
    }
    
  }

  /** 
   * @brief 回代得到细网格的解
   * 
   * 公式: finer_sol += mu * (new_coarser_sol - old_coarser_sol). 步长
   * mu 默认为 1, 引入是为了防止这一步出现 nan. 暂时把这个 mu 放到局部
   * 单元控制, 不知会不会引入不稳定因素.
   *
   * @note 本来还在想, 对这个校正量做一个重构回到细网格上再校正, 可能会
   * 更好. 但目前并没有好算法, 在 Operator::Interp 若是重构插值的话, 要
   * 求传入的解是在标准展开下的, 于是就需要对这个校正量做一个标准投影.
   * 而计算显示, 对这个校正量做标准投影几乎立马出 nan, 所以还是先采用恒
   * 等的延拓算子吧
   *
   * @param coarser_sol 其值为 new_coarser_sol - old_coarser_sol, 展开
   *                    参数系痛 old_coarser_sol, 于是也与 finer_sol 对
   *                    应单元相同.
   * @param finer_sol 返回新的解
   */
  void BackSubstitution( const sol_t& coarser_sol, sol_t& finer_sol ){
    for( size_t i=0; i<coarser_sol.size(); ++i )
    {
      double mu = 1.;
      BackSubstitution( coarser_sol[i], finer_sol[2*i], finer_sol[2*i+1], mu );
    }
  }
  /** 
   * @brief 回代
   *
   * 按算法, 此函数的初值 dis_h_l, dis_h_r 以及 dis_H 有相同的展开参数
   * 系 u_H, theta_H. 计算完成后, dis_h_l, dis_h_r 在各自的标准展开参数
   * 系下.
   * 
   * @param dis_H 
   * @param dis_h_l 
   * @param dis_h_r 
   * @param mu 默认更新步长
   */
  void BackSubstitution( const dis_t& dis_H, dis_t& dis_h_l, dis_t& dis_h_r, double& mu ) {
    _PreservingPositivityStepLength<dis_t, dis_t, PROBLEM>( mu, dis_h_l, dis_H );
    _PreservingPositivityStepLength<dis_t, dis_t, PROBLEM>( mu, dis_h_r, dis_H );
    
    dis_h_l.Add(mu, dis_H);
    dis_h_r.Add(mu, dis_H);

    dis_t tmp_l( dis_H.Center(), dis_H.Scaling(), dis_H.GetOrder() );
    dis_t tmp_r( dis_H.Center(), dis_H.Scaling(), dis_H.GetOrder() );
    tmp_l.ProjectToStdSpace( dis_h_l );
    tmp_r.ProjectToStdSpace( dis_h_r );
    std::swap( tmp_l, dis_h_l );
    std::swap( tmp_r, dis_h_r );
  }

  /** 
   * @brief 计算粗网格上的校正量
   * 
   * @param old_sol 粗网格初值
   * @param new_sol 粗网格新的解, 返回 new_sol - old_sol, 但是展开参数
   *                系是 old_sol
   */
  void CalcCorrection( const sol_t& old_sol, sol_t& new_sol ) {
    for( unsigned int i=0; i<new_sol.size(); ++i ) {
      dis_t delta_dis( old_sol[i].Center(), old_sol[i].Scaling(), old_sol[i].GetOrder() );
      Project( new_sol[i], delta_dis );
      delta_dis -= old_sol[i];
      std::swap( new_sol[i], delta_dis );
    }
  }

  /** 
   * @brief 外力的限制
   * 
   * 当前的外力接口也不算很好, 这个限制就是调用 force_t 的 CoarserFrom
   * 函数, 同时将其中的变量以及方程残量做限制到粗网格上.
   *
   * @param _level 
   */
  void RestrictionForce( const unsigned int _level ) {
    multi_force.GetDataHandlePtr( _level-1 ) -> CoarseFrom( multi_force.GetDataHandle( _level ) );    
  }
  /** 
   * @brief 粗网格上外力方程的右端项
   * 
   * @param _sol 
   * @param _level 注意是细网格的层级
   */
  void CoarseForceRHS( const sol_t& _sol, const unsigned int _level ) {
    multi_force.GetDataHandlePtr( _level-1 ) -> CoarseForceEquationRHS( _sol );
  }
  /** 
   * @brief 外力的回代
   * 
   * @param old_coarser_force 粗网格上旧的外力
   * @param _level
   */
  void BackSubstitutionForce( const force_t& old_coarser_force, 
                              const unsigned int _level ) {
    // 粗网格上校正量, 存储在自己的某些变量中
    multi_force.GetDataHandlePtr( _level-1 ) -> CalcCorrection( old_coarser_force );
    multi_force.GetDataHandlePtr( _level ) -> BackSubstitution( multi_force.GetDataHandle(_level - 1) );
  }
}; // class NonlinearMultiGrid

/**
 * @class FullNonlinearMultiGrid
 * @brief 非线性多重网格方法的 Full MG 实现
 * 
 * 顾名思义, 可知此类的方法. 与 NonlinearMultiGrid 方法相比, Full MG 只
 * 是准备一个好的初值, 因而此类继承自 NonlinearMultiGrid 类, 同时将仅修
 * 改初值准备函数 initial_value();
 */
template <typename _BASICMETHOD>
class FullNonlinearMultiGrid : public NonlinearMultiGrid<_BASICMETHOD> {
public:
  typedef NonlinearMultiGrid<_BASICMETHOD> BASE;

  typedef typename _BASICMETHOD::RESIDUAL RESIDUAL;
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

  typedef MultiLevelDataHandle<mesh_t> multi_mesh_t;
  typedef MultiLevelDataHandle<mesh_t, force_t> multi_force_t;

  //enum { dim = dis_t::dim  };

public:
  FullNonlinearMultiGrid() : BASE() {};
#if 0
  FullNonlinearMultiGrid( const ptr_mesh_t& p_mesh, 
                          const unsigned int _n_levels, 
                          const unsigned int _n_pre_smooth, 
                          const unsigned int _n_post_smooth, CycleType _vw_cycle )
    : BASE(p_mesh, _n_levels, _n_pre_smooth, _n_post_smooth, _vw_cycle ) { };

  FullNonlinearMultiGrid( const ptr_mesh_t& p_mesh, 
                          const unsigned int _n_levels, 
                          const unsigned int _n_pre_smooth, 
                          const unsigned int _n_post_smooth, CycleType _vw_cycle, 
                          const std::shared_ptr<_BASICMETHOD>& _basic_method ) 
    : BASE( p_mesh, _n_levels, _n_pre_smooth, _n_post_smooth, _vw_cycle, _basic_method ) {};
#endif
  FullNonlinearMultiGrid( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) : BASE(p_mesh, cf) { };
  FullNonlinearMultiGrid( const ptr_mesh_t& p_mesh, std::istream& is ) : BASE(p_mesh, is) { };

public:
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
    std::cerr << "FullNonlinearMultiGrid Solver: Prepare initial value ..." << std::endl;

    std::ofstream os_res( "res_step", std::ios::app );
    os_res.precision( 12 );

    os_res << "% FullNonlinearMultiGrid Solver: Prepare initial value ..." << std::endl;
    
    timeval start_time;
    gettimeofday(&start_time, NULL);
  
    (this->multi_force).Reinit( p_force_handle );
    ptr_mesh_t _p_mesh = this->p_multi_mesh->GetDataHandlePtr( 0 );
    const unsigned int ORDER = sol[0].GetOrder();
    sol_t _sol( _p_mesh->n_ele(), ORDER ), _rhs;
    PROBLEM::initial_value( _sol, *_p_mesh );

    for( unsigned int _level = 0; _level < this->p_multi_mesh->size()-1; ++_level ){
      os_res << "% level " << _level << ":" <<std::endl;

      ptr_force_t& _p_force_handle = (this->multi_force).GetDataHandlePtr( _level );
      _p_force_handle->initial_value( _sol );
      _rhs.Reinit( _sol );
      const unsigned int startup_steps = (this->fmg_level_steps < 0)
          ? 0
          : static_cast<unsigned int>(this->fmg_level_steps);
      SolveAtLevel( _sol, _rhs, _level, startup_steps );

#if 1 // 保存各层的最优解数据, 以免要用时重新算
    std::stringstream file_sol;
    file_sol << "Sol" << 1000+_level << "th.dat"; // 1000 + _level 以和
                                                  // 原来的 n_run 表示
                                                  // 各次解区别.
    SaveMacroVars( _sol, *_p_mesh, file_sol.str() );
    
    std::stringstream file_dis;
    file_dis << "DisSol" << 1000+_level << "th.dat";
    SaveSolution( _sol, *_p_mesh, file_dis.str() );

    _p_force_handle->save_results( 1000+_level );
    PROBLEM::save_current( _sol, *_p_mesh, 1000+_level );
//    _BV::save_one_over_tau( *sol, *p_mesh, *p_force_handle, n_run);
#endif

      timeval level_time;
      gettimeofday( &level_time, NULL );
      os_res << "% Elapsed time for level " << _level << " is: " << getElapsedTime( start_time, level_time ) << "s\n" <<std::endl;
      start_time = level_time;
     
      ptr_mesh_t p_next_mesh = this->p_multi_mesh->GetDataHandlePtr(_level+1);
      sol_t next_level_sol( p_next_mesh->n_ele(), ORDER );
      Operator::Interp<mesh_t,sol_t,PROBLEM>( *_p_mesh, _sol, *p_next_mesh, next_level_sol );
      std::swap(_sol, next_level_sol); 
      _p_mesh = p_next_mesh;
    }

    if (this->fmg_finest_steps >= 0) {
      const unsigned int finest_level = this->p_multi_mesh->size() - 1;
      ptr_force_t& p_finest_force = (this->multi_force).GetDataHandlePtr(finest_level);
      p_finest_force->initial_value(_sol);
      _rhs.Reinit(_sol);
      SolveAtLevel(_sol, _rhs, finest_level, static_cast<unsigned int>(this->fmg_finest_steps));
    }

    std::swap( _sol, sol );
//  this->multi_force.Reinit( p_force_handle );

    os_res.close(); // 不知这里使用会否对其他函数操作相同 os_res 产生影响.
    std::cerr << "FullNonlinearMultiGrid Solver: Prepare initial value is OK!" << std::endl;
  }

private:
  /** 
   * @brief 多重网格方法类的精确求解函数
   *
   * 本质上等价于 SingleGridMethod 中的 ExactSolve() 函数, 也等同于
   * SteadyState 类中的最外层循环 run(). 
   *
   * 前两个传入参数分别是 Boltzmann 方程的解和右端项. @p _level 表示做
   * 一个 _level 层网格求解, 当 @p _level = 0 时, 是一层网格, 即做一个
   * SingleGridMethod 求解. @p max_steps 表示求解的最大次数, 若其值为默
   * 认的 0, 则一直求解到收敛 (残量满足收敛准则).
   *
   * 因为这个函数主要是为 Full mg 准备好初值而设的, 从经济的瀑布型多重
   * 网格可以看到, 各层 max_steps 一个合适的选择可以大大改进效率. 当然
   * 在 Full mg 中, 每层迭代的步数不会如瀑布型那么多, 所以经济的瀑布多
   * 重网格的策略可以以后再行考虑. 
   *
   * 另外一个是, 在准备初值的计算中, 目前基本考虑每层都求解到稳态, 这当
   * 然也不错, 一次性就把各网格上的稳态解都得到了. 但是如果仅考虑最密网
   * 格上的解, 或许可以放宽准备初值时其余网格方法的精度. 我们已经知道,
   * 固定的 max_steps 对于 full mg 不是一个好的选择, 所以可以考虑加一个
   * 参数 @p _tol, 使得准备初值时, 只要残量小于 _tol 即停止计算, 将返回
   * 值作为下一层网格的初值. 这个也先留着之后考虑.
   *
   * @note 这个函数名不能取为 Solve, 否则将覆盖掉继承自基类的同名函数.
   * (这个不能算重载么?)
   * 
   * @param sol 
   * @param rhs 
   * @param _level
   * @param max_steps 
   */
  void SolveAtLevel( sol_t& sol, const sol_t& rhs, const unsigned int _level, 
                     const unsigned int max_steps=0 ) {
    std::cerr << "Solve by FullNonlinearMultiGrid with level " << _level << " ..." << std::endl;
    const ptr_mesh_t& p_current_mesh = this->p_multi_mesh->GetDataHandlePtr(_level);
    const ptr_force_t& p_current_force_handle = (this->multi_force).GetDataHandlePtr(_level);
    double res;
    RESIDUAL::GetResidualPL2Norm( res, sol, rhs, *p_current_mesh, *p_current_force_handle );
    unsigned int step = 0;
    
    std::ofstream os_res( "res_step", std::ios::app );
    os_res.precision( 12 );
    os_res << step << " " << res << "\n";
    const bool use_legacy_exact = (max_steps == 0 && this->fmg_level_steps < 0 && this->fmg_finest_steps < 0);
    if (!use_legacy_exact && max_steps == 0) {
      std::cerr << "stop after " << step << " steps with residual: " << res << std::endl;
      os_res.close();
      return;
    }
    do {
      ++step;
      this->MultiGridOneStep( sol, rhs, _level );
      this->PostProcessing( sol, p_current_mesh );
      RESIDUAL::GetResidualPL2Norm( res, sol, rhs, *p_current_mesh, *p_current_force_handle );
      os_res << step << " " << res << std::endl;
      std::cout << step << " " << res << std::endl;

      /**
       * 准备的初值并不需要一定是粗网格的稳态解, 所以为了效率可以放宽粗
       * 网格的迭代终止条件, 但若放的太宽, 又会使得最密网格的迭代步数增
       * 加 (于是反而更慢), 所以具体的选择要小心些. 如,
       * @code
       * while ((res >=(this->p_basic_method->reset_tol()*1e3) || step<6) && (step<max_steps || max_steps==0));
       * @endcode
       * 以上的语句要结合 max_steps > 0 使用, 对测试的 15360 网格上的器
       * 件, 从 126s 减少到 88s, 但更密的网格这个判断又未必合适了.
      */
      getControl();
    } while (res >=(this->p_basic_method->reset_tol()) &&
             (use_legacy_exact ? true : (step < max_steps)));

    /**
     * 计算完毕, 可以考虑把数据保存下来
     * 
     */

    // _level=0的话, step 要乘上某个常数
    std::cerr << "stop after " << step << " steps with residual: " << res << std::endl;
    os_res.close();
  }
}; // class FullNonlinearMultiGrid



#endif // __MULTI_GRID_H__

/**
 * end of file
 * 
 */
