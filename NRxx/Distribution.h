// vim: fdm=syntax:ts=4:sw=4:syntax=cpp.doxygen:tw=70:fo+=Mm

#ifndef __Distribution_h_
#define __Distribution_h_


#include <limits>
#include <algorithm>



#include "config.h"
#include "Vector.h"
#include "Matrix.h"
#include "CoefTable.h"

#ifdef SPECTRAL_PROJECTION
#include "SpectralForODE.h"
#endif // SPECTRAL_PROJECTION

NRXX_NAMESPACE_OPEN

/**
 * 判断条件 @p COND 是否成立，否则抛出异常 @p EXCEPTION。
 */
#define Assert(COND, EXCEPTION) \
	if (!(COND)) throw (EXCEPTION);

/**
 * 用于描述分布函数的类。一个离散了的分布函数事实上只包含离散的中心、
 * 伸缩因子和一系列谱系数。此类主要提供对分布函数的一系列运算。模板参
 * 数 @p DIM 和 @p _T 的意义参见 @ref CoefTable，模板参数 @p K 表示由
 * 于降维或多原子等原因约化后产生的额外的分布函数的数目。
 */
template <int DIM, class _T, int K> class Distribution;

/**
 * 多原子情形。其中包含了三个特殊的参数，Prandtl 数 @ref Pr、松弛碰撞
 * 数 @ref Z 以及内部自由度数 @ref delta。
 */
template <int DIM, class _T, int NR> class Distribution {
 private:
  /**
   * 单原子情形下的分布函数类。
   */
  typedef Distribution<DIM, _T, 0> _MonatomicDistribution;

  /**
   * 多原子情形事实上只需要两个完全相同的单原子分布函数类即可表
   * 示。
   */
  _MonatomicDistribution __dis[NR + 1];

  /**
   * 迭代器类型。模板参数 @p _IT 表示系数表的迭代器类型。
   */
  template <class _IT>
    class _IteratorBase {
  public:
    /**
     * 元素值类型。
     */
    typedef typename _IT::value_type value_type;

    /**
     * 元素指针类型。
     */
    typedef typename _IT::pointer pointer;

    /**
     * 元素引用类型。
     */
    typedef typename _IT::reference reference;

  private:
    /**
     * 用 @p NR+1 个系数表的迭代器来组合成多原子分布的迭代器。
     */
    _IT __it[NR + 1];

    /**
     * 该变量用于标识当前迭代器指向第几个分布函数。
     */
    unsigned int __n_dis;

    /**
     * 迭代器所指元素的阶数。
     */
    unsigned int __order;

    /**
     * 当前所指元素在分布函数中的线性位置。
     */
    unsigned int __position;

    /**
     * 获取非约化部分每阶以下矩的个数。
     */
    static std::vector<unsigned int> _GetMomentNumber() {
      std::vector<unsigned int> __n(MAX_ORDER);
      for (size_t __i = 0; __i < MAX_ORDER; __i++)
        __n[__i] = MomentNumber::_GetMomentNumber(DIM, __i);
      return __n;
    }

    /**
     * 获取约化部分每阶以下矩的个数。
     */
    static std::vector<unsigned int> _GetMomentNumber_Reduced() {
      std::vector<unsigned int> __n(MAX_ORDER, 0);
      for (size_t __i = 2; __i < MAX_ORDER; __i++)
        __n[__i] = MomentNumber::_GetMomentNumber(DIM, __i - 2);
      return __n;
    }

    /**
     * 获取 @p __order 阶以下矩的个数，@p __is_reduced 表示计算约化部
     * 分还是非约化部分。
     */
    static unsigned int _GetMomentNumber(bool __is_reduced, unsigned int __order) {
      static std::vector<unsigned int> __n = _GetMomentNumber();
      static std::vector<unsigned int> __n_reduced = _GetMomentNumber_Reduced();
      if (__is_reduced)
        return __n_reduced[__order];
      else
        return __n[__order];
    }

    /**
     * 构造函数。直接对内部成员赋值。
     */
    _IteratorBase(const _IT __it_f[NR + 1], unsigned int __nd,
        unsigned int __ord, unsigned int __pos) :
      __n_dis(__nd), __order(__ord), __position(__pos)
    {
      for (int __i = 0; __i <= NR; __i++) __it[__i] = __it_f[__i];
    }

  public:
    /**
     * 默认构造函数。
     */
    _IteratorBase() : __n_dis(0) {}

    /**
     * 前缀增量操作。
     */
    _IteratorBase& operator++() {
      ++__it[__n_dis];

      if (++__position >= _GetMomentNumber(__n_dis, __order)) {
        if (__order >= 2) __n_dis = (__n_dis + 1) % (NR + 1);
        if (!__n_dis) __order++;
        __position = _GetMomentNumber(__n_dis, __order - 1);
      }
      return *this;
    }

    /**
     * 后缀增量操作。
     */
    _IteratorBase operator++(int) {
      _IteratorBase __bak = *this;
      operator++();
      return __bak;
    }

    /**
     * 判断两个迭代器是否相等。
     */
    bool operator==(const _IteratorBase& __iter) const {
      for (int __i = 0; __i <= NR; __i++)
        if (__it[__i] != __iter.__it[__i]) return false;
      return __n_dis == __iter.__n_dis &&
        __order == __iter.__order &&
        __position == __iter.__position;
    }

    /**
     * 判断两个迭代器是否不等。
     */
    bool operator!=(const _IteratorBase& __iter) const {
      return !operator==(__iter);
    }

    /**
     * Dereference 操作。
     */
    reference operator*() {
      return __it[__n_dis].operator*();
    }

    /**
     * 取指针操作。
     */
    pointer operator->() {
      return __it[__n_dis].operator->();
    }

    /**
     * 取出元素下标。
     */
    const typename _MonatomicDistribution::Indices& GetIndices() const {
      return __it[__n_dis].GetIndices();
    }

    /**
     * 获取当前矩所在的阶数。
     */
    unsigned int Order() const {
      return __order;
    }

    /**
     * 判断当前迭代器是否位于第几个分布函数上。
     */
    unsigned int NDis() const {
      return __n_dis;
    }

    friend class Distribution<DIM, _T, NR>;
  };

 public:
  /**
   * 各系数的类型。
   */
  typedef _T value_type;

  /**
   * 速度的类型。
   */
  typedef Vector<DIM, _T> Velocity;

  /**
   * 维数。
   */
  static const int dim = DIM;

  /**
   * 约化部分分布函数的个数。
   */
  static const int n_reduced = NR;

  /**
   * 内部自由度的个数。
   */
  _T delta[NR];

  /**
   * 矩的下标类型。
   */
  typedef typename _MonatomicDistribution::Indices Indices;

  /**
   * 迭代器类型。
   */
  typedef _IteratorBase<typename CoefTable<DIM, _T>::Iterator> Iterator;

  /**
   * @p const 迭代器类型。
   */
  typedef _IteratorBase<typename CoefTable<DIM, _T>::ConstIterator> ConstIterator;

  /**
   * 默认构造函数。
   */
  Distribution() {}

  /**
   * 构造函数。指定展开中心 @p __v 和伸缩因子 @p __x。
   */
  Distribution(const Velocity& __v, const _T& __x) {
    for (int __i = 0; __i <= NR; __i++)
      __dis[__i].Reinit(__v, __x);
  }

  /**
   * 构造函数。指定展开阶数 @p __M。
   */
  explicit Distribution(unsigned int __M) {
    __dis[0].Reinit(__M);
    for (int __i = 1; __i <= NR; __i++)
      __dis[__i].Reinit(__M - 2);
  }

  /**
   * 构造函数。指定展开中心 @p __v，伸缩因子 @p __x 及展开阶数@p __M。
   */
  Distribution(const Velocity& __v, const _T& __x, unsigned int __M) {
    __dis[0].Reinit(__v, __x, __M);
    for (int __i = 1; __i <= NR; __i++)
      __dis[__i].Reinit(__v, __x, __M - 2);
  }

  /**
   * 重设展开中心 @p __v。
   */
  void Reinit(const Velocity& __v) {
    for (int __i = 0; __i <= NR; __i++)
      __dis[__i].Reinit(__v);
  }

  /**
   * 重设展开中心 @p __v 和伸缩因子 @p __x。
   */
  void Reinit(const Velocity& __v, const _T& __x) {
    for (int __i = 0; __i <= NR; __i++)
      __dis[__i].Reinit(__v, __x);
  }

  /**
   * 重设展开的阶数。
   */
  void Reinit(unsigned int __M) {
    __dis[0].Reinit(__M);
    for (int __i = 1; __i <= NR; __i++)
      __dis[__i].Reinit(__M - 2);
  }

  /**
   * 重设展开中心 @p __v 和伸缩因子 @p __x 及展开阶数 @p __M。
   */
  void Reinit(const Velocity& __v, const _T& __x, unsigned int __M) {
    __dis[0].Reinit(__v, __x, __M);
    for (int __i = 1; __i <= NR; __i++)
      __dis[__i].Reinit(__v, __x, __M - 2);
  }

  /**
   * 复制另一分布函数的所有参数。
   */
  void Reinit(const Distribution& __f) {
    __dis[0].Reinit(__f.__dis[0]);
    for (int __i = 1; __i <= NR; __i++) {
      __dis[__i].Reinit(__f.__dis[__i]);
      delta[__i-1] = __f.delta[__i-1];
    }
  }

  /**
   * 获取展开中心。
   */
  const Velocity& Center() const {
    return __dis[0].Center();
  }

  /**
   * 获取平动变量伸缩因子。
   */
  const _T& Scaling() const {
    return __dis[0].Scaling();
  }

  /**
   * 获取展开的阶数。
   */
  unsigned int GetOrder() const {
    return __dis[0].GetOrder();
  }

  /**
   * 重设展开的阶数，保持所有系数不变。
   */
  void ResetOrder(unsigned int __M) {
    __dis[0].ResetOrder(__M);
    for (int __i = 1; __i <= NR; __i++)
      __dis[__i].ResetOrder(__M - 2);
  }

  /**
   * 将分布函数置零。这个函数使得我们可以写 <tt>f = 0;</tt> 这样
   * 的语句，但为了避免歧义，参数只允许为零。
   */
  Distribution& operator=(const _T& __zero) {
    for (int __i = 0; __i <= NR; __i++) __dis[__i] = __zero;
    return *this;
  }

  /**
   * 加上一个分布函数 @p __d，为了效率，总假定两个分布函数具有相
   * 同的中心和伸缩因子。
   */
  Distribution& operator+=(const Distribution& __d) {
    for (int __i = 0; __i <= NR; __i++) __dis[__i] += __d.__dis[__i];
    return *this;
  }

  /**
   * 减去一个分布函数 @p __d，为了效率，总假定两个分布函数具有相
   * 同的中心和伸缩因子。
   */
  Distribution& operator-=(const Distribution& __d) {
    for (int __i = 0; __i <= NR; __i++) __dis[__i] -= __d.__dis[__i];
    return *this;
  }

  /**
   * 将当前分布函数乘上一个常数 @p __s。
   */
  Distribution& operator*=(const _T& __s) {
    for (int __i = 0; __i <= NR; __i++) __dis[__i] *= __s;
    return *this;
  }

  /**
   * 将当前分布函数除以一个常数 @p __s。
   */
  Distribution& operator/=(const _T& __s) {
    for (int __i = 0; __i <= NR; __i++) __dis[__i] /= __s;
    return *this;
  }

  /**
   * 将当前分布函数加上 @p __s 倍的 @p __d。
   */
  void Add(const _T& __s, const Distribution& __d) {
    for (int __i = 0; __i <= NR; __i++) __dis[__i].Add(__s, __d.__dis[__i]);
  }

 
  /**
   * 取分布函数的一些矩。@p _VECTOR 是一个向量类型，@p __v 必须
   * 至少含有 @p DIM+2 个元素。在函数运行完成后，<tt>__v[0]</tt>
   * 保存密度，<tt>__v[1]</tt> 至 <tt>__v[DIM]</tt> 保存动量，@p
   * __v 的最后一个元素保存能量。
   */
  template <class _VECTOR>
    void Moments(_VECTOR& __v) const {
    __dis[0].Moments(__v);

    _T __T = __dis[0].Scaling() * __dis[0].Scaling();
    for (int __i = 1; __i <= NR; __i++)
      __v[DIM + 1] += 0.5 * delta[__i - 1] * __T * __dis[0].front() - __dis[__i].front();
  }

  /**
   * 这个函数重载使得该成员函数对 @p const 指针可用。
   */
  void Moments(_T * const __v) const {
    _T * __tmp = __v;
    Moments<_T *>(__tmp);
  }

  /**
   * 获取对应于该分布函数的原始变量。@p _VECTOR 是一个向量类型，
   * @p __v 必须至少含有 @p DIM+2 个元素。在函数运行完成后，
   * <tt>__v[0]</tt> 保存密度，<tt>__v[1]</tt> 至
   * <tt>__v[DIM]</tt> 保存速度，@p __v 的最后一个元素保存温度。
   */
  template <class _VECTOR>
    void PrimitiveVars(_VECTOR& __v) const {
    _T __tmp = 0;
    Moments(__v);
    for (unsigned int __i = 1; __i <= DIM; __i++) {
      __v[__i] /= __v[0];
      __tmp += __v[__i] * __v[__i];
    }

    _T __delta = DIM;
    for (int __i = 0; __i < NR; __i++) __delta += delta[__i];
    __v[DIM + 1] = (2 * __v[DIM + 1] / __v[0] - __tmp) / __delta;
  }

  /**
   * 这个函数重载使得该成员函数对 @p const 指针可用。
   */
  void PrimitiveVars(_T * const __v) const {
    _T * __tmp = __v;
    PrimitiveVars<_T *>(__tmp);
  }

  /**
   * 计算宏观速度。@p _VECTOR 为一向量类型，须至少含有 @p DIM 个
   * 元素。函数运行完成后 @p __v 中保存宏观速度的各个分量。
   */
  template <class _VECTOR>
    void MacroVelocity(_VECTOR& __v) const {
    __dis[0].MacroVelocity(__v);
  }

  /**
   * 这个函数重载使得该成员函数对 @p const 指针可用。
   */
  void MacroVelocity(_T * const __v) const {
    _T * __tmp = __v;
    MacroVelocity<_T *>(__tmp);
  }

  /**
   * 取平动温度。
   */
  _T TransTemperature() const {
    return __dis[0].TransTemperature();
  }

  /**
   * 取第 @p __i 个分布函数所给出来的温度。
   */
  _T DirecTemperature(unsigned int __i) const {
    if (!__i) return TransTemperature();
    _T __T = Scaling() * Scaling();
    return __T - 2 * __dis[__i].front() / (__dis[0].front() * delta[__i - 1]);
  }

  /**
   * 获取热通量。
   */
  template <class _VECTOR>
    void HeatFlux(_VECTOR& __v) const {
    __dis[0].HeatFlux(__v);

    _T __velocity[DIM];
    MacroVelocity(__velocity);

    _MonatomicDistribution __tmp_dis;
    __tmp_dis.Reinit(__dis[0]); __tmp_dis.Reinit(1);
    _T __theta = __dis[0].Scaling() * __dis[0].Scaling();
    _T __delta = 0;
    for (int __i = 0; __i < NR; __i++) {
      __delta += delta[__i];
      __tmp_dis.Add(-1.0, __dis[__i+1]);
    }
    __tmp_dis.Add(0.5 * __delta * __theta, __dis[0]);

    const Velocity& __c = __tmp_dis.Center();
    Indices __ind;

    for (int __i = 0; __i < DIM; __i++) {
      __ind[__i] = 1;
      __v[__i] += (__c[__i] - __velocity[__i]) * __tmp_dis.front() + __tmp_dis(__ind);
      __ind[__i] = 0;
    }
  }

  /**
   * 这个函数重载使得该成员函数对 @p const 指针可用。
   */
  void HeatFlux(_T * const __v) const {
    _T * __tmp = __v;
    HeatFlux<_T *>(__tmp);
  }

  /**
   * 使用精确投影方式将分布函数 @p __f 投影到当前分布函数的中心和伸缩
   * 因子下。参数 @p __m 指定投影的最高阶数，其值为负数时表示由函数自
   * 动选择最高阶数。高于这个最高阶数的部分将被清零。
   */
  void ExactProject(const Distribution& __f, int __m = -1) {
    // 要求两个分布函数有相同的 delta 值
    for (int __i = 0; __i < NR; __i++)
      delta[__i] = __f.delta[__i];

    __dis[0].ExactProject(__f.__dis[0], __m);
    _T __T = Scaling() * Scaling();
    _T __T_f = __f.Scaling() * __f.Scaling();
    for (unsigned int __i = 1; __i <= NR; __i++) {
      __dis[__i].ExactProject(__f.__dis[__i], __m - 2);
      __dis[__i].Add(0.5 * delta[__i - 1] * (__T - __T_f), __dis[0]);
    }
  }

  /**
   * 将分布函数 @p __f 投影到当前分布函数的中心和伸缩因子下。参
   * 数 @p __m 指定投影的最高阶数，其值为负数时表示由函数自动选
   * 择最高阶数。高于这个最高阶数的部分将被清零。
   */
  void Project(const Distribution& __f, int __m = -1) {
    // 要求两个分布函数有相同的 delta 值
    for (int __i = 0; __i < NR; __i++)
      delta[__i] = __f.delta[__i];

    __dis[0].Project(__f.__dis[0], __m);
    _T __T = Scaling() * Scaling();
    _T __T_f = __f.Scaling() * __f.Scaling();
    for (unsigned int __i = 1; __i <= NR; __i++) {
      __dis[__i].Project(__f.__dis[__i], __m - 2);
      __dis[__i].Add(0.5 * delta[__i - 1] * (__T - __T_f), __dis[0]);
    }
  }

  /**
   * 把分布函数 @p __f 投影到其标准空间。这与 @ref Project 仅有
   * 的不同在于该函数重置了展开中心和伸缩因子。
   */
  void ExactProjectToStdSpace(const Distribution& __f, int __m = -1) {
    std::vector<_T> __pv(DIM + 2);
    __f.PrimitiveVars(__pv);
    Reinit(Velocity(&__pv[1]), sqrt(__pv[DIM + 1]));
    for (unsigned int __i = 0; __i < NR; __i++)
      delta[__i] = __f.delta[__i];
    ExactProject(__f, __m);
  }

  /**
   * 把分布函数 @p __f 投影到其标准空间。这与 @ref Project 仅有
   * 的不同在于该函数重置了展开中心和伸缩因子。
   */
  void ProjectToStdSpace(const Distribution& __f, int __m = -1) {
    std::vector<_T> __pv(DIM + 2);
    __f.PrimitiveVars(__pv);
    Reinit(Velocity(&__pv[1]), sqrt(__pv[DIM + 1]));
    for (unsigned int __i = 0; __i < NR; __i++)
      delta[__i] = __f.delta[__i];
    Project(__f, __m);
  }

  /**
   * 把当前分布函数乘上第 @p __i 维的速度，结果存在 @p __f 中。
   */
  void MulVelocity(unsigned int __i, Distribution& __f) const {
    for (unsigned int __k = 0; __k < NR; __k++)
      __f.delta[__k] = delta[__k];
    for (unsigned int __k = 0; __k <= NR; __k++)
      __dis[__k].MulVelocity(__i, __f.__dis[__k]);
  }

  /**
   * 生成一个麦克斯韦分布函数，以展开中心为平均速度，伸缩因子的
   * 平方为温度，@p __rho 为密度。
   */
  void MakeMaxwellian(const _T& __rho) {
    __dis[0].MakeMaxwellian(__rho);
    for (unsigned int __i = 1; __i <= NR; __i++) __dis[__i] = 0;
  }

  /**
   * 通过下标取元素。
   */
  _T& operator()(const Indices& __ind, unsigned int __n_dis) {
    return __dis[__n_dis](__ind);
  }

  /**
   * 通过下标取元素。
   */
  const _T& operator()(const Indices& __ind, unsigned int __n_dis) const {
    return __dis[__n_dis](__ind);
  }

  /**
   * 指向起始元素的迭代器。
   */
  Iterator begin() {
    typename _MonatomicDistribution::Iterator __it[NR + 1];
    for (unsigned int __i = 0; __i <= NR; __i++)
      __it[__i] = __dis[__i].begin();
    return Iterator(__it, 0, 0, 0);
  }

  /**
   * 指向起始元素的 @p const 迭代器。
   */
  ConstIterator begin() const {
    typename _MonatomicDistribution::ConstIterator __it[NR + 1];
    for (unsigned int __i = 0; __i <= NR; __i++)
      __it[__i] = __dis[__i].begin();
    return ConstIterator(__it, 0, 0, 0);
  }

  /**
   * 指向末尾元素的迭代器。
   */
  Iterator end() {
    unsigned int __M = __dis[0].GetOrder();
    typename _MonatomicDistribution::Iterator __it[NR + 1];
    for (unsigned int __i = 0; __i <= NR; __i++)
      __it[__i] = __dis[__i].end();
    return Iterator(__it, 0, __M + 1, MomentNumber::_GetMomentNumber(DIM, __M));
  }

  /**
   * 指向末尾元素的 @p const 迭代器。
   */
  ConstIterator end() const {
    unsigned int __M = __dis[0].GetOrder();
    typename _MonatomicDistribution::ConstIterator __it[NR + 1];
    for (unsigned int __i = 0; __i <= NR; __i++)
      __it[__i] = __dis[__i].end();
    return ConstIterator(__it, 0, __M + 1, MomentNumber::_GetMomentNumber(DIM, __M));
  }

  /**
   * 指向 @p __m 阶矩起始元素的迭代器。
   */
  Iterator begin(unsigned int __m) {
    typename _MonatomicDistribution::Iterator __it[NR + 1];
    if (__m == 0)
      return begin();
    else if (__m < 2) {
      __it[0] = __dis[0].begin(__m);
      for (unsigned int __i = 1; __i <= NR; __i++)
        __it[__i] = __dis[__i].begin();
      return Iterator(__it, 0, __m, MomentNumber::_GetMomentNumber(DIM, __m - 1));
    } else {
      __it[0] = __dis[0].begin(__m);
      for (unsigned int __i = 1; __i <= NR; __i++)
        __it[__i] = __dis[__i].begin(__m - 2);
      return Iterator(__it, 0, __m, MomentNumber::_GetMomentNumber(DIM, __m - 1));
    }
  }

  /**
   * 指向 @p __m 阶矩的起始元素的 @p const 迭代器。
   */
  ConstIterator begin(unsigned int __m) const {
    typename _MonatomicDistribution::ConstIterator __it[NR + 1];
    if (__m == 0)
      return begin();
    else if (__m < 2) {
      __it[0] = __dis[0].begin(__m);
      for (unsigned int __i = 1; __i <= NR; __i++)
        __it[__i] = __dis[__i].begin();
      return ConstIterator(__it, 0, __m, MomentNumber::_GetMomentNumber(DIM, __m - 1));
    } else {
      __it[0] = __dis[0].begin(__m);
      for (unsigned int __i = 1; __i <= NR; __i++)
        __it[__i] = __dis[__i].begin(__m - 2);
      return ConstIterator(__it, 0, __m, MomentNumber::_GetMomentNumber(DIM, __m - 1));
    }
  }

  /**
   * 指向 @p __m 阶矩末尾元素的迭代器。
   */
  Iterator end(unsigned int __m) {
    return begin(__m + 1);
  }

  /**
   * 指向 @p __m 阶矩的末尾元素的 @p const 迭代器。
   */
  ConstIterator end(unsigned int __m) const {
    return begin(__m + 1);
  }

  /**
   * 取第一个系数，即密度。
   */
  _T& front() { return __dis[0].front(); }

  /**
   * 取第一个系数，即密度，@p const 版本。
   */
  const _T& front() const { return __dis[0].front(); }

  /**
   * 对分布函数做 Givens 变换。设变换矩阵为 @p A，则
   * \f[
   *     A_{ii} = A_{jj} = \cos\theta, \quad
   *     A_{ij} = -A_{ji} = \sin\theta, \quad
   *     A_{kl} = \delta_{kl}, \qquad
   *     k,l \neq i,j.
   * \f]
   */
  void GivensTransform(int __i, int __j, _T __theta, const Distribution& __src)
  {
    for (int __k = 0; __k < NR; ++__k)
      delta[__k] = __src.delta[__k];
    for (int __k = 0; __k <= NR; ++__k) {
      __dis[__k].GivensTransform(__i, __j, __theta, __src.__dis[__k]);
    }
  }

  /**
   * 对分布函数做一正交变换，变换矩阵由 @p __A 给出。@p __src 为
   * 输入分布函数 f(v)，变换得到的结果 f(Av) 存储在当前分布函数
   * 中。用户需自行保证 @p __A 为一正交矩阵，程序只对 @p __A 是
   * 否为方阵进行检验。同时 @p __A 的行数 @p m 须小于或等于 @p
   * DIM，当小于时，可认为变换矩阵为
   *   \f[ \mathrm{diag}\{A, I_{D-m}\}, \f]
   * 其中 \f$D\f$ 即为 @p DIM。
   */
  void OrthTransform(const Matrix<_T>& __A, const Distribution& __src)
  {
    for (int __i = 0; __i < NR; ++__i)
      delta[__i] = __src.delta[__i];
    for (int __k = 0; __k <= NR; ++__k) {
      __dis[__k].OrthTransform(__A, __src.__dis[__k]);
    }
  }

  /**
   * 对分布函数 @p __src 做一镜面反射变换，其中 @p __normal 为镜
   * 面单位法向，结果保存在当前分布函数中。用户需自行保证 @p
   * __normal 的范数为一。
   */
  template <class _VECTOR>
    void SymmTransform(const _VECTOR& __normal, const Distribution& __src) {
    for (int __i = 0; __i < NR; ++__i)
      delta[__i] = __src.delta[__i];
    for (int __k = 0; __k <= NR; ++__k) {
      __dis[__k].SymmTransform(__normal, __src.__dis[__k]);
    }
  }

  /**
   * 对分布函数 @p __src 做截断，截平面与第 @p __k 个坐标轴垂直
   * ，且第 @p __k 维坐标为 @p __u。当 @p __sign 为零或正时取正
   * 的部分，否则取负的部分。
   */
  void HalfSpaceCutOff(unsigned int __k, _T __u,
                       const Distribution& __src, int __sign = 1) {
    for (int __i = 0; __i < NR; ++__i)
      delta[__i] = __src.delta[__i];
    for (int __i = 0; __i <= NR; ++__i) {
      __dis[__i].HalfSpaceCutOff(__k, __u, __src.__dis[__i], __sign);
    }
  }
};

template <int DIM, class _T>
  class Distribution<DIM, _T, 0> : public CoefTable<DIM, _T>
{
 private:
  /**
   * 展开时的中心位置，一般为分布函数所对应的速度。
   */
  Vector<DIM, _T> __c;

  /**
   * 平动变量的伸缩因子，一般为平动温度的开方。
   */
  _T __s;

  /**
   * 迭代器类型。其中模板参数 @p _IT 表示系数表的迭代器类型。
   */
  template <class _IT>
    class _IteratorBase : public _IT {
  public:
    /**
     * 默认构造函数。
     */
  _IteratorBase() : _IT() {}

    /**
     * 构造函数。从基类复制。
     */
  _IteratorBase(const _IT& __it) : _IT(__it) {}

    /**
     * 无约化情形下仅有一个分布函数。
     */
    static unsigned int NDis() {
      return 0;
    }
  };

  /**
   * 求一个角度 @p theta，使得
   * <tt>
   *   -x * sin(theta) + y * cos(theta) = 0.
   * </tt>
   */
  _T _GivensParameter(_T __x, _T __y) {
    return atan2(__y, __x);
  }

  /**
   * 将矩阵 @p __M 左乘一个 Givens 矩阵，这个 Givens 矩阵形如
   * \f[
   *     A_{ii} = A_{jj} = \cos\theta, \quad
   *     A_{ij} = -A_{ji} = \sin\theta, \quad
   *     A_{kl} = \delta_{kl}, \qquad k,l \neq i,j.
   * \f]
   */
  void _GivensTransform(unsigned int __i, unsigned int __j,
                        _T __theta, Matrix<_T>& __M)
  {
    _T __cos = cos(__theta);
    _T __sin = sin(__theta);
    unsigned int __N = __M.n();
    for (unsigned int __k = 0; __k < __N; __k++) {
      _T& __a = __M(__i, __k);
      _T& __b = __M(__j, __k);
      _T __c = __b;

      __b = __b * __cos - __a * __sin;
      __a = __a * __cos + __c * __sin;
    }
  }

public:
  void ProjectRK(const Distribution& __f, int __m, int __Nstep)
  {
    int __n = 0;
    _T __dt = (_T)1 / __Nstep;

    Velocity __du_dt = (__f.Center() - __c) * __dt;
    _T __half_dtheta_dt = 0.5 * (__f.__s * __f.__s - __s * __s) * __dt;

    // 采用四级四阶 Runge-Kutta 方法求解
    _Base __t[4];
    Iterator __it, __it_end;
    typename _Base::Iterator __it1[DIM], __it2[DIM], __it1_pt[DIM], __it2_pt[DIM];

    do {
      for (unsigned int __i = 0; __i < 4; __i++)
        __t[__i].Reinit(__m);

#define RK4(I, DELTA_COEF)                                              \
      if (__ind[__j] >= (I)) {                                          \
        *__it += (DELTA_COEF) * ((*__it##I[__j]) + __k[__i] * (*__it##I##_pt[__j])); \
        ++__it##I[__j], ++__it##I##_pt[__j];                            \
      }

      _T __k[4] = {0, .5, .5, 1};
      for (unsigned int __i = 0; __i < 4; __i++) {
        _Base* __pt = (__i ? &__t[__i - 1] : this);
        __it = __t[__i].begin();
        __it_end = __t[__i].end();

        for (int __j = 0; __j < DIM; __j++) {
          __it1[__j] = _Base::begin();
          __it2[__j] = _Base::begin();
          __it1_pt[__j] = __pt->begin();
          __it2_pt[__j] = __pt->begin();
        }

        for (; __it != __it_end; ++__it) {
          const Indices& __ind = __it.GetIndices();
          for (int __j = 0; __j < DIM; __j++) {
            RK4(2, __half_dtheta_dt);
            RK4(1, __du_dt[__j]);
          }
        }
      }
#undef RK4

      _T *__p = &_Base::front(), *__p_end = __p + __t[0].size();
      _T *__p_t[4] = {&__t[0].front(), &__t[1].front(), &__t[2].front(), &__t[3].front()};
      for (; __p != __p_end; ++__p)
      {
        *__p += ((*__p_t[0]) + 2 * (*__p_t[1]) + 2 * (*__p_t[2]) + (*__p_t[3])) / 6;
        ++__p_t[0]; ++__p_t[1]; ++__p_t[2]; ++__p_t[3];
      }
    } while (++ __n < __Nstep);
  }

 public:
  /**
   * 为与约化情形保持一致而预留的指针。
   */
  static const _T * delta;

  /**
   * 各系数的类型。
   */
  typedef _T value_type;

  /**
   * 维数。
   */
  static const int dim = DIM;

  /**
   * 约化部分分布函数个数，总为零。
   */
  static const int n_reduced = 0;

  /**
   * 基类类型。
   */
  typedef CoefTable<DIM, _T> _Base;

  /**
   * 基类类型的 @p const 版本。
   */
  typedef const CoefTable<DIM, _T> _CBase;

  /**
   * 速度类型。
   */
  typedef Vector<DIM, _T> Velocity;

  /**
   * 元素指标类型。
   */
  typedef typename _Base::Indices Indices;

  /**
   * 迭代器类型。
   */
  typedef _IteratorBase<typename _Base::Iterator> Iterator;

  /**
   * @p const 迭代器类型。
   */
  typedef _IteratorBase<typename _Base::ConstIterator> ConstIterator;

  /**
   * 默认构造函数。
   */
 Distribution() : _Base() {}

  /**
   * 构造函数。指定展开中心 @p __v 和伸缩因子 @p __x。
   */
 Distribution(const Velocity& __v, const _T& __x) :
  _Base(), __c(__v), __s(__x) {}

  /**
   * 构造函数。指定展开阶数 @p __M。
   */
 explicit Distribution(unsigned int __M) : _Base(__M) {}

  /**
   * 构造函数。指定展开中心 @p __v，伸缩因子 @p __x 及展开阶数
   * @p __M。
   */
 Distribution(const Velocity& __v, const _T& __x, unsigned int __M) :
  _Base(__M), __c(__v), __s(__x) {}

  /**
   * 重设展开中心 @p __v。
   */
  void Reinit(const Velocity& __v) {
    __c = __v;
  }

  /**
   * 重设展开中心 @p __v 和伸缩因子 @p __x。
   */
  void Reinit(const Velocity& __v, const _T& __x) {
    __c = __v; __s = __x;
  }

  /**
   * 重设展开阶数 @p __M。
   */
  void Reinit(unsigned int __M) {
    _Base::Reinit(__M);
  }

  /**
   * 重设展开中心 @p __v 和伸缩因子 @p __x 及展开阶数 @p __M。
   */
  void Reinit(const Velocity& __v, const _T& __x, unsigned int __M) {
    __c = __v; __s = __x;
    _Base::Reinit(__M);
  }

  /**
   * 复制另一分布函数的所有参数。
   */
  void Reinit(const Distribution& __f) {
    Reinit(__f.Center(), __f.Scaling(), __f.GetOrder());
  }

  /**
   * 获取展开中心。
   */
  const Velocity& Center() const {
    return __c;
  }

  /**
   * 获取伸缩因子。
   */
  const _T& Scaling() const {
    return __s;
  }

  /**
   * 将分布函数置零。这个函数使得我们可以写 <tt>f = 0;</tt> 这样
   * 的语句，但为了避免歧义，参数只允许为零。
   */
  Distribution& operator=(const _T& __zero) {
    _Base::operator=(__zero);
    return *this;
  }

  /**
   * 加上一个分布函数 @p __d，为了效率，总假定两个分布函数具有相
   * 同的中心和伸缩因子。
   */
  Distribution& operator+=(const Distribution& __d) {
    Iterator __it = _Base::begin();
    Iterator __it_end = _Base::end();
    ConstIterator __it_d = __d._Base::begin();
    ConstIterator __it_d_end = __d._Base::end();

    for (; __it != __it_end && __it_d != __it_d_end;
         ++__it, ++__it_d) *__it += *__it_d;
    return *this;
  }

  /**
   * 减去一个分布函数 @p __d，为了效率，总假定两个分布函数具有相
   * 同的中心和伸缩因子。
   */
  Distribution& operator-=(const Distribution& __d) {
    Iterator __it = _Base::begin();
    Iterator __it_end = _Base::end();
    ConstIterator __it_d = __d._Base::begin();
    ConstIterator __it_d_end = __d._Base::end();

    for (; __it != __it_end && __it_d != __it_d_end;
         ++__it, ++__it_d) *__it -= *__it_d;
    return *this;
  }

  /**
   * 将当前分布函数乘上一个常数 @p __s。
   */
  Distribution& operator*=(const _T& __s) {
    Iterator __it = _Base::begin();
    Iterator __it_end = _Base::end();

    for (; __it != __it_end; ++__it) *__it *= __s;
    return *this;
  }

  /**
   * 将当前分布函数除以一个常数 @p __s。
   */
  Distribution& operator/=(const _T& __s) {
    Iterator __it = _Base::begin();
    Iterator __it_end = _Base::end();

    for (; __it != __it_end; ++__it) *__it /= __s;
    return *this;
  }

  /**
   * 获取指标 @p __i 对应的元素。
   */
  const _T& operator()(const Indices& __i, int __n_dis = 0) const {
    return _CBase::operator()(__i);
  }

  /**
   * 获取指标 @p __i 对应的元素。
   */
  _T& operator()(const Indices& __i, int __n_dis = 0) {
    return _Base::operator()(__i);
  }

  /**
   * 将当前分布函数加上 @p __s 倍的 @p __d。
   */
  void Add(const _T& __s, const Distribution& __d) {
    Iterator __it = _Base::begin();
    Iterator __it_end = _Base::end();
    ConstIterator __it_d = __d._Base::begin();
    ConstIterator __it_d_end = __d._Base::end();

    for (; __it != __it_end && __it_d != __it_d_end;
         ++__it, ++__it_d) *__it += __s * *__it_d;
  } 

  /**
   * 取分布函数的一些矩。@p _VECTOR 是一个向量类型，@p __v 必须
   * 至少含有 @p DIM+2 个元素。在函数运行完成后，<tt>__v[0]</tt>
   * 保存密度，<tt>__v[1]</tt> 至 <tt>__v[DIM]</tt> 保存动量，@p
   * __v 的最后一个元素保存能量。
   */
  template <class _VECTOR>
    void Moments(_VECTOR& __v) const {
    Assert(_Base::GetOrder() >= 2, std::length_error("Moment order insufficient"));
    _T __s2 = __s * __s;
    Indices __ind;
    __v[0] = _Base::operator()(__ind);

    _T __tmp1 = 0, __tmp2 = 0, __tmp3 = 0;
    for (int __i = 0; __i < DIM; __i++) {
      __ind[__i] = 1;
      __v[1 + __i] = __c[__i] * __v[0] + _Base::operator()(__ind);
      __tmp1 += __c[__i] * __v[1 + __i];
      __tmp2 += __c[__i] * __c[__i];

      __ind[__i] = 2;
      __tmp3 += _Base::operator()(__ind);
      __ind[__i] = 0;
    }
    __v[DIM + 1] = __tmp1 + __tmp3 +
      0.5 * (DIM * __s2 - __tmp2) * __v[0];
  }

  /**
   * 这个函数重载使得该成员函数对 @p const 指针可用。
   */
  void Moments(_T * const __v) const {
    _T * __tmp = __v;
    Moments<_T *>(__tmp);
  }

  /**
   * 计算热通量。@p _VECTOR 是一个向量类型，@p __v 必须至少含有
   * @p DIM 个元素。在函数运行完成后，<tt>__v[0]</tt> 至
   * <tt>__v[DIM-1]</tt> 保存 @p DIM 个方向上的热通量。
   */
  template <class _VECTOR>
    void HeatFlux(_VECTOR& __v) const {
    Assert(_Base::GetOrder() >= 3, std::length_error("Moment order insufficient"));

    double __var[DIM + 2];
    PrimitiveVars(__var);

    double __tmp[DIM];
    for(int __i = 0; __i < DIM; ++__i)
      __tmp[__i] = __c[__i] - __var[__i + 1];

    for (int __i = 0; __i < DIM; __i++) {
      Indices __ind;
      __v[__i] = (__tmp[__i] * __tmp[__i] * __tmp[__i] + 3 * __tmp[__i] * __s * __s) * (*this)(__ind) / 2;
      __ind[__i] = 3;
      __v[__i] += 3 * (*this)(__ind);
      __ind[__i] = 2;
      __v[__i] += 3 * __tmp[__i] * (*this)(__ind);
      __ind[__i] = 1;
      __v[__i] += 1.5 * (__tmp[__i] * __tmp[__i] + __s * __s) * (*this)(__ind);

      for (int __j = 0; __j < DIM; __j++) {
        if (__j != __i) {
          __ind[__i] = 1;
          __ind[__j] = 2;
          __v[__i] += (*this)(__ind);
          __ind[__j] = 1;
          __v[__i] += __tmp[__j] * (*this)(__ind);
          __ind[__j] = 0;
          __v[__i] += (__s * __s + __tmp[__j] * __tmp[__j]) * (*this)(__ind)/2;
          __ind[__i] = 0;
          __ind[__j] = 2;
          __v[__i] += __tmp[__i] * (*this)(__ind);
          __ind[__j] = 1;
          __v[__i] += __tmp[__i] * __tmp[__j] * (*this)(__ind);
          __ind[__j] = 0;
          __v[__i] += __tmp[__i] * (__s * __s + __tmp[__j] * __tmp[__j]) * (*this)(__ind)/2;
        }
      }
    }
  }

  /**
   * 这个函数重载使得该成员函数对 @p const 指针可用。
   */
  void HeatFlux(_T * const __v) const {
    _T * __tmp = __v;
    HeatFlux<_T *>(__tmp);
  }

  /**
   * 获取对应于该分布函数的原始变量。@p _VECTOR 是一个向量类型，
   * @p __v 必须至少含有 @p DIM+2 个元素。在函数运行完成后，
   * <tt>__v[0]</tt> 保存密度，<tt>__v[1]</tt> 至
   * <tt>__v[DIM]</tt> 保存速度，@p __v 的最后一个元素保存温度。
   */
  template <class _VECTOR>
    void PrimitiveVars(_VECTOR& __v) const {
    _T __tmp = 0;
    Moments(__v);
    const double epsilon = 1.e-44;
    if (fabs(__v[0]) < epsilon) {
      //
      __v[0] = 0.0;
      for (unsigned int __i = 1; __i<=DIM;__i++)
        __v[__i] = __c[__i-1];
      __v[DIM+1] =__s*__s;
      return;
    } else {
      for (unsigned int __i = 1; __i <= DIM; __i++) {
        __v[__i] /= __v[0];
        __tmp += __v[__i] * __v[__i];
      }
      __v[DIM + 1] = (__v[DIM + 1] * 2 / __v[0] - __tmp) / DIM;
    }
  }

  /**
   * 这个函数重载使得该成员函数对 @p const 指针可用。
   */
  void PrimitiveVars(_T * const __v) const {
    _T * __tmp = __v;
    PrimitiveVars<_T *>(__tmp);
  }

  /**
   * 计算宏观速度。@p _VECTOR 为一向量类型，须至少含有 @p DIM 个
   * 元素。函数运行完成后 @p __v 中保存宏观速度的各个分量。
   */
  template <class _VECTOR>
    void MacroVelocity(_VECTOR& __v) const {
    _T __pv[DIM + 2];
    PrimitiveVars(__pv);
    for (int __i = 0; __i < DIM; __i++)
      __v[__i] = __pv[__i + 1];
  }

  /**
   * 这个函数重载使得该成员函数对 @p const 指针可用。
   */
  void MacroVelocity(_T * const __v) const {
    _T * __tmp = __v;
    MacroVelocity<_T *>(__tmp);
  }

  /**
   * 取平移温度。单原子情形下平移温度也就是平均温度。
   */
  _T TransTemperature() const {
    _T __pv[DIM + 2];
    PrimitiveVars(__pv);
    return __pv[DIM + 1];
  }

  /**
   * 为了和约化情形保持一致而保留的成员函数。
   */
  _T DirecTemperature(unsigned int) const {
    return TransTemperature();
  }

  /**
   * 将分布函数 @p __f 投影到当前分布函数的中心和伸缩因子下。参
   * 数 @p __m 指定投影的最高阶数，其值为负数时表示由函数自动选
   * 择最高阶数。
   */
  void ExactProject(const Distribution& __f, int __m = -1)
  {
    if (__m < 0)
      __m = _Base::GetOrder();
    else
      __m = std::min<int>(_Base::GetOrder(), __m);

    unsigned int __n = std::min(&(*_Base::begin(__m + 1)) - &(*_Base::begin()),
                                &(*__f._Base::end()) - &(*__f._Base::begin()));
    _Base::operator=((_T)0);
    memcpy(&(*_Base::begin()), &(*__f._Base::begin()), sizeof(_T) * __n);

    Velocity __du = __f.Center() - __c;
    _T __half_dtheta = 0.5 * (__f.__s * __f.__s - __s * __s);

    _Base *__aux1 = new _Base(__m), *__aux2 = new _Base(__m);
    unsigned int size = std::min(__aux1->size(), __f.size());
    memcpy(&__aux1->front(), &__f._Base::front(), size * sizeof(_T));
    memcpy(&_Base::front(), &__f._Base::front(), size * sizeof(_T));

    Iterator __it_self;
    typename _Base::Iterator __it1[DIM], __it2[DIM], __it, __it_end;

    for (int __k = 1; __k <= __m; __k++) {
      *__aux2 = 0;

#define RECUR(I, DELTA_COEF)                                            \
      if (__ind[__j] >= (I)) {                                          \
        *__it += (DELTA_COEF) * (*__it##I[__j]);                        \
        ++__it##I[__j];                                                 \
      }

      __it = __aux2->begin(__k);
      __it_end = __aux2->end();

      for (int __j = 0; __j < DIM; __j++) {
        __it1[__j] = __aux1->begin(__k - 1);
        __it2[__j] = __aux1->begin(__k > 1 ? __k - 2 : 0);
      }

      for (__it_self = begin(__k); __it != __it_end; ++__it, ++__it_self) {
        const Indices& __ind = __it.GetIndices();
        for (int __j = 0; __j < DIM; __j++) {
          RECUR(2, __half_dtheta);
          RECUR(1, __du[__j]);
        }

        *__it /= __k;
        *__it_self += *__it;
      }
#undef RECUR

      std::swap(__aux1, __aux2);
    }

    delete __aux1; delete __aux2;
  }

  /**
   * 将分布函数 @p __f 投影到当前分布函数的中心和伸缩因子下。参
   * 数 @p __m 指定投影的最高阶数，其值为负数时表示由函数自动选
   * 择最高阶数。高于这个最高阶数的部分将被清零。求解中使用RK4
   * 和精确求解的混合方法，当两个值非常相近时，只需要1～2步RK4
   * 即可实现；当两个值不是非常相近时则使用精确求解的办法。
   */
  void Project(const Distribution& __f, int __m = -1)
  {
    if (__m < 0)
      __m = _Base::GetOrder();
    else
      __m = std::min<int>(_Base::GetOrder(), __m);

    // 先将 @p __f 拷贝到当前分布上
    unsigned int __n = std::min(&(*_Base::begin(__m + 1)) - &(*_Base::begin()),
                                &(*__f._Base::end()) - &(*__f._Base::begin()));
    _Base::operator=((_T)0);
    memcpy(&(*_Base::begin()), &(*__f._Base::begin()), sizeof(_T) * __n);

    double __s_f = __f.Scaling();

    double __tr = __s_f / __s;
    Velocity __w = __f.Center() - __c;

    int __Nstep1 = int(__w.Length() / 0.1);
    int __Nstep2 = int((__tr - 1) / 0.01);
    int __Nstep3 = int((1 / __tr - 1) / 0.01);
    int __Nstep = std::max(std::max(__Nstep1, __Nstep2), __Nstep3) + 1;
    if(__Nstep < 4)
      ProjectRK(__f, __m, __Nstep);
    else
      ExactProject(__f, __m);
  }

  /**
   * 把分布函数 @p __f 投影到其标准空间。这与 @ref ExactProject
   * 仅有的不同在于该函数重置了 @ref __c 和 @ref __s。
   */
  void ExactProjectToStdSpace(const Distribution& __f, int __m = -1) {
      std::vector<_T> __pv(DIM + 2);
      __f.PrimitiveVars(__pv);
      Reinit(Velocity(&__pv[1]), sqrt(__pv[DIM + 1]));
      ExactProject(__f, __m);
  }

  /**
   * 把分布函数 @p __f 投影到其标准空间。这与 @ref Project 仅有
   * 的不同在于该函数重置了 @ref __c 和 @ref __s。
   */
  void ProjectToStdSpace(const Distribution& __f, int __m = -1) {
      std::vector<_T> __pv(DIM + 2);
      __f.PrimitiveVars(__pv);
      Reinit(Velocity(&__pv[1]), sqrt(__pv[DIM + 1]));
      Project(__f, __m);
  }

  /**
   * 把当前分布函数乘上第 @p __i 维的速度，结果存在 @p __f 中。
   */
  void MulVelocity(unsigned int __i, Distribution& __f) const {
    unsigned int __M = _Base::GetOrder();
    __f.Reinit(__c, __s, __M);

    _T __s2 = __s * __s;
    Iterator __it = __f.begin();

#define ADD_TO_RESULT(D, S) {                   \
      __ind[__i] += (D);                        \
      *__it += (S) * (*this)(__ind);            \
      __ind[__i] -= (D);                        \
    }

    Iterator __it_end = __f.begin(__M);
    for (; __it != __it_end; ++__it) {
      Indices __ind(__it.GetIndices());
      *__it = __c[__i] * (*this)(__ind);

      if (__ind[__i] > 0) ADD_TO_RESULT(-1, __s2);
      ADD_TO_RESULT(1, __ind[__i]);
    }

    __it_end = __f.end();
    for (; __it != __it_end; ++__it) {
      Indices __ind(__it.GetIndices());
      *__it = __c[__i] * (*this)(__ind);
      if (__ind[__i] > 0) ADD_TO_RESULT(-1, __s2);
    }
  }

  /**
   * 把当前分布函数乘上第 @p __i 维的速度，结果存在 @p __f 中。
   * 这个函数只对当前分布函数的前 @p __M 阶做乘法，且不改变 @p
   * __f 的阶数。
   */	
  void MulVelocity(unsigned int __i, Distribution& __f, unsigned int __M) const {
	__M = std::min(_Base::GetOrder(), __M);
	__f.Reinit(__c, __s); __f = 0;
	unsigned int __M_f = __f._Base::GetOrder();
	_T __s2 = __s * __s;

	for (unsigned int __j = 0; __j <= __M_f; ++__j) {
      if (__j > __M + 1) break;
      Iterator __it = __f.begin(__j);
      Iterator __it_end = __f.begin(__j + 1);

      for (; __it != __it_end; ++__it) {
        Indices __ind(__it.GetIndices());
        if (__j <= __M)
          *__it = __c[__i] * (*this)(__ind);

        if (__ind[__i] > 0) ADD_TO_RESULT(-1, __s2);
        if (__j <= __M - 1) ADD_TO_RESULT(1, __ind[__i]);
      }
	}
  }
#undef ADD_TO_RESULT

  /**
   * 生成一个麦克斯韦分布函数，以 @ref __c 为平均速度，@ref __s
   * 的平方为温度，@p __rho 为密度。
   */
  void MakeMaxwellian(const _T& __rho) {
	Assert(_Base::GetOrder() != (unsigned int)-1,
           std::length_error("Moment order insufficient"));
	Indices __ind;
	_Base::operator=((_T)0);
	_Base::operator()(__ind) = __rho;
  }

  /**
   * 指向起始元素的迭代器。
   */
  Iterator begin() {
	return Iterator(_Base::begin());
  }

  /**
   * 指向起始元素的 @p const 迭代器。
   */
  ConstIterator begin() const {
	return ConstIterator(_CBase::begin());
  }

  /**
   * 指向末尾元素的迭代器。
   */
  Iterator end() {
	return Iterator(_Base::end());
  }

  /**
   * 指向末尾元素的 @p const 迭代器。
   */
  ConstIterator end() const {
	return ConstIterator(_CBase::end());
  }

  /**
   * 指向 @p __m 阶矩起始元素的迭代器。
   */
  Iterator begin(unsigned int __m) {
	return Iterator(_Base::begin(__m));
  }

  /**
   * 指向 @p __m 阶矩起始元素的 @p const 迭代器。
   */
  ConstIterator begin(unsigned int __m) const {
	return ConstIterator(_CBase::begin(__m));
  }

  /**
   * 指向 @p __m 阶矩末尾元素的迭代器。
   */
  Iterator end(unsigned int __m) {
	return begin(__m + 1);
  }

  /**
   * 指向 @p __m 阶矩起始元素的 @p const 迭代器。
   */
  ConstIterator end(unsigned int __m) const {
	return begin(__m + 1);
  }

  /**
   * 对分布函数做 Givens 变换。设变换矩阵为 @p A，则
   * \f[
   *     A_{ii} = A_{jj} = \cos\theta, \quad
   *     A_{ij} = -A_{ji} = \sin\theta, \quad
   *     A_{kl} = \delta_{kl}, \qquad
   *     k,l \neq i,j.
   * \f]
   */
  void GivensTransform(int __i, int __j, _T __theta, const Distribution& __src)
  {
	// 把矩阵写出来方便取元素
	Matrix<_T> __A((IdentityMatrix<_T>(DIM)));
	_T __cos = cos(__theta);
	_T __sin = sin(__theta);
	__A(__i, __i) = __cos; __A(__j, __j) = __cos;
	__A(__i, __j) = __sin; __A(__j, __i) = -__sin;

	ConstIterator __it = __src.begin();
	ConstIterator __it_end = __src.end();

	// 旋转后的速度
	const Velocity& __v = __src.Center();
	Velocity __v1(__v);
	__v1[__i] = __v[__i] * __cos - __v[__j] * __sin;
	__v1[__j] = __v[__j] * __cos + __v[__i] * __sin;

	Reinit(__v1, __src.Scaling(), __src.GetOrder());

	unsigned int __n = &(*__it_end) - &(*__it);
	std::vector<_Base> __tmp(__n);

	for (unsigned int __n = 0; __it != __it_end; ++__it, __n++) {
      Indices __ind = __it.GetIndices();
      unsigned int __M = 0;
      for (unsigned int __d = 0; __d < DIM; __d++)
        __M += __ind[__d];

      unsigned int __N = __ind[__i] + __ind[__j];

      __tmp[__n].Reinit(__M);

      if (__M == 0) {
        __tmp[__n].front() = 1;
        _Base::front() = *__it;
      } else {
        for (unsigned int __x = 0; __x <= __N; ++ __x) {
          Indices __ind1 = __ind;
          __ind1[__i] = __x; __ind1[__j] = __N - __x;
          _T& __ele = __tmp[__n](__ind1);
          _T& __dst = (*this)(__ind1);

          unsigned int __d;
          for (__d = 0; __d < DIM; __d++)
            if (__ind1[__d] > 0) break;
          __ind1[__d] --;

          for (unsigned int __k = 0; __k < DIM; __k++) {
            if (__ind[__k] > 0) {
              __ind[__k] --;
              unsigned int __n = &__src(__ind) - &__src.front();
              __ind[__k] ++;

              __ele += __ind[__k] * __A(__k, __d)
                * __tmp[__n](__ind1) / (__ind1[__d] + 1);
            }
          }
          __dst += (*__it) * __ele;
        }
      }
	}
  }

  /**
   * 对分布函数做一正交变换，变换矩阵由 @p __A 给出。@p __src 为
   * 输入分布函数 f(v)，变换得到的结果 f(Av) 存储在当前分布函数
   * 中。用户需自行保证 @p __A 为一正交矩阵，程序只对 @p __A 是
   * 否为方阵进行检验。同时 @p __A 的行数 @p m 须小于或等于 @p
   * DIM，当小于时，可认为变换矩阵为
   *   \f[ \mathrm{diag}\{A, I_{D-m}\}, \f]
   * 其中 \f$D\f$ 即为 @p DIM。
   */
  void OrthTransform(const Matrix<_T>& __A, const Distribution& __src)
  {
	Matrix<_T> __B(__A);
	Assert(__A.m() == __A.n(), std::domain_error("The parameter is not a square matrix"));
	Assert(__A.m() <= DIM, std::length_error("Matrix too small"));

	int __D = std::min<int>(__A.m(), DIM);

	Distribution __tmp;
	const Distribution * __p0 = &__src;
	Distribution * __p1, * __p2;
	if (__D * (__D - 1) / 2 % 2)
	  __p1 = &__tmp, __p2 = this;
	else
	  __p1 = this, __p2 = &__tmp;

	for (unsigned int __i = 0; __i < __D - 1; __i++)
	  for (unsigned int __j = __i + 1; __j < __D; __j++) {
        _T __theta = _GivensParameter(__B(__i, __i), __B(__j, __i));
        _GivensTransform(__i, __j, __theta, __B);
        __p2->GivensTransform(__j, __i, __theta, *__p0);

        std::swap(__p1, __p2);
        __p0 = __p1;
	  }

	std::vector<bool> __neg(__D, false);
	for (unsigned int __d = 0; __d < __D; __d++)
	  if (__B(__d, __d) < 0) {
        __c[__d] = -__c[__d], __neg[__d] = true;
	  }

	for (Iterator __it = begin(), __it_end = end(); 
         __it != __it_end; ++__it)
	{
      const Indices& __ind = __it.GetIndices();
      for (unsigned int __d = 0; __d < __D; __d++)
        if (__ind[__d] % 2 && __neg[__d]) *__it = -(*__it);
	}
  }

  /**
   * 对分布函数 @p __src 做一镜面反射变换，其中 @p __normal 为镜
   * 面单位法向，结果保存在当前分布函数中。用户需自行保证 @p
   * __normal 的范数为一。
   */
  template <class _VECTOR>
    void SymmTransform(const _VECTOR& __normal, const Distribution& __src) {
	Matrix<_T> __A((IdentityMatrix<_T>(DIM)));
	for (unsigned int __i = 0; __i < DIM; __i++)
	  for (unsigned int __j = 0; __j < DIM; __j++)
		__A(__i, __j) -= 2 * __normal[__i] * __normal[__j];
	OrthTransform(__A, __src);
  }

  /**
   * 对分布函数 @p __src 做截断，截平面与第 @p __k 个坐标轴垂直
   * ，且第 @p __k 维坐标为 @p __u。当 @p __sign 为零或正时取正
   * 的部分，否则取负的部分。
   */
  void HalfSpaceCutOff(unsigned int __k, _T __u,
                       const Distribution& __src, int __sign = 1)
  {
	const Velocity& __c = __src.Center();
	const _T& __s = __src.Scaling();
	unsigned int __M = __src.GetOrder();

	// 计算 Hermite 多项式的值
	std::vector<_T> __He(__M + 1);
	_T __v_k = (__u - __c[__k]) / __s;
	__He[0] = 1; __He[1] = __v_k;
	for (unsigned int __i = 2; __i <= __M; __i++) {
      __He[__i] = __He[1] * __He[__i - 1] - (__i - 1) * __He[__i - 2];
	}

	// exp(-(xi_k - c_k)^2 / 2 theta) 的值
	_T __exp = exp(-__v_k * __v_k / 2);

	// 计算函数 K_{m,n} 的值
	std::vector<std::vector<_T> > __K(__M + 1, std::vector<_T>(__M + 1));
	__K[0][0] = 1 / sqrt(2 * M_PI);
	for (unsigned int __i = 1; __i <= __M; __i++) {
      __K[0][__i] = __He[__i - 1] * __exp / sqrt(2 * M_PI);
      __K[__i][0] = __K[__i - 1][0] / __i;
	}
	__K[0][0] = 0.5 * erfc(__v_k / M_SQRT2);

	for (unsigned int __i = 1; __i <= __M; __i++)
	  __K[__i][0] *= __He[__i - 1] * __exp;

	for (unsigned int __i = 1; __i <= __M; __i++)
	  for (unsigned int __j = 1; __j <= __M; __j++)
		__K[__i][__j] = __K[__i][0] * __He[__j] + __j * __K[__i - 1][__j - 1] / __i;

	// theta 的方幂
	std::vector<_T> __pow(2 * __M + 1);
	__pow[__M] = 1;
	for (unsigned int __i = 1; __i <= __M; __i++) {
      __pow[__M + __i] = __pow[__M + __i - 1] * __s;
      __pow[__M - __i] = __pow[__M - __i + 1] / __s;
	}

	// 计算半空间上的积分
	Reinit(__c, __s, __M);
	for (unsigned int __m = 0; __m <= __M; __m++) {
      for (Iterator __it = _Base::begin(__m),
             __it_end = _Base::begin(__m + 1);
           __it != __it_end; ++__it)
      {
        const Indices& __ind = __it.GetIndices();
        unsigned int __N = __M - (__m - __ind[__k]);

        Indices __ind1 = __ind;
        for (unsigned int __n = 0; __n <= __N; __n++) {
          __ind1[__k] = __n;
          *__it += __src(__ind1) * __pow[__ind[__k] + __M - __n] * __K[__ind[__k]][__n];
        }
      }
	}

	// 取负的部分时
	if (__sign < 0) {
      ConstIterator __it_src = __src.begin();
      for (Iterator __it = _Base::begin(), __it_end = _Base::end();
           __it != __it_end; ++__it, ++__it_src)
      {
        *__it = *__it_src - *__it;
      }
	}
  }

  /**
   * 求半空间上的麦克斯韦分布在任意一个空间中的表示。其中 @p
   * __dir 表示截面的法向，且截面总穿过麦克斯韦分布所表示的平均
   * 速度的位置。@p __pv 为麦克斯韦分布所对应的各原始变量，即密
   * 度、速度、温度，截面截得的两部分中，总选取 @p __dir 所指向
   * 的那一半。
   */
  template <class ARRAY, class VECTOR>
    void HalfMaxwellian(const ARRAY& __pv, const VECTOR& __dir) {
	// 构造一系列 Givens 变换把 __dir 变成 -e_1
	Vector<DIM, _T> __dir_bak(&__dir[0]);
	std::vector<_T> __pv_bak(&__pv[0], &__pv[DIM + 2]);
	std::vector<_T> __phi(DIM - 1);
	for (int __i = DIM - 1; __i > 0; __i--) {
      _T __rot = _GivensParameter(__dir_bak[__i - 1], __dir_bak[__i]);
      _T __cos = cos(__rot), __sin = sin(__rot);
      __phi[__i - 1] = __rot;
      __dir_bak[__i - 1] = __dir_bak[__i - 1] * __cos + __dir_bak[__i] * __sin;

      _T __u_i = -__pv_bak[__i] * __sin + __pv_bak[__i + 1] * __cos;
      __pv_bak[__i] = __pv_bak[__i] * __cos + __pv_bak[__i + 1] * __sin;
      __pv_bak[__i + 1] = __u_i;

      _T __c_i = -__c[__i - 1] * __sin + __c[__i] * __cos;
      __c[__i - 1] = __c[__i - 1] * __cos + __c[__i] * __sin;
      __c[__i] = __c_i;
	}
	if (__dir_bak[0] > 0) {
      __c[0] = -__c[0];
      __pv_bak[1] = -__pv_bak[1];
	}

	// 开始求分布函数
	unsigned int __M = _Base::GetOrder();
	std::vector<double> __K(__M, 0);

	_T __theta = __s * __s;
	__K[0] = sqrt(__pv_bak[DIM + 1] / (2 * M_PI)) *
      exp(-(__c[0] - __pv_bak[1]) * (__c[0] - __pv_bak[1]) / (2 * __pv_bak[DIM + 1]));
	for (unsigned int __i = 2; __i < __M; __i += 2)
	  __K[__i] = -__K[__i - 2] * (__i - 1) / ((__i + 1) * __i) * __theta;

	std::vector<double> __H[DIM];
	__H[0].resize(__M + 1);
	__H[0][0] = 0.5 * erfc((__pv_bak[1] - __c[0]) / sqrt(2 * __pv[DIM + 1]));
	__H[0][1] = (__pv[1] - __c[0]) * __H[0][0] - __K[0];

	for (unsigned int __i = 2; __i <= __M; __i++) {
      __H[0][__i] = ((__pv_bak[1] - __c[0]) * __H[0][__i - 1] +
                     (__pv_bak[DIM + 1] - __theta) * __H[0][__i - 2]) / __i - __K[__i - 1];
	}

	for (unsigned int __k = 1; __k < DIM; __k++) {
      __H[__k].resize(__M + 1);
      __H[__k][0] = 1; __H[__k][1] = __pv_bak[__k + 1] - __c[__k];
      for (unsigned int __i = 2; __i <= __M; __i++) 
        __H[__k][__i] = ((__pv_bak[__k + 1] - __c[__k]) * __H[__k][__i - 1] +
                         (__pv_bak[4] - __theta) * __H[__k][__i - 2]) / __i;
	}

	Distribution __tmp[DIM];
	__tmp[0].Reinit(*this);
	Iterator __it = __tmp[0].begin();
	Iterator __it_end = __tmp[0].end();
	for (; __it != __it_end; ++__it) {
      Indices __ind = __it.GetIndices();
      *__it = __pv_bak[0];
      for (unsigned int __k = 0; __k < DIM; __k++) *__it *= __H[__k][__ind[__k]];
	}

	// 最后通过旋转变换回来
	for (unsigned int __i = 0; __i < DIM - 1; __i++)
	  __tmp[__i + 1].GivensTransform(__i, __i + 1, __phi[__i], __tmp[__i]);
	*this = __tmp[DIM - 1];
	if (__dir_bak[0] > 0) {
      __c[0] = -__c[0];
      Iterator __it = begin();
      Iterator __it_end = end();
      for (; __it != __it_end; ++__it) {
        Indices __ind = __it.GetIndices();
        if (__ind[0] % 2 == 1) *__it = -(*__it);
      }
	}
  }

  /**
   * 获取分布函数在 @p __p 点处的值。
   */
  template <class _VECTOR> value_type Value(const _VECTOR& __p) const {
	  _Base __basis_value(_Base::GetOrder());
	  Velocity __v;

	  typename _Base::Iterator __it = __basis_value.begin();
	  *__it = 1;
	  for (unsigned int __i = 0; __i < dim; __i++) {
		  __v[__i] = (__p[__i] - __c[__i]) / __s;
		  *__it *= exp(-__v[__i] * __v[__i] / 2) * (M_2_SQRTPI / (2 * M_SQRT2)) / __s;
	  }

	  value_type __sqr_s = __s * __s;
	  typename _Base::Iterator __it_end = __basis_value.end();
	  for (++__it; __it != __it_end; ++__it) {
		  Indices __ind(__it.GetIndices());
		  unsigned int __i = 0;
		  while (__ind[__i] == 0) __i++;
		  __ind[__i]--;
		  *__it = __basis_value(__ind) * __v[__i] / __s;

		  if (__ind[__i] != 0) {
			  __ind[__i]--;
			  *__it -= (__ind[__i] + 1) * __basis_value(__ind) / __sqr_s;
		  }
	  }

	  return std::inner_product(__basis_value.begin(), __it_end, begin(), value_type(0));
  }
};

template <int DIM, class _T> const _T * Distribution<DIM, _T, 0>::delta = NULL;

/**
 * 投影操作。将分布函数 @p __src 投影到 @p __dst 所设置的空间中，参数
 * @p __m 若为正，则表示仅对前 @p __m 阶矩做投影。
 */
template <class _DISTRIBUTION>
void Project(const _DISTRIBUTION& __src,
             _DISTRIBUTION& __dst, int __m = -1)
{
  __dst.Project(__src, __m);
}

/**
 * 将分布函数 @p __src 乘上第 @p __i 维的速度，结果存在 @p __dst 中。
 */
template <class _DISTRIBUTION>
void MulVelocity(unsigned int __i,
                 const _DISTRIBUTION& __src, _DISTRIBUTION& __dst)
{
  __src.MulVelocity(__i, __dst);
}

NRXX_NAMESPACE_CLOSE

#endif // __Distribution_h_
