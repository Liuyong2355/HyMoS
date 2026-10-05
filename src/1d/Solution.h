/**
 * @file   Solution.h
 * @author HU Zhicheng <huzhicheng1986@gmail.com>
 * @date   Thu Dec 26 11:02:42 2013
 * 
 * @brief  封装问题解的类型
 * 
 * 感觉 NRxx 中 Grid 类和 Solution 类封装有些复杂, 等到空间高维问题再考
 * 虑用 NRxx 中的程序, 现在先给自己封装一个.
 */


#ifndef __SSP_SOLUTION_H__  // 不知道与 NRxx 中 Solution.h 定义相同的宏
#define __SSP_SOLUTION_H__  // 会不会出问题, 还是再加个标识以示区别吧

#include <vector>
/// 为了封装的类支持序列化操作 (这个概念还不甚明白)
#include <boost/archive/text_oarchive.hpp>
#include <boost/archive/text_iarchive.hpp>
/// 派生类的序列化需要用到
#include <boost/serialization/base_object.hpp>

template <typename DISTRIBUTION>
class Solution : public std::vector<DISTRIBUTION> {
public:
  typedef std::vector<DISTRIBUTION> BASE;
  typedef DISTRIBUTION dis_t;

public:
  Solution(){ };
  Solution( unsigned int n ) : BASE(n) { };
  Solution( unsigned int n, unsigned int ORDER ) : BASE(n, DISTRIBUTION(ORDER)) { };
  /** 
   * @brief 用传入的解重初始化
   * 
   * @param _sol 从这里获取重初始化的参数系
   */
  void Reinit( Solution& _sol ) {
    this->resize( _sol.size() );
    for( unsigned int i=0; i<this->size(); ++i ) {
      (*this)[i].Reinit( _sol[i] );
    }
  }

  Solution& operator+=( const Solution& sol ) {
    for( unsigned int i = 0; i<this->size(); ++i ){
      (*this)[i] += sol[i];
    }

    return *this;
  }
  Solution& operator-=( const Solution& sol ) {
    for( unsigned int i = 0; i<this->size(); ++i ){
      (*this)[i] -= sol[i];
    }

    return *this;
  }
  /**
   * 数乘
   * 
   */
  template <typename NUMBER>
  Solution& operator*=( const NUMBER& s ) {
    for( unsigned int i = 0; i<this->size(); ++i ){
      (*this)[i] *= s;
    }

    return *this;
  }
  /**
   * 数除
   * 
   */
  template <typename NUMBER>
  Solution& operator/=( const NUMBER& s ) {
    (*this) *= 1.0 / s;

    return *this;
  }


private:
  //为了能让串行化类库能够访问私有成员，所以要声明一个友元类
  friend class boost::serialization::access;
  //串行化的函数，这一个函数完成对象的保存与恢复
  template <class Archive>
  void serialize(Archive & ar, const unsigned int version) {
    /**
     * 若 a 是仅有的数据成员, 则 "ar & a;" 即可. 对于派生类, 要保存父类
     * 的数据, 不能直接调用父类的serialize函数, 而要如下操作. 注意要添
     * 加头文件 #include <boost/serialization/base_object.hpp>. 还有直
     * 接使用 "ar & *this;" 编译并未报错, 而是运行时段错误.
     */
    ar & boost::serialization::base_object<BASE>(*this);
  }
};

#endif // __SSP_SOLUTION_H__

/**
 * end of file
 * 
 */
