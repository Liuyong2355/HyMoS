/**
 * @file   Mesh1D.h
 * @author HU Zhicheng <huzhicheng1986@gmail.com>
 * @date   Fri Dec 21 10:58:53 2012
 * 
 * @brief  一维网格类
 * 
 * 为非均匀网格的使用, 还是如 YWFan 那样把网格封装一下
 */
#ifndef __MESH_1D_h__
#define __MESH_1D_h__

#include<vector>
#include<string>
#include<fstream>
#include <memory>

#include "ConfigFile.h"

template <class _T>
class Mesh1D : public std::vector<_T>{
public:
  typedef _T point_t;
  typedef std::vector<point_t> _Base;

public:
  point_t x_l;                  /**< 区间左边界 */
  point_t x_r;                  /**< 区域右边界 */

  static const int dim = 1;     /**< 空间维数 */

public:
  Mesh1D(){};
  Mesh1D( const std::string _str ){
    ReadMesh(_str);
  }
  /** 
   * @brief 构造函数
   * 
   * @param x0 左边界
   * @param x1 右边界, 当然要求 x1 > x0
   * @param N 网格单元数
   */
  Mesh1D( const point_t& x0, const point_t& x1, const size_t N ) : x_l(x0), x_r(x1) {
    Reinit( N );
  }
  Mesh1D( const ConfigFile& cf ) {
    read_config( cf );
  }
  Mesh1D( const std::shared_ptr<Mesh1D>& p_old_mesh, const unsigned int N ) 
  : x_l(p_old_mesh->x_l), x_r(p_old_mesh->x_r) {
    Reinit( N );
  }
  /// 此构造函数要求读入流数据由 dump_data() 导出
  Mesh1D( std::istream& is ) {
    unsigned int N;
    is.read((char *)&N, sizeof(unsigned int));
    this->resize( N );
    for (size_t i = 0; i < this->size(); ++ i)
      is.read((char *)&((*this)[i]), sizeof(double));
    x_l = *(this->begin());
    x_r = *(this->end()-1);
  }

  /** 
   * @brief 从文件中读入一个网格
   * 
   * @param _str 文件名
   */
  void ReadMesh( const std::string _str ){
    std::ifstream io(_str.c_str());
    unsigned int N; /// 注意 N+1 才是网格点数目
    io >> N;
    this->resize(N+1);
    for( size_t i=0; i<=N; ++i )
      io >> (*this)[i];
    x_l = *(this->begin());
    x_r = *(this->end()-1);
  }

  void read_config( const ConfigFile& cf ){
    int m_t = (int)cf.Value("Mesh", "Type");
    if( m_t == 0 ){
      int Nx = (int)cf.Value("Mesh","Nx");

      x_l = (double)cf.Value("Mesh","L0");
      x_r = x_l + (double)cf.Value("Mesh","Len");

      Reinit( Nx ); 
    }
    else { /// 读入网格程序暂放于此, 以后应该将程序read_config 和
           /// initialize两部分完全分开.
      std::string file = (std::string)cf.Value("Mesh", "File" );
      ReadMesh( file );
    }
  }

  /** 
   * @brief 初始化为均匀网格
   *
   * 要求已经知晓区域左右边界
   * 
   * @param _N 网格单元数
   */
  void Reinit( const size_t _N ){
    this->resize( _N+1 );
    _T h = (x_r-x_l) / _N;
    for( size_t i=0; i<_N; ++i )
      (*this)[i] = x_l + i * h;
    (*this)[_N] = x_r;
  }
  /** 
   * @brief 初始化为均匀网格
   * 
   * @param x0 左边界
   * @param x1 右边界, 要求 x1 > x0
   * @param _N 网格单元数
   */
  void Reinit( const point_t x0, const point_t x1, const size_t _N ){
    x_l = x0;
    x_r = x1;
    Reinit( _N );
  }

  void dump_data( std::ostream& os ) const {
    unsigned int N = this->size();
    os.write((const char *)&(N), sizeof(int));
    for (size_t i = 0; i < this->size(); ++ i)
      os.write((const char *)&((*this)[i]), sizeof(double));  
  }

  /** 
   * @brief 计算单元的重心
   * 
   * @param i 单元指标
   * 
   * @return 单元重心
   */
  point_t Center( const size_t i ) const {
    return ((*this)[i]+(*this)[i+1]) * 0.5;
  }
  /** 
   * @brief 计算单元的长度 (面积)
   * 
   * @param i 单元指标
   * 
   * @return 单元长度
   */
  _T Length( const size_t i ) const {
    return (*this)[i+1] - (*this)[i];
  }
  /** 
   * @brief 计算网格单元数
   * 
   * @return 返回网格单元数目, 比网格点数少一.
   */
  size_t n_ele() const {
    return this->size()-1;
  }

  /** 
   * @brief 粗化网格
   *
   * 将相邻两个单元拼成新网格的一个单元. 在多重网格方法中应用
   * 
   * @param _mesh 细网格, 要求细网格单元数是偶数 (此函数并不检测)
   */
  void CoarseFrom( const Mesh1D &_mesh ){
    if( _mesh.n_ele() % 2 != 0 )
    {/**
      * 虽然有这个提示, 由于此函数以及 MultiGrid 限制函数等以粗网格单
      * 元做循环, 多重网格方法还是可以运行下去的, 虽然多重网格回代的时
      * 候可能不会更改最后一个细网格单元上的解 (这个解的校正就要靠光滑
      * 步了).
      */
      std::cerr << "The mesh could not be coarse: the number of the elements is odd!" << std::endl;
    }
    this->resize( _mesh.n_ele() / 2 + 1 );
    typename _Base::iterator it=this->begin();
    const typename _Base::iterator it_end=this->end();
    typename _Base::const_iterator it_fine_mesh = _mesh.begin();
    for( ; it != it_end; ++it, it_fine_mesh += 2 ){
      *it = *it_fine_mesh;
    }
    x_l = *(this->begin());
    x_r = *(this->end()-1);
  }

  /** 
   * @brief 全局加密网格一次
   * 
   * 将 @p _mesh 中的每个单元都均分为二. 当 @p _mesh 本身是非均匀网格时,
   * 这样得到的网格可能光滑性不是很好, 但有更好又简单易实现的方式么?
   *
   * @param _mesh
   */
  void RefineFrom( const Mesh1D &_mesh ) {
    this->resize( _mesh.n_ele() * 2 + 1 );
    typename _Base::const_iterator it_coarse_mesh = _mesh->begin();
    const typename _Base::const_iterator the_coarse_end = _mesh.end()-1;
    typename _Base::iterator it = this->begin();
    for( ; it_coarse_mesh != the_coarse_end; ++ it_coarse_mesh, it += 2 ) {
      *it = *it_coarse_mesh;
      *(it+1) = 0.5 * (*it_coarse_mesh + *(it_coarse_mesh+1));
    }
    *(this->end()-1) = _mesh.x_r;
    x_l = _mesh.x_l;
    x_r = _mesh.x_r;
  }

  double DomainSize() const {
    return x_r - x_l;
  }
};

#endif // __MESH_1D_h__
