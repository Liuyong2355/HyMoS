/**
 * @file Collision.h
 * @brief BGK类型碰撞模型
 * @author Guanghan Li NUAA
 * @version 1.0
 * @date 2022-10-20
 */

# ifndef __COLLISION__H__
# define __COLLISION__H__

enum ColModel {BGK, ESBGK, Shakhov};

template <ColModel CM, typename SOLUTION, typename PROBLEM>
class Collision{
	public:
	    typedef SOLUTION sol_t;
	    typedef PROBLEM Problem;
	    typedef typename sol_t::dis_t dis_t;
	public:
		static void collision(size_t nx, size_t ny,
				const sol_t& sol, sol_t& res);
};

template <typename SOLUTION, typename PROBLEM>
class Collision<BGK, SOLUTION, PROBLEM>{
	public:
	    typedef SOLUTION sol_t;
	    typedef PROBLEM Problem;
	    typedef typename sol_t::dis_t dis_t;
	public:
		static void collision(size_t nx, size_t ny,
				const sol_t& sol, sol_t& res){
		  	const dis_t& dis_sol = sol[nx][ny];
		  	double tau = Problem::getTau(dis_sol);
			double nu = 1.0 / tau;
		  	dis_t& dis_res = res[nx][ny];
		  	typename dis_t::Iterator the_dis = dis_res.begin(2);
		  	typename dis_t::Iterator end_dis = dis_res.end();
		  	typename dis_t::ConstIterator src_dis = dis_sol.begin(2);
		  	for( ; the_dis != end_dis; ++the_dis, ++src_dis) {
		  		 *the_dis += -nu * (*src_dis);
		  	}
		}

};

template <typename SOLUTION, typename PROBLEM>
class Collision<ESBGK, SOLUTION, PROBLEM>{
	public:
	    typedef SOLUTION sol_t;
	    typedef PROBLEM Problem;
	    typedef typename sol_t::dis_t dis_t;
	public:
		static void collision(size_t nx, size_t ny,
		  						   const sol_t& sol, sol_t& res){
		  	const dis_t& dis_sol = sol[nx][ny];
		  	double tau = Problem::getTau(dis_sol);
			//ES-BGK碰撞项需要乘以Pr
		  	double nu = Problem::Pr / tau;
		  	dis_t& dis_res = res[nx][ny];
		  	typename dis_t::Iterator the_dis = dis_res.begin(2);
		  	typename dis_t::Iterator end_dis = dis_res.end();
		  	typename dis_t::ConstIterator src_dis = dis_sol.begin(2);
		  	for( ; the_dis != end_dis; ++the_dis, ++src_dis) {
		  		 *the_dis +=  -nu * (*src_dis);
		  	}

		  	dis_t dis_es;
		  	GetEquilibriumDistribution(dis_es, dis_sol, Problem::Pr);
		  	dis_es.front() = 0.0;
		  	dis_res.Add(nu, dis_es);
		}

		void GetEquilibriumDistribution(dis_t& dis_es, const dis_t& src_dis,
		  										 const double _Pr){
		  	/// 取 \f$ \sigma_{ij} = (1+ \delta_{ij}) f_{e_i+e_j}. \f$
		  	const int dim = dis_t::dim;
		  	double sigma[dim][dim];
		  	for( int i=0; i<dim; ++i )
		  	{
		  		 for( int j=0; j<dim; ++j )
		  		 {
		  			  typename dis_t::Indices ind; // 初始化 ind 是 0
		  			  ind[i] += 1;   ind[j] += 1;
		  			  sigma[i][j] = src_dis( ind );
		  			  /**
		  				* 这样的取法效率如何, 或者直接根据内部排序对应关系而使用迭代器?
		  				* 根据矩系数的排序, 二阶矩对应的分别是 \$ 0.5 \sigma_{11},
		  				* \sigma_{12}, \sigma_{13}, 0.5 \sigma_{22}, \sigma_{23}, 0.5
		  				* \sigma_{33} \$
		  				* 
		  				*/
		  		 }
		  		 sigma[i][i] *= 2;
		  	}

		  	// /**
		  	//  * 引入这个判断仅是为了确认传入的分布函数是在标准展开下的. 然后刚测
		  	//  * 试无视扰动分布函数导致标准展开参数系的改变, 仍然使用快速算法进行
		  	//  * 计算, 于是这里的这个判断需要注释掉了.
		  	//  * 
		  	//  */
		  	// if (fabs(sigma[0][0]+sigma[1][1]+sigma[2][2] )> 1e-14 )
		  	// {
		  	//   std::cerr <<"sigma error"<<std::endl;
		  	//   getchar();
		  	// }

		  	/// 计算 \f$ f_{ES} \f$ 的矩系数
		  	const double _rho = src_dis.front();
		  	const double factor = ( 1 - 1./_Pr ) / _rho;

		  	dis_es.Reinit( src_dis );
		  	dis_es.MakeMaxwellian( _rho );
		  	typename dis_t::Iterator the_dis = dis_es.begin(2);
		  	typename dis_t::Iterator end_dis = dis_es.end();
		  	for( ; the_dis != end_dis; ++ the_dis )
		  	{
		  		 typename dis_t::Indices ind = the_dis.GetIndices();
		  		 /// 若 \f$ |\alpha| \f$ 是奇数, 系数为 0.
		  		 int _sum = 0;
		  		 for( int i=0; i<dim; ++i )
		  		 _sum += ind[i];
		  		 if( _sum % 2 == 1 )   continue;

		  		 /// 确定指标 i, 使得 ind[i] > 0.
		  		 for( int i=0; i<dim; ++i )
		  		 {
		  			  if( ind[i] == 0 )   continue;
		   
		  			  double _G_alpha = 0.0;      
		  			  ind[i] -= 1;   // 此时 ind[i] >= 0;
		  			  for( int k=0; k<dim; ++k )
		  			  {
		  					if( ind[k] == 0 )   continue;  // 根据递推式需要跳过指标为负的项
		  					ind[k] -= 1;
		  					_G_alpha += sigma[i][k] * dis_es( ind );
		  					ind[k] += 1;
		  			  }
		  			  ind[i] += 1; // 回到当前指标
		  			  *the_dis = _G_alpha * factor / ind[i];
		  			  // 只需要一个指标 i 就可以了.
		  			  break;
		  		 }
		  	}
		}
};

template <typename SOLUTION, typename PROBLEM>
class Collision<Shakhov, SOLUTION, PROBLEM>{
	public:
	    typedef SOLUTION sol_t;
	    typedef PROBLEM Problem;
	    typedef typename sol_t::dis_t dis_t;
	public:
		static void collision(size_t nx, size_t ny,
		  						   const sol_t& sol, sol_t& res){
		  	const dis_t& dis_sol = sol[nx][ny];
		  	double tau = Problem::getTau(dis_sol);
			double nu = 1.0 / tau;
		  	dis_t& dis_res = res[nx][ny];
		  	typename dis_t::Iterator the_dis = dis_res.begin(2);
		  	typename dis_t::Iterator end_dis = dis_res.end();
		  	typename dis_t::ConstIterator src_dis = dis_sol.begin(2);
		  	for( ; the_dis != end_dis; ++the_dis, ++src_dis) {
		  		 *the_dis +=  -nu * (*src_dis);
		  	}

		  	const int dim = dis_t::dim;
		  	double heat_flux[dim];
		  	dis_sol.HeatFlux( heat_flux );
			//这个地方莫名其妙加了static
		  	double factor = nu * ( 1 - Problem::Pr ) / (dim+2);
		  	typename dis_t::Indices ind; // 初始化 ind 应该是 0
		  	for( int i=0; i<dim; ++i )
		  	{
		  		 ind[i] = 1;
		  		 for( int j=0; j<dim; ++j )
		  		 {
		  			  ind[j] += 2;
		  			  dis_res( ind ) += factor * heat_flux[i];
		  			  ind[j] -= 2;
		  		 }
		  		 ind[i] = 0;
		  	}
		}
};

# endif
