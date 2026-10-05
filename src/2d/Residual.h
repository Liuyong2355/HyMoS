/**
 * @file Residual.h
 * @brief 残量类
 * @author Guanghan Li NUAA
 * @version 1.0
 * @date 2022-10-12
 */

# ifndef __RESIDUAL__H__
# define __RESIDUAL__H__

# include "Collision.h"
# include "config.h"
# include <Eigen/Core>
# include <Eigen/Dense>

template <typename SOLUTION, ColModel CM, typename PROBLEM>
class Residual{
	public:
		typedef SOLUTION sol_t;
	    typedef PROBLEM Problem;
	public:
	    typedef typename sol_t::dis_t dis_t;
	    typedef typename sol_t::mesh_t mesh_t; 
	    typedef typename dis_t::value_type value_t;
	    typedef typename dis_t::Indices index_t;
	    typedef typename dis_t::Velocity velocity_t;
	    typedef typename Eigen::VectorXd Euler_sol_t;
		typedef typename std::vector<std::vector<Euler_sol_t>> Euler_sol_vec_t;	  
	public:
		static void getResidualL2Norm(value_t& l2_norm, const sol_t& sol,
		  		 const sol_t& rhs, const mesh_t& mesh){
		  	sol_t res;	
		  	getResidual(sol, rhs, res, mesh);
		  	getResidualL2Norm(l2_norm, res, mesh);
		}

		static void getResidualL2Norm(value_t& l2_norm, 
		  		 const sol_t& res, const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	l2_norm = 0;
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		    double local_norm = L2NormDis(res[nx][ny]);
		  		    l2_norm += local_norm * local_norm * mesh.Area(nx, ny);
		  		}
		  	}
		  	l2_norm = sqrt(l2_norm / mesh.DomainSize());
		}
		
		static double L2NormDis(const dis_t& dis){
		  	dis_t normalized_dis( dis.GetOrder() );
		  	NormalizedDistribution( normalized_dis, dis );
		  	double l2_norm = 0.0;
		  	typename dis_t::Iterator it_dis = normalized_dis.begin();
		  	typename dis_t::Iterator the_end = normalized_dis.end();
		  	for( ; it_dis != the_end; ++it_dis )
		  	{
		  		 l2_norm += (*it_dis) * (*it_dis);
		  	}
		  	// 归一化的时候只是做到了量纲统一, 所以要得到真正的二范数, 还要乘一个
		  	// 常数.
		  	static int dim = dis_t::dim;
		  	static double pi_c = pow( 2*M_PI, - 0.25*dim );

		  	return sqrt(l2_norm) * pi_c * pow( normalized_dis.Scaling(), -dim );
		}

		static void NormalizedDistribution(dis_t& normalized_dis, const dis_t& dis){
		  	static int dim = dis_t::dim;
		  	double one_over_s = 1. / dis.Scaling();
		  	double coe = 1;
		  	unsigned int old_order = 0;

		  	normalized_dis.Reinit(dis.Center(), dis.Scaling());
		  	typename dis_t::ConstIterator it_dis = dis.begin();
		  	//typename dis_t::ConstIterator it_end = dis.end();
		  	typename dis_t::Iterator it_normalized_dis = normalized_dis.begin();
		  	const typename dis_t::Iterator the_end = normalized_dis.end();
		  	for( ; it_normalized_dis != the_end; ++ it_dis, ++ it_normalized_dis) {
		  		unsigned int order = it_dis.Order();
		  		if ( order > old_order) {
		  		     coe *= one_over_s;
		  		     old_order = order;
		  		}
		  		double local_value = coe;
		  		const typename dis_t::Indices& ind = it_dis.GetIndices();
		  		for( int i = 0; i < dim; ++i) {
		  		     local_value *= SqrtFactorialNumber<int, double> (ind[i]);
		  		}
		  		*it_normalized_dis = *it_dis * local_value;
		  	}
		}

		template <typename INT, typename DOUBLE> 
		static double SqrtFactorialNumber(const INT i ){
		  	static std::map<INT,DOUBLE> sqrt_factorial;
		  	typename std::map<INT,DOUBLE>::iterator it=sqrt_factorial.find(i);
		  	if ( it != sqrt_factorial.end() )
		  		 return it->second;

		  	if( i > 1 )
		  	{
		  	sqrt_factorial.insert( typename std::map<INT,DOUBLE>::
		  			  value_type( i, sqrt(i) * SqrtFactorialNumber<INT, DOUBLE>( i-1 )) );
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

		static void getResidualPL2Norm( value_t& pl2_norm,
				const sol_t& sol, const sol_t& rhs,
				const mesh_t& mesh){
			sol_t res;
			getResidual(sol, rhs, res, mesh);
			getResidualPL2Norm(pl2_norm, res, mesh);
		}

		static void getResidualPL2Norm(value_t& pl2_norm,
				const sol_t& res, 
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
			pl2_norm = 0.0;
			for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
					double local_norm = PL2NormDis(res[nx][ny]);
					//std::cout << local_norm << " ";
			  		pl2_norm += local_norm * local_norm * mesh.Area(nx,ny);
		  		}
		  	}
			/**
			 * 对这全局 l2_norm 做个归一化, 或者说是求平均的 l2_norm.
			 * 
			 */    
			pl2_norm = sqrt( pl2_norm / mesh.DomainSize() );
		}	

		static double PL2NormDis( const dis_t& dis ) {
			// 下面这个构造函数参数必须是 unsigned int
			dis_t normalized_dis( std::min( 3u, dis.GetOrder() ) );
			NormalizedDistribution( normalized_dis, dis );
			return L2NormSDis( normalized_dis );
		}

		static double L2NormSDis( const dis_t& normalized_dis ) {
			double l2_norm = 0.0;
			typename dis_t::ConstIterator it_dis = normalized_dis.begin();
			typename dis_t::ConstIterator the_end = normalized_dis.end();
			for( ; it_dis != the_end; ++it_dis )
			{
			  l2_norm += (*it_dis) * (*it_dis);
			}
			// 归一化的时候只是做到了量纲统一, 所以要得到真正的二范数, 还要乘一个
			// 常数.
			static int dim = dis_t::dim;
			static double pi_c = pow( 2*M_PI, - 0.25*dim );

			return sqrt(l2_norm) * pi_c * pow( normalized_dis.Scaling(), -dim );  
		}

		static void getResidual(const sol_t& sol,
				const sol_t& rhs, 
				sol_t& res,
				const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	res.Reinit(sol);
		  	#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		    getResidual(sol, rhs, mesh, res, nx, ny);
		  		}
		  	}
		}

		static void getResidual(const sol_t& sol,
				const sol_t& rhs,
				const mesh_t& mesh,
				sol_t& res,
				size_t nx, size_t ny){
			Project(rhs[nx][ny], res[nx][ny]);
			convection(sol, mesh, res, nx, ny);
			Collision<CM, SOLUTION, PROBLEM>::collision(nx, ny, sol, res);
		}

		static void getResidual(const sol_t& sol,
				sol_t& res,
				const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	res.Reinit(sol);
		  	#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		    convection(sol, mesh, res, nx, ny);
		  		    Collision<CM, SOLUTION, PROBLEM>::collision(nx, ny, sol, res);
		  		}
		  	}
		}
        
		static void convection(const sol_t& sol,
				sol_t& res,
		  		const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		    convection(sol, mesh, res, nx, ny);
		  		}
		  	}
		}

		static void convection(const sol_t& sol,
				const mesh_t& mesh,
				sol_t& res,
				size_t nx, size_t ny){
			#ifdef RECONSTRUCT
			secondOrderDisPri(sol, mesh, res, nx, ny);
			//secondOrderDisCon(sol, mesh, res, nx, ny);
			#else	
			firstOrderDis(sol, mesh, res, nx, ny);
			#endif
		}
        
		static void firstOrderDis(const sol_t& sol,
		  		const mesh_t& mesh,
				sol_t& res, 
				size_t nx, size_t ny){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
		  	dis_t dis_l, dis_m, dis_r, flux_l2r, flux_r2l;
		  	std::vector<double> normal{1,0,0};
		  	double dx = mesh.dx(nx, ny);
		  	dis_t tmp(res[nx][ny].Center(), 
		  			res[nx][ny].Scaling(),
		  			res[nx][ny].GetOrder());
		  	if(nx == 0){
		  		PROBLEM::leftBoundaryCondition(sol[nx][ny], dis_l);
		  		dis_m = sol[nx][ny];
				dis_r = sol[nx+1][ny];
		  	}else if(nx == (Nx-1)){
		  		PROBLEM::rightBoundaryCondition(sol[nx][ny], dis_r);
				dis_l = sol[nx-1][ny];
		  		dis_m = sol[nx][ny];
		  	}else{
		  		dis_l = sol[nx-1][ny];
		  		dis_m = sol[nx][ny];
				dis_r = sol[nx+1][ny];
		  	}

			HLLFlux(dis_l, dis_m, normal, flux_l2r, flux_r2l);
		  	Project(flux_r2l, tmp);
		  	res[nx][ny].Add(-1.0/dx, tmp);

		  	HLLFlux(dis_m, dis_r, normal, flux_l2r, flux_r2l);
		  	Project(flux_l2r, tmp);
		  	res[nx][ny].Add(-1.0/dx, tmp);	

			normal = {0,1,0};
		  	double dy = mesh.dy(nx, ny);
			if(ny == 0){
		  		PROBLEM::bottomBoundaryCondition(sol[nx][ny], dis_l);
		  		dis_m = sol[nx][ny];
				dis_r = sol[nx][ny+1];
		  	}else if(ny == (Ny-1)){
		  		PROBLEM::topBoundaryCondition(sol[nx][ny], dis_r);
		  		dis_l = sol[nx][ny-1];
				dis_m = sol[nx][ny];
		  	}else{
				dis_l = sol[nx][ny-1];
				dis_m = sol[nx][ny];
				dis_r = sol[nx][ny+1];
		  	}

			HLLFlux(dis_l, dis_m, normal, flux_l2r, flux_r2l);
		  	Project(flux_r2l, tmp);
		  	res[nx][ny].Add(-1.0/dy, tmp);

		  	HLLFlux(dis_m, dis_r, normal, flux_l2r, flux_r2l);
		  	Project(flux_l2r, tmp);
		  	res[nx][ny].Add(-1.0/dy, tmp);
		}

		//对展开参数和展开系数线性重构	
		static void secondOrderDisPri(const sol_t& sol,
		  		const mesh_t& mesh,
				sol_t& res, 
				size_t nx, size_t ny){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
		  	dis_t re_lr, re_ml, re_mr, re_rl, dummy, flux_l2r, flux_r2l, flux;
			dis_t p_dis_lr, p_dis_ml, p_dis_mr, p_dis_rl;
		  	std::vector<double> normal{1.0,0,0};
		  	double dx = mesh.dx(nx, ny);
		  	dis_t tmp(res[nx][ny].Center(), 
		  			res[nx][ny].Scaling(),
		  			res[nx][ny].GetOrder());
		  	if(nx == 0){
				reconHorPri(sol, mesh, nx, ny, re_ml, re_mr);
				reconHorPri(sol, mesh, nx+1, ny, re_rl, dummy);
		  		PROBLEM::leftBoundaryCondition(re_ml, re_lr);
		  	}else if(nx == (Nx-1)){
				reconHorPri(sol, mesh, nx-1, ny, dummy, re_lr);
				reconHorPri(sol, mesh, nx, ny, re_ml, re_mr);
		  		PROBLEM::rightBoundaryCondition(re_mr, re_rl);
		  	}else{
				reconHorPri(sol, mesh, nx-1, ny, dummy, re_lr);
				reconHorPri(sol, mesh, nx, ny, re_ml, re_mr);
				reconHorPri(sol, mesh, nx+1, ny, re_rl, dummy);
		  	}

		 	p_dis_lr.Reinit( sol[nx][ny] ); Project( re_lr, p_dis_lr );
			p_dis_ml.Reinit( sol[nx][ny] ); Project( re_ml, p_dis_ml );
			p_dis_mr.Reinit( sol[nx][ny] ); Project( re_mr, p_dis_mr );
			p_dis_rl.Reinit( sol[nx][ny] ); Project( re_rl, p_dis_rl );
			
			NumericalFluxHor(flux, p_dis_lr, p_dis_ml);
		  	res[nx][ny].Add(1.0/dx, flux);
			  
			NumericalFluxHor(flux, p_dis_mr, p_dis_rl);
		  	res[nx][ny].Add(-1.0/dx, flux);
				
			normal = {0,1,0};
		  	double dy = mesh.dy(nx, ny);
			if(ny == 0){
				reconVerPri(sol, mesh, nx, ny, re_ml, re_mr);
				reconVerPri(sol, mesh, nx, ny+1, re_rl, dummy);
		  		PROBLEM::bottomBoundaryCondition(re_ml, re_lr);
		  	}else if(ny == (Ny-1)){
				reconVerPri(sol, mesh, nx, ny-1, dummy, re_lr);
				reconVerPri(sol, mesh, nx, ny, re_ml, re_mr);
		  		PROBLEM::topBoundaryCondition(re_mr, re_rl);
		  	}else{
				reconVerPri(sol, mesh, nx, ny-1, dummy, re_lr);
				reconVerPri(sol, mesh, nx, ny, re_ml, re_mr);
				reconVerPri(sol, mesh, nx, ny+1, re_rl, dummy);
		  	}

			p_dis_lr.Reinit( sol[nx][ny] ); Project( re_lr, p_dis_lr );
		   p_dis_ml.Reinit( sol[nx][ny] ); Project( re_ml, p_dis_ml );
		   p_dis_mr.Reinit( sol[nx][ny] ); Project( re_mr, p_dis_mr );
		 	p_dis_rl.Reinit( sol[nx][ny] ); Project( re_rl, p_dis_rl );

			NumericalFluxVer(flux, p_dis_lr, p_dis_ml);
		  	res[nx][ny].Add(1.0/dy, flux);
			  
			NumericalFluxVer(flux, p_dis_mr, p_dis_rl);
		  	res[nx][ny].Add(-1.0/dy, flux);
		}

		//对守恒变量和展开系数线性重构	
		static void secondOrderDisCon(const sol_t& sol,
		  		const mesh_t& mesh,
				sol_t& res, 
				size_t nx, size_t ny){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
		  	dis_t re_lr, re_ml, re_mr, re_rl, dummy, flux_l2r, flux_r2l, flux;
			dis_t p_dis_lr, p_dis_ml, p_dis_mr, p_dis_rl;
		  	std::vector<double> normal{1.0,0,0};
		  	double dx = mesh.dx(nx, ny);
		  	dis_t tmp(res[nx][ny].Center(), 
		  			res[nx][ny].Scaling(),
		  			res[nx][ny].GetOrder());
		  	if(nx == 0){
				reconHorCon(sol, mesh, nx, ny, re_ml, re_mr);
				reconHorCon(sol, mesh, nx+1, ny, re_rl, dummy);
		  		PROBLEM::leftBoundaryCondition(re_ml, re_lr);
		  	}else if(nx == (Nx-1)){
				reconHorCon(sol, mesh, nx-1, ny, dummy, re_lr);
				reconHorCon(sol, mesh, nx, ny, re_ml, re_mr);
		  		PROBLEM::rightBoundaryCondition(re_mr, re_rl);
		  	}else{
				reconHorCon(sol, mesh, nx-1, ny, dummy, re_lr);
				reconHorCon(sol, mesh, nx, ny, re_ml, re_mr);
				reconHorCon(sol, mesh, nx+1, ny, re_rl, dummy);
		  	}

		 	p_dis_lr.Reinit( sol[nx][ny] ); Project( re_lr, p_dis_lr );
			p_dis_ml.Reinit( sol[nx][ny] ); Project( re_ml, p_dis_ml );
			p_dis_mr.Reinit( sol[nx][ny] ); Project( re_mr, p_dis_mr );
			p_dis_rl.Reinit( sol[nx][ny] ); Project( re_rl, p_dis_rl );
			
			NumericalFluxHor(flux, p_dis_lr, p_dis_ml);
		  	res[nx][ny].Add(1.0/dx, flux);
			  
			NumericalFluxHor(flux, p_dis_mr, p_dis_rl);
		  	res[nx][ny].Add(-1.0/dx, flux);
				
			normal = {0,1,0};
		  	double dy = mesh.dy(nx, ny);
			if(ny == 0){
				reconVerCon(sol, mesh, nx, ny, re_ml, re_mr);
				reconVerCon(sol, mesh, nx, ny+1, re_rl, dummy);
		  		PROBLEM::bottomBoundaryCondition(re_ml, re_lr);
		  	}else if(ny == (Ny-1)){
				reconVerCon(sol, mesh, nx, ny-1, dummy, re_lr);
				reconVerCon(sol, mesh, nx, ny, re_ml, re_mr);
		  		PROBLEM::topBoundaryCondition(re_mr, re_rl);
		  	}else{
				reconVerCon(sol, mesh, nx, ny-1, dummy, re_lr);
				reconVerCon(sol, mesh, nx, ny, re_ml, re_mr);
				reconVerCon(sol, mesh, nx, ny+1, re_rl, dummy);
		  	}

			p_dis_lr.Reinit( sol[nx][ny] ); Project( re_lr, p_dis_lr );
		   p_dis_ml.Reinit( sol[nx][ny] ); Project( re_ml, p_dis_ml );
		   p_dis_mr.Reinit( sol[nx][ny] ); Project( re_mr, p_dis_mr );
		 	p_dis_rl.Reinit( sol[nx][ny] ); Project( re_rl, p_dis_rl );

			NumericalFluxVer(flux, p_dis_lr, p_dis_ml);
		  	res[nx][ny].Add(1.0/dy, flux);
			  
			NumericalFluxVer(flux, p_dis_mr, p_dis_rl);
		  	res[nx][ny].Add(-1.0/dy, flux);
		}

		//在一维的情形，胡老师在做线性重构时并未调用最原始的HLLFlux()函数
		//而是重新写了NumericalFlux()，这两个函数在功能上都是计算边界上的数值通量
		//但是下面的函数更加简洁一些
		static void NumericalFluxHor( dis_t& flux,
										  const dis_t& dis_l, const dis_t& dis_r){
			dis_t v_dis_l, v_dis_r;

			MulVelocity( 0, dis_l, v_dis_l );
			MulVelocity( 0, dis_r, v_dis_r );

			value_t pv_l[dis_t::dim + 2], pv_r[dis_t::dim + 2];
			dis_l.PrimitiveVars( pv_l );
			dis_r.PrimitiveVars( pv_r );
			
			unsigned int M = dis_l.GetOrder();
			value_t C = MaxRootOfHermitePolynomial<value_t>( M+1 );
			value_t lambda_ll = pv_l[1] - C * sqrt( pv_l[dis_t::dim + 1] ),
			  		lambda_lr = pv_r[1] - C * sqrt( pv_r[dis_t::dim + 1] ),
			  		lambda_rl = pv_l[1] + C * sqrt( pv_l[dis_t::dim + 1] ),
			  		lambda_rr = pv_r[1] + C * sqrt( pv_r[dis_t::dim + 1] );
			value_t lambda_l = std::min<value_t>( lambda_ll, lambda_lr ),
			  		lambda_r = std::max<value_t>( lambda_rl, lambda_rr );

			if ( lambda_l > 0 ) {
			  MulVelocity( 0, dis_l, flux );
			} else if ( lambda_r < 0 ) {
			  MulVelocity( 0, dis_r, flux );
			} else {
			  value_t c1 = lambda_r / (lambda_r - lambda_l);
			  value_t c2 = lambda_l / (lambda_r - lambda_l);
			  value_t c3 = lambda_l * lambda_r / (lambda_r - lambda_l);

			  dis_t v_dis_l, v_dis_r;
			  MulVelocity( 0, dis_l, v_dis_l );
			  MulVelocity( 0, dis_r, v_dis_r );
			  flux = dis_l; flux *= -c3;
			  flux.Add( c3, dis_r );
			  flux.Add( c1, v_dis_l );
			  flux.Add( -c2, v_dis_r );
			}
		}

		static void NumericalFluxVer( dis_t& flux,
										  const dis_t& dis_l, const dis_t& dis_r){
			dis_t v_dis_l, v_dis_r;
			MulVelocity( 1, dis_l, v_dis_l );
			MulVelocity( 1, dis_r, v_dis_r );

			value_t pv_l[dis_t::dim + 2], pv_r[dis_t::dim + 2];
			dis_l.PrimitiveVars( pv_l );
			dis_r.PrimitiveVars( pv_r );
			
			unsigned int M = dis_l.GetOrder();
			value_t C = MaxRootOfHermitePolynomial<value_t>( M+1 );
			value_t lambda_ll = pv_l[2] - C * sqrt( pv_l[dis_t::dim + 1] ),
			  		lambda_lr = pv_r[2] - C * sqrt( pv_r[dis_t::dim + 1] ),
			  		lambda_rl = pv_l[2] + C * sqrt( pv_l[dis_t::dim + 1] ),
			  		lambda_rr = pv_r[2] + C * sqrt( pv_r[dis_t::dim + 1] );
			value_t lambda_l = std::min<value_t>( lambda_ll, lambda_lr ),
			  		lambda_r = std::max<value_t>( lambda_rl, lambda_rr );

			if ( lambda_l > 0 ) {
			  MulVelocity( 1, dis_l, flux );
			} else if ( lambda_r < 0 ) {
			  MulVelocity( 1, dis_r, flux );
			} else {
			  value_t c1 = lambda_r / (lambda_r - lambda_l);
			  value_t c2 = lambda_l / (lambda_r - lambda_l);
			  value_t c3 = lambda_l * lambda_r / (lambda_r - lambda_l);

			  dis_t v_dis_l, v_dis_r;
			  MulVelocity( 1, dis_l, v_dis_l );
			  MulVelocity( 1, dis_r, v_dis_r );
			  flux = dis_l; flux *= -c3;
			  flux.Add( c3, dis_r );
			  flux.Add( c1, v_dis_l );
			  flux.Add( -c2, v_dis_r );
			}
		}

		//对包含密度、速度和温度的守恒向量进行重构
		static void reconHorCon(const sol_t& sol,
				const mesh_t& mesh,
				size_t nx, size_t ny,
				dis_t& re_l,
				dis_t& re_r){
			size_t Nx = mesh.n_ele_x();
			Eigen::VectorXd cv_l, cv_m, cv_r, re_cv_l, re_cv_r;
			if(nx == 0){
				re_l = sol[nx][ny]; re_l *= value_t(1.5); re_l.Add(-0.5, sol[nx+1][ny]);
			  	re_r = sol[nx][ny]; re_r *= value_t(0.5); re_r.Add( 0.5, sol[nx+1][ny]);

				//对守恒变量进行重构
				getcv(sol[nx][ny], cv_l);
				getcv(sol[nx+1][ny], cv_r);
				re_cv_l = 1.5 * cv_l - 0.5 * cv_r;
				re_cv_r = 0.5 * cv_l + 0.5 * cv_r;
		  	}else if(nx == (Nx-1)){
				re_l = sol[nx][ny]; re_l *= value_t(0.5); re_l.Add( 0.5, sol[nx-1][ny]);
			  	re_r = sol[nx][ny]; re_r *= value_t(1.5); re_r.Add(-0.5, sol[nx-1][ny]);

				//对守恒变量进行重构
				getcv(sol[nx-1][ny], cv_l);
				getcv(sol[nx][ny], cv_r);
				re_cv_l = 0.5 * cv_r + 0.5 * cv_l;
				re_cv_r = 1.5 * cv_r - 0.5 * cv_l;
		  	}else{
				getcv(sol[nx-1][ny], cv_l);
				getcv(sol[nx][ny], cv_m);
				getcv(sol[nx+1][ny], cv_r);
				re_cv_l = cv_m - 0.25 * (cv_r - cv_l); 
				re_cv_r = cv_m + 0.25 * (cv_r - cv_l);

				re_l.Reinit(sol[nx][ny]), re_r.Reinit(sol[nx][ny]);
				typename dis_t::Iterator it_re_l = re_l.begin(), it_re_r = re_r.begin();
				typename dis_t::ConstIterator it_l = sol[nx-1][ny].begin(), it_r = sol[nx+1][ny].begin();
				for(typename dis_t::ConstIterator it = sol[nx][ny].begin(), it_end = sol[nx][ny].end();
				    it != it_end; ++it, ++it_l, ++it_r, ++it_re_l, ++it_re_r){
					value_t g = (*it_r - *it_l) / 2.0;
					*it_re_l = *it - 0.5 * g;
					*it_re_r = *it + 0.5 * g;
				}	
		  	}

			velocity_t center;
			value_t scale;
			getPrim(re_cv_l, center, scale);
			re_l.Reinit(center, scale);
			getPrim(re_cv_r, center, scale);
			re_r.Reinit(center, scale);
		}

		static void reconVerCon(const sol_t& sol,
				const mesh_t& mesh,
				size_t nx, size_t ny,
				dis_t& re_l,
				dis_t& re_r){
			size_t Ny = mesh.n_ele_y();
			Eigen::VectorXd cv_l, cv_m, cv_r, re_cv_l, re_cv_r;
			if(ny == 0){
				re_l = sol[nx][ny]; re_l *= value_t(1.5); re_l.Add(-0.5, sol[nx][ny+1]);
			  	re_r = sol[nx][ny]; re_r *= value_t(0.5); re_r.Add( 0.5, sol[nx][ny+1]);

				//对守恒变量进行重构
				getcv(sol[nx][ny], cv_l);
				getcv(sol[nx][ny+1], cv_r);
				re_cv_l = 1.5 * cv_l - 0.5 * cv_r;
				re_cv_r = 0.5 * cv_l + 0.5 * cv_r;

		  	}else if(ny == (Ny-1)){
				re_l = sol[nx][ny]; re_l *= value_t(0.5); re_l.Add( 0.5, sol[nx][ny-1]);
			  	re_r = sol[nx][ny]; re_r *= value_t(1.5); re_r.Add(-0.5, sol[nx][ny-1]);

				//对守恒变量进行重构
				getcv(sol[nx][ny-1], cv_l);
				getcv(sol[nx][ny], cv_r);
				re_cv_l = 0.5 * cv_r + 0.5 * cv_l;
				re_cv_r = 1.5 * cv_r - 0.5 * cv_l;
				
		  	}else{
				getcv(sol[nx][ny-1], cv_l);
				getcv(sol[nx][ny], cv_m);
				getcv(sol[nx][ny+1], cv_r);
				re_cv_l = cv_m - 0.25 * (cv_r - cv_l); 
				re_cv_r = cv_m + 0.25 * (cv_r - cv_l);

				re_l.Reinit(sol[nx][ny]), re_r.Reinit(sol[nx][ny]);
				typename dis_t::Iterator it_re_l = re_l.begin(), it_re_r = re_r.begin();
				typename dis_t::ConstIterator it_l = sol[nx][ny-1].begin(), it_r = sol[nx][ny+1].begin();
				for(typename dis_t::ConstIterator it = sol[nx][ny].begin(), it_end = sol[nx][ny].end();
				    it != it_end; ++it, ++it_l, ++it_r, ++it_re_l, ++it_re_r){
					value_t g = (*it_r - *it_l) / 2;
					*it_re_l = *it - 0.5 * g;
					*it_re_r = *it + 0.5 * g;
				}	
		  	}

			velocity_t center;
			value_t scale;
			getPrim(re_cv_l, center, scale);
			re_l.Reinit(center, scale);
			getPrim(re_cv_r, center, scale);
			re_r.Reinit(center, scale);
		}

		static void getcv(const dis_t& dis,
				Eigen::VectorXd& cv){
			double pv[5];
			dis.PrimitiveVars(pv);
			double rho = pv[0];
			double u1 = pv[1];
			double u2 = pv[2];
			double u3 = pv[3];
			double theta = pv[4];
			cv.resize(5);
			cv[0] = rho;
			cv[1] = rho * u1;
			cv[2] = rho * u2;
			cv[3] = rho * u3;
			cv[4] = rho * (1.5 * theta + (u1*u1 + u2*u2 + u3*u3)/2);
		}

		static void getPrim(Eigen::VectorXd& cv,
		  	  velocity_t& center,
		  	  value_t& scale){
		    double rho = cv[0];
		    double u1 = cv[1] / rho;
		    double u2 = cv[2] / rho;
			double u3 = cv[3] / rho;
		  	double theta = (cv[4] / rho - (u1*u1 + u2*u2 + u3*u3)/2) / 1.5;
		  	center = {u1, u2, u3};
		  	scale = sqrt(theta);
		}
		
		//对展开式参数和矩系数作重构
		static void reconHorPri(const sol_t& sol,
				const mesh_t& mesh,
				size_t nx, size_t ny,
				dis_t& re_l,
				dis_t& re_r){
			size_t Nx = mesh.n_ele_x();
			if(nx == 0){
				re_l = sol[nx][ny]; re_l *= value_t(1.5); re_l.Add(-0.5, sol[nx+1][ny]);
			  	re_r = sol[nx][ny]; re_r *= value_t(0.5); re_r.Add( 0.5, sol[nx+1][ny]);

				re_l.Reinit( sol[nx][ny].Center() * 1.5 - sol[nx+1][ny].Center() * 0.5,
                   sol[nx][ny].Scaling() * 1.5 - sol[nx+1][ny].Scaling() * 0.5 );
      		re_r.Reinit( sol[nx][ny].Center() * 0.5 + sol[nx+1][ny].Center() * 0.5,
                   sol[nx][ny].Scaling() * 0.5 + sol[nx+1][ny].Scaling() * 0.5 );

		  	}else if(nx == (Nx-1)){
				re_l = sol[nx][ny]; re_l *= value_t(0.5); re_l.Add( 0.5, sol[nx-1][ny]);
			  	re_r = sol[nx][ny]; re_r *= value_t(1.5); re_r.Add(-0.5, sol[nx-1][ny]);

				re_l.Reinit( sol[nx][ny].Center() * 0.5 + sol[nx-1][ny].Center() * 0.5,
                   sol[nx][ny].Scaling() * 0.5 + sol[nx-1][ny].Scaling() * 0.5 );
      		re_r.Reinit( sol[nx][ny].Center() * 1.5 - sol[nx-1][ny].Center() * 0.5,
                   sol[nx][ny].Scaling() * 1.5 - sol[nx-1][ny].Scaling() * 0.5 );
				
		  	}else{
				velocity_t center_diff = value_t( 0.25 ) * (sol[nx+1][ny].Center() - sol[nx-1][ny].Center());
				velocity_t center_l( sol[nx][ny].Center() - center_diff );
				velocity_t center_r( sol[nx][ny].Center() + center_diff );

				value_t scale_diff = 0.25 * ( sol[nx+1][ny].Scaling() - sol[nx-1][ny].Scaling() );
				re_l.Reinit( center_l, sol[nx][ny].Scaling() - scale_diff, sol[nx][ny].GetOrder() );
				re_r.Reinit( center_r, sol[nx][ny].Scaling() + scale_diff, sol[nx][ny].GetOrder() );

				typename dis_t::Iterator it_re_l = re_l.begin(), it_re_r = re_r.begin();
				typename dis_t::ConstIterator it_l = sol[nx-1][ny].begin(), it_r = sol[nx+1][ny].begin();
				for(typename dis_t::ConstIterator it = sol[nx][ny].begin(), it_end = sol[nx][ny].end();
				    it != it_end; ++it, ++it_l, ++it_r, ++it_re_l, ++it_re_r){
					value_t g = (*it_r - *it_l) / 2.0;
					*it_re_l = *it - 0.5 * g;
					*it_re_r = *it + 0.5 * g;
				}	
		  	}
		}

		static void reconVerPri(const sol_t& sol,
				const mesh_t& mesh,
				size_t nx, size_t ny,
				dis_t& re_l,
				dis_t& re_r){
			size_t Ny = mesh.n_ele_y();
			if(ny == 0){
				re_l = sol[nx][ny]; re_l *= value_t(1.5); re_l.Add(-0.5, sol[nx][ny+1]);
			  	re_r = sol[nx][ny]; re_r *= value_t(0.5); re_r.Add( 0.5, sol[nx][ny+1]);

				re_l.Reinit( sol[nx][ny].Center() * 1.5 - sol[nx][ny+1].Center() * 0.5,
                   sol[nx][ny].Scaling() * 1.5 - sol[nx][ny+1].Scaling() * 0.5 );
      		re_r.Reinit( sol[nx][ny].Center() * 0.5 + sol[nx][ny+1].Center() * 0.5,
                   sol[nx][ny].Scaling() * 0.5 + sol[nx][ny+1].Scaling() * 0.5 );

		  	}else if(ny == (Ny-1)){
				re_l = sol[nx][ny]; re_l *= value_t(0.5); re_l.Add( 0.5, sol[nx][ny-1]);
			  	re_r = sol[nx][ny]; re_r *= value_t(1.5); re_r.Add(-0.5, sol[nx][ny-1]);

				re_l.Reinit( sol[nx][ny].Center() * 0.5 + sol[nx][ny-1].Center() * 0.5,
                   sol[nx][ny].Scaling() * 0.5 + sol[nx][ny-1].Scaling() * 0.5 );
      		re_r.Reinit( sol[nx][ny].Center() * 1.5 - sol[nx][ny-1].Center() * 0.5,
                   sol[nx][ny].Scaling() * 1.5 - sol[nx][ny-1].Scaling() * 0.5 );
				
		  	}else{
				velocity_t center_diff = value_t( 0.25 ) * (sol[nx][ny+1].Center() - sol[nx][ny-1].Center());
				velocity_t center_l( sol[nx][ny].Center() - center_diff );
				velocity_t center_r( sol[nx][ny].Center() + center_diff );

				value_t scale_diff = 0.25 * ( sol[nx][ny+1].Scaling() - sol[nx][ny-1].Scaling() );
				re_l.Reinit( center_l, sol[nx][ny].Scaling() - scale_diff, sol[nx][ny].GetOrder() );
				re_r.Reinit( center_r, sol[nx][ny].Scaling() + scale_diff, sol[nx][ny].GetOrder() );

				typename dis_t::Iterator it_re_l = re_l.begin(), it_re_r = re_r.begin();
				typename dis_t::ConstIterator it_l = sol[nx][ny-1].begin(), it_r = sol[nx][ny+1].begin();
				for(typename dis_t::ConstIterator it = sol[nx][ny].begin(), it_end = sol[nx][ny].end();
				    it != it_end; ++it, ++it_l, ++it_r, ++it_re_l, ++it_re_r){
					value_t g = (*it_r - *it_l) / 2;
					*it_re_l = *it - 0.5 * g;
					*it_re_r = *it + 0.5 * g;
				}	
		  	}
		}
};

# endif
