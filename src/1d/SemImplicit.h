/**
 * @file SemImplicit.h
 * @brief 碰撞项作隐式处理的格式
 * @author Guanghan Li NUAA
 * @version 1
 * @date 2022-11-22
 */

# ifndef __SEMIMPLICIT__H__
# define __SEMIMPLICIT__H__

# include <omp.h>

template <typename _RESIDUAL>
class SemImpShakhov{
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
	private:
	    double CFL;
	    double dt;
	    double tolerance;
	public:
		  
		SemImpShakhov() : CFL(0.5), tolerance(1e-8) {}

		SemImpShakhov(const ptr_mesh_t& p_mesh, const ConfigFile& cf){
		  	read_config(p_mesh, cf);
		}

		SemImpShakhov(const ptr_mesh_t& p_mesh, std::istream& is){
		  	is.read((char *)&CFL, sizeof(double));
		  	is.read((char *)&tolerance, sizeof(double));
		}

		void read_config( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) {
		  	CFL = (double)cf.Value( "Time", "CFL" );
		  	tolerance = (double)cf.Value( "Error", "Tol" );
		}

		void initial_value( sol_t& sol, const ptr_mesh_t& p_mesh, 
                    const ptr_force_t& p_force_handle ) {
		  	PROBLEM::initial_value( sol, *p_mesh );
		}

		void PostProcessing( sol_t& sol, const ptr_mesh_t& p_mesh ) {
			PROBLEM::PostProcessing( sol, *p_mesh );
		};

		void reset_tol( const double _tol ) { 
		  	tolerance = _tol;
		}

		void dump_data( std::ostream& os ) const {
		  	os.write((const char *)&CFL, sizeof(double));
		  	os.write((const char *)&tolerance, sizeof(double));
		}

		void print_config() const {
		  	std::cout << "SemImpShakhov: CFL = " << CFL << std::endl;
		  	std::cout << "SemImpShakhov: Tol = " << tolerance << std::endl;
		}

		void Reinit(const ptr_mesh_t& p_mesh) {};
		
		void ExactSolve(sol_t& sol, const sol_t& rhs, 
				const ptr_mesh_t& p_mesh, 
				const ptr_force_t& p_force_handle, 
		  				  const unsigned int max_steps=0 ) {
		  	double initial_res;
		  	RESIDUAL::GetResidualPL2Norm( initial_res, sol, rhs, *p_mesh, *p_force_handle );
		  	// 以下的设置是不是效率又防止最后残量降不下来.
		  	const double _Tol = (initial_res > 10*tolerance) ? tolerance : 0.1 * tolerance;

		  	double res;
		  	unsigned int step = 0;
		  	do{
		  		++step;
		  		Solve( sol, rhs, p_mesh, p_force_handle, 1 );
		   		//RESIDUAL::GetResidualMaxNorm( res, sol, rhs, *p_mesh, *p_force_handle );
		  		RESIDUAL::GetResidualPL2Norm( res, sol, rhs, *p_mesh, *p_force_handle );
			} while (res >= _Tol && (step<max_steps || max_steps==0) );
		   		//std::cout<< "coarsest steps: = " << step << ", res = " << res << std::endl;
		}
		 
		void Solve( sol_t& sol, const sol_t& rhs, 
				const ptr_mesh_t& p_mesh, 
              	const ptr_force_t& p_force_handle,
              	const unsigned int steps = 1 ) {
		 	sol_t res(sol.size());
		 	for(unsigned int i = 0; i < steps; i++){
		 		TimeStep(sol, *p_mesh);
		 		GetConvection(sol, rhs, res, *p_mesh);
		 		Update(sol, res, p_force_handle);
		 	}
		}    

		void TimeStep(const sol_t& sol, const mesh_t& mesh){
		     int n_ele = sol.size();
		     dt = 10;
#pragma omp parallel for reduction(min:dt)
		     for(int i = 0; i < n_ele; i++){
		   		dt = std::min(dt, TimeStep(sol, mesh, i));
		     }
		}

		double TimeStep(const sol_t& sol, const mesh_t& mesh, int i){
		     value_t C_M =
		   		Transform<value_t>::MaxRootOfHermitePolynomial(sol[i].GetOrder()
		   				  + 1);
		     double lambda = fabs(sol[i].Center()[0]) + C_M * sol[i].Scaling();
		     return CFL * mesh.Length(i) / lambda;
		}

		void GetConvection(const sol_t& sol,
				const sol_t& rhs,
				sol_t& res,
				const mesh_t& mesh) {
		     int n_ele = sol.size();
#pragma omp parallel for
		     for(int i = 0; i < n_ele; i++){
		   		res[i].Reinit(sol[i]);
					Project( rhs[i], res[i] );			
		   		RESIDUAL::Convection(res[i], sol[i], sol, mesh, i);
		     }
		}

		void Update(sol_t& sol,
				const sol_t& res,
				const ptr_force_t& p_force_handle) {
		    int n_ele = sol.size();
#pragma omp parallel for
		    for(int i = 0; i < n_ele; i++){
		   	  	Update(sol[i], res[i], p_force_handle, i);
		    }
		}

		void Update(dis_t& dis_sol,
				const dis_t& dis_res,
				const ptr_force_t& p_force_handle,
				unsigned int i){
		  	dis_t new_sol;
		  	new_sol.Reinit(dis_sol);
			//1.考虑对流项
		  	dis_sol.Add(dt, dis_res);
			//2.标准投影，获得新时刻的物理量
			new_sol.ProjectToStdSpace(dis_sol);
		  	std::swap(dis_sol, new_sol);
			double nu = PROBLEM::OneOverTau(dis_sol, *p_force_handle, i);
			const int dim = dis_t::dim;
			double heat_flux[dim];
		  	dis_sol.HeatFlux( heat_flux );
			//3.作系数调整，得到新时刻分布函数的展开式，首先处理BGK的部分
		  	typename dis_t::Iterator the_dis = dis_sol.begin(2);
		  	typename dis_t::Iterator end_dis = dis_sol.end();
		  	for(; the_dis != end_dis; the_dis++){
				*the_dis *= 1.0 / ( 1.0 + dt * nu);
		  	}
			//再处理alpha=3的部分系数
			double factor = dt * nu * ( 1 - PROBLEM::Pr ) / (dim+2);
			typename dis_t::Indices ind; // 初始化 ind 应该是 0
		  	for(int j = 0; j < dim; j++){
				ind[j] = 1;
				//新时刻的热通量预测值
				heat_flux[j] /= (1 + dt * nu * PROBLEM::Pr);
		  		for(int i = 0; i < dim; i++){
					ind[i] += 2;
		  		    dis_sol(ind) +=  1 / (1.0 + dt * nu) * factor * heat_flux[j];
		  		    ind[i] -= 2;
		  		}
		  		ind[j] = 0;
		  	}
		}
};

template <typename _RESIDUAL>
class SemImpShakhovSGS{
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
	private:
	    double CFL;
	    double dt;
	    double tolerance;
	public:
		  
		SemImpShakhovSGS() : CFL(0.5), tolerance(1e-8) {}

		SemImpShakhovSGS(const ptr_mesh_t& p_mesh, const ConfigFile& cf){
		  	read_config(p_mesh, cf);
		}

		SemImpShakhovSGS(const ptr_mesh_t& p_mesh, std::istream& is){
		  	is.read((char *)&CFL, sizeof(double));
		  	is.read((char *)&tolerance, sizeof(double));
		}

		void read_config( const ptr_mesh_t& p_mesh, const ConfigFile& cf ) {
		  	CFL = (double)cf.Value( "Time", "CFL" );
		  	tolerance = (double)cf.Value( "Error", "Tol" );
		}

		void initial_value( sol_t& sol, const ptr_mesh_t& p_mesh, 
                    const ptr_force_t& p_force_handle ) {
		  	PROBLEM::initial_value( sol, *p_mesh );
		}

		void PostProcessing( sol_t& sol, const ptr_mesh_t& p_mesh ) {
			//这里需要作守恒性校正
			PROBLEM::PostProcessing( sol, *p_mesh );
		};

		void reset_tol( const double _tol ) { 
		  	tolerance = _tol;
		}

		void dump_data( std::ostream& os ) const {
		  	os.write((const char *)&CFL, sizeof(double));
		  	os.write((const char *)&tolerance, sizeof(double));
		}

		void print_config() const {
		  	std::cout << "SemImpShakhovSGS: CFL = " << CFL << std::endl;
		  	std::cout << "SemImpShakhovSGS: Tol = " << tolerance << std::endl;
		}

		void Reinit(const ptr_mesh_t& p_mesh) {};
		
		void ExactSolve(sol_t& sol, const sol_t& rhs, 
				const ptr_mesh_t& p_mesh, 
				const ptr_force_t& p_force_handle, 
		  				  const unsigned int max_steps=0 ) {
		  	double initial_res;
		  	RESIDUAL::GetResidualPL2Norm( initial_res, sol, rhs, *p_mesh, *p_force_handle );
		  	// 以下的设置是不是效率又防止最后残量降不下来.
		  	const double _Tol = (initial_res > 10*tolerance) ? tolerance : 0.1 * tolerance;

		  	double res;
		  	unsigned int step = 0;
		  	do{
		  		++step;
		  		Solve( sol, rhs, p_mesh, p_force_handle, 1 );
		   		//RESIDUAL::GetResidualMaxNorm( res, sol, rhs, *p_mesh, *p_force_handle );
		  		RESIDUAL::GetResidualPL2Norm( res, sol, rhs, *p_mesh, *p_force_handle );
			} while (res >= _Tol && (step<max_steps || max_steps==0) );
		   		//std::cout<< "coarsest steps: = " << step << ", res = " << res << std::endl;
		}
		 
		void Solve( sol_t& sol, const sol_t& rhs, 
				const ptr_mesh_t& p_mesh, 
            	const ptr_force_t& p_force_handle,
            	const unsigned int steps = 1 ) {
		 	for(unsigned int i = 0; i < steps; i++){
		 		LeftToRight(sol, rhs, p_mesh, p_force_handle);
		  		RightToLeft(sol, rhs, p_mesh, p_force_handle);
		 	}
		}

		void LeftToRight(sol_t& sol, const sol_t& rhs, 
            const ptr_mesh_t& p_mesh, 
            const ptr_force_t& p_force_handle){
		   	int n_ele = sol.size();
			sol_t res(n_ele);
		  	for(int nx = 0; nx < n_ele; nx++){
		  		LocalCellOneStep(sol, res, rhs, p_mesh, p_force_handle, nx);
		  	}
		}

		void RightToLeft(sol_t& sol, const sol_t& rhs, 
            const ptr_mesh_t& p_mesh, 
            const ptr_force_t& p_force_handle){
		    int n_ele = sol.size();
			sol_t res(n_ele);
		  	for(int nx = n_ele; nx > 0; nx--){
		  		LocalCellOneStep(sol, res, rhs, p_mesh, p_force_handle, nx-1);
		  	}
		}

		void LocalCellOneStep(sol_t& sol, sol_t& res,
				const sol_t& rhs, 
				const ptr_mesh_t& p_mesh, 
            	const ptr_force_t& p_force_handle,
		  	  	int nx){
		    dt = TimeStep(sol, *p_mesh, nx);
			GetConvection(sol, res, rhs, *p_mesh, nx);
		  	//dis_t res;
		  	//res.Reinit(sol[nx]);
		  	//Project( rhs[nx], res);					  
		  	//RESIDUAL::Convection(res, sol[nx], sol, *p_mesh, nx);
		  	update(sol[nx], res[nx], p_force_handle, nx);
		}

		void TimeStep(const sol_t& sol, const mesh_t& mesh){
		     int n_ele = sol.size();
		     dt = 10;
		     for(int i = 0; i < n_ele; i++){
		   		dt = std::min(dt, TimeStep(sol, mesh, i));
		     }
		}

		double TimeStep(const sol_t& sol, const mesh_t& mesh, int i){
		     value_t C_M =
		   		Transform<value_t>::MaxRootOfHermitePolynomial(sol[i].GetOrder()
		   				  + 1);
		     double lambda = fabs(sol[i].Center()[0]) + C_M * sol[i].Scaling();
		     return CFL * mesh.Length(i) / lambda;
		}

		void GetConvection(const sol_t& sol, sol_t& res,
				const sol_t& rhs,
				const mesh_t& mesh,
				int nx) {
			res[nx].Reinit(sol[nx]);
			Project( rhs[nx], res[nx] );			
		   RESIDUAL::Convection(res[nx], sol[nx], sol, mesh, nx);
		}

		void update(dis_t& dis_sol,
				const dis_t& dis_res,
				const ptr_force_t& p_force_handle,
				unsigned int i){
		  	dis_t new_sol;
		  	new_sol.Reinit(dis_sol);
			//1.考虑对流项
		  	dis_sol.Add(dt, dis_res);
			//2.标准投影，获得新时刻的物理量
			new_sol.ProjectToStdSpace(dis_sol);
		  	std::swap(dis_sol, new_sol);
			double nu = PROBLEM::OneOverTau(dis_sol, *p_force_handle, i);
			const int dim = dis_t::dim;
			double heat_flux[dim];
		  	dis_sol.HeatFlux( heat_flux );
			//3.作系数调整，得到新时刻分布函数的展开式，首先处理BGK的部分
		  	typename dis_t::Iterator the_dis = dis_sol.begin(2);
		  	typename dis_t::Iterator end_dis = dis_sol.end();
		  	for(; the_dis != end_dis; the_dis++){
				*the_dis *= 1.0 / ( 1.0 + dt * nu);
		  	}

			//再处理alpha=3的部分系数
			double factor = dt * nu * ( 1 - PROBLEM::Pr ) / (dim+2);
			typename dis_t::Indices ind; // 初始化 ind 应该是 0
		  	for(int j = 0; j < dim; j++){
				ind[j] = 1;
				//新时刻的热通量预测值
				heat_flux[j] /= (1 + dt * nu * PROBLEM::Pr);
		  		for(int i = 0; i < dim; i++){
					ind[i] += 2;
		  		    dis_sol(ind) +=  1 / (1.0 + dt * nu) * factor * heat_flux[j];
		  		    ind[i] -= 2;
		  		}
		  		ind[j] = 0;
		  	}
		}
};

# endif
