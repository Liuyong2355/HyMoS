/**
 * @file SingleGridMethod.h
 * @brief 单层网格求解器类
 * @author Guanghan Li NUAA
 * @version 1.0
 * @date 2022-10-12
 */

# ifndef __SINGLEGRIDMETHOD__H__
# define __SINGLEGRIDMETHOD__H__

# include "EulerSolver.h"
# include "config.h"
# include <limits>

template <typename RESIDUAL>
class SingleGridMethod{
	public:
	    typedef RESIDUAL Residual;
	    typedef typename Residual::Problem Problem;
	    typedef typename Residual::sol_t sol_t;
	    typedef typename sol_t::dis_t dis_t;
	    typedef typename sol_t::mesh_t mesh_t; 
	    typedef typename dis_t::value_type value_t;
	    typedef typename dis_t::Indices index_t;
	    typedef typename dis_t::Velocity velocity_t;
	    typedef std::shared_ptr<mesh_t> ptr_mesh_t;
	public:
		double dt;
		double CFL;
		unsigned int EulerEqsSolvingSteps;
    public:
		SingleGridMethod() : dt(1.0), CFL(0.5), EulerEqsSolvingSteps(1) {};

		SingleGridMethod(const ptr_mesh_t& p_mesh, std::istream& is){
			is.read((char *)&CFL, sizeof(double));
		}

		void initialize(sol_t& sol, mesh_t& mesh){
		  	Problem::initialize(sol, mesh);
		}

		void read_config(const ConfigFile& cf){
		  	dt = (double)cf.Value("Time", "delta_t");
		  	CFL = (double)cf.Value("Time", "CFL");
		  	EulerEqsSolvingSteps = (int)cf.Value("SingleGridMethod",
		  			  "EulerEqsSolvingSteps");
		}

		void dump_data( std::ostream& os ) const {
		  	os.write((const char *)&CFL, sizeof(double));
		}

		void load_data (std::istream& is){
			is.read((char *)&CFL, sizeof(double));			
		}

		void timeStep(const sol_t& sol, const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
		  	double dt_new = std::numeric_limits<double>::max();
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		    double dt_local = timeStep(nx, ny, sol, mesh);
		  		    dt_new = (dt_local < dt_new) ? dt_local : dt_new;
		  		}
		  	}
		  	dt = dt_new;
		}

		double timeStep(size_t nx, size_t ny,
				const sol_t& sol, const mesh_t& mesh){
		  	double dx = mesh.dx(nx, ny);
		  	double dy = mesh.dy(nx, ny);
		  	const dis_t& dis = sol[nx][ny];
		  	typename dis_t::Velocity v = dis.Center();
		  	double sqrt_theta = dis.Scaling();
		  	unsigned int M = dis.GetOrder();
		  	double C = MaxRootOfHermitePolynomial<double>(M + 1);
		  	return CFL / ((fabs(v[0]) + C * sqrt_theta) / dx +
		  			        (fabs(v[1]) + C * sqrt_theta) / dy);
		}
		
		void update(sol_t& sol, sol_t& res, const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
		  	#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		    update(nx, ny, sol, res);
		  		}
		  	}
		}

		void update(size_t nx, size_t ny, 
				sol_t& sol, sol_t& res){
		  	dis_t& dis_sol = sol[nx][ny];
		  	const dis_t& dis_res = res[nx][ny];
		  	dis_t dis_tmp(dis_sol); 
		  	Project(dis_res, dis_tmp);
		  	dis_sol.Add(dt, dis_tmp);
		  	dis_t tmp_dis = dis_sol;
		  	tmp_dis.ProjectToStdSpace(dis_sol);
		  	std::swap(tmp_dis,dis_sol);			
		}

		void postProcessing( sol_t& sol, const mesh_t& mesh ) {
		  	Problem::postProcessing( sol, mesh );
		}
        
		virtual void solve(sol_t& sol, sol_t& rhs, const mesh_t& mesh){}
		
		virtual void print_config(){}

};

template <typename RESIDUAL>
class ExplicitTime : public SingleGridMethod<RESIDUAL>{
	public:
	    typedef RESIDUAL Residual;
	    typedef SingleGridMethod<Residual> BASE;
	    typedef typename Residual::sol_t sol_t;
	    typedef typename sol_t::dis_t dis_t;
	    typedef typename sol_t::mesh_t mesh_t; 
	    typedef std::shared_ptr<mesh_t> ptr_mesh_t;
	private:
	    sol_t res;
	public:
		ExplicitTime(){};

		ExplicitTime(const ptr_mesh_t& p_mesh, std::istream& is){
		  	is.read((char *)&(BASE::CFL), sizeof(double));
		}
		  	
		void solve(sol_t& sol, const sol_t& rhs,
				const mesh_t& mesh,
				const size_t steps = 1){
		  	for(size_t i = 0; i < steps; i++){
		  		BASE::timeStep(sol, mesh);
		  		Residual::getResidual(sol, rhs, res, mesh);
		  		BASE::update(sol, res, mesh);
		  	}
		}
		
		void print_config(){
		  	std::cout << "ExplicitTime: CFL = " << BASE::CFL << std::endl;
		}
		
		void dump_data( std::ostream& os ) const {
		  	os.write((const char *)&(BASE::CFL), sizeof(double));
		}
};

template <typename RESIDUAL>
class GaussSeidel : public SingleGridMethod<RESIDUAL>{
	public:
	    typedef RESIDUAL Residual;
	    typedef SingleGridMethod<Residual> BASE;
	    typedef typename Residual::sol_t sol_t;
	    typedef typename sol_t::dis_t dis_t;
	    typedef typename sol_t::mesh_t mesh_t; 
	    typedef typename dis_t::value_type value_t;
	    typedef typename dis_t::Indices index_t;
	    typedef typename dis_t::Velocity velocity_t;
	    typedef std::shared_ptr<mesh_t> ptr_mesh_t;
	private:
	    sol_t res;
	public:
		GaussSeidel(){};

		GaussSeidel(const ptr_mesh_t& p_mesh, std::istream& is){
		  	is.read((char *)&(BASE::CFL), sizeof(double));
		}

		void solve(sol_t& sol, const sol_t& rhs, const mesh_t& mesh,
		  		    size_t steps = 1){
		  	for(size_t i = 0; i < steps; i++){
		  		symSweeping(sol, rhs, mesh);	
		  	}
		}

		void symSweeping(sol_t& sol, const sol_t& rhs, const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
			res.Reinit(sol);
			#pragma omp parallel for
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
				   BASE::dt = BASE::timeStep(nx, ny, sol, mesh);
				   Residual::getResidual(sol, rhs, mesh, res, nx, ny);
				   BASE::update(nx, ny, sol, res);
				}
			}
		}

		void print_config(){
		  	std::cout << "GaussSeidel: CFL = " << BASE::CFL << std::endl;
		}
		
		void dump_data( std::ostream& os ) const {
		  	os.write((const char *)&(BASE::CFL), sizeof(double));
		}
};

template <typename RESIDUAL>
class SymGaussSeidel : public SingleGridMethod<RESIDUAL>{
	public:
	    typedef RESIDUAL Residual;
	    typedef SingleGridMethod<Residual> BASE;
	    typedef typename Residual::sol_t sol_t;
	    typedef typename sol_t::dis_t dis_t;
	    typedef typename sol_t::mesh_t mesh_t; 
	    typedef typename dis_t::value_type value_t;
	    typedef typename dis_t::Indices index_t;
	    typedef typename dis_t::Velocity velocity_t;
	    typedef std::shared_ptr<mesh_t> ptr_mesh_t;
	private:
	    sol_t res;
	public:
		SymGaussSeidel(){};

		SymGaussSeidel(const ptr_mesh_t& p_mesh, std::istream& is){
		  	is.read((char *)&(BASE::CFL), sizeof(double));
		}

		void solve(sol_t& sol, const sol_t& rhs, const mesh_t& mesh,
		  		    size_t steps = 1){
		  	for(size_t i = 0; i < steps; i++){
		  		symSweeping(sol, rhs, mesh);	
		  	}
		}

		void symSweeping(sol_t& sol, const sol_t& rhs, const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();

			res.Reinit(sol);
			#pragma omp parallel for
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
				   BASE::dt = BASE::timeStep(nx, ny, sol, mesh);
				   Residual::getResidual(sol, rhs, mesh, res, nx, ny);
				   BASE::update(nx, ny, sol, res);
				}
			}
			
			res.Reinit(sol);
			#pragma omp parallel for
			for(size_t ny = Ny; ny > 0; ny--){
				for(size_t nx = Nx; nx > 0; nx--){
				   BASE::dt = BASE::timeStep(nx-1, ny-1, sol, mesh);
				   Residual::getResidual(sol, rhs, mesh, res, nx-1, ny-1);
				   BASE::update(nx-1, ny-1, sol, res);
				}
			}
		}

		void print_config(){
		  	std::cout << "SymGaussSeidel: CFL = " << BASE::CFL << std::endl;
		}
		
		void dump_data( std::ostream& os ) const {
		  	os.write((const char *)&(BASE::CFL), sizeof(double));
		}
};

template <typename RESIDUAL>
class FastSweeping : public SingleGridMethod<RESIDUAL>{
	public:
	    typedef RESIDUAL Residual;
	    typedef SingleGridMethod<Residual> BASE;
	    typedef typename Residual::sol_t sol_t;
	    typedef typename sol_t::dis_t dis_t;
	    typedef typename sol_t::mesh_t mesh_t; 
	    typedef std::shared_ptr<mesh_t> ptr_mesh_t;
	private:
	    sol_t res;
	public:
		FastSweeping(){};

		FastSweeping(const ptr_mesh_t& p_mesh, std::istream& is){
		  	is.read((char *)&(BASE::CFL), sizeof(double));
		}
		
		void solve(sol_t& sol, const sol_t& rhs, const mesh_t& mesh, 
		  		    const size_t steps = 1){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	res.Reinit(sol);
		  	for(size_t i = 0; i < steps; i++){
		  		//D1
		  		#pragma omp parallel for
		  		for(size_t ny = 0; ny < Ny; ny++){
		  		    for(size_t nx = 0; nx < Nx; nx++){
		  		   	   	BASE::dt = BASE::timeStep(nx, ny, sol, mesh);
		  		   	   	Residual::getResidual(sol, rhs, mesh, res, nx, ny);
		  		   	   	BASE::update(nx, ny, sol, res);
		  		    }
		  		}

		  		//D2
		  		#pragma omp parallel for
		  		for(size_t ny = 0; ny < Ny; ny++){
		  		    for(size_t nx = Nx; nx > 0; nx--){
		  		   	   	BASE::dt = BASE::timeStep(nx-1, ny, sol,mesh);
		  		   	   	Residual::getResidual(sol, rhs, mesh, res, nx-1, ny);
		  		   	   	BASE::update(nx-1, ny, sol, res);
		  		    }
		  		}

		  		//D3
		  		#pragma omp parallel for
		  		for(size_t ny = Ny; ny > 0; ny--){
		  		    for(size_t nx = Nx; nx > 0; nx--){
		  		   	   BASE::dt = BASE::timeStep(nx-1, ny-1, sol,mesh);
		  		   	   Residual::getResidual(sol, rhs, mesh, res, nx-1, ny-1);
		  		   	   BASE::update(nx-1, ny-1, sol, res);
		  		    }
		  		}

		  		//D4
		  		#pragma omp parallel for
		  		for(size_t ny = Ny; ny > 0; ny--){
		  		    for(size_t nx = 0; nx < Nx; nx++){
		  		   	   	BASE::dt = BASE::timeStep(nx, ny-1, sol,mesh);
		  		   	   	Residual::getResidual(sol, rhs, mesh, res, nx, ny-1);
		  		   	   	BASE::update(nx, ny-1, sol, res);
		  		    }
		  		}
		  	}
		}

		void print_config(){
		  	std::cout << "FastSweeping: CFL = " << BASE::CFL << std::endl;
		}

		void dump_data( std::ostream& os ) const {
		  	os.write((const char *)&(BASE::CFL), sizeof(double));
		}
};

template <typename RESIDUAL>
class SISShakhovTime : public SingleGridMethod<RESIDUAL>{
	public:
	    typedef RESIDUAL Residual;
	    typedef SingleGridMethod<Residual> BASE;
	    typedef typename BASE::Problem Problem;
	    typedef typename Residual::sol_t sol_t;
	    typedef typename sol_t::dis_t dis_t;
		typedef typename dis_t::Indices index_t;
	    typedef typename sol_t::mesh_t mesh_t; 
	    typedef std::shared_ptr<mesh_t> ptr_mesh_t;
	private:
	    sol_t res;
	public:
		SISShakhovTime(){};

		SISShakhovTime(const ptr_mesh_t& p_mesh, std::istream& is){
		  	is.read((char *)&(BASE::CFL), sizeof(double));
		}
		
		void solve(sol_t& sol, const sol_t& rhs, const mesh_t& mesh, 
		 		    const size_t steps = 1){
		 	res.Reinit(sol);
		 	for(size_t i = 0; i < steps; i++){
		 		BASE::timeStep(sol, mesh);
		 		ProjectRHS(rhs, res, mesh);
		 		Residual::convection(sol, res, mesh); 
		 		update(sol, res, mesh);
		 	}
		}

		void update(sol_t& sol, const sol_t& res, const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
			#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		    dis_t& dis_sol = sol[nx][ny];
		  		    const dis_t& dis_res = res[nx][ny];
		  		    update(dis_sol, dis_res, BASE::dt);
		  		}
		  	}
		}

		void update(dis_t& dis_sol,
				const dis_t& dis_res,
				double dt){
		  	dis_t new_sol;
		  	new_sol.Reinit(dis_sol);
			//1.考虑对流项
		  	dis_sol.Add(dt, dis_res);
			//2.标准投影，获得新时刻的物理量
			new_sol.ProjectToStdSpace(dis_sol);
		  	std::swap(dis_sol, new_sol);
			double tau = Problem::getTau(dis_sol);
			double nu = 1.0/tau;
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
			double factor = dt * nu * ( 1 - Problem::Pr ) / (dim+2);
			typename dis_t::Indices ind; // 初始化 ind 应该是 0
		  	for(int j = 0; j < dim; j++){
				ind[j] = 1;
				//新时刻的热通量预测值
				heat_flux[j] /= (1 + dt * nu * Problem::Pr);
		  		for(int i = 0; i < dim; i++){
					ind[i] += 2;
		  		    dis_sol(ind) +=  1 / (1.0 + dt * nu) * factor * heat_flux[j];
		  		    ind[i] -= 2;
		  		}
		  		ind[j] = 0;
		  	}
		}
		
		void ProjectRHS(const sol_t& rhs, sol_t& res, const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		    dis_t& dis_res = res[nx][ny];
		  		    const dis_t& dis_rhs = rhs[nx][ny];
		  		    Project(dis_rhs, dis_res);
		  		}
		  	}
		}

	    void print_config(){
		  	std::cout << "SISShakhovTime: CFL = " << BASE::CFL << std::endl;
		}

		void dump_data( std::ostream& os ) const {
		  	os.write((const char *)&(BASE::CFL), sizeof(double));
		}
};

template <typename RESIDUAL>
class SISShakhovSymGS : public SingleGridMethod<RESIDUAL>{
	public:
	    typedef RESIDUAL Residual;
	    typedef SingleGridMethod<Residual> BASE;
	    typedef typename BASE::Problem Problem;
	    typedef typename Residual::sol_t sol_t;
	    typedef typename sol_t::dis_t dis_t;
		typedef typename dis_t::Indices index_t;
	    typedef typename sol_t::mesh_t mesh_t; 
	    typedef std::shared_ptr<mesh_t> ptr_mesh_t;
	private:
	    sol_t res;
	public:
		SISShakhovSymGS(){};

		SISShakhovSymGS(const ptr_mesh_t& p_mesh, std::istream& is){
		  	is.read((char *)&(BASE::CFL), sizeof(double));
		}
		
		void solve(sol_t& sol, const sol_t& rhs, const mesh_t& mesh, 
		 		    const size_t steps = 1){
		 	for(size_t i = 0; i < steps; i++){
		 		symSweeping(sol, rhs, mesh);
		 	}
		}

		void symSweeping(sol_t& sol, const sol_t& rhs, const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x();
		 	size_t Ny = mesh.n_ele_y();
			res.Reinit(sol);
		 	ProjectRHS(rhs, res, mesh);
			#pragma omp parallel for				 
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
				   double dt = BASE::timeStep(nx, ny, sol, mesh);
				   Residual::convection(sol, mesh, res, nx, ny);
				   update(sol[nx][ny], res[nx][ny], dt);
				}
			}
			
			res.Reinit(sol);
		 	ProjectRHS(rhs, res, mesh);
			#pragma omp parallel for				 
			for(size_t ny = Ny; ny > 0; ny--){
				for(size_t nx = Nx; nx > 0; nx--){
				   double dt = BASE::timeStep(nx-1, ny-1, sol, mesh);
				   Residual::convection(sol, mesh, res, nx-1, ny-1);
				   update(sol[nx-1][ny-1], res[nx-1][ny-1], dt);
				}
			}
		}

		void update(dis_t& dis_sol,
				const dis_t& dis_res,
				double dt){
		  	dis_t new_sol;
		  	new_sol.Reinit(dis_sol);
			//1.考虑对流项
		  	dis_sol.Add(dt, dis_res);
			//2.标准投影，获得新时刻的物理量
			new_sol.ProjectToStdSpace(dis_sol);
		  	std::swap(dis_sol, new_sol);
			double tau = Problem::getTau(dis_sol);
			double nu = 1.0/tau;
			//取出数值通量
			const int dim = dis_t::dim;
			double heat_flux[dim];
		  	dis_sol.HeatFlux( heat_flux );
			//3.作系数调整，得到新时刻分布函数的展开式，首先处理BGK的部分
		  	typename dis_t::Iterator the_dis = dis_sol.begin(2);
		  	typename dis_t::Iterator end_dis = dis_sol.end();
		  	for(; the_dis != end_dis; the_dis++){
				*the_dis *= 1.0 / ( 1.0 + dt * nu);
		  	}
			//处理alpha=3的部分系数
			double factor = nu * ( 1 - Problem::Pr ) / (dim+2);
			typename dis_t::Indices ind; // 初始化 ind 应该是 0
		  	for(int j = 0; j < dim; j++){
				ind[j] = 1;
				//新时刻的热通量预测值
				heat_flux[j] = heat_flux[j] / (1 + dt * nu * Problem::Pr);
		  		for(int i = 0; i < dim; i++){
					ind[i] += 2;
		  		    dis_sol(ind) += dt /(1.0 + dt*nu) * factor * heat_flux[j];
		  		    ind[i] -= 2;
		  		}
		  		ind[j] = 0;
		  	}

			////另一种实现方式，但效率似乎差不多
			//typename dis_t::Iterator the_dis = dis_sol.begin(2);
		  	//typename dis_t::Iterator end_dis = dis_sol.begin(3);
		  	//for(; the_dis != end_dis; the_dis++){
			//	*the_dis *= 1.0 / ( 1.0 + dt * nu);
		  	//}
			////alpha=3的部分系数组成线性方程组并求解
			//double factor = (1 - Problem::Pr) / (dim+2);
			//double coe = dt * nu;
			//Eigen::MatrixXd A = Eigen::MatrixXd::Zero(3,3); 
			//A << 1-coe*(factor*3-1), -coe*factor, -coe*factor,
			//  -coe*factor*3, 1-coe*(factor-1), -coe*factor,
			//  -coe*factor*3, -coe*factor, 1-coe*(factor-1);
			//Eigen::MatrixXd A_inv = A.inverse();
			//
			//typename dis_t::Indices ind11(3,0,0), ind12(1,2,0), ind13(1,0,2);
			//Eigen::Vector3d b(dis_sol(ind11),dis_sol(ind12),dis_sol(ind13));
			//Eigen::Vector3d x = A_inv * b;
			//dis_sol(ind11) = x[0], dis_sol(ind12) = x[1], dis_sol(ind13) = x[2];
			//
			//typename dis_t::Indices ind21(0,3,0), ind22(2,1,0), ind23(0,1,2);
			//b << dis_sol(ind21), dis_sol(ind22), dis_sol(ind23);
			//x = A_inv * b;
			//dis_sol(ind21) = x[0], dis_sol(ind22) = x[1], dis_sol(ind23) = x[2];

			//typename dis_t::Indices ind31(0,0,3), ind32(2,0,1), ind33(0,2,1);
			//b << dis_sol(ind31), dis_sol(ind32), dis_sol(ind33);
			//x = A_inv * b;
			//dis_sol(ind31) = x[0], dis_sol(ind32) = x[1], dis_sol(ind33) = x[2];

			//typename dis_t::Indices ind(1,1,1);
			//dis_sol(ind) *= 1.0 / ( 1.0 + dt * nu);
			//
			//the_dis = dis_sol.begin(4);
		  	//end_dis = dis_sol.end();
		  	//for(; the_dis != end_dis; the_dis++){
			//	*the_dis *= 1.0 / ( 1.0 + dt * nu);
		  	//}
		}
		
		void ProjectRHS(const sol_t& rhs, sol_t& res, const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		    dis_t& dis_res = res[nx][ny];
		  		    const dis_t& dis_rhs = rhs[nx][ny];
		  		    Project(dis_rhs, dis_res);
		  		}
		  	}
		}

	    void print_config(){
		  	std::cout << "SISShakhovSymGS: CFL = " << BASE::CFL << std::endl;
		}

		void dump_data( std::ostream& os ) const {
		  	os.write((const char *)&(BASE::CFL), sizeof(double));
		}

};

template <typename RESIDUAL>
class EulerBased : public SingleGridMethod<RESIDUAL>{
	public:
	   typedef RESIDUAL Residual;
	   typedef SingleGridMethod<Residual> BASE;
	   typedef typename Residual::sol_t sol_t;
	   typedef typename sol_t::dis_t dis_t;
		typedef typename dis_t::Velocity velocity_t;
		typedef typename dis_t::value_type value_t;
	   typedef typename sol_t::mesh_t mesh_t;
	   typedef typename Residual::Euler_sol_t Euler_sol_t;
	   typedef typename Residual::Euler_sol_vec_t Euler_sol_vec_t;
	   typedef std::shared_ptr<mesh_t> ptr_mesh_t;
	private:
		enum FIMType { FIM1 = 1, FIM2 = 2, FIM3 = 3 };
		int fim_type;

		ExplicitTime<RESIDUAL> high_order_fim1_;
		SISShakhovTime<RESIDUAL> high_order_fim2_;
		SISShakhovSymGS<RESIDUAL> high_order_fim3_;

		ExplicitTimeMacEqs<RESIDUAL> euler_fim1_;
		ExplicitTimeMacEqs<RESIDUAL> euler_fim2_;
		SymGSMacEqs<RESIDUAL> euler_fim3_;
		
	public:

		EulerBased() : fim_type(FIM3) {}

		EulerBased(const ptr_mesh_t& p_mesh, std::istream& is) : fim_type(FIM3) {
			is.read((char *)&fim_type, sizeof(int));
			if(fim_type == FIM1){
				high_order_fim1_.load_data(is);
				euler_fim1_.load_data(is);
			}else if(fim_type == FIM2){
				high_order_fim2_.load_data(is);
				euler_fim2_.load_data(is);
			}else{
				fim_type = FIM3;
				high_order_fim3_.load_data(is);
				euler_fim3_.load_data(is);
			}
		}

		void read_config(const ConfigFile& cf){
			fim_type = readFIMType(cf);
			if(fim_type == FIM1){
				high_order_fim1_.read_config(cf);
				euler_fim1_.read_config(cf);
			}else if(fim_type == FIM2){
				high_order_fim2_.read_config(cf);
				euler_fim2_.read_config(cf);
			}else{
				fim_type = FIM3;
				high_order_fim3_.read_config(cf);
				euler_fim3_.read_config(cf);
			}
		}
		
		void dump_data( std::ostream& os ) const {
			os.write((const char *)&fim_type, sizeof(int));
			if(fim_type == FIM1){
				high_order_fim1_.dump_data(os);
				euler_fim1_.dump_data(os);
			}else if(fim_type == FIM2){
				high_order_fim2_.dump_data(os);
				euler_fim2_.dump_data(os);
			}else{
				high_order_fim3_.dump_data(os);
				euler_fim3_.dump_data(os);
			}
		}

		//这是一个更新高阶矩方程组的解的函数，矩方程的解只能在这个函数中发生变化，
		//在宏观方程的演化中应该避免出现矩方程组的解，以免在任何可能的地方对其作出
		//改变. 
		void Euler2HighOrder(sol_t& sol, 
				const Euler_sol_vec_t& Euler_sol_vec,
		  		const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
				 	//Euler_sol_t Euler_sol = Euler_sol_vec[nx][ny];
					//double rho = Euler_sol[0];
					//double u1 = Euler_sol[1] / rho;
					//double u2 = Euler_sol[2] / rho;
					//double u3 = Euler_sol[3] / rho;
					//double theta = EULERSOLVER.getTemperature(Euler_sol);
					////传递密度、速度、温度
					//sol[nx][ny].front() = rho;
					//typename dis_t::Velocity Center{u1, u2, u3};
					//value_t Scaling = sqrt(theta);
					//sol[nx][ny].Reinit(Center, Scaling);					

					//这里提供第二种宏观量的使用方式
					//对前后两步的宏观量做线性组合生成下一时刻的量
					//线性组合的权重是否会对算法的效率产生影响,目前没有系统测试过
					linearCombination(sol[nx][ny], Euler_sol_vec[nx][ny]);
		  		}
		  	}
		}

		void linearCombination(dis_t& dis,
				const Euler_sol_t& Euler_sol){
			//取中间步的宏观量
			value_t rho_old = dis.front();
			velocity_t v_old = dis.Center();
			value_t T_old = dis.TransTemperature();

			//取流体方程更新后的宏观量
			double rho = Euler_sol[0];
			double u1 = Euler_sol[1] / rho;
			double u2 = Euler_sol[2] / rho;
			double u3 = Euler_sol[3] / rho;
			double theta = getTemperature(Euler_sol);

			//做线性组合
			double tau = 1.0;
			rho = tau * rho + rho_old * (1-tau);
			u1 = tau * u1 + v_old[0] * (1-tau);
			u2 = tau * u2 + v_old[1] * (1-tau);
			u3 = tau * u3 + v_old[2] * (1-tau);
			theta = tau * theta + T_old * (1-tau);

			//重设分布函数
			dis.front() = rho;
			velocity_t Center{u1, u2, u3};
			value_t Scaling = sqrt(theta);
			dis.Reinit(Center, Scaling);
		}

		void HighOrder2Euler(const sol_t& sol,
				Euler_sol_vec_t& Euler_sol_vec,
		  		Euler_sol_vec_t& closure_vec, 
		  		const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		    const dis_t& dis = sol[nx][ny];
		  		    Euler_sol_t& Euler_sol = Euler_sol_vec[nx][ny];
		  		    Euler_sol_t& closure = closure_vec[nx][ny];
					getMacroVar(dis, Euler_sol, closure);
		  		}
		  	}
		}

		void getRHS(Euler_sol_vec_t& rhs_vec,
				const sol_t& rhs,
		  		const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	//注意，程序能够并行计算，与变量的范围密切相关
		  	//pv若定义在for循环之外，那么并行计算就会构成竞争关系，程序崩溃
		  	#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		    double pv[5];
		  		    const dis_t& dis_rhs = rhs[nx][ny];
		  		    dis_rhs.PrimitiveVars(pv);
		  		    Euler_sol_t& Euler_rhs = rhs_vec[nx][ny];
		  		    Euler_rhs.resize(5);
		  		    Euler_rhs[0] = pv[0];
		  		    Euler_rhs[1] = pv[0] * pv[1];
		  		    Euler_rhs[2] = pv[0] * pv[2];
		  		    Euler_rhs[3] = pv[0] * pv[3];
		  		    Euler_rhs[4] = pv[0] * (3.0/2.0 * pv[4] + 
		  		   			 (pv[1]*pv[1] + pv[2]*pv[2] + pv[3]*pv[3])/2.0 );
		  		}
		  	}
		}

		void solve(sol_t& sol, const sol_t& rhs, const mesh_t& mesh, 
		  		    const size_t steps = 1){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	Euler_sol_vec_t Euler_sol_vec(Nx, std::vector<Euler_sol_t>(Ny));
		  	Euler_sol_vec_t closure_vec(Nx, std::vector<Euler_sol_t>(Ny));
		  	Euler_sol_vec_t rhs_vec(Nx, std::vector<Euler_sol_t>(Ny));
		  	getRHS(rhs_vec, rhs, mesh);
		  	for(size_t i = 0; i < steps; i++){
				solveHighOrder(sol, rhs, mesh);
		  		HighOrder2Euler(sol, Euler_sol_vec, closure_vec, mesh);
				solveEuler(sol, Euler_sol_vec, closure_vec, rhs_vec, mesh);
		  		Euler2HighOrder(sol, Euler_sol_vec, mesh);
		  	}
		}

		void print_config(){
			std::cout << "FIM mode = FIM-" << fim_type << std::endl;
			if(fim_type == FIM1){
				std::cout << "High-order solver = ExplicitTime" << std::endl;
				std::cout << "Hydrodynamic solver = ExplicitTimeMacEqs" << std::endl;
				high_order_fim1_.print_config();
				euler_fim1_.print_config();
			}else if(fim_type == FIM2){
				std::cout << "High-order solver = SISShakhovTime" << std::endl;
				std::cout << "Hydrodynamic solver = ExplicitTimeMacEqs" << std::endl;
				high_order_fim2_.print_config();
				euler_fim2_.print_config();
			}else{
				std::cout << "High-order solver = SISShakhovSymGS" << std::endl;
				std::cout << "Hydrodynamic solver = SymGSMacEqs" << std::endl;
				high_order_fim3_.print_config();
				euler_fim3_.print_config();
			}
		}

	private:
		static int readFIMType(const ConfigFile& cf){
			try{
				int value = (int)cf.Value("Solver", "FIM");
				if(value >= FIM1 && value <= FIM3){
					return value;
				}
				throw std::string("Solver/FIM must be 1, 2, or 3");
			}catch(const std::string& e){
				if(e == "Solver/FIM must be 1, 2, or 3"){
					throw;
				}
			}
			return FIM3;
		}

		void solveHighOrder(sol_t& sol, const sol_t& rhs, const mesh_t& mesh){
			if(fim_type == FIM1){
				high_order_fim1_.solve(sol, rhs, mesh);
			}else if(fim_type == FIM2){
				high_order_fim2_.solve(sol, rhs, mesh);
			}else{
				high_order_fim3_.solve(sol, rhs, mesh);
			}
		}

		void solveEuler(sol_t& sol,
				Euler_sol_vec_t& Euler_sol_vec,
				const Euler_sol_vec_t& closure_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			if(fim_type == FIM1){
				euler_fim1_.solve(sol, Euler_sol_vec, closure_vec, rhs_vec, mesh);
			}else if(fim_type == FIM2){
				euler_fim2_.solve(sol, Euler_sol_vec, closure_vec, rhs_vec, mesh);
			}else{
				euler_fim3_.solve(sol, Euler_sol_vec, closure_vec, rhs_vec, mesh);
			}
		}

		void getMacroVar(const dis_t& dis, Euler_sol_t& Euler_sol, Euler_sol_t& closure){
			if(fim_type == FIM1){
				euler_fim1_.getMacroVar(dis, Euler_sol, closure);
			}else if(fim_type == FIM2){
				euler_fim2_.getMacroVar(dis, Euler_sol, closure);
			}else{
				euler_fim3_.getMacroVar(dis, Euler_sol, closure);
			}
		}

		double getTemperature(const Euler_sol_t& Euler_sol){
			if(fim_type == FIM1){
				return euler_fim1_.getTemperature(Euler_sol);
			}else if(fim_type == FIM2){
				return euler_fim2_.getTemperature(Euler_sol);
			}else{
				return euler_fim3_.getTemperature(Euler_sol);
			}
		}
};

# endif
