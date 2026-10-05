/**
 * @file EulerSolver.h
 * @brief Euler方程求解器类
 * @author Guanghan Li NUAA
 * @version 1
 * @date 2022-10-29
 */

# ifndef __EULERSOLVER__H__
# define __EULERSOLVER__H__

# include <vector> 
# include <cmath>
# include <NRxx/Transform.h>
# include "ConfigFile.h"
# include "config.h"
# include <Eigen/Core>
# include <Eigen/Dense>
# include <Eigen/StdVector>
# include <chrono>
# include <fstream>
# include <omp.h>

template <typename RESIDUAL>
class EulerSolverBase{
	public:
	    typedef typename RESIDUAL::mesh_t mesh_t;
	    typedef typename RESIDUAL::PROBLEM PROBLEM;
	    typedef typename RESIDUAL::sol_t sol_t;
	    typedef typename sol_t::dis_t dis_t;
	    typedef typename dis_t::value_type value_t;
	    typedef typename dis_t::Indices index_t;
	    typedef std::shared_ptr<mesh_t> ptr_mesh_t;
	    typedef typename Eigen::VectorXd Euler_sol_t;
	    typedef typename std::vector<Euler_sol_t> Euler_sol_vec_t;
	public:
	    double dt, CFL;
	    int maxSteps, M;
	    double res_norm, diff, macTol;
	public:
		    EulerSolverBase() : dt(1000),
		                        CFL(0.8),
		                        maxSteps(1),
		                        M(0),
		                        res_norm(0.),
		                        diff(1.),
		                        macTol(1e-8){};
		    EulerSolverBase(const ptr_mesh_t& /*p_mesh*/, std::istream& is) : EulerSolverBase() {
				is.read((char *)&CFL, sizeof(double));
				is.read((char *)&maxSteps, sizeof(int));
				is.read((char *)&macTol, sizeof(double));
	    }
	    ~EulerSolverBase() {};
	public:
		void read_config(const ptr_mesh_t& p_mesh, const ConfigFile& cf){
		  	CFL = (double)cf.Value( "Time", "CFL_HEs" );
		  	macTol = (double)cf.Value( "Error", "Tol" );
		  	maxSteps = (int)cf.Value("SingleGridMethod", "EulerEqsSolvingSteps");
		}

		void dump_data(std::ostream& os) const {
			os.write((const char *)&CFL, sizeof(double));
			os.write((const char *)&maxSteps, sizeof(int));
			os.write((const char *)&macTol, sizeof(double));
		}

		void print_config() const {
			std::cout << "Euler solver: CFL = " << CFL << std::endl;
			std::cout << "Euler solver: maxSteps = " << maxSteps << std::endl;
			std::cout << "Euler solver: macTol = " << macTol << std::endl;
		}

		void timeStep(const Euler_sol_vec_t& Euler_sol_vec, 
		  			    const mesh_t& mesh){
		  	size_t n_ele = mesh.n_ele();
			dt = 1000;
#pragma omp parallel for reduction(min:dt)
		  	for(size_t i = 0; i < n_ele; i++){
		  		 double dt_local = timeStep(Euler_sol_vec[i], mesh, i);
		  		 dt = (dt_local < dt) ? dt_local : dt;
		  	}
		}

		double timeStep(const Euler_sol_t& Euler_sol, 
		  					const mesh_t& mesh,
		  					unsigned int i){
		  	double u1 = Euler_sol[1] / Euler_sol[0];
		  	double theta = getTemperature(Euler_sol);
		  	double C_M = Transform<double>::MaxRootOfHermitePolynomial(M+1);
		  	return CFL * mesh.Length(i) / (fabs(u1) + C_M * sqrt(theta));
		}

		double getTemperature(const Euler_sol_t& Euler_sol){
		  double rho = Euler_sol[0];
		  double u1 = Euler_sol[1] / rho;
		  double u2 = Euler_sol[2] / rho;
		  double u3 = Euler_sol[3] / rho;
		  return ( Euler_sol[4] / rho - (u1*u1 + u2*u2 + u3*u3) / 2.0 )
		  	  * (2.0 / 3.0);
		}

		void getResidual(const sol_t& sol,
				const Euler_sol_vec_t& Euler_sol_vec,
		  		const Euler_sol_vec_t& closure_vec, 
		  		const mesh_t& mesh,
		  		Euler_sol_vec_t& res_vec){
		  	size_t n_ele = mesh.n_ele();
#pragma omp parallel for
		  	for(size_t nx = 0; nx < n_ele; nx++){
		  		 getResidual(sol, Euler_sol_vec, closure_vec, mesh, nx, res_vec[nx]);
		  	}
		}

		void getResidual(const sol_t& sol,
		  		 const Euler_sol_vec_t& Euler_sol_vec,
		  		 const Euler_sol_vec_t& closure_vec, 
		  		 const mesh_t& mesh,
		  		 size_t nx,
		  		 Euler_sol_t& res){
		  	#ifdef RECONSTRUCT
			secondOrderDisPri(sol, Euler_sol_vec, closure_vec, mesh, nx, res);
			//secondOrderDisCon(sol, Euler_sol_vec, closure_vec, mesh, nx, res);
			#else
		   firstOrderDis(sol, Euler_sol_vec, closure_vec, mesh, nx, res);
			#endif
		}

		void firstOrderDis(const sol_t& sol,
		  		 const Euler_sol_vec_t& Euler_sol_vec,
		  		 const Euler_sol_vec_t& closure_vec, 
				 const mesh_t& mesh,
		  		 size_t nx,
		  		 Euler_sol_t& res){
		    unsigned int n_ele = sol.size();
			 Euler_sol_t flux_l, flux_r;
		    Euler_sol_t sol_l, sol_m, sol_r, Eulercloure_l, Eulercloure_m, Eulercloure_r, temp;
		    if(nx == 0){
		  	  leftBoundaryCondition(sol[nx], Euler_sol_vec[nx], sol_l, Eulercloure_l);
		  	  sol_m = Euler_sol_vec[nx], Eulercloure_m = closure_vec[nx];
		  	  sol_r = Euler_sol_vec[nx+1], Eulercloure_r = closure_vec[nx+1];
		    }else if( nx == (n_ele - 1)){
		  	  rightBoundaryCondition(sol[nx], Euler_sol_vec[nx], sol_r, Eulercloure_r);
		  	  sol_m = Euler_sol_vec[nx], Eulercloure_m = closure_vec[nx];
		  	  sol_l = Euler_sol_vec[nx-1], Eulercloure_l = closure_vec[nx-1];
		    }else{
		  	 sol_l = Euler_sol_vec[nx-1], Eulercloure_l = closure_vec[nx-1];
		  	 sol_m = Euler_sol_vec[nx], Eulercloure_m = closure_vec[nx];
		  	 sol_r = Euler_sol_vec[nx+1], Eulercloure_r = closure_vec[nx+1];
		    }
		    
		    hllFlux(sol_l, sol_m, Eulercloure_l, Eulercloure_m, temp);
		    flux_l = temp;

		    hllFlux(sol_m, sol_r, Eulercloure_m, Eulercloure_r, temp);
		    flux_r = temp;

			//这里返回-R，注意不要随意更改符号，在多重网格算法中
		    //粗网格上残量的定义为r-R，我们这里又用了+=
		    res += -(flux_r - flux_l) / mesh.Length(nx);
		}

		void secondOrderDisPri(const sol_t& sol,
		  		 const Euler_sol_vec_t& Euler_sol_vec,
		  		 const Euler_sol_vec_t& closure_vec, 
		  		 const mesh_t& mesh,
		  		 size_t nx,
		  		 Euler_sol_t& res){
		    size_t n_ele = mesh.n_ele();
			Euler_sol_t flux_l, flux_r;
		    //临时变量，用来存边界处守恒向量的极限
		    Euler_sol_t lim_lr, lim_ml, lim_mr, lim_rl, dummy, temp;
		    //临时变量，用来存对应网格上的封闭项
		    Euler_sol_t Eulercloure_lr, Eulercloure_ml, Eulercloure_mr, Eulercloure_rl;
		    if(nx == 0){
			  reconstructPrimi(Euler_sol_vec, mesh, nx, lim_ml, lim_mr);
		  	  reconstruct(closure_vec, mesh, nx, Eulercloure_ml, Eulercloure_mr);

			  reconstructPrimi(Euler_sol_vec, mesh, nx+1, lim_rl, dummy);
		  	  reconstruct(closure_vec, mesh, nx+1, Eulercloure_rl, dummy);

		  	  dis_t dis_ml;
		  	  dis_ml = sol[0]; dis_ml *= value_t( 1.5 ); dis_ml.Add( -0.5, sol[1] );
		  	  leftBoundaryCondition(dis_ml, lim_ml, lim_lr, Eulercloure_lr);
		    }else if(nx == (n_ele - 1)){
			  reconstructPrimi(Euler_sol_vec, mesh, nx, lim_ml, lim_mr);
		  	  reconstruct(closure_vec, mesh, nx, Eulercloure_ml, Eulercloure_mr);

			  reconstructPrimi(Euler_sol_vec, mesh, nx-1, dummy, lim_lr);
		  	  reconstruct(closure_vec, mesh, nx-1, dummy, Eulercloure_lr);
		  	  
		  	  dis_t dis_mr;
		  	  dis_mr = sol[nx]; dis_mr *= value_t( 1.5 ); dis_mr.Add( -0.5, sol[nx-1] );
		  	  rightBoundaryCondition(dis_mr, lim_mr, lim_rl, Eulercloure_rl);
		    }else{
			  reconstructPrimi(Euler_sol_vec, mesh, nx-1, dummy, lim_lr);
		  	  reconstruct(closure_vec, mesh, nx-1, dummy, Eulercloure_lr);

			  reconstructPrimi(Euler_sol_vec, mesh, nx, lim_ml, lim_mr);
		  	  reconstruct(closure_vec, mesh, nx, Eulercloure_ml, Eulercloure_mr);

			  reconstructPrimi(Euler_sol_vec, mesh, nx+1, lim_rl, dummy);
		  	  reconstruct(closure_vec, mesh, nx+1, Eulercloure_rl, dummy);  
		    }
		    
		    hllFlux(lim_lr, lim_ml, Eulercloure_lr, Eulercloure_ml, flux_l);
		    hllFlux(lim_mr, lim_rl, Eulercloure_mr, Eulercloure_rl, flux_r);

			res += -(flux_r - flux_l) / mesh.Length(nx);
		}

		void reconstructPrimi(const Euler_sol_vec_t& Euler_sol_vec,
		  	  const mesh_t& mesh,
		  	  size_t nx,
		  	  Euler_sol_t& leftLim,
		  	  Euler_sol_t& rightLim){
			size_t n_ele = mesh.n_ele();
			Euler_sol_t primitiveLeft(5), primitiveMid(5), primitiveRight(5), slope(5);
			Euler_sol_t primiLeftLim(5), primiRightLim(5);
			if(nx == 0){
				conserToPrimi(Euler_sol_vec[nx], primitiveMid);
				conserToPrimi(Euler_sol_vec[nx+1], primitiveRight);
				primiLeftLim = 1.5 * primitiveMid - 0.5 * primitiveRight;
				primiRightLim = 0.5 * primitiveMid + 0.5 * primitiveRight;
			}else if(nx == n_ele-1){
				conserToPrimi(Euler_sol_vec[nx-1], primitiveLeft);
				conserToPrimi(Euler_sol_vec[nx], primitiveMid);
				primiLeftLim = 0.5 * primitiveMid + 0.5 * primitiveLeft;
				primiRightLim = 1.5 * primitiveMid - 0.5 * primitiveLeft;
			}else{
				double dx = mesh.Length(nx);
				conserToPrimi(Euler_sol_vec[nx-1], primitiveLeft);
				conserToPrimi(Euler_sol_vec[nx], primitiveMid);
				conserToPrimi(Euler_sol_vec[nx+1], primitiveRight);
				slope = (primitiveRight - primitiveLeft) / (2 * dx);
				primiLeftLim = primitiveMid - (0.5 * dx) * slope;
			   primiRightLim = primitiveMid + (0.5 * dx) * slope;
			}
			primiToConser(primiLeftLim, leftLim);
			primiToConser(primiRightLim, rightLim);
		}

		void conserToPrimi(const Euler_sol_t& Euler_sol,
				Euler_sol_t& primi){
		   double rho = Euler_sol[0];
		   double u1 = Euler_sol[1] / rho;
		   double u2 = Euler_sol[2] / rho;
		   double u3 = Euler_sol[3] / rho;
		   double theta = (Euler_sol[4] /rho - (u1*u1 + u2*u2 + u3*u3) / 2.0) / 1.5;
			primi.resize(5);
			primi[0] = rho;
			primi[1] = u1;
			primi[2] = u2;
			primi[3] = u3;
			primi[4] = sqrt(theta);
		}

		void primiToConser(const Euler_sol_t& primi,
				Euler_sol_t& conser){
			double rho = primi[0];
		   double u1 = primi[1];
		   double u2 = primi[2];
		   double u3 = primi[3];
		   double theta = primi[4] * primi[4]; 
			conser.resize(5);
			conser[0] = rho;
		  	conser[1] = rho * u1;
		  	conser[2] = rho * u2;
		  	conser[3] = rho * u3;
		  	conser[4] = rho * (1.5 * theta + (u1*u1 + u2*u2 + u3*u3) / 2.0);
		}

		void secondOrderDisCon(const sol_t& sol,
		  		 const Euler_sol_vec_t& Euler_sol_vec,
		  		 const Euler_sol_vec_t& closure_vec, 
		  		 const mesh_t& mesh,
		  		 size_t nx,
		  		 Euler_sol_t& res){
		    size_t n_ele = mesh.n_ele();
			Euler_sol_t flux_l, flux_r;
		    //临时变量，用来存边界处守恒向量的极限
		    Euler_sol_t lim_lr, lim_ml, lim_mr, lim_rl, dummy, temp;
		    //临时变量，用来存对应网格上的封闭项
		    Euler_sol_t Eulercloure_lr, Eulercloure_ml, Eulercloure_mr, Eulercloure_rl;
		    if(nx == 0){
		  	  reconstruct(Euler_sol_vec, mesh, nx, lim_ml, lim_mr);
		  	  reconstruct(closure_vec, mesh, nx, Eulercloure_ml, Eulercloure_mr);

		  	  reconstruct(Euler_sol_vec, mesh, nx+1, lim_rl, dummy);
		  	  reconstruct(closure_vec, mesh, nx+1, Eulercloure_rl, dummy);

		  	  dis_t dis_ml;
		  	  dis_ml = sol[0]; dis_ml *= value_t( 1.5 ); dis_ml.Add( -0.5, sol[1] );
		  	  leftBoundaryCondition(dis_ml, lim_ml, lim_lr, Eulercloure_lr);
		    }else if(nx == (n_ele - 1)){
		  	  reconstruct(Euler_sol_vec, mesh, nx, lim_ml, lim_mr);
		  	  reconstruct(closure_vec, mesh, nx, Eulercloure_ml, Eulercloure_mr);

		  	  reconstruct(Euler_sol_vec, mesh, nx-1, dummy, lim_lr);
		  	  reconstruct(closure_vec, mesh, nx-1, dummy, Eulercloure_lr);
		  	  
		  	  dis_t dis_mr;
		  	  dis_mr = sol[nx]; dis_mr *= value_t( 1.5 ); dis_mr.Add( -0.5, sol[nx-1] );
		  	  rightBoundaryCondition(dis_mr, lim_mr, lim_rl, Eulercloure_rl);
		    }else{
		  	  reconstruct(Euler_sol_vec, mesh, nx-1, dummy, lim_lr);
		  	  reconstruct(closure_vec, mesh, nx-1, dummy, Eulercloure_lr);

		  	  reconstruct(Euler_sol_vec, mesh, nx, lim_ml, lim_mr);
		  	  reconstruct(closure_vec, mesh, nx, Eulercloure_ml, Eulercloure_mr);

		  	  reconstruct(Euler_sol_vec, mesh, nx+1, lim_rl, dummy);
		  	  reconstruct(closure_vec, mesh, nx+1, Eulercloure_rl, dummy);  
		    }
		    
		    hllFlux(lim_lr, lim_ml, Eulercloure_lr, Eulercloure_ml, flux_l);
		    hllFlux(lim_mr, lim_rl, Eulercloure_mr, Eulercloure_rl, flux_r);

			res += -(flux_r - flux_l) / mesh.Length(nx);
		}

		void reconstruct(const Euler_sol_vec_t& Euler_sol_vec,
		  	  const mesh_t& mesh,
		  	  size_t nx,
		  	  Euler_sol_t& leftLim,
		  	  Euler_sol_t& rightLim){
		    size_t n_ele = mesh.n_ele();
		    Euler_sol_t slope;
		    if(nx == 0){
		        leftLim = 1.5 * Euler_sol_vec[nx] - 0.5 * Euler_sol_vec[nx+1];
		        rightLim = 0.5 * Euler_sol_vec[nx] + 0.5 * Euler_sol_vec[nx+1];
		    }else if(nx == n_ele-1){
		        leftLim = 0.5 * Euler_sol_vec[nx] + 0.5 * Euler_sol_vec[nx-1];
		        rightLim = 1.5 * Euler_sol_vec[nx] - 0.5 * Euler_sol_vec[nx-1];
		    }else{
		        double dx = mesh.Length(nx);
		        slope = (Euler_sol_vec[nx+1] - Euler_sol_vec[nx-1]) / (2 * dx );
		        leftLim = Euler_sol_vec[nx] - (0.5 * dx) * slope;
		        rightLim = Euler_sol_vec[nx] + (0.5 * dx) * slope;
		    }
		}

		void leftBoundaryCondition(const dis_t& dis,
		  	  const Euler_sol_t& Euler_sol,
		  	  Euler_sol_t& Euler_sol_ghost,
		  	  Euler_sol_t& closure_ghost){
		    dis_t dis_tmp = dis;
		    //用与边界邻接单元上的宏观量更新分布函数
		    //目前测试表明，是否更新，对最终迭代影响不大，所以暂时注释在这里
		    updateMacroVar(dis_tmp, Euler_sol);
		    dis_t dis_ghost;
		    PROBLEM::BoundaryValue(dis_ghost, dis_tmp, false);
		    getMacroVar(dis_ghost, Euler_sol_ghost, closure_ghost);

			//Kn=0.001时使用NS方程进行求解，使用稳态解的边界条件，固定值
			//Euler_sol_ghost.resize(5);
			//Euler_sol_ghost << 1.0510268654089585,0.00011071639981621082, 
			//	-0.65660726302748074, 0, 1.7835309699433728 ;
			//closure_ghost.resize(5);
			//closure_ghost << 1.0522867139046852, 3.5773017161983617e-05, 
			//	-0.0041310224623859419, -0, -0.0025675442011927284 ;
		}

		void rightBoundaryCondition(const dis_t& dis,
		  	  const Euler_sol_t& Euler_sol,
		  	  Euler_sol_t& Euler_sol_ghost,
		  	  Euler_sol_t& closure_ghost){
			dis_t dis_tmp = dis;
		   updateMacroVar(dis_tmp, Euler_sol);
		   dis_t dis_ghost;
		   PROBLEM::BoundaryValue(dis_ghost, dis_tmp, true);
		   getMacroVar(dis_ghost, Euler_sol_ghost, closure_ghost);	

			//Kn=0.001时使用NS方程进行求解，使用稳态解的边界条件，固定值
			//Euler_sol_ghost.resize(5);
			//Euler_sol_ghost << 1.0510268658060797, -0.00011071638288025833,
			//	0.65660726486101872, 0, 1.7835309710853635;
			//closure_ghost.resize(5);
			//closure_ghost << 1.0522867139540339, 3.5773023392540446e-05 , 
			//	-0.0041310209485781336, 0, 0.0025675434895811753;
		}

		void updateMacroVar(dis_t& dis,
		  	  const Euler_sol_t& Euler_sol){
		    double rho = Euler_sol[0];
		    double u1 = Euler_sol[1] / rho;
		    double u2 = Euler_sol[2] / rho;
		    double u3 = Euler_sol[3] / rho;
		    double theta = (Euler_sol[4] /rho - (u1*u1 + u2*u2 + u3*u3) / 2.0) / 1.5;

		    dis.front() = rho;
		    typename dis_t::Velocity Center{u1, u2, u3};
		    value_t Scaling = sqrt(theta);
		    dis.Reinit(Center, Scaling);
		}

		void getMacroVar(const dis_t& dis,
		  		 Euler_sol_t& Euler_sol,
		  		 Euler_sol_t& closure){
		  	static index_t ind_0(2,0,0), ind_1(1,1,0), ind_2(1,0,1), 
		  			  ind_3(3,0,0), ind_4(1,2,0), ind_5(1,0,2);
		  	Euler_sol.resize(5), closure.resize(5);
		  	//从分布函数中取出原始变量
		  	double pv[5];  
		  	dis.PrimitiveVars(pv);
			//使用原始变量生成守恒向量
		  	Euler_sol[0] = pv[0];
		  	Euler_sol[1] = pv[0] * pv[1];
		  	Euler_sol[2] = pv[0] * pv[2];
		  	Euler_sol[3] = pv[0] * pv[3];
		  	Euler_sol[4] = pv[0] * (3.0/2.0 * pv[4] + 
		  			  (pv[1]*pv[1] + pv[2]*pv[2] + pv[3]*pv[3])/2.0 );

		  	//closure terms: pressure, stress tensor, heat flux
		  	closure[0] = pv[0] * pv[4]; 
		  	closure[1] = 2 * dis(ind_0);
		  	closure[2] = dis(ind_1);
		  	closure[3] = dis(ind_2);
		  	closure[4] = 3 * dis(ind_3) + dis(ind_4) + dis(ind_5);
		}

		void hllFlux(const Euler_sol_t& sol_l, 
				const Euler_sol_t& sol_r,
		  		const Euler_sol_t& closure_l, 
		  		const Euler_sol_t& closure_r,
		  		Euler_sol_t& flux){
		  	double S_L, S_R;
		  	flux.resize(5);
		  	Eigen::VectorXd flux_1(5), flux_2(5);			
		  	flux_1.resize(5), flux_2.resize(5);
		  	double u_l = sol_l[1] / sol_l[0];
		  	double u_r = sol_r[1] / sol_r[0];
		  	double theta_l = getTemperature(sol_l);
		  	double theta_r = getTemperature(sol_r);
		  	double C_M = Transform<double>::MaxRootOfHermitePolynomial( M+1 );
		  	S_L = std::min(u_l - C_M*sqrt(theta_l), u_r - C_M*sqrt(theta_r));
		  	S_R = std::max(u_l + C_M*sqrt(theta_l), u_r + C_M*sqrt(theta_r));
		  	fluxFunction(sol_l, closure_l, flux_1);
		  	fluxFunction(sol_r, closure_r, flux_2);
		  	if(S_L >= 0){
				flux = flux_1;
		  	}else if (S_L < 0 && 0 < S_R){
		  		flux = (S_R * flux_1 - S_L * flux_2 + S_L * S_R * (sol_r - sol_l)) / (S_R - S_L);				
		  	}else if (S_R <= 0){
		  		flux = flux_2;
		  	}
		}

		void fluxFunction(const Euler_sol_t& Euler_sol, 
				const Euler_sol_t& closure, 
		  		Euler_sol_t& F){
		  	double rho = Euler_sol[0];
		  	double u1 = Euler_sol[1] / rho;
		  	double u2 = Euler_sol[2] / rho;
		  	double u3 = Euler_sol[3] / rho;
		  	double theta = getTemperature(Euler_sol);
		  	double p = rho * theta;
		  	F[0] = rho * u1;
		  	F[1] = rho * u1 * u1 + p + closure[1];
		  	F[2] = rho * u1 * u2 + closure[2];
		  	F[3] = rho * u1 * u3 + closure[3];
		  	F[4] = u1 * (Euler_sol[4] + p) + u1*closure[1] + u2*closure[2] + 
		  		  u3*closure[3] + closure[4];
		}

		void updateSolution(Euler_sol_vec_t& Euler_sol_vec, 
		  	    const Euler_sol_vec_t& res_vec,
		  		 const mesh_t& mesh){
		  	size_t n_ele = mesh.n_ele();
#pragma omp parallel for
		  	for(size_t nx = 0; nx < n_ele; nx++){
		  		 updateSolution(Euler_sol_vec[nx], res_vec[nx]);
		  	}
		}
		  	
		void updateSolution(Euler_sol_t& Euler_sol,
		  	   const Euler_sol_t& res){
		    Euler_sol += dt * res;
		}

		double getL1Norm(const Euler_sol_vec_t& Euler_sol_vec){
		     int n_ele = Euler_sol_vec.size();
		  	double normL1 = 0;
#pragma omp parallel for reduction(+:normL1)
		  	for(int nx = 0; nx < n_ele; nx++){
		  		double temp = getL1Norm(Euler_sol_vec[nx]);
		  		normL1 += temp; 
		  	}
		  	return normL1;
		}

		double getL1Norm(const Euler_sol_t& Euler_sol){
		     int n_ele = Euler_sol.size();
		     double normL1 = 0;
		  	for(int k = 0; k < n_ele; k++){
		  		normL1 += abs(Euler_sol[k]);
		  	}
		  	return normL1;
		}

	    double getL2Norm(const Euler_sol_vec_t& Euler_sol_vec){
		     int n_ele = Euler_sol_vec.size();
		  	double normL2 = 0;
#pragma omp parallel for reduction(+:normL2)
		  	for(int nx = 0; nx < n_ele; nx++){
		  		 double temp = getL2Norm(Euler_sol_vec[nx]);
		  		 normL2 += temp; 
		  	}
		  	return sqrt(normL2/(5*n_ele));
		}

		double getL2Norm(const Euler_sol_t& Euler_sol){
		     int n_ele = Euler_sol.size();
		     double normL2 = 0;
		  	for(int k = 0; k < n_ele; k++){
		  		 normL2 += pow(Euler_sol[k], 2);
		  	}
		  	return normL2;
		}

		double getInfNorm(const Euler_sol_vec_t& Euler_sol_vec){
		     int n_ele = Euler_sol_vec.size();
		  	double normInf = 0;
#pragma omp parallel
				{
					double thread_max = 0;
#pragma omp for nowait
					for(int nx = 0; nx < n_ele; nx++){
						double temp = getInfNorm(Euler_sol_vec[nx]);
						thread_max = (thread_max < temp) ? temp : thread_max;
					}
#pragma omp critical
					{
						normInf = (normInf < thread_max) ? thread_max : normInf;
					}
				}
		  	return normInf;
		}

		double getInfNorm(const Euler_sol_t& Euler_sol){
		     int n_ele = Euler_sol.size();
		     double normInf = 0;
		  	for(int k = 0; k < n_ele; k++){
		  		 double temp = abs(Euler_sol[k]);
		  		 normInf = (normInf < temp) ? temp : normInf; 
		  	}
		  	return normInf;
		}
		  	
		double getDifNorm(const Euler_sol_vec_t& Euler_sol_vec_old,
		  		const Euler_sol_vec_t& Euler_sol_vec_new,
		  		const mesh_t& mesh){
		     int n_ele = mesh.n_ele();
		     Euler_sol_vec_t Euler_sol_vec_dif(n_ele);
#pragma omp parallel for
		  	for(int nx = 0; nx < n_ele; nx++){
		  		 int n_var = Euler_sol_vec_old[nx].size();
		  		 Euler_sol_vec_dif[nx].resize(n_var);
		  		 for(int k = 0; k < n_var; k++){
		  			  Euler_sol_vec_dif[nx][k] = (Euler_sol_vec_old[nx][k] -
		  				Euler_sol_vec_new[nx][k]);
		  		 }
		  	}
		  	//return getInfNorm(Euler_sol_vec_dif);
		  	return getL1Norm(Euler_sol_vec_dif);
		  	//return getL2Norm(Euler_sol_vec_dif);
		}

		double getResNorm(const sol_t& sol,
		  	    const Euler_sol_vec_t& Euler_sol_vec,
		  		 const Euler_sol_vec_t& closure_vec, 
		  		 const mesh_t& mesh){
		     int n_ele = mesh.n_ele();
		     Euler_sol_vec_t res_vec(n_ele, Eigen::VectorXd::Zero(5));
		     //res_vec.resize(n_ele, Eigen::VectorXd::Zero(5));
		     //res_vec.resize(n_ele, std::vector<double>(5,0));
		     getResidual(sol, Euler_sol_vec, closure_vec, mesh, res_vec);
		     return getL2Norm(res_vec);
		}

		void printMacro(const Euler_sol_t& Euler_sol){
		  	double rho = Euler_sol[0];
		  	double u1 = Euler_sol[1] / rho;
		  	double u2 = Euler_sol[2] / rho;
		  	double u3 = Euler_sol[3] / rho;
		  	double theta = getTemperature(Euler_sol);
		  	std::cout << "rho = " << rho << std::endl
		  		 << "u1 = " << u1 << std::endl
		  		 << "u2 = " << u2 << std::endl
		  		 << "u3 = " << u3 << std::endl	
		  		 << "theta = " << theta << std::endl;
		}

		void printClosure(const Euler_sol_t& closure){
		  	std::cout << "p = " << closure[0] << std::endl
		  		 << "sigma11 = " << closure[1] << std::endl
		  		 << "sigma12 = " << closure[2] << std::endl
		  		 << "sigma13 = " << closure[3] << std::endl	
		  		 << "q1 = " << closure[4] << std::endl;
		}

		void ghostMacro(Euler_sol_t& Euler_sol_ghost, 
		  		 const Euler_sol_t& Euler_sol){
		  	double rho = Euler_sol[0];
		  	double u1 = -Euler_sol[1] /rho;
		  	double u2 = Euler_sol[2] /rho;
		  	double u3 = Euler_sol[3] /rho;
		  	double theta = getTemperature(Euler_sol);
		  	Euler_sol_ghost[0] = rho;
		  	Euler_sol_ghost[1] = rho * u1;
		  	Euler_sol_ghost[2] = rho * u2;
		  	Euler_sol_ghost[3] = rho * u3;
		  	Euler_sol_ghost[4] = rho * (3.0/2.0 * theta 
		  			  + (u1*u1 + u2*u2 + u3*u3) / 2.0);
		}

		virtual void solve(const sol_t& sol,
					 Euler_sol_vec_t& Euler_sol_vec, 
					 const Euler_sol_vec_t& closure_vec, 
					 const Euler_sol_vec_t& rhs_vec,
					 const mesh_t& mesh,
					 const size_t steps = 1) {}
};

template <typename RESIDUAL>
class EulerTimeForward : public EulerSolverBase<RESIDUAL>{
	 public:
		  typedef EulerSolverBase<RESIDUAL> Base;
		  typedef typename Base::sol_t sol_t;
		  typedef typename Base::dis_t dis_t;
		  typedef typename Base::mesh_t mesh_t;
		  typedef typename Base::PROBLEM PROBLEM;
		  typedef typename Base::Euler_sol_t Euler_sol_t;
		  typedef typename Base::Euler_sol_vec_t Euler_sol_vec_t;
	 public:
		  EulerTimeForward() {}
		  EulerTimeForward(const typename Base::ptr_mesh_t& p_mesh, std::istream& is) : Base(p_mesh, is) {}
		  void solve(const sol_t& sol,
					 Euler_sol_vec_t& Euler_sol_vec, 
					 const Euler_sol_vec_t& closure_vec, 
					 const Euler_sol_vec_t& rhs_vec,
					 const mesh_t& mesh){
				Base::M = sol[0].GetOrder();
				int steps = 0; double dif = 1.0;
				//设置宏观方程组的终止条件，迭代到最大步数，或者相邻两步解相差很小
				//迭代步数通过一些测试获得，相差的小量与矩方程组Tol一致
				//由于两步解的差的范数是单调下降的，可以想到前期基本都是由最大步数来
				//控制，在最后几步的时候宏观方程只求解一步即可满足条件
				//如果能有一个随着迭代进行自适应调整宏观方程组求解步数的策略会更好
				while((steps < Base::maxSteps) && (dif > Base::macTol)){
					 Euler_sol_vec_t Euler_sol_vec_old = Euler_sol_vec;
					 timeOneStep(sol, Euler_sol_vec, closure_vec, rhs_vec, mesh);
					 steps++;
					 dif = Base::getDifNorm(Euler_sol_vec_old,
							 Euler_sol_vec, mesh);
					}
			  }
		  
		  void timeOneStep(const sol_t& sol,
					 Euler_sol_vec_t& Euler_sol_vec, 
					 const Euler_sol_vec_t& closure_vec, 
					 const Euler_sol_vec_t& rhs_vec,
					 const mesh_t& mesh){
			  Euler_sol_vec_t res_vec = rhs_vec;
			  Base::timeStep(Euler_sol_vec, mesh);
			  Base::getResidual(sol, Euler_sol_vec, closure_vec, mesh, res_vec);
			  Base::updateSolution(Euler_sol_vec, res_vec, mesh);
		  }
};

template <typename RESIDUAL>
class EulerGSIteration : public EulerSolverBase<RESIDUAL>{
	 public:
		  typedef EulerSolverBase<RESIDUAL> Base;
		  typedef typename Base::sol_t sol_t;
		  typedef typename Base::dis_t dis_t;
		  typedef typename Base::mesh_t mesh_t;
		  typedef typename Base::PROBLEM PROBLEM;
		  typedef typename Base::Euler_sol_t Euler_sol_t;
		  typedef typename Base::Euler_sol_vec_t Euler_sol_vec_t;
	 public:
		  EulerGSIteration() {}
		  EulerGSIteration(const typename Base::ptr_mesh_t& p_mesh, std::istream& is) : Base(p_mesh, is) {}
		  void solve(const sol_t& sol,
		   		 Euler_sol_vec_t& Euler_sol_vec, 
		   		 const Euler_sol_vec_t& closure_vec, 
		   		 const Euler_sol_vec_t& rhs_vec,
		   		 const mesh_t& mesh,
				 const size_t step = 1){
		   	Base::M = sol[0].GetOrder();
		   	int steps = 0; double dif = 1.0;
		   	while((steps < Base::maxSteps) && (dif > Base::macTol)){
		   		Euler_sol_vec_t Euler_sol_vec_old = Euler_sol_vec;
		   		gsOneStep(sol, Euler_sol_vec, closure_vec, rhs_vec, mesh);
		   		steps++;
		   		dif = Base::getDifNorm(Euler_sol_vec_old, Euler_sol_vec, mesh);
		   	}
		  }

		  void gsOneStep(const sol_t& sol,
		   		 Euler_sol_vec_t& Euler_sol_vec, 
		   		 const Euler_sol_vec_t& closure_vec, 
		   		 const Euler_sol_vec_t& rhs_vec,
		   		 const mesh_t& mesh){
			  leftToRight(sol, Euler_sol_vec, closure_vec, rhs_vec, mesh);
			  rightToLeft(sol, Euler_sol_vec, closure_vec, rhs_vec, mesh);
		  }

		  void leftToRight(const sol_t& sol,
		   		 Euler_sol_vec_t& Euler_sol_vec, 
		   		 const Euler_sol_vec_t& closure_vec, 
		   		 const Euler_sol_vec_t& rhs_vec,
		   		 const mesh_t& mesh){
		      size_t n_ele = mesh.n_ele();
			  //已经有好几次在残量这个地方出错了，有的是残量没有正确赋值，导致算法
			  //发散。还有残量没有初始化为右端项，导致作为多重网格的光滑子时残量降
			  //不下来。这个地方特别注意
		      Euler_sol_vec_t res_vec = rhs_vec;
		      for(size_t nx = 0; nx < n_ele; nx++){
		   		 localCellOneStep(sol, Euler_sol_vec, closure_vec, rhs_vec, res_vec, mesh, nx);
			  }
		  }

		  void rightToLeft(const sol_t& sol,
		   		 Euler_sol_vec_t& Euler_sol_vec, 
		   		 const Euler_sol_vec_t& closure_vec, 
		   		 const Euler_sol_vec_t& rhs_vec,
		   		 const mesh_t& mesh){
		      size_t n_ele = mesh.n_ele();
		      Euler_sol_vec_t res_vec = rhs_vec;
		      for(size_t nx = n_ele; nx > 0; nx--){
		   		 localCellOneStep(sol, Euler_sol_vec, closure_vec, rhs_vec, res_vec, mesh, nx-1);
			  }
		  }

		  void localCellOneStep(const sol_t& sol,
		   		 Euler_sol_vec_t& Euler_sol_vec, 
		   		 const Euler_sol_vec_t& closure_vec, 
		   		 const Euler_sol_vec_t& rhs_vec,
		   		 Euler_sol_vec_t& res_vec,
		   		 const mesh_t& mesh,
		   		 size_t nx){
		      Base::dt = Base::timeStep(Euler_sol_vec[nx], mesh, nx);
		   	  Base::getResidual(sol, Euler_sol_vec, closure_vec, mesh, nx,
		   			res_vec[nx]);
		   	  Base::updateSolution(Euler_sol_vec[nx], res_vec[nx]);
		  }
};

template <typename RESIDUAL>
class NavierStokes : public EulerSolverBase<RESIDUAL>{
	public:
		typedef EulerSolverBase<RESIDUAL> Base;
		typedef typename Base::sol_t sol_t;
		typedef typename Base::dis_t dis_t;
		typedef typename Base::mesh_t mesh_t;
		typedef typename Base::PROBLEM PROBLEM;
		typedef typename Base::Euler_sol_t Euler_sol_t;
		typedef typename Base::Euler_sol_vec_t Euler_sol_vec_t;
	public:
		void solve(const sol_t& sol,
				Euler_sol_vec_t& Euler_sol_vec, 
				Euler_sol_vec_t& closure_vec, 
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh, 
				const int steps = 1){
			size_t n_ele = mesh.n_ele();
			Base::M = sol[0].GetOrder();
			//Euler_sol_vec_t res_vec(n_ele, Eigen::VectorXd::Zero(5));
			int NSsteps = 0;
			double resNorm;
			//当NS方程解的残量下降到1e-8时停止迭代
			do{
				 NSsteps++;
				 Base::timeStep(Euler_sol_vec, mesh);
				 NSclosure(Euler_sol_vec, closure_vec, mesh);
				 Euler_sol_vec_t res_vec(n_ele, Eigen::VectorXd::Zero(5));
				 Base::getResidual(sol, Euler_sol_vec, closure_vec, mesh, res_vec);
				 Base::updateSolution(Euler_sol_vec, res_vec, mesh);
				 resNorm = Base::getL2Norm(res_vec);
				 std::cout << "steps = " << NSsteps << '\t' << "res_norm = " << resNorm << std::endl;
			}while(resNorm>1e-8);
			outputSolution(Euler_sol_vec, closure_vec, mesh);
			exit(0);
		}
		  
		void NSclosure(const Euler_sol_vec_t& Euler_sol_vec,
		  	  Euler_sol_vec_t& closure_vec,
		  	  const mesh_t& mesh){
		    size_t n_ele = mesh.n_ele();
#pragma omp parallel for
		    for(size_t i = 0; i < n_ele; i++){
				Euler_sol_t Euler_sol_l, Euler_sol_m, Euler_sol_r;
				if(i == 0){
					//ghost单元守恒量，Kn=0.001
					Euler_sol_t left(5);
					left << 1.0510268654089585,0.00011071639981621082, 
					-0.65660726302748074, 0, 1.7835309699433728;
					Euler_sol_l = left;
					Euler_sol_m = Euler_sol_vec[i];
					Euler_sol_r = Euler_sol_vec[i+1];
		  	  	}else if(i == n_ele-1){
					Euler_sol_t right(5);
					right << 1.0510268658060797, -0.00011071638288025833,
					0.65660726486101872, 0, 1.7835309710853635;

					Euler_sol_l = Euler_sol_vec[i-1];
					Euler_sol_m = Euler_sol_vec[i];
					Euler_sol_r = right;
		  	  	}else{
					Euler_sol_l = Euler_sol_vec[i-1];
					Euler_sol_m = Euler_sol_vec[i];
					Euler_sol_r = Euler_sol_vec[i+1];
				}
				getNSClosure(closure_vec[i], Euler_sol_l, Euler_sol_m, Euler_sol_r, mesh);
		    }
		}
		
		//中心差分
		void getNSClosure(Euler_sol_t& closure,
		  	  const Euler_sol_t& solL, 
		  	  const Euler_sol_t& solM,
		  	  const Euler_sol_t& solR,
		  	  const mesh_t& mesh){
		    Euler_sol_t macL(5), macM(5), macR(5);
		    getMac(macL, solL);
		 	getMac(macM, solM);
		    getMac(macR, solR);
		    //网格单元长度，均匀网格
		    double h = mesh.Length(0);
		    //求压力
		    double p = macM[0] * macM[4]; 
		    //粘性指数
		    double w = 0.81;
		    //平均碰撞频率
		    double nu = macM[0] * pow( macM[4], 1-w ) / PROBLEM::Kn *
		  		 sqrt(M_PI/2);
		    //p
		    closure[0] = p;
		    //sigma_{11}
		    closure[1] = -(4.0/3.0) * (p/nu) * (macR[1] - macL[1]) / (2*h);
		    //sigma_{12}
		    closure[2] = - (p/nu) * (macR[2] - macL[2]) / (2*h);
		    //sigma_{13}
		    closure[3] = - (p/nu) * (macR[3] - macL[3]) / (2*h);
		    //q_{1}
		    closure[4] = -(5.0/2.0) * (p/nu) * (macR[4] - macL[4]) / (2*h);
		}

		void getMac(Euler_sol_t& mac, const Euler_sol_t& Euler_sol){
			mac[0] = Euler_sol[0];
			mac[1] = Euler_sol[1] / mac[0];
			mac[2] = Euler_sol[2] / mac[0];
			mac[3] = Euler_sol[3] / mac[0];
			mac[4] = Base::getTemperature(Euler_sol);
		}

		void outputSolution(const Euler_sol_vec_t& Euler_sol_vec,
				const Euler_sol_vec_t& closure_vec,
				const mesh_t& mesh){
			size_t n_ele = mesh.n_ele();
			std::ofstream os("NS.dat");
			os.precision(12);
			os << "% variables= \"x\", \"rho\", "
				<< " \"u\", \"v\", \"w\", \"T\", \"p\" " 
			   << " \"sigma11\", \"sigma12\", \"sigma13\", \"q1\" " 
				<< std::endl;
			Euler_sol_t mac(5);
			for(size_t nx = 0; nx < n_ele; nx++){
				 os << mesh.Center(nx) << '\t';
				 //宏观量
				 getMac(mac, Euler_sol_vec[nx]);
				 os << mac[0] << '\t';
				 os << mac[1] << '\t';
				 os << mac[2] << '\t';
				 os << mac[3] << '\t';
				 os << mac[4] << '\t';
				 //封闭项
				 os << closure_vec[nx][0] << '\t';
				 os << closure_vec[nx][1] << '\t';
				 os << closure_vec[nx][2] << '\t';
				 os << closure_vec[nx][3] << '\t';
				 os << closure_vec[nx][4] << std::endl;
			}
		}
};

template <typename RESIDUAL>
class NewtonIte : public EulerSolverBase<RESIDUAL>{
	public:
	   typedef EulerSolverBase<RESIDUAL> Base;
	   typedef typename Base::sol_t sol_t;
	   typedef typename Base::dis_t dis_t;
	   typedef typename Base::mesh_t mesh_t;
		typedef std::shared_ptr<mesh_t> ptr_mesh_t;
	   typedef typename Base::PROBLEM PROBLEM;
	   typedef typename Base::Euler_sol_t Euler_sol_t;
	   typedef typename Base::Euler_sol_vec_t Euler_sol_vec_t;
		typedef typename Eigen::MatrixXd mat_t;
		typedef typename std::vector<mat_t> mat_vec_t;
	private:
		double perturbation_ratio;    /**< 计算 Jacobian 矩阵数值微分扰动比例 */
  		double lambda;                /**< 正则化 Jacobian 矩阵的比例系数 */
  		double mu;                    /**< 更新解步长, 准确的 Newton 迭代应为 1 */
		size_t s1, s2;            /**< 线性多重网格光滑步数，最粗网格迭代步数 */
		size_t maxMGSteps;    		  /**< 线性多重网格最大迭代步数 */
		size_t coarsestGrid;		  /**< 线性多重网格最粗网格规模 */
	public:
	    NewtonIte() : perturbation_ratio(1e-6), lambda(0.0), mu(1.0),
		s1(2), s2(2), coarsestGrid(8), maxMGSteps(1) {};
	    ~NewtonIte() {};
	public:
		class JacMatrix{
			public:
				mat_vec_t cent, left, right, ll, rr;
			public:
				JacMatrix(){}

				JacMatrix(size_t n, const mesh_t& mesh){
					size_t n_ele = mesh.n_ele();
					if(n == 3){
						cent.resize(n_ele, Eigen::MatrixXd::Zero(5,5)),
						left.resize(n_ele, Eigen::MatrixXd::Zero(5,5)),
						right.resize(n_ele, Eigen::MatrixXd::Zero(5,5));
					}else if(n == 5){
						cent.resize(n_ele, Eigen::MatrixXd::Zero(5,5)),
						left.resize(n_ele, Eigen::MatrixXd::Zero(5,5)),
						right.resize(n_ele, Eigen::MatrixXd::Zero(5,5)),
						ll.resize(n_ele, Eigen::MatrixXd::Zero(5,5)),
						rr.resize(n_ele, Eigen::MatrixXd::Zero(5,5));
					}
				}		

				void resize(size_t n, const mesh_t& mesh){
					size_t n_ele = mesh.n_ele();
					if(n == 3){
						cent.resize(n_ele, Eigen::MatrixXd::Zero(5,5)),
						left.resize(n_ele, Eigen::MatrixXd::Zero(5,5)),
						right.resize(n_ele, Eigen::MatrixXd::Zero(5,5));
					}else if(n == 5){
						cent.resize(n_ele, Eigen::MatrixXd::Zero(5,5)),
						left.resize(n_ele, Eigen::MatrixXd::Zero(5,5)),
						right.resize(n_ele, Eigen::MatrixXd::Zero(5,5)),
						ll.resize(n_ele, Eigen::MatrixXd::Zero(5,5)),
						rr.resize(n_ele, Eigen::MatrixXd::Zero(5,5));
					}
				}
		};

		void read_config(const ptr_mesh_t& p_mesh, const ConfigFile& cf){
		  	Base::read_config(p_mesh, cf);
			perturbation_ratio = (double)cf.Value( "SingleGridMethod", "Perturbation_ratio" );
    		lambda = (double)cf.Value( "SingleGridMethod", "RegularizeJacobiMatrixScaling" );
    		mu = (double)cf.Value( "SingleGridMethod", "UpdatingDisStepScaling");
    		s1 = (int)cf.Value("SingleGridMethod", "PreSmoothingSteps");
    		s2 = (int)cf.Value("SingleGridMethod", "PostSmoothingSteps");
			coarsestGrid = (int)cf.Value("SingleGridMethod", "CoarsestGrid");
			maxMGSteps = (int)cf.Value("SingleGridMethod", "MaxMGsteps");
		}

	public:
		void solve(const sol_t& sol,
		  		 Euler_sol_vec_t& Euler_sol_vec, 
		  		 const Euler_sol_vec_t& closure_vec, 
		  		 const Euler_sol_vec_t& rhs_vec,
		  		 const mesh_t& mesh){
		    size_t n_ele = mesh.n_ele();
		    Base::M = sol[0].GetOrder();
		    int steps = 0; double norm = 1;
		    //设置Newton迭代终止条件，最大步数，或解的变化足够小
			 //测试牛顿迭代所需的时间
			 //double matTime = 0;
			 //auto start = std::chrono::high_resolution_clock::now();
		    while((steps < Base::maxSteps) && (norm > 1e-12)){
				//这里注意要将delU置为0
				Euler_sol_vec_t delU(n_ele, Eigen::VectorXd::Zero(5));
				Euler_sol_vec_t res_vec = rhs_vec;
				Base::getResidual(sol, Euler_sol_vec, closure_vec, mesh, res_vec);
				#ifdef RECONSTRUCT
				JacMatrix jacMat(5, mesh);
				getJacMatSecDis(sol, Euler_sol_vec, closure_vec, res_vec, rhs_vec, mesh, jacMat);
				solveLinEqsSecDis(jacMat, delU, res_vec, mesh);
				#else
				JacMatrix jacMat(3, mesh);
				getJacMatFirDis(sol, Euler_sol_vec, closure_vec, res_vec, rhs_vec, mesh, jacMat);
				solveLinEqsFirDis(jacMat, delU, res_vec, mesh);
				#endif
				update(delU, mesh, Euler_sol_vec);
				steps++;
		    }
			 //auto end = std::chrono::high_resolution_clock::now();
			 //matTime = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			 //std::cout << "Elapsed time: " << matTime << " ms\n";
		}

		void getJacMatFirDis(const sol_t& sol,
		  		 const Euler_sol_vec_t& Euler_sol_vec, 
		  		 const Euler_sol_vec_t& closure_vec, 
		  		 const Euler_sol_vec_t& res_vec,
				 const Euler_sol_vec_t& rhs_vec,
		  		 const mesh_t& mesh,
				 JacMatrix& jacMat){
			size_t n_ele = mesh.n_ele();
			Euler_sol_vec_t temp = Euler_sol_vec;
		   for(size_t i = 0; i < n_ele; i++){
				//传入当前考虑的单元位置和需要扰动的自变量的位置
				getMatrix(sol, temp, closure_vec, rhs_vec, mesh, i, i, res_vec[i], jacMat.cent[i]);
				if(i == 0){
				 	getMatrix(sol, temp, closure_vec, rhs_vec, mesh, i, i+1, res_vec[i], jacMat.right[i]);
				}else if(i == n_ele-1){
				 	getMatrix(sol, temp, closure_vec, rhs_vec, mesh, i, i-1, res_vec[i], jacMat.left[i]);
				}else{
				 	getMatrix(sol, temp, closure_vec, rhs_vec, mesh, i, i-1, res_vec[i], jacMat.left[i]);
				 	getMatrix(sol, temp, closure_vec, rhs_vec, mesh, i, i+1, res_vec[i], jacMat.right[i]);
				}
			}
		}

		//==2025.07.13(星期日)==
		//不应该在局部对一个大的向量申请和释放空间
		//Euler_sol_vec_t temp = Euler_sol_vec;
		void getMatrix(const sol_t& sol,
				 Euler_sol_vec_t& temp,
		  		 const Euler_sol_vec_t& closure_vec, 
				 const Euler_sol_vec_t& rhs_vec,
		  		 const mesh_t& mesh,
		  		 const int nx,
		  		 const int ny,
		  		 const Euler_sol_t& res,
		  		 Eigen::MatrixXd& jacMatrix){
			for(size_t i = 0; i < 5; i++){
				double turb = perturbation_ratio * temp[ny][i] + 1e-8;
		  		//对编号ny的单元上第i个自变量做扰动
		  		temp[ny][i] += turb;
		  		//为扰动后的残量申请空间，要在外层多重网格方法中作为光滑子，这里
				//应该初始化为右端项
				Euler_sol_t resTurb = rhs_vec[nx];
		  		//求扰动后的残量
		  		Base::getResidual(sol, temp, closure_vec, mesh, nx, resTurb);
		  		//生成雅可比矩阵的第i列
		  		for(size_t j = 0; j < 5; j++){
					jacMatrix(j,i) = -(resTurb[j] - res[j]) / turb;
		  		   if(nx == ny && i == j){
						//使用局部残量范数和Knudsen数对矩阵做正则化
						jacMatrix(j, i) += lambda * (res.norm() + PROBLEM::Kn);
		  		   }
		  		}
		  		//把扰动去掉
		  		temp[ny][i] -= turb;
		  	}
		}

		void solveLinEqsFirDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			double res = 100;
			int steps = 0;
			//while(res > 1e-8 && steps < 8){
			//	mgOneStepFirDis(jacMat, sol_vec, rhs_vec, mesh);
			//	getResNormFirDis(jacMat, sol_vec, rhs_vec, res);
			//	steps++;
			//}
			for(size_t i = 0; i < maxMGSteps; i++){
				mgOneStepFirDis(jacMat, sol_vec, rhs_vec, mesh);
			}
		}

		void update(const Euler_sol_vec_t& delU,
		  	  const mesh_t& mesh,
		  	  Euler_sol_vec_t& Euler_sol_vec){
		    size_t n_ele = mesh.n_ele();
		    for(size_t i = 0; i < n_ele; i++){
				 Euler_sol_vec[i] += mu * delU[i];
		    }
		}		

		void symmetrySweepingFirDis(const JacMatrix& jacMat,
				 Euler_sol_vec_t& delU,
		  		 const Euler_sol_vec_t& rhs_vec,
		  		 const mesh_t& mesh, 
				 const size_t steps){
		    size_t n_ele = mesh.n_ele();
		    for(size_t l = 0; l < steps; l++){
				for(size_t i = 0; i < n_ele; i++){
					size_t nx = i;
					sweepingFirDis(jacMat, rhs_vec, mesh, delU, nx);
				}
				for(size_t i = n_ele; i > 0; i--){
					size_t nx = i-1;
					sweepingFirDis(jacMat, rhs_vec, mesh, delU, nx);
				}
			}
		}

		void sweepingFirDis(const JacMatrix& jacMat,
		  		 const Euler_sol_vec_t& rhs_vec,
		  		 const mesh_t& mesh, 
		  		 Euler_sol_vec_t& delU,
		  		 const size_t nx){
		    size_t n_ele = mesh.n_ele();
		    Euler_sol_t right(5);
			right = rhs_vec[nx];
		    if(nx == 0){
		  	  	right += - jacMat.right[nx] * delU[nx+1];
			}else if( nx == n_ele -1){
		  	  	right += - jacMat.left[nx] * delU[nx-1];
			}else{
		  	  	right += - jacMat.left[nx] * delU[nx-1] 
						 - jacMat.right[nx] * delU[nx+1];
		    }
		    delU[nx] = jacMat.cent[nx].inverse() * right;
		}

		void preSmoothFirDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			symmetrySweepingFirDis(jacMat, sol_vec, rhs_vec, mesh, s1);
		}
		
		void postSmoothFirDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec, 
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			symmetrySweepingFirDis(jacMat, sol_vec, rhs_vec, mesh, s2);
		}

	  	void coarsestSolveFirDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			//最粗网格上矩阵规模小，精确求解代价也不大
			exactSolveFirDis(jacMat, sol_vec, rhs_vec, mesh);
		}

		void exactSolveFirDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			size_t n_ele = mesh.n_ele();	
			mat_t fullJacMat = Eigen::MatrixXd::Zero(5*n_ele, 5*n_ele);
			Eigen::VectorXd fullRHS = Eigen::VectorXd::Zero(5*n_ele);
			getFullMatFir(fullJacMat, jacMat, mesh);
			getRHS(fullRHS, rhs_vec, mesh);
			Eigen::VectorXd fullSOL = Eigen::VectorXd::Zero(5*n_ele);
			fullSOL = fullJacMat.inverse() * fullRHS;
			for(size_t i = 0; i < n_ele; i++){
				sol_vec[i] = fullSOL.segment(i*5,5);
			}
		}

		void getFullMatFir(mat_t& fullJacMat,
				const JacMatrix& jacMat,
				const mesh_t& mesh){
			size_t n_ele = mesh.n_ele();
			for(size_t i = 0; i < n_ele; i++){
				if(i == 0){
					fullJacMat.block<5,5>(i,i) = jacMat.cent[i];
					fullJacMat.block<5,5>(i,(i+1)*5) = jacMat.right[i];
				}else if(i == n_ele-1){
					fullJacMat.block<5,5>(i*5,(i-1)*5) = jacMat.left[i];
					fullJacMat.block<5,5>(i*5,i*5) = jacMat.cent[i];
				}else{
					fullJacMat.block<5,5>(i*5,(i-1)*5) = jacMat.left[i];
					fullJacMat.block<5,5>(i*5,i*5) = jacMat.cent[i];
					fullJacMat.block<5,5>(i*5,(i+1)*5) = jacMat.right[i];
				}
			}
		}

		void getRHS(Eigen::VectorXd& fullRHS,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			size_t n_ele = mesh.n_ele();
			for(size_t i = 0; i < n_ele; i++){
				fullRHS.segment(i*5,5) = rhs_vec[i];
			}
		}

		void restriction(const Euler_sol_vec_t& sol_vec,
				Euler_sol_vec_t& coarser_sol_vec){
			size_t n_ele = sol_vec.size();
			size_t n_ele_coarser = n_ele / 2;
			coarser_sol_vec.resize(n_ele_coarser);
			for(size_t i = 0; i < n_ele_coarser; i++){
				coarser_sol_vec[i] = (sol_vec[2*i] + sol_vec[2*i+1]);
			}
		}

		void restrictFirDis(const JacMatrix& Mat,
				JacMatrix& coarMat,
				const mesh_t& mesh){
			size_t n_ele = mesh.n_ele();
			coarMat.resize(3, mesh);
			//按粗网格进行循环，对应关系要写清楚
			for(size_t i = 0; i < n_ele; i++){
				coarMat.cent[i] = (Mat.left[2*i+1] + Mat.right[2*i] + Mat.cent[2*i] + Mat.cent[2*i+1]);
				if(i == 0){
					coarMat.right[i] = Mat.right[2*i+1];
				}else if(i == n_ele-1){
					coarMat.left[i] = Mat.left[2*i];
				}else{
					coarMat.left[i] = Mat.left[2*i];
					coarMat.right[i] = Mat.right[2*i+1];
				}
			}
		}

		void backSubstitution(Euler_sol_vec_t& correction,
				Euler_sol_vec_t& sol_vec){
			size_t n_ele = correction.size();
			for(size_t i = 0; i < n_ele; i++){
				sol_vec[2*i] +=  correction[i];
				sol_vec[2*i+1] +=  correction[i];
			}
		}
	
		void getResFirDis(const JacMatrix& jacMat,
				const Euler_sol_vec_t& sol_vec, 
				const Euler_sol_vec_t& rhs_vec,
				Euler_sol_vec_t& res_vec){
			size_t n_ele = sol_vec.size();
			res_vec.resize(n_ele);
			Euler_sol_t right;
			for(size_t i = 0; i < n_ele; i++){
				right = rhs_vec[i];
				if(i == 0){
					res_vec[i] = -(jacMat.cent[i] * sol_vec[i]
					+ jacMat.right[i] * sol_vec[i+1]) + right ;
				}else if(i == n_ele-1){
					res_vec[i] = -(jacMat.left[i] * sol_vec[i-1] 
					+ jacMat.cent[i] * sol_vec[i]) + right;
				}else{
					res_vec[i] = -(jacMat.left[i] * sol_vec[i-1]
					+ jacMat.cent[i] * sol_vec[i] 
					+ jacMat.right[i] * sol_vec[i+1]) + right;
				}
			}		
		}

		void getResNormFirDis(const JacMatrix& jacMat,
				const Euler_sol_vec_t& sol_vec, 
				const Euler_sol_vec_t& rhs_vec,
				double& res){
			size_t n_ele = sol_vec.size();
			Euler_sol_vec_t res_vec;
			getResFirDis(jacMat, sol_vec, rhs_vec, res_vec);
			res = 0;
			for(size_t i = 0; i < n_ele; i++){
				res += res_vec[i].norm();
			}
			res = res / n_ele;
		}

		void subGridIteFirDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			mgOneStepFirDis(jacMat, sol_vec, rhs_vec, mesh);
		}

		void mgOneStepFirDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			size_t n_ele = sol_vec.size();
			size_t n_ele_coarser = n_ele / 2;

			//if it is the coarsest grid
			if(n_ele == coarsestGrid){
				coarsestSolveFirDis(jacMat, sol_vec, rhs_vec, mesh);
				return;
			}

			preSmoothFirDis(jacMat, sol_vec, rhs_vec, mesh);

			{
				mesh_t coarser_mesh;
				coarser_mesh.CoarseFrom(mesh);
				//粗网格上的解向量
				Euler_sol_vec_t coarser_sol_vec(n_ele_coarser, Eigen::VectorXd::Zero(5)); 

				//calculate the right hand side of the coarser grid problem 
				Euler_sol_vec_t res_vec; 
				getResFirDis(jacMat, sol_vec, rhs_vec, res_vec);
				Euler_sol_vec_t coarser_rhs_vec;
				restriction(res_vec, coarser_rhs_vec);

				JacMatrix coarMat;
				restrictFirDis(jacMat, coarMat, coarser_mesh);

				//recursively call NMG iteration
				subGridIteFirDis(coarMat, coarser_sol_vec, coarser_rhs_vec, coarser_mesh);

				//update fine grid solution
				backSubstitution(coarser_sol_vec, sol_vec);
			}

			postSmoothFirDis(jacMat, sol_vec, rhs_vec, mesh);
		}

		//二阶离散	
		void getJacMatSecDis(const sol_t& sol,
		  		 const Euler_sol_vec_t& Euler_sol_vec, 
		  		 const Euler_sol_vec_t& closure_vec, 
		  		 const Euler_sol_vec_t& res_vec,
				 const Euler_sol_vec_t& rhs_vec,
		  		 const mesh_t& mesh,
				 JacMatrix& jacMat){
			size_t n_ele = mesh.n_ele();
			getJacMatFirDis(sol, Euler_sol_vec, closure_vec, res_vec, rhs_vec, mesh, jacMat);
			Euler_sol_vec_t temp = Euler_sol_vec;
		   for(size_t i = 0; i < n_ele; i++){
			  if(i == 0 || i == 1){
			  		getMatrix(sol, temp, closure_vec, rhs_vec, mesh, i, i+2, res_vec[i], jacMat.rr[i]);
			  }else if(i == n_ele-1 || i == n_ele-2){
			      getMatrix(sol, temp, closure_vec, rhs_vec, mesh, i, i-2, res_vec[i], jacMat.ll[i]);
			  }else{
			  		getMatrix(sol, temp, closure_vec, rhs_vec, mesh, i, i-2, res_vec[i], jacMat.ll[i]);
			  		getMatrix(sol, temp, closure_vec, rhs_vec, mesh, i, i+2, res_vec[i], jacMat.rr[i]);
			  }
		   }
		}

		void solveLinEqsSecDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			for(size_t i = 0; i < maxMGSteps; i++){
				mgOneStepSecDis(jacMat, sol_vec, rhs_vec, mesh);
			}
		}

		void symmetrySweepingSecDis(const JacMatrix& jacMat,
				 Euler_sol_vec_t& delU,
		  		 const Euler_sol_vec_t& rhs_vec,
		  		 const mesh_t& mesh, 
				 const size_t steps){
		    size_t n_ele = mesh.n_ele();
		    for(size_t l = 0; l < steps; l++){
				for(size_t i = 0; i < n_ele; i++){
					size_t nx = i;
					sweepingSecDis(jacMat, rhs_vec, mesh, delU, nx);
				}
				for(size_t i = n_ele; i > 0; i--){
					size_t nx = i-1;
					sweepingSecDis(jacMat, rhs_vec, mesh, delU, nx);
				}
			}
		}

		void sweepingSecDis(const JacMatrix& jacMat,
		  		 const Euler_sol_vec_t& rhs_vec,
		  		 const mesh_t& mesh, 
		  		 Euler_sol_vec_t& delU,
		  		 const size_t nx){
		    size_t n_ele = mesh.n_ele();
		    Euler_sol_t right(5);
			right = rhs_vec[nx];
		    if(nx == 0){
		  	  	right += - jacMat.right[nx] * delU[nx+1]
						 - jacMat.rr[nx] * delU[nx+2];
			}else if(nx == 1){
				right += - jacMat.left[nx] * delU[nx-1]
						 - jacMat.right[nx] * delU[nx+1]
						 - jacMat.rr[nx] * delU[nx+2];
			}else if(nx == n_ele-2){
		  	  	right += - jacMat.ll[nx] * delU[nx-2]
						 - jacMat.left[nx] * delU[nx-1]
						 - jacMat.right[nx] * delU[nx+1];
			}else if(nx == n_ele-1){
		  	  	right += - jacMat.ll[nx] * delU[nx-2] 
						 - jacMat.left[nx] * delU[nx-1];
		    }else{
				right += - jacMat.ll[nx] * delU[nx-2] 
						 - jacMat.left[nx] * delU[nx-1]
						 - jacMat.right[nx] * delU[nx+1]
						 - jacMat.rr[nx] * delU[nx+2];
			}
		    delU[nx] = jacMat.cent[nx].inverse() * right;
		}

		void preSmoothSecDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			symmetrySweepingSecDis(jacMat, sol_vec, rhs_vec, mesh, s1);
		}
		
		void postSmoothSecDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec, 
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			symmetrySweepingSecDis(jacMat, sol_vec, rhs_vec, mesh, s2);
		}

		void coarsestSolveSecDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec, 
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			exactSolveSecDis(jacMat, sol_vec, rhs_vec, mesh);
		}

		void exactSolveSecDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			size_t n_ele = mesh.n_ele();	
			mat_t fullJacMat = Eigen::MatrixXd::Zero(5*n_ele, 5*n_ele);
			Eigen::VectorXd fullRHS = Eigen::VectorXd::Zero(5*n_ele);
			getFullMatSec(fullJacMat, jacMat, mesh);
			getRHS(fullRHS, rhs_vec, mesh);
			Eigen::VectorXd fullSOL = Eigen::VectorXd::Zero(5*n_ele);
			fullSOL = fullJacMat.inverse() * fullRHS;
			for(size_t i = 0; i < n_ele; i++){
				sol_vec[i] = fullSOL.segment(i*5,5);
			}
		}

		void getFullMatSec(mat_t& fullJacMat,
				const JacMatrix& jacMat,
				const mesh_t& mesh){
			size_t n_ele = mesh.n_ele();
			getFullMatFir(fullJacMat, jacMat, mesh);
			for(size_t i = 0; i < n_ele; i++){
				if(i == 0 || i == 1){
					fullJacMat.block<5,5>(i,(i+2)*5) = jacMat.rr[i];
				}else if(i == n_ele-2 || i == n_ele-1){
					fullJacMat.block<5,5>(i*5,(i-2)*5) = jacMat.ll[i];
				}else{
					fullJacMat.block<5,5>(i*5,(i-2)*5) = jacMat.ll[i];
					fullJacMat.block<5,5>(i*5,(i+2)*5) = jacMat.rr[i];
				}
			}
		}

		void restrictSecDis(const JacMatrix& Mat,
				JacMatrix& coarMat,
				const mesh_t& mesh){
			size_t n_ele = mesh.n_ele();
			coarMat.resize(3, mesh);
			//按粗网格进行循环，对应关系要写清楚
			for(size_t i = 0; i < n_ele; i++){
				coarMat.cent[i] = (Mat.left[2*i+1] + Mat.right[2*i] + Mat.cent[2*i] + Mat.cent[2*i+1]);
				if(i == 0){
					coarMat.right[i] = (Mat.rr[2*i] + Mat.right[2*i+1] + Mat.rr[2*i+1]) ;
				}else if(i == n_ele-1){
					coarMat.left[i] = (Mat.ll[2*i] + Mat.left[2*i] + Mat.ll[2*i+1]);
				}else{
					coarMat.left[i] = (Mat.ll[2*i] + Mat.left[2*i] + Mat.ll[2*i+1]);
					coarMat.right[i] = (Mat.rr[2*i] + Mat.right[2*i+1] + Mat.rr[2*i+1]);
				}
			}
		}
	
		void getResSecDis(const JacMatrix& jacMat,
				const Euler_sol_vec_t& sol_vec, 
				const Euler_sol_vec_t& rhs_vec,
				Euler_sol_vec_t& res_vec){
			size_t n_ele = sol_vec.size();
			res_vec.resize(n_ele);
			Euler_sol_t right;
			for(size_t i = 0; i < n_ele; i++){
				right = rhs_vec[i];
				if(i == 0){
					res_vec[i] = -(jacMat.cent[i] * sol_vec[i]
					+ jacMat.right[i] * sol_vec[i+1]
					+ jacMat.rr[i] * sol_vec[i+2]) + right ;
				}else if(i == 1){
					res_vec[i] = -(jacMat.left[i] * sol_vec[i-1]
					+ jacMat.cent[i] * sol_vec[i]
					+ jacMat.right[i] * sol_vec[i+1]
					+ jacMat.rr[i] * sol_vec[i+2]) + right ;
				}else if(i == n_ele-2){
					res_vec[i] = -(jacMat.ll[i] * sol_vec[i-2]
					+ jacMat.left[i] * sol_vec[i-1]
					+ jacMat.cent[i] * sol_vec[i]
					+ jacMat.right[i] * sol_vec[i+1]) + right ;
				}
				else if(i == n_ele-1){
					res_vec[i] = -(jacMat.ll[i] * sol_vec[i-2]
					+ jacMat.left[i] * sol_vec[i-1]
					+ jacMat.cent[i] * sol_vec[i]) + right;
				}else{
					res_vec[i] = -(jacMat.ll[i] * sol_vec[i-2]
					+ jacMat.left[i] * sol_vec[i-1]
					+ jacMat.cent[i] * sol_vec[i] 
					+ jacMat.right[i] * sol_vec[i+1]
					+ jacMat.rr[i] * sol_vec[i+2]) + right;
				}
			}		
		}

		void subGridIteSecDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			mgOneStepFirDis(jacMat, sol_vec, rhs_vec, mesh);
		}

		void mgOneStepSecDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			size_t n_ele = sol_vec.size();
			size_t n_ele_coarser = n_ele / 2;

			//if it is the coarsest grid
			if(n_ele == coarsestGrid){
				coarsestSolveSecDis(jacMat, sol_vec, rhs_vec, mesh);
				return;
			}

			preSmoothSecDis(jacMat, sol_vec, rhs_vec, mesh);

			{
				mesh_t coarser_mesh;
				coarser_mesh.CoarseFrom(mesh);
				//粗网格上的解向量
				Euler_sol_vec_t coarser_sol_vec(n_ele_coarser, Eigen::VectorXd::Zero(5)); 

				//calculate the right hand side of the coarser grid problem 
				Euler_sol_vec_t res_vec; 
				getResSecDis(jacMat, sol_vec, rhs_vec, res_vec);
				Euler_sol_vec_t coarser_rhs_vec;
				restriction(res_vec, coarser_rhs_vec);

				JacMatrix coarMat;
				restrictSecDis(jacMat, coarMat, coarser_mesh);

				//recursively call NMG iteration
				subGridIteSecDis(coarMat, coarser_sol_vec, coarser_rhs_vec, coarser_mesh);

				//update fine grid solution
				backSubstitution(coarser_sol_vec, sol_vec);
			}

			postSmoothSecDis(jacMat, sol_vec, rhs_vec, mesh);
		}
};

# endif
