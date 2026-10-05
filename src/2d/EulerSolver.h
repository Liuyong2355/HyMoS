/**
 * @file EulerSolver.h
 * @brief 二维Euler方程求解器类
 * @author Guanghan Li NUAA
 * @version 1
 * @date 2022-10-30
 */

# ifndef __EULERSOLVER2D__H__
# define __EULERSOLVER2D__H__

# include<time.h>
# include <sys/time.h>
# include <vector>
# include <cmath>
# include <limits>
# include <NRxx/Transform.h>
# include "config.h"
# include <Eigen/Core>
# include <Eigen/Dense>
# include <chrono>
# include <fstream>

template <typename RESIDUAL>
class EulerSolverBase{
	public:
	 	typedef typename RESIDUAL::dis_t dis_t;
	  	typedef typename RESIDUAL::sol_t sol_t;
	  	typedef typename sol_t::mesh_t mesh_t;
	   typedef typename dis_t::Velocity velocity_t;
	   typedef typename dis_t::Indices index_t;
	   typedef typename dis_t::value_type value_t;
	   typedef typename RESIDUAL::Problem Problem;
	   typedef typename RESIDUAL::Euler_sol_t Euler_sol_t;
	   typedef typename RESIDUAL::Euler_sol_vec_t Euler_sol_vec_t;
	public:
		double dt, CFL;
		int M;
		double tol, res;
		int maxSteps;
	public:
		EulerSolverBase() : dt(1000.0), CFL(0.49), tol(1e-8), maxSteps(1) {}
		~EulerSolverBase() {};
	public:
	   void read_config(const ConfigFile& cf){
		   CFL = (double)cf.Value( "Time", "CFL_HEs" );
		   tol = (value_t)cf.Value("Error", "Tol" );
		   maxSteps = (int)cf.Value("SingleGridMethod", "EulerEqsSolvingSteps");
		}

		void dump_data( std::ostream& os ) const {
			os.write((const char *)&CFL, sizeof(double));
			os.write((const char *)&tol, sizeof(value_t));
			os.write((const char *)&maxSteps, sizeof(int));
		}

		void load_data(std::istream& is){
		  	is.read((char *)&CFL, sizeof(double));
		  	is.read((char *)&tol, sizeof(value_t));
		  	is.read((char *)&maxSteps, sizeof(int));
		}

		void print_config(){
			std::cout << "CFL number (HEs): " << CFL << std::endl;
			std::cout << "tol (HEs): " << tol << std::endl;
			std::cout << "gamma_2: " << maxSteps << std::endl;
		}

		void timeStep(const Euler_sol_vec_t& Euler_sol_vec,
						  const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
			double dt_new = std::numeric_limits<double>::max();
			for(size_t ny = 0; ny < Ny; ny++){
			    for(size_t nx = 0; nx < Nx; nx++){
			   	   double dt_local = timeStep(Euler_sol_vec[nx][ny],
			   	   		  mesh, nx, ny);
			   	   dt_new = (dt_local < dt_new) ? dt_local : dt_new;
			    }
			}
			dt = dt_new;
		}

		double timeStep(const Euler_sol_t& Euler_sol,
				const mesh_t& mesh,
				size_t nx, size_t ny){
			double u1 = Euler_sol[1]/Euler_sol[0];
		   double u2 = Euler_sol[2]/Euler_sol[0];
		   double theta = getTemperature(Euler_sol);
			double C_M = Transform<double>::MaxRootOfHermitePolynomial( M+1 );
		   return CFL / ( (fabs(u1) + C_M * sqrt(theta)) / mesh.dx(nx, ny) +
						(fabs(u2) + C_M * sqrt(theta)) / mesh.dy(nx, ny) );
		}

		void getResidual(const sol_t& sol, 
				const Euler_sol_vec_t& Euler_sol_vec, 
				const Euler_sol_vec_t& closure_vec,
				Euler_sol_vec_t& res_vec,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
			#pragma omp parallel for
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					getResidual(sol, Euler_sol_vec, closure_vec, res_vec[nx][ny], mesh, nx, ny);
				}
			}
		}

		void getResidual(const sol_t& sol, 
				  const Euler_sol_vec_t& Euler_sol_vec, 
				  const Euler_sol_vec_t& closure_vec,
				  Euler_sol_t& res,
				  const mesh_t& mesh, 
				  size_t nx, size_t ny){
			#ifdef RECONSTRUCT
			secondOrderDisPri(sol, Euler_sol_vec, closure_vec, res, mesh, nx, ny);
			//secondOrderDisCon(sol, Euler_sol_vec, closure_vec, res, mesh, nx, ny);
			#else
			firstOrderDis(sol, Euler_sol_vec, closure_vec, res, mesh, nx, ny);
     		#endif
		}

		void firstOrderDis(const sol_t& sol, 
				  const Euler_sol_vec_t& Euler_sol_vec, 
				  const Euler_sol_vec_t& closure_vec,
				  Euler_sol_t& res,
				  const mesh_t& mesh, 
				  size_t nx, size_t ny){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
			//下面这几个变量应该定义为局部变量，而非之前那样，作为类的私有成员
			//作为类成员，会造成在并行时的竞争关系，导致算法发散
			Euler_sol_t flux_l, flux_r, flux_b, flux_t,
			   		 Euler_sol_l, Euler_sol_m, Euler_sol_r, 
					 closure_l, closure_m, closure_r, temp;
			if(nx == 0){
			    leftBoundaryCondition(sol[nx][ny], Euler_sol_vec[nx][ny], Euler_sol_l,
			   							  closure_l);
			    Euler_sol_m = Euler_sol_vec[nx][ny], closure_m = closure_vec[nx][ny];
			    Euler_sol_r = Euler_sol_vec[nx+1][ny], closure_r = closure_vec[nx+1][ny];
			}else if(nx == (Nx -1)){
			    rightBoundaryCondition(sol[nx][ny], Euler_sol_vec[nx][ny], Euler_sol_r,
			   								closure_r);
			    Euler_sol_m = Euler_sol_vec[nx][ny], closure_m = closure_vec[nx][ny];
			    Euler_sol_l = Euler_sol_vec[nx-1][ny], closure_l = closure_vec[nx-1][ny];
			}else{
			    Euler_sol_l = Euler_sol_vec[nx-1][ny], closure_l = closure_vec[nx-1][ny];
			    Euler_sol_m = Euler_sol_vec[nx][ny], closure_m = closure_vec[nx][ny];
			    Euler_sol_r = Euler_sol_vec[nx+1][ny], closure_r = closure_vec[nx+1][ny];
			}

			hllFluxHor(Euler_sol_l, Euler_sol_m, closure_l, closure_m, flux_l);
			hllFluxHor(Euler_sol_m, Euler_sol_r, closure_m, closure_r, flux_r);

			Euler_sol_t Euler_sol_b, Euler_sol_top, closure_b, closure_t;
			if(ny == 0){
			    bottomBoundaryCondition(sol[nx][ny], Euler_sol_vec[nx][ny], Euler_sol_b,
			   								 closure_b);
			    Euler_sol_m = Euler_sol_vec[nx][ny], closure_m = closure_vec[nx][ny];
			    Euler_sol_top = Euler_sol_vec[nx][ny+1], closure_t = closure_vec[nx][ny+1];
			}else if(ny == (Ny - 1)){
			    topBoundaryCondition(sol[nx][ny], Euler_sol_vec[nx][ny], Euler_sol_top,
			   								 closure_t);
			    Euler_sol_m = Euler_sol_vec[nx][ny], closure_m = closure_vec[nx][ny];
			    Euler_sol_b = Euler_sol_vec[nx][ny-1], closure_b = closure_vec[nx][ny-1];
			}else{
			    Euler_sol_b = Euler_sol_vec[nx][ny-1], closure_b = closure_vec[nx][ny-1];
			    Euler_sol_m = Euler_sol_vec[nx][ny], closure_m = closure_vec[nx][ny];
			    Euler_sol_top = Euler_sol_vec[nx][ny+1], closure_t = closure_vec[nx][ny+1];
			}

			hllFluxVer(Euler_sol_b, Euler_sol_m, closure_b, closure_m, flux_b);
			hllFluxVer(Euler_sol_m, Euler_sol_top, closure_m, closure_t, flux_t);

			res +=  -((flux_r - flux_l) / mesh.dx(nx,ny) 
					+ (flux_t - flux_b) / mesh.dy(nx,ny));
		}

		//对原始变量做重构
		void secondOrderDisPri(const sol_t& sol,
				const Euler_sol_vec_t& Euler_sol_vec, 
				const Euler_sol_vec_t& closure_vec,
				Euler_sol_t& res,
				const mesh_t& mesh, 
				size_t nx, size_t ny){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
			//存矩形单元四条边上计算出来的数值通量
			Euler_sol_t flux_l, flux_r, flux_b, flux_t;
			//存守恒变量在边界处的极限值
			Euler_sol_t re_lr, re_ml, re_mr, re_rl, dummy;
			//存封闭项在边界处的极限值
			Euler_sol_t closure_lr, closure_ml, closure_mr, closure_rl;
			if(nx == 0){
				reconstructPrimiHor(Euler_sol_vec, mesh, nx, ny, re_ml, re_mr);
				reconstructHor(closure_vec, mesh, nx, ny, closure_ml, closure_mr);
				
				reconstructPrimiHor(Euler_sol_vec, mesh, nx+1, ny, re_rl, dummy);
				reconstructHor(closure_vec, mesh, nx+1, ny, closure_rl, dummy);

				dis_t dis_ml;
			   dis_ml = sol[nx][ny]; dis_ml *= value_t(1.5); dis_ml.Add(-0.5, sol[nx+1][ny]);
				leftBoundaryCondition(dis_ml, re_ml, re_lr, closure_lr);
			}else if(nx == (Nx-1)){
				reconstructPrimiHor(Euler_sol_vec, mesh, nx, ny, re_ml, re_mr);
				reconstructHor(closure_vec, mesh, nx, ny, closure_ml, closure_mr);

				reconstructPrimiHor(Euler_sol_vec, mesh, nx-1, ny, dummy, re_lr);
				reconstructHor(closure_vec, mesh, nx-1, ny, dummy, closure_lr);
				
				dis_t dis_mr;
			   dis_mr = sol[nx][ny]; dis_mr *= value_t(1.5); dis_mr.Add(-0.5, sol[nx-1][ny]);
				rightBoundaryCondition(dis_mr, re_mr, re_rl, closure_rl);
			}else{
				reconstructPrimiHor(Euler_sol_vec, mesh, nx-1, ny, dummy, re_lr);
				reconstructHor(closure_vec, mesh, nx-1, ny, dummy, closure_lr);

				reconstructPrimiHor(Euler_sol_vec, mesh, nx, ny, re_ml, re_mr);
				reconstructHor(closure_vec, mesh, nx, ny, closure_ml, closure_mr);
				
				reconstructPrimiHor(Euler_sol_vec, mesh, nx+1, ny, re_rl, dummy);
				reconstructHor(closure_vec, mesh, nx+1, ny, closure_rl, dummy);
			}

				// The right interface flux uses the right reconstructed state.
				hllFluxHor(re_lr, re_ml, closure_lr, closure_ml, flux_l);
			hllFluxHor(re_mr, re_rl, closure_mr, closure_rl, flux_r);
			res += -(flux_r - flux_l) / mesh.dx(nx,ny);
				
			if(ny == 0){
				reconstructPrimiVer(Euler_sol_vec, mesh, nx, ny, re_ml, re_mr);
				reconstructVer(closure_vec, mesh, nx, ny, closure_ml, closure_mr);
				
				reconstructPrimiVer(Euler_sol_vec, mesh, nx, ny+1, re_rl, dummy);
				reconstructVer(closure_vec, mesh, nx, ny+1, closure_rl, dummy);
				
				dis_t dis_ml;
			   dis_ml = sol[nx][ny]; dis_ml *= value_t(1.5); dis_ml.Add(-0.5, sol[nx][ny+1]);
				bottomBoundaryCondition(dis_ml, re_ml, re_lr, closure_lr);
			}else if(ny == (Ny-1)){
				reconstructPrimiVer(Euler_sol_vec, mesh, nx, ny, re_ml, re_mr);
				reconstructVer(closure_vec, mesh, nx, ny, closure_ml, closure_mr);

				reconstructPrimiVer(Euler_sol_vec, mesh, nx, ny-1, dummy, re_lr);
				reconstructVer(closure_vec, mesh, nx, ny-1, dummy, closure_lr);

				dis_t dis_mr;
			   dis_mr = sol[nx][ny]; dis_mr *= value_t(1.5); dis_mr.Add(-0.5, sol[nx][ny-1]);
				topBoundaryCondition(dis_mr, re_mr, re_rl, closure_rl);
			}else{
				reconstructPrimiVer(Euler_sol_vec, mesh, nx, ny-1, dummy, re_lr);
				reconstructVer(closure_vec, mesh, nx, ny-1, dummy, closure_lr);
				
				reconstructPrimiVer(Euler_sol_vec, mesh, nx, ny, re_ml, re_mr);
				reconstructVer(closure_vec, mesh, nx, ny, closure_ml, closure_mr);
				
				reconstructPrimiVer(Euler_sol_vec, mesh, nx, ny+1, re_rl, dummy);
				reconstructVer(closure_vec, mesh, nx, ny+1, closure_rl, dummy);
			}
			
			hllFluxVer(re_lr, re_ml, closure_lr, closure_ml, flux_b);
			hllFluxVer(re_mr, re_rl, closure_mr, closure_rl, flux_t);
			
			res += - (flux_t - flux_b) / mesh.dy(nx,ny);
		}
		
		void reconstructPrimiHor(const Euler_sol_vec_t& Euler_sol_vec,
				const mesh_t& mesh,
				size_t nx, size_t ny,
				Euler_sol_t& leftLim,
				Euler_sol_t& rightLim){
			size_t Nx = mesh.n_ele_x();
			Euler_sol_t primitiveLeft(5), primitiveMid(5), primitiveRight(5), slope(5);
			Euler_sol_t primiLeftLim(5), primiRightLim(5);
			if(nx == 0){
				conserToPrimi(Euler_sol_vec[nx][ny], primitiveMid);
				conserToPrimi(Euler_sol_vec[nx+1][ny], primitiveRight);
				primiLeftLim = 1.5 * primitiveMid - 0.5 * primitiveRight;
				primiRightLim = 0.5 * primitiveMid + 0.5 * primitiveRight;
			}else if(nx == (Nx-1)){
				conserToPrimi(Euler_sol_vec[nx-1][ny], primitiveLeft);
				conserToPrimi(Euler_sol_vec[nx][ny], primitiveMid);
				primiLeftLim = 0.5 * primitiveMid + 0.5 * primitiveLeft;
				primiRightLim = 1.5 * primitiveMid - 0.5 * primitiveLeft;
			}else{
				double dx = mesh.dx(nx,ny);
				conserToPrimi(Euler_sol_vec[nx-1][ny], primitiveLeft);
				conserToPrimi(Euler_sol_vec[nx][ny], primitiveMid);
				conserToPrimi(Euler_sol_vec[nx+1][ny], primitiveRight);
				slope = (primitiveRight - primitiveLeft) / (2);
				primiLeftLim = primitiveMid - (0.5 ) * slope;
			   primiRightLim = primitiveMid + (0.5 ) * slope;
			}
			primiToConser(primiLeftLim, leftLim);
			primiToConser(primiRightLim, rightLim);
		}		

		void reconstructPrimiVer(const Euler_sol_vec_t& Euler_sol_vec,
				const mesh_t& mesh,
				size_t nx, size_t ny,
				Euler_sol_t& leftLim,
				Euler_sol_t& rightLim){
			size_t Ny = mesh.n_ele_y();
			Euler_sol_t primitiveLeft(5), primitiveMid(5), primitiveRight(5), slope(5);
			Euler_sol_t primiLeftLim(5), primiRightLim(5);
			if(ny == 0){
				conserToPrimi(Euler_sol_vec[nx][ny], primitiveMid);
				conserToPrimi(Euler_sol_vec[nx][ny+1], primitiveRight);
				primiLeftLim = 1.5 * primitiveMid - 0.5 * primitiveRight;
				primiRightLim = 0.5 * primitiveMid + 0.5 * primitiveRight;
			}else if(ny == (Ny-1)){
				conserToPrimi(Euler_sol_vec[nx][ny-1], primitiveLeft);
				conserToPrimi(Euler_sol_vec[nx][ny], primitiveMid);
				primiLeftLim = 0.5 * primitiveMid + 0.5 * primitiveLeft;
				primiRightLim = 1.5 * primitiveMid - 0.5 * primitiveLeft;
			}else{
				double dy = mesh.dy(nx,ny);
				conserToPrimi(Euler_sol_vec[nx][ny-1], primitiveLeft);
				conserToPrimi(Euler_sol_vec[nx][ny], primitiveMid);
				conserToPrimi(Euler_sol_vec[nx][ny+1], primitiveRight);
				slope = (primitiveRight - primitiveLeft) / (2);
				primiLeftLim = primitiveMid - (0.5 ) * slope;
			   primiRightLim = primitiveMid + (0.5 ) * slope;
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

		//对守恒变量做重构
		void secondOrderDisCon(const sol_t& sol,
				const Euler_sol_vec_t& Euler_sol_vec, 
				const Euler_sol_vec_t& closure_vec,
				Euler_sol_t& res,
				const mesh_t& mesh, 
				size_t nx, size_t ny){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
			//存矩形单元四条边上计算出来的数值通量
			Euler_sol_t flux_l, flux_r, flux_b, flux_t;
			//存守恒变量在边界处的极限值
			Euler_sol_t re_lr, re_ml, re_mr, re_rl, dummy;
			//存封闭项在边界处的极限值
			Euler_sol_t closure_lr, closure_ml, closure_mr, closure_rl;
		
			if(nx == 0){
				reconstructHor(Euler_sol_vec, mesh, nx, ny, re_ml, re_mr);
				reconstructHor(closure_vec, mesh, nx, ny, closure_ml, closure_mr);
				
				reconstructHor(Euler_sol_vec, mesh, nx+1, ny, re_rl, dummy);
				reconstructHor(closure_vec, mesh, nx+1, ny, closure_rl, dummy);

				dis_t dis_ml;
			   dis_ml = sol[nx][ny]; dis_ml *= value_t(1.5); dis_ml.Add(-0.5, sol[nx+1][ny]);
				leftBoundaryCondition(dis_ml, re_ml, re_lr, closure_lr);
			}else if(nx == (Nx-1)){
				reconstructHor(Euler_sol_vec, mesh, nx, ny, re_ml, re_mr);
				reconstructHor(closure_vec, mesh, nx, ny, closure_ml, closure_mr);

				reconstructHor(Euler_sol_vec, mesh, nx-1, ny, dummy, re_lr);
				reconstructHor(closure_vec, mesh, nx-1, ny, dummy, closure_lr);
				
				dis_t dis_mr;
			   dis_mr = sol[nx][ny]; dis_mr *= value_t(1.5); dis_mr.Add(-0.5, sol[nx-1][ny]);
				rightBoundaryCondition(dis_mr, re_mr, re_rl, closure_rl);
			}else{
				reconstructHor(Euler_sol_vec, mesh, nx-1, ny, dummy, re_lr);
				reconstructHor(closure_vec, mesh, nx-1, ny, dummy, closure_lr);

				reconstructHor(Euler_sol_vec, mesh, nx, ny, re_ml, re_mr);
				reconstructHor(closure_vec, mesh, nx, ny, closure_ml, closure_mr);
				
				reconstructHor(Euler_sol_vec, mesh, nx+1, ny, re_rl, dummy);
				reconstructHor(closure_vec, mesh, nx+1, ny, closure_rl, dummy);
			}

				// The right interface flux uses the right reconstructed state.
				hllFluxHor(re_lr, re_ml, closure_lr, closure_ml, flux_l);
			hllFluxHor(re_mr, re_rl, closure_mr, closure_rl, flux_r);
			res += -(flux_r - flux_l) / mesh.dx(nx,ny);
			
			if(ny == 0){
				reconstructVer(Euler_sol_vec, mesh, nx, ny, re_ml, re_mr);
				reconstructVer(closure_vec, mesh, nx, ny, closure_ml, closure_mr);
				
				reconstructVer(Euler_sol_vec, mesh, nx, ny+1, re_rl, dummy);
				reconstructVer(closure_vec, mesh, nx, ny+1, closure_rl, dummy);
				
				dis_t dis_ml;
			   dis_ml = sol[nx][ny]; dis_ml *= value_t(1.5); dis_ml.Add(-0.5, sol[nx][ny+1]);
				bottomBoundaryCondition(dis_ml, re_ml, re_lr, closure_lr);
			}else if(ny == (Ny-1)){
				reconstructVer(Euler_sol_vec, mesh, nx, ny, re_ml, re_mr);
				reconstructVer(closure_vec, mesh, nx, ny, closure_ml, closure_mr);

				reconstructVer(Euler_sol_vec, mesh, nx, ny-1, dummy, re_lr);
				reconstructVer(closure_vec, mesh, nx, ny-1, dummy, closure_lr);

				dis_t dis_mr;
			   dis_mr = sol[nx][ny]; dis_mr *= value_t(1.5); dis_mr.Add(-0.5, sol[nx][ny-1]);
				topBoundaryCondition(dis_mr, re_mr, re_rl, closure_rl);
			}else{
				reconstructVer(Euler_sol_vec, mesh, nx, ny-1, dummy, re_lr);
				reconstructVer(closure_vec, mesh, nx, ny-1, dummy, closure_lr);
				
				reconstructVer(Euler_sol_vec, mesh, nx, ny, re_ml, re_mr);
				reconstructVer(closure_vec, mesh, nx, ny, closure_ml, closure_mr);
				
				reconstructVer(Euler_sol_vec, mesh, nx, ny+1, re_rl, dummy);
				reconstructVer(closure_vec, mesh, nx, ny+1, closure_rl, dummy);
			}
			
			hllFluxVer(re_lr, re_ml, closure_lr, closure_ml, flux_b);
			hllFluxVer(re_mr, re_rl, closure_mr, closure_rl, flux_t);
			
			res += - (flux_t - flux_b) / mesh.dy(nx,ny);
		}
		
		void reconstructHor(const Euler_sol_vec_t& Euler_sol_vec,
				const mesh_t& mesh,
				size_t nx, size_t ny,
				Euler_sol_t& leftLim,
				Euler_sol_t& rightLim){
			size_t Nx = mesh.n_ele_x();
			Euler_sol_t slope;
			if(nx == 0){
				leftLim = 1.5 * Euler_sol_vec[nx][ny] - 0.5 * Euler_sol_vec[nx+1][ny];
			    rightLim = 0.5 * Euler_sol_vec[nx][ny] + 0.5 * Euler_sol_vec[nx+1][ny];
			}else if(nx == (Nx-1)){
				leftLim = 0.5 * Euler_sol_vec[nx][ny] + 0.5 * Euler_sol_vec[nx-1][ny];
			    rightLim = 1.5 * Euler_sol_vec[nx][ny] - 0.5 * Euler_sol_vec[nx-1][ny];
			}else{
				double dx = mesh.dx(nx,ny);
			    slope = (Euler_sol_vec[nx+1][ny] - Euler_sol_vec[nx-1][ny]) / (2 * dx );
			    leftLim = Euler_sol_vec[nx][ny] - (0.5 * dx) * slope;
			    rightLim = Euler_sol_vec[nx][ny] + (0.5 * dx) * slope;
			}
		}

		void reconstructVer(const Euler_sol_vec_t& Euler_sol_vec,
				const mesh_t& mesh,
				size_t nx, size_t ny,
				Euler_sol_t& leftLim,
				Euler_sol_t& rightLim){
			size_t Ny = mesh.n_ele_y();
			Euler_sol_t slope;
			if(ny == 0){
				leftLim = 1.5 * Euler_sol_vec[nx][ny] - 0.5 * Euler_sol_vec[nx][ny+1];
			   rightLim = 0.5 * Euler_sol_vec[nx][ny] + 0.5 * Euler_sol_vec[nx][ny+1];
			}else if(ny == (Ny-1)){
				leftLim = 0.5 * Euler_sol_vec[nx][ny] + 0.5 * Euler_sol_vec[nx][ny-1];
			   rightLim = 1.5 * Euler_sol_vec[nx][ny] - 0.5 * Euler_sol_vec[nx][ny-1];
			}else{
				double dy = mesh.dy(nx, ny);
			   slope = (Euler_sol_vec[nx][ny+1] - Euler_sol_vec[nx][ny-1]) / (2 * dy) ;
			   leftLim = Euler_sol_vec[nx][ny] - (0.5 * dy) * slope;
			   rightLim = Euler_sol_vec[nx][ny] + (0.5 * dy) * slope;
			}
		}

		void hllFluxHor(const Euler_sol_t& Euler_sol_l,
				const Euler_sol_t& Euler_sol_r,
				const Euler_sol_t& closure_l,
				const Euler_sol_t& closure_r,
				Euler_sol_t& flux){
			double S_L, S_R;
			Euler_sol_t flux_1, flux_2;
			flux_1.resize(5), flux_2.resize(5), flux.resize(5);
			getSHorizontal(S_L, S_R, Euler_sol_l, Euler_sol_r);
			fluxFunHorizontal(Euler_sol_l, closure_l, flux_1);
			fluxFunHorizontal(Euler_sol_r, closure_r, flux_2);
		    if(S_L >= 0){
				flux = flux_1;
		    }else if (S_L < 0 && 0 < S_R){
				flux = (S_R * flux_1 - S_L * flux_2 + S_L * S_R
				* (Euler_sol_r - Euler_sol_l) ) / (S_R - S_L);
		    }else if (S_R <= 0){
				flux = flux_2;
		    }
		}

		void hllFluxVer(const Euler_sol_t& Euler_sol_b,
				const Euler_sol_t& Euler_sol_m,
				const Euler_sol_t& closure_b,
				const Euler_sol_t& closure_m,
				Euler_sol_t& flux){
		    double S_B, S_T;
			Euler_sol_t flux_1, flux_2;
			flux_1.resize(5), flux_2.resize(5), flux.resize(5);
			getSVertical(S_B, S_T, Euler_sol_b, Euler_sol_m);
			fluxFunVertical(Euler_sol_b, closure_b, flux_1);
			fluxFunVertical(Euler_sol_m, closure_m, flux_2);
		    if(S_B >= 0){
				flux = flux_1;
		    }else if (S_B < 0 && 0 < S_T){
				flux = (S_T * flux_1 - S_B * flux_2 + S_B * S_T
				* (Euler_sol_m - Euler_sol_b)) / (S_T - S_B);
		    }else if (S_T <= 0){
				flux = flux_2;
		    }
		}

		void getSHorizontal(double& S_L, double& S_R,
				const Euler_sol_t& Euler_sol_l,
				const Euler_sol_t& Euler_sol_r){
		    double u_l = Euler_sol_l[1] / Euler_sol_l[0];
		    double u_r = Euler_sol_r[1] / Euler_sol_r[0];
		    double theta_l = getTemperature(Euler_sol_l);
		    double theta_r = getTemperature(Euler_sol_r);
		    double C_M = Transform<double>::MaxRootOfHermitePolynomial( M+1 );
		    S_L = std::min(u_l - C_M*sqrt(theta_l), u_r - C_M*sqrt(theta_r));
		    S_R = std::max(u_l + C_M*sqrt(theta_l), u_r + C_M*sqrt(theta_r));
		}

		void getSVertical(double& S_B, double& S_T,
				const Euler_sol_t& Euler_sol_b,
				const Euler_sol_t& Euler_sol_top){
		    double u_b = Euler_sol_b[2] / Euler_sol_b[0];
		    double u_t = Euler_sol_top[2] / Euler_sol_top[0];
		    double theta_b = getTemperature(Euler_sol_b);
		    double theta_t = getTemperature(Euler_sol_top);
		    double C_M = Transform<double>::MaxRootOfHermitePolynomial( M+1 );
		    S_B = std::min(u_b - C_M*sqrt(theta_b), u_t - C_M*sqrt(theta_t));
		    S_T = std::max(u_b + C_M*sqrt(theta_b), u_t + C_M*sqrt(theta_t));
		}
		
		double getTemperature(const Euler_sol_t& Euler_sol){
		    double rho = Euler_sol[0];
		    double u1 = Euler_sol[1] / rho;
		    double u2 = Euler_sol[2] / rho;
		    double u3 = Euler_sol[3] / rho;
		    return ( Euler_sol[4] / rho - (u1*u1 + u2*u2 + u3*u3) / 2.0 ) * (2.0 / 3.0);
		}

		void updateSolution(Euler_sol_vec_t& Euler_sol_vec,
		 	   const Euler_sol_vec_t& res_vec,
			   const mesh_t& mesh){
		    size_t Nx = mesh.n_ele_x();
		    size_t Ny = mesh.n_ele_y();
		    #pragma omp parallel for
		    for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					updateSolution(Euler_sol_vec[nx][ny], res_vec[nx][ny]);
				}
		    }
		}

		void updateSolution(Euler_sol_t& Euler_sol,
		 	   const Euler_sol_t& res){
			Euler_sol += dt * res;
		}

		void printMac(const dis_t& dis){
			double pv[5];
			dis.PrimitiveVars(pv);
			double rho = pv[0];
			double u1 = pv[1];
			double u2 = pv[2];
			double u3 = pv[3];
			double theta = pv[4];
			std::cout << "rho = " << rho << std::endl;
			std::cout << "u = " << u1 << " " << u2 << " " << u3 << std::endl;
			std::cout << "theta = " << theta << std::endl;
		}

		void printMac(const Euler_sol_t& Euler_sol){
			double rho = Euler_sol[0];
			double u1 = Euler_sol[1] / rho;
			double u2 = Euler_sol[2] / rho;
			double u3 = Euler_sol[3] / rho;
			double theta = getTemperature(Euler_sol);
			std::cout << "rho = " << rho << std::endl;
			std::cout << "u = " << u1 << " " << u2 << " " << u3 << std::endl;
			std::cout << "theta = " << theta << std::endl;
		}

		void leftBoundaryCondition(const dis_t& dis,
		 	   const Euler_sol_t& Euler_sol,
		 	   Euler_sol_t& Euler_sol_l,
		 	   Euler_sol_t& closure_l){
		    dis_t temp = dis;
		    updateMacroVar(temp, Euler_sol);
		    dis_t dis_ghost;
		    Problem::leftBoundaryCondition(temp, dis_ghost);
		    getMacroVar(dis_ghost, Euler_sol_l, closure_l);
		}

		void rightBoundaryCondition(const dis_t& dis,
		 		const Euler_sol_t& Euler_sol,
		 		Euler_sol_t& Euler_sol_r,
		 		Euler_sol_t& closure_r){
		    dis_t temp = dis;
		    updateMacroVar(temp, Euler_sol);
		    dis_t dis_ghost;
		    Problem::rightBoundaryCondition(temp, dis_ghost); 
		    getMacroVar(dis_ghost, Euler_sol_r, closure_r);
		}

		void topBoundaryCondition(const dis_t& dis, 
				  const Euler_sol_t& Euler_sol,
				  Euler_sol_t& Euler_sol_top,
				  Euler_sol_t& closure_t){
		    dis_t temp = dis;
		    updateMacroVar(temp, Euler_sol);
		    dis_t dis_ghost;
		    Problem::topBoundaryCondition(temp, dis_ghost); 
		    getMacroVar(dis_ghost, Euler_sol_top, closure_t);
		}

		void bottomBoundaryCondition(const dis_t& dis,
		   	 const Euler_sol_t& Euler_sol,
				 Euler_sol_t& Euler_sol_b,
				 Euler_sol_t& closure_b){
			dis_t temp = dis;
			updateMacroVar(temp, Euler_sol);
			dis_t dis_ghost;
			Problem::bottomBoundaryCondition(temp, dis_ghost); 
			getMacroVar(dis_ghost, Euler_sol_b, closure_b);
		}

		void getMacroVar(const dis_t& dis, 
				  Euler_sol_t& ghost, 
				  Euler_sol_t& closure){
			static index_t sigma11(2,0,0), sigma12(1,1,0), sigma13(1,0,1),
					sigma22(0,2,0), sigma23(0,1,1), sigma33(0,0,2);
			ghost.resize(5), closure.resize(10);
			//守恒量
			double pv[5];
			dis.PrimitiveVars(pv);
			ghost[0] = pv[0];
			ghost[1] = pv[0] * pv[1];
			ghost[2] = pv[0] * pv[2];
			ghost[3] = pv[0] * pv[3];
			ghost[4] = pv[0] * (1.5 * pv[4] + (pv[1]*pv[1] + pv[2]*pv[2] +
						  pv[3]*pv[3])/2.0 );

			//压力，应力张量
			closure[0] = pv[0] * pv[4];
			closure[1] = 2 * dis(sigma11);
			closure[2] = dis(sigma12);
			closure[3] = dis(sigma13);
			closure[4] = 2 * dis(sigma22);
			closure[5] = dis(sigma23);
			closure[6] = 2 * dis(sigma33);

			//热通量
			double heatflux[3];
			dis.HeatFlux(heatflux);
			closure[7] = heatflux[0];
			closure[8] = heatflux[1];
			closure[9] = heatflux[2];
		}

		//这个地方允许对传入的分布函数作出改变，是存在风险的
		//只是为了在边界条件的函数中对临时的分布函数作出改变
		//不应该在这个文件之外调用
		void updateMacroVar(dis_t& dis,
				  const Euler_sol_t& Euler_sol){
			double rho = Euler_sol[0];
			double u1 = Euler_sol[1] / rho;
			double u2 = Euler_sol[2] / rho;
			double u3 = Euler_sol[3] / rho;
			double theta = getTemperature(Euler_sol);
			//传递密度、速度、温度
			dis.front() = rho;
			typename dis_t::Velocity Center{u1, u2, u3};
			value_t Scaling = sqrt(theta);
			dis.Reinit(Center, Scaling);
		}

		void fluxFunHorizontal(const Euler_sol_t& Euler_sol,
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
					  u3*closure[3] + closure[7];
		}

		void fluxFunVertical(const Euler_sol_t& Euler_sol,
				  const Euler_sol_t& closure, 
				  Euler_sol_t& F){
			double rho = Euler_sol[0];
			double u1 = Euler_sol[1] / rho;
			double u2 = Euler_sol[2] / rho;
			double u3 = Euler_sol[3] / rho;
			double theta = getTemperature(Euler_sol);
			double p = rho * theta;
			F[0] = rho * u2;
			F[1] = rho * u2 * u1 + closure[2];
			F[2] = rho * u2 * u2 + p + closure[4];
			F[3] = rho * u2 * u3 + closure[5];
			F[4] = u2 * (Euler_sol[4] + p) + u1*closure[2] + u2*closure[4] + 
					  u3*closure[5] + closure[8];
		}

		double getL1Norm(const Euler_sol_vec_t& Euler_sol_vec,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
			double normL1 = 0;
			for(size_t ny = 0; ny < Ny; ny++){
			   for(size_t nx = 0; nx < Nx; nx++){
				   double temp = getL1Norm(Euler_sol_vec[nx][ny]);
				   normL1 += temp;
			   }
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

		double getL2Norm(const Euler_sol_vec_t& Euler_sol_vec,
		 					  const mesh_t& mesh){
		 	 size_t Nx = mesh.n_ele_x();
		 	 size_t Ny = mesh.n_ele_y();
		 	 double normL2 = 0;
		 	 for(size_t ny = 0; ny < Ny; ny++){
		 		  for(size_t nx = 0; nx < Nx; nx++){
		 			   double temp = getL2Norm(Euler_sol_vec[nx][ny]);
		 				normL2 += temp;
		 		  }
		 	 }
		 	 return sqrt(normL2);
		}

		double getL2Norm(const Euler_sol_t& Euler_sol){
		  	int n_ele = Euler_sol.size();
		  	double normL2 = 0;
		  	for(int k = 0; k < n_ele; k++){
		  		 normL2 += pow(Euler_sol[k], 2);
		  	}
		  	return normL2;
		}

		double getInfNorm(const Euler_sol_vec_t& Euler_sol_vec,
		  					  const mesh_t& mesh){
		     size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	double normInf = 0;
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		  for(size_t nx = 0; nx < Nx; nx++){
		  			   double temp = getInfNorm(Euler_sol_vec[nx][ny]);
		  				normInf = (normInf < temp) ? temp : normInf;
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
		     size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	Euler_sol_vec_t Euler_sol_vec_dif(Nx, std::vector<Euler_sol_t>(Ny));
		  	#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		 for(size_t nx = 0; nx < Nx; nx++){
		  			  int n_var = Euler_sol_vec_old[nx][ny].size();
		  			  Euler_sol_vec_dif[nx][ny].resize(n_var);
		  			  for(int k = 0; k < n_var; k++){
		  				   Euler_sol_vec_dif[nx][ny][k] = (Euler_sol_vec_old[nx][ny][k] -
		  					Euler_sol_vec_new[nx][ny][k]);
		  			  }
		  		 }
		  	}
		  	//有三种范数可供选择
		  	//return getInfNorm(Euler_sol_vec_dif, mesh);
		  	return getL1Norm(Euler_sol_vec_dif, mesh);
		  	//return getL2Norm(Euler_sol_vec_dif, mesh);
		}

		double getResidualNorm(const sol_t& sol, 
				  const Euler_sol_vec_t& Euler_sol_vec,
				  const Euler_sol_vec_t& closure_vec,
				  const Euler_sol_vec_t& rhs_vec,
				  const mesh_t& mesh){
		    size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
			Euler_sol_vec_t res_vec = rhs_vec;
		  	getResidual(sol, Euler_sol_vec, closure_vec, res_vec, mesh);
			res = 0;
		    for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					res += res_vec[nx][ny].norm();
					//res += res_vec[nx][ny].lpNorm<1>();
					//res += res_vec[nx][ny].lpNorm<Infinity>();
				}
			}
			res = res / (Nx * Ny);
			return res;
		}

		void printMat(const Eigen::MatrixXd& mat){
		    std::cout << mat << std::endl;
		}

		void printVec(const Eigen::VectorXd& vec){
		    std::cout << vec << std::endl;
		}
		
		//时间测试代码，使用时将其放置在需要测试的代码前后
		//引入头文件chrono
		//void time(){
		//	double Time = 0;
		//	auto start = std::chrono::high_resolution_clock::now();
		//	//需要测试的代码块
		//	auto end = std::chrono::high_resolution_clock::now();
		//	Time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
		//	std::cout << "Elapsed time: " << Time << " ms\n";
		//}

		virtual void solve(const sol_t& sol, 
				  Euler_sol_vec_t& Euler_sol_vec, 
				  const Euler_sol_vec_t& closure_vec,
				  const Euler_sol_vec_t& rhs_vec,
				  const mesh_t& mesh,
				  unsigned int steps = 1){}
};

template <typename RESIDUAL>
class ExplicitTimeMacEqs : public EulerSolverBase<RESIDUAL>{
	public:
	    typedef EulerSolverBase<RESIDUAL> Base;
	    typedef typename RESIDUAL::dis_t dis_t;
	    typedef typename RESIDUAL::sol_t sol_t;
	    typedef typename sol_t::mesh_t mesh_t;
	    typedef typename dis_t::Velocity velocity_t;
	    typedef typename dis_t::Indices index_t;
	    typedef typename dis_t::value_type value_t;
	    typedef typename RESIDUAL::Problem Problem;
	    typedef typename RESIDUAL::Euler_sol_t Euler_sol_t;
	    typedef typename RESIDUAL::Euler_sol_vec_t Euler_sol_vec_t;
	public:
		void solve(const sol_t& sol, 
		  	       Euler_sol_vec_t& Euler_sol_vec, 
		  			 const Euler_sol_vec_t& closure_vec,
		  			 const Euler_sol_vec_t& rhs_vec,
		  			 const mesh_t& mesh){
		  	Base::M = sol[0][0].GetOrder();
		  	int steps = 0; double dif = 1.0;
			while((steps < Base::maxSteps) && (dif > Base::tol)){
		  		Euler_sol_vec_t Euler_sol_vec_old = Euler_sol_vec;
		  		Euler_sol_vec_t res_vec = rhs_vec;
		  		Base::timeStep(Euler_sol_vec, mesh);
		  		Base::getResidual(sol, Euler_sol_vec, closure_vec, res_vec, mesh);
		  		Base::updateSolution(Euler_sol_vec, res_vec, mesh);
		  		steps++;
		  		dif = Base::getDifNorm(Euler_sol_vec_old, Euler_sol_vec, mesh);
		  	}
		}
};

template <typename RESIDUAL>
class GaussSeidelMacEqs : public EulerSolverBase<RESIDUAL>{
	public:
	    typedef EulerSolverBase<RESIDUAL> Base;
	    typedef typename RESIDUAL::dis_t dis_t;
	    typedef typename RESIDUAL::sol_t sol_t;
	    typedef typename sol_t::mesh_t mesh_t;
	    typedef typename dis_t::Velocity velocity_t;
	    typedef typename dis_t::Indices index_t;
	    typedef typename dis_t::value_type value_t;
	    typedef typename RESIDUAL::Problem Problem;
	    typedef typename RESIDUAL::Euler_sol_t Euler_sol_t;
	    typedef typename RESIDUAL::Euler_sol_vec_t Euler_sol_vec_t;
	public:
		void solve(const sol_t& sol, 
		  	  Euler_sol_vec_t& Euler_sol_vec, 
		  	  const Euler_sol_vec_t& closure_vec,
		  	  const Euler_sol_vec_t& rhs_vec,
		  	  const mesh_t& mesh){
		  	Base::M = sol[0][0].GetOrder();
		  	int steps = 0; double dif = 1.0;
		    while((steps < Base::maxSteps) && (dif > Base::tol)){
		  		Euler_sol_vec_t Euler_sol_vec_old = Euler_sol_vec;
		  		symSweeping(sol, Euler_sol_vec, closure_vec, rhs_vec, mesh);
		  		steps++;
				dif = Base::getDifNorm(Euler_sol_vec_old, Euler_sol_vec, mesh);
		  	}
		}

		void symSweeping(const sol_t& sol, 
		  	  Euler_sol_vec_t& Euler_sol_vec, 
		  	  const Euler_sol_vec_t& closure_vec,
		  	  const Euler_sol_vec_t& rhs_vec,
		  	  const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	//res_vec这个量在使用的时候，单层网格上要保证每次计算开始时是0
		  	//多重网格上计算开始时是方程右端项，这样，后面的累加操作才是正确的.
		  	Euler_sol_vec_t res_vec = rhs_vec;
		  	#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		   Base::dt = Base::timeStep(Euler_sol_vec[nx][ny], mesh, nx, ny);
		  		   Base::getResidual(sol, Euler_sol_vec, closure_vec, res_vec[nx][ny], mesh, nx, ny);
		  		   Base::updateSolution(Euler_sol_vec[nx][ny], res_vec[nx][ny]);
		  		}
		  	}
		}
};

template <typename RESIDUAL>
class SymGSMacEqs : public EulerSolverBase<RESIDUAL>{
	public:
	    typedef EulerSolverBase<RESIDUAL> Base;
	    typedef typename RESIDUAL::dis_t dis_t;
	    typedef typename RESIDUAL::sol_t sol_t;
	    typedef typename sol_t::mesh_t mesh_t;
	    typedef typename dis_t::Velocity velocity_t;
	    typedef typename dis_t::Indices index_t;
	    typedef typename dis_t::value_type value_t;
	    typedef typename RESIDUAL::Problem Problem;
	    typedef typename RESIDUAL::Euler_sol_t Euler_sol_t;
	    typedef typename RESIDUAL::Euler_sol_vec_t Euler_sol_vec_t;
	public:
		void solve(const sol_t& sol, 
		  	  Euler_sol_vec_t& Euler_sol_vec, 
		  	  const Euler_sol_vec_t& closure_vec,
		  	  const Euler_sol_vec_t& rhs_vec,
		  	  const mesh_t& mesh){
		  	Base::M = sol[0][0].GetOrder();
		  	int steps = 0; double dif = 1.0;
		   while((steps < Base::maxSteps) && (dif > Base::tol)){
		  		Euler_sol_vec_t Euler_sol_vec_old = Euler_sol_vec;
		  		symSweeping(sol, Euler_sol_vec, closure_vec, rhs_vec, mesh);
		  		steps++;
				dif = Base::getDifNorm(Euler_sol_vec_old, Euler_sol_vec, mesh);
		  	}
		}

		void symSweeping(const sol_t& sol, 
		  	  Euler_sol_vec_t& Euler_sol_vec, 
		  	  const Euler_sol_vec_t& closure_vec,
		  	  const Euler_sol_vec_t& rhs_vec,
		  	  const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	//res_vec这个量在使用的时候，单层网格上要保证每次计算开始时是0
		  	//多重网格上计算开始时是方程右端项，这样，后面的累加操作才是正确的.
		  	Euler_sol_vec_t res_vec = rhs_vec;
		  	#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		   Base::dt = Base::timeStep(Euler_sol_vec[nx][ny], mesh, nx, ny);
		  		   Base::getResidual(sol, Euler_sol_vec, closure_vec, res_vec[nx][ny], mesh, nx, ny);
		  		   Base::updateSolution(Euler_sol_vec[nx][ny], res_vec[nx][ny]);
		  		}
		  	}
			
			res_vec = rhs_vec;
			#pragma omp parallel for
		  	for(size_t ny = Ny; ny > 0; ny--){
		  		for(size_t nx = Nx; nx > 0; nx--){
		  		    Base::dt = Base::timeStep(Euler_sol_vec[nx-1][ny-1], mesh, nx-1, ny-1);
		  		    Base::getResidual(sol, Euler_sol_vec, closure_vec, res_vec[nx-1][ny-1], mesh, nx-1, ny-1);
		  		    Base::updateSolution(Euler_sol_vec[nx-1][ny-1], res_vec[nx-1][ny-1]);
		  		}
		  	}
		}
};

template <typename RESIDUAL>
class FastSwpMacEqs : public EulerSolverBase<RESIDUAL>{
	public:
	    typedef EulerSolverBase<RESIDUAL> Base;
	    typedef typename RESIDUAL::dis_t dis_t;
	    typedef typename RESIDUAL::sol_t sol_t;
	    typedef typename sol_t::mesh_t mesh_t;
	    typedef typename dis_t::Velocity velocity_t;
	    typedef typename dis_t::Indices index_t;
	    typedef typename dis_t::value_type value_t;
	    typedef typename RESIDUAL::Problem Problem;
	    typedef typename RESIDUAL::Euler_sol_t Euler_sol_t;
	    typedef typename RESIDUAL::Euler_sol_vec_t Euler_sol_vec_t;
	public:
		void solve(const sol_t& sol, 
		  	  Euler_sol_vec_t& Euler_sol_vec, 
		  	  const Euler_sol_vec_t& closure_vec,
		  	  const Euler_sol_vec_t& rhs_vec,
		  	  const mesh_t& mesh){
		  	Base::M = sol[0][0].GetOrder();
		  	int steps = 0; double dif = 1.0;
		    while((steps < Base::maxSteps) && (dif > Base::tol)){
		  		Euler_sol_vec_t Euler_sol_vec_old = Euler_sol_vec;
				symSweeping(sol, Euler_sol_vec, closure_vec, rhs_vec, mesh);
		  		steps++;
				dif = Base::getDifNorm(Euler_sol_vec_old, Euler_sol_vec, mesh);
		  	}
		}

		void symSweeping(const sol_t& sol, 
		  	  Euler_sol_vec_t& Euler_sol_vec, 
		  	  const Euler_sol_vec_t& closure_vec,
		  	  const Euler_sol_vec_t& rhs_vec,
		  	  const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	//res_vec这个量在使用的时候，单层网格上要保证每次计算开始时是0
		  	//多重网格上计算开始时是方程右端项，这样，后面的累加操作才是正确的.
		  	Euler_sol_vec_t res_vec = rhs_vec;
		  	#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		    Base::dt = Base::timeStep(Euler_sol_vec[nx][ny], mesh, nx, ny);
		  		    Base::getResidual(sol, Euler_sol_vec, closure_vec, res_vec[nx][ny], mesh, nx, ny);
		  		    Base::updateSolution(Euler_sol_vec[nx][ny], res_vec[nx][ny]);
		  		}
		  	}

			res_vec = rhs_vec;
			#pragma omp parallel for
		  	for(size_t ny = Ny; ny > 0; ny--){
		  		for(size_t nx = Nx; nx > 0; nx--){
		  		    Base::dt = Base::timeStep(Euler_sol_vec[nx-1][ny-1], mesh, nx-1, ny-1);
		  		    Base::getResidual(sol, Euler_sol_vec, closure_vec, res_vec[nx-1][ny-1], mesh, nx-1, ny-1);
		  		    Base::updateSolution(Euler_sol_vec[nx-1][ny-1], res_vec[nx-1][ny-1]);
		  		}
		  	}

			res_vec = rhs_vec;
			#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = Nx; nx > 0; nx--){
		  		    Base::dt = Base::timeStep(Euler_sol_vec[nx-1][ny], mesh, nx-1, ny);
		  		    Base::getResidual(sol, Euler_sol_vec, closure_vec, res_vec[nx-1][ny], mesh, nx-1, ny);
		  		    Base::updateSolution(Euler_sol_vec[nx-1][ny], res_vec[nx-1][ny]);
		  		}
		  	}

			res_vec = rhs_vec;
			#pragma omp parallel for
		  	for(size_t ny = Ny; ny > 0; ny--){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		    Base::dt = Base::timeStep(Euler_sol_vec[nx][ny-1], mesh, nx, ny-1);
		  		    Base::getResidual(sol, Euler_sol_vec, closure_vec, res_vec[nx][ny-1], mesh, nx, ny-1);
		  		    Base::updateSolution(Euler_sol_vec[nx][ny-1], res_vec[nx][ny-1]);
		  		}
		  	}
		}
};

template <typename RESIDUAL>
class NavierStokes : public EulerSolverBase<RESIDUAL>{
	public:
	    typedef EulerSolverBase<RESIDUAL> Base;
	    typedef typename Base::sol_t sol_t;
	    typedef typename Base::dis_t dis_t;
	    typedef typename Base::mesh_t mesh_t;
	    typedef typename Base::Problem Problem;
	    typedef typename Base::Euler_sol_t Euler_sol_t;
	    typedef typename Base::Euler_sol_vec_t Euler_sol_vec_t;
	public:
		void solve(const sol_t& sol,
				Euler_sol_vec_t& Euler_sol_vec, 
		  	    Euler_sol_vec_t& closure_vec, 
		  		const Euler_sol_vec_t& rhs_vec,
		  		const mesh_t& mesh, 
		  		const int steps = 1){
			size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
			//设置初值
		  	for(size_t ny = 0; ny < Ny; ny++){
		       for(size_t nx = 0; nx < Nx; nx++){
				   Euler_sol_vec[nx][ny] << 1, 0, 0, 0, 1.5;
				   closure_vec[nx][ny].resize(10);
		       }
			}
			Base::M = sol[0][0].GetOrder();
			int NSsteps = 0;
			double res_norm=1;

			timeval start_time_NS;
		  	gettimeofday( &start_time_NS, NULL );

			std::ofstream os_res( "res_step_NS", std::ios::app );
		  	os_res.precision( 12 );

			//当NS方程解的残量下降到1e-8时停止迭代
			do{
			   //Base::res_vec = rhs_vec;
				Euler_sol_vec_t res_vec = rhs_vec;
			   NSsteps++;
				//计算时间步长
			   Base::timeStep(Euler_sol_vec, mesh);
				//求出应力张量和热通量
			   NSclosure(sol, Euler_sol_vec, closure_vec, mesh);
				//求残量
			   Base::getResidual(sol, Euler_sol_vec, closure_vec,
			  	 	 res_vec, mesh);
				//更新解
			   Base::updateSolution(Euler_sol_vec, res_vec, mesh);
				//解更新后求残量范数
			   res_norm = Base::getResidualNorm(sol, Euler_sol_vec, closure_vec,
			  	 	 res_vec, mesh );
				//显示到屏幕
			   std::cout << "steps = " << NSsteps << '\t' 
			  			  << "res_norm = " << res_norm << std::endl;
				timeval end_time_NS;
				gettimeofday( &end_time_NS, NULL );
				//保存步数，残量，时间
				os_res << NSsteps << " "
		  			    << res_norm << " "
		  				<< getElapsedTime(start_time_NS, end_time_NS) << std::endl ;
			}while(res_norm>1e-8);
			os_res << std::endl;
		  	os_res.close();
			outputSolution(Euler_sol_vec, closure_vec, mesh);
			exit(0);
		}
		
		double getElapsedTime(const timeval& _start, const timeval& _end){
			 return _end.tv_sec - _start.tv_sec + 
					 (_end.tv_usec - _start.tv_usec) / 1000000.;
		}
		
		void NSclosure(const sol_t& sol,
		  			     const Euler_sol_vec_t& Euler_sol_vec, 
		  			     Euler_sol_vec_t& closure_vec,
		  				  const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	#pragma omp parallel for
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		   NSclosure(sol, Euler_sol_vec, closure_vec, mesh, nx, ny);
		  		}
		  	}
		}

		void NSclosure(const sol_t& sol,
				const Euler_sol_vec_t& Euler_sol_vec, 
		  		Euler_sol_vec_t& closure_vec,
		  	  	const mesh_t& mesh, 
		  	  	size_t nx,
		  	  	size_t ny){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	Euler_sol_t Euler_sol_l, Euler_sol_m, Euler_sol_r,
		  					Euler_sol_top, Euler_sol_b, closure;
		  	if(nx == 0){
		  		 Base::leftBoundaryCondition(sol[nx][ny], Euler_sol_vec[nx][ny], Euler_sol_l, closure);
		  		 Euler_sol_m = Euler_sol_vec[nx][ny];
		  		 Euler_sol_r = Euler_sol_vec[nx+1][ny];
		  	}else if(nx == (Nx-1)){
		  		 Base::rightBoundaryCondition(sol[nx][ny], Euler_sol_vec[nx][ny], Euler_sol_r, closure);
		  		 Euler_sol_m = Euler_sol_vec[nx][ny];
		  		 Euler_sol_l = Euler_sol_vec[nx-1][ny];
		  	}else{
		  		 Euler_sol_l = Euler_sol_vec[nx-1][ny];
		  		 Euler_sol_m = Euler_sol_vec[nx][ny];
		  		 Euler_sol_r = Euler_sol_vec[nx+1][ny];
		  	}

		  	if(ny == 0){
		  		 Base::bottomBoundaryCondition(sol[nx][ny], Euler_sol_vec[nx][ny], Euler_sol_b, closure);
		  		 Euler_sol_top = Euler_sol_vec[ny][ny+1];
		  	}else if(ny == (Ny-1)){
		  		 Base::topBoundaryCondition(sol[nx][ny], Euler_sol_vec[nx][ny], Euler_sol_top, closure);
		  		 Euler_sol_b = Euler_sol_vec[nx][ny-1];
		  	}else{
		  		 Euler_sol_top = Euler_sol_vec[nx][ny+1];
		  		 Euler_sol_b = Euler_sol_vec[nx][ny-1];
		  	}

		  	getNSClosure(closure_vec[ny][nx], Euler_sol_m, Euler_sol_l, Euler_sol_r,
		  					Euler_sol_b, Euler_sol_top, mesh);
		}
		
		//中心差分
		void getNSClosure(Euler_sol_t& closure,
				const Euler_sol_t& solM, 
				const Euler_sol_t& solL,
		 		const Euler_sol_t& solR,
		  		const Euler_sol_t& solB,
		 		const Euler_sol_t& solT,
		 		const mesh_t& mesh){
		 	Euler_sol_t macL(5), macM(5), macR(5), macB(5), macT(5);
		 	getMac(macL, solL);
		 	getMac(macM, solM);
		 	getMac(macR, solR);
		  	getMac(macB, solB);
		  	getMac(macT, solT);
		 	//网格单元长度，均匀网格
		 	double dx = mesh.dx(0,0);
		  	double dy = mesh.dy(0,0);
		 	//求压力 p = rho * theta;
		 	double p = macM[0] * macM[4]; 
		 	//粘性指数
		 	double w = 0.81;
		 	//平均碰撞频率
		  	double nu = macM[0] * pow( macM[4], 1-w ) / Problem::Kn *
		  		 sqrt(M_PI/2);
		 	//p
		 	closure[0] = p;
		 	//sigma_{11}
		  	double uR = macR[1], uL = macL[1], vT = macT[2], vB = macB[2];
		 	closure[1] = -(4.0/3.0) * (p/nu) * (uR - uL) / (2*dx) + 
		  					  (2.0/3.0) * (p/nu) * (vT - vB) / (2*dy);
		 	//sigma_{12} = sigma_{21}
		  	double vR = macR[2], vL = macL[2], uT = macT[1], uB = macB[1];
		 	closure[2] = - (p/nu) * ( (vR - vL) / (2*dx) + (uT - uB) / (2*dy));
		 	//sigma_{13} = sigma_{31}
		  	double wR = macR[3], wL = macL[3];
		 	closure[3] = - (p/nu) * (wR - wL) / (2*dx);
		  	//sigma_{22}
		  	//double uR = macR[1], uL = macL[1], vT = macT[2], vB = macB[2];
		  	closure[4] = -(4.0/3.0) * (p/nu) * (vT - vB) / (2*dy) + 
		  					  (2.0/3.0) * (p/nu) * (uR - uL) / (2*dx);
		  	//sigma_{23} = sigma_{32}
		  	double wT = macT[3], wB = macB[3];
		  	closure[5] = - (p/nu) * (wT - wB) / (2*dy);
		 	//q_{1}
		  	double TR = macR[4], TL = macL[4];
		 	closure[7] = -(5.0/2.0) * (p/nu) * (TR - TL) / (2*dx);
		  	//q_{2}
		  	double TT = macT[4], TB = macB[4];
		 	closure[8] = -(5.0/2.0) * (p/nu) * (TT - TB) / (2*dy);
		}
	  
		void getMac(Euler_sol_t& mac, const Euler_sol_t& Euler_sol){
		  	//从守恒量中取出rho,u,v,w,theta
		  	mac[0] = Euler_sol[0];
		  	mac[1] = Euler_sol[1] / mac[0];
		  	mac[2] = Euler_sol[2] / mac[0];
		  	mac[3] = Euler_sol[3] / mac[0];
		  	mac[4] = Base::getTemperature(Euler_sol);
		}

		void outputSolution(const Euler_sol_vec_t& Euler_sol_vec,
		  						 const Euler_sol_vec_t& closure_vec,
		  						 const mesh_t& mesh){
		  	size_t Nx = mesh.n_ele_x();
		  	size_t Ny = mesh.n_ele_y();
		  	std::ofstream os("NS.dat");
		  	os.precision(12);
		  	os << "% variables= \"x\", \"y\",\"rho\", "
		  		<< " \"u\", \"v\", \"w\", \"T\", \"p\" " 
		  	   << " \"sigma11\", \"sigma12\", \"sigma13\" "
		  		<< " \"sigma22\", \"sigma23\", \"sigma33\" "
		  		<< " \"q1\", \"q2\", \"q3\" " 
		  		<< std::endl;
		  	os << "zone i = " << Nx << ",j=" << Ny 
		       << ", datapacking = \"point\" " << std::endl ;
		  	//std::vector<double> mac(5);
		  	std::vector<double> center(2);
			Euler_sol_t mac(5);
			//Euler_sol_t center(2);
		  	for(size_t ny = 0; ny < Ny; ny++){
		  		for(size_t nx = 0; nx < Nx; nx++){
		  		    mesh.center(center, nx, ny);
		  		    os << center[0] << '\t';
		  		    os << center[1] << '\t';
		  		    //宏观量
		  		    getMac(mac, Euler_sol_vec[ny][nx]);
		  		    for(size_t i = 0; i < 5; i++)
		  		   	   os << mac[i] << '\t';
		  		    //封闭项
		  		    for(size_t i = 0; i < 10; i++)
		  		   	   os << closure_vec[ny][nx][i] << '\t';
		  		    os << std::endl;
		  		}
		  	}
		}
};

template <typename RESIDUAL>
class NewtonIteMacEqs : public EulerSolverBase<RESIDUAL>{
	public:
		typedef RESIDUAL Residual;
		typedef EulerSolverBase<RESIDUAL> Base;
		typedef typename Base::sol_t sol_t;
		typedef typename Base::dis_t dis_t;
		typedef typename Base::mesh_t mesh_t;
		typedef typename Base::Problem Problem;
		typedef typename Base::Euler_sol_t Euler_sol_t;
		typedef typename Base::Euler_sol_vec_t Euler_sol_vec_t;
		typedef typename Eigen::MatrixXd mat_t;
		typedef typename std::vector<std::vector<mat_t>> mat_vec_t;
		//在这个地方又浪费时间。Eigen的文档已经说清楚了，高版本的编译器会
		//自动处理和Eigen数据结构相关的内存分配，编译的错误也指向
		//sol[0].GetOrder()，不仔细看提示，一直在这里纠结
		//typedef typename std::vector<mat_t, Eigen::aligned_allocator<mat_t>> mat_vec_t;
	private:
		double perturbation_ratio;    /**< 计算 Jacobian 矩阵数值微分扰动比例 */
  		double lambda;                /**< 正则化 Jacobian 矩阵的比例系数 */
  		double tau;                    /**< 更新解步长, 准确的 Newton 迭代应为 1 */
		size_t s1, s2;            /**< 线性多重网格光滑步数 */
		size_t maxMGSteps;    		  /**< 线性多重网格最大迭代步数 */
		size_t coarsestGrid;		  /**< 线性多重网格最粗网格规模 */
	public:
	    NewtonIteMacEqs() : perturbation_ratio(1e-6), lambda(4.0), tau(1.0),
			s1(2), s2(2), maxMGSteps(1), coarsestGrid(8) {};
	    ~NewtonIteMacEqs() {};
	public:
		 void dump_data( std::ostream& os ) const {
			 Base::dump_data(os);
			 os.write((const char *)&perturbation_ratio, sizeof(double));
			 os.write((const char *)&lambda, sizeof(double));
			 os.write((const char *)&tau, sizeof(double));
			 os.write((const char *)&s1, sizeof(size_t));
			 os.write((const char *)&s2, sizeof(size_t));
			 os.write((const char *)&maxMGSteps, sizeof(size_t));
			 os.write((const char *)&coarsestGrid, sizeof(size_t));
		 }

		 void load_data(std::istream& is){
			 Base::load_data(is);
			 is.read((char *)&perturbation_ratio, sizeof(double));
		  	 is.read((char *)&lambda, sizeof(double));
			 is.read((char *)&tau, sizeof(double));
			 is.read((char *)&s1, sizeof(size_t));
			 is.read((char *)&s2, sizeof(size_t));
		  	 is.read((char *)&maxMGSteps, sizeof(size_t));
			 is.read((char *)&coarsestGrid, sizeof(size_t));
		 }

		 void print_config(){
			 Base::print_config();
			 std::cout << "perturbation_ratio: " << perturbation_ratio << std::endl;
			 std::cout << "lambda: " << lambda << std::endl;
			 std::cout << "tau: " << tau << std::endl;
			 std::cout << "s1: " << s1 << std::endl;
			 std::cout << "s2: " << s2 << std::endl;
			 std::cout << "maxMGSteps: " << maxMGSteps << std::endl;
			 std::cout << "coarsestGrid: " << coarsestGrid << std::endl;
		 }
	public:
		class JacMatrix{
			public:
				mat_vec_t cent, left, right, below, up, ll, rr, bb, uu;
			public:
				JacMatrix(){}

				JacMatrix(size_t n, const mesh_t& mesh){
					size_t Nx = mesh.n_ele_x();
					size_t Ny = mesh.n_ele_y();
					if(n == 5){
						cent.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						left.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						right.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						below.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						up.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5)));
					}else if(n == 9){
						cent.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						left.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						right.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						below.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						up.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5)));
						ll.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						rr.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						bb.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						uu.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5)));
					}
				}		

				void resize(size_t n, const mesh_t& mesh){
					size_t Nx = mesh.n_ele_x();
					size_t Ny = mesh.n_ele_y();
					if(n == 5){
						cent.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						left.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						right.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						below.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						up.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5)));
					}else if(n == 9){
						cent.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						left.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						right.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						below.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						up.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5)));
						ll.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						rr.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						bb.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5))),
						uu.resize(Nx, std::vector<mat_t>(Ny, mat_t::Zero(5,5)));
					}
				}
		};

		void read_config(const ConfigFile& cf){
		  	Base::read_config(cf);
			perturbation_ratio = (double)cf.Value("SingleGridMethod", "Perturbation_ratio");
    		lambda = (double)cf.Value("SingleGridMethod", "RegularizeJacobiMatrixScaling");
    		tau = (double)cf.Value("SingleGridMethod", "UpdatingDisStepScaling");
    		s1 = (int)cf.Value("SingleGridMethod", "PreSmoothingSteps");
    		s2 = (int)cf.Value("SingleGridMethod", "PostSmoothingSteps");
			//s3 = (int)cf.Value("SingleGridMethod", "CoarsestGridSteps");
			coarsestGrid = (int)cf.Value("SingleGridMethod", "CoarsestGridSize");
			maxMGSteps = (int)cf.Value("SingleGridMethod", "MaxMGSteps");
		}
	
	public:
		void solve(const sol_t& sol,
		  		 Euler_sol_vec_t& Euler_sol_vec, 
		  		 const Euler_sol_vec_t& closure_vec, 
		  		 const Euler_sol_vec_t& rhs_vec,
		  		 const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
			Base::M = sol[0][0].GetOrder();
			int steps = 0; double norm = 1;
			//设置Newton迭代终止条件，最大步数，或解的变化足够小
			//double Time = 0;
			//auto start = std::chrono::high_resolution_clock::now();
			while((steps < Base::maxSteps) && (norm > 1e-12)){
				//这里注意要将delU置为0
				Euler_sol_vec_t delU(Nx, std::vector<Euler_sol_t>(Ny, Euler_sol_t::Zero(5)));
				Euler_sol_vec_t res_vec(Nx, std::vector<Euler_sol_t>(Ny, Euler_sol_t::Zero(5)));
				res_vec = rhs_vec;
				Base::getResidual(sol, Euler_sol_vec, closure_vec, res_vec, mesh);
				#ifdef RECONSTRUCT
				JacMatrix jacMat(9, mesh);
				getJacMatSecDis(sol, Euler_sol_vec, closure_vec, res_vec, rhs_vec, mesh, jacMat);
				solveLinEqsSecDis(jacMat, delU, res_vec, mesh);
				#else
				JacMatrix jacMat(5, mesh);
				getJacMatFirDis(sol, Euler_sol_vec, closure_vec, res_vec, rhs_vec, mesh, jacMat);
				solveLinEqsFirDis(jacMat, delU, res_vec, mesh);
				#endif
				update(delU, mesh, Euler_sol_vec);
				steps++;
		   }
			//auto end = std::chrono::high_resolution_clock::now();
			//Time = std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count();
			//std::cout << "Elapsed time: " << Time << " ms\n";
		}

		void getJacMatFirDis(const sol_t& sol,
		  		 const Euler_sol_vec_t& Euler_sol_vec, 
		  		 const Euler_sol_vec_t& closure_vec, 
		  		 const Euler_sol_vec_t& res_vec,
				 const Euler_sol_vec_t& rhs_vec,
		  		 const mesh_t& mesh, 
		  		 JacMatrix& jacMat){
			size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
			Euler_sol_vec_t temp = Euler_sol_vec;
			#pragma omp parallel for 
		   for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					//传入当前考虑的单元位置和需要扰动的自变量的位置
					Eigen::Vector2i p1, p2;
					p1 << nx, ny;
					getMatrix(sol, temp, closure_vec, mesh, p1, p1, res_vec, rhs_vec, jacMat.cent);
					//左边界	
					if(nx == 0){
						p2 << nx+1, ny;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.right);
					}
					//右边界
					else if (nx == Nx-1){
						p2 << nx-1, ny;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.left);
					}else{
						p2 << nx-1, ny;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.left);
						p2 << nx+1, ny;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.right);
					}

					//下边界
					if (ny == 0){
						p2 << nx, ny+1;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.up);
					}
					//上边界
					else if (ny == Ny-1){
						p2 << nx, ny-1;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.below);
					}
					else{
						p2 << nx, ny-1;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.below);
						p2 << nx, ny+1;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.up);
					}
				}
			}
		}

		void getMatrix(const sol_t& sol,
				 Euler_sol_vec_t& temp,
		  		 const Euler_sol_vec_t& closure_vec, 
		  		 const mesh_t& mesh,
		  		 const Eigen::Vector2i& p1,
		  		 const Eigen::Vector2i& p2,
		  		 const Euler_sol_vec_t& res_vec,
				 const Euler_sol_vec_t& rhs_vec,
				 mat_vec_t& mat_vec){
			size_t nx, ny;
		  	for(size_t i = 0; i < 5; i++){
		  		//对编号p2的单元上第i个自变量做扰动
				nx = p2[0], ny = p2[1];
				double turb = perturbation_ratio * temp[nx][ny][i] + 1e-8;
		  		temp[nx][ny][i] += turb;
		  		//为扰动后的残量申请空间
				nx = p1[0], ny = p1[1];
		  		Euler_sol_t resTurb = rhs_vec[nx][ny];
		  		//求扰动后的残量
		  		Base::getResidual(sol, temp, closure_vec, resTurb, mesh, nx, ny);
		  		//生成雅可比矩阵的第i列
		  		for(size_t j = 0; j < 5; j++){
		  		    mat_vec[nx][ny](j,i) = -(resTurb[j] - res_vec[nx][ny][j]) / turb;
		  		    if(p1 == p2 && i == j){
						 //用局部残量的范数和Knudsen数对矩阵作正则化 
						 mat_vec[nx][ny](j,i) += lambda * (Problem::Kn + res_vec[nx][ny].norm());
		  		    }
		  		}
		  		//把扰动去掉
				nx = p2[0], ny = p2[1];
		  		temp[nx][ny][i] -= turb;
		  	}
		}

		void solveLinEqsFirDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			for(size_t s = 0; s < maxMGSteps; s++){
				mgOneStepFirDis(jacMat, sol_vec, rhs_vec, mesh);
			}
		}

		void update(const Euler_sol_vec_t& delU,
		  	  const mesh_t& mesh,
		  	  Euler_sol_vec_t& Euler_sol_vec){
			size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
			#pragma omp parallel for
		    for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
		  	  		Euler_sol_vec[nx][ny] += tau * delU[nx][ny];
				}
		    }
		}

		void symmetrySweepingFirDis(const JacMatrix& jacMat,
				 Euler_sol_vec_t& delU,
		  		 const Euler_sol_vec_t& rhs_vec,
		  		 const mesh_t& mesh, 
				 const size_t steps = 1){
			size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
		    for(size_t s = 0; s < steps; s++){
				//在线性多重网格方法的各个步骤使用了并行计算．
				//由于网格的粗化，这时候需要注意粗网格规模与线程数的关系
				//为了保证计算稳定，线程数应该要小于最粗网格的网格数
				//现在常用的设置是最粗网格8*8，线程数4.
				//所以使用该求解器，不能一味地增加并行规模
				#pragma omp parallel for
				for(size_t ny = 0; ny < Ny; ny++){
					for(size_t nx = 0; nx < Nx; nx++){
						Eigen::Vector2i p;
						p << nx, ny;
						sweepingFirDis(jacMat, rhs_vec, mesh, delU, p);
					}
				}

				#pragma omp parallel for
				for(size_t ny = Ny; ny > 0; ny--){
					for(size_t nx = Nx; nx > 0; nx--){
						Eigen::Vector2i p;
						p << nx-1, ny-1;
						sweepingFirDis(jacMat, rhs_vec, mesh, delU, p);
					}
				}
			}
		}

		void sweepingFirDis(const JacMatrix& jacMat,
				const Euler_sol_vec_t& rhs_vec,
		  		const mesh_t& mesh, 
		  		Euler_sol_vec_t& delU,
		  		const Eigen::Vector2i& p){
		   	size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
			size_t nx = p[0], ny = p[1];
		   	Euler_sol_t right(5);
			right = rhs_vec[nx][ny];
			if(nx == 0){
			 	right += - jacMat.right[nx][ny] * delU[nx+1][ny];
			}else if(nx == Nx-1){
			 	right += - jacMat.left[nx][ny] * delU[nx-1][ny];
			}else{
			 	right += - jacMat.left[nx][ny] * delU[nx-1][ny] 
						 - jacMat.right[nx][ny] * delU[nx+1][ny]; 
			}

			if(ny == 0){
		  	  	right += - jacMat.up[nx][ny] * delU[nx][ny+1];
		   	}else if(ny == Ny-1){
		  	 	right += - jacMat.below[nx][ny] * delU[nx][ny-1];
		   	}else{
		  	 	right += - jacMat.up[nx][ny] * delU[nx][ny+1]
						 - jacMat.below[nx][ny] * delU[nx][ny-1]; 
		   	}

		   	delU[nx][ny] = jacMat.cent[nx][ny].inverse() * right;
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
			//生成满矩阵并求逆
			exactSolveFirDis(jacMat, sol_vec, rhs_vec, mesh);
			//迭代法求解s3次
			//symmetrySweepingFirDis(jacMat, sol_vec, rhs_vec, mesh, s3);
		}

		void exactSolveFirDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
			size_t n_ele = Nx * Ny;
			mat_t fullJacMat = Eigen::MatrixXd::Zero(5*n_ele, 5*n_ele);
			Eigen::VectorXd fullRHS = Eigen::VectorXd::Zero(5*n_ele);
			getFullMatFir(fullJacMat, jacMat, mesh);
			getRHS(fullRHS, rhs_vec, mesh);
			Eigen::VectorXd fullSOL = Eigen::VectorXd::Zero(5*n_ele);
			fullSOL = fullJacMat.inverse() * fullRHS;
			#pragma omp parallel for
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					size_t ind = ny * Nx + nx;
					sol_vec[nx][ny] =  fullSOL.segment(ind,5);
				}
			}
		}	

		void getFullMatFir(mat_t& fullJacMat,
				const JacMatrix& jacMat,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
			#pragma omp parallel for
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					size_t row = ny * Nx + nx, col = ny * Nx + nx;
					if(nx == 0){
						fullJacMat.block<5,5>(row*5,col*5) = jacMat.cent[nx][ny];
						fullJacMat.block<5,5>(row*5,(col+1)*5) = jacMat.right[nx][ny];
					}else if(nx == Nx-1){
						fullJacMat.block<5,5>(row*5,(col-1)*5) = jacMat.left[nx][ny];
						fullJacMat.block<5,5>(row*5,col*5) = jacMat.cent[nx][ny];
					}else{
						fullJacMat.block<5,5>(row*5,(col-1)*5) = jacMat.left[nx][ny];
						fullJacMat.block<5,5>(row*5,col*5) = jacMat.cent[nx][ny];
						fullJacMat.block<5,5>(row*5,(col+1)*5) = jacMat.right[nx][ny];
					}
					
					if(ny == 0){
						col = (ny+1)*Nx + nx;
						fullJacMat.block<5,5>(row*5,col*5) = jacMat.up[nx][ny];
					}else if(ny == Ny-1){
						col = (ny-1)*Nx + nx;
						fullJacMat.block<5,5>(row*5,col*5) = jacMat.below[nx][ny];
					}else{
						col = (ny+1)*Nx + nx;
						fullJacMat.block<5,5>(row*5,col*5) = jacMat.up[nx][ny];
						col = (ny-1)*Nx + nx;
						fullJacMat.block<5,5>(row*5,col*5) = jacMat.below[nx][ny];
					}
				}
			}
		}

		void getRHS(Eigen::VectorXd& fullRHS,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
			#pragma omp parallel for
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					size_t row = ny * Nx + nx, col = ny * Nx + nx;
					fullRHS.segment(row*5,5) = rhs_vec[nx][ny];
				}
			}
		}

		void subGridIteration(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			mgOneStepFirDis(jacMat, sol_vec, rhs_vec, mesh);
		}

		void restrictionFirDis(const JacMatrix& Mat,
				JacMatrix& coarMat,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
			coarMat.resize(5, mesh);	
			#pragma omp parallel for
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					coarMat.cent[nx][ny] =  
					(Mat.cent[2*nx][2*ny] + Mat.right[2*nx][2*ny] + Mat.up[2*nx][2*ny] +
					 Mat.cent[2*nx+1][2*ny] + Mat.left[2*nx+1][2*ny] + Mat.up[2*nx+1][2*ny] + 
					 Mat.cent[2*nx][2*ny+1] + Mat.below[2*nx][2*ny+1] + Mat.right[2*nx][2*ny+1] + 
					 Mat.cent[2*nx+1][2*ny+1] + Mat.left[2*nx+1][2*ny+1] + Mat.below[2*nx+1][2*ny+1]);

					if(nx == 0){
						coarMat.right[nx][ny] = Mat.right[2*nx+1][2*ny] + Mat.right[2*nx+1][2*ny+1];
					}
					else if(nx == Nx-1){
						coarMat.left[nx][ny] = Mat.left[2*nx][2*ny] + Mat.left[2*nx][2*ny+1];
					}
					else{
						coarMat.left[nx][ny] = Mat.left[2*nx][2*ny] + Mat.left[2*nx][2*ny+1];
						coarMat.right[nx][ny] = Mat.right[2*nx+1][2*ny] + Mat.right[2*nx+1][2*ny+1];
					}

					if(ny == 0){
						coarMat.up[nx][ny] = Mat.up[2*nx][2*ny+1] + Mat.up[2*nx+1][2*ny+1];
					}
					else if(ny == Ny-1){
						coarMat.below[nx][ny] = Mat.below[2*nx][2*ny] + Mat.below[2*nx+1][2*ny];
					}
					else{
						coarMat.below[nx][ny] = Mat.below[2*nx][2*ny] + Mat.below[2*nx+1][2*ny];
						coarMat.up[nx][ny] = Mat.up[2*nx][2*ny+1] + Mat.up[2*nx+1][2*ny+1];
					}
				}
			}
		}
		
		void restriction(const Euler_sol_vec_t& sol_vec,
				Euler_sol_vec_t& coarser_sol_vec,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
			coarser_sol_vec.resize(Nx, std::vector<Euler_sol_t>(Ny));
			#pragma omp parallel for
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					coarser_sol_vec[nx][ny] = (sol_vec[2*nx][2*ny] 
					+ sol_vec[2*nx+1][2*ny] + sol_vec[2*nx][2*ny+1]
					+ sol_vec[2*nx+1][2*ny+1]);
				}
			}
		}

		void backSubstitution(Euler_sol_vec_t& correction,
				Euler_sol_vec_t& sol_vec,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
			double relax = 1.0;
			#pragma omp parallel for
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					sol_vec[2*nx][2*ny] += relax * correction[nx][ny];
					sol_vec[2*nx+1][2*ny] += relax * correction[nx][ny];
					sol_vec[2*nx][2*ny+1] += relax * correction[nx][ny];
					sol_vec[2*nx+1][2*ny+1] += relax * correction[nx][ny];
				}
			}
		}

		void mgOneStepFirDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
			//if it is the coarsest grid
			if(Nx == coarsestGrid){
				coarsestSolveFirDis(jacMat, sol_vec, rhs_vec, mesh);
				return;
			}

			preSmoothFirDis(jacMat, sol_vec, rhs_vec, mesh);

			{
				//calculate the right hand side of the coarser grid problem 
				Euler_sol_vec_t res_vec(Nx, std::vector<Euler_sol_t>(Ny, Euler_sol_t::Zero(5))); 
				getResFirDis(jacMat, sol_vec, rhs_vec, res_vec, mesh);
				mesh_t coar_mesh = mesh;
				coar_mesh.Coarser();
				Nx = coar_mesh.n_ele_x(), Ny = coar_mesh.n_ele_y();
				Euler_sol_vec_t coarser_rhs_vec(Nx, std::vector<Euler_sol_t>(Ny, Euler_sol_t::Zero(5)));
				restriction(res_vec, coarser_rhs_vec, coar_mesh);

				JacMatrix coarMat;
				restrictionFirDis(jacMat, coarMat, coar_mesh);

				//recursively call NMG iteration
				//粗网格上的解向量
				Euler_sol_vec_t coarser_sol_vec(Nx, std::vector<Euler_sol_t>(Ny, Euler_sol_t::Zero(5)));
				subGridIteration(coarMat, coarser_sol_vec, coarser_rhs_vec, coar_mesh);

				//update fine grid solution
				backSubstitution(coarser_sol_vec, sol_vec, coar_mesh);
			}

			postSmoothFirDis(jacMat, sol_vec, rhs_vec, mesh);
		}

		void getResFirDis(const JacMatrix& jacMat,
				const Euler_sol_vec_t& sol_vec, 
				const Euler_sol_vec_t& rhs_vec,
				Euler_sol_vec_t& res_vec,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
			Euler_sol_t right;
			#pragma omp parallel for
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					res_vec[nx][ny] = rhs_vec[nx][ny] 
									- jacMat.cent[nx][ny] * sol_vec[nx][ny];
					if(nx == 0){
						res_vec[nx][ny] += - jacMat.right[nx][ny] * sol_vec[nx+1][ny];						
					}
					else if(nx == Nx-1){
						res_vec[nx][ny] += - jacMat.left[nx][ny] * sol_vec[nx-1][ny];
					}
					else{
						res_vec[nx][ny] += - (jacMat.left[nx][ny] * sol_vec[nx-1][ny]
										   +  jacMat.right[nx][ny] * sol_vec[nx+1][ny]);
					}

					if(ny == 0){
						res_vec[nx][ny] += - jacMat.up[nx][ny] * sol_vec[nx][ny+1];
					}
					else if(ny == Ny-1){
						res_vec[nx][ny] += - jacMat.below[nx][ny] * sol_vec[nx][ny-1];
					}
					else{
						res_vec[nx][ny] += - (jacMat.below[nx][ny] * sol_vec[nx][ny-1]
										   +  jacMat.up[nx][ny] * sol_vec[nx][ny+1]);
					}
				}
			}
		}

		//二阶离散
		void getJacMatSecDis(const sol_t& sol,
		  		 const Euler_sol_vec_t& Euler_sol_vec, 
		  		 const Euler_sol_vec_t& closure_vec, 
		  		 const Euler_sol_vec_t& res_vec,
				 const Euler_sol_vec_t& rhs_vec,
		  		 const mesh_t& mesh, 
		  		 JacMatrix& jacMat){
			size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
			getJacMatFirDis(sol, Euler_sol_vec, closure_vec, res_vec, rhs_vec, mesh, jacMat);
			Euler_sol_vec_t temp = Euler_sol_vec;
			#pragma omp parallel for 
		   for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					//传入当前考虑的单元位置和需要扰动的自变量的位置
				   Eigen::Vector2i p1, p2;
					p1 << nx, ny;
					//左边界	
					if(nx == 0 || nx == 1){
						p2 << nx+2, ny;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.rr);
					}
					//右边界
					else if (nx == Nx-1 || nx == Nx-2){
						p2 << nx-2, ny;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.ll);
					}else{
						p2 << nx-2, ny;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.ll);
						p2 << nx+2, ny;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.rr);
					}

					//下边界
					if (ny == 0 || ny == 1){
						p2 << nx, ny+2;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.uu);
					}
					//上边界
					else if (ny == Ny-1 || ny == Ny-2){
						p2 << nx, ny-2;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.bb);
					}
					else{
						p2 << nx, ny-2;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.bb);
						p2 << nx, ny+2;
						getMatrix(sol, temp, closure_vec, mesh, p1, p2, res_vec, rhs_vec, jacMat.uu);
					}
				}
		   }
		}

		void solveLinEqsSecDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			for(size_t s = 0; s < maxMGSteps; s++){
				mgOneStepSecDis(jacMat, sol_vec, rhs_vec, mesh);
				//symmetrySweepingSecDis(jacMat, sol_vec, rhs_vec, mesh);
			}
		}

		void symmetrySweepingSecDis(const JacMatrix& jacMat,
				 Euler_sol_vec_t& delU,
		  		 const Euler_sol_vec_t& rhs_vec,
		  		 const mesh_t& mesh, 
				 const size_t steps = 1){
			size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
		    for(size_t s = 0; s < steps; s++){
				#pragma omp parallel for
				for(size_t ny = 0; ny < Ny; ny++){
					for(size_t nx = 0; nx < Nx; nx++){
						Eigen::Vector2i p;
						p << nx, ny;
						sweepingSecDis(jacMat, rhs_vec, mesh, delU, p);
					}
				}

				#pragma omp parallel for
				for(size_t ny = Ny; ny > 0; ny--){
					for(size_t nx = Nx; nx > 0; nx--){
						Eigen::Vector2i p;
						p << nx-1, ny-1;
						sweepingSecDis(jacMat, rhs_vec, mesh, delU, p);
					}
				}
			}
		}

		void sweepingSecDis(const JacMatrix& jacMat,
				const Euler_sol_vec_t& rhs_vec,
		  		const mesh_t& mesh, 
		  		Euler_sol_vec_t& delU,
		  		const Eigen::Vector2i& p){
		   	size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
			size_t nx = p[0], ny = p[1];
		   	Euler_sol_t right(5);
			right = rhs_vec[nx][ny];
			if(nx == 0){
			 	right += - jacMat.right[nx][ny] * delU[nx+1][ny]
						 - jacMat.rr[nx][ny] * delU[nx+2][ny];
			}else if(nx == 1){
				right += - jacMat.left[nx][ny] * delU[nx-1][ny]
						 - jacMat.right[nx][ny] * delU[nx+1][ny]
						 - jacMat.rr[nx][ny] * delU[nx+2][ny];
			}else if(nx == Nx-2){
			 	right += - jacMat.ll[nx][ny] * delU[nx-2][ny]
						 - jacMat.left[nx][ny] * delU[nx-1][ny]
						 - jacMat.right[nx][ny] * delU[nx+1][ny];
			}else if(nx == Nx-1){
			 	right += - jacMat.ll[nx][ny] * delU[nx-2][ny]
						 - jacMat.left[nx][ny] * delU[nx-1][ny];
			}else{
			 	right += - jacMat.ll[nx][ny] * delU[nx-2][ny]
						 - jacMat.left[nx][ny] * delU[nx-1][ny] 
						 - jacMat.right[nx][ny] * delU[nx+1][ny]
						 - jacMat.rr[nx][ny] * delU[nx+2][ny];
			}

			if(ny == 0){
		  	  	right += - jacMat.uu[nx][ny] * delU[nx][ny+2]
						 - jacMat.up[nx][ny] * delU[nx][ny+1];
		   	}else if(ny == 1){
		  	 	right += - jacMat.uu[nx][ny] * delU[nx][ny+2]
						 - jacMat.up[nx][ny] * delU[nx][ny+1];
						 - jacMat.below[nx][ny] * delU[nx][ny-1];
		   	}else if(ny == Ny-2){
		  	 	right += - jacMat.up[nx][ny] * delU[nx][ny+1];
						 - jacMat.below[nx][ny] * delU[nx][ny-1]
						 - jacMat.bb[nx][ny] * delU[nx][ny-2];
		   	}else if(ny == Ny-1){
		  	 	right += - jacMat.below[nx][ny] * delU[nx][ny-1]
						 - jacMat.bb[nx][ny] * delU[nx][ny-2];
		   	}else{
		  	 	right += - jacMat.bb[nx][ny] * delU[nx][ny-2]
						 - jacMat.below[nx][ny] * delU[nx][ny-1]
						 - jacMat.up[nx][ny] * delU[nx][ny+1]
						 - jacMat.uu[nx][ny] * delU[nx][ny+2];
		   	}

		   	delU[nx][ny] = jacMat.cent[nx][ny].inverse() * right;
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
			//对最粗网格上的系统作精确求解
			exactSolveSecDis(jacMat, sol_vec, rhs_vec, mesh);
			//使用迭代法求解s3次
			//symmetrySweepingFirDis(jacMat, sol_vec, rhs_vec, mesh, s3);
		}

		void exactSolveSecDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
			size_t n_ele = Nx * Ny;
			mat_t fullJacMat = Eigen::MatrixXd::Zero(5*n_ele, 5*n_ele);
			Eigen::VectorXd fullRHS = Eigen::VectorXd::Zero(5*n_ele);
			getFullMatSec(fullJacMat, jacMat, mesh);
			getRHS(fullRHS, rhs_vec, mesh);
			Eigen::VectorXd fullSOL = Eigen::VectorXd::Zero(5*n_ele);
			fullSOL = fullJacMat.inverse() * fullRHS;
			#pragma omp parallel for
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					size_t ind = ny * Nx + nx;
					sol_vec[nx][ny] =  fullSOL.segment(ind,5);
				}
			}
		}

		void getFullMatSec(mat_t& fullJacMat,
				const JacMatrix& jacMat,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
			getFullMatFir(fullJacMat, jacMat, mesh);
			#pragma omp parallel for
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					size_t row = ny * Nx + nx, col = ny * Nx + nx;
					if(nx == 0 || nx == 1){
						fullJacMat.block<5,5>(row,(col+2)*5) = jacMat.rr[nx][ny];
					}else if(nx == Nx-2 || nx == Nx-1){
						fullJacMat.block<5,5>(row*5,(col-2)*5) = jacMat.ll[nx][ny];
					}else{
						fullJacMat.block<5,5>(row*5,(col-2)*5) = jacMat.ll[nx][ny];
						fullJacMat.block<5,5>(row*5,(col+2)*5) = jacMat.rr[nx][ny];
					}

					if(ny == 0 || ny == 1){
						col = (ny+2)*Nx + nx;
						fullJacMat.block<5,5>(row,col*5) = jacMat.uu[nx][ny];
					}else if(ny == Ny-2 || ny == Ny-1){
						col = (ny-2)*Nx + nx;
						fullJacMat.block<5,5>(row*5,col*5) = jacMat.bb[nx][ny];
					}else{
						col = (ny+2)*Nx + nx;
						fullJacMat.block<5,5>(row*5,col*5) = jacMat.uu[nx][ny];
						col = (ny-2)*Nx + nx;
						fullJacMat.block<5,5>(row*5,col*5) = jacMat.bb[nx][ny];
					}
				}
			}
		}		

		void subGridIteSecDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			mgOneStepFirDis(jacMat, sol_vec, rhs_vec, mesh);
		}

		void restrictionSecDis(const JacMatrix& Mat,
				JacMatrix& coarMat,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x();
			size_t Ny = mesh.n_ele_y();
			coarMat.resize(5, mesh);	
			#pragma omp parallel for
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					coarMat.cent[nx][ny] =  
					(Mat.cent[2*nx][2*ny] + Mat.right[2*nx][2*ny] + Mat.up[2*nx][2*ny] +
					 Mat.cent[2*nx+1][2*ny] + Mat.left[2*nx+1][2*ny] + Mat.up[2*nx+1][2*ny] + 
					 Mat.cent[2*nx][2*ny+1] + Mat.below[2*nx][2*ny+1] + Mat.right[2*nx][2*ny+1] + 
					 Mat.cent[2*nx+1][2*ny+1] + Mat.left[2*nx+1][2*ny+1] + Mat.below[2*nx+1][2*ny+1]);

					if(nx == 0){
						coarMat.right[nx][ny] = 
							(Mat.right[2*nx+1][2*ny] + Mat.right[2*nx+1][2*ny+1] 
							 + Mat.rr[2*nx+1][2*ny] + Mat.rr[2*nx+1][2*ny+1]
							 + Mat.rr[2*nx][2*ny] + Mat.rr[2*nx][2*ny+1]);
					}
					else if(nx == Nx-1){
						coarMat.left[nx][ny] = 
							(Mat.left[2*nx][2*ny] + Mat.left[2*nx][2*ny+1] 
							 + Mat.ll[2*nx][2*ny] + Mat.ll[2*nx][2*ny+1]
							 + Mat.ll[2*nx+1][2*ny] + Mat.ll[2*nx+1][2*ny+1]);
					}
					else{
						coarMat.left[nx][ny] = (
							Mat.left[2*nx][2*ny] + Mat.left[2*nx][2*ny+1] 
							+ Mat.ll[2*nx][2*ny] + Mat.ll[2*nx][2*ny+1]
							+ Mat.ll[2*nx+1][2*ny] + Mat.ll[2*nx+1][2*ny+1]);
						coarMat.right[nx][ny] = (
							Mat.right[2*nx+1][2*ny] + Mat.right[2*nx+1][2*ny+1]
							+ Mat.rr[2*nx+1][2*ny] + Mat.rr[2*nx+1][2*ny+1]
							+ Mat.rr[2*nx][2*ny] + Mat.rr[2*nx][2*ny+1]);
					}

					if(ny == 0){
						coarMat.up[nx][ny] = (
							Mat.up[2*nx][2*ny+1] + Mat.up[2*nx+1][2*ny+1]
							+ Mat.uu[2*nx][2*ny+1] + Mat.uu[2*nx+1][2*ny+1]
							+ Mat.uu[2*nx][2*ny] + Mat.uu[2*nx+1][2*ny]);
					}
					else if(ny == Ny-1){
						coarMat.below[nx][ny] = (
							Mat.below[2*nx][2*ny] + Mat.below[2*nx+1][2*ny]
							+ Mat.bb[2*nx][2*ny] + Mat.bb[2*nx+1][2*ny]
							+ Mat.bb[2*nx][2*ny+1] + Mat.bb[2*nx+1][2*ny+1]);

					}else{
						coarMat.below[nx][ny] = (
							Mat.below[2*nx][2*ny] + Mat.below[2*nx+1][2*ny]
							+ Mat.bb[2*nx][2*ny] + Mat.bb[2*nx+1][2*ny]
							+ Mat.bb[2*nx][2*ny+1] + Mat.bb[2*nx+1][2*ny+1]);
						coarMat.up[nx][ny] = (
							Mat.up[2*nx][2*ny+1] + Mat.up[2*nx+1][2*ny+1]
							+ Mat.uu[2*nx][2*ny+1] + Mat.uu[2*nx+1][2*ny+1]
							+ Mat.uu[2*nx][2*ny] + Mat.uu[2*nx+1][2*ny]);
					}
				}
			}
		}
		
		void mgOneStepSecDis(const JacMatrix& jacMat,
				Euler_sol_vec_t& sol_vec,
				const Euler_sol_vec_t& rhs_vec,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
			//if it is the coarsest grid
			if(Nx == coarsestGrid){
				coarsestSolveSecDis(jacMat, sol_vec, rhs_vec, mesh);
				return;
			}

			preSmoothSecDis(jacMat, sol_vec, rhs_vec, mesh);

			{
				//calculate the right hand side of the coarser grid problem 
				Euler_sol_vec_t res_vec(Nx, std::vector<Euler_sol_t>(Ny, Euler_sol_t::Zero(5))); 
				getResSecDis(jacMat, sol_vec, rhs_vec, res_vec, mesh);
				mesh_t coar_mesh = mesh;
				coar_mesh.Coarser();
				Nx = coar_mesh.n_ele_x(), Ny = coar_mesh.n_ele_y();
				Euler_sol_vec_t coarser_rhs_vec(Nx, std::vector<Euler_sol_t>(Ny, Euler_sol_t::Zero(5)));
				restriction(res_vec, coarser_rhs_vec, coar_mesh);

				JacMatrix coarMat;
				restrictionSecDis(jacMat, coarMat, coar_mesh);

				//recursively call NMG iteration
				//粗网格上的解向量
				Euler_sol_vec_t coarser_sol_vec(Nx, std::vector<Euler_sol_t>(Ny, Euler_sol_t::Zero(5)));
				subGridIteSecDis(coarMat, coarser_sol_vec, coarser_rhs_vec, coar_mesh);

				//update fine grid solution
				backSubstitution(coarser_sol_vec, sol_vec, coar_mesh);
			}

			postSmoothSecDis(jacMat, sol_vec, rhs_vec, mesh);
		}

		void getResSecDis(const JacMatrix& jacMat,
				const Euler_sol_vec_t& sol_vec, 
				const Euler_sol_vec_t& rhs_vec,
				Euler_sol_vec_t& res_vec,
				const mesh_t& mesh){
			size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
			Euler_sol_t right;
			#pragma omp parallel for
			for(size_t ny = 0; ny < Ny; ny++){
				for(size_t nx = 0; nx < Nx; nx++){
					res_vec[nx][ny] = rhs_vec[nx][ny] - jacMat.cent[nx][ny] * sol_vec[nx][ny];
					if(nx == 0){
						res_vec[nx][ny] += - jacMat.right[nx][ny] * sol_vec[nx+1][ny]
							- jacMat.rr[nx][ny] * sol_vec[nx+2][ny];						
					}else if(nx == 1){
						res_vec[nx][ny] += - jacMat.left[nx][ny] * sol_vec[nx-1][ny]
							- jacMat.right[nx][ny] * sol_vec[nx+1][ny]
							- jacMat.rr[nx][ny] * sol_vec[nx+2][ny];
					}else if(nx == Nx-2){
						res_vec[nx][ny] += - jacMat.ll[nx][ny] * sol_vec[nx-2][ny]
							- jacMat.left[nx][ny] * sol_vec[nx-1][ny]
							- jacMat.right[nx][ny] * sol_vec[nx+1][ny];
					}else if(nx == Nx-1){
						res_vec[nx][ny] += - jacMat.ll[nx][ny] * sol_vec[nx-2][ny]
							- jacMat.left[nx][ny] * sol_vec[nx-1][ny];
					}
					else{
						res_vec[nx][ny] += - jacMat.ll[nx][ny] * sol_vec[nx-2][ny]
							- jacMat.left[nx][ny] * sol_vec[nx-1][ny]
							- jacMat.right[nx][ny] * sol_vec[nx+1][ny]
							- jacMat.rr[nx][ny] * sol_vec[nx+2][ny];
					}

					if(ny == 0){
						res_vec[nx][ny] += - jacMat.up[nx][ny] * sol_vec[nx][ny+1]
							- jacMat.uu[nx][ny] * sol_vec[nx][ny+2];
					}else if(ny == 1){
						res_vec[nx][ny] += - jacMat.below[nx][ny] * sol_vec[nx][ny-1]
							- jacMat.up[nx][ny] * sol_vec[nx][ny+1]
							- jacMat.uu[nx][ny] * sol_vec[nx][ny+2];
					}else if(ny == Ny-2){
						res_vec[nx][ny] += - jacMat.bb[nx][ny] * sol_vec[nx][ny-2]
							- jacMat.below[nx][ny] * sol_vec[nx][ny-1]
							- jacMat.up[nx][ny] * sol_vec[nx][ny+1];
					}else if(ny == Ny-1){
						res_vec[nx][ny] += - jacMat.bb[nx][ny] * sol_vec[nx][ny-2]
							- jacMat.below[nx][ny] * sol_vec[nx][ny-1];
					}
					else{
						res_vec[nx][ny] += - jacMat.bb[nx][ny] * sol_vec[nx][ny-2]
							- jacMat.below[nx][ny] * sol_vec[nx][ny-1]
							- jacMat.up[nx][ny] * sol_vec[nx][ny+1]
							- jacMat.uu[nx][ny] * sol_vec[nx][ny+2];
					}
				}
			}
		}
}; 

# endif
