/**
 * @file MultigridMethod.h
 * @brief 多重网格求解器
 * @author Guanghan Li NUAA
 * @version 1.0
 * @date 2022-10-12
 */

# ifndef __MULTIGRID__H__
# define __MULTIGRID__H__

# include <algorithm>
# include "OperatorInterpolate.h"

template <typename SINGLEGRIDSOLVER>
class NonlinearMultiGrid{
	 public:
		  typedef SINGLEGRIDSOLVER SingleGridSolver;
		  typedef typename SingleGridSolver::Residual Residual;
		  typedef typename Residual::sol_t sol_t;
		  typedef typename Residual::Problem Problem;
		  typedef typename sol_t::dis_t dis_t;
		  typedef typename sol_t::mesh_t mesh_t;
		  typedef typename dis_t::value_type value_t;
		  typedef typename dis_t::Indices index_t;
		  typedef typename dis_t::Velocity velocity_t;
		  typedef std::shared_ptr<mesh_t> ptr_mesh_t;
		  typedef hymos::twod::operator_interpolate::Transfer<sol_t> transfer_t;
	 protected:
		  int CycleType;
		  int PreSteps, PostSteps, CoarsestSteps;
		  int CoarsestGrid;
		  int FMGLevelSteps;
		  int FMGFinestSteps;
		  int FMGProlongationMode;
		  int FMGPreservedOrder;
		  double tolerance;
		  std::shared_ptr<SingleGridSolver> p_basic_method;
	 public:
		  NonlinearMultiGrid() : CycleType(1),
							PreSteps(2),
							PostSteps(2),
							CoarsestSteps(4),
							CoarsestGrid(8),
							FMGLevelSteps(1),
							FMGFinestSteps(3),
							FMGProlongationMode(0),
							FMGPreservedOrder(0),
							tolerance(1e-8),
							p_basic_method(new SingleGridSolver()){};

		  NonlinearMultiGrid(const ptr_mesh_t& p_mesh, std::istream& is){
				is.read((char *)&PreSteps, sizeof(int));
				is.read((char *)&PostSteps, sizeof(int));
				is.read((char *)&CoarsestSteps, sizeof(int));
				is.read((char *)&CoarsestGrid, sizeof(int));
				is.read((char *)&CycleType, sizeof(int));
				is.read((char *)&FMGLevelSteps, sizeof(int));
				is.read((char *)&FMGFinestSteps, sizeof(int));
				is.read((char *)&FMGProlongationMode, sizeof(int));
				is.read((char *)&FMGPreservedOrder, sizeof(int));
				is.read((char *)&tolerance, sizeof(double));
				p_basic_method.reset(new SingleGridSolver(p_mesh, is));
		  }

		  bool is_coarsestgrid(const mesh_t& mesh){
				return CoarsestGrid == (int)mesh.n_ele_x()
					&& CoarsestGrid == (int)mesh.n_ele_y();
		  }

		  void initialize(sol_t& sol, mesh_t& mesh){
				Problem::initialize(sol, mesh);
		  }

		  void read_config(const ConfigFile& cf){
				p_basic_method->read_config(cf);
				PreSteps = (int)cf.Value("MultiGridMethod", "PreSmoothingSteps");
				PostSteps = (int)cf.Value("MultiGridMethod", "PostSmoothingSteps");
				CoarsestSteps = (int)cf.Value("MultiGridMethod",
						  "CoarsestGridSmoothingSteps");
				CoarsestGrid = (int)cf.Value("MultiGridMethod", "CoarsestGridSize");
				tolerance = (double)cf.Value("Error", "Tol");
				try{
					CycleType = (int)cf.Value("MultiGridMethod", "CycleType");
				}catch(const std::string&){
					CycleType = 1;
				}
				if(CycleType < 1 || CycleType > 3){
					throw std::string("MultiGridMethod/CycleType must be 1 (V-cycle), 2 (W-cycle), or 3 (F-cycle)");
				}
				if(CoarsestGrid <= 0){
					throw std::string("MultiGridMethod/CoarsestGridSize is invalid");
				}
				int nx = (int)cf.Value("Mesh", "Nx");
				int ny = (int)cf.Value("Mesh", "Ny");
				if(nx < CoarsestGrid || ny < CoarsestGrid){
					throw std::string("Mesh/Nx, Mesh/Ny and MultiGridMethod/CoarsestGridSize are incompatible");
				}
				while(nx > CoarsestGrid || ny > CoarsestGrid){
					if(nx <= CoarsestGrid || ny <= CoarsestGrid || nx % 2 != 0 || ny % 2 != 0){
						throw std::string("Mesh/Nx and Mesh/Ny must stay even while coarsening toward MultiGridMethod/CoarsestGridSize");
					}
					nx /= 2;
					ny /= 2;
				}
				if(nx != CoarsestGrid || ny != CoarsestGrid){
					throw std::string("Mesh/Nx, Mesh/Ny and MultiGridMethod/CoarsestGridSize are incompatible");
				}
				try{
					FMGLevelSteps = std::max(0, (int)cf.Value("MultiGridMethod", "FMGLevelStartupSteps"));
				}catch(const std::string&){
					FMGLevelSteps = 1;
				}
				try{
					FMGFinestSteps = std::max(0, (int)cf.Value("MultiGridMethod", "FMGFinestStartupSteps"));
				}catch(const std::string&){
					FMGFinestSteps = 3;
				}
				try{
					FMGProlongationMode = std::max(0, (int)cf.Value("MultiGridMethod", "FMGProlongationMode"));
				}catch(const std::string&){
					FMGProlongationMode = 0;
				}
				try{
					FMGPreservedOrder = std::max(0, (int)cf.Value("MultiGridMethod", "FMGPreservedOrder"));
				}catch(const std::string&){
					FMGPreservedOrder = 0;
				}
		  }

		  void print_config(){
				std::cout << " NMG: CycleType = " << CycleType << std::endl;
				std::cout << " NMG: PreSmoothingSteps = " << PreSteps << std::endl;
				std::cout << " NMG: PostSmoothingSteps = " << PostSteps << std::endl;
				std::cout << " NMG: CoarsestSteps = " << CoarsestSteps << std::endl;
				std::cout << " NMG: CoarsestGrid = " << CoarsestGrid << "*" <<
					 CoarsestGrid << std::endl;
				std::cout << " NMG: Tol = " << tolerance << std::endl;
				std::cout << " FMG: LevelStartupSteps = " << FMGLevelSteps << std::endl;
				std::cout << " FMG: FinestStartupSteps = " << FMGFinestSteps << std::endl;
				std::cout << " FMG: ProlongationMode = " << FMGProlongationMode << std::endl;
				std::cout << " FMG: PreservedOrder = " << FMGPreservedOrder << std::endl;
				p_basic_method->print_config();
		  }

		  void dump_data(std::ostream& os) const{
				os.write((const char *)&PreSteps, sizeof(int));
				os.write((const char *)&PostSteps, sizeof(int));
				os.write((const char *)&CoarsestSteps, sizeof(int));
				os.write((const char *)&CoarsestGrid, sizeof(int));
				os.write((const char *)&CycleType, sizeof(int));
				os.write((const char *)&FMGLevelSteps, sizeof(int));
				os.write((const char *)&FMGFinestSteps, sizeof(int));
				os.write((const char *)&FMGProlongationMode, sizeof(int));
				os.write((const char *)&FMGPreservedOrder, sizeof(int));
				os.write((const char *)&tolerance, sizeof(double));
				p_basic_method->dump_data(os);
		  }

		  void load_data(std::istream& is) const{
				is.read((char *)&PreSteps, sizeof(int));
				is.read((char *)&PostSteps, sizeof(int));
				is.read((char *)&CoarsestSteps, sizeof(int));
				is.read((char *)&CoarsestGrid, sizeof(int));
				is.read((char *)&CycleType, sizeof(int));
				is.read((char *)&FMGLevelSteps, sizeof(int));
				is.read((char *)&FMGFinestSteps, sizeof(int));
				is.read((char *)&FMGProlongationMode, sizeof(int));
				is.read((char *)&FMGPreservedOrder, sizeof(int));
				is.read((char *)&tolerance, sizeof(double));
				p_basic_method->load_data(is);
		  }

		  void postProcessing( sol_t& sol, const mesh_t& mesh ) {
				Problem::postProcessing( sol, mesh );
        }

		  void coarsestSolve(sol_t& sol, const sol_t& rhs, const mesh_t& mesh){
				p_basic_method->solve(sol, rhs, mesh, CoarsestSteps);
		  }

		  void preSmoothing(sol_t& sol, const sol_t& rhs, const mesh_t& mesh){
				p_basic_method->solve(sol, rhs, mesh, PreSteps);
		  }

		  void postSmoothing(sol_t& sol, const sol_t& rhs, const mesh_t& mesh){
				p_basic_method->solve(sol, rhs, mesh, PostSteps);
		  }

		  // 回代后修复物理约束：若某单元 ρ 或 θ 变负（大修正量导致），
		  // 用当前参考系参数重建 Maxwellian，防止后续 FIM 触发 NaN。
		  void repairPhysicalConstraints(sol_t& sol, const mesh_t& mesh) {
				size_t Nx = mesh.n_ele_x();
				size_t Ny = mesh.n_ele_y();
				for (size_t ny = 0; ny < Ny; ny++) {
					for (size_t nx = 0; nx < Nx; nx++) {
						dis_t& dis = sol[nx][ny];
						double pv[dis_t::dim + 2];
						dis.PrimitiveVars(pv);
						bool ok = (pv[0] > 1e-12) && (pv[dis_t::dim + 1] > 1e-12);
						for (int k = 0; ok && k < (int)(dis_t::dim + 2); ++k)
							ok = std::isfinite(pv[k]);
						if (!ok) {
							const double rho = (std::isfinite(pv[0]) && pv[0] > 1e-12)
								? pv[0] : 1.0;
							dis.Reinit(dis.Center(), dis.Scaling(), dis.GetOrder());
							dis.MakeMaxwellian(rho);
						}
					}
				}
		  }

		  void restriction(sol_t& finer_sol, sol_t& coarser_sol,
								 const mesh_t& mesh){
				transfer_t::RestrictSolution(finer_sol, coarser_sol, mesh);
		  }

		  void restriction(sol_t& coarser_rhs, sol_t& finer_res, sol_t&
					 coarser_sol, const mesh_t& mesh){
				transfer_t::RestrictResidual(coarser_rhs, finer_res, coarser_sol, mesh);
		  }

		  void restriction(dis_t& dis_H, const dis_t& dis_h){
				dis_t tmp(dis_H.Center(), dis_H.Scaling(), dis_H.GetOrder());
				Project(dis_h, tmp);
				typename dis_t::Iterator it_dis = dis_H.begin();
            	typename dis_t::Iterator the_end = dis_H.end();
				typename dis_t::Iterator it_dis_h = tmp.begin();
				for(; it_dis != the_end; ++it_dis, ++it_dis_h){
					 *it_dis += *it_dis_h / 4;
				}
		  }

		 	  void subGridIteration(sol_t& sol, const sol_t& rhs, mesh_t& mesh,
						const int cycle_type){
					switch(cycle_type){
					case 1:
						multigridOneStep(sol, rhs, mesh, 1);
						break;
					case 2:
						multigridOneStep(sol, rhs, mesh, 2);
						multigridOneStep(sol, rhs, mesh, 2);
						break;
					case 3:
						multigridOneStep(sol, rhs, mesh, 3);
						multigridOneStep(sol, rhs, mesh, 1);
						break;
					default:
						throw std::string("MultiGridMethod/CycleType must be 1 (V-cycle), 2 (W-cycle), or 3 (F-cycle)");
					}
			  }

		  void calcCorrection(const sol_t& old_sol,
				  sol_t& new_sol,
				  const mesh_t& mesh){
				size_t Nx = mesh.n_ele_x();
				size_t Ny = mesh.n_ele_y();
				#pragma omp parallel for
				for(size_t ny = 0; ny < Ny; ny++){
					for(size_t nx = 0; nx < Nx; nx++){
					    dis_t delta_dis(old_sol[nx][ny].Center(),
					      	  old_sol[nx][ny].Scaling(),
					      	  old_sol[nx][ny].GetOrder());
					    Project(new_sol[nx][ny], delta_dis);
					    delta_dis -= old_sol[nx][ny];
					    std::swap(new_sol[nx][ny], delta_dis);
					}
				}

		  }

		  void backSubstitution(const sol_t& coarser_sol,
				  sol_t& finer_sol,
				  const mesh_t& mesh){
				transfer_t::ProlongCorrection(coarser_sol, finer_sol, mesh);
		  }

		  void backSubstitution(const dis_t& dis_H, dis_t& dis_h){
				double mu = 1.0;
				dis_h.Add(mu, dis_H);
				dis_t tmp(dis_H.Center(), dis_H.Scaling(), dis_H.GetOrder());
				tmp.ProjectToStdSpace(dis_h);
				std::swap(tmp, dis_h);
		  }

			  void multigridOneStep(sol_t& sol, const sol_t& rhs, mesh_t& mesh,
						const int cycle_type){
					if(is_coarsestgrid(mesh)){
						 coarsestSolve(sol, rhs, mesh);
						 return;
				}

				preSmoothing(sol, rhs, mesh);
				{
					sol_t finer_res;
					Residual::getResidual(sol, rhs, finer_res, mesh);

					mesh.Coarser();
					sol_t coarser_sol(mesh.n_ele_x(), mesh.n_ele_y(), sol.GetOrder());
					restriction(sol, coarser_sol, mesh);
					sol_t init_coarser_sol;
					init_coarser_sol = coarser_sol;

					sol_t coarser_res(mesh.n_ele_x(), mesh.n_ele_y(), sol.GetOrder());
					Residual::getResidual(coarser_sol, coarser_res, mesh);
					sol_t coarser_rhs(mesh.n_ele_x(), mesh.n_ele_y(), sol.GetOrder());
					restriction(coarser_rhs, finer_res, coarser_sol, mesh);
					//程序和论文相差一个负号
					coarser_rhs -= coarser_res ;

						subGridIteration(coarser_sol, coarser_rhs, mesh, cycle_type);
					calcCorrection(init_coarser_sol, coarser_sol, mesh);
					backSubstitution(coarser_sol, sol, mesh);
					mesh.Finer();
					repairPhysicalConstraints(sol, mesh);
				}
				postSmoothing(sol, rhs, mesh);
		  }

			  void solve(sol_t& sol, const sol_t& rhs, const mesh_t& mesh){
				   mesh_t Mesh = mesh;
					multigridOneStep(sol, rhs, Mesh, CycleType);
			  }

};

template <typename SINGLEGRIDSOLVER>
class FullNonlinearMultiGrid : public NonlinearMultiGrid<SINGLEGRIDSOLVER>{
	public:
		typedef NonlinearMultiGrid<SINGLEGRIDSOLVER> BASE;
		typedef typename BASE::SingleGridSolver SingleGridSolver;
		typedef typename BASE::Residual Residual;
		typedef typename BASE::sol_t sol_t;
		typedef typename BASE::Problem Problem;
		typedef typename BASE::dis_t dis_t;
		typedef typename BASE::mesh_t mesh_t;
		typedef typename BASE::value_t value_t;
		typedef typename BASE::ptr_mesh_t ptr_mesh_t;

		FullNonlinearMultiGrid() : BASE() {}

		FullNonlinearMultiGrid(const ptr_mesh_t& p_mesh, std::istream& is) : BASE(p_mesh, is) {}

		void initialize(sol_t& sol, mesh_t& mesh){
			std::cerr << "FullNonlinearMultiGrid Solver: Prepare initial value ..." << std::endl;
			std::ofstream os_res("res_step", std::ios::app);
			os_res.precision(12);
			os_res << "% FullNonlinearMultiGrid Solver: Prepare initial value ..." << std::endl;

			mesh_t current_mesh = mesh;
			while(current_mesh.n_ele_x() > (size_t)this->CoarsestGrid){
				current_mesh.Coarser();
			}

			sol_t current_sol(current_mesh.n_ele_x(), current_mesh.n_ele_y(), sol.GetOrder());
			Problem::initialize(current_sol, current_mesh);
			std::cerr << "FMG level 0 mesh = "
					  << current_mesh.n_ele_x() << "x" << current_mesh.n_ele_y() << std::endl;
			os_res << "% level 0 mesh = "
			       << current_mesh.n_ele_x() << "x" << current_mesh.n_ele_y() << std::endl;
			sol_t current_rhs;
			current_rhs.Reinit(current_sol);
			solveAtLevel(current_sol, current_rhs, current_mesh, os_res, this->CoarsestSteps);

			unsigned int level = 0;
			while(current_mesh.n_ele_x() < mesh.n_ele_x()){
				mesh_t finer_mesh = current_mesh;
				finer_mesh.Finer();
				sol_t finer_sol(finer_mesh.n_ele_x(), finer_mesh.n_ele_y(), sol.GetOrder());
				prolongation(current_sol, finer_sol, current_mesh);
				++level;
				std::cerr << "FMG prolongation to level " << level << " mesh = "
						  << finer_mesh.n_ele_x() << "x" << finer_mesh.n_ele_y() << std::endl;
				os_res << "% level " << level << " mesh = "
				       << finer_mesh.n_ele_x() << "x" << finer_mesh.n_ele_y() << std::endl;
				current_rhs.Reinit(finer_sol);
				const unsigned int startup_steps =
					(finer_mesh.n_ele_x() == mesh.n_ele_x() && finer_mesh.n_ele_y() == mesh.n_ele_y())
						? (unsigned int)this->FMGFinestSteps
						: (unsigned int)this->FMGLevelSteps;
				solveAtLevel(finer_sol, current_rhs, finer_mesh, os_res, startup_steps);
				std::swap(current_sol, finer_sol);
				current_mesh = finer_mesh;
			}

			std::swap(sol, current_sol);
			os_res << "% FullNonlinearMultiGrid Solver: Prepare initial value is OK!" << std::endl;
			os_res.close();
			std::cerr << "FullNonlinearMultiGrid Solver: Prepare initial value is OK!" << std::endl;
		}

		void print_config(){
			std::cout << " SolverType = fmg" << std::endl;
			BASE::print_config();
		}

	private:
		void solveAtLevel(sol_t& sol, const sol_t& rhs, const mesh_t& mesh, std::ostream& os_res, unsigned int max_steps){
			double res_norm = 0.0;
			Residual::getResidualPL2Norm(res_norm, sol, rhs, mesh);
			unsigned int steps = 0;
			os_res << "%   initial residual = " << res_norm << std::endl;
			sol_t safe_sol = sol;
			while(steps < max_steps){
				++steps;
				BASE::solve(sol, rhs, mesh);
				BASE::postProcessing(sol, mesh);
				Residual::getResidualPL2Norm(res_norm, sol, rhs, mesh);
				os_res << "%   step " << steps << " residual = " << res_norm << std::endl;
				if(res_norm != res_norm){
					sol = safe_sol;
					os_res << "%   NaN detected, reverted to last good state" << std::endl;
					std::cerr << "FMG level solve: NaN detected at step " << steps
							  << ", reverted to last good state" << std::endl;
					break;
				}
				safe_sol = sol;
			}
			std::cerr << "FMG level solve stopped after " << steps
					  << " steps with residual = " << res_norm << std::endl;
		}

		void prolongation(const sol_t& coarser_sol, sol_t& finer_sol, const mesh_t& mesh){
			if (this->FMGProlongationMode == 1) {
				BASE::transfer_t::ProlongInitialGuessLinear(coarser_sol, finer_sol, mesh);
			} else if (this->FMGProlongationMode == 2) {
				BASE::transfer_t::ProlongInitialGuessWithRepair(coarser_sol, finer_sol, mesh,
					(unsigned int)this->FMGPreservedOrder);
			} else {
				BASE::transfer_t::ProlongInitialGuessMaxwellian(coarser_sol, finer_sol, mesh);
			}
		}
};

# endif
