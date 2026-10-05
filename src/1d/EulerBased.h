#ifndef __IMPLICIT_H_
#define __IMPLICIT_H_

#include "EulerSolver.h"
#include "SemImplicit.h"
#include "SingleGridMethod.h"
#include <omp.h>

template <typename _RESIDUAL>
class EulerBased {
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

  typedef typename Eigen::VectorXd Euler_sol_t;
  typedef typename std::vector<Euler_sol_t> Euler_sol_vec_t;

private:
  enum FIMType { FIM1 = 1, FIM2 = 2, FIM3 = 3 };

  double tolerance;
  int fim_type;

  ExplicitTimeForward<RESIDUAL> high_order_fim1_;
  SemImpShakhov<RESIDUAL> high_order_fim2_;
  SemImpShakhovSGS<RESIDUAL> high_order_fim3_;

  EulerTimeForward<RESIDUAL> euler_fim1_;
  EulerTimeForward<RESIDUAL> euler_fim2_;
  EulerGSIteration<RESIDUAL> euler_fim3_;

public:
  EulerBased() : tolerance(1e-8), fim_type(FIM3) {}

  EulerBased(const ptr_mesh_t& p_mesh, const ConfigFile& cf) : tolerance(1e-8), fim_type(FIM3) {
    read_config(p_mesh, cf);
  }

  EulerBased(const ptr_mesh_t& p_mesh, std::istream& is) : tolerance(1e-8), fim_type(FIM3) {
    is.read((char*)&fim_type, sizeof(int));
    is.read((char*)&tolerance, sizeof(double));

    if (fim_type == FIM1) {
      high_order_fim1_ = ExplicitTimeForward<RESIDUAL>(p_mesh, is);
      euler_fim1_ = EulerTimeForward<RESIDUAL>(p_mesh, is);
    } else if (fim_type == FIM2) {
      high_order_fim2_ = SemImpShakhov<RESIDUAL>(p_mesh, is);
      euler_fim2_ = EulerTimeForward<RESIDUAL>(p_mesh, is);
    } else {
      fim_type = FIM3;
      high_order_fim3_ = SemImpShakhovSGS<RESIDUAL>(p_mesh, is);
      euler_fim3_ = EulerGSIteration<RESIDUAL>(p_mesh, is);
    }
  }

  void initial_value(sol_t& sol, const ptr_mesh_t& p_mesh, const ptr_force_t& /*p_force_handle*/) {
    PROBLEM::initial_value(sol, *p_mesh);
  }

  void ExactSolve(sol_t& sol,
                  const sol_t& rhs,
                  const ptr_mesh_t& p_mesh,
                  const ptr_force_t& p_force_handle,
                  const unsigned int max_steps = 0) {
    double initial_res;
    RESIDUAL::GetResidualPL2Norm(initial_res, sol, rhs, *p_mesh, *p_force_handle);
    const double _Tol = (initial_res > 10 * tolerance) ? tolerance : 0.1 * tolerance;

    double res;
    unsigned int step = 0;
    do {
      ++step;
      Solve(sol, rhs, p_mesh, p_force_handle, 1);
      RESIDUAL::GetResidualPL2Norm(res, sol, rhs, *p_mesh, *p_force_handle);
    } while (res >= _Tol && (step < max_steps || max_steps == 0));
  }

  void PostProcessing(sol_t& sol, const ptr_mesh_t& p_mesh) {
    PROBLEM::PostProcessing(sol, *p_mesh);
  }

  void Reinit(const ptr_mesh_t& /*p_mesh*/) {}

  void read_config(const ptr_mesh_t& p_mesh, const ConfigFile& cf) {
    fim_type = ReadFIMType(cf);
    tolerance = (double)cf.Value("Error", "Tol");

    if (fim_type == FIM1) {
      high_order_fim1_.read_config(p_mesh, cf);
      euler_fim1_.read_config(p_mesh, cf);
    } else if (fim_type == FIM2) {
      high_order_fim2_.read_config(p_mesh, cf);
      euler_fim2_.read_config(p_mesh, cf);
    } else {
      fim_type = FIM3;
      high_order_fim3_.read_config(p_mesh, cf);
      euler_fim3_.read_config(p_mesh, cf);
    }
  }

  void print_config() const {
    std::cout << "FIM mode = FIM-" << fim_type << std::endl;
    if (fim_type == FIM1) {
      std::cout << "High-order solver = ExplicitTimeForward" << std::endl;
      std::cout << "Hydrodynamic solver = EulerTimeForward" << std::endl;
      high_order_fim1_.print_config();
      euler_fim1_.print_config();
    } else if (fim_type == FIM2) {
      std::cout << "High-order solver = SemImpShakhov" << std::endl;
      std::cout << "Hydrodynamic solver = EulerTimeForward" << std::endl;
      high_order_fim2_.print_config();
      euler_fim2_.print_config();
    } else {
      std::cout << "High-order solver = SemImpShakhovSGS" << std::endl;
      std::cout << "Hydrodynamic solver = EulerGSIteration" << std::endl;
      high_order_fim3_.print_config();
      euler_fim3_.print_config();
    }
  }

  void dump_data(std::ostream& os) const {
    os.write((const char*)&fim_type, sizeof(int));
    os.write((const char*)&tolerance, sizeof(double));

    if (fim_type == FIM1) {
      high_order_fim1_.dump_data(os);
      euler_fim1_.dump_data(os);
    } else if (fim_type == FIM2) {
      high_order_fim2_.dump_data(os);
      euler_fim2_.dump_data(os);
    } else {
      high_order_fim3_.dump_data(os);
      euler_fim3_.dump_data(os);
    }
  }

  void reset_tol(const double _tol) {
    tolerance = _tol;
    if (fim_type == FIM1) {
      high_order_fim1_.reset_tol(_tol);
    } else if (fim_type == FIM2) {
      high_order_fim2_.reset_tol(_tol);
    } else {
      high_order_fim3_.reset_tol(_tol);
    }
  }

  double& reset_tol() { return tolerance; }
  const double& reset_tol() const { return tolerance; }

  void EulerToHighOrder(sol_t& sol, const Euler_sol_vec_t& Euler_sol_vec) {
#pragma omp parallel for
    for (int i = 0; i < (int)sol.size(); i++) {
      linearCombination(sol[i], Euler_sol_vec[i]);
    }
  }

  void linearCombination(dis_t& dis, const Euler_sol_t& Euler_sol) {
    value_t rho_old = dis.front();
    velocity_t v_old = dis.Center();
    value_t T_old = dis.TransTemperature();

    double rho = Euler_sol[0];
    double u1 = Euler_sol[1] / rho;
    double u2 = Euler_sol[2] / rho;
    double u3 = Euler_sol[3] / rho;
    double theta = GetTemperature(Euler_sol);

    double tau = 1.0;
    rho = tau * rho + rho_old * (1 - tau);
    u1 = tau * u1 + v_old[0] * (1 - tau);
    u2 = tau * u2 + v_old[1] * (1 - tau);
    u3 = tau * u3 + v_old[2] * (1 - tau);
    theta = tau * theta + T_old * (1 - tau);

    dis.front() = rho;
    velocity_t Center{u1, u2, u3};
    value_t Scaling = sqrt(theta);
    dis.Reinit(Center, Scaling);
  }

  void HighOrderToEuler(Euler_sol_vec_t& Euler_sol_vec, Euler_sol_vec_t& closure_vec, const sol_t& sol) {
    unsigned int n_ele = sol.size();
    Euler_sol_vec.resize(n_ele);
    closure_vec.resize(n_ele);
 #pragma omp parallel for
    for (int i = 0; i < (int)n_ele; i++) {
      GetMacroVar(sol[i], Euler_sol_vec[i], closure_vec[i]);
    }
  }

  void getRHS(Euler_sol_vec_t& rhs_vec, const sol_t& rhs) {
    unsigned int n_ele = rhs.size();
    rhs_vec.resize(n_ele);
    double pv[5];
#pragma omp parallel for private(pv)
    for (int i = 0; i < (int)n_ele; i++) {
      rhs[i].PrimitiveVars(pv);
      rhs_vec[i].resize(5);
      rhs_vec[i][0] = pv[0];
      rhs_vec[i][1] = pv[0] * pv[1];
      rhs_vec[i][2] = pv[0] * pv[2];
      rhs_vec[i][3] = pv[0] * pv[3];
      rhs_vec[i][4] = pv[0] * (3.0 / 2.0 * pv[4] + (pv[1] * pv[1] + pv[2] * pv[2] + pv[3] * pv[3]) / 2.0);
    }
  }

  void Solve(sol_t& sol,
             const sol_t& rhs,
             const ptr_mesh_t& p_mesh,
             const ptr_force_t& p_force_handle,
             const unsigned int steps = 1) {
    Euler_sol_vec_t Euler_sol_vec, closure_vec, rhs_vec;
    getRHS(rhs_vec, rhs);

    for (unsigned int i = 0; i < steps; ++i) {
      SolveHighOrder(sol, rhs, p_mesh, p_force_handle);
      HighOrderToEuler(Euler_sol_vec, closure_vec, sol);
      SolveEuler(sol, Euler_sol_vec, closure_vec, rhs_vec, *p_mesh);
      EulerToHighOrder(sol, Euler_sol_vec);
    }
  }

private:
  static int ReadFIMType(const ConfigFile& cf) {
    try {
      const int value = (int)cf.Value("Solver", "FIM");
      if (value >= FIM1 && value <= FIM3) {
        return value;
      }
      throw std::string("Solver/FIM must be 1, 2, or 3");
    } catch (const std::string& e) {
      if (e == "Solver/FIM must be 1, 2, or 3") {
        throw;
      }
    }
    return FIM3;
  }

  void SolveHighOrder(sol_t& sol, const sol_t& rhs, const ptr_mesh_t& p_mesh, const ptr_force_t& p_force_handle) {
    if (fim_type == FIM1) {
      high_order_fim1_.Solve(sol, rhs, p_mesh, p_force_handle);
    } else if (fim_type == FIM2) {
      high_order_fim2_.Solve(sol, rhs, p_mesh, p_force_handle);
    } else {
      high_order_fim3_.Solve(sol, rhs, p_mesh, p_force_handle);
    }
  }

  void SolveEuler(sol_t& sol,
                  Euler_sol_vec_t& Euler_sol_vec,
                  const Euler_sol_vec_t& closure_vec,
                  const Euler_sol_vec_t& rhs_vec,
                  const mesh_t& mesh) {
    if (fim_type == FIM1) {
      euler_fim1_.solve(sol, Euler_sol_vec, closure_vec, rhs_vec, mesh);
    } else if (fim_type == FIM2) {
      euler_fim2_.solve(sol, Euler_sol_vec, closure_vec, rhs_vec, mesh);
    } else {
      euler_fim3_.solve(sol, Euler_sol_vec, closure_vec, rhs_vec, mesh);
    }
  }

  void GetMacroVar(const dis_t& dis, Euler_sol_t& Euler_sol, Euler_sol_t& closure) {
    if (fim_type == FIM1) {
      euler_fim1_.getMacroVar(dis, Euler_sol, closure);
    } else if (fim_type == FIM2) {
      euler_fim2_.getMacroVar(dis, Euler_sol, closure);
    } else {
      euler_fim3_.getMacroVar(dis, Euler_sol, closure);
    }
  }

  double GetTemperature(const Euler_sol_t& Euler_sol) {
    if (fim_type == FIM1) {
      return euler_fim1_.getTemperature(Euler_sol);
    } else if (fim_type == FIM2) {
      return euler_fim2_.getTemperature(Euler_sol);
    } else {
      return euler_fim3_.getTemperature(Euler_sol);
    }
  }
};

#endif
