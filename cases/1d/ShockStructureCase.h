// Usage:
// - This file defines the shock-structure problem setup for the FIM backend.
// - Register it in `src/engine/CaseRegistry.cpp` to expose it to CLI or GUI.

#ifndef HYMOS_CASES_1D_SHOCKSTRUCTURECASE_H
#define HYMOS_CASES_1D_SHOCKSTRUCTURECASE_H

#define MAX_ORDER 50
#define MAX_VECTOR_DIMENSION 3
#include <fstream>
#include <memory>
#include <sstream>
#include <vector>

#include <NRxx/BoundaryCondition.h>
#include <NRxx/Distribution.h>

#include "src/common/ConfigFile.h"
#include "src/1d/ConfiguredSolver.h"
#include "src/1d/EulerBased.h"
#include "src/1d/ExternalForce.h"
#include "src/1d/MultiGrid.h"
#include "src/1d/Residual.h"
#include "src/1d/Solution.h"
#include "src/1d/SteadyState.h"
#include "src/1d/Mesh1D.h"

namespace hymos {
namespace cases {
namespace oned {
namespace shock_structure {

class Case;

typedef Mesh1D<double> mesh_t;
typedef NRXX::Distribution<MAX_VECTOR_DIMENSION, double, 0> dis_t;
typedef Solution<dis_t> sol_t;
typedef NoForce<mesh_t, sol_t, Case> force_t;
typedef Residual<mesh_t, sol_t, HLL, force_t, SHAKHOV, Case> residual_t;
typedef EulerBased<residual_t> single_grid_solver_t;
typedef NonlinearMultiGrid<single_grid_solver_t> nmg_solver_t;
typedef FullNonlinearMultiGrid<single_grid_solver_t> fmg_solver_t;
typedef ConfiguredSolver1D<single_grid_solver_t, nmg_solver_t, fmg_solver_t> solver_t;

typedef SteadyState<mesh_t, sol_t, force_t, solver_t, Case> app_t;

class Case {
public:
  typedef dis_t::Velocity velocity_t;
  typedef std::shared_ptr<mesh_t> ptr_mesh_t;

  static double Temperature;
  static double Kn;
  static double Pr;
  static double Ma;

  static void Reinit(const ptr_mesh_t& /*p_mesh*/) {}

  static double OneOverTau(const dis_t& dis, const force_t& /*force_handle*/, const unsigned int /*i_ele*/) {
    static const double w = 0.72;
    return dis.front() * pow(dis.Scaling(), 2 - 2 * w) * (5 - 2 * w) * (7 - 2 * w) / 15 / Kn * sqrt(2 / M_PI);
  }

  static void initial_value(sol_t& sol, mesh_t& mesh) {
    for (unsigned int i = 0; i < sol.size(); ++i) {
      initial_value(sol[i], mesh.Center(i));
    }
  }

  static void initial_value(dis_t& dis, const double x) {
    if (x < 0) {
      BoundaryValue(dis, dis, false);
    } else {
      BoundaryValue(dis, dis, true);
    }
  }

  static void BoundaryValue(dis_t& bv, const dis_t& dis, const bool& flag) {
    static double rho_l = 1.0;
    static double u_w_l = sqrt(5.0 / 3.0) * Ma;
    static double theta_l = 1.0;

    static double rho_r = 4 * rho_l * Ma * Ma / (Ma * Ma + 3);
    static double u_w_r = sqrt(5.0 / 3.0) * (Ma * Ma + 3) / 4 / Ma;
    static double theta_r = (5 * Ma * Ma - 1) * (Ma * Ma + 3) / 16 / Ma / Ma;

    if (flag) {
      bv.Reinit(velocity_t(u_w_r, 0, 0), sqrt(theta_r), dis.GetOrder());
      bv.MakeMaxwellian(rho_r);
    } else {
      bv.Reinit(velocity_t(u_w_l, 0, 0), sqrt(theta_l), dis.GetOrder());
      bv.MakeMaxwellian(rho_l);
    }
  }

  static void PostProcessing(sol_t& /*sol*/, const mesh_t& /*mesh*/) {}

  static void read_config(const ptr_mesh_t& p_mesh, const ConfigFile& cf) {
    Temperature = (double)cf.Value("System", "Temperature");
    Kn = (double)cf.Value("Const", "Kn");
    Pr = (double)cf.Value("Const", "Pr");
    Ma = (double)cf.Value("Const", "Ma");
    Reinit(p_mesh);
  }

  static void print_config() {
    std::cout << "ShockStructure: ref Temperature = " << Temperature << std::endl;
    std::cout << "ShockStructure: Knudsen number = " << Kn << std::endl;
    std::cout << "ShockStructure: Prandtl number = " << Pr << std::endl;
    std::cout << "ShockStructure: Mach number = " << Ma << std::endl;
  }

  static void load_data(const ptr_mesh_t& p_mesh, std::istream& is) {
    is.read((char*)&Temperature, sizeof(double));
    is.read((char*)&Kn, sizeof(double));
    is.read((char*)&Pr, sizeof(double));
    is.read((char*)&Ma, sizeof(double));
    Reinit(p_mesh);
  }

  static void dump_data(std::ostream& os) {
    os.write((const char*)&Temperature, sizeof(double));
    os.write((const char*)&Kn, sizeof(double));
    os.write((const char*)&Pr, sizeof(double));
    os.write((const char*)&Ma, sizeof(double));
  }

  static void save_one_over_tau(const sol_t& sol, const mesh_t& mesh, const force_t& force_handle, const unsigned int n_run) {
    std::stringstream file_tau;
    file_tau << "OOTau" << n_run << "th.dat";
    std::ofstream os_ootau(file_tau.str().c_str(), std::ios::out);
    os_ootau.precision(15);

    for (unsigned int i = 0; i < sol.size(); ++i) {
      os_ootau << mesh.Center(i) << " " << OneOverTau(sol[i], force_handle, i) << std::endl;
    }
  }

  static void save_current(const sol_t& /*sol*/, const mesh_t& /*mesh*/, const unsigned int /*n_run*/) {}

  static void NumericalFlux(const dis_t& dis0, const dis_t& dis1, const std::vector<double>& normal, dis_t& flux_l2r, dis_t& flux_r2l) {
    HLLFlux(dis0, dis1, normal, flux_l2r, flux_r2l);
  }

  static void GetDeltaFlux(const dis_t& dis0, const dis_t& dis1, const std::vector<double>& normal, const double df, const typename dis_t::Indices& ind, dis_t& dflux_l2r) {
    DeltaHLLFlux(dis0, dis1, normal, df, ind, dflux_l2r);
  }
};

double Case::Temperature = 1.0;
double Case::Kn = 0.0;
double Case::Pr = 2.0 / 3.0;
double Case::Ma = 1.0;

}  // namespace shock_structure
}  // namespace oned
}  // namespace cases
}  // namespace hymos

#endif
