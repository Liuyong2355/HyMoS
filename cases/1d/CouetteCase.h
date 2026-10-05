// Usage:
// - This file defines the Couette problem setup for the FIM backend.
// - Register it in `src/engine/CaseRegistry.cpp` to expose it to CLI or GUI.

#ifndef HYMOS_CASES_1D_COUETTECASE_H
#define HYMOS_CASES_1D_COUETTECASE_H

#define MAX_ORDER 50
#define MAX_VECTOR_DIMENSION 3
#define SINGLE_GRID_SOLVER

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
namespace couette {

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
  static double LeftWallTemperature;
  static double LeftWallU2;
  static double RightWallTemperature;
  static double RightWallU2;

  static const OrderReductionStrategy ors = ORS1;

  static double OptionalValue(const ConfigFile& cf, const char* section,
                              const char* key, double fallback) {
    try {
      return (double)cf.Value(section, key);
    } catch (const std::string&) {
      return fallback;
    }
  }

  static void Reinit(const ptr_mesh_t& /*p_mesh*/) {}

  static double OneOverTau(const dis_t& dis, const force_t& /*force_handle*/, const unsigned int /*i_ele*/) {
    static const double w = 0.81;
    return dis.front() * pow(dis.Scaling(), 2 - 2 * w) / Kn * sqrt(M_PI / 2);
  }

  static void initial_value(sol_t& sol, mesh_t& mesh) {
    for (unsigned int i = 0; i < sol.size(); ++i) {
      initial_value(sol[i], mesh.Center(i));
    }
  }

  static void initial_value(dis_t& dis, const double /*x*/) {
    static double u = 0.0;
    static double s = 1.0;
    static double rho = 1.0;
    dis.Reinit(velocity_t(0, u, 0), s);
    dis.MakeMaxwellian(rho);
  }

  static void BoundaryValue(dis_t& bv, const dis_t& dis, const bool& flag) {
    static double chi = 1.0;

    if (flag) {
      static std::vector<double> normal(1, 1.0);
      Reflection::BoundaryValue(dis, velocity_t(0, RightWallU2, 0),
                                RightWallTemperature, normal, chi, bv);
    } else {
      static std::vector<double> normal(1, -1.0);
      Reflection::BoundaryValue(dis, velocity_t(0, LeftWallU2, 0),
                                LeftWallTemperature, normal, chi, bv);
    }
  }

  static void PostProcessing(sol_t& sol, const mesh_t& mesh) {
    double total_mass = 0.0;
    for (size_t i = 0; i < sol.size(); ++i) {
      total_mass += sol[i].front() * mesh.Length(i);
    }
    sol *= mesh.DomainSize() / total_mass;
  }

  static void read_config(const ptr_mesh_t& p_mesh, const ConfigFile& cf) {
    Temperature = (double)cf.Value("System", "Temperature");
    Kn = (double)cf.Value("Const", "Kn");
    Pr = (double)cf.Value("Const", "Pr");
    Ma = OptionalValue(cf, "Const", "Ma", 1.0);
    LeftWallTemperature = OptionalValue(cf, "LeftWall", "Temperature", 1.0);
    LeftWallU2 = OptionalValue(cf, "LeftWall", "TangentialVelocity", -0.62885);
    RightWallTemperature = OptionalValue(cf, "RightWall", "Temperature", 1.0);
    RightWallU2 = OptionalValue(cf, "RightWall", "TangentialVelocity", 0.62885);
    Reinit(p_mesh);
  }

  static void print_config() {
    std::cout << "Couette flow: ref Temperature = " << Temperature << std::endl;
    std::cout << "Couette flow: Knudsen number = " << Kn << std::endl;
    std::cout << "Couette flow: Prandtl number = " << Pr << std::endl;
    std::cout << "Couette flow: Mach number = " << Ma << std::endl;
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
double Case::LeftWallTemperature = 1.0;
double Case::LeftWallU2 = -0.62885;
double Case::RightWallTemperature = 1.0;
double Case::RightWallU2 = 0.62885;

}  // namespace couette
}  // namespace oned
}  // namespace cases
}  // namespace hymos

#endif
