// Usage:
// - This file defines the cavity problem setup for the FIM backend.
// - Register it in `src/engine/CaseRegistry.cpp` to expose it to CLI or GUI.

#ifndef HYMOS_CASES_2D_CAVITYCASE_H
#define HYMOS_CASES_2D_CAVITYCASE_H

#include <fstream>
#include <memory>
#include <vector>

#include <NRxx/BoundaryCondition.h>
#include <NRxx/CoefTable.h>
#include <NRxx/Distribution.h>
#include <NRxx/NumericalFlux.h>
#include <NRxx/Transform.h>

#include "src/common/ConfigFile.h"
#include "src/2d/ConfiguredSolver.h"
#include "src/2d/Collision.h"
#include "src/2d/Mesh2D.h"
#include "src/2d/MultiGrid.h"
#include "src/2d/Residual.h"
#include "src/2d/SingleGridMethod.h"
#include "src/2d/Solution.h"
#include "src/2d/SteadyState.h"

namespace hymos {
namespace cases {
namespace twod {
namespace cavity {

class Case;

typedef NRXX::Distribution<3, double, 0> dis_t;
typedef Mesh2D<double> mesh_t;
typedef Solution<dis_t, mesh_t> sol_t;
typedef Residual<sol_t, Shakhov, Case> residual_t;
typedef EulerBased<residual_t> single_grid_solver_t;
typedef NonlinearMultiGrid<single_grid_solver_t> nmg_solver_t;
typedef FullNonlinearMultiGrid<single_grid_solver_t> fmg_solver_t;
typedef ConfiguredSolver2D<single_grid_solver_t, nmg_solver_t, fmg_solver_t> solver_t;

typedef SteadyState<solver_t> app_t;

class Case {
public:
  typedef dis_t::Velocity velocity_t;
  typedef std::shared_ptr<mesh_t> ptr_mesh_t;

  static double Temperature;
  static double Kn;
  static double Pr;
  static double Ma;
  static double chi;
  static const double u_l;
  static const double u_r;
  static const double u_t;
  static const double u_b;
  static const double T_l;
  static const double T_r;
  static const double T_t;
  static const double T_b;

  static void read_config(const ptr_mesh_t& /*p_mesh*/, const ConfigFile& cf) {
    Kn = (double)cf.Value("Const", "Kn");
    Pr = (double)cf.Value("Const", "Pr");
  }

  static void print_config() {
    std::cerr << " Kn = " << Kn << std::endl;
    std::cerr << " Pr = " << Pr << std::endl;
  }

  static void dump_data(std::ostream& os) {
    os.write((const char*)&Temperature, sizeof(double));
    os.write((const char*)&Kn, sizeof(double));
    os.write((const char*)&Pr, sizeof(double));
    os.write((const char*)&Ma, sizeof(double));
  }

  static void load_data(std::istream& is) {
    is.read((char*)&Temperature, sizeof(double));
    is.read((char*)&Kn, sizeof(double));
    is.read((char*)&Pr, sizeof(double));
    is.read((char*)&Ma, sizeof(double));
  }

  static void initialize(sol_t& sol, mesh_t& mesh) {
    size_t Nx = mesh.n_ele_x();
    size_t Ny = mesh.n_ele_y();
    static double u = 0.0;
    static double s = 1.0;
    static double rho = 1.0;
    for (size_t ny = 0; ny < Ny; ny++) {
      for (size_t nx = 0; nx < Nx; nx++) {
        sol[nx][ny].Reinit(velocity_t(0, u, 0), s);
        sol[nx][ny].MakeMaxwellian(rho);
      }
    }
  }

  static double getTau(const dis_t& dis) {
    const double rho = dis.front();
    const double theta = dis.TransTemperature();
    const double omega = 0.81;
    return sqrt(2 / M_PI) * Kn * pow(theta, omega - 1) / rho;
  }

  static void postProcessing(sol_t& sol, const mesh_t& mesh) {
    size_t Nx = mesh.n_ele_x();
    size_t Ny = mesh.n_ele_y();
    double total_mass = 0.0;
    for (size_t ny = 0; ny < Ny; ny++) {
      for (size_t nx = 0; nx < Nx; nx++) {
        total_mass += sol[nx][ny].front() * mesh.Area(nx, ny);
      }
    }
    sol *= mesh.DomainSize() / total_mass;
  }

  static void leftBoundaryCondition(const dis_t& dis, dis_t& ghost_dis) {
    std::vector<double> normal{-1, 0};
    velocity_t u(0, u_l, 0);
    Reflection::BoundaryValue(dis, u, T_l, normal, chi, ghost_dis);
  }

  static void rightBoundaryCondition(const dis_t& dis, dis_t& ghost_dis) {
    std::vector<double> normal{1, 0};
    velocity_t u(0, u_r, 0);
    Reflection::BoundaryValue(dis, u, T_r, normal, chi, ghost_dis);
  }

  static void topBoundaryCondition(const dis_t& dis, dis_t& ghost_dis) {
    std::vector<double> normal{0, 1};
    velocity_t u(u_t, 0.0, 0.0);
    Reflection::BoundaryValue(dis, u, T_t, normal, chi, ghost_dis);
  }

  static void bottomBoundaryCondition(const dis_t& dis, dis_t& ghost_dis) {
    std::vector<double> normal{0, -1};
    velocity_t u(u_b, 0.0, 0.0);
    Reflection::BoundaryValue(dis, u, T_b, normal, chi, ghost_dis);
  }

  static void outputSolution(const sol_t& sol, const mesh_t& mesh) {
    size_t Nx = mesh.n_ele_x(), Ny = mesh.n_ele_y();
    std::ofstream os("data.dat");
    os.precision(12);
    os << " variables= \"x\", \"y\", \"density\", "
       << " \"u\", \"v\", \"w\", "
       << " \"T\", \"tauxx\", \"tauxy\", \"tauyy\", "
       << " \"qx\",\"qy\" " << std::endl;
    os << "zone i = " << Nx << ",j=" << Ny << ", datapacking = \"point\" " << std::endl;
    std::vector<double> center(2);
    for (size_t ny = 0; ny < Ny; ny++) {
      for (size_t nx = 0; nx < Nx; nx++) {
        mesh.center(center, nx, ny);
        os << center[0] << '\t' << center[1] << '\t';
        const dis_t& dis = sol[nx][ny];
        dis_t::Velocity v = dis.Center();
        os << dis.front() << '\t';
        os << v[0] << '\t' << v[1] << '\t' << v[2] << '\t';
        os << dis.TransTemperature() << '\t';
        dis_t::Indices ind_sigma_xx(2, 0, 0);
        dis_t::Indices ind_sigma_xy(1, 1, 0);
        dis_t::Indices ind_sigma_yy(0, 2, 0);
        os << dis(ind_sigma_xx, 0) << '\t' << dis(ind_sigma_xy, 0) << '\t' << dis(ind_sigma_yy, 0) << '\t';
        double heat_flux[3];
        dis.HeatFlux(heat_flux);
        os << heat_flux[0] << '\t' << heat_flux[1] << std::endl;
      }
    }
  }
};

double Case::Temperature = 1.0;
double Case::Kn = 0.1;
double Case::Pr = 2.0 / 3.0;
double Case::Ma = 1.0;
double Case::chi = 1.0;
const double Case::u_l = 0;
const double Case::u_r = 0;
const double Case::u_b = 0;
const double Case::u_t = 50 / sqrt(273 * 1.38e-23 / 6.63e-26);
const double Case::T_l = 1.0;
const double Case::T_r = 1.0;
const double Case::T_t = 1.0;
const double Case::T_b = 1.0;

}  // namespace cavity
}  // namespace twod
}  // namespace cases
}  // namespace hymos

#endif
