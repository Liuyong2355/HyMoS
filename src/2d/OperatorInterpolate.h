#ifndef HYMOS_SRC_2D_OPERATORINTERPOLATE_H
#define HYMOS_SRC_2D_OPERATORINTERPOLATE_H

#include <cmath>

namespace hymos {
namespace twod {
namespace operator_interpolate {

inline double Minmod(const double a, const double b) {
  if (a * b <= 0.0) {
    return 0.0;
  }
  return (fabs(a) < fabs(b)) ? a : b;
}

template <typename SOLUTION>
class Transfer {
public:
  typedef SOLUTION sol_t;
  typedef typename sol_t::dis_t dis_t;
  typedef typename sol_t::mesh_t mesh_t;
  typedef typename dis_t::Velocity velocity_t;

  static void RestrictSolution(const sol_t& finer_sol, sol_t& coarser_sol, const mesh_t& coarser_mesh) {
    size_t Nx = coarser_mesh.n_ele_x();
    size_t Ny = coarser_mesh.n_ele_y();
    #pragma omp parallel for
    for (size_t ny = 0; ny < Ny; ny++) {
      for (size_t nx = 0; nx < Nx; nx++) {
        const dis_t& dis_h_l = finer_sol[2 * nx][2 * ny];
        const dis_t& dis_h_r = finer_sol[2 * nx + 1][2 * ny];
        const dis_t& dis_h_t = finer_sol[2 * nx][2 * ny + 1];
        const dis_t& dis_h_b = finer_sol[2 * nx + 1][2 * ny + 1];
        static int dim = 3;
        const int Order = const_cast<sol_t&>(finer_sol).GetOrder();
        double pv_l[dim + 2], pv_r[dim + 2], pv_b[dim + 2], pv_t[dim + 2];
        dis_h_l.Moments(pv_l);
        dis_h_r.Moments(pv_r);
        dis_h_b.Moments(pv_b);
        dis_h_t.Moments(pv_t);

        double rho_H = (pv_l[0] + pv_r[0] + pv_b[0] + pv_t[0]) / 4;
        velocity_t u_H;
        double tmp = 0;
        for (int i = 0; i < dim; i++) {
          u_H[i] = (pv_l[i + 1] + pv_r[i + 1] + pv_b[i + 1] + pv_t[i + 1]) / 4 / rho_H;
          tmp += u_H[i] * u_H[i];
        }
        double theta_H = (pv_l[dim + 1] + pv_r[dim + 1] + pv_b[dim + 1] + pv_t[dim + 1]) / 4 / rho_H - 0.5 * tmp;
        theta_H *= 2. / dim;
        double s_H = sqrt(theta_H);
        dis_t& dis_H = coarser_sol[nx][ny];
        dis_H.Reinit(u_H, s_H, Order);

        RestrictDistribution(dis_H, dis_h_l);
        RestrictDistribution(dis_H, dis_h_r);
        RestrictDistribution(dis_H, dis_h_t);
        RestrictDistribution(dis_H, dis_h_b);
      }
    }
  }

  static void RestrictResidual(sol_t& coarser_rhs, const sol_t& finer_res, const sol_t& coarser_sol, const mesh_t& coarser_mesh) {
    size_t Nx = coarser_mesh.n_ele_x();
    size_t Ny = coarser_mesh.n_ele_y();
    #pragma omp parallel for
    for (size_t ny = 0; ny < Ny; ny++) {
      for (size_t nx = 0; nx < Nx; nx++) {
        coarser_rhs[nx][ny].Reinit(coarser_sol[nx][ny].Center(),
                                   coarser_sol[nx][ny].Scaling(),
                                   coarser_sol[nx][ny].GetOrder());
        RestrictDistribution(coarser_rhs[nx][ny], finer_res[2 * nx][2 * ny]);
        RestrictDistribution(coarser_rhs[nx][ny], finer_res[2 * nx + 1][2 * ny]);
        RestrictDistribution(coarser_rhs[nx][ny], finer_res[2 * nx][2 * ny + 1]);
        RestrictDistribution(coarser_rhs[nx][ny], finer_res[2 * nx + 1][2 * ny + 1]);
      }
    }
  }

  static void ProlongCorrection(const sol_t& coarser_sol, sol_t& finer_sol, const mesh_t& coarser_mesh) {
    size_t Nx = coarser_mesh.n_ele_x();
    size_t Ny = coarser_mesh.n_ele_y();
    #pragma omp parallel for
    for (size_t ny = 0; ny < Ny; ny++) {
      for (size_t nx = 0; nx < Nx; nx++) {
        ProlongCorrection(coarser_sol[nx][ny], finer_sol[2 * nx][2 * ny]);
        ProlongCorrection(coarser_sol[nx][ny], finer_sol[2 * nx + 1][2 * ny]);
        ProlongCorrection(coarser_sol[nx][ny], finer_sol[2 * nx][2 * ny + 1]);
        ProlongCorrection(coarser_sol[nx][ny], finer_sol[2 * nx + 1][2 * ny + 1]);
      }
    }
  }

  static void ProlongInitialGuess(const sol_t& coarser_sol, sol_t& finer_sol, const mesh_t& coarser_mesh) {
    ProlongInitialGuessWithRepair(coarser_sol, finer_sol, coarser_mesh, 0);
  }

  static void ProlongInitialGuessWithRepair(const sol_t& coarser_sol,
                                           sol_t& finer_sol,
                                           const mesh_t& coarser_mesh,
                                           unsigned int preserved_order) {
    ProlongInitialGuessLinear(coarser_sol, finer_sol, coarser_mesh);

    size_t Nx = coarser_mesh.n_ele_x();
    size_t Ny = coarser_mesh.n_ele_y();
    #pragma omp parallel for
    for (size_t ny = 0; ny < Ny; ny++) {
      for (size_t nx = 0; nx < Nx; nx++) {
        const dis_t& parent = coarser_sol[nx][ny];
        RepairInitialGuess(parent, finer_sol[2 * nx][2 * ny], preserved_order);
        RepairInitialGuess(parent, finer_sol[2 * nx + 1][2 * ny], preserved_order);
        RepairInitialGuess(parent, finer_sol[2 * nx][2 * ny + 1], preserved_order);
        RepairInitialGuess(parent, finer_sol[2 * nx + 1][2 * ny + 1], preserved_order);
      }
    }
  }

  static void ProlongInitialGuessMaxwellian(const sol_t& coarser_sol,
                                            sol_t& finer_sol,
                                            const mesh_t& coarser_mesh) {
    size_t Nx = coarser_mesh.n_ele_x();
    size_t Ny = coarser_mesh.n_ele_y();
    #pragma omp parallel for
    for (size_t ny = 0; ny < Ny; ny++) {
      for (size_t nx = 0; nx < Nx; nx++) {
        const dis_t& parent = coarser_sol[nx][ny];
        MakeMaxwellianFrom(parent, finer_sol[2 * nx][2 * ny]);
        MakeMaxwellianFrom(parent, finer_sol[2 * nx + 1][2 * ny]);
        MakeMaxwellianFrom(parent, finer_sol[2 * nx][2 * ny + 1]);
        MakeMaxwellianFrom(parent, finer_sol[2 * nx + 1][2 * ny + 1]);
      }
    }
  }

  static void ProlongInitialGuessLinear(const sol_t& coarser_sol, sol_t& finer_sol, const mesh_t& coarser_mesh) {
    size_t Nx = coarser_mesh.n_ele_x();
    size_t Ny = coarser_mesh.n_ele_y();
    #pragma omp parallel for
    for (size_t ny = 0; ny < Ny; ny++) {
      for (size_t nx = 0; nx < Nx; nx++) {
        const dis_t& dis = coarser_sol[nx][ny];
        const dis_t& dis_l = (nx > 0) ? coarser_sol[nx - 1][ny] : dis;
        const dis_t& dis_r = (nx + 1 < Nx) ? coarser_sol[nx + 1][ny] : dis;
        const dis_t& dis_b = (ny > 0) ? coarser_sol[nx][ny - 1] : dis;
        const dis_t& dis_top = (ny + 1 < Ny) ? coarser_sol[nx][ny + 1] : dis;

        dis_t pro_dis_l(dis.Center(), dis.Scaling(), dis.GetOrder());
        dis_t pro_dis_r(dis.Center(), dis.Scaling(), dis.GetOrder());
        dis_t pro_dis_b(dis.Center(), dis.Scaling(), dis.GetOrder());
        dis_t pro_dis_t(dis.Center(), dis.Scaling(), dis.GetOrder());
        Project(dis_l, pro_dis_l);
        Project(dis_r, pro_dis_r);
        Project(dis_b, pro_dis_b);
        Project(dis_top, pro_dis_t);

        dis_t recon_lb(dis.Center(), dis.Scaling(), dis.GetOrder());
        dis_t recon_rb(dis.Center(), dis.Scaling(), dis.GetOrder());
        dis_t recon_lt(dis.Center(), dis.Scaling(), dis.GetOrder());
        dis_t recon_rt(dis.Center(), dis.Scaling(), dis.GetOrder());

        const double dx_l = (nx > 0) ? 0.5 * (coarser_mesh.dx(nx - 1, ny) + coarser_mesh.dx(nx, ny)) : coarser_mesh.dx(nx, ny);
        const double dx_r = (nx + 1 < Nx) ? 0.5 * (coarser_mesh.dx(nx + 1, ny) + coarser_mesh.dx(nx, ny)) : coarser_mesh.dx(nx, ny);
        const double dy_b = (ny > 0) ? 0.5 * (coarser_mesh.dy(nx, ny - 1) + coarser_mesh.dy(nx, ny)) : coarser_mesh.dy(nx, ny);
        const double dy_t = (ny + 1 < Ny) ? 0.5 * (coarser_mesh.dy(nx, ny + 1) + coarser_mesh.dy(nx, ny)) : coarser_mesh.dy(nx, ny);
        const double half_child_dx = 0.25 * coarser_mesh.dx(nx, ny);
        const double half_child_dy = 0.25 * coarser_mesh.dy(nx, ny);

        typename dis_t::ConstIterator it_dis = dis.begin();
        typename dis_t::ConstIterator it_end = dis.end();
        typename dis_t::Iterator it_pro_l = pro_dis_l.begin();
        typename dis_t::Iterator it_pro_r = pro_dis_r.begin();
        typename dis_t::Iterator it_pro_b = pro_dis_b.begin();
        typename dis_t::Iterator it_pro_t = pro_dis_t.begin();
        typename dis_t::Iterator it_lb = recon_lb.begin();
        typename dis_t::Iterator it_rb = recon_rb.begin();
        typename dis_t::Iterator it_lt = recon_lt.begin();
        typename dis_t::Iterator it_rt = recon_rt.begin();

        for (; it_dis != it_end;
             ++it_dis, ++it_pro_l, ++it_pro_r, ++it_pro_b, ++it_pro_t,
             ++it_lb, ++it_rb, ++it_lt, ++it_rt) {
          const double slope_x_l = (*it_dis - *it_pro_l) / dx_l;
          const double slope_x_r = (*it_pro_r - *it_dis) / dx_r;
          const double slope_y_b = (*it_dis - *it_pro_b) / dy_b;
          const double slope_y_t = (*it_pro_t - *it_dis) / dy_t;

          const double slope_x = Minmod(slope_x_l, slope_x_r);
          const double slope_y = Minmod(slope_y_b, slope_y_t);

          *it_lb = *it_dis - slope_x * half_child_dx - slope_y * half_child_dy;
          *it_rb = *it_dis + slope_x * half_child_dx - slope_y * half_child_dy;
          *it_lt = *it_dis - slope_x * half_child_dx + slope_y * half_child_dy;
          *it_rt = *it_dis + slope_x * half_child_dx + slope_y * half_child_dy;
        }

        dis_t& new_dis_lb = finer_sol[2 * nx][2 * ny];
        dis_t& new_dis_rb = finer_sol[2 * nx + 1][2 * ny];
        dis_t& new_dis_lt = finer_sol[2 * nx][2 * ny + 1];
        dis_t& new_dis_rt = finer_sol[2 * nx + 1][2 * ny + 1];
        new_dis_lb.Reinit(dis.Center(), dis.Scaling(), dis.GetOrder());
        new_dis_rb.Reinit(dis.Center(), dis.Scaling(), dis.GetOrder());
        new_dis_lt.Reinit(dis.Center(), dis.Scaling(), dis.GetOrder());
        new_dis_rt.Reinit(dis.Center(), dis.Scaling(), dis.GetOrder());
        new_dis_lb.ProjectToStdSpace(recon_lb);
        new_dis_rb.ProjectToStdSpace(recon_rb);
        new_dis_lt.ProjectToStdSpace(recon_lt);
        new_dis_rt.ProjectToStdSpace(recon_rt);
      }
    }
  }

private:
  static bool PrimitiveVarsAreAdmissible(const dis_t& dis) {
    double pv[dis_t::dim + 2];
    dis.PrimitiveVars(pv);
    if (!(pv[0] > 1e-12) || !(pv[dis_t::dim + 1] > 1e-12)) {
      return false;
    }
    for (int i = 0; i < dis_t::dim + 2; ++i) {
      if (!std::isfinite(pv[i])) {
        return false;
      }
    }
    return true;
  }

  static void MakeMaxwellianFrom(const dis_t& src, dis_t& dst) {
    double pv[dis_t::dim + 2];
    src.PrimitiveVars(pv);
    typename dis_t::Velocity center(&pv[1]);
    const double rho = std::max(pv[0], 1e-12);
    const double theta = std::max(pv[dis_t::dim + 1], 1e-12);
    dst.Reinit(center, sqrt(theta), src.GetOrder());
    dst.MakeMaxwellian(rho);
  }

  static void RepairInitialGuess(const dis_t& parent, dis_t& child, unsigned int preserved_order) {
    if (!PrimitiveVarsAreAdmissible(child)) {
      MakeMaxwellianFrom(parent, child);
      return;
    }

    dis_t filtered(child.Center(), child.Scaling(), child.GetOrder());
    filtered.ProjectToStdSpace(child, preserved_order);
    if (!PrimitiveVarsAreAdmissible(filtered)) {
      MakeMaxwellianFrom(parent, child);
      return;
    }
    std::swap(filtered, child);
  }

  static void RestrictDistribution(dis_t& dis_H, const dis_t& dis_h) {
    dis_t tmp(dis_H.Center(), dis_H.Scaling(), dis_H.GetOrder());
    Project(dis_h, tmp);
    typename dis_t::Iterator it_dis = dis_H.begin();
    typename dis_t::Iterator the_end = dis_H.end();
    typename dis_t::Iterator it_dis_h = tmp.begin();
    for (; it_dis != the_end; ++it_dis, ++it_dis_h) {
      *it_dis += *it_dis_h / 4;
    }
  }

  static void ProlongCorrection(const dis_t& dis_H, dis_t& dis_h) {
    double mu = 1.0;
    dis_h.Add(mu, dis_H);
    dis_t tmp(dis_H.Center(), dis_H.Scaling(), dis_H.GetOrder());
    tmp.ProjectToStdSpace(dis_h);
    std::swap(tmp, dis_h);
  }

};

}  // namespace operator_interpolate
}  // namespace twod
}  // namespace hymos

#endif
