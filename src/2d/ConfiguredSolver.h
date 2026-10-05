#ifndef HYMOS_SRC_2D_CONFIGUREDSOLVER_H
#define HYMOS_SRC_2D_CONFIGUREDSOLVER_H

#include <memory>

#include "src/common/SolverType.h"

template <typename SINGLE_SOLVER, typename NMG_SOLVER, typename FMG_SOLVER>
class ConfiguredSolver2D {
public:
  typedef typename SINGLE_SOLVER::Residual Residual;
  typedef typename Residual::sol_t sol_t;
  typedef typename sol_t::mesh_t mesh_t;
  typedef std::shared_ptr<mesh_t> ptr_mesh_t;

  ConfiguredSolver2D() : solver_type_(hymos::common::kSolverSingle) {
    single_solver_.reset(new SINGLE_SOLVER());
  }

  ConfiguredSolver2D(const ptr_mesh_t& p_mesh, std::istream& is)
      : solver_type_(hymos::common::kSolverSingle) {
    int stored_type = 0;
    is.read((char*)&stored_type, sizeof(int));
    if (stored_type == (int)hymos::common::kSolverNMG) {
      solver_type_ = hymos::common::kSolverNMG;
      nmg_solver_.reset(new NMG_SOLVER(p_mesh, is));
    } else if (stored_type == (int)hymos::common::kSolverFMG) {
      solver_type_ = hymos::common::kSolverFMG;
      fmg_solver_.reset(new FMG_SOLVER(p_mesh, is));
    } else {
      solver_type_ = hymos::common::kSolverSingle;
      single_solver_.reset(new SINGLE_SOLVER(p_mesh, is));
    }
  }

  void initialize(sol_t& sol, mesh_t& mesh) {
    if (solver_type_ == hymos::common::kSolverNMG) {
      EnsureNMG()->initialize(sol, mesh);
    } else if (solver_type_ == hymos::common::kSolverFMG) {
      EnsureFMG()->initialize(sol, mesh);
    } else {
      EnsureSingle()->initialize(sol, mesh);
    }
  }

  void read_config(const ConfigFile& cf) {
    solver_type_ = hymos::common::ParseRuntimeSolverType(cf, true);
    if (solver_type_ == hymos::common::kSolverNMG) {
      EnsureNMG()->read_config(cf);
    } else if (solver_type_ == hymos::common::kSolverFMG) {
      EnsureFMG()->read_config(cf);
    } else {
      solver_type_ = hymos::common::kSolverSingle;
      EnsureSingle()->read_config(cf);
    }
  }

  void dump_data(std::ostream& os) const {
    const int stored_type = (int)solver_type_;
    os.write((const char*)&stored_type, sizeof(int));
    if (solver_type_ == hymos::common::kSolverNMG) {
      nmg_solver_->dump_data(os);
    } else if (solver_type_ == hymos::common::kSolverFMG) {
      fmg_solver_->dump_data(os);
    } else {
      single_solver_->dump_data(os);
    }
  }

  void print_config() {
    std::cout << "SolverType = " << hymos::common::RuntimeSolverTypeName(solver_type_) << std::endl;
    if (solver_type_ == hymos::common::kSolverNMG) {
      EnsureNMG()->print_config();
    } else if (solver_type_ == hymos::common::kSolverFMG) {
      EnsureFMG()->print_config();
    } else {
      EnsureSingle()->print_config();
    }
  }

  void postProcessing(sol_t& sol, const mesh_t& mesh) {
    if (solver_type_ == hymos::common::kSolverNMG) {
      EnsureNMG()->postProcessing(sol, mesh);
    } else if (solver_type_ == hymos::common::kSolverFMG) {
      EnsureFMG()->postProcessing(sol, mesh);
    } else {
      EnsureSingle()->postProcessing(sol, mesh);
    }
  }

  void solve(sol_t& sol, const sol_t& rhs, const mesh_t& mesh) {
    if (solver_type_ == hymos::common::kSolverNMG) {
      EnsureNMG()->solve(sol, rhs, mesh);
    } else if (solver_type_ == hymos::common::kSolverFMG) {
      EnsureFMG()->solve(sol, rhs, mesh);
    } else {
      EnsureSingle()->solve(sol, rhs, mesh);
    }
  }

private:
  SINGLE_SOLVER* EnsureSingle() {
    if (single_solver_.get() == NULL) {
      single_solver_.reset(new SINGLE_SOLVER());
    }
    return single_solver_.get();
  }

  NMG_SOLVER* EnsureNMG() {
    if (nmg_solver_.get() == NULL) {
      nmg_solver_.reset(new NMG_SOLVER());
    }
    return nmg_solver_.get();
  }

  FMG_SOLVER* EnsureFMG() {
    if (fmg_solver_.get() == NULL) {
      fmg_solver_.reset(new FMG_SOLVER());
    }
    return fmg_solver_.get();
  }

  hymos::common::RuntimeSolverType solver_type_;
  std::unique_ptr<SINGLE_SOLVER> single_solver_;
  std::unique_ptr<NMG_SOLVER> nmg_solver_;
  std::unique_ptr<FMG_SOLVER> fmg_solver_;
};

#endif
