#ifndef HYMOS_SRC_1D_CONFIGUREDSOLVER_H
#define HYMOS_SRC_1D_CONFIGUREDSOLVER_H

#include <memory>

#include "src/common/SolverType.h"

template <typename SINGLE_SOLVER, typename NMG_SOLVER, typename FMG_SOLVER>
class ConfiguredSolver1D {
public:
  typedef typename SINGLE_SOLVER::RESIDUAL RESIDUAL;
  typedef typename RESIDUAL::mesh_t mesh_t;
  typedef typename RESIDUAL::sol_t sol_t;
  typedef typename RESIDUAL::force_t force_t;
  typedef std::shared_ptr<mesh_t> ptr_mesh_t;
  typedef std::shared_ptr<force_t> ptr_force_t;

  ConfiguredSolver1D() : solver_type_(hymos::common::kSolverSingle) {
    single_solver_.reset(new SINGLE_SOLVER());
  }

  ConfiguredSolver1D(const ptr_mesh_t& p_mesh, const ConfigFile& cf)
      : solver_type_(hymos::common::kSolverSingle) {
    read_config(p_mesh, cf);
  }

  ConfiguredSolver1D(const ptr_mesh_t& p_mesh, std::istream& is)
      : solver_type_(hymos::common::kSolverSingle) {
    int stored_type = 0;
    is.read((char*)&stored_type, sizeof(int));
    solver_type_ = NormalizeStoredType(stored_type);
    if (solver_type_ == hymos::common::kSolverNMG) {
      nmg_solver_.reset(new NMG_SOLVER(p_mesh, is));
    } else if (solver_type_ == hymos::common::kSolverFMG) {
      fmg_solver_.reset(new FMG_SOLVER(p_mesh, is));
    } else {
      solver_type_ = hymos::common::kSolverSingle;
      single_solver_.reset(new SINGLE_SOLVER(p_mesh, is));
    }
  }

  void Reinit(const ptr_mesh_t& p_mesh) {
    if (solver_type_ == hymos::common::kSolverNMG) {
      EnsureNMG()->Reinit(p_mesh);
    } else if (solver_type_ == hymos::common::kSolverFMG) {
      EnsureFMG()->Reinit(p_mesh);
    } else {
      EnsureSingle()->Reinit(p_mesh);
    }
  }

  void read_config(const ptr_mesh_t& p_mesh, const ConfigFile& cf) {
    solver_type_ = hymos::common::ParseRuntimeSolverType(cf, true);
    if (solver_type_ == hymos::common::kSolverNMG) {
      EnsureNMG()->read_config(p_mesh, cf);
    } else if (solver_type_ == hymos::common::kSolverFMG) {
      EnsureFMG()->read_config(p_mesh, cf);
    } else {
      solver_type_ = hymos::common::kSolverSingle;
      EnsureSingle()->read_config(p_mesh, cf);
    }
  }

  void print_config() const {
    std::cout << "SolverType = " << hymos::common::RuntimeSolverTypeName(solver_type_) << std::endl;
    if (solver_type_ == hymos::common::kSolverNMG) {
      nmg_solver_->print_config();
    } else if (solver_type_ == hymos::common::kSolverFMG) {
      fmg_solver_->print_config();
    } else {
      single_solver_->print_config();
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

  void reset_tol(const double tol) {
    if (solver_type_ == hymos::common::kSolverNMG) {
      EnsureNMG()->reset_tol(tol);
    } else if (solver_type_ == hymos::common::kSolverFMG) {
      EnsureFMG()->reset_tol(tol);
    } else {
      EnsureSingle()->reset_tol(tol);
    }
  }

  void initial_value(sol_t& sol, const ptr_mesh_t& p_mesh, const ptr_force_t& p_force_handle) {
    if (solver_type_ == hymos::common::kSolverNMG) {
      EnsureNMG()->initial_value(sol, p_mesh, p_force_handle);
    } else if (solver_type_ == hymos::common::kSolverFMG) {
      EnsureFMG()->initial_value(sol, p_mesh, p_force_handle);
    } else {
      EnsureSingle()->initial_value(sol, p_mesh, p_force_handle);
    }
  }

  void Solve(sol_t& sol, const sol_t& rhs, const ptr_mesh_t& p_mesh, const ptr_force_t& p_force_handle, const unsigned int steps = 1) {
    if (solver_type_ == hymos::common::kSolverNMG) {
      EnsureNMG()->Solve(sol, rhs, p_mesh, p_force_handle, steps);
    } else if (solver_type_ == hymos::common::kSolverFMG) {
      EnsureFMG()->Solve(sol, rhs, p_mesh, p_force_handle, steps);
    } else {
      EnsureSingle()->Solve(sol, rhs, p_mesh, p_force_handle, steps);
    }
  }

  void PostProcessing(sol_t& sol, const ptr_mesh_t& p_mesh) {
    if (solver_type_ == hymos::common::kSolverNMG) {
      EnsureNMG()->PostProcessing(sol, p_mesh);
    } else if (solver_type_ == hymos::common::kSolverFMG) {
      EnsureFMG()->PostProcessing(sol, p_mesh);
    } else {
      EnsureSingle()->PostProcessing(sol, p_mesh);
    }
  }

private:
  static hymos::common::RuntimeSolverType NormalizeStoredType(int value) {
    if (value == (int)hymos::common::kSolverNMG) {
      return hymos::common::kSolverNMG;
    }
    if (value == (int)hymos::common::kSolverFMG) {
      return hymos::common::kSolverFMG;
    }
    return hymos::common::kSolverSingle;
  }

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
