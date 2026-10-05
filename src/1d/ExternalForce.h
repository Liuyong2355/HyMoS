/**
 * @file ExternalForce.h
 * @brief Minimal force support for the current 1D gas-flow cases
 */

#ifndef __EXTERNAL_FORCE_H__
#define __EXTERNAL_FORCE_H__

#include <memory>
#include <ostream>

#include "ConfigFile.h"

template <typename MESH, typename SOLUTION, typename _BV>
class NoForce {
public:
  typedef NoForce<MESH, SOLUTION, _BV> force_t;
  typedef MESH mesh_t;
  typedef std::shared_ptr<mesh_t> ptr_mesh_t;
  typedef SOLUTION sol_t;
  typedef typename sol_t::dis_t dis_t;
  typedef typename dis_t::Indices index_t;

  NoForce() {}
  NoForce(const ptr_mesh_t& /*mesh*/) {}
  NoForce(const ptr_mesh_t& mesh, const ConfigFile& cf) { read_config(mesh, cf); }
  NoForce(const ptr_mesh_t& /*mesh*/, std::istream& /*is*/) {}

  void Reinit(const ptr_mesh_t& /*p_mesh*/) {}
  void read_config(const ptr_mesh_t& /*mesh*/, const ConfigFile& /*cf*/) {}
  void print_config() const {}
  void dump_data(std::ostream& /*os*/) const {}
  void save_results(unsigned int /*n_run*/) const {}

  void GetAcceleration(dis_t& /*rhs*/, const dis_t& /*dis*/, const unsigned int /*i_ele*/) const {}
  void GetDeltaAcceleration(dis_t& /*dest_dis*/,
                            const dis_t& /*src_dis*/,
                            const double /*df*/,
                            const index_t& /*ind*/,
                            const unsigned int /*i_ele*/) const {}

  void initial_value(const sol_t& /*sol*/) {}
  void CalcExternalForce(const sol_t& /*sol*/) {}
  void CalcExternalForce(const sol_t& /*sol*/, const unsigned int /*i_ele*/) {}
  void CalcReverseExternalForce(const sol_t& /*sol*/, const unsigned int /*i_ele*/) {}
  void CoarseFrom(const force_t& /*force*/) {}
  void CoarseForceEquationRHS(const sol_t& /*sol*/) {}
  void BackSubstitution(const force_t& /*coarse_force*/) {}
  void CalcCorrection(const force_t& /*old_force*/) {}
};

#endif
