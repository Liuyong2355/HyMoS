/**
 * @file scat.h
 * @brief Collision models used by the current 1D FIM solver path
 */

#ifndef __SCAT_H__
#define __SCAT_H__

#include <NRxx/Distribution.h>

enum ColModel { BGK, ESBGK, SHAKHOV };

template <ColModel CM, class DISTRIBUTION, class PROBLEM>
class Scat;

template <ColModel CM, typename DISTRIBUTION, typename FORCE, typename _BV>
inline void Collision(DISTRIBUTION& rhs,
                      const DISTRIBUTION& dis,
                      const FORCE& force_handle,
                      const unsigned int i_ele) {
  double one_over_tau = _BV::OneOverTau(dis, force_handle, i_ele);
  Scat<CM, DISTRIBUTION, _BV>::Collision(rhs, dis, one_over_tau);
}

template <ColModel CM, typename DISTRIBUTION, typename FORCE, typename _BV>
inline void GetDeltaCollision(DISTRIBUTION& dest_dis,
                              const DISTRIBUTION& src_dis,
                              const double df,
                              const typename DISTRIBUTION::Indices& ind,
                              const FORCE& force_handle,
                              const unsigned int i_ele) {
  double one_over_tau = _BV::OneOverTau(src_dis, force_handle, i_ele);
  Scat<CM, DISTRIBUTION, _BV>::GetDeltaCollision(dest_dis, src_dis, df, ind, one_over_tau);
}

template <typename DISTRIBUTION>
inline void _DeltaLinearCollision(DISTRIBUTION& dest_dis,
                                  const double df,
                                  const typename DISTRIBUTION::Indices& ind,
                                  const double one_over_tau) {
  dest_dis(ind) += df * one_over_tau;
}

template <ColModel CM, class DISTRIBUTION, class PROBLEM>
class Scat {
public:
  static void Collision(DISTRIBUTION& dest_dis, const DISTRIBUTION& src_dis, const double one_over_tau);
  static void GetDeltaCollision(DISTRIBUTION& dest_dis,
                                const DISTRIBUTION& src_dis,
                                const double df,
                                const typename DISTRIBUTION::Indices& ind,
                                const double one_over_tau);
};

template <class DISTRIBUTION, class PROBLEM>
class Scat<BGK, DISTRIBUTION, PROBLEM> {
public:
  static void Collision(DISTRIBUTION& dest_dis, const DISTRIBUTION& src_dis, const double one_over_tau) {
    typename DISTRIBUTION::Iterator the_dis = dest_dis.begin(2);
    typename DISTRIBUTION::Iterator end_dis = dest_dis.end();
    typename DISTRIBUTION::ConstIterator the_src = src_dis.begin(2);
    for (; the_dis != end_dis; ++the_dis, ++the_src) {
      *the_dis -= one_over_tau * (*the_src);
    }
  }

  static void GetDeltaCollision(DISTRIBUTION& dest_dis,
                                const DISTRIBUTION& /*src_dis*/,
                                const double df,
                                const typename DISTRIBUTION::Indices& ind,
                                const double one_over_tau) {
    _DeltaLinearCollision(dest_dis, df, ind, one_over_tau);
  }
};

template <class DISTRIBUTION, class PROBLEM>
class Scat<ESBGK, DISTRIBUTION, PROBLEM> {
public:
  static void Collision(DISTRIBUTION& dest_dis, const DISTRIBUTION& src_dis, const double one_over_tau) {
    Scat<BGK, DISTRIBUTION, PROBLEM>::Collision(dest_dis, src_dis, one_over_tau * PROBLEM::Pr);

    const int dim = DISTRIBUTION::dim;
    double sigma[dim][dim];
    for (int i = 0; i < dim; ++i) {
      for (int j = 0; j < dim; ++j) {
        typename DISTRIBUTION::Indices ind;
        ind[i] += 1;
        ind[j] += 1;
        sigma[i][j] = src_dis(ind);
      }
      sigma[i][i] *= 2;
    }

    const double rho = src_dis.front();
    const double factor = (1 - 1.0 / PROBLEM::Pr) * one_over_tau * PROBLEM::Pr / rho;

    DISTRIBUTION dis_es;
    dis_es.Reinit(src_dis);
    dis_es.MakeMaxwellian(rho);

    typename DISTRIBUTION::Iterator the_dis = dis_es.begin(2);
    typename DISTRIBUTION::Iterator end_dis = dis_es.end();
    for (; the_dis != end_dis; ++the_dis) {
      typename DISTRIBUTION::Indices ind = the_dis.GetIndices();
      int sum = 0;
      for (int i = 0; i < dim; ++i) {
        sum += ind[i];
      }
      if (sum % 2 == 1) {
        continue;
      }

      for (int i = 0; i < dim; ++i) {
        if (ind[i] == 0) {
          continue;
        }
        double G_alpha = 0.0;
        ind[i] -= 1;
        for (int k = 0; k < dim; ++k) {
          if (ind[k] == 0) {
            continue;
          }
          ind[k] -= 1;
          G_alpha += sigma[i][k] * dis_es(ind);
          ind[k] += 1;
        }
        ind[i] += 1;
        *the_dis = G_alpha * factor / ind[i];
        break;
      }
    }

    typename DISTRIBUTION::Iterator dest_it = dest_dis.begin(2);
    typename DISTRIBUTION::ConstIterator src_it = dis_es.begin(2);
    for (; dest_it != dest_dis.end(); ++dest_it, ++src_it) {
      *dest_it += *src_it;
    }
  }

  static void GetDeltaCollision(DISTRIBUTION& dest_dis,
                                const DISTRIBUTION& /*src_dis*/,
                                const double df,
                                const typename DISTRIBUTION::Indices& ind,
                                const double one_over_tau) {
    _DeltaLinearCollision(dest_dis, df, ind, one_over_tau * PROBLEM::Pr);
  }
};

template <class DISTRIBUTION, class PROBLEM>
class Scat<SHAKHOV, DISTRIBUTION, PROBLEM> {
public:
  static void Collision(DISTRIBUTION& dest_dis, const DISTRIBUTION& src_dis, const double one_over_tau) {
    Scat<BGK, DISTRIBUTION, PROBLEM>::Collision(dest_dis, src_dis, one_over_tau);

    const int dim = DISTRIBUTION::dim;
    double heat_flux[dim];
    src_dis.HeatFlux(heat_flux);
    double factor = one_over_tau * (1 - PROBLEM::Pr) / (dim + 2);
    typename DISTRIBUTION::Indices ind;
    for (int i = 0; i < dim; ++i) {
      ind[i] = 1;
      for (int j = 0; j < dim; ++j) {
        ind[j] += 2;
        dest_dis(ind) += factor * heat_flux[i];
        ind[j] -= 2;
      }
      ind[i] = 0;
    }
  }

  static void GetDeltaCollision(DISTRIBUTION& dest_dis,
                                const DISTRIBUTION& /*src_dis*/,
                                const double df,
                                const typename DISTRIBUTION::Indices& ind,
                                const double one_over_tau) {
    _DeltaLinearCollision(dest_dis, df, ind, one_over_tau * PROBLEM::Pr);
  }
};

#endif
