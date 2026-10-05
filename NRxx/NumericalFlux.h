// vim: fdm=syntax:ts=4:sw=4:syntax=cpp.doxygen:tw=70:fo+=Mm

#ifndef _NumericalFlux_h_
#define _NumericalFlux_h_

#include <vector>
#include <numeric>
#include <algorithm>
#include "Vector.h"
#include "Transform.h"

template <class T>
inline T simpson(T l, T m, T r) {
     return (l + 4 * m + r) / 6;
}

template <class DISTRIBUTION>
void HLLFlux(const DISTRIBUTION& dis0, const DISTRIBUTION& dis1,
        const std::vector<double>& normal, DISTRIBUTION& flux_l2r, DISTRIBUTION& flux_r2l)
{
    typedef typename DISTRIBUTION::value_type T;

    unsigned int M = dis0.GetOrder();
    T C = MaxRootOfHermitePolynomial<T>(M + 1);

    const typename DISTRIBUTION::Velocity& c0 = dis0.Center();
    const typename DISTRIBUTION::Velocity& c1 = dis1.Center();
    const T& s0 = dis0.Scaling();
    const T& s1 = dis1.Scaling();
    T u0 = std::inner_product(normal.begin(), normal.end(), &c0[0], (T)0);
    T u1 = std::inner_product(normal.begin(), normal.end(), &c1[0], (T)0);

    T lambda_l = std::min<T>(u0 - C * s0, u1 - C * s1);
    T lambda_r = std::max<T>(u0 + C * s0, u1 + C * s1);
    static const int dim = DISTRIBUTION::dim;

    NRXX::Vector<dim + 2, T> pv_l, pv_r;
    dis0.PrimitiveVars(pv_l); 
    dis1.PrimitiveVars(pv_r);
    NRXX::Vector<dim + 2, T> pv_diff = pv_r - pv_l;

#ifndef NO_REG
    DISTRIBUTION nc_flux(M);
    typename DISTRIBUTION::Iterator it = nc_flux.begin(M);
    typename DISTRIBUTION::Iterator it_end = nc_flux.end(M);

    for (; it != it_end; ++it) {
        typename DISTRIBUTION::Indices ind = it.GetIndices();
        for (unsigned int j = 0; j < normal.size(); j++) {
            T tmp = 0;
            ind[j] += 1;
            int n_dis = it.NDis();

            for (int d = 0; d < dim; d++) {
                if (ind[d] > 0) {
                    ind[d] -= 1;
                    const T& f_l = dis0(ind, n_dis);
                    const T& f_r = dis1(ind, n_dis);

                    tmp -= simpson(f_l, (f_l / pv_l[0] + f_r / pv_r[0]) * (pv_l[0] + pv_r[0]) / 4, f_r) * pv_diff[d + 1];
                    ind[d] += 1;
                }

                if (ind[d] > 1) {
                    ind[d] -= 2;
                    const T& f_l = dis0(ind, n_dis);
                    const T& f_r = dis1(ind, n_dis);
                    tmp -= 0.5 * simpson(f_l, (f_l / pv_l[0] + f_r / pv_r[0]) * (pv_l[0] + pv_r[0]) / 4, f_r) * pv_diff[dim + 1];
                    ind[d] += 2;
                }
            }

            if (n_dis) {
                const T& f_l = -dis0.delta[n_dis-1] * dis0(ind, 0);
                const T& f_r = -dis1.delta[n_dis-1] * dis1(ind, 0);
                tmp -= 0.5 * simpson(f_l, (f_l / pv_l[0] + f_r / pv_r[0]) * (pv_l[0] + pv_r[0]) / 4, f_r) * pv_diff[dim + 1];
            }

            tmp *= ind[j];
            *it += tmp * normal[j];

            ind[j] -= 1;
        }
    }
#endif // NO_REG

    if (lambda_l > 0) {
        flux_l2r.Reinit(dis0);
        flux_r2l.Reinit(dis0);
        for (size_t i = 0, n = normal.size(); i < n; i++) {
            DISTRIBUTION tmp;
            MulVelocity(i, dis0, tmp);
            flux_l2r.Add(normal[i], tmp);
        }
        flux_r2l.Add(-1.0, flux_l2r);
#ifndef NO_REG
        flux_r2l += nc_flux;
#endif // NO_REG
    } else if (lambda_r < 0) {
        flux_l2r.Reinit(dis1);
        flux_r2l.Reinit(dis1);
        for (size_t i = 0, n = normal.size(); i < n; i++) {
            DISTRIBUTION tmp;
            MulVelocity(i, dis1, tmp);
            flux_l2r.Add(normal[i], tmp);
        }
        flux_r2l.Add(-1.0, flux_l2r);
#ifndef NO_REG
        flux_l2r += nc_flux;
#endif // NO_REG
    } else {
        T c1 = lambda_r / (lambda_r - lambda_l);
        T c2 = lambda_l / (lambda_r - lambda_l);
        T c3 = lambda_l * lambda_r / (lambda_r - lambda_l);

        flux_l2r.Reinit(dis0);
        for (size_t i = 0, n = normal.size(); i < n; i++) {
            DISTRIBUTION tmp;
            MulVelocity(i, dis0, tmp);
            flux_l2r.Add(c1 * normal[i], tmp);
        }

        DISTRIBUTION tmp1;
        tmp1.Reinit(dis0);
        tmp1.ResetOrder(M + 1);
        Project(dis1, tmp1);

        for (size_t i = 0, n = normal.size(); i < n; i++) {
            DISTRIBUTION tmp;
            MulVelocity(i, tmp1, tmp);
            flux_l2r.Add(-c2 * normal[i], tmp);
        }

        flux_l2r.Add(c3, tmp1);
        flux_l2r.Add(-c3, dis0);
        flux_r2l = flux_l2r; flux_r2l *= -1.0;
#ifndef NO_REG
        flux_l2r.Add(-c2, nc_flux);
        flux_r2l.Add(c1, nc_flux);
#endif // NO_REG
    }
}

#endif // _NumericalFlux_h_
