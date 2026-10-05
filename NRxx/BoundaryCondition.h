// vim: fdm=syntax:ts=4:sw=4:syntax=cpp.doxygen:tw=70:fo+=Mm
#ifndef NRXX_BOUNDARYCONDITION_H
#define NRXX_BOUNDARYCONDITION_H

#include "Distribution.h"
#include "Transform.h"

namespace Reflection {

template <class VECTOR, class DISTRIBUTION>
void HalfMaxwellian(const VECTOR& pv, DISTRIBUTION& dis) {
	typedef typename DISTRIBUTION::value_type T;
	static const int DIM = DISTRIBUTION::dim;

	unsigned int M = dis.GetOrder();
	std::vector<double> K(M, 0);

	const typename DISTRIBUTION::Velocity& u = dis.Center();
	T theta = dis.Scaling() * dis.Scaling();
	K[0] = sqrt(pv[DIM + 1] / (2 * M_PI)) * exp(-(u[0] - pv[1]) * (u[0] - pv[1]) / (2 * pv[DIM + 1]));
	for (unsigned int i = 2; i < M; i += 2)
		K[i] = -K[i - 2] * (i - 1) / ((i + 1) * i) * theta;

	std::vector<double> H[DIM];
	H[0].resize(M + 1);
	H[0][0] = 0.5 * erfc((pv[1] - u[0]) / sqrt(2 * pv[DIM + 1]));
	H[0][1] = (pv[1] - u[0]) * H[0][0] - K[0];

	for (unsigned int i = 2; i <= M; i++) {
		H[0][i] = ((pv[1] - u[0]) * H[0][i - 1] + (pv[DIM + 1] - theta) * H[0][i - 2]) / i - K[i - 1];
	}

	for (int k = 1; k < DIM; k++) {
		H[k].resize(M + 1);
		H[k][0] = 1; H[k][1] = pv[k + 1] - u[k];
		for (unsigned int i = 2; i <= M; i++) 
			H[k][i] = ((pv[k + 1] - u[k]) * H[k][i - 1] + (pv[DIM + 1] - theta) * H[k][i - 2]) / i;
	}

	typename DISTRIBUTION::Iterator it = dis.begin();
	typename DISTRIBUTION::Iterator it_end = dis.end();
	for (; it != it_end; ++it) {
        const typename DISTRIBUTION::Indices& ind = it.GetIndices();
		if (it.NDis()) {
			*it = 0.5 * dis.delta[it.NDis()-1] * dis(ind, 0) * (theta - pv[DIM + 1]);
        } else {
			*it = pv[0];
			for (int k = 0; k < DIM; k++) *it *= H[k][ind[k]];
		}
	}
}

template <class DISTRIBUTION>
void Flip(DISTRIBUTION& dis) {
    typename DISTRIBUTION::Velocity c = dis.Center();
    c[0] = -c[0];
    dis.Reinit(c);

    typename DISTRIBUTION::Iterator it = dis.begin();
    typename DISTRIBUTION::Iterator it_end = dis.end();
    for (; it != it_end; ++it) {
        const typename DISTRIBUTION::Indices& ind = it.GetIndices();
        if (ind[0] % 2) *it = -*it;
    }
}

/**
 * 实现麦克斯韦边界条件。@p dis 在靠近边界的网格上的分布函数，@p u 是
 * 固壁的速度，@p theta 是固壁的温度，@p normal 是指向固壁内部的法向，
 * 即计算区域的外法向，用户需自行保证 @p normal 为一单位向量，@p chi
 * 是固壁的适应系数。函数完成后，ghost 网格单元上的分布函数存于 @p
 * result 中。
 */
template <class DISTRIBUTION>
void BoundaryValue(const DISTRIBUTION& dis,
        const typename DISTRIBUTION::Velocity& u,
        typename DISTRIBUTION::value_type theta,
        const std::vector<double>& normal,
        typename DISTRIBUTION::value_type chi,
        DISTRIBUTION& result)
{
    typedef typename DISTRIBUTION::value_type T;
	static const int DIM = DISTRIBUTION::dim;

    // 一系列 Givens 变换，把 normal 变成 e_1
    unsigned int n = normal.size();
    std::vector<T> psi(n - 1);
    std::vector<double> v = normal;
    typename DISTRIBUTION::Velocity u_rot(u);
    std::vector<DISTRIBUTION> tmp_dis(n);
    tmp_dis[n - 1] = dis;
    for (int i = n - 1; i > 0; i--) {
        psi[i - 1] = atan2(v[i], v[i - 1]);
        T cos_psi = cos(psi[i - 1]), sin_psi = sin(psi[i - 1]);
        v[i - 1] = cos_psi * v[i - 1] + sin_psi * v[i];

        T u_i = -u_rot[i - 1] * sin_psi + u_rot[i] * cos_psi;
        u_rot[i - 1] = u_rot[i - 1] * cos_psi + u_rot[i] * sin_psi;
        u_rot[i] = u_i;

        tmp_dis[i - 1].GivensTransform(i - 1, i, -psi[i - 1], tmp_dis[i]);
    }

    if (v[0] < 0) {
        u_rot[0] = -u_rot[0];
        Flip(tmp_dis[0]);
    }

    // 半空间 Maxwellian 部分
    std::vector<T> pv(DIM + 2);
    pv[0] = 1; pv[DIM + 1] = theta;
    for (int i = 0; i < DIM; i++) pv[i + 1] = u_rot[i];

    DISTRIBUTION half_Maxwellian;
    half_Maxwellian.Reinit(tmp_dis[0]);
    for (int i = 1; i < DIM; i++) u_rot[i] = tmp_dis[0].Center()[i];
    half_Maxwellian.Reinit(u_rot);  // 注意要放在以 u_rot 为速度的空间中做
    HalfMaxwellian(pv, half_Maxwellian);

    // 半分布函数部分
    std::vector<DISTRIBUTION> rotated_dis(n);
    rotated_dis[0].HalfSpaceCutOff(0, tmp_dis[0].Center()[0], tmp_dis[0]);
    rotated_dis[0].Reinit(u_rot);   // 切一半时只要切的位置与速度即可，但
                                    // 为了获得正确的动量，还必须把中心移
                                    // 至 u_rot 处

    // 组合在一起时需要保证边界通量为零
    std::vector<T> mnts_hm(DIM + 2), mnts_rd(DIM + 2);
    half_Maxwellian.Moments(mnts_hm);
    rotated_dis[0].Moments(mnts_rd);
    typename DISTRIBUTION::Indices ind(1,0,0);
    half_Maxwellian *= (0.5 * tmp_dis[0](ind) - mnts_rd[1]) / mnts_hm[1];
//half_Maxwellian *= -mnts_rd[1] / mnts_hm[1];

    // 边界外侧的分布函数
    typename DISTRIBUTION::Iterator it = tmp_dis[0].begin();
    typename DISTRIBUTION::Iterator it_end = tmp_dis[0].end();
    typename DISTRIBUTION::Iterator it_rot = rotated_dis[0].begin();
    typename DISTRIBUTION::Iterator it_hm = half_Maxwellian.begin();

    for (; it != it_end; ++it, ++it_rot, ++it_hm) {
        const typename DISTRIBUTION::Indices& ind = it.GetIndices();
        if (ind[0] % 2 == 1)
            *it_rot = (4 * chi * (*it_rot + *it_hm) - (2 + chi) * (*it)) / (2 - chi);
        else
            *it_rot = *it;
    }

    u_rot[0] = u_rot[0] * 2.0 - tmp_dis[0].Center()[0];
    rotated_dis[0].Reinit(u_rot);

    // 旋转回来
    if (v[0] < 0) Flip(rotated_dis[0]);
    for (unsigned int i = 0; i < n - 1; i++) {
        rotated_dis[i + 1].GivensTransform(i, i + 1, psi[i], rotated_dis[i]);
    }

    // 最终结果
    result = rotated_dis[n - 1];
}

}  // namespace Reflection

#endif
