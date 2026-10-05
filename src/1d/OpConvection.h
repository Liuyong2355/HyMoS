/**
 * @file   OpConvection.h
 * @author CAI Zhenning
 * @date   Wed Sep 24 09:56:16 2014
 * 
 * @brief  用于计算对流
 * 
 * 
 */

#ifndef __OP_CONVECTION_H__
#define __OP_CONVECTION_H__

#include <NRxx/NumericalFlux.h>
# include <Eigen/Core>

template <typename MESH, typename SOLUTION, typename _BV>
class OpConvection {
public:
  typedef MESH mesh_t;
  typedef SOLUTION sol_t;

  typedef typename sol_t::dis_t dis_t;
  typedef typename dis_t::value_type value_t;
  typedef typename dis_t::Velocity velocity_t;

private:
  //对展开参数和展开系数进行重构
  static void Reconstruct( dis_t& re_l, dis_t& re_r,
                           const sol_t& sol, unsigned int i )
  {
    if ( i == 0 ) {
      re_l = sol[0]; re_l *= value_t( 1.5 ); re_l.Add( -0.5, sol[1] );
      re_r = sol[0]; re_r *= value_t( 0.5 ); re_r.Add( 0.5, sol[1] );
      re_l.Reinit( sol[0].Center() * 1.5 - sol[1].Center() * 0.5,
                   sol[0].Scaling() * 1.5 - sol[1].Scaling() * 0.5 );
      re_r.Reinit( sol[0].Center() * 0.5 + sol[1].Center() * 0.5,
                   sol[0].Scaling() * 0.5 + sol[1].Scaling() * 0.5 );
    } else if ( i == sol.size() - 1 ) {
      re_l = sol[i]; re_l *= value_t(0.5); re_l.Add( 0.5, sol[i-1] );
      re_r = sol[i]; re_r *= value_t(1.5); re_r.Add( -0.5, sol[i-1] );
      re_l.Reinit( sol[i].Center() * 0.5 + sol[i-1].Center() * 0.5,
                   sol[i].Scaling() * 0.5 + sol[i-1].Scaling() * 0.5 );
      re_r.Reinit( sol[i].Center() * 1.5 - sol[i-1].Center() * 0.5,
                   sol[i].Scaling() * 1.5 - sol[i-1].Scaling() * 0.5 );
    } else {
      velocity_t center_diff = value_t( 0.25 ) * (sol[i+1].Center() - sol[i-1].Center());
      velocity_t center_l( sol[i].Center() - center_diff );
      velocity_t center_r( sol[i].Center() + center_diff );

      value_t scale_diff = 0.25 * ( sol[i+1].Scaling() - sol[i-1].Scaling() );
      re_l.Reinit( center_l, sol[i].Scaling() - scale_diff, sol[i].GetOrder() );
      re_r.Reinit( center_r, sol[i].Scaling() + scale_diff, sol[i].GetOrder() );

      typename dis_t::Iterator it_re_l = re_l.begin(), it_re_r = re_r.begin();
      typename dis_t::ConstIterator it_l = sol[i-1].begin(), it_r = sol[i+1].begin();
      for (typename dis_t::ConstIterator it = sol[i].begin(), it_end = sol[i].end();
          it != it_end; ++it, ++it_l, ++it_r, ++it_re_l, ++it_re_r)
      {
        value_t g = (*it_r - *it_l) / 2;
        *it_re_l = *it - 0.5 * g;
        *it_re_r = *it + 0.5 * g;
      }
    }
  }
  
  //==2024.06.15(星期六)==
  //对守恒变量和展开系数做重构
  /*
  static void Reconstruct( dis_t& re_l, dis_t& re_r,
                           const sol_t& sol, unsigned int i )
  {
		Eigen::VectorXd cv_l, cv_m, cv_r, re_cv_l, re_cv_r;
		 if ( i == 0 ) {
			//对矩系数进行重构
			re_l = sol[0]; re_l *= value_t( 1.5 ); re_l.Add( -0.5, sol[1] );
			re_r = sol[0]; re_r *= value_t( 0.5 ); re_r.Add( 0.5, sol[1] );

			//对守恒变量进行重构
		  getcv(sol[0], cv_l);
		  getcv(sol[1], cv_r);
		  re_cv_l = 1.5 * cv_l - 0.5 * cv_r;
		  re_cv_r = 0.5 * cv_l + 0.5 * cv_r;
		 } else if ( i == sol.size() - 1 ) {
			re_l = sol[i]; re_l *= value_t(0.5); re_l.Add( 0.5, sol[i-1] );
			re_r = sol[i]; re_r *= value_t(1.5); re_r.Add( -0.5, sol[i-1] );
		
		  getcv(sol[i-1], cv_l);
		  getcv(sol[i], cv_r);
		  re_cv_l = 0.5 * cv_r + 0.5 * cv_l;
		  re_cv_r = 1.5 * cv_r - 0.5 * cv_l;
		 } else {
		  getcv(sol[i-1], cv_l);
		  getcv(sol[i], cv_m);
		  getcv(sol[i+1], cv_r);
		  re_cv_l = cv_m - 0.25 * (cv_r - cv_l); 
		  re_cv_r = cv_m + 0.25 * (cv_r - cv_l);

		  re_l.Reinit(sol[i]), re_r.Reinit(sol[i]);
			typename dis_t::Iterator it_re_l = re_l.begin(), it_re_r = re_r.begin();
			typename dis_t::ConstIterator it_l = sol[i-1].begin(), it_r = sol[i+1].begin();
			for (typename dis_t::ConstIterator it = sol[i].begin(), it_end = sol[i].end();
				 it != it_end; ++it, ++it_l, ++it_r, ++it_re_l, ++it_re_r)
			{
			  value_t g = (*it_r - *it_l) / 2;
			  *it_re_l = *it - 0.5 * g;
			  *it_re_r = *it + 0.5 * g;
			}
		 }

		velocity_t center;
		value_t scale;
		getPrim(re_cv_l, center, scale);
		 re_l.Reinit(center, scale);
		getPrim(re_cv_r, center, scale);
		 re_r.Reinit(center, scale);
  }

  static void getcv(const dis_t& dis, Eigen::VectorXd& cv){
		double pv[5];
		dis.PrimitiveVars(pv);
		double rho = pv[0];
		double u1 = pv[1];
		double u2 = pv[2];
		double u3 = pv[3];
		double theta = pv[4];
		cv.resize(5);
		cv[0] = rho;
		cv[1] = rho * u1;
		cv[2] = rho * u2;
		cv[3] = rho * u3;
		cv[4] = rho * (1.5 * theta + (u1*u1 + u2*u2 + u3*u3)/2);
  }

  static void getPrim(Eigen::VectorXd& cv,
		  velocity_t& center,
		  value_t& scale){
		double rho = cv[0];
		double u1 = cv[1] / rho;
		double u2 = cv[2] / rho;
		double u3 = cv[3] / rho;
		double theta = (cv[4] / rho - (u1*u1 + u2*u2 + u3*u3)/2) / 1.5;
		center = {u1, u2, u3};
		scale = sqrt(theta);
  }*/

  static void NumericalFlux( dis_t& flux,
                             const dis_t& dis_l, const dis_t& dis_r)
  {
    dis_t v_dis_l, v_dis_r;
    MulVelocity( 0, dis_l, v_dis_l );
    MulVelocity( 0, dis_r, v_dis_r );

    value_t pv_l[dis_t::dim + 2], pv_r[dis_t::dim + 2];
    dis_l.PrimitiveVars( pv_l );
    dis_r.PrimitiveVars( pv_r );
    
    unsigned int M = dis_l.GetOrder();
    value_t C = MaxRootOfHermitePolynomial<value_t>( M+1 );
    value_t lambda_ll = pv_l[1] - C * sqrt( pv_l[dis_t::dim + 1] ),
            lambda_lr = pv_r[1] - C * sqrt( pv_r[dis_t::dim + 1] ),
            lambda_rl = pv_l[1] + C * sqrt( pv_l[dis_t::dim + 1] ),
            lambda_rr = pv_r[1] + C * sqrt( pv_r[dis_t::dim + 1] );
    value_t lambda_l = std::min<value_t>( lambda_ll, lambda_lr ),
            lambda_r = std::max<value_t>( lambda_rl, lambda_rr );

    if ( lambda_l > 0 ) {
      MulVelocity( 0, dis_l, flux );
    } else if ( lambda_r < 0 ) {
      MulVelocity( 0, dis_r, flux );
    } else {
      value_t c1 = lambda_r / (lambda_r - lambda_l);
      value_t c2 = lambda_l / (lambda_r - lambda_l);
      value_t c3 = lambda_l * lambda_r / (lambda_r - lambda_l);

      dis_t v_dis_l, v_dis_r;
      MulVelocity( 0, dis_l, v_dis_l );
      MulVelocity( 0, dis_r, v_dis_r );
      flux = dis_l; flux *= -c3;
      flux.Add( c3, dis_r );
      flux.Add( c1, v_dis_l );
      flux.Add( -c2, v_dis_r );
    }
  }

public:
  static void Convection( dis_t& res_dis, const dis_t& dis, 
                          const sol_t& sol, const mesh_t& mesh,
                          const unsigned int i_ele )
  {
    dis_t dis_lr, dis_ml, dis_mr, dis_rl, dummy;
    if ( i_ele == 0 ) {
      Reconstruct( dis_ml, dis_mr, sol, i_ele );
      Reconstruct( dis_rl, dummy, sol, i_ele+1 );
      _BV::BoundaryValue( dis_lr, dis_ml, false );
    } else if ( i_ele == sol.size() - 1 ) {
      Reconstruct( dummy, dis_lr, sol, i_ele-1 );
      Reconstruct( dis_ml, dis_mr, sol, i_ele );
      _BV::BoundaryValue( dis_rl, dis_mr, true );
    } else {
      Reconstruct( dummy, dis_lr, sol, i_ele-1 );
      Reconstruct( dis_ml, dis_mr, sol, i_ele );
      Reconstruct( dis_rl, dummy, sol, i_ele+1 );
    }

    dis_t p_dis_lr, p_dis_ml, p_dis_mr, p_dis_rl;
    p_dis_lr.Reinit( dis ); Project( dis_lr, p_dis_lr );
    p_dis_ml.Reinit( dis ); Project( dis_ml, p_dis_ml );
    p_dis_mr.Reinit( dis ); Project( dis_mr, p_dis_mr );
    p_dis_rl.Reinit( dis ); Project( dis_rl, p_dis_rl );
    
    dis_t flux;
    value_t dx = mesh.Length( i_ele );
    NumericalFlux( flux, p_dis_lr, p_dis_ml );
    res_dis.Add( 1./dx, flux );

    NumericalFlux( flux, p_dis_mr, p_dis_rl );
    res_dis.Add( -1./dx, flux );
  }
};

#endif // __OP_CONVECTION_H__
