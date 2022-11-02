#include "clucu_power_utils.h"
/*-----------------------------------------------------------------------------*/
/*                                       插值结构体                                  */
/*-----------------------------------------------------------------------------*/
clucu_f2d_t *clucu_f2d_t_new(int nz,double *z_arr,
                         int nk,double *lk_arr,
                         double *fkz_arr,
                         double *fk_arr,
                         double *fz_arr,
                         int is_factorizable,
                         int extrap_order_lok,
                         int extrap_order_hik,
                         clucu_f2d_extrap_growth_t extrap_linear_growth,
                         int is_fkz_log,
                         double growth_factor_0,
                         int growth_exponent,
                         clucu_f2d_interp_t interp_type,
                         int *status) {
  int s2dstatus=0;
  clucu_f2d_t *f2d = malloc(sizeof(clucu_f2d_t));
  if (f2d == NULL)
    *status = CLUCU_ERROR_MEMORY;

  if (*status == 0) {
    is_factorizable = is_factorizable || (z_arr == NULL) || (lk_arr == NULL) || (fkz_arr == NULL);
    f2d->is_factorizable = is_factorizable;
    f2d->is_k_constant = ((lk_arr == NULL) || ((fkz_arr == NULL) && (fk_arr == NULL)));
    f2d->is_z_constant = ((z_arr == NULL) || ((fkz_arr == NULL) && (fz_arr == NULL)));
    f2d->extrap_order_lok = extrap_order_lok;
    f2d->extrap_order_hik = extrap_order_hik;
    f2d->extrap_linear_growth = extrap_linear_growth;
    f2d->is_log = is_fkz_log;
    f2d->growth_factor_0 = growth_factor_0;
    f2d->growth_exponent = growth_exponent;
    f2d->fkz = NULL;
    f2d->fk = NULL;
    f2d->fz = NULL;

    if (!(f2d->is_k_constant)) { //If it's not constant
      f2d->lkmin = lk_arr[0];
      f2d->lkmax = lk_arr[nk-1];
    }
    if (!(f2d->is_z_constant)) {
      f2d->zmin = z_arr[0];
      f2d->zmax = z_arr[nz-1];
    }
  }

  if ((extrap_order_lok > 2) || (extrap_order_lok < 0) || (extrap_order_hik > 2) || (extrap_order_hik < 0))
    *status = CLUCU_ERROR_INCONSISTENT;

  if ((extrap_linear_growth != clucu_f2d_clucugrowth) &&
      (extrap_linear_growth != clucu_f2d_constantgrowth) &&
      (extrap_linear_growth != clucu_f2d_no_extrapol))
    *status = CLUCU_ERROR_INCONSISTENT;

if(*status == 0) {
    switch(interp_type) 
    {
        case(clucu_f2d_3):
            //如果T(k),D(z)可以分开，比如EH98和BBKS
            //我们就分别设置fk,fz
            if (f2d->is_factorizable) 
            {
                // Do not allocate spline if constant
                if(f2d->is_k_constant)
                f2d->fk = NULL;
                else { //Otherwise allocate and check
                f2d->fk = gsl_spline_alloc(gsl_interp_cspline, nk);
                if(f2d->fk == NULL)
                *status = CLUCU_ERROR_MEMORY;
                }

                // Do not allocate spline if constant
                if (f2d->is_z_constant)
                    f2d->fz = NULL;
                else
                { //Otherwise allocate and check
                    f2d->fz = gsl_spline_alloc(gsl_interp_cspline, nz);
                    if (f2d->fz == NULL)
                    *status = CLUCU_ERROR_MEMORY;
                }
            }
            //如果T(k),D(z)不能分开，比如有中微子的影响D(k,z)，用CLASS输出的功率谱P(k,z)
            //我们就整体设置fkz
            else 
            {
                // Do not allocate spline if constant
                if ((f2d->is_k_constant) || (f2d->is_z_constant))
                    f2d->fkz = NULL;
                else
                { //Otherwise allocate and check
                    f2d->fkz = gsl_spline2d_alloc(gsl_interp2d_bicubic, nk, nz);
                    if (f2d->fkz == NULL)
                        *status = CLUCU_ERROR_MEMORY;
                }
            }
            break;
        default:
            f2d->fk = NULL;
            f2d->fz = NULL;
            f2d->fkz = NULL;
    }
}

  if (*status == 0) 
  {
    if (f2d->is_factorizable) 
    {
      if (f2d->fk != NULL)
        s2dstatus |= gsl_spline_init(f2d->fk, lk_arr, fk_arr, nk);
      if (f2d->fz != NULL)
        s2dstatus |= gsl_spline_init(f2d->fz, z_arr, fz_arr, nz);
    }
    else {
      if (f2d->fkz != NULL)
        s2dstatus=gsl_spline2d_init(f2d->fkz, lk_arr, z_arr, fkz_arr, nk, nz);
    }
    if (s2dstatus)
      *status = CLUCU_ERROR_SPLINE;
  }

  return f2d;
}
//输入log_pkz, logk, z。 输出pkz。不管插值时是不是log，输出不带log
double clucu_f2d_t_eval(clucu_f2d_t *f2d,double lk,double z,void *cosmo, int *status) 
{
    //z是否超出区间
    int is_hiz, is_loz;
    double z_ev = z;
    if (f2d->is_z_constant) 
    {
        is_hiz = 0;
        is_loz = 0;
    }
    else 
    {
        is_hiz = z > f2d->zmax;
        is_loz = z < f2d->zmin;
        if (is_loz) { // Are we above the interpolation range in a?
            if (f2d->extrap_linear_growth == clucu_f2d_no_extrapol) {
            *status=CLUCU_ERROR_SPLINE_EV;
            return NAN;
            }
            z_ev = f2d->zmin;
        }
        else if (is_hiz) { // Are we below the interpolation range in a?
            if (f2d->extrap_linear_growth == clucu_f2d_no_extrapol) {
            *status=CLUCU_ERROR_SPLINE_EV;
            return NAN;
            }
            z_ev = f2d->zmax;
        }
    }

    //k是否超出区间
    int is_hik, is_lok;
    double fkz_pre, fkz_post;
    double lk_ev = lk;
    if (f2d->is_k_constant) 
    {
        is_hik = 0;
        is_lok = 0;
    }
    else 
    {
        is_hik = lk > f2d->lkmax;
        is_lok = lk < f2d->lkmin;
        if (is_hik) // Are we above the interpolation range in k?
            lk_ev = f2d->lkmax;
        else if (is_lok) // Are we below the interpolation range in k?
            lk_ev = f2d->lkmin;
    }

    //先求ev的值
    //如果没有超出边界，则直接return ev的值
    //如果超出了边界，则ev是上边界/下边界，先求边界的值
    // Evaluate spline
    int spstatus=0;
    if (f2d->is_factorizable) 
    {
        double fk, fz;
        if (f2d->fk == NULL) 
        {
            if (f2d->is_log)
                fk = 0;
            else
                fk = 1;
        }
        else
            spstatus |= gsl_spline_eval_e(f2d->fk, lk_ev, NULL, &fk);


        if (f2d->fz == NULL) 
        {
            if (f2d->is_log)
            fz = 0;
            else
            fz = 1;
        }
        else
            spstatus |= gsl_spline_eval_e(f2d->fz, z_ev, NULL, &fz);
        
        if (f2d->is_log)
            fkz_pre = fk+fz;
        else
            fkz_pre = fk*fz;
    }
    else 
    {
        if (f2d->fkz == NULL) 
        {
            if (f2d->is_log)
                fkz_pre = 0;
            else
                fkz_pre = 1;
        }
        else
            spstatus = gsl_spline2d_eval_e(f2d->fkz, lk_ev, z_ev, NULL, NULL, &fkz_pre);
    }

    if (spstatus) 
    {
        CLUCU_RAISE_GSL_WARNING(spstatus,"k=%.2e,z=%.2e\n",exp(lk_ev),z_ev);
        *status = CLUCU_ERROR_SPLINE_EV;
        return NAN;
    }

    //开始外推
    // Now extrapolate in k if needed
    if (is_hik) 
    {
        fkz_post = fkz_pre;
        //logk线性外推(边界处一阶导数不变)
        if (f2d->extrap_order_hik > 0) 
        {   
            double pd;
            double dlk = lk-lk_ev;
            if (f2d->is_factorizable)
                spstatus = gsl_spline_eval_deriv_e(f2d->fk, lk_ev, NULL, &pd);
            else
                spstatus = gsl_spline2d_eval_deriv_x_e(f2d->fkz, lk_ev, z_ev, NULL, NULL, &pd);
            if (spstatus) 
            {
                CLUCU_RAISE_GSL_WARNING(spstatus,"k=%.2e,z=%.2e\n",exp(lk_ev),z_ev);
                *status = CLUCU_ERROR_SPLINE_EV;
                return NAN;
            }
            fkz_post += pd*dlk;
            if (f2d->extrap_order_hik > 1) 
            //logk外推(边界处 2 阶导数不变)
            {
                if (f2d->is_factorizable)
                    spstatus = gsl_spline_eval_deriv2_e(f2d->fk, lk_ev, NULL, &pd);
                else
                    spstatus = gsl_spline2d_eval_deriv_xx_e(f2d->fkz, lk_ev, z_ev, NULL, NULL, &pd);
                if (spstatus) 
                {
                    CLUCU_RAISE_GSL_WARNING(spstatus,"k=%.2e,z=%.2e\n",exp(lk_ev),z_ev);
                    *status=CLUCU_ERROR_SPLINE_EV;
                    return NAN;
                }
                fkz_post += pd*dlk*dlk*0.5;
            }
        }
    }
    else if (is_lok) 
    {
        fkz_post = fkz_pre;
        if (f2d->extrap_order_lok > 0) 
        {
            double pd;
            double dlk = lk-lk_ev;
            if (f2d->is_factorizable)
                spstatus = gsl_spline_eval_deriv_e(f2d->fk, lk_ev, NULL, &pd);
            else
                spstatus = gsl_spline2d_eval_deriv_x_e(f2d->fkz, lk_ev, z_ev, NULL, NULL, &pd);
            if (spstatus) 
            {
                CLUCU_RAISE_GSL_WARNING(spstatus,"k=%.2e,z=%.2e\n",exp(lk_ev),z_ev);
                *status = CLUCU_ERROR_SPLINE_EV;
                return NAN;
            }
            fkz_post += pd*dlk;

            if (f2d->extrap_order_lok > 1) 
            {
                if (f2d->is_factorizable)
                    spstatus = gsl_spline_eval_deriv2_e(f2d->fk, lk_ev, NULL, &pd);
                else
                    spstatus = gsl_spline2d_eval_deriv_xx_e(f2d->fkz, lk_ev, z_ev, NULL, NULL, &pd);
                if (spstatus) 
                {
                    CLUCU_RAISE_GSL_WARNING(spstatus,"k=%.2e,z=%.2e\n",exp(lk_ev),z_ev);
                    *status = CLUCU_ERROR_SPLINE_EV;
                    return NAN;
                }
                fkz_post += pd*dlk*dlk*0.5;
            }
        }
    }
    else
    fkz_post = fkz_pre;

    // Exponentiate if needed
    if (f2d->is_log)
    fkz_post = exp(fkz_post);

    // Extrapolate in a if needed
    if (is_hiz) 
    {
        double gz;
        if (f2d->extrap_linear_growth == clucu_f2d_clucugrowth) 
        { // Use CLUCU's growth function
            clucu_cosmology *csm = (clucu_cosmology *)cosmo;
            if (!csm->computed_growth) 
            {
                *status = CLUCU_ERROR_GROWTH_INIT;
                CLUCU_RAISE_WARNING(*status,"growth factor splines have not been precomputed!");
                clucu_cosmology_set_status_message(
                    csm,
                    "clucu_f2d.c: clucu_f2d_t_eval(): growth factor splines have not been precomputed!");
                return NAN;
            }
            gz = (
            clucu_growth_factor(csm, z,status) /
            clucu_growth_factor(csm, z_ev,status));
        }
        else if( f2d->extrap_linear_growth == clucu_f2d_constantgrowth)
         // Use constant growth factor
            gz = f2d->growth_factor_0;

        fkz_post *= pow(gz, f2d->growth_exponent);
    }

    return fkz_post;
}

void clucu_f2d_t_free(clucu_f2d_t *f2d)
{
  if(f2d != NULL) {
    if(f2d->fkz != NULL)
      gsl_spline2d_free(f2d->fkz);
    if(f2d->fk != NULL)
      gsl_spline_free(f2d->fk);
    if(f2d->fz != NULL)
      gsl_spline_free(f2d->fz);
    free(f2d);
  }
}