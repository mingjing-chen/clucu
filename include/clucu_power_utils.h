#pragma once

#include "clucu_background.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>
/**
 * Create a clucu_f2d_t structure.
 * @param na number of elements in a_arr.
 * @param a_arr array of scale factor values at which the function is defined. The array should be ordered.
 * @param nk number of elements of lk_arr.
 * @param lk_arr array of logarithmic wavenumbers at which the function is defined (i.e. this array contains ln(k), NOT k). The array should be ordered.
 * @param fka_arr array of size na * nk containing the 2D function. The 2D ordering is such that fka_arr[ia*nk+ik] = f(k=exp(lk_arr[ik]),a=a_arr[ia]).
 * @param fk_arr array of size nk containing the k-dependent part of the function. Only relevant if is_factorizable is true.
 * @param fa_arr array of size na containing the a-dependent part of the function. Only relevant if is_factorizable is true.
 * @param is_factorizable if not 0, fk_arr and fa_arr will be used as 1-D arrays to construct a factorizable 2D function.
 * @param extrap_order_lok Order of the polynomial that extrapolates on wavenumbers smaller than the minimum of lk_arr. Allowed values: 0 (constant), 1 (linear extrapolation) and 2 (quadratic extrapolation). Extrapolation happens in ln(k).
 * @param extrap_order_hik Order of the polynomial that extrapolates on wavenumbers larger than the maximum of lk_arr. Allowed values: 0 (constant), 1 (linear extrapolation) and 2 (quadratic extrapolation). Extrapolation happens in ln(k).
 * @param extrap_linear_growth: clucu_f2d_extrap_growth_t value defining how the function with scale factors below the interpolation range. Allowed values: clucu_f2d_clucugrowth (scale with the CLUCU linear growth factor), clucu_f2d_constantgrowth (scale by multiplying the function at the earliest available scale factor by a constant number, defined by `growth_factor_0`), clucu_f2d_no_extrapol (throw an error if the function is ever evaluated outside the interpolation range in a). Note that, above the interpolation range (i.e. for low redshifts), the function will be assumed constant.
 * @param is_fka_log: if not zero, `fka_arr` contains ln(f(k,a)) instead of f(k,a). If the function is factorizable, then `fk_arr` holds ln(K(k)) and `fa_arr` holds ln(A(a)), where f(k,a)=K(k)*A(a).
 * @param growth_factor_0: custom growth function. Irrelevant if extrap_linear_growth!=clucu_f2d_constantgrowth.
 * @param growth_exponent: power to which the extrapolating growth factor should be exponentiated when extrapolating (e.g. usually 2 for linear power spectra).
 * @param interp_type: 2D interpolation method. Currently only clucu_f2d_3 is implemented (bicubic interpolation).
 * @param status Status flag. 0 if there are no errors, nonzero otherwise.
 */
clucu_f2d_t *clucu_f2d_t_new(int nz,double *z_arr,
			 int nk,double *lk_arr,
			 double *fkz_arr,
			 double *fk_arr,
			 double *fz_arr,
			 int is_factorizable,
			 int extrap_order_lok,
			 int extrap_order_hik,
			 clucu_f2d_extrap_growth_t extrap_linear_growth,
			 int is_fka_log,
			 double growth_factor_0,
			 int growth_exponent,
			 clucu_f2d_interp_t interp_type,
			 int *status);

/**
 * Evaluate 2D function of k and a defined by clucu_f2d_t structure.
 * @param fka clucu_f2d_t structure defining f(k,a).
 * @param lk Natural logarithm of the wavenumber.
 * @param a Scale factor.
 * @param cosmo clucu_cosmology structure, only needed if evaluating f(k,a) at small scale factors outside the interpolation range, and if fka was initialized with extrap_linear_growth = clucu_f2d_clucugrowth.
 * @param status Status flag. 0 if there are no errors, nonzero otherwise.
 */
double clucu_f2d_t_eval(clucu_f2d_t *fka,double lk,double a,void *cosmo,
		      int *status);

/**
 * F2D structure destructor.
 * Frees up all memory associated with a f2d structure.
 * @param fka Structure to be freed.
 */
void clucu_f2d_t_free(clucu_f2d_t *fka);
