#pragma once
#include "clucu_error.h"
#include "clucu_utils.h"
#include "clucu_config.h"
#include "clucu_core.h"
#include "clucu_background.h"
#include "clucu_probe.h"

#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp2d.h>

void clucu_compute_kernel_n(clucu_cosmology *cosmo,int *status);
void clucu_compute_kernel_nb(clucu_cosmology *cosmo,int *status);
double kernel_n(clucu_cosmology *cosmo,double ln_mass,double z,int *status);
double kernel_nb(clucu_cosmology *cosmo,double ln_mass,double z,int *status);

double weight_h(clucu_cosmology *cosmo,double ln_mass1,double ln_mass2,double z,int *status);

double clucu_barn(clucu_cosmology *cosmo,double ln_mass1,double ln_mass2,double z1,double z2,int *status);
void clucu_compute_barn(clucu_cosmology *cosmo,int *status);
void clucu_compute_cl_hh(clucu_cosmology *cosmo,int *status);
double clucu_cl_kk(clucu_cosmology *cosmo,double ell,double chi1,double chi2,int *status);
void clucu_compute_cl_kk(clucu_cosmology *cosmo,int *status);

double clucu_cl_hk2(clucu_cosmology *cosmo,double ell,double ln_mass1,double ln_mass2,double chi1,double chi2,int *status);
void clucu_compute_cl_hk2(clucu_cosmology *cosmo,int CEN_massob,int CEN_zob,int *status);
double clucu_cl_hk1(clucu_cosmology *cosmo,double ell,double ln_mass1,double ln_mass2,double chi1,double chi2,int *status);
void clucu_compute_cl_hk1(clucu_cosmology *cosmo,int CEN_massob,int CEN_zob,int *status);
double cl_kk_integrand(double chi,void *params);
//double cl_hh_integrand_ln_mass(double ln_mass,void *params);

//double cl_hk1_integrand_ln_mass(double ln_mass,void *params);


void SaveClz_hh(clucu_cosmology *cosmo);
void SaveClz_hh_hat(clucu_cosmology *cosmo);
void SaveClz_kk(clucu_cosmology *cosmo);
void SaveClz_kk_hat(clucu_cosmology *cosmo);
void SaveClz_hk1(clucu_cosmology *cosmo);
void SaveClz_hk2(clucu_cosmology *cosmo);
void SaveClz_hk(clucu_cosmology *cosmo);
