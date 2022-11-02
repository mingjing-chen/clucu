#pragma once
#include "clucu_power.h"

double MassFuncPress1974(double sigma,int *status);
double MassFuncSheth1999(clucu_cosmology *cosmo,double sigma,int *status);
double MassFuncJenkins2001(clucu_cosmology *cosmo,double sigma,int *status);
double MassFuncTinker2008(clucu_cosmology *cosmo,double sigma, double z,int *status);
double MassFuncBocquet2016(clucu_cosmology *cosmo,double sigma,double mass, double z,int *status);

double clucu_fitting_function(clucu_cosmology *cosmo,double sigma, double mass,double z,int *status);
double clucu_halo_dn_dlnm(clucu_cosmology *cosmo,double mass,double z,int *status);
double clucu_halo_bias(clucu_cosmology *cosmo,double mass,double z,int *status);
double clucu_halo_dndlnm_bias(clucu_cosmology *cosmo,double mass,double z,int *status);
double clucu_halo_dndlnm_bias_k(clucu_cosmology *cosmo,double mass,double z,double k,int *status);

void clucu_compute_sigmaS3(clucu_cosmology *cosmo,int *status);
double NGhmf_factor(clucu_cosmology *cosmo,double sigma, double mass,double z,int *status);
