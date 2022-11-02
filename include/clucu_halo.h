#pragma once
#include "clucu_error.h"
#include "clucu_utils.h"
#include "clucu_config.h"
#include "clucu_core.h"
#include "clucu_background.h"
#include "clucu_power.h"
#include "clucu_halo_bias.h"
#include "clucu_halo_massfunction.h"
#include "clucu_halo_concentration.h"


double nfw_mc(double c);
void clucu_convert_concentration(clucu_cosmology *cosmo,
			       double delta_old, int nc, double c_old[],
			       double delta_new, double c_new[],int *status);

double clucu_Mold_to_Mnew(clucu_cosmology *cosmo,halo_define_type Told,halo_define_type Tnew,double Mold,double z,int *status);
double clucu_M200m_to_M500c(clucu_cosmology *cosmo,double M200m,double z,int *status);
double clucu_M200c_to_M180m(clucu_cosmology *cosmo,double M200c,double z,int *status);
double clucu_M200m_to_M200c(clucu_cosmology *cosmo,double M200m,double z,int *status);
double clucu_M200c_to_M200m(clucu_cosmology *cosmo,double M200c,double z,int *status);

void clucu_compute_massconvetr(clucu_cosmology *cosmo,halo_define_type Told,halo_define_type T200,int *status);
double clucu_Mold_to_M200(clucu_cosmology *cosmo,halo_define_type Told,halo_define_type T200, double Mold, double z,int *status);


double nfw_scale(clucu_cosmology *cosmo,double mass,double z,double *out,int *status);
double nfw_uk(clucu_cosmology *cosmo,double k,double mass,double z,int *status);
