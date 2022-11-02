#pragma once
#include "clucu_power.h"


double ConcentrationDuffy2008(clucu_cosmology *cosmo,halo_define_type type,double mass,double z,int *status);

double ConcentrationBhattacharya2013(clucu_cosmology *cosmo,halo_define_type type,double mass,double z,int *status);

double ConcentrationDiemer2015(clucu_cosmology *cosmo,halo_define_type type,double mass,double z,int *status);

double clucu_concentration(clucu_cosmology *cosmo,halo_define_type type,double mass, double z,int *status);
