#pragma once
#include "clucu_power.h"
double GISDB_factor(clucu_cosmology *cosmo, double k, double z);
double HaloBiasPress1974(clucu_cosmology *cosmo,double sigma);
double HaloBiasSheth1999(clucu_cosmology *cosmo,double sigma);
double HaloBiasSheth2001(clucu_cosmology *cosmo,double sigma);
double HaloBiasTinker2010(clucu_cosmology *cosmo,double sigma);
double HaloBiasBhattacharya2011(clucu_cosmology *cosmo,double sigma);
double clucu_halo_bias_sigma(clucu_cosmology *cosmo,double sigma);
