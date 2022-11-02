#pragma once
#include "clucu_error.h"
#include "clucu_utils.h"
#include "clucu_config.h"
#include "clucu_core.h"
#include "clucu_background.h"
#include "clucu_probe.h"
#include "clucu_cls.h"

#include <gsl_sf_bessel.h>

double SamplingVariance_integrand_chi(double chi,void *params);
void clucu_compute_samplingvariance(clucu_cosmology *cosmo,int *status);
void clucu_compute_cov_NN(clucu_cosmology *cosmo,int *status);
void clucu_compute_cov_hhhh(clucu_cosmology *cosmo,int *status);
void clucu_compute_cov_hkhh(clucu_cosmology *cosmo,int *status);
void clucu_compute_cov_hkhk(clucu_cosmology *cosmo,int *status);

void SaveCov_NN(clucu_cosmology *cosmo);
void SaveCov_hhhh(clucu_cosmology *cosmo);
void SaveCov_hkhk(clucu_cosmology *cosmo);

