#pragma once
#include "clucu_error.h"
#include "clucu_utils.h"
#include "clucu_config.h"
#include "clucu_core.h"
#include "clucu_background.h"
#include "clucu_halo.h"
#include "clucu_survey.h"
#include "clucu_probe_utils.h"




double numbercounts_integrand_ln_mass_true(double ln_mass_true, void *params);
void clucu_compute_cluster_numbercounts(clucu_cosmology *cosmo,int CEN_massob,int CEN_zob,int *status);

double average_hmf_integrand_ln_mass_true(double ln_mass_true, void *params);
double average_hmfb_integrand_ln_mass_true(double ln_mass_true, void *params);
void clucu_compute_cluster_averagebias(clucu_cosmology *cosmo,int *status);
void clucu_compute_cluster_averagebias_k(clucu_cosmology *cosmo,int *status);
double clucu_averageBias_k(clucu_cosmology *cosmo, double k, double z,int *status);


void clucu_compute_cluster_volumeeffect(clucu_cosmology *cosmo,int *status);
void clucu_compute_cluster_powerspectrum(clucu_cosmology *cosmo,int *status);
