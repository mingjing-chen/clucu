#pragma once
#include "clucu_error.h"
#include "clucu_utils.h"
#include "clucu_config.h"
#include "clucu_core.h"
#include "clucu_background.h"
#include "clucu_halo.h"

double zob_to_ztrue(clucu_cosmology *cosmo,double z_ob,double *out);
static double expective_ln_mass_ob1(clucu_cosmology *cosmo, double ln_mass_true,double z_true);
static double expective_ln_mass_ob2(clucu_cosmology *cosmo, double ln_mass_true,double z_true);
static double sigma_ln_mass_ob1(clucu_cosmology *cosmo,double ln_mass_true,double z_true);
double sigma_ln_mass_ob2(clucu_cosmology *cosmo,double ln_mass_true,double z_true);
double mass_ob_probability(clucu_cosmology *cosmo,double ln_mass_ob, double ln_mass_true,double z_true);
double redshift_probability(clucu_cosmology *cosmo,double z_ob,double z_true);

void clucu_compute_average_one_over_chis(clucu_cosmology *cosmo);


static double NFW_fx_x3(double x);
static double Mf3(double x, void *params);
double M200cz_to_M180mz(clucu_cosmology *cosmo,double M200c, double z);
static double Mf4(double x, void *params);
double M180mz_to_M200cz(clucu_cosmology *cosmo,double M180m, double z);

double M200mz_to_M500cz(clucu_cosmology *cosmo,double M200m, double z);

double M500cz_to_M200mz(clucu_cosmology *cosmo,double M500c, double z);

static double Mmin_SZE180(clucu_cosmology *cosmo,double z);
static double Mmin_xray180(clucu_cosmology *cosmo,double z);
double Mmin_limit(clucu_cosmology *cosmo,double z);