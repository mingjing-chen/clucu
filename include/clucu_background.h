#pragma once
#include "clucu_error.h"
#include "clucu_utils.h"
#include "clucu_config.h"
#include "clucu_core.h"

#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_roots.h>

//species_x labels
typedef enum clucu_species_x_label {
  clucu_species_crit_label=0,
  clucu_species_m_label=1,
  clucu_species_l_label=2,
  clucu_species_g_label=3,
  clucu_species_k_label=4,
  clucu_species_ur_label=5,
  clucu_species_nu_label=6,

  clucu_species_cb_label=7,//transfer function额外的选项
  clucu_species_nu1_label=8,
  clucu_species_nu2_label=9,
  clucu_species_nu3_label=10,
} clucu_species_x_label;

typedef enum clucu_background_label {
  clucu_background_cla_label=0,
  clucu_background_class_label=1,
} clucu_background_label;

void RunClass(clucu_cosmology *cosmo);

double clucu_omega_x(clucu_cosmology *cosmo, double z, clucu_species_x_label label, int *status);

double clucu_rho_x(clucu_cosmology *cosmo, double z, clucu_species_x_label label, int is_comoving, int *status);

void clucu_compute_background(clucu_cosmology *cosmo,clucu_background_label label,int *status);

double clucu_h_over_h0(clucu_cosmology *cosmo,double z,int *status);

////哈勃参数，单位1/Mpc
double clucu_hubble_parameter(clucu_cosmology *cosmo,double z,int *status);

double clucu_comoving_radial_distance(clucu_cosmology *cosmo,double z,int *status);

double clucu_sinn(clucu_cosmology *cosmo, double chi,int *status);

double clucu_comoving_angular_distance(clucu_cosmology *cosmo,double z,int *status);

double clucu_angular_diameter_distance(clucu_cosmology *cosmo,double z,int *status);

double clucu_luminosity_distance(clucu_cosmology *cosmo,double z,int *status);

double clucu_growth_factor(clucu_cosmology *cosmo,double z,int *status);

double clucu_growth_rate(clucu_cosmology *cosmo,double z,int *status);

double clucu_chi_to_z(clucu_cosmology *cosmo,double chi,int *status);

double clucu_volume_element(clucu_cosmology *cosmo,double z,int *status);

