#pragma once
#include "clucu_error.h"
#include "clucu_utils.h"
#include "clucu_config.h"
#include "clucu_core.h"
#include "clucu_background.h"
#include "clucu_power.h"

#define delta_short_constant_centroid 0.015
#define boost_initial_conditions 1.0
#define precision_scale 10
#define precision_normalization 4
#define z_collapse  0.7
#define z_initial 200.
#define z_collapse_min 0.001
#define N_delta_long 2
#define delta_long_max (delta_short_constant_centroid/10.0)

#define do_clustering 0
#define debug_mode 1



//class改为输出tk
//cosmo->data部分，改为插值tk
//对每个mass(Nm),k(Nk),delta_L(2个)，都要计算塌缩方程，求delta_i

int clucu_collapse(clucu_cosmology *cosmo,int *status);
double clucu_SoverL(clucu_cosmology *cosmo,double mass,double k,int *status);
double clucu_delta_crit(clucu_cosmology *cosmo,double mass,double k,int *status);


double find_z_collapse_nothing
(clucu_cosmology * const cosmo, const double Ri, const double Rpi, const double delta_long,const double k_long,const double Mhalo);
double find_z_collapse_1nu
(clucu_cosmology * const cosmo, const double Ri, const double Rpi, const double delta_long,const double k_long,const double Mhalo,
const long Nz_solution, double *R_solution, double *Mnu_solution);
double find_z_collapse_2nu
(clucu_cosmology * const cosmo, const double Ri, const double Rpi, const double delta_long,const double k_long,const double Mhalo,
const long Nz_solution, double *R_solution, double *Mnu1_solution, double *Mnu2_solution);
double find_z_collapse_3nu
(clucu_cosmology * const cosmo, const double Ri, const double Rpi, const double delta_long,const double k_long,const double Mhalo,
const long Nz_solution, double *R_solution, double *Mnu1_solution, double *Mnu2_solution, double *Mnu3_solution);