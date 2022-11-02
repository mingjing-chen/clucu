#pragma once
#include "clucu_error.h"
#include "clucu_utils.h"
#include "clucu_config.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_const_mksa.h>
#include <gsl/gsl_roots.h>

// maximum number of species
#define CLUCU_MAX_NU_SPECIES 3
// limits for the precomputed spline of phase
// space diagram in MNU/T
#define CLUCU_NU_MNUT_MIN 1e-4
#define CLUCU_NU_MNUT_MAX 5000
// and number of points
#define CLUCU_NU_MNUT_N 2000

// The combination of constants required in Omeganuh2
//=(4*STBOLTZ/c^3) / rho_crit, 单位h^2
#define NU_CONST ( \
  8. * pow(M_PI,5) *pow((clucu_constants.KBOLTZ/ clucu_constants.HPLANCK),3)* \
  clucu_constants.KBOLTZ/(15. *pow( clucu_constants.CLIGHT,3))* \
  (8. * M_PI * clucu_constants.GNEWT) / \
  (3. * 100.*100.*1000.*1000. /clucu_constants.MPC_TO_METER /clucu_constants.MPC_TO_METER  * clucu_constants.CLIGHT * clucu_constants.CLIGHT))
//3，4行是 c^2 * rho_crit的倒数，rho_crit的单位是h^-2 kg/m^3。
//这里总共5个c，刚好使单位归一化。最终单位h^2
//1，2行：STBOLTZ/c = (2*pi^5*KBOLTZ^4)  / (15*HPLANCK^3) / CLIGHT^3
//最后结果 = (4*STBOLTZ/c) / (c^2*rho_crit) =  (4*STBOLTZ/c^3) / rho_crit。分子分母单位均为kg m^-3

// STBOLTZ,单位 W m^-2 K^-4 = kg s^-3 K^-4
// STBOLTZ/c^3, 单位 kg m^-3 K^-4
// T温度，单位K
// STBOLTZ/c^3*T^4，单位 kg m^-3，即为质量密度的单位
//前面还得乘 60/pi**2，因为STBOLTZ的定义，参考热统8.4节，STBOLTZ根据kb,hbar,c推出

typedef enum clucu_neutrino_mass_split_label{
    clucu_nu_splited = 0,
    clucu_nu_single = 1,
    clucu_nu_normal = 2,
    clucu_nu_inverted = 3,
    clucu_nu_equal = 4,
} clucu_neutrino_mass_split_label;

/**
 * Returns density of one neutrino species at a scale factor a.
 * Users are encouraged to access this quantity via the function clucu_omega_x.
 * @param a Scale factor
 * @param Neff The effective number of species with neutrino mass mnu.
 * @param mnu Pointer to array containing neutrino mass (can be 0).
 * @param T_CMB Temperature of the CMB
 * @param status Status flag. 0 if there are no errors, nonzero otherwise.
 * For specific cases see documentation for clucu_error.c
 * @return OmNuh2 Fractional energy density of neutrions with mass mnu, multiplied by h squared.
 */
//double clucu_Omeganuh2(double z, int N_nu_mass, double* mnu, double T_CMB, int * status);

/**
 * Returns mass of one neutrino species at a scale factor a.
 * @param a Scale factor
 * @param Neff The effective number of species with neutrino mass mnu.
 * @param OmNuh2 Fractional energy density of neutrions with mass mnu, multiplied by h squared. (can be 0).
 * @param T_CMB Temperature of the CMB
 * @param status Status flag. 0 if there are no errors, nonzero otherwise.
 * For specific cases see documentation for clucu_error.c
 * @return Mnu Neutrino mass [eV].
 */
//double* clucu_nu_masses(double OmNuh2, clucu_neutrino_mass_split_label mass_split, double T_CMB, int * status);

double nu_phasespace_intg(double mnuOT, int* status);

double density_WDM(double mass, double Temp);
double pressure_WDM(double mass, double Temp);
double EoS_WDM(double mass, double Temp);

double clucu_neutrino_mass_split(double sumnu, clucu_neutrino_mass_split_label label, double *mnu, int *N_nu_mass, int *status);

