#include "clucu_neutrino.h"
//TODO:改为对分布函数的积分


// Global variable to hold the neutrino phase-space spline
gsl_spline* density_spline = NULL;
gsl_spline* pressure_spline = NULL;

// these are NOT adjustable
// this phase space integral is only done once and the following is accurate
// enough according to tests done by the devs
/**
 * Absolute precision in neutrino root finding
 */
#define GSL_EPSABS_NU 1E-7

/**
 * Relative precision in neutrino root finding
 */
#define GSL_EPSREL_NU 1E-7

/**
 * Number of iterations for neutrino root finding
 */
#define GSL_N_ITERATION_NU 1000



/* ------- ROUTINE: nu_integrand ------
INPUTS: x: dimensionless momentum, *r: pointer to a dimensionless mass / temperature
TASK: Integrand of phase-space massive neutrino integral
*/
//对p积分，单位Tnu。m的单位也是Tnu。积分结果需要乘T^4，变成T的单位K
static double density_integrand(double p, void *moverT) {
    double m = *((double*)(moverT));
    double fermi=1./(exp(p)+1.); //p in units of Tnu
    double integrand = sqrt(p*p+m*m);
    return fermi * integrand *p*p *2/(2.*M_PI*M_PI);//2 because 2 degrees of freedom
}
static double pressure_integrand(double p, void *moverT) {
    double m = *((double*)(moverT));
    double fermi=1./(exp(p)+1.); //p in units of Tnu
    double integrand = p*p/(3.*sqrt(p*p+m*m));
    return fermi * integrand *p*p *2/(2.*M_PI*M_PI);
}
//中微子质量不变，但温度随红移变化，所以积分结果会变。
static gsl_spline* calculate_density_spline() {
    int N = CLUCU_NU_MNUT_N;
    double *mnut = NULL;
    double *y = NULL;
    gsl_spline* spl = NULL;
    gsl_integration_cquad_workspace * workspace = NULL;
    int stat = 0, gslstatus;
    gsl_function F;

    mnut = clucu_linear_spacing(log(CLUCU_NU_MNUT_MIN), log(CLUCU_NU_MNUT_MAX), N);
    y = malloc(sizeof(double)*CLUCU_NU_MNUT_N);
    if ((y == NULL) || (mnut == NULL))
        CLUCU_RAISE_WARNING(CLUCU_ERROR_NU_INT,"..")

    workspace = gsl_integration_cquad_workspace_alloc(GSL_N_ITERATION_NU);
    if (workspace == NULL)
        CLUCU_RAISE_WARNING(CLUCU_ERROR_NU_INT,"..")

    F.function = &density_integrand;
    for (int i=0; i < CLUCU_NU_MNUT_N; i++) {
        double mnut_ = exp(mnut[i]);
        F.params = &(mnut_);
        gslstatus = gsl_integration_cquad(&F, 0, 1000.0,
                                        GSL_EPSABS_NU,
                                        GSL_EPSREL_NU,
                                        workspace, &y[i], NULL, NULL);
        if (gslstatus != GSL_SUCCESS)
            CLUCU_RAISE_GSL_WARNING(gslstatus,"..");
    }

    spl = gsl_spline_alloc(gsl_interp_akima, CLUCU_NU_MNUT_N);
    if (spl == NULL)
        CLUCU_RAISE_WARNING(CLUCU_ERROR_NU_INT,"..");

    stat |= gsl_spline_init(spl, mnut, y, CLUCU_NU_MNUT_N);
    if (stat)
        CLUCU_RAISE_GSL_WARNING(gslstatus,"...");

  gsl_integration_cquad_workspace_free(workspace);
  free(mnut);
  free(y);

  return spl;
}
static gsl_spline* calculate_pressure_spline() {
    int N = CLUCU_NU_MNUT_N;
    double *mnut = NULL;
    double *y = NULL;
    gsl_spline* spl = NULL;
    gsl_integration_cquad_workspace * workspace = NULL;
    int stat = 0, gslstatus;
    gsl_function F;

    mnut = clucu_linear_spacing(log(CLUCU_NU_MNUT_MIN), log(CLUCU_NU_MNUT_MAX), N);
    y = malloc(sizeof(double)*CLUCU_NU_MNUT_N);
    if ((y == NULL) || (mnut == NULL))
        CLUCU_RAISE_WARNING(CLUCU_ERROR_NU_INT,"..")

    workspace = gsl_integration_cquad_workspace_alloc(GSL_N_ITERATION_NU);
    if (workspace == NULL)
        CLUCU_RAISE_WARNING(CLUCU_ERROR_NU_INT,"..")

    F.function = &pressure_integrand;
    for (int i=0; i < CLUCU_NU_MNUT_N; i++) {
        double mnut_ = exp(mnut[i]);
        F.params = &(mnut_);
        gslstatus = gsl_integration_cquad(&F, 0, 1000.0,
                                        GSL_EPSABS_NU,
                                        GSL_EPSREL_NU,
                                        workspace, &y[i], NULL, NULL);
        if (gslstatus != GSL_SUCCESS)
            CLUCU_RAISE_GSL_WARNING(gslstatus,"..");
    }

    spl = gsl_spline_alloc(gsl_interp_akima, CLUCU_NU_MNUT_N);
    if (spl == NULL)
        CLUCU_RAISE_WARNING(CLUCU_ERROR_NU_INT,"..");

    stat |= gsl_spline_init(spl, mnut, y, CLUCU_NU_MNUT_N);
    if (stat)
        CLUCU_RAISE_GSL_WARNING(gslstatus,"...");

  gsl_integration_cquad_workspace_free(workspace);
  free(mnut);
  free(y);

  return spl;
}
//mass单位eV，Temp单位K，result单位=M_sun/Mpc^3
double density_WDM(double mass, double Temp){
//calculates rhobar(mass,Temp) in Kelvin^4 as in Eq. (19) of Loverde14.
//for a WDM particle of mass mass (in eV) and temperature Temp (in K)
    // Check if the global variable for the phasespace spline has been defined yet:
    if (density_spline == NULL)
        density_spline = calculate_density_spline();

	//double moverT=mass/(Temp*clucu_constants.K_TO_EV); //mnu in units of Tnu.
    double moverT=mass / Temp * (clucu_constants.EV_IN_J / (clucu_constants.KBOLTZ));

    double result;
    int gslstatus = gsl_spline_eval_e(density_spline, log(moverT), NULL, &result);
    if (gslstatus != GSL_SUCCESS)
        CLUCU_RAISE_GSL_WARNING(gslstatus,"nu_phasespace_intg():");

    result *= pow(Temp,4.); //in K^4
    result *= 60./M_PI/M_PI* clucu_constants.STBOLTZ/pow(clucu_constants.CLIGHT,3.); //in kg/m^3
    result *= 1./clucu_constants.SOLAR_MASS * pow(clucu_constants.MPC_TO_METER,3.);//in M_sun/Mpc^3
	return result;
}
//mass单位eV，Temp单位K，result单位=M_sun/Mpc^3
double pressure_WDM(double mass, double Temp){
//calculates rhobar(mass,Temp) in Kelvin^4 as in Eq. (19) of Loverde14.
//for a WDM particle of mass mass (in eV) and temperature Temp (in K)

	if (pressure_spline == NULL)
        pressure_spline = calculate_pressure_spline();

	double moverT=mass/(Temp*clucu_constants.K_TO_EV); //mnu in units of Tnu.
    
    double result;
    int gslstatus = gsl_spline_eval_e(pressure_spline, log(moverT), NULL, &result);
    if (gslstatus != GSL_SUCCESS)
        CLUCU_RAISE_GSL_WARNING(gslstatus,"nu_phasespace_intg():");

    result *= pow(Temp,4.); //in K^4
    result *= 60./M_PI/M_PI* clucu_constants.STBOLTZ/pow(clucu_constants.CLIGHT,3.); //in kg/m^3
    result *= 1./clucu_constants.SOLAR_MASS * pow(clucu_constants.MPC_TO_METER,3.);//in M_sun/Mpc^3
	return result;
}
double EoS_WDM(double mass, double Temp){
//calculates w, Equation of state, of WDM with mass and Temp, (in eV and K, respectively)

	double w = pressure_WDM(mass,Temp)/density_WDM(mass,Temp);

	return w;
}

double clucu_neutrino_mass_split(double sumnu, clucu_neutrino_mass_split_label label, double *mnu, int *N_nu_mass, int *status){
    if(label == clucu_nu_normal){
        *N_nu_mass = 3;
        mnu[0] = (
        2./3.* sumnu - 1./6. * pow(-6. * clucu_constants.DELTAM12_sq + 12. * clucu_constants.DELTAM13_sq_pos + 4. * sumnu*sumnu, 0.5)
        - 0.25 * clucu_constants.DELTAM12_sq / (2./3.* sumnu - 1./6. * pow(-6. * clucu_constants.DELTAM12_sq + 12. *
        clucu_constants.DELTAM13_sq_pos + 4. * sumnu*sumnu, 0.5)));
      mnu[1] = (
        2./3.* sumnu - 1./6. * pow(-6. * clucu_constants.DELTAM12_sq + 12. * clucu_constants.DELTAM13_sq_pos + 4. * sumnu*sumnu, 0.5)
        + 0.25 * clucu_constants.DELTAM12_sq / (2./3.* sumnu - 1./6. * pow(-6. * clucu_constants.DELTAM12_sq + 12. *
        clucu_constants.DELTAM13_sq_pos + 4. * sumnu*sumnu, 0.5)));
      mnu[2] = (
        -1./3. * sumnu + 1./3 * pow(-6. * clucu_constants.DELTAM12_sq +
          12. * clucu_constants.DELTAM13_sq_pos + 4. * sumnu*sumnu, 0.5));

      if (mnu[0] < 0 || mnu[1] < 0 || mnu[2] < 0) {
        // The user has provided a sum that is below the physical limit.
        if (sumnu < 1e-14) {
          mnu[0] = 0.;
          mnu[1] = 0.;
          mnu[2] = 0.;
        }
        else {
          *status = CLUCU_ERROR_MNU_UNPHYSICAL;
          CLUCU_RAISE_WARNING(*status,"Sum of neutrinos masses for this Omeganu "
                    "value is incompatible with the normal mass hierarchy.");
        }
      }
    }
    else if (label == clucu_nu_inverted){
        *N_nu_mass = 3;
        mnu[0] = (
        2./3.* sumnu - 1./6. * pow(-6. * clucu_constants.DELTAM12_sq + 12. * clucu_constants.DELTAM13_sq_neg + 4. * sumnu * sumnu, 0.5)
        - 0.25 * clucu_constants.DELTAM12_sq / (2./3.* sumnu - 1./6. * pow(-6. * clucu_constants.DELTAM12_sq + 12. *
        clucu_constants.DELTAM13_sq_neg + 4. * sumnu * sumnu, 0.5)));
      mnu[1] = (
        2./3.* sumnu - 1./6. * pow(-6. * clucu_constants.DELTAM12_sq + 12. * clucu_constants.DELTAM13_sq_neg + 4. * sumnu*sumnu, 0.5)
        + 0.25 * clucu_constants.DELTAM12_sq / (2./3.* sumnu - 1./6. * pow(-6. *
        clucu_constants.DELTAM12_sq + 12. * clucu_constants.DELTAM13_sq_neg + 4. * sumnu*sumnu, 0.5)));
      mnu[2] = (
        -1./3. * sumnu + 1./3 * pow(-6. * clucu_constants.DELTAM12_sq +
          12. * clucu_constants.DELTAM13_sq_neg + 4. * sumnu*sumnu, 0.5));

      if (mnu[0] < 0 || mnu[1] < 0 || mnu[2] < 0) {
        // The user has provided a sum that is below the physical limit.
        if (sumnu < 1e-14) {
          mnu[0] = 0.;
          mnu[1] = 0.;
          mnu[2] = 0.;;
        }
        else {
          *status = CLUCU_ERROR_MNU_UNPHYSICAL;
          CLUCU_RAISE_WARNING(*status,"Sum of neutrinos masses for this Omeganu "
                    "value is incompatible with the inverted mass hierarchy.");
        }
      }
    }
    else if (label == clucu_nu_equal){
        *N_nu_mass = 3;
        mnu[0] = sumnu/3.;
        mnu[1] = sumnu/3.;
        mnu[2] = sumnu/3.;
    }
    else if (label == clucu_nu_single){
        *N_nu_mass = 1;
        mnu[0] = sumnu;
    }
    else
        CLUCU_RAISE_WARNING(1,"wrong in neutrino");

}

//定义为rho_nu(z)/rho_crit(0), 单位h^2
// double clucu_Omeganuh2(double mass, double Temp)
// {
//     double rho_crit0 = clucu_constants.RHO_CRITICAL;//h^-2 M_sun/Mpc^3
//     double rho_nuz = density_WDM(mass, Temp);//M_sun/Mpc^3

//     return rho_nuz/rho_crit0;
// }

#undef GSL_EPSABS_NU
#undef GSL_EPSREL_NU
#undef GSL_N_ITERATION_NU
