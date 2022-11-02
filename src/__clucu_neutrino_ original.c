#include "clucu_neutrino.h"
//TODO:改为对分布函数的积分


// Global variable to hold the neutrino phase-space spline
gsl_spline* nu_spline = NULL;

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
static double nu_integrand(double x, void *r) {
  double rat = *((double*)(r));
  double x2 = x*x;
  return sqrt(x2 + rat*rat) / (exp(x)+1.0) * x2;
}

/* ------- ROUTINE: clucu_calculate_nu_phasespace_spline ------
TASK: Get the spline of the result of the phase-space integral required for massive neutrinos.
*/

static gsl_spline* calculate_nu_phasespace_spline() {
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

    F.function = &nu_integrand;
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

    //y[0]=(pi^4/15) * (7/8)
    //即为 int dx x^3 / (exp(x)+1.0)的积分结果
    double renorm = 1./y[0];
    for (int i=0; i < CLUCU_NU_MNUT_N; i++)
        y[i] *= renorm;


    spl = gsl_spline_alloc(gsl_interp_akima, CLUCU_NU_MNUT_N);
    if (spl == NULL)
        CLUCU_RAISE_WARNING(CLUCU_ERROR_NU_INT,"..");

    stat |= gsl_spline_init(spl, mnut, y, CLUCU_NU_MNUT_N);
    if (stat)
        CLUCU_RAISE_GSL_WARNING(gslstatus," calculate_nu_phasespace_spline():");

  gsl_integration_cquad_workspace_free(workspace);
  free(mnut);
  free(y);

  return spl;
}

/* ------- ROUTINE: clucu_nu_phasespace_intg ------
INPUTS: mnuOT: the dimensionless mass / temperature of a single massive neutrino
TASK: Get the value of the phase space integral at mnuOT
*/
double nu_phasespace_intg(double mnuOT, int *status)
{
  // Check if the global variable for the phasespace spline has been defined yet:
  if (nu_spline == NULL)
    nu_spline = calculate_nu_phasespace_spline();

  double integral_value = 0.;

  // First check the cases where we are in the limits.
  if (mnuOT < CLUCU_NU_MNUT_MIN)
    return 7./8.;
  else if (mnuOT > CLUCU_NU_MNUT_MAX)
    return 0.2776566337 * mnuOT;

  int gslstatus = gsl_spline_eval_e(nu_spline, log(mnuOT), NULL, &integral_value);
  if (gslstatus != GSL_SUCCESS)
    CLUCU_RAISE_GSL_WARNING(gslstatus,"nu_phasespace_intg():");

    return integral_value * 7./8.;
}
//TODO：由于0.00017这个限制我没看懂，所以我现在mnu只有一个元素，N_nu_mass只为1
/* -------- ROUTINE: Omeganuh2 ---------
INPUTS: a: scale factor, Nnumass: number of massive neutrino species,
        mnu: total mass in eV of neutrinos, T_CMB: CMB temperature, status: pointer to status integer.
TASK: Compute Omeganu * h^2 as a function of time.
!! To all practical purposes, Neff is simply N_nu_mass !!
*/
double clucu_Omeganuh2(double z, int N_nu_mass, double* mnu, double T_CMB, int *status)
{
  double a=1./(1.+z);
  double Tnu, a4, prefix_massless, OmNuh2;
  double Tnu_eff, mnuOT, intval, prefix_massive;

  // First check if N_nu_mass is 0
  if (N_nu_mass == 0) return 0.0;

  Tnu = T_CMB*pow(4./11.,1./3.);
  a4 = a*a*a*a;
  
  // Tnu_eff is used in the massive case because CLASS uses an effective
  // temperature of nonLCDM components to match to mnu / Omeganu =93.14eV. Tnu_eff = T_ncdm * T_CMB = 0.71611 * T_CMB
  Tnu_eff = Tnu * clucu_constants.TNCDM / (pow(4./11.,1./3.));

  // Define the prefix using the effective temperature (to get mnu / Omega = 93.14 eV) for the massive case:
  prefix_massive = NU_CONST * Tnu_eff * Tnu_eff * Tnu_eff * Tnu_eff;

  OmNuh2 = 0.; // Initialize to 0 - we add to this for each massive neutrino species.
  for(int i=0; i < N_nu_mass; i++) {
      
      //TODO: 这个函数不是在求 有质量中微子的密度分数吗？为什么出现无质量中微子
      //mnu的元素应该>0.00017啊，这不应该报错吗
    // Check whether this species is effectively massless
    // In this case, invoke the analytic massless limit:
    if (mnu[i] < 0.00017) {  // Limit taken from Lesgourgues et al. 2012
        prefix_massless = NU_CONST  * Tnu * Tnu * Tnu * Tnu;	
        OmNuh2 = N_nu_mass*prefix_massless*7./8./a4 + OmNuh2;	
        //跟下面的比起来，即
        //Tnu_eff换成了Tnu
        //intval（积分结果*7./8.）换成了1*7./8.
        //求和符号，直接*3

        //循环3次，每次乘3，这不就算重复了吗，这是干嘛

        //当mnuOT < CLUCU_NU_MNUT_MIN时，积分结果为1
    } else {  
       // For the true massive case:  
       // Get mass over T (mass (eV) / ((kb eV/s/K) Tnu_eff (K))
       // This returns the density normalized so that we get nuh2 at a=0
       mnuOT = mnu[i] / (Tnu_eff/a) * (clucu_constants.EV_IN_J / (clucu_constants.KBOLTZ));

       // Get the value of the phase-space integral
       intval = nu_phasespace_intg(mnuOT,status);
       OmNuh2 = intval*prefix_massive/a4 + OmNuh2;
    }
  }

  return OmNuh2;
}

/* -------- ROUTINE: Omeganuh2_to_Mnu ---------
INPUTS: OmNuh2: neutrino mass density today Omeganu * h^2,
        label: how you want to split up the masses, see clucu_neutrinos.h for options,
        T_CMB: CMB temperature, status: pointer to status integer.
TASK: Given Omeganuh2 today, the method of splitting into masses, and the temperature of the
      CMB, output a pointer to the array of neutrino masses (may be length 1 if label asks for sum)
*/
double* clucu_nu_masses(double OmNuh2, clucu_neutrino_mass_splits mass_split,
                      double T_CMB, int *status) {
  double sumnu;
  double *mnu = NULL;

  sumnu = 93.14 * OmNuh2;

    // Now split the sum up into three masses depending on the label given:
    if (mass_split == clucu_nu_normal){
        mnu = malloc(3*sizeof(double));
        if (mnu == NULL)
            CLUCU_RAISE_WARNING(CLUCU_ERROR_MEMORY,"..")

        // See CLUCU note for how we get these expressions for the neutrino masses in
        // normal and inverted hierarchy.
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
        else 
            CLUCU_RAISE_WARNING(CLUCU_ERROR_MNU_UNPHYSICAL,
                                "Sum of neutrinos masses for this Omeganu "
                                "value is incompatible with the requested mass hierarchy.");
        }
    }
    else if (mass_split == clucu_nu_inverted) {
        
        mnu = malloc(3*sizeof(double));
        if (mnu == NULL)
            CLUCU_RAISE_WARNING(CLUCU_ERROR_MEMORY,"..")

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
        else 
            CLUCU_RAISE_WARNING(CLUCU_ERROR_MNU_UNPHYSICAL,
                                "Sum of neutrinos masses for this Omeganu "
                                "value is incompatible with the requested mass hierarchy.");
        }
    }
    else if (mass_split == clucu_nu_equal) {
        mnu = malloc(3*sizeof(double));
        if (mnu == NULL)
            CLUCU_RAISE_WARNING(CLUCU_ERROR_MEMORY,"..")

        mnu[0] = sumnu/3.;
        mnu[1] = sumnu/3.;
        mnu[2] = sumnu/3.;
    }
    else if ((mass_split == clucu_nu_sum) || (mass_split == clucu_nu_single)) {
        mnu = malloc(sizeof(double));
        if (mnu == NULL)
            CLUCU_RAISE_WARNING(CLUCU_ERROR_MEMORY,"..")

        mnu[0] = sumnu;
    }
    else {
        CLUCU_RAISE_WARNING(CLUCU_ERROR_MNU_UNPHYSICAL,
                    "mass option = %d not yet supported!", mass_split);
    }

    return mnu;
}
//在python包实际输出时，mnuz只会有非零值，即长度只能是3('normal', 'inverted', 'equal')，1('sum', 'single')



#undef GSL_EPSABS_NU
#undef GSL_EPSREL_NU
#undef GSL_N_ITERATION_NU
