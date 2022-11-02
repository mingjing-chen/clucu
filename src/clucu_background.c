#include "clucu_background.h"

/*-------------------------------------------------------------------------------------------------------*/
/*                                               计算组分的密度                                            */
/*-------------------------------------------------------------------------------------------------------*/
/* --------- ROUTINE: h_over_h0 ---------
INPUT: scale factor, cosmology
TASK: Compute E(a)=H(a)/H0
*/
static double h_over_h0(double a, clucu_cosmology *cosmo, int *status)
{
  return sqrt(
    (cosmo->param.Omega_m +
     cosmo->param.Omega_l * pow(a,-3*(cosmo->param.w0+cosmo->param.wa)) *exp(3*cosmo->param.wa*(a-1)) +
     cosmo->param.Omega_k * a +
     (cosmo->param.Omega_g + cosmo->param.Omega_nu_rel) / a 
     ) / (a*a*a));
}
/* --------- ROUTINE: clucu_omega_x ---------
INPUT: cosmology object, scale factor, species label
TASK: Compute the density relative to critical, Omega(a) for a given species.
Possible values for "label":
clucu_species_crit_label <- critical (physical)
clucu_species_m_label <- matter
clucu_species_l_label <- DE
clucu_species_g_label <- radiation
clucu_species_k_label <- curvature
clucu_species_ur_label <- massless neutrinos
clucu_species_nu_label <- massive neutrinos
*/
double clucu_omega_x(clucu_cosmology *cosmo, double z, clucu_species_x_label label, int *status)
{
    double a=1./(1.+z);

    double hnorm = h_over_h0(a, cosmo, status);

    switch(label)
    {
        case clucu_species_crit_label :
            return 1.;
        case clucu_species_cb_label :
            return (cosmo->param.Omega_c+cosmo->param.Omega_b)/ (a*a*a) / hnorm / hnorm ;
        case clucu_species_m_label :
            return cosmo->param.Omega_m/ (a*a*a) / hnorm / hnorm ;
        case clucu_species_l_label :
            return
            cosmo->param.Omega_l *
            pow(a,-3 * (1 + cosmo->param.w0 + cosmo->param.wa)) *
            exp(3 * cosmo->param.wa * (a-1)) / hnorm / hnorm;
        case clucu_species_g_label :
            return cosmo->param.Omega_g / (a*a*a*a) / hnorm / hnorm;
        case clucu_species_k_label :
            return cosmo->param.Omega_k / (a*a) / hnorm / hnorm;
        case clucu_species_ur_label :
            return cosmo->param.Omega_nu_rel / (a*a*a*a) / hnorm / hnorm;
        case clucu_species_nu_label :
            return cosmo->param.Omega_nu_mass_sum / (a*a*a) / hnorm / hnorm;
        default:
            *status = CLUCU_ERROR_PARAMETERS;
            CLUCU_RAISE_WARNING(*status,"Species %d not supported", label);
            clucu_cosmology_set_status_message(
            cosmo, "clucu_background.c: clucu_omega_x(): Species %d not supported", label);
            return NAN;
    }
}

/* --------- ROUTINE: clucu_rho_x ---------
INPUT: cosmology object, scale factor, species label
TASK: Compute rho_x(a), with x defined by species label.
Possible values for "label":
clucu_species_crit_label <- critical (physical)
clucu_species_m_label <- matter (physical)
clucu_species_l_label <- DE (physical)
clucu_species_g_label <- radiation (physical)
clucu_species_k_label <- curvature (physical)
clucu_species_ur_label <- massless neutrinos (physical)
clucu_species_nu_label <- massive neutrinos (physical)
*/
//单位: M_sun/Mpc^3
double clucu_rho_x(clucu_cosmology *cosmo, double z, clucu_species_x_label label, int is_comoving, int *status)
{
    double a=1./(1.+z);
    double comfac;
    if (is_comoving)
        comfac = a*a*a;
    else
        comfac = 1.0;
    double hnorm = h_over_h0(a, cosmo, status);
    //rho_cz，不是共动的，我们求的是物理的
    double rhocrit =
    clucu_constants.RHO_CRITICAL *
    (cosmo->param.h) *
    (cosmo->param.h) * hnorm * hnorm * comfac;
    //单位: M_sun/Mpc^3

    return rhocrit * clucu_omega_x(cosmo, z, label, status);
}

/*-------------------------------------------------------------------------------------------------------*/
/*                                               计算chi并插值                                           */
/*-------------------------------------------------------------------------------------------------------*/
// Structure to hold parameters of chi_integrand
typedef struct {
    clucu_cosmology *cosmo;
    int * status;
} chipar;
/* --------- ROUTINE: chi_integrand ---------
INPUT: scale factor
TASK: compute the integrand of the comoving distance
*/
//单位Mpc/h
static double chi_integrand(double a, void * param_void)
{
    clucu_cosmology * cosmo = ((chipar *)param_void)->cosmo;
    int *status = ((chipar *)param_void)->status;

    return clucu_constants.CLIGHT_HMPC/(a*a*h_over_h0(a, cosmo, status));
}

/* --------- ROUTINE: compute_chi ---------
INPUT: scale factor, cosmology
OUTPUT: chi -> radial comoving distance
TASK: compute radial comoving distance at a
*/
//单位Mpc/h
void compute_chi(double a, clucu_cosmology *cosmo, double * chi, int * stat)
{
    int gslstatus;
    double result;
    chipar p;

    p.cosmo=cosmo;
    p.status=stat;

    gsl_integration_cquad_workspace * workspace = NULL;

    gsl_function F;
    F.function = &chi_integrand;
    F.params = &p;

    workspace = gsl_integration_cquad_workspace_alloc(cosmo->gsl_param.N_ITERATION);

    if (workspace == NULL)
    {
        *stat = CLUCU_ERROR_MEMORY;
        CLUCU_RAISE_WARNING(*stat,NULL);
    }
    else
    {
        //TODO: CQUAD is great, but slower than other methods. This could be sped up if it becomes an issue.
        gslstatus=gsl_integration_cquad(
            &F, a, 1.0, 0.0, cosmo->gsl_param.INTEGRATION_DISTANCE_EPSREL, workspace, &result, NULL, NULL);
        *chi=result/cosmo->param.h;
        //单位Mpc/h

        if (gslstatus != GSL_SUCCESS)
        {
            *stat = CLUCU_ERROR_COMPUTECHI;
            CLUCU_RAISE_GSL_WARNING(gslstatus,NULL);
        }
    }

    gsl_integration_cquad_workspace_free(workspace);
}
/*----------------------*/
/*        a(chi)      */
/*----------------------*/
//Root finding for a(chi)
typedef struct {
    double chi;
    clucu_cosmology *cosmo;
    int * status;
} Fpar;

static double fzero(double a,void *param)
{
    double chi,chia,a_use=a;

    chi=((Fpar *)param)->chi;
    compute_chi(a_use,((Fpar *)param)->cosmo,&chia, ((Fpar *)param)->status);

    return chi-chia;
}

static double dfzero(double a,void *param)
{
    clucu_cosmology *cosmo=((Fpar *)param)->cosmo;
    int *stat = ((Fpar *)param)->status;

    chipar p;
    p.cosmo=cosmo;
    p.status=stat;
    //单位Mpc
    return chi_integrand(a,&p)/cosmo->param.h;
}

static void fdfzero(double a,void *param,double *f,double *df)
{
    *f=fzero(a,param);
    *df=dfzero(a,param);
}

/* --------- ROUTINE: a_of_chi ---------
INPUT: comoving distance chi, cosmology, stat, a_old, gsl_root_fdfsolver
OUTPUT: scale factor
TASK: compute the scale factor that corresponds to a given comoving distance chi
Note: This routine uses a root solver to find an a such that compute_chi(a) = chi.
The root solver uses the derivative of compute_chi (which is chi_integrand) and
the value itself.
*/
static void a_of_chi(double chi, clucu_cosmology *cosmo, int* stat, double *a_old, gsl_root_fdfsolver *s)
{
    if(chi==0)
        *a_old=1;
    else
    {
        Fpar p;
        gsl_function_fdf FDF;
        double a_previous,a_current=*a_old;

        p.cosmo=cosmo;
        p.chi=chi;
        p.status=stat;
        FDF.f=&fzero;
        FDF.df=&dfzero;
        FDF.fdf=&fdfzero;
        FDF.params=&p;
        gsl_root_fdfsolver_set(s,&FDF,a_current);

        int iter=0, gslstatus;
        do
        {
            iter++;
            gslstatus=gsl_root_fdfsolver_iterate(s);
            if(gslstatus!=GSL_SUCCESS) CLUCU_RAISE_GSL_WARNING(gslstatus,NULL);
            a_previous=a_current;
            a_current=gsl_root_fdfsolver_root(s);
            gslstatus=gsl_root_test_delta(a_current, a_previous, 0, cosmo->gsl_param.ROOT_EPSREL);
        } while(gslstatus==GSL_CONTINUE && iter <= cosmo->gsl_param.ROOT_N_ITERATION);

        *a_old=a_current;

        // Allows us to pass a status to h_over_h0 for the neutrino integral calculation.
        if(gslstatus==GSL_SUCCESS)
            *stat = *(p.status);
        else 
        {
            *stat = CLUCU_ERROR_COMPUTECHI;
            CLUCU_RAISE_WARNING(*stat,NULL);
        }
    }
}
/*----------------------------*/
/*        E(z),chi(z)插值      */
/*----------------------------*/
/* ----- ROUTINE: clucu_cosmology_compute_distances ------
INPUT: cosmology
TASK: if not already there, make a table of comoving distances and of E(a)
*/
void clucu_cosmology_compute_distances(clucu_cosmology *cosmo, int *status)
{
    //Do nothing if everything is computed already
    if(cosmo->computed_distances)
    return;

    // Create logarithmically and then linearly-spaced values of the scale factor
    int nz = cosmo->spline_param.Z_SPLINE_NLIN_BG+cosmo->spline_param.Z_SPLINE_NLOG_BG-1;
    double *z = clucu_linlog_spacing(cosmo->spline_param.Z_SPLINE_MIN,
                                    cosmo->spline_param.Z_SPLINE_MINLOG_BG,
                                    cosmo->spline_param.Z_SPLINE_MAX_BG,
                                    cosmo->spline_param.Z_SPLINE_NLIN_BG,
                                    cosmo->spline_param.Z_SPLINE_NLOG_BG);

    // Allocate arrays for all three of E(a), chi(a), and a(chi)
    double *E_z = malloc(sizeof(double)*nz);
    double *chi_z = malloc(sizeof(double)*nz);
    // Allocate E(z) and chi(z) and a(chi) splines
    gsl_spline * E = gsl_spline_alloc(gsl_interp_cspline, nz);
    gsl_spline * chi = gsl_spline_alloc(gsl_interp_cspline, nz);
    gsl_spline * zchi = gsl_spline_alloc(gsl_interp_cspline, nz);

    //Check for too little memory
    if (z == NULL || E_z == NULL || chi_z == NULL || E == NULL || chi == NULL || zchi == NULL)
    {
        *status=CLUCU_ERROR_MEMORY;
        CLUCU_RAISE_WARNING(*status,"ran out of memory");
        clucu_cosmology_set_status_message(
            cosmo, "clucu_background.c: clucu_cosmology_compute_distances(): ran out of memory");
    }


    // Fill in E(z) - note, this step cannot change the status variable
    if (!*status)
    {
        for (int i=0; i<nz; i++)
            E_z[i] = h_over_h0(1./(1.+z[i]), cosmo, status);
    }

    // Create a E(z) spline
    if (!*status)
    {
        if (gsl_spline_init(E, z, E_z, nz))
        {
            *status = CLUCU_ERROR_SPLINE;
            CLUCU_RAISE_WARNING(*status,"Error creating  E(z) spline");
            clucu_cosmology_set_status_message(
            cosmo, "clucu_background.c: clucu_cosmology_compute_distances(): Error creating  E(z) spline");
        }
    }

    // Compute chi(z)
    if (!*status)
    {
        for (int i=0; i<nz; i++)
            compute_chi(1./(1.+z[i]), cosmo, &chi_z[i], status);
        if (*status)
        {
            *status = CLUCU_ERROR_INTEG;
            CLUCU_RAISE_WARNING(*status,"chi(z) integration error ");
            clucu_cosmology_set_status_message(
            cosmo, "clucu_background.c: clucu_cosmology_compute_distances(): chi(z) integration error ");
        }
    }

    // Initialize chi(z) spline
    if (!*status)
    {
        if (gsl_spline_init(chi, z, chi_z, nz))
        {//in Mpc
            *status = CLUCU_ERROR_SPLINE;
            CLUCU_RAISE_WARNING(*status,"Error creating  chi(z) spline");
            clucu_cosmology_set_status_message(
            cosmo, "clucu_background.c: clucu_cosmology_compute_distances(): Error creating  chi(z) spline");
        }
    }

    // Initialize z(chi) spline
    if (!*status)
    {
        if (gsl_spline_init(zchi, chi_z, z, nz))
        {//in Mpc
            *status = CLUCU_ERROR_SPLINE;
            CLUCU_RAISE_WARNING(*status,"Error creating  z(chi) spline");
            clucu_cosmology_set_status_message(
            cosmo, "clucu_background.c: clucu_cosmology_compute_distances(): Error creating  z(chi) spline");
        }
    }

    if (*status)
    { //If there was an error, free the GSL splines and return
        gsl_spline_free(E); // Note: you are allowed to call gsl_free() on NULL
        gsl_spline_free(chi);
        gsl_spline_free(zchi);
        E = NULL;
        chi = NULL;
        zchi = NULL;
    }

    //TODO: The interval in chi (5. Mpc) should be made a macro
    free(z); 
    free(E_z);
    free(chi_z);
    z = NULL;
    E_z = NULL;
    chi_z = NULL;
    //Note: you are allowed to call free() on NULL

    if (*status == 0)
    {
        //If there were no errors, attach the splines to the cosmo struct and end the function.
        cosmo->data.E             = E;
        cosmo->data.chi           = chi;
        cosmo->data.zchi          = zchi;
        cosmo->computed_distances = true;
    }
}

/*-------------------------------------------------------------------------------------------------------*/
/*                                               计算D(z)并插值                                           */
/*-------------------------------------------------------------------------------------------------------*/
/* --------- ROUTINE: growth_ode_system ---------
INPUT: scale factor
TASK: Define the ODE system to be solved in order to compute the growth (of the density)
*/
static int growth_ode_system(double a,const double y[],double dydt[],void *param)
{
    int status = 0;
    clucu_cosmology * cosmo = param;

    double hnorm=h_over_h0(a,cosmo, &status);
    double om=clucu_omega_x(cosmo, a, clucu_species_m_label, &status);

    dydt[1]=1.5*hnorm*a*om*y[0];
    dydt[0]=y[1]/(a*a*a*hnorm);

    return status;
}

/* --------- ROUTINE: growth_factor_and_growth_rate ---------
INPUT: scale factor, cosmology
TASK: compute the growth (D(z)) and the growth rate, logarithmic derivative (f?)
*/
static int growth_factor_and_growth_rate(double a, double *gf, double *fg, clucu_cosmology *cosmo, int *stat)
{
    if(a < cosmo->gsl_param.EPS_SCALEFAC_GROWTH)
    {
        *gf = a;
        *fg = 1;
        return 0;
    }
    else
    {
        int gslstatus;
        double y[2];
        double ainit = cosmo->gsl_param.EPS_SCALEFAC_GROWTH;


        gsl_odeiv2_system sys = {growth_ode_system, NULL, 2, cosmo};

        gsl_odeiv2_driver *d = gsl_odeiv2_driver_alloc_y_new(
        &sys, gsl_odeiv2_step_rkck,
        0.1*cosmo->gsl_param.EPS_SCALEFAC_GROWTH, 0, cosmo->gsl_param.ODE_GROWTH_EPSREL);

        if (d == NULL)
            return CLUCU_ERROR_MEMORY;

        y[0] = cosmo->gsl_param.EPS_SCALEFAC_GROWTH;
        y[1] = (
        cosmo->gsl_param.EPS_SCALEFAC_GROWTH *
        cosmo->gsl_param.EPS_SCALEFAC_GROWTH *
        cosmo->gsl_param.EPS_SCALEFAC_GROWTH*
        h_over_h0(cosmo->gsl_param.EPS_SCALEFAC_GROWTH, cosmo, stat));

        gslstatus = gsl_odeiv2_driver_apply(d, &ainit, a, y);
        gsl_odeiv2_driver_free(d);

        if(gslstatus != GSL_SUCCESS)
        {
            CLUCU_RAISE_GSL_WARNING(gslstatus,NULL);
            return gslstatus;
        }

        *gf = y[0];
        *fg = y[1]/(a*a*h_over_h0(a, cosmo, stat)*y[0]);

        return 0;
    }
}
/* ----- ROUTINE: clucu_cosmology_compute_growth ------
INPUT: cosmology
TASK: if not already there, make a table of growth function and growth rate
      normalize growth to input parameter growth0
*/
void clucu_cosmology_compute_growth(clucu_cosmology *cosmo, int* status)
{
    if (cosmo->computed_growth)
    return;

    // Create logarithmically and then linearly-spaced values of the scale factor
    int chistatus = 0;
    int nz = cosmo->spline_param.Z_SPLINE_NLIN_BG+cosmo->spline_param.Z_SPLINE_NLOG_BG-1;
    double *z = clucu_linlog_spacing(cosmo->spline_param.Z_SPLINE_MIN,
                                    cosmo->spline_param.Z_SPLINE_MINLOG_BG,
                                    cosmo->spline_param.Z_SPLINE_MAX_BG,
                                    cosmo->spline_param.Z_SPLINE_NLIN_BG,
                                    cosmo->spline_param.Z_SPLINE_NLOG_BG);
    gsl_integration_cquad_workspace * workspace = NULL;
    gsl_function F;
    double growth0, fgrowth0;
    double *y = NULL;
    double *y2 = NULL;
    double df, integ;


    // allocate space for y, which will be all three
    // of E(a), chi(a), D(a) and f(a) in turn.
    if (*status == 0)
    {
        y = malloc(sizeof(double)*nz);
        if (y == NULL)
        {
            *status = CLUCU_ERROR_MEMORY;
            CLUCU_RAISE_WARNING(*status,"ran out of memory");
            clucu_cosmology_set_status_message(
                cosmo, "clucu_background.c: clucu_cosmology_compute_growth(): ran out of memory");
        }
    }

    if (*status == 0)
    {
        y2 = malloc(sizeof(double)*nz);
        if (y2 == NULL) {
            *status = CLUCU_ERROR_MEMORY;
            CLUCU_RAISE_WARNING(*status,"ran out of memory");
            clucu_cosmology_set_status_message(
                cosmo, "clucu_background.c: clucu_cosmology_compute_growth(): ran out of memory");
        }
    }

    if (*status == 0)
    {
        // Get the growth factor and growth rate at z=0
        chistatus |= growth_factor_and_growth_rate(1., &growth0, &fgrowth0, cosmo, status);

        // Get the growth factor and growth rate at other redshifts
        for(int i=0; i<nz; i++)
        {
            chistatus |= growth_factor_and_growth_rate(1./(1.+z[i]), &(y[i]), &(y2[i]), cosmo, status);
            // Normalizing to the growth factor to the growth today
            y[i] /= growth0;
        }

        if (chistatus  || *status)
        {
            if (chistatus)
            {
                *status = CLUCU_ERROR_INTEG;
                CLUCU_RAISE_WARNING(*status,"integral for linear growth factor didn't converge");
                clucu_cosmology_set_status_message(
                    cosmo, "clucu_background.c: clucu_cosmology_compute_growth(): integral for linear growth factor didn't converge");
            }
        }
    }

    gsl_spline *growth = NULL;
    gsl_spline *fgrowth = NULL;

    if (*status == 0)
    {
        //*(cosmo->spline_param.A_SPLINE_TYPE)
        growth = gsl_spline_alloc(gsl_interp_cspline, nz);
        fgrowth = gsl_spline_alloc(gsl_interp_cspline, nz);

        if (growth == NULL || fgrowth == NULL)
        {
            *status = CLUCU_ERROR_MEMORY;
            CLUCU_RAISE_WARNING(*status,"ran out of memory");
            clucu_cosmology_set_status_message(
                cosmo, "clucu_background.c: clucu_cosmology_compute_growth(): ran out of memory");
        }
    }

    if (*status == 0)
    {
        chistatus = gsl_spline_init(growth, z, y, nz);

        if (chistatus)
        {
            *status = CLUCU_ERROR_SPLINE;
            CLUCU_RAISE_WARNING(*status,"Error creating D(a) spline");
            clucu_cosmology_set_status_message(
            cosmo, "clucu_background.c: clucu_cosmology_compute_growth(): Error creating D(a) spline");
        }
    }

    if (*status == 0)
    {
        chistatus = gsl_spline_init(fgrowth, z, y2, nz);
        if (chistatus)
        {
            *status = CLUCU_ERROR_SPLINE;
            CLUCU_RAISE_WARNING(*status,"Error creating f(a) spline");
            clucu_cosmology_set_status_message(
            cosmo, "clucu_background.c: clucu_cosmology_compute_growth(): Error creating f(a) spline");
        }
    }

    if (*status == 0)
    {
        // assign all the splines we've just made to the structure.
        cosmo->data.growth = growth;
        cosmo->data.fgrowth = fgrowth;
        cosmo->data.growth0 = growth0;
        cosmo->computed_growth = true;
    }
    else
    {
        gsl_spline_free(growth);
        gsl_spline_free(fgrowth);
    }

    free(z);
    free(y);
    free(y2);
    gsl_integration_cquad_workspace_free(workspace);
}
/*-------------------------------------------------------------------------------------------------------*/
/*                            根据插值后的chi(z),D(z) ，计算其他背景量DA,DL,f                                 */
/*-------------------------------------------------------------------------------------------------------*/
void clucu_background_cla(clucu_cosmology *cosmo,int *status)
{
    clucu_cosmology_compute_distances(cosmo,status);
    CLUCU_CHECK_ERROR( );
    clucu_cosmology_compute_growth(cosmo,status);
    CLUCU_CHECK_ERROR( );
}
/*-------------------------------------------------------------------------------------------------------*/
/*                            CLASS                                 */
/*-------------------------------------------------------------------------------------------------------*/
//运行class，得到很多个z的P(k)，和一个背景表格。只读取H(z)。
void RunClass(clucu_cosmology *cosmo)
{
    if(cosmo->runned_class == true) return;
    if(cosmo->class_id != cosmo->cosmo_id){
        cosmo->runned_class = true;
        return;  
    }
	char namedir[50];
    char nameclassini[50];
    char cmd[100];
	sprintf(namedir,"../data/%s", cosmo->classname);
	mkdir("../data", 0755);
	mkdir(namedir, 0755);
    sprintf(nameclassini,"../data/%s/p%d.ini", cosmo->classname, cosmo->class_id);
    sprintf(cmd," ../src/class  ../data/%s/p%d.ini", cosmo->classname, cosmo->class_id);


    //创建ini文件
    FILE* fp1=fopen(nameclassini,"w");
    if(! fp1) 
    {
        cosmo->runned_class = false;
        return;
    }
    
    //基本参数
//	fprintf(fp1,"h = %.15lf  \nT_cmb =%lf  \nomega_b = %lf   \nomega_cdm = %lf    \n",
//                    cosmo->param.h,cosmo->param.T_CMB,  cosmo->param.Obh2,  cosmo->param.Och2);
    fprintf(fp1,"h = %.15lf  \nOmega_b = %lf   \nOmega_cdm = %lf    \n",
                    cosmo->param.h,  cosmo->param.Omega_b,  cosmo->param.Omega_c);
    fprintf(fp1,"Omega_k = %lf    \nT_cmb =%lf  \n",cosmo->param.Omega_k,cosmo->param.T_CMB);
    fprintf(fp1,"\n");
    
    //暗能量参数
    if(isnormal(cosmo->param.w0))
        fprintf(fp1,"Omega_Lambda = 0   \nfluid_equation_of_state = CLP     \nw0_fld = %lf    \nwa_fld = %lf    \n",cosmo->param.w0,cosmo->param.wa);
    else
        fprintf(fp1,"Omega_Lambda = %lf     \n",cosmo->param.Omega_l);
    fprintf(fp1,"\n");
    
    //功率谱参数
    fprintf(fp1,"k_pivot = %lf  \nn_s = %lf    \nalpha_s = %lf    \n",cosmo->param.k_pivot,cosmo->param.n_s,cosmo->param.alpha_s);
    if(!isnormal(cosmo->param.sigma8) && isnormal(cosmo->param.A_s))
        fprintf(fp1,"A_s = %le    \n",cosmo->param.A_s);
    else if(isnormal(cosmo->param.sigma8) && !isnormal(cosmo->param.A_s))
        fprintf(fp1,"sigma8 = %lf    \n",cosmo->param.sigma8);
    fprintf(fp1,"\n");
    
    //中微子参数
    if(cosmo->param.N_nu_mass > 0)
    {
        fprintf(fp1,"N_ur = %lf    \nN_ncdm = %d        \n",cosmo->param.N_nu_rel,cosmo->param.N_nu_mass);
        fprintf(fp1,"m_ncdm = ");
        for(int i=0;i<cosmo->param.N_nu_mass-1;i++)
            fprintf(fp1,"%lf,",cosmo->param.m_nu[i]);
        fprintf(fp1,"%lf\n",cosmo->param.m_nu[cosmo->param.N_nu_mass-1]);
    }
    fprintf(fp1,"\n");
    
    if(cosmo->config.matter_power_spectrum_method == clucu_boltzmann_class_Pk)
        fprintf(fp1,"\noutput = mPk    \n");
    else if(cosmo->config.matter_power_spectrum_method == clucu_boltzmann_class_Tk)
        fprintf(fp1,"\noutput = mTk    \n");
    fprintf(fp1,"P_k_max_1/Mpc = 100.    \nk_per_decade_for_pk = 1000    \nroot = ../data/%s/p%d_    \nwrite background = yes\n",cosmo->classname, cosmo->class_id);
    fprintf(fp1,"z_pk=");
    for(int i=0;i<cosmo->data.spline_Nz-1;i++)
        fprintf(fp1,"%.4f,",cosmo->data.spline_z[i]);
    fprintf(fp1,"%.4f",cosmo->data.spline_z[cosmo->data.spline_Nz-1]);

    fclose(fp1);

    //调用class，得到输出文件。这里只读取H(z)
    system(cmd);
    cosmo->runned_class = true;
    return;
}
//Mpc
void clucu_background_class(clucu_cosmology *cosmo,int *status)
//H,dC,dA,dL,
{   
    if(cosmo->runned_class==false)
    {
        RunClass(cosmo);
    }
    char nameclassbg[150];
    sprintf(nameclassbg,"../data/%s/p%d_background.dat", cosmo->classname, cosmo->class_id);

    FILE* fp1=fopen(nameclassbg,"r");
    if(! fp1)
    {
        *status = CLUCU_ERROR_FILE;
        CLUCU_RAISE_WARNING(*status,"failed in opening \"%s\"",nameclassbg);
        cosmo->computed_distances = false;
        cosmo->computed_growth = false;
        return;
    }
    int nz=1000;
    double *z=Create1Grid(nz);
    double *H_z=Create1Grid(nz);
    double *E_z=Create1Grid(nz);
    double *chi_z=Create1Grid(nz);
    double *growth_z=Create1Grid(nz);
    double *fgrowth_z=Create1Grid(nz);
    //Check for too little memory
    if (z == NULL || H_z == NULL || E_z == NULL || chi_z == NULL || growth_z == NULL || fgrowth_z == NULL)
    {
        *status=CLUCU_ERROR_MEMORY;
        CLUCU_RAISE_WARNING(*status,"ran out of memory");
    }
    for(int table_i=0;table_i<4624;table_i++)
    {
        char abandon[1000];
        int i=nz-(table_i-4624+nz)-1;
        if(table_i<(4624-nz))
            fgets(abandon,1000,fp1);
        else
        {
            if(!isnormal(cosmo->param.w0))//无动态暗能量
            {
                if(cosmo->param.N_nu_mass == 0)
                    fscanf(fp1,"%lf  %*lf  %*lf  %lf  %lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %lf  %lf",&z[i],&H_z[i],&chi_z[i],&growth_z[i],&fgrowth_z[i]);
                else if(cosmo->param.N_nu_mass == 1)
                    fscanf(fp1,"%lf  %*lf  %*lf  %lf  %lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %*lf  %lf  %lf",&z[i],&H_z[i],&chi_z[i],&growth_z[i],&fgrowth_z[i]);
                else if(cosmo->param.N_nu_mass == 2)
                    fscanf(fp1,"%lf  %*lf  %*lf  %lf  %lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %*lf  %*lf  %*lf    %lf  %lf",&z[i],&H_z[i],&chi_z[i],&growth_z[i],&fgrowth_z[i]);
                else if(cosmo->param.N_nu_mass == 3)
                    fscanf(fp1,"%lf  %*lf  %*lf  %lf  %lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %lf  %lf",&z[i],&H_z[i],&chi_z[i],&growth_z[i],&fgrowth_z[i]);
            }
            else//有动态暗能量
            {
                if(cosmo->param.N_nu_mass == 0)
                    fscanf(fp1,"%lf  %*lf  %*lf  %lf  %lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %*lf  %*lf  %*lf    %lf  %lf",&z[i],&H_z[i],&chi_z[i],&growth_z[i],&fgrowth_z[i]);
                else if(cosmo->param.N_nu_mass == 1)
                    fscanf(fp1,"%lf  %*lf  %*lf  %lf  %lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %lf  %lf",&z[i],&H_z[i],&chi_z[i],&growth_z[i],&fgrowth_z[i]);
                else if(cosmo->param.N_nu_mass == 2)
                    fscanf(fp1,"%lf  %*lf  %*lf  %lf  %lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %*lf  %lf  %lf",&z[i],&H_z[i],&chi_z[i],&growth_z[i],&fgrowth_z[i]);
                else if(cosmo->param.N_nu_mass == 3)
                    fscanf(fp1,"%lf  %*lf  %*lf  %lf  %lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %*lf  %*lf  %*lf    %*lf  %*lf  %*lf  %lf  %lf",&z[i],&H_z[i],&chi_z[i],&growth_z[i],&fgrowth_z[i]);
            }
        }
    }
    if( z[0]!=0. || chi_z[0]!=0. || growth_z[0]!=1. )
    {
        *status=CLUCU_ERROR_FILE;
        CLUCU_RAISE_WARNING(*status,"failed in loading background: z[0]=%e,chi[0]=%e,D[0]=%e;z[1]=%e,chi[1]=%e,D[1]=%e",z[0],chi_z[0],growth_z[0],z[1],chi_z[1],growth_z[1]);
    }
    fclose(fp1);
    for(int i=0;i<nz;i++) E_z[i] = H_z[i]/H_z[0];
    double growth0=growth_z[0];

    gsl_spline *E = gsl_spline_alloc(gsl_interp_cspline, nz);
    gsl_spline *chi = gsl_spline_alloc(gsl_interp_cspline, nz);
    gsl_spline *zchi = gsl_spline_alloc(gsl_interp_cspline, nz);
    gsl_spline *growth = gsl_spline_alloc(gsl_interp_cspline, nz);
    gsl_spline *fgrowth = gsl_spline_alloc(gsl_interp_cspline, nz);
    if (E == NULL || chi == NULL || zchi == NULL || growth == NULL || fgrowth == NULL)
    {
        *status=CLUCU_ERROR_MEMORY;
        CLUCU_RAISE_WARNING(*status,"ran out of memory");
    }
 
    // Create a E(z) spline
    if (!*status){
        if (gsl_spline_init(E, z, E_z, nz)){
            *status = CLUCU_ERROR_SPLINE;
            CLUCU_RAISE_WARNING(*status,"Error creating  E(z) spline");
            clucu_cosmology_set_status_message(
            cosmo, "clucu_background_class.c: LoadBackGround(): Error creating  E(z) spline");
        }
    }
    // Initialize chi(z) spline
    if (!*status){
        if (gsl_spline_init(chi, z, chi_z, nz)){//in Mpc
            *status = CLUCU_ERROR_SPLINE;
            CLUCU_RAISE_WARNING(*status,"Error creating  chi(z) spline");
            clucu_cosmology_set_status_message(
            cosmo, "clucu_background.c: clucu_cosmology_compute_distances(): Error creating  chi(z) spline");
        }
    }

    // Initialize z(chi) spline
    if (!*status){
        if (gsl_spline_init(zchi, chi_z, z, nz)){//in Mpc
            *status = CLUCU_ERROR_SPLINE;
            CLUCU_RAISE_WARNING(*status,"Error creating  z(chi) spline");
            clucu_cosmology_set_status_message(
            cosmo, "clucu_background.c: clucu_cosmology_compute_distances(): Error creating  z(chi) spline");
        }
    }

    // Initialize D(z) spline
    if (!*status){
        if (gsl_spline_init(growth, z, growth_z, nz)){//in Mpc
            *status = CLUCU_ERROR_SPLINE;
            CLUCU_RAISE_WARNING(*status,"Error creating  growth(z) spline");
            clucu_cosmology_set_status_message(
            cosmo, "clucu_background.c: clucu_cosmology_compute_distances(): Error creating  z(chi) spline");
        }
    }

    // Initialize f(z) spline
    if (!*status){
        if (gsl_spline_init(fgrowth, z, fgrowth_z, nz)){//in Mpc
            *status = CLUCU_ERROR_SPLINE;
            CLUCU_RAISE_WARNING(*status,"Error creating  fgrowth(z) spline");
            clucu_cosmology_set_status_message(
            cosmo, "clucu_background.c: clucu_cosmology_compute_distances(): Error creating  z(chi) spline");
        }
    }

    if (*status){ //If there was an error, free the GSL splines and return
        gsl_spline_free(E); // Note: you are allowed to call gsl_free() on NULL
        gsl_spline_free(chi);
        gsl_spline_free(zchi);
        gsl_spline_free(growth);
        gsl_spline_free(fgrowth);
        E = NULL;
        chi = NULL;
        zchi = NULL;
        growth = NULL;
        fgrowth = NULL;
    }
    free(z); 
    free(H_z);
    free(E_z);
    free(chi_z);
    free(growth_z);
    free(fgrowth_z);
    z = NULL;
    E_z = NULL;
    chi_z = NULL;
    growth_z = NULL;
    fgrowth_z = NULL;
    if (*status == 0) {
        //If there were no errors, attach the splines to the cosmo struct and end the function.
        cosmo->data.E             = E;
        cosmo->data.chi           = chi;
        cosmo->data.zchi          = zchi;
        cosmo->computed_distances = true;
        cosmo->data.growth = growth;
        cosmo->data.fgrowth = fgrowth;
        cosmo->data.growth0 = growth0;
        cosmo->computed_growth = true;
    }
}

/*-------------------------------------------------------------------------------------------------------*/
/*                            最终                                 */
/*-------------------------------------------------------------------------------------------------------*/

void clucu_compute_background(clucu_cosmology *cosmo,clucu_background_label label,int *status)
{
    switch(label)
    {
        case clucu_background_cla_label: clucu_background_cla(cosmo,status);break;
        case clucu_background_class_label: clucu_background_class(cosmo,status);break;
    }
    CLUCU_CHECK_ERROR( );
}

double clucu_h_over_h0(clucu_cosmology *cosmo,double z,int *status)
{
    if(!cosmo->computed_distances){
        *status = CLUCU_ERROR_DISTANCES_INIT;
        CLUCU_RAISE_WARNING(*status,"distance splines have not been precomputed!");
        clucu_cosmology_set_status_message(
            cosmo,
            "clucu_background.c: clucu_h_over_h0(): distance splines have not been precomputed!");
        return NAN;
    }

    double h_over_h0;
    int gslstatus = gsl_spline_eval_e(cosmo->data.E, z, NULL, &h_over_h0);
    if(gslstatus != GSL_SUCCESS){
        CLUCU_RAISE_GSL_WARNING(gslstatus,"z=%.2e is outside interpolation range [%.2e,%.2e]",z,cosmo->data.E->interp->xmin,cosmo->data.E->interp->xmax);
        *status = gslstatus;
        clucu_cosmology_set_status_message(cosmo, "clucu_background.c: clucu_h_over_h0(): Scale factor outside interpolation range.");
        return NAN;
    }

    return h_over_h0;
}
double clucu_hubble_parameter(clucu_cosmology *cosmo,double z,int *status)
{
    double H = cosmo->param.H0 * 1000.0 / clucu_constants.CLIGHT;	   //单位：1/Mpc
    double Hz = H * clucu_h_over_h0(cosmo,z,status);
    CLUCU_CHECK_ERROR(NAN);

    return Hz;
}
//chi Mpc
double clucu_comoving_radial_distance(clucu_cosmology *cosmo,double z,int *status)
{
    double a=1./(1.+z);
    if((a > (1.0 - 1.e-8)) && (a<=1.0)) {
    return 0.;
    }
    else if(a>1.){
        *status = CLUCU_ERROR_COMPUTECHI;
        CLUCU_RAISE_WARNING(*status,"redshift cannot be lower than 0. z=%e",z);
        clucu_cosmology_set_status_message(cosmo, "clucu_background.c: scale factor cannot be larger than 1.");
        return NAN;
    }
    else{
        if(!cosmo->computed_distances){
            *status = CLUCU_ERROR_DISTANCES_INIT;
            CLUCU_RAISE_WARNING(*status,"distance splines have not been precomputed!");
            clucu_cosmology_set_status_message(
            cosmo,
            "clucu_background.c: clucu_comoving_radial_distance(): distance splines have not been precomputed!");
            return NAN;
        }
    }

    double crd;
    int gslstatus = gsl_spline_eval_e(cosmo->data.chi, z, NULL, &crd);
    if(gslstatus != GSL_SUCCESS){
        CLUCU_RAISE_GSL_WARNING(gslstatus,"z=%.2e is outside interpolation range [%.2e,%.2e]",z,cosmo->data.chi->interp->xmin,cosmo->data.chi->interp->xmax);
        *status = gslstatus;
        clucu_cosmology_set_status_message(
        cosmo, "clucu_background.c: clucu_comoving_radial_distance(): Scale factor outside interpolation range.");
        return NAN;
    }
    return crd;
}
double clucu_sinn(clucu_cosmology *cosmo, double chi,int *status)
{
    //////
    //         { sin(x)  , if k==1
    // sinn(x)={  x      , if k==0
    //         { sinh(x) , if k==-1
    switch(cosmo->param.k_sign){
    case -1:
        return sinh(cosmo->param.sqrtk * chi) / cosmo->param.sqrtk;
    case 1:
        return sin(cosmo->param.sqrtk*chi) / cosmo->param.sqrtk;
    case 0:
        return chi;
    default:
        *status = CLUCU_ERROR_PARAMETERS;
        CLUCU_RAISE_WARNING(*status,"ill-defined cosmo->params.k_sign = %d",cosmo->param.k_sign);
        clucu_cosmology_set_status_message(cosmo, "clucu_background.c: clucu_sinn: ill-defined cosmo->params.k_sign = %d",
                                        cosmo->param.k_sign);
        return NAN;

    }
}
//dC
double clucu_comoving_angular_distance(clucu_cosmology *cosmo,double z,int *status)
{
    double crd = clucu_comoving_radial_distance(cosmo,z,status);
    CLUCU_CHECK_ERROR(NAN);

    double cad = clucu_sinn(cosmo,crd,status);
    CLUCU_CHECK_ERROR(NAN);

    return cad;
}
//dA
double clucu_angular_diameter_distance(clucu_cosmology *cosmo,double z,int *status)
{
    double cad = clucu_comoving_angular_distance(cosmo,z,status);
    CLUCU_CHECK_ERROR(NAN);
    
    return cad/(1.+z);
}
//dL
double clucu_luminosity_distance(clucu_cosmology *cosmo,double z,int *status)
{
    double cad = clucu_comoving_angular_distance(cosmo,z,status);
    CLUCU_CHECK_ERROR(NAN);
    
    return cad*(1.+z);
}
double clucu_growth_factor(clucu_cosmology *cosmo,double z,int *status)
{
    double a=1./(1.+z);
    if(a==1.){
        return 1.;
    }
    else if(a>1.){
        *status = CLUCU_ERROR_COMPUTECHI;
        CLUCU_RAISE_WARNING(*status,"redshift cannot be lower than 0. z=%e",z);
        clucu_cosmology_set_status_message(
            cosmo, "clucu_background.c: scale factor cannot be larger than 1.");
        return NAN;
    }
    else{
        if(!cosmo->computed_growth){
            *status = CLUCU_ERROR_GROWTH_INIT;
            CLUCU_RAISE_WARNING(*status,"growth factor splines have not been precomputed!");
            clucu_cosmology_set_status_message(
            cosmo,
            "clucu_background.c: clucu_growth_factor(): growth factor splines have not been precomputed!");
            return NAN;
        }
        if (*status != CLUCU_ERROR_NOT_IMPLEMENTED) {
            double D;
            int gslstatus = gsl_spline_eval_e(cosmo->data.growth, z, NULL, &D);
            if(gslstatus != GSL_SUCCESS){
                CLUCU_RAISE_GSL_WARNING(gslstatus,"z=%.2e is outside interpolation range [%.2e,%.2e]",z,cosmo->data.growth->interp->xmin,cosmo->data.growth->interp->xmax);
                *status |= gslstatus;
                clucu_cosmology_set_status_message(
                    cosmo, "clucu_background.c: clucu_growth_factor(): Scale factor outside interpolation range.");
                return NAN;
            }
            return D;
        }
        else{
            return NAN;
        }
    }
}
double clucu_growth_rate(clucu_cosmology *cosmo,double z,int *status)
{
    double a=1./(1.+z);
    if(a>1.){
        *status = CLUCU_ERROR_COMPUTECHI;
        CLUCU_RAISE_WARNING(*status,"redshift cannot be lower than 0. z=%e",z);
        clucu_cosmology_set_status_message(cosmo, "clucu_background.c: scale factor cannot be larger than 1.");
        return NAN;
    }
    else{
        if (!cosmo->computed_growth){
            *status = CLUCU_ERROR_GROWTH_INIT;
            CLUCU_RAISE_WARNING(*status,"growth factor splines have not been precomputed!");
            clucu_cosmology_set_status_message(
            cosmo,
            "clucu_background.c: clucu_growth_rate(): growth rate splines have not been precomputed!");
        }
        if(*status != CLUCU_ERROR_NOT_IMPLEMENTED){
            double g;
            int gslstatus = gsl_spline_eval_e(cosmo->data.fgrowth, z, NULL ,&g);
            if(gslstatus != GSL_SUCCESS){
                CLUCU_RAISE_GSL_WARNING(gslstatus,"z=%.2e is outside interpolation range [%.2e,%.2e]",z,cosmo->data.fgrowth->interp->xmin,cosmo->data.fgrowth->interp->xmax);
                *status |= gslstatus;
                clucu_cosmology_set_status_message(cosmo, "clucu_background.c: clucu_growth_rate(): Scale factor outside interpolation range.");
                return NAN;
            }
            return g;
        } else {
            return NAN;
        }
    }
}
double clucu_chi_to_z(clucu_cosmology *cosmo,double chi,int *status)
{
    if(chi<0.){
        *status = CLUCU_ERROR_COMPUTECHI;
        CLUCU_RAISE_WARNING(*status,"comoving distance cannot be lower than 0. chi=%e",chi);
        clucu_cosmology_set_status_message(cosmo, "clucu_background.c: scale factor cannot be larger than 1.");
        return NAN;
    }
    else{
        if (!cosmo->computed_growth){
            *status = CLUCU_ERROR_GROWTH_INIT;
            CLUCU_RAISE_WARNING(*status,"growth factor splines have not been precomputed!");
            clucu_cosmology_set_status_message(
            cosmo,
            "clucu_background.c: clucu_growth_rate(): growth rate splines have not been precomputed!");
        }
        if(*status != CLUCU_ERROR_NOT_IMPLEMENTED){
            double z;
            int gslstatus = gsl_spline_eval_e(cosmo->data.zchi, chi, NULL ,&z);
            if(gslstatus != GSL_SUCCESS){
                CLUCU_RAISE_GSL_WARNING(gslstatus,"chi=%.2e is outside interpolation range [%.2e,%.2e]",chi,cosmo->data.zchi->interp->xmin,cosmo->data.zchi->interp->xmax);
                *status |= gslstatus;
                clucu_cosmology_set_status_message(cosmo, "clucu_background.c: clucu_growth_rate(): Scale factor outside interpolation range.");
                return NAN;
            }
            return z;
        } else {
            return NAN;
        }
    }
}
//共动体积，单位Mpc^3
double clucu_volume_element(clucu_cosmology *cosmo,double z,int *status)
{
    double dC =  clucu_comoving_angular_distance(cosmo,z,status);//共动距离，单位Mpc
    CLUCU_CHECK_ERROR(NAN);
    

    double H = clucu_hubble_parameter(cosmo,z,status);//哈勃参数，单位1/Mpc
    CLUCU_CHECK_ERROR(NAN);

    return dC*dC/H;//共动体积，单位Mpc^3
}