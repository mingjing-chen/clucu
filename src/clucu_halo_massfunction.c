#include "clucu_halo_massfunction.h"
#include "clucu_halo_bias.h"
/*----------------------*/
/*        拟合函数       */
/*----------------------*/
//决定了halo_define_method
//但concentration可以灵活选择halo_define，以转换质量
//拟合函数，PS
double MassFuncPress1974(double sigma,int *status)
{
    double delta_c = 1.68647;
    double nu = delta_c/sigma;
	return sqrt(2.0 / M_PI) * nu * exp(-0.5 * nu*nu);
}
//178m
double MassFuncSheth1999(clucu_cosmology *cosmo,double sigma,int *status)
{
    double A = 0.3222;
    double a = 0.707;
    double p = 0.3;
    double delta_crit =1.68647;
    double nu = delta_crit / sigma;
    double anu2 = a * nu*nu;
	return A  *  sqrt(2.0*a/M_PI)  *  (1.0+pow(anu2,-p)) * nu * exp(-anu2/2.0);
}
//拟合函数，Jenkins 2001。M180m.这篇文章的主要结果是fof0.2，但我没收录
double MassFuncJenkins2001(clucu_cosmology *cosmo,double sigma,int *status)
{
    if(cosmo->config.halo_define_method != clucu_delta_180m)
        printf("wrong type in MassFuncJenkins2001");
    return 0.301 * exp(-pow(fabs(log(1.0 / sigma) + 0.64), 3.82));
}
//拟合函数，Tinker08。z=0-2.5
double MassFuncTinker2008(clucu_cosmology *cosmo,double sigma, double z,int *status)
{
//    if(cosmo->config.halo_define_method != clucu_delta_200m)
//        printf("wrong type in MassFuncTinker2008");
	double delta = 200.0; //SOmean，DELTA取200的意思
    double alpha = pow(10.0, -pow(0.75 / log10(delta / 75.0), 1.2)  );
	double A = 0.186 * pow(1.0 + z, -0.14);
	double a = 1.47 * pow(1.0 + z, -0.06);
	double b = 2.57 * pow(1.0 + z, -alpha);
	double c = 1.19;
	return A * (pow(sigma / b, -a) + 1.0) * exp(-c / pow(sigma, 2.0));
}
//fof 0.2
double MassFuncBhattacharya2011(clucu_cosmology *cosmo,double sigma, double z,int *status)
{
    double A = 0.333 * pow(1.+z,-0.11);
    double a = 0.788 * pow(1.+z,-0.01);
    double p = 0.807;
    double q = 1.795;
    double delta_crit =1.68647;
    double nu = delta_crit / sigma;
    double anu2 = a * nu*nu;
    double term1 = A * sqrt(2./M_PI);
    double term2 = exp( -0.5*anu2 );
    double term3 = 1.+pow(anu2,-p);
    double term4 = pow(   nu * sqrt(a),   q);//这里之前写错啦，大错误！
    return term1 * term2 *term3 * term4;
}
//200m,200c,500c,默认包括重子效应(table2的hydro)
double MassFuncBocquet2016(clucu_cosmology *cosmo,double sigma,double mass, double z,int *status)
{
    double A0,a0,b0,c0,Az,az,bz,cz;
    switch(cosmo->config.halo_define_method)
    {
        case clucu_delta_200m: A0=0.228; a0=2.15; b0=1.69; c0=1.30; Az=0.285; az=-0.058; bz=-0.366; cz=-0.045;break;
        case clucu_delta_200c: A0=0.202; a0=2.21; b0=2.00; c0=1.57; Az=1.147; az=0.375; bz=-1.074; cz=-0.196;break;
        case clucu_delta_500c: A0=0.180; a0=2.29; b0=2.44; c0=1.97; Az=1.088; az=0.150; bz=-1.008; cz=-0.322;break;
        default: printf("wrong type in MassFuncBocquet2016");
    }
    double A = A0 * pow(1.+z,Az);
    double a = a0 * pow(1.+z,az);
    double b = b0 * pow(1.+z,bz);
    double c = c0 * pow(1.+z,cz);
    double fsigma = A * ( pow(sigma/b,-a) +1. ) * exp( -c /sigma /sigma);//200m
    if(cosmo->config.halo_define_method==clucu_delta_200c)
    {
        double Omega_m = clucu_omega_x(cosmo,z,clucu_species_m_label,status);
        double lnM = log(mass);//Msun,不带h
        double gamma0 = 3.54E-2 + pow(Omega_m,0.09);
        double gamma1 = 4.56E-2 + 2.68E-2 / Omega_m;
        double gamma2 = 0.721 + 3.50E-2 / Omega_m;
        double gamma3 = 0.628 + 0.164 / Omega_m;
        double delta0 = -1.67E-2 + 2.18E-2 * Omega_m;
        double delta1 = 6.52E-3 - 6.86E-3 * Omega_m;
        double gamma = gamma0 + gamma1 * exp(-pow((gamma2 - z) / gamma3,2));
        double delta = delta0 + delta1 * z;
        double M200c_M200m = gamma + delta * lnM;
        double fsigma = fsigma * M200c_M200m;
    }
    else if(cosmo->config.halo_define_method==clucu_delta_500c)
    {
        double Omega_m = clucu_omega_x(cosmo,z,clucu_species_m_label,status);
        double lnM = log(mass);//Msun,不带h
        double alpha0 = 0.880 + 0.329 * Omega_m;
        double alpha1 = 1.00 + 4.31E-2 / Omega_m;
        double alpha2 = -0.365 + 0.254 / Omega_m;
        double alpha = alpha0 * (alpha1 * z + alpha2) / (z + alpha2);
        double beta = -1.7E-2 + 3.74E-3 * Omega_m;
        double M500c_M200m = alpha + beta * lnM;
        double fsigma = fsigma * M500c_M200m;
    }
    return fsigma;
}

/*----------------------*/
/*        非高斯性       */
/*----------------------*/
//计算sigma * S3 (M),对log10M插值. 质量单位不带h
double SigmaS3(clucu_cosmology *cosmo,double mass)
{
    double Om = cosmo->param.Omega_m;
    double Ob = cosmo->param.Omega_b;
    double h = cosmo->param.h;
    double fnl = cosmo->param.f_nl;
    double sigma8 = cosmo->param.sigma8;

    double Gamma = Om*h * exp( -Ob*( 1.+sqrt(2.*h)/Om ) );
    double m10 = mass*h/pow(10.,10)  *pow(Gamma,3) / cosmo->param.Omh2;
    double index = -0.0272 -0.11*(cosmo->param.n_s-0.96) -0.0008*log10(m10);
    
    return 8.66e-5*fnl*Om*pow(Gamma,-1.4)*sigma8 *pow(m10,index);
}
void clucu_compute_sigmaS3(clucu_cosmology *cosmo,int *status)
{
    int nm = cosmo->spline_param.LOG10M_SPLINE_NM;
    double *m = NULL;
    double *y = Create1Grid(nm);

    // create linearly-spaced values of log10-mass.
    m = clucu_linear_spacing(cosmo->spline_param.LOG10M_SPLINE_MIN,
                            cosmo->spline_param.LOG10M_SPLINE_MAX, nm);
    if(m == NULL ||
    (fabs(m[0]-cosmo->spline_param.LOG10M_SPLINE_MIN)>1e-5) ||
    (fabs(m[nm-1]-cosmo->spline_param.LOG10M_SPLINE_MAX)>1e-5) ||
    (m[nm-1]>10E17)) 
    {
        *status = CLUCU_ERROR_MEMORY;
        CLUCU_RAISE_WARNING(*status,"Error creating linear spacing in m");
        clucu_cosmology_set_status_message(cosmo,
                                        "clucu_massfunc.c: clucu_compute_logsigma(): "
                                        "Error creating linear spacing in m\n");
    }
    //set y
    for(int i=0;i<nm;i++)
        y[i] = SigmaS3(cosmo,pow(10.,m[i]));
    
    //alloc
    gsl_spline *spline_sigmaS3 = NULL;
    if (*status == 0)
    {
        spline_sigmaS3 = gsl_spline_alloc(gsl_interp_cspline, nm);
        if (spline_sigmaS3 == NULL)
        {
            *status = CLUCU_ERROR_MEMORY;
            CLUCU_RAISE_WARNING(*status,"ran out of memory");
        }
    }

    //init
    if (*status == 0)
    {
        int gslstatus = gsl_spline_init(spline_sigmaS3, m, y, nm);

        if (gslstatus)
        {
            *status = CLUCU_ERROR_SPLINE;
            CLUCU_RAISE_WARNING(*status,"Error creating sigmaS3 spline");
        }
    }

    if (*status == 0)
    {
        cosmo->data.spline_sigmaS3 = spline_sigmaS3;
        cosmo->computed_sigmaS3 = true;
    }
    else
    {
        gsl_spline_free(spline_sigmaS3);
    }

    free(m);
    free(y);
}
double NGhmf_factor(clucu_cosmology *cosmo,double sigma, double mass,double z,int *status)
{
    if(fabs(cosmo->param.f_nl)<1.0e-15)
        return 1.;
    if(!cosmo->computed_sigmaS3){
        *status = CLUCU_ERROR_GROWTH_INIT;
        CLUCU_RAISE_WARNING(*status,"sigmaS3 splines have not been precomputed!");
        return NAN;
    }
    double dlnsigma_dlog10M=clucu_dlnsigma_dlnm(cosmo, mass, z, status);
    double sigmaS3;
    int gslstatus = gsl_spline_eval_e(cosmo->data.spline_sigmaS3,log10(mass),NULL, &sigmaS3);
    double dsigmaS3_dlog10M;
    gslstatus = gsl_spline_eval_deriv_e(cosmo->data.spline_sigmaS3,log10(mass),NULL, &dsigmaS3_dlog10M);
    if(gslstatus != GSL_SUCCESS){
        CLUCU_RAISE_GSL_WARNING(gslstatus,"m=%.2e is outside interpolation range [%.2e,%.2e]",mass,pow(10.,cosmo->data.spline_sigmaS3->interp->xmin),pow(10.,cosmo->data.spline_sigmaS3->interp->xmax));
        *status |= gslstatus;
        return NAN;
    }
    double nu = 1.686/sigma;
    double dsigmaS3_dlnnu = -dsigmaS3_dlog10M/dlnsigma_dlog10M;
    return 1. + (pow(nu,3.)-3.*nu)*sigmaS3/6. -(nu-1./nu)*dsigmaS3_dlnnu/6.;
}

double clucu_fitting_function(clucu_cosmology *cosmo,double sigma, double mass,double z,int *status)
{
    double result;
	switch(cosmo->config.mass_function_method)
	{
		case clucu_massfunc_press1974:          result = MassFuncPress1974(sigma,status);break;
        case clucu_massfunc_sheth1999:          result = MassFuncSheth1999(cosmo,sigma,status);break;
        case clucu_massfunc_jenkins2001:        result = MassFuncJenkins2001(cosmo,sigma,status);break;
		case clucu_massfunc_tinker2008:         result = MassFuncTinker2008(cosmo,sigma, z,status);break;
        case clucu_massfunc_bhattacharya2011:   result = MassFuncBhattacharya2011(cosmo,sigma,z,status);break;
        case clucu_massfunc_bocquet2016:        result = MassFuncBocquet2016(cosmo,sigma,mass,z,status);break;
        default: printf("wrong type in clucu_concentration\n");
	}
    return result * NGhmf_factor(cosmo,sigma,mass, z,status);
}

//dn/dlnm，输入一个m,z,输出dn/dlnm(m,z)。单位 1/Mpc^3
double clucu_halo_dn_dlnm(clucu_cosmology *cosmo,double mass,double z,int *status)
{
    double rho_m = clucu_rho_x(cosmo, 0., clucu_species_m_label, 1, status);
    double sigma = clucu_sigma(cosmo,mass,z,status);
    double abs_dlnsigma_dlnm = clucu_dlnsigma_dlnm(cosmo,mass,z,status);
    double f_sigma = clucu_fitting_function(cosmo,sigma,mass, z,status);
    double dn_dlnm = rho_m*f_sigma*abs_dlnsigma_dlnm/mass;
    return dn_dlnm;
}
double clucu_halo_bias(clucu_cosmology *cosmo,double mass,double z,int *status)
{
    double sigma = clucu_sigma(cosmo,mass,z,status);
    double bias = clucu_halo_bias_sigma(cosmo,sigma);
    return bias;
}
// double clucu_halo_bias_k(clucu_cosmology *cosmo,double mass,double z,double k,int *status)
// {
//     double rho_m = clucu_rho_x(cosmo, 0., clucu_species_m_label, 1,status);
//     double sigma = clucu_sigma(cosmo,mass,z,status);
//     double bias = clucu_halo_bias_sigma(cosmo,sigma);
//     double SoverL = 
//     return (bias-1.)*SoverL;
//     double GISDB = GISDB_factor(cosmo,k,z);
//     return ((bias-1.)*GISDB+1.);
// }

double clucu_halo_dndlnm_bias(clucu_cosmology *cosmo,double mass,double z,int *status)
{
    double rho_m = clucu_rho_x(cosmo, 0., clucu_species_m_label, 1,status);
    double sigma = clucu_sigma(cosmo,mass,z,status);
    double abs_dlnsigma_dlnm = clucu_dlnsigma_dlnm(cosmo,mass,z,status);
    double f_sigma = clucu_fitting_function(cosmo,sigma,mass, z,status);
    double dn_dlnm = rho_m*f_sigma*abs_dlnsigma_dlnm/mass;
    double bias = clucu_halo_bias_sigma(cosmo,sigma);
    return dn_dlnm*bias;
}
double clucu_halo_dndlnm_bias_k(clucu_cosmology *cosmo,double mass,double z,double k,int *status)
{
    double rho_m = clucu_rho_x(cosmo, 0., clucu_species_m_label, 1,status);
    double sigma = clucu_sigma(cosmo,mass,z,status);
    double abs_dlnsigma_dlnm = clucu_dlnsigma_dlnm(cosmo,mass,z,status);
    double f_sigma = clucu_fitting_function(cosmo,sigma,mass, z,status);
    double dn_dlnm = rho_m*f_sigma*abs_dlnsigma_dlnm/mass;
    double bias = clucu_halo_bias_sigma(cosmo,sigma);
    double GISDB = GISDB_factor(cosmo,k,z);
    return dn_dlnm*((bias-1.)*GISDB+1.);
}