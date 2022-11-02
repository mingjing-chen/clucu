//根据nfw profile
//根据c_old计算c_new
#include "clucu_halo.h"
double dc_NakamuraSuto(clucu_cosmology *cosmo, double z,int *status)
{
    double Om_mz = clucu_omega_x(cosmo, z, clucu_species_m_label, status);
    double dc0 = (3./20.)*pow(12.*M_PI,2./3.);
    double dc = dc0*(1.+0.012299*log10(Om_mz));
    return dc;
}
double Dv_NakamuraSuto(clucu_cosmology *cosmo, double z,int *status)
{
    double Om_mz = clucu_omega_x(cosmo, z, clucu_species_m_label, status);

    return 18.*M_PI*M_PI *pow(Om_mz,-0.573);
}
/*
 * convert_concentration_single finds the concentration c'
 * for a mass definition with overdensity Delta' given
 * the concentration c for a mass definition with overdensity
 * Delta assuming an NFW density profile.
 *
 * To do so, it solves the following equation numerically:
 *    c^3 * f(c) = Delta' f(c')
 * where f(x) = x^3 * (x + 1) / ((1+x) * log(x + 1) - x).
 * 
 * The equation is solved using a Newton-Raphson approach.
 * The functions nfw_fx, nfw_f, nfw_df and nfw_fdf implement
 * the function whose zero we try to find and its derivative.
 *
 * clucu_convert_concentration does the same thing for an
 * array of input concentrations.
 */
/*--------------------------*/
/*        nfw profile       */
/*--------------------------*/
double nfw_mc(double c)
{
    return log(1.0 + c) - c/(1.+c);
}

/*---------------------------*/
/*        mass convert       */
/*---------------------------*/
//升级算法，使用ccl的导数寻根算法
//并且200c,180m的信息包含在d_factor中
//这里的fx，是我之前的程序中的x^3/fx，即星系团的密度rho分之一
static double nfw_fx(double x)
{
  if(x>0.01) {
    double xp1=1+x;
    double lxp1=log(xp1);
    double den=1./(xp1*lxp1-x);
    return xp1*x*x*x*den;
  }
  else {
    return 2*x;
  }
}

static double nfw_f(double x,void *params)
{
  double offset = *((double *)params);
  return nfw_fx(x)-offset;
}

static double nfw_df(double x,void *params)
{
  if(x>0.01) {
    double xp1=1+x;
    double lxp1=log(xp1);
    double den=1./(xp1*lxp1-x);
    return x*x*(3*xp1*xp1*lxp1-x*(4*x+3))*den*den;
  }
  else {
    return 2;
  }
}

static void nfw_fdf(double x,void *params,
		    double *y, double *dy)
{
  double offset = *((double *)params);
  if(x>0.01) {
    double xp1=1+x;
    double lxp1=log(xp1);
    double den=1./(xp1*lxp1-x);
    *y=xp1*x*x*x*den-offset;
    *dy=x*x*(3*xp1*xp1*lxp1-x*(4*x+3))*den*den;
  }
  else {
    *y=2*x-offset;
    *dy=2;
  }
}

static int convert_concentration_single(double d_factor, double c_old,
					double *c_new, double c_start)
{
  double c0, offset = d_factor * nfw_fx(c_old);
  int status, iter=0, max_iter=100;
  gsl_function_fdf FDF;
  gsl_root_fdfsolver *s=NULL;//创建求解器
  FDF.f = &nfw_f;
  FDF.df = &nfw_df;
  FDF.fdf = &nfw_fdf;
  FDF.params = &offset;
  
  s=gsl_root_fdfsolver_alloc(gsl_root_fdfsolver_newton);//求解方法
  if (s==NULL)
    return CLUCU_ERROR_MEMORY;
  
  gsl_root_fdfsolver_set (s, &FDF, c_start);//初始化求解器
  *c_new = c_start;
  do
    {
      iter++;
      c0 = *c_new;
      status = gsl_root_fdfsolver_iterate (s);
      *c_new = gsl_root_fdfsolver_root (s);
      status = gsl_root_test_delta (*c_new, c0, 0, 1e-4);

    }
  while (status == GSL_CONTINUE && iter < max_iter);
  gsl_root_fdfsolver_free (s);

  return status;
}
/**
 * @brief 
 * 例：M180cz -> M200mz: 
 * delta_old = 180*rho_cz = 180 * E(z)^2 * rho_c0
 * delta_new = 200*rho_mz = 200 * (1+z)^3 * Omega_m0 * rho_c0
 * 
 * @param cosmo 
 * @param delta_old 
 * @param nc 
 * @param c_old 
 * @param delta_new 
 * @param c_new 
 * @param status 
 */

void clucu_convert_concentration(clucu_cosmology *cosmo,
			       double delta_old, int nc, double c_old[],
			       double delta_new, double c_new[],int *status)
{
  if(nc<=0)
    return;

  int ii,st=0;
  double d_factor = delta_old/delta_new;
  for(ii=0;ii<nc;ii++) {
    st=convert_concentration_single(d_factor, c_old[ii], &(c_new[ii]), c_old[ii]);
    if(st!=GSL_SUCCESS) {
      *status=CLUCU_ERROR_ROOT;
      clucu_cosmology_set_status_message(cosmo,
        "clucu_mass_conversion.c: clucu_convert_concentration(): "
        "NR solver failed to find a root\n");
      return;
    }
  }
}
//单位: M_sun/Mpc^3
double clucu_rho_delta(clucu_cosmology *cosmo,halo_define_type type,double z,int *status)
{
    switch(type)
    {
        case clucu_delta_200m: return 200. * clucu_rho_x(cosmo,z,clucu_species_m_label,0,status);break;
        case clucu_delta_200c: return 200. * clucu_rho_x(cosmo,z,clucu_species_crit_label,0,status);break;
        case clucu_delta_500m: return 500. * clucu_rho_x(cosmo,z,clucu_species_m_label,0,status);break;
        case clucu_delta_500c: return 500. * clucu_rho_x(cosmo,z,clucu_species_crit_label,0,status);break;
        case clucu_delta_180m: return 180. * clucu_rho_x(cosmo,z,clucu_species_m_label,0,status);break;
        case clucu_delta_vir: return Dv_NakamuraSuto(cosmo,z,status) *  clucu_rho_x(cosmo,z,clucu_species_m_label,0,status);break;
    }
    printf("wrong type in clucu_rho_delta");
    return 1;
}
// 例：M180cz -> M200mz: 
// delta_old = 180*rho_cz = 180 * E(z)^2 * rho_c0
// delta_new = 200*rho_mz = 200 * (1+z)^3 * Omega_m0 * rho_c0
double clucu_Mold_to_Mnew(clucu_cosmology *cosmo,halo_define_type Told,halo_define_type Tnew,double Mold,double z,int *status)
{
    double delta_old = clucu_rho_delta(cosmo,Told,z,status);
    double delta_new = clucu_rho_delta(cosmo,Tnew,z,status);
    double d_factor = delta_old/delta_new;
    double c_old=clucu_concentration(cosmo,Told,Mold,z,status);
    double c_new;
    convert_concentration_single( d_factor,  c_old,  &c_new, c_old);
    return nfw_mc(c_new) / nfw_mc(c_old) * Mold;
}
double clucu_M200m_to_M500c(clucu_cosmology *cosmo,double M200m,double z,int *status)
{
    return clucu_Mold_to_Mnew(cosmo, clucu_delta_200m, clucu_delta_500c, M200m, z,status);
}
double clucu_M200c_to_M180m(clucu_cosmology *cosmo,double M200c,double z,int *status)
{
    return clucu_Mold_to_Mnew(cosmo, clucu_delta_200c, clucu_delta_180m, M200c, z,status);
}
double clucu_M200m_to_M200c(clucu_cosmology *cosmo,double M200m,double z,int *status)
{
    return clucu_Mold_to_Mnew(cosmo, clucu_delta_200m, clucu_delta_200c, M200m, z,status);
}
double clucu_M200c_to_M200m(clucu_cosmology *cosmo,double M200c,double z,int *status)
{
    return clucu_Mold_to_Mnew(cosmo, clucu_delta_200c, clucu_delta_200m, M200c, z,status);
}
//如果要算M500c_to_M200m怎么办呢，没有c500c的定义
//但这些concentration都有200c，我们用200c算到目标type，然后二维插值（log10 mass,a），和sigma的插值一样
//Tnew must be 200c or 200m
void clucu_compute_massconvetr(clucu_cosmology *cosmo,halo_define_type Told,halo_define_type T200,int *status)
{
    int nz=100;
    double *z=clucu_linear_spacing(0.,2.,nz);
    int nm = 100;
    double *log10m200 = clucu_linear_spacing(13.,18,nm);
    double *log10m_new = NULL;
    double *y = NULL;

    log10m_new = malloc(sizeof(double)*nm);
    y = malloc(sizeof(double)*nm*nz);
    int number_z,number_m;
    for(int number = 0;number<nz*nm;number+=1)
    {
        number_z = number/nm;
        number_m = number % nm;
        //由M200算目标质量
        log10m_new[number_m]=log10(clucu_Mold_to_Mnew(cosmo, T200, Told,pow( 10.,log10m200[number_m]), z[number_z],status));
        //y得设置成m200
        y[number_z*nm + number_m] = log10m200[number_m];
    }
    gsl_spline2d *Mold_to_M200 = NULL;
    Mold_to_M200 = gsl_spline2d_alloc(gsl_interp2d_bicubic, nm, nz);
    gsl_spline2d_init(Mold_to_M200, log10m_new, z, y, nm, nz);

    cosmo->computed_massconvert = true;
    cosmo->data.spline_massconvert = Mold_to_M200;

    free(z);
    free(log10m200);
    free(log10m_new);
    free(y);

}
double clucu_Mold_to_M200(clucu_cosmology *cosmo,halo_define_type Told,halo_define_type T200, double Mold, double z,int *status)
{
    
    if(log10(Mold) < cosmo->data.spline_massconvert->interp_object.xmin
    || log10(Mold) > cosmo->data.spline_massconvert->interp_object.xmax
    || z < cosmo->data.spline_massconvert->interp_object.ymin
    || z > cosmo->data.spline_massconvert->interp_object.ymax){
        *status = CLUCU_ERROR_SPLINE_EV;
        CLUCU_RAISE_WARNING(*status,"mass=%.2e,z=%.2e is outside interpolation range [%.2e,%.2e][%.2e,%.2e]. [cen-min,max-cen]: ln_mass:[%.2e,%.2e] z:[%.2e,%.2e]",
                                        Mold,z,
                                        pow(10.,cosmo->data.spline_massconvert->interp_object.xmin),
                                        pow(10.,cosmo->data.spline_massconvert->interp_object.xmax),
                                        cosmo->data.spline_massconvert->interp_object.ymin,
                                        cosmo->data.spline_massconvert->interp_object.ymax,

                                        log10(Mold)-cosmo->data.spline_massconvert->interp_object.xmin,
                                        cosmo->data.spline_massconvert->interp_object.xmax-log10(Mold),
                                        z-cosmo->data.spline_massconvert->interp_object.ymin,
                                        cosmo->data.spline_massconvert->interp_object.ymax-z);
        }
    double log10M200;
    int gslstatus = gsl_spline2d_eval_e(cosmo->data.spline_massconvert, log10(Mold),
                                        z, NULL, NULL, &log10M200);
    if(gslstatus != GSL_SUCCESS)
    {
        CLUCU_RAISE_GSL_WARNING(gslstatus,NULL);
        *status|= gslstatus;
    }
    return pow(10.,log10M200);
}

/*------------------------*/
/*        nfw for k       */
/*------------------------*/
double nfw_scale(clucu_cosmology *cosmo,double mass,double z,double *out,int *status)
{
    double m=mass;
    //这里mass的质量定义，得和concentration的质量定义一样
    double c=clucu_concentration(cosmo,cosmo->config.halo_define_method, mass, z,status);
    double mc=nfw_mc(c);
    double rho =  clucu_rho_delta(cosmo,cosmo->config.halo_define_method,z,status);//单位: M_sun/Mpc^3
    out[0] = c*c*c / (3.*mc) * rho;//单位: M_sun/Mpc^3
    out[1] = pow( m / (4.* M_PI * (out[0]) * mc) ,1./3.);//单位:Mpc
}
//nfw profile uk
//单位都不带h
double nfw_uk(clucu_cosmology *cosmo,double k,double mass,double z,int *status)
{
    double c=clucu_concentration(cosmo,cosmo->config.halo_define_method, mass, z,status);
    double out[2];
    nfw_scale(cosmo,mass,z, out,status);
    double rhos=out[0];
    double rs=out[1];
    double x=(1.+z)*k*rs;
    //double A=4.*M_PI*rhos*pow(rs,3.)/mass;
    double A = 1./nfw_mc(c);
    double SSSi=Si((1.+c)*x)-Si(x);
    double CCCi=Ci((1.+c)*x)-Ci(x);
    return A*(   sin(x) * SSSi + cos(x) * CCCi - sin(c*x)/((1.+c)*x)  );
}