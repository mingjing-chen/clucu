#include "clucu_survey.h"

static double expective_ln_mass_ob1(clucu_cosmology *cosmo, double ln_mass_true,double z_true)
{
    double M_pivot = 3.0e14/cosmo->param.h;  //M_sun/h改为Msun
	double z_pivot = 0.6;
	double M = exp(ln_mass_true);
	double z = z_true;
	double mean_lnN = cosmo->survey.mass_ob_A + cosmo->survey.mass_ob_B*log(M/M_pivot) + cosmo->survey.mass_ob_Bz*log((1.0+z)/(1.0+z_pivot)) + cosmo->survey.mass_ob_Cz*pow( log((1.0+z)/(1.0+z_pivot)), 2.) ;
    return mean_lnN;
}
static double expective_ln_mass_ob2(clucu_cosmology *cosmo, double ln_mass_true,double z_true)
{
    double ln_mass_bias = cosmo->survey.lnM_b0 + cosmo->survey.s_b[0] * log(1.+z_true);
    return ln_mass_true + ln_mass_bias;
}
static double expective_ln_mass_ob3(clucu_cosmology *cosmo, double ln_mass_true,double z_true)
{
    double ln_mass_pivot = log(3.0e14/cosmo->param.h);  //M_sun/h改为Msun
    double ln_mass_bias = cosmo->survey.lnM_b0;
    for(int i=0;i<3;i++){
        ln_mass_bias += cosmo->survey.q_b[i] * pow(ln_mass_true-ln_mass_pivot,i+1);
        ln_mass_bias += cosmo->survey.s_b[i] * pow(z_true,i+1);
    }

    return ln_mass_true + ln_mass_bias;
}
static double sigma_ln_mass_ob1(clucu_cosmology *cosmo,double ln_mass_true,double z_true)
{
    double M_pivot = 3.0e14/cosmo->param.h;  //M_sun/h改为Msun
	double z_pivot = 0.6;
	double M = exp(ln_mass_true);
	double z = z_true;
    double sigma = cosmo->survey.mass_ob_sigma0 + cosmo->survey.mass_ob_q*log(M/M_pivot) + cosmo->survey.mass_ob_qz*log((1.0+z)/(1.0+z_pivot)) + cosmo->survey.mass_ob_pz*pow( log((1.0+z)/(1.0+z_pivot)), 2.) ;
    return sigma;
}
double sigma_ln_mass_ob2(clucu_cosmology *cosmo,double ln_mass_true,double z_true)
{
    double M200m=exp(ln_mass_true);
    M200m=max(M200m,cosmo->survey.mass_ob_min);
    double M500c=M200mz_to_M500cz(cosmo,M200m,z_true);
    double M=M500c;
    double Mp=3.0e14/cosmo->param.h;  //M_sun/h改为Msun
    double zp=0.45;
    double Al=cosmo->survey.mass_ob_A;//79.8
    double Bl=cosmo->survey.mass_ob_B;//0.93
    double Cl=cosmo->survey.mass_ob_Bz;//-0.49
    double Dl=cosmo->survey.mass_ob_sigma0;//0.217
    double kappa=cosmo->survey.mass_ob_q;//0.01
    double beta=cosmo->survey.mass_ob_qz;//1
    double lambda_ob = Al * pow(M/Mp,Bl) * pow((1.+z_true)/(1.+zp),Cl);
    double sigma_lnl = sqrt(Dl*Dl + 1./lambda_ob);
    double sigma_lnM0 = sigma_lnl/Bl;
    double sigma_lnM = sqrt( sigma_lnM0*sigma_lnM0 + kappa*pow(1.+z_true, 2.*beta) );

    return sigma_lnM;
}
static double sigma_ln_mass_ob3(clucu_cosmology *cosmo,double ln_mass_true,double z_true)
{
    double ln_mass_pivot = log(3.0e14/cosmo->param.h);  //M_sun/h改为Msun
    double sigma_lnM = cosmo->survey.sigma_lnM0;
    for(int i=0;i<3;i++){
        sigma_lnM += cosmo->survey.q_sigma_lnm[i] * pow(ln_mass_true-ln_mass_pivot,i+1);
        sigma_lnM += cosmo->survey.s_sigma_lnm[i] * pow(z_true,i+1);
    }
    // if(sigma_lnM<0.)
    //     CLUCU_RAISE_WARNING(23421323,"q=[%.1e,%.1e,%.1e],s=[%.1e,%.1e,%.1e],mass_true=%.1e,z_true=%.1e,sigma_lnM=%.2f",
    //                         cosmo->survey.q_sigma_lnm[0],
    //                         cosmo->survey.q_sigma_lnm[1],
    //                         cosmo->survey.q_sigma_lnm[2],
    //                         cosmo->survey.s_sigma_lnm[0],
    //                         cosmo->survey.s_sigma_lnm[1],
    //                         cosmo->survey.s_sigma_lnm[2],
    //                         exp(ln_mass_true),
    //                         z_true,
    //                         sigma_lnM
    //                         )
    return sigma_lnM;
}
//公式里M_pivot带h，所以这里的质量 M_sun/h
double mass_ob_probability(clucu_cosmology *cosmo,double ln_mass_ob, double ln_mass_true,double z_true)
{
	double mean_lnN,sigma;
    switch(cosmo->config.mass_observable_method)
    {
        case clucu_murrata2019: mean_lnN = expective_ln_mass_ob1(cosmo,ln_mass_true,z_true);sigma = sigma_ln_mass_ob1(cosmo,ln_mass_true,z_true);break;
        case clucu_my_csst: mean_lnN = expective_ln_mass_ob2(cosmo,ln_mass_true,z_true);sigma = sigma_ln_mass_ob2(cosmo,ln_mass_true,z_true);break;
        case clucu_oguri2011: mean_lnN = expective_ln_mass_ob3(cosmo,ln_mass_true,z_true);sigma = sigma_ln_mass_ob3(cosmo,ln_mass_true,z_true);break;
    }
    double x = ln_mass_ob - mean_lnN;
	double probability;
	//if sigma=0,p=1(x=0),p=0(x!=0)
	if(fabs(sigma)<=1.0e-15)
		if(fabs(x)<=1.0e-15)
			probability = 1.0;
		else 
			probability = 0.0;
	else 
		probability =  exp(-0.5*pow(x/sigma,2)) /  ( sqrt(2.0*M_PI)*sigma );

	return probability;
}
double redshift_probability(clucu_cosmology *cosmo,double z_ob,double z_true)
{
    double sigma0 = cosmo->survey.redshift_sigma0;
	double sigma = sigma0 *(1.0+z_true);
	z_ob = z_ob;
	z_true = z_true;
	double x = z_ob - z_true;
	return GaussianDistribution(x,sigma);
}
double zob_to_ztrue(clucu_cosmology *cosmo,double z_ob,double *out)
{
	out[0] = z_ob-3.5*cosmo->survey.redshift_sigma0*(1.0+z_ob);
    out[1] = z_ob-3.5*cosmo->survey.redshift_sigma0*(1.0+z_ob);
    out[2] = z_ob-2.5*cosmo->survey.redshift_sigma0*(1.0+z_ob);
    out[3] = z_ob+2.5*cosmo->survey.redshift_sigma0*(1.0+z_ob);
    out[4] = z_ob+3.5*cosmo->survey.redshift_sigma0*(1.0+z_ob);
    out[5] = z_ob+3.5*cosmo->survey.redshift_sigma0*(1.0+z_ob);
    for(int i=0;i<=5;i++)
    {
        if (out[i]<=0.)
            out[i]=0.;
    }

	return 1;
}

double galaxy_redshift_distribution(double z_source,void *params)
{
    clucu_cosmology *cosmo=(clucu_cosmology *)params;
    double z0=cosmo->survey.source_redshift/3.;
    double x=z_source/z0;
    return 0.5 * x*x *exp(-x) / z0;
}
double average_one_over_chis_integration(double z_source,void *params)
{
    clucu_cosmology *cosmo=(clucu_cosmology *)params;
    int status=0.;
    double chi=clucu_comoving_angular_distance(cosmo,z_source,&status);
    return galaxy_redshift_distribution(z_source,cosmo)/chi;//单位 1/Mpc
}
void clucu_compute_average_one_over_chis(clucu_cosmology *cosmo)
{
    double result1,result2;
    double a=1.5,b=5.;
    gsl_integration_cquad_workspace *workspace =  NULL;
    workspace = gsl_integration_cquad_workspace_alloc(cosmo->gsl_param.N_ITERATION);
    gsl_function F;
    F.function=&average_one_over_chis_integration;
    F.params=cosmo;
    gsl_integration_cquad(&F,
                            a,
                            b,
                            0.0, 1.0e-5,
                            workspace,&result1,NULL,NULL);
    F.function=&galaxy_redshift_distribution;
    gsl_integration_cquad(&F,
                            a,
                            b,
                            0.0, 1.0e-5,
                            workspace,&result2,NULL,NULL);
    gsl_integration_cquad_workspace_free(workspace);
    cosmo->data.average_one_over_chis=result1/result2;
    cosmo->computed_average_one_over_chis =true;
}

/*----------------------*/
/*        质量下限       */
/*----------------------*/

static double NFW_fx_x3(double x)
{
	return nfw_mc(x) / pow(x, 3.0);
}
struct NFW_params
{
	double z;
	clucu_cosmology *cosmo;
}NFW_params;
//M200cz,to,M180mz
static double Mf3(double x, void *params)
{
    
	//将指向空类型的指针params，转化成指向结构体Mf3_params的指针  (struct Mf3_params *)params
	double z = ((struct NFW_params *)params)->z;
	clucu_cosmology *cosmo = ((struct NFW_params *)params)->cosmo;
    //TODO需要改为concentrate
	double x200c = 5.0;
    int status=0;
    double Omz = clucu_omega_x(cosmo,z,clucu_species_m_label,&status);
	return NFW_fx_x3(x) / NFW_fx_x3(x200c) - (180.0 / 200.0) * Omz; //rho_m(z), rho_crit(z)
}
//M200c,to,M180m
double M200cz_to_M180mz(clucu_cosmology *cosmo,double M200c, double z)
{

	int status;
	int iter = 0, max_iter = 40; //迭代计数
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 0.0; //root的初始值
	double x_lo = 0.1, x_hi = 40.0;
	gsl_function F;
	struct NFW_params params = {z, cosmo};
	F.function = &Mf3;
	F.params = &params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_lo, x_hi,0, 1.0e-10);

	} while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free(s);

	double M180m = nfw_mc(r) / nfw_mc(5.0) * M200c;

	return M180m;
}
static double Mf4(double x, void *params)
{
	//将指向空类型的指针params，转化成指向结构体Mf3_params的指针  (struct Mf3_params *)params
	double z = ((struct NFW_params *)params)->z;
	clucu_cosmology *cosmo = ((struct NFW_params *)params)->cosmo;
	double x180 = 5.0*180.0/200.0;
    int status=0;
    double Omz = clucu_omega_x(cosmo,z,clucu_species_m_label,&status);
	return NFW_fx_x3(x) / NFW_fx_x3(x180) - 200.0 / 180.0 / Omz; //rho_m(z), rho_crit(z)
}
//M180mz,to,M200cz
double M180mz_to_M200cz(clucu_cosmology *cosmo,double M180m, double z)
{

	int status;
	int iter = 0, max_iter = 100; //迭代计数
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 0.0; //root的初始值
	double x_lo = 0.1, x_hi = 40.0;
	gsl_function F;
	struct NFW_params params = {z, cosmo};
	F.function = &Mf4;
	F.params = &params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_lo, x_hi,
										0, 0.001);

	} while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free(s);

	double M200c = nfw_mc(r) / nfw_mc(5.0*180.0/200.0) * M180m;

	return M200c;
}
/**
 * @brief M200m,to,M500c
 * wang sheng 的文章都是带z,M200cz,to,M180mz
 * tinker2008明确表明是200mz
 * 有z=0的定义吗？
 * @param cosmo 
 * @param M200m 
 * @param z 
 * @return double 
 */
static double Mf5(double x, void *params)
{
    //将指向空类型的指针params，转化成指向结构体Mf3_params的指针  (struct Mf3_params *)params
	double z = ((struct NFW_params *)params)->z;
	clucu_cosmology *cosmo = ((struct NFW_params *)params)->cosmo;
    //TODO需要改为concentrate
	//double x200m = 11.915226265364348;
    double x200m = 10.;
    //11.915226265364348
    int status=0;
    double Omz = clucu_omega_x(cosmo,z,clucu_species_m_label,&status);
	return NFW_fx_x3(x) / NFW_fx_x3(x200m) - (500. / 200.0) / Omz; //rho_m(z), rho_crit(z)
}
double M200mz_to_M500cz(clucu_cosmology *cosmo,double M200m, double z)
{
    double x200m = 11.915226265364348;
    int status;
	int iter = 0, max_iter = 40; //迭代计数
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 0.0; //root的初始值
	double x_lo = 0.1, x_hi = 40.0;
	gsl_function F;
	struct NFW_params params = {z, cosmo};
	F.function = &Mf5;
	F.params = &params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_lo, x_hi,0, 1.0e-10);

	} while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free(s);

	double M500 = nfw_mc(r) / nfw_mc(x200m) * M200m;

	return M500;
}
static double Mf6(double x, void *params)
{
    //将指向空类型的指针params，转化成指向结构体Mf3_params的指针  (struct Mf3_params *)params
	double z = ((struct NFW_params *)params)->z;
	clucu_cosmology *cosmo = ((struct NFW_params *)params)->cosmo;
    //TODO需要改为concentrate
	double x500c = 5.;
    int status=0;
    double Omz = clucu_omega_x(cosmo,z,clucu_species_m_label,&status);
	return NFW_fx_x3(x) / NFW_fx_x3(x500c) - (200. / 500.0) * Omz; //rho_m(z), rho_crit(z)
}
double M500cz_to_M200mz(clucu_cosmology *cosmo,double M500c, double z)
{
    //double x500c =4.1265466796234529;
    double x500c = 5.;
    int status;
	int iter = 0, max_iter = 40; //迭代计数
	const gsl_root_fsolver_type *T;
	gsl_root_fsolver *s;
	double r = 0.0; //root的初始值
	double x_lo = 0.1, x_hi = 40.0;
	gsl_function F;
	struct NFW_params params = {z, cosmo};
	F.function = &Mf6;
	F.params = &params;

	T = gsl_root_fsolver_brent;
	s = gsl_root_fsolver_alloc(T);
	gsl_root_fsolver_set(s, &F, x_lo, x_hi);

	do
	{
		iter++;
		status = gsl_root_fsolver_iterate(s);
		r = gsl_root_fsolver_root(s);
		x_lo = gsl_root_fsolver_x_lower(s);
		x_hi = gsl_root_fsolver_x_upper(s);
		status = gsl_root_test_interval(x_lo, x_hi,0, 1.0e-10);

	} while (status == GSL_CONTINUE && iter < max_iter);

	gsl_root_fsolver_free(s);

	double M200m = nfw_mc(r) / nfw_mc(x500c) * M500c;

	return M200m;
}
static double Mmin_SZE180(clucu_cosmology *cosmo,double z)
{
	double log10ASZ=cosmo->survey.log10ASZ;
	double betaSZ=cosmo->survey.betaSZ;
	double gammaSZ=cosmo->survey.gammaSZ;
	double fSZ=cosmo->survey.flux;
	double M0 = 1.0e15; //单位 M_sun
	double ASZ = pow(10.0, log10ASZ);
int status=0;
	double L = fSZ * pow(clucu_angular_diameter_distance(cosmo,z,&status), 2.0); //fSZ单位mJy，角直径距离单位Mpc
	double fnu = 0.9521;								 //the frequency dependence of the SZE distortion
	double fICM = 0.12;
	double Mbeta = L / (ASZ * fnu * fICM * pow(clucu_h_over_h0(cosmo,z,&status), 2.0 / 3.0) * pow(1 + z, gammaSZ));
	double M200 = pow(Mbeta, 1.0 / betaSZ) * M0; //输出质量，单位M_sun

	double M180 = M200cz_to_M180mz(cosmo,M200, z);
	if (M180 < 1.0e14 / cosmo->param.h)
		M180 = 1.0e14 / cosmo->param.h;
	return M180; //输出质量，单位M_sun
}
static double Mmin_xray180(clucu_cosmology *cosmo,double z)
{
    double log10Ax=cosmo->survey.log10ASZ;
	double betax=cosmo->survey.betaSZ;
	double gammax=cosmo->survey.gammaSZ;
	double fx=cosmo->survey.flux;

	double M0 = 1.0e15; //单位 M_sun
	double Ax = pow(10.0, log10Ax);
int status=0;
	double L = fx * 4.0 * M_PI * pow(clucu_luminosity_distance(cosmo,z,&status), 2.0);
	double Mbeta = L / (Ax * pow(clucu_h_over_h0(cosmo,z,&status), 2.0) * pow(1.0 + z, gammax));
	double M200 = pow(Mbeta, 1.0 / betax) * M0; //输出质量，单位M_sun
	double M180 = M200cz_to_M180mz(cosmo,M200, z);
	if (M180 < 1.0e14 / cosmo->param.h)
		M180 = 1.0e14 / cosmo->param.h;//质量Msun 和 1.0e14 h^-1 Msun比较
	return M180; //输出质量，单位M_sun
}
double Mmin_limit(clucu_cosmology *cosmo,double z)
{
	double M;
	if (cosmo->survey.log10ASZ < 0)
		M = Mmin_xray180(cosmo,z);
	else
		M = Mmin_SZE180(cosmo,z);
	return M;
}
