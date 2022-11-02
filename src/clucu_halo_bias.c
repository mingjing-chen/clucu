
#include "clucu_halo_bias.h"
//k: Mpc−1
double GISDB_factor(clucu_cosmology *cosmo, double k, double z)
{
	if(cosmo->GISDB==0) return 1.0;
	double alpha = 4.0;
	double Delta_LCDM = 4.8e-3;
	double theta_cmb = cosmo->param.T_CMB/2.7;
	double k_eq = 0.0746*cosmo->param.Omh2 / pow(theta_cmb,2.) ; // = 0.01Mpc-1 = 0.015  h Mpc−1
	
	double k_fs = (0.08/sqrt(1.0+z)) * (cosmo->param.m_nu_sum/0.1) * cosmo->param.h; // Mpc−1
//	double f_nu = cosmo->param.Onuh2_mass / (cosmo->param.Omh2 + cosmo->param.Onuh2_mass);
    double f_nu = cosmo->param.Onuh2_mass_sum / cosmo->param.Omh2;
	double q = 5.0*k/k_fs;
	double Delta_L = 0.6*f_nu;
	double Delta_q = 1.6;

	double R_LCDM = 1.0 + Delta_LCDM * tanh(alpha*k/k_eq);
	double R_nu = 1.0 + 0.5 * Delta_L *( tanh(log(q)/Delta_q)+1.0 );

	return R_LCDM*R_nu;
}
double NGbias_factor(clucu_cosmology *cosmo, double k, double z)
{
    int status=0;
	double H0 = clucu_hubble_parameter(cosmo,0.,&status); //1/Mpc
    double Om = cosmo->param.Omega_m;
    double fnl = cosmo->param.f_nl;
    double delta_c = 1.686;
    double Dz = clucu_growth_factor(cosmo,z,&status);
//TODO: 有没有T（k）?
    double Tk = clucu_transfer(cosmo, k,&status);
    return fnl * delta_c * 3.*Om*H0*H0 /Tk /Dz /k /k;
}
double HaloBiasPress1974(clucu_cosmology *cosmo,double sigma)
{
    double delta_crit=1.686;
	double nu2 = pow(delta_crit / sigma, 2.0);
	double bias = 1.0 + (nu2-1.) / delta_crit;

	return bias;
}

//178m
double HaloBiasSheth1999(clucu_cosmology *cosmo,double sigma)
{
	double a = 0.707;
	double p = 0.3;
    double delta_crit = 1.68647;
	double anu2 = a * pow(delta_crit / sigma, 2.0);
	double bias = 1.0 + (anu2 - 1.0) / delta_crit + 2.0 * p / delta_crit / (1.0 + pow(anu2, p));

	return bias;
}

//这个应该是fof0.2, SMT2001
double HaloBiasSheth2001(clucu_cosmology *cosmo,double sigma)
{
	double a = 0.707;
    double b=0.5;
    double c=0.6;
    double delta_crit = 1.68647;
    double nu = delta_crit/sigma;
    double anu2 = a * nu*nu;
    double term1 = sqrt(a)*anu2;
    double term2 = sqrt(a)*b*pow(anu2,1.0-c);
    double term3 = -pow(anu2,c) / ( pow(anu2,c)+b*(1.0-c)*(1.0-0.5*c)  );
	double bias = 1.0 + (term1+term2+term3) / (sqrt(a)*delta_crit);

	return bias;
}

//200m,300m,400m,600m,800m..
double HaloBiasTinker2010(clucu_cosmology *cosmo,double sigma)
{
    if(cosmo->config.halo_define_method != clucu_delta_200m)
        printf("wrong type in HaloBiasTinker10");
    double dc = 1.68647;
    double delta = 200;
    double ld=log10(delta);
    double xp = exp( -pow(4./ld,4.) );
    double A = 1.0 + 0.24 * ld * xp;
    double a = 0.44 * ld - 0.88;
    double B = 0.183;
    double b = 1.5;
    double C = 0.019 + 0.107 * ld + 0.19*xp;
    double c = 2.4;
    
    double nu = dc/sigma;
    double nua = pow(nu,a);

    return 1. - A * nua / (nua + pow(dc,a)) + B * pow(nu,b) + C * pow(nu,c);
}
//fof 0.2
double HaloBiasBhattacharya2011(clucu_cosmology *cosmo,double sigma)
{
    double a = 0.788;
    double az = 0.01;
    double p = 0.807;
    double q = 1.795;
    double dc = 1.68647;
    double nu = dc/sigma;
    a=a*pow(a,az);
    double anu2 = a*nu*nu;
    double bias = 1.0 + (anu2 - q) / dc + 2.0 * p / dc / (1.0 + pow(anu2, p));
	return bias;
}
double clucu_halo_bias_sigma(clucu_cosmology *cosmo,double sigma)
{
    switch(cosmo->config.bias_function_method)
    {
        case clucu_bias_press1974:          return HaloBiasPress1974(cosmo,sigma);break;
        case clucu_bias_sheth1999:          return HaloBiasSheth1999(cosmo,sigma);break;
        case clucu_bias_sheth2001:          return HaloBiasSheth2001(cosmo,sigma);break;
        case clucu_bias_tinker2010:         return HaloBiasTinker2010(cosmo,sigma);break;
        case clucu_bias_bhattacharya2011:   return HaloBiasBhattacharya2011(cosmo,sigma);break;
    }
}