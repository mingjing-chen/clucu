
#include "clucu_halo_concentration.h"
/**
 * 对于M–c relation，McClintock说用Bhattacharya et al. 2013
 * 或者 Diemer & Kravtsov 2015 没什么区别；
 * 
 * */
double ConcentrationDuffy2008(clucu_cosmology *cosmo,halo_define_type type,double mass,double z,int *status)
{
    double A,B,C;
    switch(type)
    {
        //case clucu_delta_vir: A=7.85;B=-0.081;C=-0.71;break;
        case clucu_delta_vir: 
            if(cosmo->survey.A_vir>0){
                A=cosmo->survey.A_vir;
                B=cosmo->survey.B_vir;
                C=cosmo->survey.C_vir;
            }
            break;
        case clucu_delta_200m: A=10.14;B=-0.081;C=-1.01;break;
        case clucu_delta_200c: A=5.71;B=-0.084;C=-0.47;break;
        default: printf("wrong type in ConcentrationDuffy2008\n");
    }
    double Mp = 2.0e12 /cosmo->param.h; //Msun/h改为Msun
    return A * pow(mass/Mp,B) * pow(1.+z,C);
}
double ConcentrationBhattacharya2013(clucu_cosmology *cosmo,halo_define_type type,double mass,double z,int *status)
{
    double A,B,C;
    switch(type)
    {
        case clucu_delta_vir: A=7.7;B=0.9;C=-0.29;break;
        case clucu_delta_200m: A=9.0;B=1.15;C=-0.29;break;
        case clucu_delta_200c: A=5.9;B=0.54;C=-0.35;break;
        default: printf("wrong type in ConcentrationBhattacharya2013\n");
    }    
    double delta_c = 1.686;
    double gz = clucu_growth_factor(cosmo,z,status);
    double sigma = clucu_sigma(cosmo, mass, z,status);//和质量定义无关
    double nu = delta_c / sigma;
    return A * pow(gz,B) * pow(nu,C);
}
double ConcentrationDiemer2015(clucu_cosmology *cosmo,halo_define_type type,double mass,double z,int *status)
{
    if(type != clucu_delta_200c)
        printf("wrong type in ConcentrationDiemer2015\n");
    double kappa = 1.0;
    double phi_0 = 6.58;
    double phi_1 = 1.27;
    double eta_0 = 7.28;
    double eta_1 = 1.56;
    double alpha = 1.08;
    double beta = 1.77;
    double delta_c = 1.686;
    //计算物质功率谱的斜率
    double R = M_to_R( cosmo, mass,status);
    double lk_R = log(2.0 * M_PI / R * kappa);
        // Using central finite differences
    double lk_hi = lk_R + 0.005;
    double lk_lo = lk_R - 0.005;
    double dlpk = log(clucu_power(cosmo, exp(lk_hi), z,status) /
                      clucu_power(cosmo, exp(lk_lo), z,status));
    double dlk = lk_hi - lk_lo;
    double n = dlpk / dlk;
    double sigma = clucu_sigma(cosmo, mass, z,status);//和质量定义无关
    double nu = delta_c / sigma;
    double c_min = phi_0 + n * phi_1;
    double nu_min = eta_0 + n * eta_1;
    return 0.5 * c_min * ( pow(nu_min / nu,alpha) + pow(nu / nu_min,beta) );
}

double clucu_concentration(clucu_cosmology *cosmo,halo_define_type type,double mass, double z,int *status)
{
	switch(cosmo->config.halo_concentration_method)
	{
		case clucu_concentration_duffy2008:return ConcentrationDuffy2008(cosmo,type,mass,z,status);break;
		case clucu_concentration_bhattacharya2013:return ConcentrationBhattacharya2013(cosmo,type,mass,z,status);break;
		case clucu_concentration_diemer2015:return ConcentrationDiemer2015(cosmo,type,mass,z,status);break;
        default: printf("wrong type in clucu_concentration\n");
    }
}
