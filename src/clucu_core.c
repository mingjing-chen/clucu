#include "clucu_core.h"

/** @file cosmo.c 
 *
 * @author Mingjing Chen
 * @description: 创建cosmology，初始化
 *
 */
const clucu_configuration config_default = 
{
    .matter_power_spectrum_method = clucu_boltzmann_class_Pk,
    .mass_function_method =         clucu_massfunc_tinker2008,
    .bias_function_method =         clucu_bias_sheth2001,
    .halo_concentration_method =    clucu_concentration_duffy2008,
    .mass_observable_method =       clucu_murrata2019,
    .halo_define_method =           clucu_delta_200m,
};
const clucu_configuration config_oguri = 
{
    .matter_power_spectrum_method = clucu_eisenstein_hu,
    .mass_function_method =         clucu_massfunc_bhattacharya2011,//不看，这个好像有问题，和tinker差好多
    .bias_function_method =         clucu_bias_bhattacharya2011,//不看
    .halo_concentration_method =    clucu_concentration_duffy2008,//看Delta类型
    .mass_observable_method =       clucu_oguri2011,
    .halo_define_method =           clucu_delta_vir,
};




//Precision parameters
/**
 * Default relative precision if not otherwise specified
 */
#define GSL_EPSREL 1E-4

/**
 * Default number of iterations for integration and root-finding if not otherwise
 * specified
 */
#define GSL_N_ITERATION 1000

/**
 * Default number of Gauss-Kronrod points in QAG integration if not otherwise
 * specified
 */
#define GSL_INTEGRATION_GAUSS_KRONROD_POINTS GSL_INTEG_GAUSS41

/**
 * Relative precision in sigma_R calculations
 */
#define GSL_EPSREL_SIGMAR 1E-7

/**
 * Relative precision in k_NL calculations
 */
#define GSL_EPSREL_KNL 1E-5

/**
 * Relative precision in distance calculations
 */
#define GSL_EPSREL_DIST 1E-6

/**
 * Relative precision in growth calculations
 */
#define GSL_EPSREL_GROWTH 1E-6

/**
 * Relative precision in dNdz calculations
 */
#define GSL_EPSREL_DNDZ 1E-6
const clucu_gsl_param gsl_param_default = 
{
    GSL_N_ITERATION,                     // N_ITERATION
    GSL_INTEGRATION_GAUSS_KRONROD_POINTS,// INTEGRATION_GAUSS_KRONROD_POINTS
    GSL_EPSREL,                          // INTEGRATION_EPSREL
    GSL_INTEGRATION_GAUSS_KRONROD_POINTS,// INTEGRATION_LIMBER_GAUSS_KRONROD_POINTS
    GSL_EPSREL,                          // INTEGRATION_LIMBER_EPSREL
    GSL_EPSREL_DIST,                     // INTEGRATION_DISTANCE_EPSREL
    GSL_EPSREL_SIGMAR,                   // INTEGRATION_SIGMAR_EPSREL
    GSL_EPSREL_KNL,                      // INTEGRATION_KNL_EPSREL
    GSL_EPSREL,                          // ROOT_EPSREL
    GSL_N_ITERATION,                     // ROOT_N_ITERATION
    GSL_EPSREL_GROWTH,                   // ODE_GROWTH_EPSREL
    1E-6,                                // EPS_SCALEFAC_GROWTH, 当a<该值时，增长因子D(a)=a
    1E7,                                 // HM_MMIN
    1E18,                                // HM_MMAX
    0.0,                                 // HM_EPSABS
    1E-4,                                // HM_EPSREL
    1000,                                // HM_LIMIT
    GSL_INTEG_GAUSS41,                   // HM_INT_METHOD
    true,                                // NZ_NORM_SPLINE_INTEGRATION
    true,                                 // LENSING_KERNEL_SPLINE_INTEGRATION
    .INTEGRATION_MASS_TRUE_EPSREL=1.0e-8,
    .INTEGRATION_Z_TRUE_EPSREL=1.0e-5,
    .INTEGRATION_MASS_OB_EPSREL=1.0e-5,
    .INTEGRATION_Z_OB_EPSREL=1.0e-5,
};

#undef GSL_EPSREL
#undef GSL_N_ITERATION
#undef GSL_INTEGRATION_GAUSS_KRONROD_POINTS
#undef GSL_EPSREL_SIGMAR
#undef GSL_EPSREL_KNL
#undef GSL_EPSREL_DIST
#undef GSL_EPSREL_GROWTH
#undef GSL_EPSREL_DNDZ


const clucu_spline_param spline_param_default = {

    .Z_SPLINE_MIN=0.,
    //功率谱插值
    .Z_SPLINE_MINLOG_PK=2.,
    .Z_SPLINE_MAX_PK=3.,
    .Z_SPLINE_NLIN_PK=100,
    .Z_SPLINE_NLOG_PK=2,
    //sigma插值
    .Z_SPLINE_MINLOG_SM=2.,
    .Z_SPLINE_MAX_SM=3.,
    .Z_SPLINE_NLIN_SM=100,
    .Z_SPLINE_NLOG_SM=2,
    //对chi(z),D(z)插值，背景插值
    .Z_SPLINE_MINLOG_BG=2.,
    .Z_SPLINE_MAX_BG=10000.,
    .Z_SPLINE_NLIN_BG=250,
    .Z_SPLINE_NLOG_BG=250,
  
    // mass splines
    .LOG10M_SPLINE_MIN=5,  // LOG10M_SPLINE_MIN 6
    .LOG10M_SPLINE_MAX=17,  // LOG10M_SPLINE_MAX 17
    .LOG10M_SPLINE_NM=100,  // LOG10M_SPLINE_NM 50：sigmar=0.937088；100:sigmar=0.937088；1000:sigmar=0.937088

  
  

  // k-splines and integrals
  .K_MAX_SPLINE=50,  // K_MAX_SPLINE
  .K_MAX=1E3,  // K_MAX 功率谱插值的最大值
  .K_MIN=5E-5,  // K_MIN 功率谱插值的最小值
  .DLOGK_INTEGRATION=0.025,  // DLOGK_INTEGRATION
  .DCHI_INTEGRATION=5.,  // DCHI_INTEGRATION。对chi积分时，使用插值积分法时
  .N_K=167,  // N_K 每个log10k 内，有插值几个
  .N_K_3DCOR=100000,  // N_K_3DCOR

  // correlation function parameters
  .ELL_MIN_CORR=0.01,  // ELL_MIN_CORR
  .ELL_MAX_CORR=60000,  // ELL_MAX_CORR
  .N_ELL_CORR=5000,  // N_ELL_CORR

  //Spline types
  //static const gsl_interp_type cspline_type 定义一个静态结构体cspline_type
  //const gsl_interp_type * gsl_interp_cspline = &cspline_type;将结构体的地址赋值给指针
  .Z_SPLINE_TYPE=&gsl_interp_cspline,
  .K_SPLINE_TYPE=NULL,
  .M_SPLINE_TYPE=NULL,
  .D_SPLINE_TYPE=NULL,
  .PNL_SPLINE_TYPE=NULL,
  .PLIN_SPLINE_TYPE=&gsl_interp2d_bicubic,
  .CORR_SPLINE_TYPE=NULL,
};
const clucu_parameter param_default={
    .Neff = 3.046,
    .N_nu_mass =1,
    .Onuh2_mass_sum = 6.442e-4,
    .neutrino_mass_split_label = clucu_nu_single,

    .H0 = NAN,
    .h = NAN,
    .Obh2 = 0.022,
    .Omh2 = 0.1430,
    .Omega_l = 0.6847,
    .w0 = -1.0,
    .wa = 0.0,
    .Omega_k = 0.0,

    .n_s = 0.9665,
    .alpha_s = 0.,
    .k_pivot = 0.05,
    .sigma8 = 0.8102,
    .delta_zeta = NAN,
    .A_s = NAN, 
    //.sigma8 = NAN,
    //.delta_zeta = 4.551473387816302e-05,
};
const clucu_parameter param_planck2018={
    .Neff = 3.046,
    .N_nu_mass =0,
    .Onuh2_mass = 0.,
    //.m_nu = {0.02,0.02,0.02},

    .H0 = NAN,
    .h = NAN,
    .Obh2 = 0.0224,
    .Och2 = 0.120,
    .Omega_l = 0.6847,
    .w0 = -1.0,
    .wa = 0.0,
    .Omega_k = 0.0,

    .n_s = 0.965,
    .alpha_s = 0.,
    .k_pivot = 0.05,
    .sigma8 = 0.811,
    .delta_zeta = NAN,
    .A_s = NAN, 
    //.sigma8 = NAN,
    //.delta_zeta = 4.551473387816302e-05,
};
const clucu_parameter param_wmap9={
    .Neff = 3.046,
    .N_nu_mass =0,
    .Onuh2_mass = 0.,
    //.m_nu = {0.02,0.02,0.02},

    .H0 = NAN,
    .h = NAN,
    .Obh2 = 0.02264,
    .Och2 = 0.1138,
    .Omega_l = 0.721,
    .w0 = -1.0,
    .wa = 0.0,
    .Omega_k = 0.0,

    .n_s = 0.972,
    .alpha_s = 0.,
    .k_pivot = 0.05,
    .sigma8 = 0.821,
    .delta_zeta = NAN,
    .A_s = NAN, 
    //.sigma8 = NAN,
    //.delta_zeta = 4.551473387816302e-05,
};
const clucu_parameter param_WMAP5={
    .Neff = 3.046,
    .N_nu_mass =0,
    .Onuh2_mass = 0.,
//    .m_nu = {0.},

    .H0 = NAN,
    .h = 0.7,
    .Omega_b = 0.0462,
    .Omega_m = 0.279,

    .w0 = -1.0,
    .wa = 0.0,
    .Omega_k = 0.0,

    .n_s = 0.96,
    .alpha_s = 0.,
    .k_pivot = 0.02,
    .sigma8 = NAN,
    .delta_zeta = sqrt(2.21e-9),
    //对应的CLASS算出来的.sigma8 = 0.844482,

    .tau_reio = 0.089,
};
const clucu_parameter param_oguri={
    .Neff = 3.046,
    .N_nu_mass =0,
    .Onuh2_mass = 0.,
//    .m_nu = {0.},

    .H0 = NAN,
    .h = NAN,
    .Obh2 = 0.0226,
    .Omh2 = 0.134,
    .Omega_l = 0.734,
    .w0 = -1.0,
    .wa = 0.0,
    .Omega_k = 0.0,
    

    .n_s = 0.963,
    .alpha_s = 0.,
    .k_pivot = 0.002,
    .sigma8 = NAN,
    .delta_zeta = 4.89e-5,
    //对应的CLASS算出来的.sigma8 = 0.844482,

    .tau_reio = 0.089,

};

const clucu_survey survey_sze={
    .flux=5.0,
	.log10ASZ=8.9,
	.betaSZ=1.68,
	.gammaSZ=0.0,
    .survey_area = 4000.,
    .redshift_sigma0 = 0.00,
    .mass_ob_min = NAN,
    .spline_z_minlog = 2.0,
    .spline_z_max = 2.2,
    .spline_z_nlin = 80,
    .spline_z_nlog = 3,
    .zob_min = 0.0,
    .zob_max = 2.0,
    .zob_bin = 0.05,
    .zobP_min = 0.0,
    .zobP_max = 2.0,
    .zobP_bin = 0.2,
    .k1_min = 0.005,
    .k1_max = 0.15,
    .k1_bin = 0.005,
    .k2_min = 0.005,
    .k2_max = 0.15,
    .k2_bin = 0.005,
};

const clucu_survey survey_xrayDeep={
    .flux=2.25e-14,
	.log10ASZ=-4.159,
	.betaSZ=1.807,
	.gammaSZ=0.0,
    .survey_area = 150.0,
    .redshift_sigma0 = 0.00,
    .mass_ob_min = NAN,
    .spline_z_minlog = 2.0,
    .spline_z_max = 2.2,
    .spline_z_nlin = 80,
    .spline_z_nlog = 3,
    .zob_min = 0.0,
    .zob_max = 2.0,
    .zob_bin = 0.05,
    .zobP_min = 0.0,
    .zobP_max = 2.0,
    .zobP_bin = 0.2,
    .k1_min = 0.005,
    .k1_max = 0.15,
    .k1_bin = 0.005,
    .k2_min = 0.005,
    .k2_max = 0.15,
    .k2_bin = 0.005,
};

const clucu_survey survey_xrayWide={
    .flux=1.75e-13,
	.log10ASZ=-4.159,
	.betaSZ=1.807,
	.gammaSZ=0.0,
    .survey_area = 6000.0,
    .redshift_sigma0 = 0.00,
    .mass_ob_min = NAN,
    .spline_z_minlog = 2.0,
    .spline_z_max = 2.2,
    .spline_z_nlin = 80,
    .spline_z_nlog = 3,
    .zob_min = 0.0,
    .zob_max = 2.0,
    .zob_bin = 0.05,
    .zobP_min = 0.0,
    .zobP_max = 2.0,
    .zobP_bin = 0.2,
    .k1_min = 0.005,
    .k1_max = 0.15,
    .k1_bin = 0.005,
    .k2_min = 0.005,
    .k2_max = 0.15,
    .k2_bin = 0.005,
};
const clucu_survey survey_csst={
    .mass_ob_A = 3.15,
    .mass_ob_B = 0.86, 
    .mass_ob_Bz = -0.21, 
    .mass_ob_Cz = 0.0, 
    .mass_ob_sigma0 = 0.32, 
    .mass_ob_q = -0.08,
    .mass_ob_qz = 0.03,
    .mass_ob_pz = 0.00,
    .survey_area = 17500.0,
    .redshift_sigma0 = 0.001,
    .mass_ob_min = 15,
    .mass_ob_max = 300,
    .dln_mass_ob = 0.3,
    .dlog10_mass_ob = NAN,
    .spline_z_minlog = 1.52,
    .spline_z_max = 1.7,
    .spline_z_nlin = 80,
    .spline_z_nlog = 3,
    .zob_min = 0.0,
    .zob_max = 1.5,
    .zob_bin = 0.05,
    .zobP_min = 0.0,
    .zobP_max = 1.4,
    .zobP_bin = 0.2,
    .k1_min = 0.005,
    .k1_max = 0.15,
    .k1_bin = 0.005,
    .k2_min = 0.005,
    .k2_max = 0.15,
    .k2_bin = 0.005,
};
const clucu_survey survey_csst2={
    .mass_ob_A = 79.8,
    .mass_ob_B = 0.93, 
    .mass_ob_Bz = -0.49, 
    .mass_ob_sigma0 = 0.217, 
    .mass_ob_q = 0.01,//kappa
    .mass_ob_qz = 1.,//beta=1
    .lnM_b0 = 0.,//B_M0
    .s_b = {0.,0.,0.},//alpha
    .survey_area = 17500.0,
    .redshift_sigma0 = 0.001,
    .mass_ob_min = NAN,
    .mass_ob_max = NAN,
    .mass_ob_min_h = 1.0e14,
    .mass_ob_max_h = 1.0e16,
    .dln_mass_ob = NAN,
    .dlog10_mass_ob = 0.2,

    .spline_z_minlog = 1.52,
    .spline_z_max = 1.7,
    .spline_z_nlin = 80,
    .spline_z_nlog = 3,

    .zob_min = 0.0,
    .zob_max = 1.5,
    .zob_bin = 0.05,
    .zobP_min = 0.0,
    .zobP_max = 1.5,
    .zobP_bin = 0.1,
    .k1_min = 0.005,
    .k1_max = 0.15,
    .k1_bin = 0.005,
    .k2_min = 0.005,
    .k2_max = 0.15,
    .k2_bin = 0.005,
};
const clucu_survey survey_oguri={
    //A=7.85;B=-0.081;C=-0.71
    .A_vir = 7.85,
    .B_vir = -0.081,
    .C_vir = -0.71,
    .f_cen0 = 0.75,
    .p_cenM = 0.05,
    .p_cenz = 0.,
    .sigma_s0 = 0.42,
    .p_sigmaM = 0.,
    .p_sigmaz = 0.,
    .lnM_b0 = 0.,
    .q_b = {0.,0.,0.},
    .s_b = {0.,0.,0.},
    .sigma_lnM0 = 0.3,
    .q_sigma_lnm = {0.,0.,0.},
    .s_sigma_lnm = {0.,0.,0.},
    .source_redshift = 1.,


    .mass_ob_A = 0.,
    .mass_ob_sigma0 = 0.3, //这个不能删啊，还得用这个判断有没有pm
    .survey_area = 2000.0,
    .redshift_sigma0 = 0.0,
    .mass_ob_min = NAN,
    .mass_ob_max = NAN,
    .mass_ob_min_h = 1.0e14,
    .mass_ob_max_h = 1.0e16,
    .dln_mass_ob = NAN,
    .dlog10_mass_ob = 0.2,
    
    .spline_z_minlog = 1.52,
    .spline_z_max = 1.7,
    .spline_z_nlin = 80,
    .spline_z_nlog = 3,
    
    .zob_min = 0.,
    .zob_max = 1.4,
    .zob_bin = 0.1,
    .zobP_min = 0.,
    .zobP_max = 1.4,
    .zobP_bin = 0.1,

    .k1_min = 0.005,
    .k1_max = 0.15,
    .k1_bin = 0.005,
    .k2_min = 0.005,
    .k2_max = 0.15,
    .k2_bin = 0.005,
};

void clucu_parameters_initial_knowH(clucu_cosmology *cosmo,int *status)
{
    //曲率
    if(fabs(cosmo->param.Omega_k)<1E-6)
        cosmo->param.k_sign=0;
    else if(cosmo->param.Omega_k>0)
        cosmo->param.k_sign=-1;
    else
        cosmo->param.k_sign=1;
    cosmo->param.sqrtk=sqrt(fabs(cosmo->param.Omega_k))*cosmo->param.h/clucu_constants.CLIGHT_HMPC;

    //功率谱
    if(isnormal(cosmo->param.A_s) && !isnormal(cosmo->param.delta_zeta))
        cosmo->param.delta_zeta = sqrt(cosmo->param.A_s);
    else if(!isnormal(cosmo->param.A_s) && isnormal(cosmo->param.delta_zeta))
        cosmo->param.A_s = cosmo->param.delta_zeta * cosmo->param.delta_zeta;

    //判断是已知h还是已知H0
    if(isnormal(cosmo->param.h) && !isnormal(cosmo->param.H0))
        cosmo->param.H0 = cosmo->param.h *100.;
    else if(!isnormal(cosmo->param.h) && isnormal(cosmo->param.H0))
        cosmo->param.h = cosmo->param.H0 /100.;


    //光子
    cosmo->param.T_CMB =  clucu_constants.T_CMB;
    double rho_g = 4. * clucu_constants.STBOLTZ / pow(clucu_constants.CLIGHT, 3) * pow(cosmo->param.T_CMB, 4);// kg / m^3
    double rho_crit = clucu_constants.RHO_CRITICAL * clucu_constants.SOLAR_MASS/pow(clucu_constants.MPC_TO_METER, 3) *  pow(cosmo->param.h, 2);// kg / m^3
    cosmo->param.Omega_g = rho_g/rho_crit;
    

    //中微子
    //非相对中微子
    if(cosmo->param.neutrino_mass_split_label == clucu_nu_splited){
        if(isnormal(cosmo->param.m_nu[0]) && !isnormal(cosmo->param.Omega_nu_mass[0]) && !isnormal(cosmo->param.Onuh2_mass[0])){
            for(int i=0;i<3;i++){
                if(cosmo->param.m_nu[i]>0){
                    cosmo->param.N_nu_mass++;
                    cosmo->param.Onuh2_mass[i] = density_WDM(cosmo->param.m_nu[i],cosmo->param.T_CMB*clucu_constants.TNCDM)/clucu_constants.RHO_CRITICAL;
                    cosmo->param.Omega_nu_mass[i] = cosmo->param.Onuh2_mass[i] / cosmo->param.h/cosmo->param.h;
                    cosmo->param.m_nu_sum += cosmo->param.m_nu[i];
                    cosmo->param.Omega_nu_mass_sum += cosmo->param.Omega_nu_mass[i];
                    cosmo->param.Onuh2_mass_sum += cosmo->param.Onuh2_mass[i];
                }
            }
        }
        else if(!isnormal(cosmo->param.m_nu[0]) && isnormal(cosmo->param.Omega_nu_mass[0]) && !isnormal(cosmo->param.Onuh2_mass[0])){
            for(int i=0;i<3;i++){
                if(cosmo->param.Omega_nu_mass[i]>0){
                    cosmo->param.N_nu_mass++;
                    cosmo->param.Onuh2_mass[i] = cosmo->param.Omega_nu_mass[i] *cosmo->param.h *cosmo->param.h;
                    cosmo->param.m_nu[i] = cosmo->param.Onuh2_mass[i]/93.14;
                    cosmo->param.m_nu_sum += cosmo->param.m_nu[i];
                    cosmo->param.Omega_nu_mass_sum += cosmo->param.Omega_nu_mass[i];
                    cosmo->param.Onuh2_mass_sum += cosmo->param.Onuh2_mass[i];
                }
            }
        } 
        else if(!isnormal(cosmo->param.m_nu[0]) && !isnormal(cosmo->param.Omega_nu_mass[0]) && isnormal(cosmo->param.Onuh2_mass[0])){
            for(int i=0;i<3;i++){
                if(cosmo->param.Onuh2_mass[i]>0){
                    cosmo->param.N_nu_mass++;
                    cosmo->param.Omega_nu_mass[i] = cosmo->param.Onuh2_mass[i] / cosmo->param.h/cosmo->param.h;
                    cosmo->param.m_nu[i] = cosmo->param.Onuh2_mass[i]/93.14;
                    cosmo->param.m_nu_sum += cosmo->param.m_nu[i];
                    cosmo->param.Omega_nu_mass_sum += cosmo->param.Omega_nu_mass[i];
                    cosmo->param.Onuh2_mass_sum += cosmo->param.Onuh2_mass[i];
                }
            }
        }
    }
    else{
        if(isnormal(cosmo->param.m_nu_sum) && !isnormal(cosmo->param.Omega_nu_mass_sum) && !isnormal(cosmo->param.Onuh2_mass_sum)){
            cosmo->param.Omega_nu_mass_sum = cosmo->param.m_nu_sum / 93.14 /cosmo->param.h /cosmo->param.h;
            cosmo->param.Onuh2_mass_sum = cosmo->param.m_nu_sum / 93.14;
        }
        else if(!isnormal(cosmo->param.m_nu_sum) && isnormal(cosmo->param.Omega_nu_mass_sum) && !isnormal(cosmo->param.Onuh2_mass_sum)){
            cosmo->param.m_nu_sum = cosmo->param.Omega_nu_mass_sum *93.14 *cosmo->param.h *cosmo->param.h;
            cosmo->param.Onuh2_mass_sum = cosmo->param.Omega_nu_mass_sum *cosmo->param.h *cosmo->param.h;
        }
        else if(!isnormal(cosmo->param.m_nu_sum) && !isnormal(cosmo->param.Omega_nu_mass_sum) && isnormal(cosmo->param.Onuh2_mass_sum)){
            cosmo->param.m_nu_sum = cosmo->param.Onuh2_mass_sum *93.14;
            cosmo->param.Omega_nu_mass_sum = cosmo->param.Onuh2_mass_sum /cosmo->param.h /cosmo->param.h;
        }
        clucu_neutrino_mass_split(cosmo->param.m_nu_sum, cosmo->param.neutrino_mass_split_label, cosmo->param.m_nu, &(cosmo->param.N_nu_mass), status);
        for(int i=0;i<cosmo->param.N_nu_mass;i++){
            cosmo->param.Onuh2_mass[i] = cosmo->param.m_nu[i] / 93.14;
            cosmo->param.Omega_nu_mass[i] = cosmo->param.Onuh2_mass[i] / cosmo->param.h/cosmo->param.h;
        }
    }
    //相对论中微子
    cosmo->param.N_nu_rel = cosmo->param.Neff - cosmo->param.N_nu_mass * pow(clucu_constants.TNCDM, 4) / pow(4./11.,4./3.);
    double T_nu= (cosmo->param.T_CMB) * pow(4./11.,1./3.);
    double rho_nu_rel = cosmo->param.N_nu_rel* 7.0/8.0 * 4. * clucu_constants.STBOLTZ / pow(clucu_constants.CLIGHT, 3) * pow(T_nu, 4);// in kg / m^3，相对论中微子的密度
    cosmo->param.Omega_nu_rel = rho_nu_rel/rho_crit;
    cosmo->param.Onuh2_rel = cosmo->param.Omega_nu_rel * cosmo->param.h * cosmo->param.h;

//判断是已知Om，还是已知Omh2。改为Om
    if( (!isnormal(cosmo->param.Omega_b) && isnormal(cosmo->param.Obh2)) || (!isnormal(cosmo->param.Omega_c) && isnormal(cosmo->param.Och2)) || (!isnormal(cosmo->param.Omega_m) && isnormal(cosmo->param.Omh2)) )
    {
        cosmo->param.Omega_c = cosmo->param.Och2 / cosmo->param.h / cosmo->param.h;
        cosmo->param.Omega_b = cosmo->param.Obh2 / cosmo->param.h / cosmo->param.h;
        cosmo->param.Omega_m = cosmo->param.Omh2 / cosmo->param.h / cosmo->param.h;
    }
    //根据b,c,m中的任意两个，求第三个：
    if(isnormal(cosmo->param.Omega_b) && isnormal(cosmo->param.Omega_c))
        cosmo->param.Omega_m = cosmo->param.Omega_b + cosmo->param. Omega_c + cosmo->param.Omega_nu_mass_sum;
    else if(isnormal(cosmo->param.Omega_m) && isnormal(cosmo->param.Omega_c))
        cosmo->param.Omega_b = cosmo->param.Omega_m - cosmo->param. Omega_c -cosmo->param.Omega_nu_mass_sum;
    else if(isnormal(cosmo->param.Omega_m) && isnormal(cosmo->param.Omega_b))
        cosmo->param.Omega_c = cosmo->param.Omega_m - cosmo->param. Omega_b -cosmo->param.Omega_nu_mass_sum;
    cosmo->param.Ogh2 = cosmo->param.Omega_g * cosmo->param.h * cosmo->param.h;
    cosmo->param.Och2 = cosmo->param.Omega_c * cosmo->param.h * cosmo->param.h;
    cosmo->param.Obh2 = cosmo->param.Omega_b * cosmo->param.h * cosmo->param.h;
    cosmo->param.Omh2 = cosmo->param.Omega_m * cosmo->param.h * cosmo->param.h;

    cosmo->param.Omega_l = 1.0 - cosmo->param.Omega_m - cosmo->param.Omega_g - cosmo->param.Omega_nu_rel - cosmo->param.Omega_k;
}
//已知Omega_L，已知Omh2,中微子Onuh2或mnu，未知h
void clucu_parameters_initial_knowOL(clucu_cosmology *cosmo,int *status)
{
    //曲率
    if(fabs(cosmo->param.Omega_k)<1E-6)
        cosmo->param.k_sign=0;
    else if(cosmo->param.Omega_k>0)
        cosmo->param.k_sign=-1;
    else
        cosmo->param.k_sign=1;
    cosmo->param.sqrtk=sqrt(fabs(cosmo->param.Omega_k))*cosmo->param.h/clucu_constants.CLIGHT_HMPC;
    
    //功率谱
    if(isnormal(cosmo->param.A_s) && !isnormal(cosmo->param.delta_zeta))
        cosmo->param.delta_zeta = sqrt(cosmo->param.A_s);
    else if(!isnormal(cosmo->param.A_s) && isnormal(cosmo->param.delta_zeta))
        cosmo->param.A_s = cosmo->param.delta_zeta * cosmo->param.delta_zeta;

    //光子
    cosmo->param.T_CMB =  clucu_constants.T_CMB;
    double rho_g = 4. * clucu_constants.STBOLTZ / pow(clucu_constants.CLIGHT, 3) * pow(cosmo->param.T_CMB, 4);// kg / m^3
    double rho_crit = clucu_constants.RHO_CRITICAL * clucu_constants.SOLAR_MASS/pow(clucu_constants.MPC_TO_METER, 3);//h^2 kg / m^3 
    cosmo->param.Ogh2 = rho_g/rho_crit;

    //中微子
    //非相对中微子
    if(cosmo->param.neutrino_mass_split_label == clucu_nu_splited){
        if(isnormal(cosmo->param.m_nu[0]) && !isnormal(cosmo->param.Omega_nu_mass[0]) && !isnormal(cosmo->param.Onuh2_mass[0])){
            for(int i=0;i<3;i++){
                if(cosmo->param.m_nu[i]>0){
                    cosmo->param.N_nu_mass++;
                    cosmo->param.Onuh2_mass[i] = density_WDM(cosmo->param.m_nu[i],cosmo->param.T_CMB*clucu_constants.TNCDM)/clucu_constants.RHO_CRITICAL;
//                    cosmo->param.Omega_nu_mass[i] = cosmo->param.Onuh2_mass[i] / cosmo->param.h/cosmo->param.h;
                    cosmo->param.m_nu_sum += cosmo->param.m_nu[i];
//                    cosmo->param.Omega_nu_mass_sum += cosmo->param.Omega_nu_mass[i];
                    cosmo->param.Onuh2_mass_sum += cosmo->param.Onuh2_mass[i];
                }
            }
        }
        else if(!isnormal(cosmo->param.m_nu[0]) && !isnormal(cosmo->param.Omega_nu_mass[0]) && isnormal(cosmo->param.Onuh2_mass[0])){
            for(int i=0;i<3;i++){
                if(cosmo->param.Onuh2_mass[i]>0){
                    cosmo->param.N_nu_mass++;
//                    cosmo->param.Omega_nu_mass[i] = cosmo->param.Onuh2_mass[i] / cosmo->param.h/cosmo->param.h;
                    cosmo->param.m_nu[i] = cosmo->param.Onuh2_mass[i]/93.14;
                    cosmo->param.m_nu_sum += cosmo->param.m_nu[i];
//                    cosmo->param.Omega_nu_mass_sum += cosmo->param.Omega_nu_mass[i];
                    cosmo->param.Onuh2_mass_sum += cosmo->param.Onuh2_mass[i];
                }
            }
        }
    }
    else{
        if(isnormal(cosmo->param.m_nu_sum) && !isnormal(cosmo->param.Omega_nu_mass_sum) && !isnormal(cosmo->param.Onuh2_mass_sum)){
//            cosmo->param.Omega_nu_mass_sum = cosmo->param.m_nu_sum / 93.14 /cosmo->param.h /cosmo->param.h;
            cosmo->param.Onuh2_mass_sum = cosmo->param.m_nu_sum / 93.14;
        }
        else if(!isnormal(cosmo->param.m_nu_sum) && !isnormal(cosmo->param.Omega_nu_mass_sum) && isnormal(cosmo->param.Onuh2_mass_sum)){
            cosmo->param.m_nu_sum = cosmo->param.Onuh2_mass_sum *93.14;
//            cosmo->param.Omega_nu_mass_sum = cosmo->param.Onuh2_mass_sum /cosmo->param.h /cosmo->param.h;
        }
        else
            CLUCU_RAISE_WARNING(1,"wrong in neutrino");
        clucu_neutrino_mass_split(cosmo->param.m_nu_sum, cosmo->param.neutrino_mass_split_label, cosmo->param.m_nu, &(cosmo->param.N_nu_mass), status);
        for(int i=0;i<cosmo->param.N_nu_mass;i++){
            cosmo->param.Onuh2_mass[i] = cosmo->param.m_nu[i] / 93.14;
//            cosmo->param.Omega_nu_mass[i] = cosmo->param.Onuh2_mass[i] / cosmo->param.h/cosmo->param.h;
        } 
    }
    //相对论中微子
    cosmo->param.N_nu_rel = cosmo->param.Neff - cosmo->param.N_nu_mass * pow(clucu_constants.TNCDM, 4) / pow(4./11.,4./3.);
    double T_nu= (cosmo->param.T_CMB) * pow(4./11.,1./3.);
    double rho_nu_rel = cosmo->param.N_nu_rel* 7.0/8.0 * 4. * clucu_constants.STBOLTZ / pow(clucu_constants.CLIGHT, 3) * pow(T_nu, 4);// in kg / m^3
    cosmo->param.Onuh2_rel = rho_nu_rel/rho_crit;


    //根据b,c,m中的任意两个，求第三个：
    if(isnormal(cosmo->param.Obh2) && isnormal(cosmo->param.Och2))
        cosmo->param.Omh2 = cosmo->param.Obh2 + cosmo->param.Och2 + cosmo->param.Onuh2_mass_sum;
    else if(isnormal(cosmo->param.Omh2) && isnormal(cosmo->param.Och2))
        cosmo->param.Obh2 = cosmo->param.Omh2 - cosmo->param.Och2 -cosmo->param.Onuh2_mass_sum;
    else if(isnormal(cosmo->param.Omh2) && isnormal(cosmo->param.Obh2))
        cosmo->param.Och2 = cosmo->param.Omh2 - cosmo->param.Obh2 -cosmo->param.Onuh2_mass_sum;
    else
        CLUCU_RAISE_WARNING(1,"wrong in b,c,m");


    cosmo->param.h = sqrt((cosmo->param.Omh2 + cosmo->param.Ogh2 + cosmo->param.Onuh2_rel ) / (1 - cosmo->param.Omega_l - cosmo->param.Omega_k));
    cosmo->param.H0 = 100.0 * cosmo->param.h;   
    for(int i=0;i<cosmo->param.N_nu_mass;i++)
        cosmo->param.Omega_nu_mass[i] = cosmo->param.Onuh2_mass[i] / pow(cosmo->param.h, 2.0);
    cosmo->param.Omega_nu_mass_sum = cosmo->param.Onuh2_mass_sum / pow(cosmo->param.h, 2.0);
    cosmo->param.Omega_nu_rel = cosmo->param.Onuh2_rel / pow(cosmo->param.h, 2.0);
    cosmo->param.Omega_g = cosmo->param.Ogh2 / pow(cosmo->param.h, 2.0);
    cosmo->param.Omega_b = cosmo->param.Obh2 / pow(cosmo->param.h, 2.0);
    cosmo->param.Omega_m = cosmo->param.Omh2 / pow(cosmo->param.h, 2.0);
    cosmo->param.Omega_c = cosmo->param.Och2 / pow(cosmo->param.h, 2.0);
    

    
}
//survey初始化
void InitializeCosmo(clucu_cosmology *cosmo)
{
    //param初始化
    if(!isnormal(cosmo->param.h) && !isnormal(cosmo->param.H0) && isnormal(cosmo->param.Omega_l))
        clucu_parameters_initial_knowOL(cosmo,&cosmo->status);
    else if( (isnormal(cosmo->param.h) || isnormal(cosmo->param.H0)) && !isnormal(cosmo->param.Omega_l))
        clucu_parameters_initial_knowH(cosmo,&cosmo->status);
     
    if(!isnormal(cosmo->survey.mass_ob_min) && isnormal(cosmo->survey.mass_ob_min_h) )
    {
        cosmo->survey.mass_ob_min = cosmo->survey.mass_ob_min_h / cosmo->param.h;
        cosmo->survey.mass_ob_max = cosmo->survey.mass_ob_max_h / cosmo->param.h;
        cosmo->survey.mass_ob_maxc = cosmo->survey.mass_ob_maxc_h / cosmo->param.h;
    }

    cosmo->data.spline_Nz = cosmo->survey.spline_z_nlin + cosmo->survey.spline_z_nlog - 1;
    cosmo->data.spline_z = clucu_linlog_spacing(0.,
                                                cosmo->survey.spline_z_minlog,
                                                cosmo->survey.spline_z_max,
                                                cosmo->survey.spline_z_nlin,
                                                cosmo->survey.spline_z_nlog);
    

    // cosmo->data.Nz_spline = cosmo->spline_param.Z_SPLINE_NLIN_PK + cosmo->spline_param.Z_SPLINE_NLOG_PK - 1;
    // cosmo->data.z_spline = clucu_linlog_spacing(cosmo->spline_param.Z_SPLINE_MIN,
    //                                             cosmo->spline_param.Z_SPLINE_MINLOG_PK,
    //                                             cosmo->spline_param.Z_SPLINE_MAX_PK,
    //                                             cosmo->spline_param.Z_SPLINE_NLIN_PK,
    //                                             cosmo->spline_param.Z_SPLINE_NLOG_PK);

    
    //根据survey初始化data
    cosmo->survey.Delta_Omega = cosmo->survey.survey_area*pow(M_PI/180,2);
//    cosmo->data.Nz_spline=(int) ( (cosmo->survey.z_max-cosmo->survey.z_min)/cosmo->survey.z_bin  +0.5)+1; //红移区间的个数
    cosmo->data.Nzob=(int) ( (cosmo->survey.zob_max-cosmo->survey.zob_min)/cosmo->survey.zob_bin +0.5 )+1; //红移区间的个数
	cosmo->data.NzobP=(int) ( (cosmo->survey.zobP_max-cosmo->survey.zobP_min)/cosmo->survey.zobP_bin  +0.5)+1; //红移区间的个数
    cosmo->data.Nk1=(int) ( (cosmo->survey.k1_max-cosmo->survey.k1_min)/cosmo->survey.k1_bin +0.5 )+1; //输入k垂直的个数
    cosmo->data.Nk2=(int) ( (cosmo->survey.k2_max-cosmo->survey.k2_min)/cosmo->survey.k2_bin +0.5 )+1; //输入k平行的个数
//    cosmo->data.z = clucu_linear_spacing(cosmo->survey.z_min, cosmo->survey.z_max, cosmo->data.Nz_spline);//红移数组;
    cosmo->data.zob = clucu_linear_spacing(cosmo->survey.zob_min, cosmo->survey.zob_max, cosmo->data.Nzob);//红移数组;
	cosmo->data.zobP = clucu_linear_spacing(cosmo->survey.zobP_min, cosmo->survey.zobP_max, cosmo->data.NzobP);//红移数组;
    cosmo->data.k1 = clucu_linear_spacing(cosmo->survey.k1_min, cosmo->survey.k1_max, cosmo->data.Nk1);//k垂直的数组
    cosmo->data.k2 = clucu_linear_spacing(cosmo->survey.k2_min, cosmo->survey.k2_max, cosmo->data.Nk2);//k平行的数
    cosmo->data.Nell=(int) ( (4.-0.)/0.1 +0.5 )+1; //红移区间的个数
    cosmo->data.ell = clucu_log10_spacing(1.,1.0e4,cosmo->data.Nell);
    if(cosmo->survey.mass_ob_min != NAN)//csst
    {
        if(cosmo->survey.dln_mass_ob >1.0e-5)
        {
            cosmo->data.Nmob=(int) ( ( log(cosmo->survey.mass_ob_max)-log(cosmo->survey.mass_ob_min) )/cosmo->survey.dln_mass_ob +0.5 )+1; //mass_ob区间的个数
            cosmo->data.ln_mass_ob = clucu_linear_spacing(log(cosmo->survey.mass_ob_min), log(cosmo->survey.mass_ob_max), cosmo->data.Nmob);//mass_ob的数组;
            cosmo->data.mass_ob = Create1Grid(cosmo->data.Nmob);
            cosmo->data.log10_mass_ob = Create1Grid(cosmo->data.Nmob);
            for(int number_mob=0;number_mob<cosmo->data.Nmob;number_mob++)
            {
                cosmo->data.mass_ob[number_mob] = exp(cosmo->data.ln_mass_ob[number_mob]);  //单位M_sun
                cosmo->data.log10_mass_ob[number_mob] = log10(cosmo->data.mass_ob[number_mob]);
            }
        }
        else if(cosmo->survey.dlog10_mass_ob >1.0e-5)
        {
            if(cosmo->survey.mass_ob_maxc > 1.0e-5)
            {
                cosmo->data.Nmob = (int) ( ( log10(cosmo->survey.mass_ob_maxc)-log10(cosmo->survey.mass_ob_min) )/cosmo->survey.dlog10_mass_ob +0.5 )+1 +1;
                cosmo->data.mass_ob = Create1Grid(cosmo->data.Nmob);
                cosmo->data.ln_mass_ob = Create1Grid(cosmo->data.Nmob);
                cosmo->data.log10_mass_ob = Create1Grid(cosmo->data.Nmob);
                for(int i=0;i<cosmo->data.Nmob-1;i++){
                    cosmo->data.log10_mass_ob[i] = log10(cosmo->survey.mass_ob_min) + cosmo->survey.dlog10_mass_ob *i;
                }
                cosmo->data.log10_mass_ob[cosmo->data.Nmob-1] = log10(cosmo->survey.mass_ob_max);
            }
            else{
                cosmo->data.Nmob=(int) ( ( log10(cosmo->survey.mass_ob_max)-log10(cosmo->survey.mass_ob_min) )/cosmo->survey.dlog10_mass_ob +0.5 )+1; //mass_ob区间的个数
                cosmo->data.log10_mass_ob = clucu_linear_spacing(log10(cosmo->survey.mass_ob_min), log10(cosmo->survey.mass_ob_max), cosmo->data.Nmob);//mass_ob的数组;
                cosmo->data.mass_ob = Create1Grid(cosmo->data.Nmob);
                cosmo->data.ln_mass_ob = Create1Grid(cosmo->data.Nmob);
            }
            
            for(int number_mob=0;number_mob<cosmo->data.Nmob;number_mob++)
            {
                cosmo->data.mass_ob[number_mob] = pow(10.,cosmo->data.log10_mass_ob[number_mob]);  //单位M_sun
                cosmo->data.ln_mass_ob[number_mob] = log(cosmo->data.mass_ob[number_mob]);
            }
        }
        else
        {
            cosmo->status = CLUCU_ERROR_MEMORY;
            clucu_cosmology_set_status_message(cosmo,
                                          "clucu_core.c: InitializeCosmo(): "
                                          "Error creating mass_ob array\n");
        }
    }

    //初始化计算量
    cosmo->data.Nbin_zob = Create1Grid(cosmo->data.Nzob);

    cosmo->data.shiftperp = Create1Grid(cosmo->data.NzobP);
    cosmo->data.shiftpara = Create1Grid(cosmo->data.NzobP);
    if(cosmo->survey.mass_ob_min != NAN)//csst
    {
        cosmo->data.Nbin_mzob = Create2Grid(cosmo->data.Nmob,cosmo->data.Nzob);
    }
    //这俩初始化为1。如果有好几个cosmo param model，就可以根据fiducial param model计算这俩。
    for(int i=0;i<cosmo->data.NzobP;i++)
    {
        cosmo->data.shiftperp[i] = 1.;
        cosmo->data.shiftpara[i]= 1.;
    }
/*
    cosmo->data.sumHMF = Create1Grid(cosmo->data.NzobP);
    cosmo->data.averageBias = Create1Grid(cosmo->data.NzobP);
    cosmo->data.Nbin_zobP = Create1Grid(cosmo->data.NzobP);

    cosmo->data.beta = Create1Grid(cosmo->data.Nzob);
    cosmo->data.V_eff = Create3Grid(cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);
    cosmo->data.V_k = Create2Grid(cosmo->data.Nk1,cosmo->data.Nk2);
    cosmo->data.bias = Create3Grid(cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);
    cosmo->data.pc = Create3Grid(cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);
    cosmo->data.pc0 = Create3Grid(cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);
    cosmo->data.lnkPc = Create3Grid(cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);

    cosmo->data.barn_mzobP=Create2Grid(cosmo->data.Nmob,cosmo->data.NzobP);
    cosmo->data.cl_hh = Create4Grid(cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.Nell);
    cosmo->data.cl_hh_hat = Create4Grid(cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.Nell);
    cosmo->data.cl_kk = Create1Grid(cosmo->data.Nell);
    cosmo->data.cl_kk_hat = Create1Grid(cosmo->data.Nell);
    cosmo->data.cl_hk1 = Create3Grid(cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.Nell);
    cosmo->data.cl_hk2 = Create3Grid(cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.Nell);
    cosmo->data.cl_hk = Create3Grid(cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.Nell);

    cosmo->data.CovNN_SV = Create3Grid(cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.Nzob);
    cosmo->data.CovNN = Create4Grid(cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.Nzob,cosmo->data.Nzob);
    cosmo->data.Covhhhh = Create7Grid(cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.NzobP,cosmo->data.Nell);
    cosmo->data.Covhkhk = Create5Grid(cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.NzobP,cosmo->data.Nell);
    cosmo->data.Covhkhh = Create6Grid(cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.NzobP,cosmo->data.Nell);
*/
}
void FreeCosmo(clucu_cosmology *cosmo)
{
    //free data
    free(cosmo->data.table_k);
    Free2Grid(cosmo->data.table_p,cosmo->data.spline_Nz,cosmo->data.table_Nk);

    //free data
    gsl_spline2d_free(cosmo->data.spline_Pkz);
    free(cosmo->data.Pkz1D);


    free(cosmo->data.sumHMF);
    free(cosmo->data.averageBias);
    free(cosmo->data.Nbin_zob);
    free(cosmo->data.Nbin_zobP);
    if(cosmo->survey.mass_ob_min != NAN)//csst 
        Free2Grid(cosmo->data.Nbin_mzob,cosmo->data.Nmob,cosmo->data.Nzob);
    free(cosmo->data.shiftperp);
    free(cosmo->data.shiftpara);
    free(cosmo->data.beta);

    Free3Grid(cosmo->data.V_eff,cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);
    Free2Grid(cosmo->data.V_k,cosmo->data.Nk1,cosmo->data.Nk2);
    Free3Grid(cosmo->data.bias,cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);
    Free3Grid(cosmo->data.pc,cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);
    Free3Grid(cosmo->data.pc0,cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);
    Free3Grid(cosmo->data.lnkPc,cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);
}




clucu_cosmology clucu_cosmology_create(clucu_configuration config,clucu_parameter param,clucu_survey survey)
{
    clucu_cosmology cosmo;
    cosmo.config = config;
    cosmo.param = param;
    cosmo.survey = survey;
    cosmo.spline_param = spline_param_default;
    cosmo.gsl_param = gsl_param_default;

    cosmo.runned_class = false;
    cosmo.computed_distances = false;
    cosmo.computed_growth = false;
    cosmo.computed_power = false;
    cosmo.computed_sigma = false;
    cosmo.computed_cluster_numbercounts = false;
    cosmo.computed_cluster_averagebias = false;
    cosmo.computed_cluster_volumeeffect = false;

    cosmo.computed_average_one_over_chis = false;
    cosmo.computed_weight_halo = false;
    cosmo.computed_weight_kappa = false;
    cosmo.computed_cl_hh = false;
    cosmo.computed_cl_kk = false;
    cosmo.computed_cl_hk1 = false;
    cosmo.computed_cl_hk2 = false;
    cosmo.computed_cl_hk = false;

  

    strcpy(cosmo.name,"cosmo_name");
    strcpy(cosmo.classname,"cosmo_classname");

    cosmo.cosmo_id=0;
    cosmo.class_id=0;
    cosmo.status=0;
    // Initialise as 0-length string
    //cosmo.status_message[0] = '\0';
//    *cosmo.status_message=NULL;
    cosmo.thread_num=1;

    cosmo.GISDB=1;
    cosmo.computed_transfer = false;
    cosmo.eh = NULL;

    return cosmo;
}


/* ------- ROUTINE: clucu_cosmology_set_status_message --------
INPUT: clucu_cosmology struct, status_string
TASK: set the status message safely.
*/
void clucu_cosmology_set_status_message(clucu_cosmology * cosmo, const char * message, ...)
{
    const int buffer = 480; /* must be < 500 - 4 */

    va_list va;
    va_start(va, message);
    // critical指定某一区域的代码，每次只能同时被一个线程执行。
    #pragma omp critical
    {
/*
        if(strlen(cosmo->status_message) != 0) 
        {
            clucu_raise_warning(CLUCU_ERROR_OVERWRITE, "Status message being overwritten:");
            fprintf(stderr, "STATUS: %d. %s\n", cosmo->status, cosmo->status_message);
        }
        vsnprintf(cosmo->status_message, buffer, message, va);
*/
        /* if bufferation happens, message[buffer - 1] is not NULL, ... will show up. */
        strcpy(&cosmo->status_message[buffer], "...");
    }
    //clucu_raise_warning(cosmo->status, cosmo->status_message);

    va_end(va);
}

void SaveNmzob(clucu_cosmology cosmo)
{
	char namedir[50];
    char nameNlzob[100];
	sprintf(namedir,"../output/%s", cosmo.name);
	sprintf(nameNlzob,"../output/%s/Nlzob%d.dat", cosmo.name,cosmo.cosmo_id);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameNlzob,"w");
    save_matrix(fp,cosmo.data.Nbin_mzob,cosmo.data.log10_mass_ob,cosmo.data.zob,cosmo.data.Nmob,cosmo.data.Nzob);
	fclose(fp);
}