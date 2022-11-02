#pragma once
#include "clucu_error.h"
#include "clucu_utils.h"
#include "clucu_config.h"
#include "clucu_neutrino.h"

#include <gsl/gsl_interp.h>

#define THREAD_RATE 0.5
#define CHECK 0
/**
 * Struct containing the parameters defining a cosmology
 */
typedef struct clucu_parameter {
    //--------hubble-------------
    double h;
	double H0;
    //--------photon-------------
    double T_CMB;
    double Omega_g;
    double Ogh2;
    //--------neutrino------------- 
    //输入计算方式，mass表示输入单个质量，
    //如果已经分配好了，输入splited
    //如果输入sum，无论是m_nu_sum,Omega_nu_mass_sum,Onuh2_mass_sum，需要指定分配方式，single,normal,inverted,equal
    clucu_neutrino_mass_split_label neutrino_mass_split_label;
    //默认输入m_nu[3] 
    //或者输入Onuh2_mass_sum，以及分配方式
    double m_nu[3];//0.06eV 
    double Omega_nu_mass[3]; // Omega_nu for MASSIVE neutrinos
    double Onuh2_mass[3]; 
    double m_nu_sum; // sum of the neutrino masses.
    double Omega_nu_mass_sum;
    double Onuh2_mass_sum;
    // Effective number of relativistic neutrino species in the early universe.
    double Neff; //3.046
    // Number of species of neutrinos which are nonrelativistic today
    int N_nu_mass; //1个0.06eV，根据m_nu计算
    // Number of species of neutrinos which are relativistic  today
    double N_nu_rel;//3.046 - 1*TNCDM^4/ (4./11.)^(4./3.)=2.0328，根据N_nu_mass计算
    // total mass of massive neutrinos (This is a pointer so that it can hold multiple masses.)
    double Omega_nu_rel; // Omega_nu for MASSLESS neutrinos
    double Onuh2_rel;
    //--------cdm-------------
    double Omega_c;
    double Och2;
    //--------baryon-------------
    double Omega_b;
    double Obh2; 
    //--------mater-------------    
    double Omh2;
    double Omega_m;
    //--------DE-------------
    double Omega_l;
	double w0;
	double wa;
    //--------geometry-------------
    double Omega_k;
    int k_sign;
    double sqrtk;
    //--------Pk-------------
    double n_s;
    double alpha_s;
	double sigma8;//sigma8和delta_zeta有一个
    double delta_zeta;
    double k_pivot;//Mpc-1
    double A_s;//=delta_zeta^2
    
    //--------其他-------------
    double tau_reio;//用在CMB
    double f_nl;//non-gaussianity

    
} clucu_parameter;

/**
 * @brief 
 * 
 */
typedef struct clucu_survey {
  // These are all functions of the redshift z.
  /*--------需要更改的望远镜参数（3个）-------------*/
	double log10ASZ;/*!< Detailed description of the member var1 */
	double betaSZ;
	double gammaSZ;
	double flux;
    /*----------concentration 参数（3个）---------------*/
    double A_vir;
    double B_vir;
    double C_vir;
    /*----------halo off centering 参数（3个）---------------*/
    double f_cen0;
    double p_cenM;
    double p_cenz;
    double sigma_s0;
    double p_sigmaM;
    double p_sigmaz;
    /*----------mass-ob 参数（14个）---------------*/
    double lnM_b0;
    double q_b[3];
    double s_b[3];
    double sigma_lnM0;
    double q_sigma_lnm[3];
    double s_sigma_lnm[3];
    /*----------source redshift 参数（1个）---------------*/
    double source_redshift;
	/*----------cluster mass ob参数（8个）---------------*/
	double mass_ob_A;
    double mass_ob_B; 
    double mass_ob_Bz; 
    double mass_ob_Cz; 
    double mass_ob_sigma0; 
    double mass_ob_q;
    double mass_ob_qz;
    double mass_ob_pz;
    /*-------------初始参数设置-------------*/
    double redshift_sigma0;
    double survey_area;
    double Delta_Omega;
    double mass_ob_min ;//输入mass_ob的最小值. 如果该值为NAN，则用mass_limit函数
    double mass_ob_maxc ;
    double mass_ob_max ;//输入mass_ob的最大值
    double mass_ob_min_h ;//输入mass_ob的最小值. 如果该值为NAN，则用mass_limit函数
    double mass_ob_maxc_h ;
    double mass_ob_max_h ;//输入mass_ob的最大值
    double dln_mass_ob ;//输入mass_ob bin。
    double dlog10_mass_ob ;//输入mass_ob bin。
//运行class的z, sigma插值的z
    double spline_z_minlog ;//输入红移的最小log
    double spline_z_max ;//输入红移的最大值
    double spline_z_nlin ;//输入红移个数，lin
    double spline_z_nlog ;//输入红移个数, log

    double zob_min ;//输入实际观测红移的最小值
    double zob_max ;//输入实际观测红移的最大值
    double zob_bin ;//输入实际观测红移bin【注意，星系团计数是0.05，星系团功率谱是0.2】
	double zobP_min ;//输入实际观测红移的最小值
    double zobP_max ;//输入实际观测红移的最大值
    double zobP_bin ;//输入实际观测红移bin【注意，星系团计数是0.05，星系团功率谱是0.2】
    double k1_min ;//计算星系团功率谱时的k垂直，单位1/Mpc
    double k1_max ;//单位1/Mpc
    double k1_bin ;//单位1/Mpc
    double k2_min ;//计算星系团功率谱时的k平行，单位1/Mpc
    double k2_max ;//单位1/Mpc
    double k2_bin ;//单位1/Mpc
 
} clucu_survey;

typedef struct clucu_data {
    int spline_Nz;
    int Nzob;
    int NzobP;
    int Nk1;
    int Nk2;
    double *spline_z;
    double *zob;
	double *zobP;
    double *k1;
    double *k2;
    
    int Nmob;
    double *ln_mass_ob;
    double *log10_mass_ob;
    double *mass_ob;
    // Distances are defined in Mpc
    gsl_spline *chi;//clucu_cosmology_compute_distances
    gsl_spline *E;
    gsl_spline *zchi;
    double growth0;//clucu_cosmology_compute_growth
    gsl_spline *growth;
    gsl_spline *fgrowth;
    //pkz
    clucu_f2d_t *log_pkz;//lnk,ln pkz
    clucu_f2d_t *transfer_g;
    clucu_f2d_t *transfer_ur;
    clucu_f2d_t *transfer_cb;
    clucu_f2d_t *transfer_nu1;
    clucu_f2d_t *transfer_nu2;
    clucu_f2d_t *transfer_nu3;
    /*----------CLASS的输出数组---------------*/
	//double *table_z,double *table_H, double *table_dC,
	//double *table_dA,double *table_dL,double *table_D,double *table_f,int table_Nz
	int table_Nk;
	double *table_k;
	double **table_p;
	double int_k_min;
	double int_k_max;
    /*----------由输出数组得到插值表---------------*/
	gsl_spline2d *spline_Pkz;
	double *Pkz1D;
    /*----------一些计算量---------------*/
    gsl_spline2d *spline_massconvert;
    gsl_spline2d *spline_sigma;//取一组log10 m，a,算一组log sigma，并插值
    gsl_spline2d *spline_kernel_n;//取一组ln m，z,算一组kernel_n，并插值
    gsl_spline2d *spline_kernel_nb;//取一组ln m，z,算一组kernel_nb，并插值
    gsl_spline *spline_sigmaS3;//sigmaS3 - log10M
    gsl_spline *spline_delta_crit;//sigmaS3 - log10M
    gsl_spline2d *spline_SoverL;
    /*-----------星系团数变量---------------*/
    double *sumHMF;
    double *averageBias;
    gsl_spline2d *spline_averageBias_k;//如果有k，则对其插值，(k,z)
	double *Nbin_zob;//不加mass-richiness 的星系团数，但有Mlimit限制
    double *Nbin_zobP;//计算bias时一起计算，作为参考
	double **Nbin_mzob;//加mass-richiness 的星系团数
    /*-----------星系团功率谱变量---------------*/
    double *shiftperp;
    double *shiftpara;
	double *beta;
	double ***V_eff;
	double **V_k;
	double ***bias;
	double ***pc;
	double ***pc0;
    double ***lnkPc;
    /*-----------角功率谱变量---------------*/
    int Nell;
    double *ell;
    double **barn_mzobP;
    double ****cl_hh;
    double *cl_kk;
    double ***cl_hk2;
    double ***cl_hk1;
    double ***cl_hk;
    void **spline_Wh;
    gsl_spline *spline_bz;
    double average_one_over_chis;
    gsl_spline *spline_wkappa_z;
    /*-----------covariance---------------*/
    double ***CovNN_SV;
    double ****CovNN;
    double ****cl_hh_hat;
    double *cl_kk_hat;
    double *****Covhkhk;
    double *******Covhhhh;
    double ******Covhkhh;
    /*-----------cmb---------------*/
    double *CTT;
    double *CTT_hat;
    double *CEE;
    double *CEE_hat;
    double *CTE;
    double *CovTT;
    double *CovEE;
    double *CovTE;
}clucu_data;

/**
 * Struct that contains all the parameters needed to create certain splines.
 * This includes splines for the scale factor, masses, and power spectra.
 */
typedef struct clucu_spline_param {
    double Z_SPLINE_MIN;
    //功率谱插值
    double Z_SPLINE_MINLOG_PK;
    double Z_SPLINE_MAX_PK;
    int Z_SPLINE_NLIN_PK;
    int Z_SPLINE_NLOG_PK;
    //sigma插值
    double Z_SPLINE_MINLOG_SM;//对sigma(m,z)插值时，z的最小log值（log插值）
    double Z_SPLINE_MAX_SM;
    int Z_SPLINE_NLIN_SM;
    int Z_SPLINE_NLOG_SM;
    //对chi(z),D(z)插值，背景插值
    double Z_SPLINE_MINLOG_BG;
    double Z_SPLINE_MAX_BG;
    int Z_SPLINE_NLIN_BG;
    int Z_SPLINE_NLOG_BG;
  
    // mass splines
    double LOG10M_SPLINE_MIN;//对sigma(m,z)插值时，log10 m的最小值
    double LOG10M_SPLINE_MAX;//对sigma(m,z)插值时，log10 m的最大值
    int LOG10M_SPLINE_NM;//对sigma(m,z)插值时，log10 m的个数

  //k-splines and integrals
  double K_MAX_SPLINE;
  double K_MAX;//功率谱积分的最大值
  double K_MIN;//功率谱积分的最小值
  double DLOGK_INTEGRATION;
  double DCHI_INTEGRATION;
  int N_K;
  int N_K_3DCOR;

  //Correlation function parameters
  double ELL_MIN_CORR;
  double ELL_MAX_CORR;
  int N_ELL_CORR;

  // interpolation types
  const gsl_interp_type **Z_SPLINE_TYPE;
  const gsl_interp_type **K_SPLINE_TYPE;
  const gsl_interp_type **M_SPLINE_TYPE;
  const gsl_interp_type **D_SPLINE_TYPE;
  const gsl_interp2d_type **PNL_SPLINE_TYPE;
  const gsl_interp2d_type **PLIN_SPLINE_TYPE;
  const gsl_interp_type **CORR_SPLINE_TYPE;
} clucu_spline_param;


/**
 * Struct that contains parameters that control the accuracy of various GSL
 * routines.
 */
typedef struct clucu_gsl_param {
  // General parameters
  //积分workspace的N
  size_t N_ITERATION;

  // Integration
  int INTEGRATION_GAUSS_KRONROD_POINTS;
  double INTEGRATION_EPSREL;
  // Limber integration
  int INTEGRATION_LIMBER_GAUSS_KRONROD_POINTS;
  double INTEGRATION_LIMBER_EPSREL;
  // Distance integrals
  double INTEGRATION_DISTANCE_EPSREL;
  // sigma_R integral
  double INTEGRATION_SIGMAR_EPSREL;//sigma(r)积分的相对误差
  // k_NL integral
  double INTEGRATION_KNL_EPSREL;

  // Root finding
  double ROOT_EPSREL;
  int ROOT_N_ITERATION;

  // ODE
  double ODE_GROWTH_EPSREL;

  // growth
  double EPS_SCALEFAC_GROWTH;

  // halo model
  double HM_MMIN;
  double HM_MMAX;
  double HM_EPSABS;
  double HM_EPSREL;
  size_t HM_LIMIT;
  int HM_INT_METHOD;

  // Flags for using spline integration
  bool NZ_NORM_SPLINE_INTEGRATION;
  bool LENSING_KERNEL_SPLINE_INTEGRATION;

  // numbercounts integral
  double INTEGRATION_MASS_TRUE_EPSREL;//dlnm_true积分的相对误差
  double INTEGRATION_Z_TRUE_EPSREL;//dz_true积分的相对误差
  double INTEGRATION_MASS_OB_EPSREL;//dlnm_true积分的相对误差
  double INTEGRATION_Z_OB_EPSREL;//dz_true积分的相对误差

} clucu_gsl_param;



/**
 * Sturct containing references to instances of the above structs, and boolean flags of precomputed values.
 */
typedef struct clucu_cosmology {
    clucu_configuration         config;
    clucu_parameter            param;
    clucu_survey            survey;
    clucu_spline_param spline_param;
    clucu_gsl_param    gsl_param;

    clucu_data             data;
    
    //不一定需要的
    bool computed_massconvert;
    
    bool runned_class;
    bool computed_distances;
    bool computed_growth;
    bool computed_power;
    bool computed_sigmaS3;

    bool computed_sigma;
        bool computed_cluster_numbercounts;
        bool computed_cluster_averagebias;
        bool computed_cluster_clusterpowerspectrum;
        bool computed_cluster_volumeeffect;
    bool computed_average_one_over_chis;
    bool computed_weight_halo;
    bool computed_weight_kappa;
    bool computed_cl_hh;
    bool computed_cl_kk;
    bool computed_cl_hk1;
    bool computed_cl_hk2;
    bool computed_cl_hk;

    bool computed_kernel_n;
    bool computed_kernel_nb;
    bool computed_cluster_barn;

    bool computed_CovNN_SV;
    bool computed_CovNN;
    bool computed_Covhhhh;
    bool computed_Covhkhk;
    bool computed_Covhkhh;


    int status;
    //this is optional - less tedious than tracking all numerical values for status in error handler function
    char status_message[500];

    int GISDB;//1表示考虑了GISDB，0表示没有考虑。

    // other flags?
    int cosmo_id;//0是fiducial，1&2是改变第一个参数。。
    int class_id;//有的模型（比如改变mass-richiness参数）不用重新跑class,这部分复制cosmo0即可
    int thread_num;
    char name[50];
    char classname[50];
    
    int other_flag_1;
    int other_flag_2;
    bool computed_transfer;
    void *eh;

} clucu_cosmology;




extern const clucu_gsl_param gsl_param_default;
extern const clucu_spline_param spline_param_default;

extern const clucu_parameter param_default;
extern const clucu_parameter param_planck2018;
extern const clucu_parameter param_wmap9;
extern const clucu_parameter param_WMAP5;
extern const clucu_parameter param_oguri;
extern const clucu_survey survey_sze;
extern const clucu_survey survey_xrayDeep;
extern const clucu_survey survey_xrayWide;
extern const clucu_survey survey_csst;
extern const clucu_survey survey_csst2;
extern const clucu_survey survey_oguri;
//初始化cosmo，1.初始化宇宙学模型，2.并给数组动态分配空间，3.线性分配z,lambda,k1,k2
void InitializeCosmo(clucu_cosmology *cosmo);
//释放cosmo
void FreeCosmo(clucu_cosmology *cosmo);

clucu_cosmology clucu_cosmology_create(clucu_configuration config,clucu_parameter param,clucu_survey survey);

void clucu_cosmology_set_status_message(clucu_cosmology * cosmo, const char * message, ...);

void SaveNmzob(clucu_cosmology cosmo);
