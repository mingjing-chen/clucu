#pragma once

/**
 * Matter power spectrum typedef.
 * Contains all information that describes a specific
 * matter power spectrum. This inclues whether we
 * want the linear power spectrum, whether we use
 * halofit, and what halo model is being used.
 */
typedef enum matter_power_spectrum_type
{
    clucu_bbks = 1,
    clucu_eisenstein_hu = 2,
    clucu_boltzmann_class_Pk = 3,
    clucu_boltzmann_class_Tk = 4,
    
    
}matter_power_spectrum_type;

typedef enum mass_function_type
{
    clucu_massfunc_press1974 =1,
    clucu_massfunc_sheth1999 =2,
    clucu_massfunc_jenkins2001 =3,
    clucu_massfunc_tinker2008 =4,
    clucu_massfunc_bhattacharya2011 =5,
    clucu_massfunc_bocquet2016 =6,
} mass_function_type;

typedef enum bias_function_type
{
    clucu_bias_press1974 =1,
    clucu_bias_sheth1999 =2,
    clucu_bias_sheth2001 =3,
    clucu_bias_tinker2010 =4,
    clucu_bias_bhattacharya2011 =5,
} bias_function_type;

typedef enum halo_concentration_type
{
  clucu_concentration_duffy2008 = 1,
  clucu_concentration_bhattacharya2013 = 2,
  clucu_concentration_diemer2015 = 3,
} halo_concentration_type;


typedef enum mass_observable_type
{
  clucu_murrata2019 = 1,
  clucu_my_csst = 2,
  clucu_oguri2011 = 3,
} mass_observable_type;

typedef enum halo_define_type
{
  clucu_delta_vir = 1,
  clucu_delta_200c = 2,
  clucu_delta_200m = 3,
  clucu_delta_500m = 4,
  clucu_delta_500c =5,
  clucu_delta_180m =6,
} halo_define_type;
//这个没包括进大结构体中，所以给他单独加了个前缀clucu
/**
 * Configuration typedef.
 * This contains the transfer function,
 * matter power spectrum, and mass function
 * that is being used currently.
 */
typedef struct clucu_configuration {
  matter_power_spectrum_type  matter_power_spectrum_method;
  mass_function_type          mass_function_method;
  bias_function_type          bias_function_method;
  halo_concentration_type     halo_concentration_method;
  mass_observable_type        mass_observable_method;
  halo_define_type            halo_define_method;
  int GISDB;
} clucu_configuration;


// The default configuration object
extern const clucu_configuration config_default;
extern const clucu_configuration config_oguri;

/*
int cosmo_id;//0是fiducial，1&2是改变第一个参数。。
	int class_id;//有的模型（比如改变mass-richiness参数）不用重新跑class,这部分复制cosmo0即可
	int runClass;//1跑，0不跑
	int thread_num;
	char name[50];
	char classname[50];
	double *delta_p;
	int Np;
	int fft;//选择拟合函数，1:PS，2:Tinker2008，3:Jenkins2001, 4:SMT
	int flag;
    int biasM_flag;//1:ST, 2:SMT, 3:Seljak
	int GISDB;//1表示考虑了GISDB，0表示没有考虑。
    */