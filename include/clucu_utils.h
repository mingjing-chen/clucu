#pragma once
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdio.h>
#include <sys/stat.h>
#include <sys/types.h>
#include <omp.h>
#include <stdbool.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_roots.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_spline2d.h>
#include <gsl/gsl_integration.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_deriv.h>
#include <gsl/gsl_linalg.h>



#include "clucu_error.h"
#define kronecker(i,j) (i==j?1:0)
#define  max(d,f)  (d>f?d:f)
#define  min(d,f)  (d<f?d:f)

double * clucu_linear_spacing(double xmin, double xmax, int N);
double * clucu_loglin_spacing(double xminlog, double xmin, double xmax, int Nlog, int Nlin);
double * clucu_linlog_spacing(double xmin, double xminlog, double xmax, int Nlin, int Nlog);
double * clucu_log_spacing(double xmin, double xmax, int N);
double * clucu_log10_spacing(double xmin, double xmax, int N);

double GaussianDistribution(double x,double sigma);
double Ci(double x);
double Si(double x);

double *Create1Grid(int m);
double **Create2Grid(int m, int n);
double ***Create3Grid(int m, int n, int t);
double ****Create4Grid(int m, int n, int t, int q);
double *****Create5Grid(int m, int n, int t, int q, int x);
double ******Create6Grid(int m, int n, int t, int q, int x, int y);
double *******Create7Grid(int m, int n, int t, int q, int x, int y, int z);
double ********Create8Grid(int m, int n, int t, int q, int x, int y, int z, int v);

void Free1Grid(double *tt, int m);
void Free2Grid(double **tt, int m, int n);
void Free3Grid(double ***tt, int m, int n, int t);
void Free4Grid(double ****tt, int m, int n, int t, int p);
void Free5Grid(double *****tt, int m, int n, int t, int p, int x);
void Free6Grid(double ******tt, int m, int n, int t, int p, int x, int y);
void Free7Grid(double *******tt, int m, int n, int t, int p, int x, int y, int z);
void Free8Grid(double ********tt, int m, int n, int t, int p, int x, int y, int z, int v);

void Copy1Grid(double *from, double *to, int m);
void Copy2Grid(double **from, double **to, int m, int n);
void Copy3Grid(double ***from, double ***to, int m, int n, int t);
void Copy4Grid(double ****from, double ****to, int m, int n, int t, int q);
void PrintMatrix(double **m,int n1,int n2);
void InverseMatrix(double *m, int n);
void save_matrix(FILE *fp,double **matrix,double *x,double *y,int nx,int ny);

typedef struct clucu_physical_constants {
    /**
     * Lightspeed / H0 in units of Mpc/h (from CODATA 2014)
     */
    double CLIGHT_HMPC;

    /**
     * Newton's gravitational constant in units of m^3/Kg/s^2
     */
    double GNEWT;

    /**
     * Solar mass in units of kg (from GSL)
     */
    double SOLAR_MASS;

    /**
     * Mpc to meters (from PDG 2013)
     */
    double MPC_TO_METER;

    /**
     * pc to meters (from PDG 2013)
     */
    double PC_TO_METER;

    /**
     * Rho critical in units of M_sun/h / (Mpc/h)^3
     */
    double RHO_CRITICAL;

    /**
     * Boltzmann constant in units of J/K
     */
    double KBOLTZ;

    /**
     * Stefan-Boltzmann constant in units of kg/s^3 / K^4
     */
    double STBOLTZ;

    /**
     * Planck's constant in units kg m^2 / s
     */
    double HPLANCK;

    /**
     * The speed of light in m/s
     */
    double CLIGHT;

    /**
     * Electron volt to Joules convestion
     */
    double EV_IN_J;

    /**
     * Temperature of the CMB in K
     */
    double T_CMB;

    /**
     * T_ncdm, as taken from CLASS, explanatory.ini
     */
    double TNCDM;

    /**
     * Kelvin to ev
     */
    double K_TO_EV;
    //8.617e-5

    /**
     * neutrino mass splitting differences
     * See Lesgourgues and Pastor, 2012 for these values.
     * Adv. High Energy Phys. 2012 (2012) 608515,
     * arXiv:1212.6154, page 13
     */
    double DELTAM12_sq;
    double DELTAM13_sq_pos;
    double DELTAM13_sq_neg;
} clucu_physical_constants;
extern const clucu_physical_constants clucu_constants;
//extern表明变量或者函数是定义在其他其他文件中的

//f2d extrapolation types for early times
typedef enum clucu_f2d_extrap_growth_t
{
  clucu_f2d_clucugrowth = 401, //Use CLUCU's linear growth
  clucu_f2d_constantgrowth = 403, //Use a constant growth factor
  clucu_f2d_no_extrapol = 404, //Do not extrapolate, just throw an exception
} clucu_f2d_extrap_growth_t;
//f2d interpolation types 插值方法
typedef enum clucu_f2d_interp_t
{
  clucu_f2d_3 = 303, //Bicubic interpolation
} clucu_f2d_interp_t;
//定义一个二维插值结构体，带边界条件。计算功率谱
typedef struct {
  double lkmin,lkmax; /**< Edges in log(k)*/
  double zmin,zmax; /**< Edges in a*/
  int is_factorizable; /**< Is this factorizable into k- and a-dependent functions? */
  int is_k_constant; /**< no k-dependence, just return 1*/
  int is_z_constant; /**< no a-dependence, just return 1*/
  int extrap_order_lok; /** 小k怎么外推 < Order of extrapolating polynomial in log(k) for low k (0, 1 or 2)*/ 
  int extrap_order_hik; /** 大k怎么外推 < Order of extrapolating polynomial in log(k) for high k (0, 1 or 2)*/
  clucu_f2d_extrap_growth_t extrap_linear_growth;  /** 大z怎么外推 < Extrapolation type at high redshifts*/
  int is_log; /**< Do I hold the values of log(f(k,a))?*/
  double growth_factor_0; /**< Constant extrapolating growth factor*/
  int growth_exponent; /**< Power to which growth should be exponentiated*/
  gsl_spline *fk; /**< Spline holding the values of the k-dependent factor*/
  gsl_spline *fz; /**< Spline holding the values of the a-dependent factor*/
  gsl_spline2d *fkz; /**< Spline holding the values of f(k,a)*/
} clucu_f2d_t;
//k外推方法。0（边界不变），1（边界1阶导不变），2（边界2阶导不变）
//z外推方法，clucu_f2d_clucugrowth（使用公式计算D(z)），else（边界固定为growth_factor_0 ** growth_exponent）
//is_log表示，你输入的数组是不是log。 clucu_f2d_t_eval函数输出的肯定不是log


gsl_matrix *gsl_matrix_dot(gsl_matrix *a,gsl_matrix *b);
gsl_matrix *gsl_matrix_inv(gsl_matrix *a);
void clucu_matrix_inv(int n,double **mat,double**inv_mat);
double gsl_matrix_trace(gsl_matrix *a);

double clucu_integ_spline(int nx,double *x,double *y,
                      double a, double b,
                      const gsl_interp_type *T);