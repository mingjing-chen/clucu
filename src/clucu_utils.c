#include "clucu_utils.h"
#include <math.h>

const clucu_physical_constants clucu_constants = {
    /**
     * Lightspeed / H0 in units of Mpc/h (from CODATA 2014)
     */
    .CLIGHT_HMPC=2997.92458,

    /**
     * Newton's gravitational constant in units of m^3/Kg/s^2
     */
    //6.6738e-11,  /(from PDG 2013) in m^3/Kg/s^2
    //6.67428e-11, // CLASS VALUE
    .GNEWT=6.67408e-11, // from CODATA 2014

    /**
     * Solar mass in units of kg (from GSL)
     */
    //GSL_CONST_MKSA_SOLAR_MASS,
    //1.9885e30, //(from PDG 2015) in Kg
    .SOLAR_MASS=1.9884754153381438E+30, //from IAU 2015

    /**
     * Mpc to meters (from PDG 2016 and using M_PI)
     */
    .MPC_TO_METER=3.085677581491367399198952281E+22,

    /**
     * pc to meters (from PDG 2016 and using M_PI)
     */
    .PC_TO_METER=3.085677581491367399198952281E+16,

    /**
     * Rho critical in units of M_sun/h / (Mpc/h)^3
     */
    .RHO_CRITICAL=((3*100*100)/(8*M_PI*6.67408e-11)) * (1000*1000*3.085677581491367399198952281E+22/1.9884754153381438E+30),

    /**
     * Boltzmann constant in units of J/K
     */
    //GSL_CONST_MKSA_BOLTZMANN,
    .KBOLTZ=1.38064852E-23, //from CODATA 2014

    /**
     * Stefan-Boltzmann constant in units of kg/s^3 / K^4
     */
    //GSL_CONST_MKSA_STEFAN_BOLTZMANN_CONSTANT,
    .STBOLTZ=5.670367E-8, //from CODATA 2014

    /**
     * Planck's constant in units kg m^2 / s
     */
    //GSL_CONST_MKSA_PLANCKS_CONSTANT_H,
    .HPLANCK=6.626070040E-34, //from CODATA 2014
    //h_bar = h /2pi

    /**
     * The speed of light in m/s
     */
    //GSL_CONST_MKSA_SPEED_OF_LIGHT,
    .CLIGHT=299792458.0, //from CODATA 2014

    /**
     * Electron volt to Joules convestion
     */
    //GSL_CONST_MKSA_ELECTRON_VOLT,
    .EV_IN_J=1.6021766208e-19,  //from CODATA 2014

    /**
     * Temperature of the CMB in K
     */
    .T_CMB=2.7255,
    //2.7255, // CLASS value

    /**
     * T_ncdm, as taken from CLASS, explanatory.ini
     */
    .TNCDM=0.71611,

    /**
     * Kelvin to ev
     */
    .K_TO_EV=8.617e-5,

    /**
     * neutrino mass splitting differences
     * See Lesgourgues and Pastor, 2012 for these values.
     * Adv. High Energy Phys. 2012 (2012) 608515,
     * arXiv:1212.6154, page 13
     */
    .DELTAM12_sq=7.62E-5,
    .DELTAM13_sq_pos=2.55E-3,
    .DELTAM13_sq_neg=-2.43E-3
};

/*----------------------*/
/*      动态分布内存      */
/*----------------------*/

/* ------- ROUTINE: clucu_linear spacing ------
INPUTS: [xmin,xmax] of the interval to be divided in N points, N-1 interval
OUTPUT: bin edges in range [xmin,xmax]
包括端点，x=[0,1,2,3],dx=1,N=4
*/
double * clucu_linear_spacing(double xmin, double xmax, int N)
{
    if (N<2) 
    {
        CLUCU_RAISE_WARNING(CLUCU_ERROR_LINSPACE,"Cannot make log-spaced array with %d points - need at least 2\n", N);
    }
    double dx = (xmax-xmin)/(N -1.);

    double *x = malloc(sizeof(double)*N);
    if (x==NULL) 
    {
        CLUCU_RAISE_WARNING(
            CLUCU_ERROR_MEMORY,
            "Could not allocate memory for linear-spaced array (N=%d)\n", N);
        return x;
    }

    for (int i=0; i<N; i++)
        x[i] = xmin + dx*i;

    x[0]=xmin; //Make sure roundoff errors don't spoil edges
    x[N-1]=xmax; //Make sure roundoff errors don't spoil edges

    return x;
}
/* ------- ROUTINE: clucu_log spacing ------
INPUTS: [xmin,xmax] of the interval to be divided logarithmically in N bins
TASK: divide an interval in N logarithmic bins
OUTPUT: bin edges in range [xmin,xmax]
*/
double * clucu_log_spacing(double xmin, double xmax, int N)
{
    if (N<2)
    {
        CLUCU_RAISE_WARNING(
            CLUCU_ERROR_LOGSPACE,
            "Cannot make log-spaced array with %d points - need at least 2\n", N);
        return NULL;
    }

    if (!(xmin>0 && xmax>0)) 
    {
        CLUCU_RAISE_WARNING(
            CLUCU_ERROR_LOGSPACE,
            "Cannot make log-spaced array xmax or xmax non-positive (had %le, %le)\n", xmin, xmax);
        return NULL;
    }

    double log_xmax = log(xmax);
    double log_xmin = log(xmin);
    double dlog_x = (log_xmax - log_xmin) /  (N-1.);

    double *x = malloc(sizeof(double)*N);
    if (x==NULL) 
    {
        CLUCU_RAISE_WARNING(
            CLUCU_ERROR_MEMORY,
            "Could not allocate memory for log-spaced array (N=%d)\n", N);
        return x;
    }

    double xratio = exp(dlog_x);
    x[0] = xmin; //Make sure roundoff errors don't spoil edges
    for (int i=1; i<N-1; i++)
        x[i] = x[i-1] * xratio;
    x[N-1]=xmax; //Make sure roundoff errors don't spoil edges

    return x;
}
double * clucu_log10_spacing(double xmin, double xmax, int N)
{
    if (N<2) 
    {
        CLUCU_RAISE_WARNING(
            CLUCU_ERROR_LOG10SPACE,
            "Cannot make log10-spaced array with %d points - need at least 2\n", N);
        return NULL;
    }

    if (!(xmin>0 && xmax>0)) 
    {
        CLUCU_RAISE_WARNING(
            CLUCU_ERROR_LOG10SPACE,
            "Cannot make log10-spaced array xmax or xmax non-positive (had %le, %le)\n", xmin, xmax);
        return NULL;
    }

    double log_xmax = log10(xmax);
    double log_xmin = log10(xmin);
    double dlog_x = (log_xmax - log_xmin) /  (N-1.);

    double *x = malloc(sizeof(double)*N);
    if (x==NULL) 
    {
        CLUCU_RAISE_WARNING(
            CLUCU_ERROR_MEMORY,
            "Could not allocate memory for log-spaced array (N=%d)\n", N);
        return x;
    }

    double xratio = pow(10.,dlog_x);
    x[0] = xmin; //Make sure roundoff errors don't spoil edges
    for (int i=1; i<N-1; i++)
        x[i] = x[i-1] * xratio;
    x[N-1]=xmax; //Make sure roundoff errors don't spoil edges

    return x;
}
/* ------- ROUTINE: clucu_linlog spacing ------
 * INPUTS: [xminlog,xmax] of the interval to be divided in bins
 *         xmin when linear spacing starts
 *         Nlog number of logarithmically spaced bins
 *         Nlin number of linearly spaced bins
 * OUTPUT: bin edges in range [xminlog,xmax]
 * [xminlog,xmin,xmax]，前面按log分，后面按lin分
 * */
double * clucu_loglin_spacing(double xminlog, double xmin, double xmax, int Nlog, int Nlin)
{
    if (Nlog<2) 
    {
        CLUCU_RAISE_WARNING(
            CLUCU_ERROR_LINLOGSPACE,
            "Cannot make log-spaced array with %d points - need at least 2\n", Nlog);
        return NULL;
    }

    if (!(xminlog>0 && xmin>0)) 
    {
        CLUCU_RAISE_WARNING(
            CLUCU_ERROR_LINLOGSPACE,
            "Cannot make log-spaced array xminlog or xmin non-positive (had %le, %le)\n", xminlog, xmin);
        return NULL;
    }

    if (xminlog>xmin)
    {
        CLUCU_RAISE_WARNING(CLUCU_ERROR_LINLOGSPACE, "ERROR: xminlog must be smaller as xmin");
        return NULL;
    }

    if (xmin>xmax)
    {
        CLUCU_RAISE_WARNING(CLUCU_ERROR_LINLOGSPACE, "ERROR: xmin must be smaller as xmax");
        return NULL;
    }

    double *x = malloc(sizeof(double)*(Nlin+Nlog-1));
    if (x==NULL) 
    {
        CLUCU_RAISE_WARNING(
            CLUCU_ERROR_MEMORY,
            "Could not allocate memory for array of size (Nlin+Nlog-1)=%d)\n", (Nlin+Nlog-1));
        return x;
    }

    double dx = (xmax-xmin)/(Nlin -1.);
    double log_xchange = log(xmin);
    double log_xmin = log(xminlog);
    double dlog_x = (log_xchange - log_xmin) /  (Nlog-1.);

    for (int i=0; i<Nlin+Nlog-1; i++) 
    {
        if (i<Nlog)
            x[i] = exp(log_xmin + dlog_x*i);
        if (i>=Nlog)
            x[i] = xmin + dx*(i-Nlog+1);
    }

    x[0]=xminlog; //Make sure roundoff errors don't spoil edges
    x[Nlog-1]=xmin; //Make sure roundoff errors don't spoil edges
    x[Nlin+Nlog-2]=xmax; //Make sure roundoff errors don't spoil edges

    return x;
}
/* ------- ROUTINE: clucu_linlog spacing ------
 * INPUTS: [xminlog,xmax] of the interval to be divided in bins
 *         xmin when linear spacing starts
 *         Nlog number of logarithmically spaced bins
 *         Nlin number of linearly spaced bins
 * OUTPUT: bin edges in range [xminlog,xmax]
 * [xmin,xminlog,xmax]，前面按lin分，后面按log分
 * */
double * clucu_linlog_spacing(double xmin, double xminlog, double xmax, int Nlin, int Nlog)
{
    if (Nlin<2) 
    {
        CLUCU_RAISE_WARNING(
            CLUCU_ERROR_LINLOGSPACE,
            "Cannot make lin-spaced array with %d points - need at least 2\n", Nlin);
        return NULL;
    }

    if (xminlog<=0) 
    {
        CLUCU_RAISE_WARNING(
            CLUCU_ERROR_LINLOGSPACE,
            "Cannot make lin-spaced array xminlog non-positive (had %le)\n", xminlog);
        return NULL;
    }

    if (xmin>xminlog)
    {
        CLUCU_RAISE_WARNING(CLUCU_ERROR_LINLOGSPACE, "ERROR: xmin must be smaller as xminlog");
        return NULL;
    }

    if (xminlog>xmax)
    {
        CLUCU_RAISE_WARNING(CLUCU_ERROR_LINLOGSPACE, "ERROR: xminlog must be smaller as xmax");
        return NULL;
    }

    double *x = malloc(sizeof(double)*(Nlin+Nlog-1));
    if (x==NULL) 
    {
        CLUCU_RAISE_WARNING(
            CLUCU_ERROR_MEMORY,
            "Could not allocate memory for array of size (Nlin+Nlog-1)=%d)\n", (Nlin+Nlog-1));
        return x;
    }

    double dx = (xminlog-xmin)/(Nlin -1.);
    double log_xmax = log(xmax);
    double log_xmin = log(xminlog);
    double dlog_x = (log_xmax - log_xmin) /  (Nlog-1.);

    for (int i=0; i<Nlin+Nlog-1; i++) 
    {
        if (i<Nlin)
            x[i] = xmin + dx*i;
        if (i>=Nlin)
            x[i] = exp(log_xmin + dlog_x*(i-Nlin+1));
    }

    x[0]=xmin; //Make sure roundoff errors don't spoil edges
    x[Nlin-1]=xminlog; //Make sure roundoff errors don't spoil edges
    x[Nlin+Nlog-2]=xmax; //Make sure roundoff errors don't spoil edges

    return x;
}


/* ------- ROUTINE: GaussianDistribution ------
 * INPUTS: x,sigma
 * TASK: gaussian distribution
 * OUTPUT: p(x,sigma)
 * */
double GaussianDistribution(double x,double sigma)
{
    if(fabs(sigma)<=1.0e-15)
		if(fabs(x)<=1.0e-15)
			return 1.0;
		else 
			return 0.0;
	else 
		return  exp(-0.5*pow(x/sigma,2)) /  ( sqrt(2.0*M_PI)*sigma );
}

/* ------- ROUTINE: Ci,Si ------
 * INPUTS: x
 * TASK: Ci,Si
 * OUTPUT: Ci,Si
 * */
static double Ci_integration(double t,void *param)
{
    return -cos(t)/t;
}
static double Si_integration(double t,void *param)
{
    if(t>1.0e-10)
        return sin(t)/t;
    else 
        return 1;
}
double Ci(double x)
{
    gsl_integration_cquad_workspace *workspace =  NULL;
    gsl_function F;
    F.function=&Ci_integration;

    workspace = gsl_integration_cquad_workspace_alloc(1000);
    double a=x;
    double b=a+2.*M_PI;
    double y;
    double result=0.;
	for(int i=0;i<200;i++)
	{
        b=a+2.*M_PI;
        int gslstatus = gsl_integration_cquad(&F,
                                            a,
                                            b,
                                            0.0, 1.0e-5,
                                            workspace,&y,NULL,NULL);
        if(gslstatus != GSL_SUCCESS) 
            CLUCU_RAISE_GSL_WARNING(gslstatus,NULL);
        result+=y;
        a=b;
    }
    gsl_integration_cquad_workspace_free(workspace);
    return result;
}
double Si(double x)
{
    gsl_integration_cquad_workspace *workspace =  NULL;
    gsl_function F;
    F.function=&Si_integration;

    workspace = gsl_integration_cquad_workspace_alloc(1000);
    double a=0;
    double b=0;
    double y;
    double result=0.;
    int N=(int)(x/(2.0*M_PI));
	for(int i=0;i<N;i++)
	{
        b=a+2.*M_PI;
        int gslstatus = gsl_integration_cquad(&F,
                                            a,
                                            b,
                                            0.0, 1.0e-5,
                                            workspace,&y,NULL,NULL);
        if(gslstatus != GSL_SUCCESS) 
            CLUCU_RAISE_GSL_WARNING(gslstatus,NULL);
        result+=y;
        a=b;
    }
    int gslstatus = gsl_integration_cquad(&F,
                                            b,
                                            x,
                                            0.0, 1.0e-5,
                                            workspace,&y,NULL,NULL);
    if(gslstatus != GSL_SUCCESS) 
            CLUCU_RAISE_GSL_WARNING(gslstatus,NULL);
    result+=y;                                        

    gsl_integration_cquad_workspace_free(workspace);
    return result;
}
/* ------- ROUTINE: 矩阵操作 ------
 * INPUTS: 
 * TASK: 
 * OUTPUT: 
 * */
/*----------------------*/
/*      创建动态数组      */
/*----------------------*/
//一边申请一边分配，地址不是连续的
double *Create1Grid(int m)
{
	double *tt =(double *)malloc(sizeof(double) * m);
	return tt;
}
double **Create2Grid(int m, int n)
{
	double **tt =(double **)malloc(sizeof(double) * m);
    for (int i = 0; i < m; i++)
        tt[i] = Create1Grid(n);
	return tt;
}
double ***Create3Grid(int m, int n, int t)
{
	double ***tt = (double ***)malloc(sizeof(double) * m);
    for (int i = 0; i < m; i++)
        tt[i] = Create2Grid(n,t);
	return tt;
}
double ****Create4Grid(int m, int n, int t, int q)
{
	double ****tt = (double ****)malloc(sizeof(double) * m);
	for (int i = 0; i < m; i++)
        tt[i] = Create3Grid(n,t,q);
	return tt;
}
double *****Create5Grid(int m, int n, int t, int q, int x)
{
	double *****tt = (double *****)malloc(sizeof(double) * m);
	for (int i = 0; i < m; i++)
        tt[i] = Create4Grid(n,t,q,x);
	return tt;
}
double ******Create6Grid(int m, int n, int t, int q, int x, int y)
{
	double ******tt = (double ******)malloc(sizeof(double) * m);
	for (int i = 0; i < m; i++)
        tt[i] = Create5Grid(n,t,q,x,y);
	return tt;
}
double *******Create7Grid(int m, int n, int t, int q, int x, int y, int z)
{
	double *******tt = (double *******)malloc(sizeof(double) * m);
	for (int i = 0; i < m; i++)
        tt[i] = Create6Grid(n,t,q,x,y,z);
	return tt;
}
double ********Create8Grid(int m, int n, int t, int q, int x, int y, int z, int v)
{
	double ********tt = (double ********)malloc(sizeof(double) * m);
	for (int i = 0; i < m; i++)
        tt[i] = Create7Grid(n,t,q,x,y,z,v);
	return tt;
}
/*----------------------*/
/*      释放动态数组      */
/*----------------------*/
void Free1Grid(double *tt, int m)
{
	if (tt != NULL)
		free(tt);
}
void Free2Grid(double **tt, int m, int n)
{
	if (tt != NULL)
		for (int i = 0; i < m; i++)
			free(tt[i]);
}
void Free3Grid(double ***tt, int m, int n, int t)
{
	if (tt != NULL)
		for (int i =0; i < m; i++)
			Free2Grid(tt[i],n,t);
}
void Free4Grid(double ****tt, int m, int n, int t, int p)
{
	if (tt != NULL)
		for (int i =0; i < m; i++)
			Free3Grid(tt[i],n,t,p);
}
void Free5Grid(double *****tt, int m, int n, int t, int p, int x)
{
	if (tt != NULL)
		for (int i =0; i < m; i++)
			Free4Grid(tt[i],n,t,p,x);
}
void Free6Grid(double ******tt, int m, int n, int t, int p, int x, int y)
{
	if (tt != NULL)
		for (int i =0; i < m; i++)
			Free5Grid(tt[i],n,t,p,x,y);
}
void Free7Grid(double *******tt, int m, int n, int t, int p, int x, int y, int z)
{
	if (tt != NULL)
	{
		for (int i =0; i < m; i++)
		{
			Free6Grid(tt[i],n,t,p,x,y,z);
		}
	}
}
void Free8Grid(double ********tt, int m, int n, int t, int p, int x, int y, int z, int v)
{
	if (tt != NULL)
	{
		for (int i =0; i < m; i++)
		{
			Free7Grid(tt[i],n,t,p,x,y,z,v);
		}
	}
}
/*----------------------*/
/*        拷贝数组       */
/*----------------------*/
void Copy1Grid(double *from, double *to, int m)
{
    memcpy(to, from,m * sizeof(double));
}

void Copy2Grid(double **from, double **to, int m, int n)
{
	for (int i = 0; i < m; i++)
		Copy1Grid(from[i], to[i], n);
}

void Copy3Grid(double ***from, double ***to, int m, int n, int t)
{
	for (int i = 0; i < m; i++)
		Copy2Grid(from[i], to[i], n, t);
}
void Copy4Grid(double ****from, double ****to, int m, int n, int t, int q)
{
	for (int i = 0; i < m; i++)
		Copy3Grid(from[i], to[i], n, t, q);
}
/*----------------------*/
/*        矩阵计算       */
/*----------------------*/
//翻转向量
void InverseMatrix(double *m, int n)
{
    double *inversed_m=Create1Grid(n);
    for(int i=0;i<n;i++)
        inversed_m[i]=m[n-i-1];
    for(int i=0;i<n;i++)
        m[i]=inversed_m[i];
    Free1Grid(inversed_m,n);
}
//矩阵乘法c=a*b
gsl_matrix *gsl_matrix_dot(gsl_matrix *a,gsl_matrix *b)
{
    if(a->size2!=b->size1)
        return NULL;
    gsl_matrix *c=gsl_matrix_calloc(a->size1,b->size2);
	for (size_t i=0;i<a->size1;i++)
	{
		for (size_t j=0;j<b->size2;j++)
		{
			double sum=0.0;
			for (size_t k=0;k<b->size1;k++)
			{
				sum+=gsl_matrix_get(a,i,k)*gsl_matrix_get(b,k,j);
			}
			gsl_matrix_set(c,i,j,sum);
		}
	}
    return c;
}
//矩阵求逆
gsl_matrix *gsl_matrix_inv(gsl_matrix *a)
{
	size_t n=a->size1;
	size_t m=a->size2;
 
	gsl_matrix *temp1=gsl_matrix_calloc(n,n);
	gsl_matrix_memcpy(temp1,a);
 
	gsl_permutation *p=gsl_permutation_calloc(n);
	int sign=0;
	gsl_linalg_LU_decomp(temp1,p,&sign);
	gsl_matrix *inverse=gsl_matrix_calloc(n,n);
 
	gsl_linalg_LU_invert(temp1,p,inverse);
	
	gsl_permutation_free(p);
	gsl_matrix_free(temp1);
    return inverse;
}
void clucu_matrix_inv(int n,double **mat,double**inv_mat)
{
    gsl_matrix *a = gsl_matrix_alloc(n,n);
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            gsl_matrix_set (a, i, j, mat[i][j]);
    gsl_matrix *inv_a = gsl_matrix_inv(a);
    for(int i=0;i<n;i++)
        for(int j=0;j<n;j++)
            inv_mat[i][j]=gsl_matrix_get(inv_a, i, j);
    gsl_matrix_free(a);
    gsl_matrix_free(inv_a);
}
//矩阵取迹
double gsl_matrix_trace(gsl_matrix *a)
{
	size_t n=a->size1;
	size_t m=a->size2;
    if(n!=m)
        return NAN;
    
    double sum=0.0;
	for (size_t i=0;i<a->size1;i++)
    {
        sum+=gsl_matrix_get(a,i,i);
	}
    return sum;
}

void PrintMatrix(double **m,int n1,int n2)
{
	for (int i = 0; i < n1; i++)
	{
		for (int j = 0; j < n2-1; j++)
			printf("%.15le   ", m[i][j]);
        printf("%.15le\n", m[i][n2-1]);
	}
	printf("\n");
}


/**
 * @brief 
 * 贝塞尔函数
 * Spherical Bessel function of order l (adapted from CAMB)
 */

#define CLUCU_GAMMA1 2.6789385347077476336556 //Gamma(1/3)
#define CLUCU_GAMMA2 1.3541179394264004169452 //Gamma(2/3)
#define CLUCU_ROOTPI12 21.269446210866192327578 //12*sqrt(pi)
//球贝塞尔
double clucu_spherical_bessel(int l,double x)
{
  double jl;
  double ax=fabs(x);
  double ax2=x*x;

  if(l<7) {
    if(l==0) {
      if(ax<0.1) jl=1-ax2*(1-ax2/20.)/6.;
      else jl=sin(x)/x;
    }
    else if(l==1) {
      if(ax<0.2) jl=ax*(1-ax2*(1-ax2/28)/10)/3;
      else jl=(sin(x)/ax-cos(x))/ax;
    }
    else if(l==2) {
      if(ax<0.3) jl=ax2*(1-ax2*(1-ax2/36)/14)/15;
      else jl=(-3*cos(x)/ax-sin(x)*(1-3/ax2))/ax;
    }
    else if(l==3) {
      if(ax<0.4)
	jl=ax*ax2*(1-ax2*(1-ax2/44)/18)/105;
      else
	jl=(cos(x)*(1-15/ax2)-sin(x)*(6-15/ax2)/ax)/ax;
    }
    else if(l==4) {
      if(ax<0.6)
	jl=ax2*ax2*(1-ax2*(1-ax2/52)/22)/945;
      else
	jl=(sin(x)*(1-(45-105/ax2)/ax2)+cos(x)*(10-105/ax2)/ax)/ax;
    }
    else if(l==5) {
      if(ax<1.0)
	jl=ax2*ax2*ax*(1-ax2*(1-ax2/60)/26)/10395;
      else {
	jl=(sin(x)*(15-(420-945/ax2)/ax2)/ax-
	    cos(x)*(1-(105-945/ax2)/ax2))/ax;
      }
    }
    else {
      if(ax<1.0)
	jl=ax2*ax2*ax2*(1-ax2*(1-ax2/68)/30)/135135;
      else {
	jl=(sin(x)*(-1+(210-(4725-10395/ax2)/ax2)/ax2)+
	    cos(x)*(-21+(1260-10395/ax2)/ax2)/ax)/ax;
      }
    }
  }
  else {
    double nu=l+0.5;
    double nu2=nu*nu;

    if(ax<1.0E-40) jl=0;
    else if((ax2/l)<0.5) {
      jl=(exp(l*log(ax/nu)-M_LN2+nu*(1-M_LN2)-(1-(1-3.5/nu2)/(30*nu2))/(12*nu))/nu)*
	(1-ax2/(4*nu+4)*(1-ax2/(8*nu+16)*(1-ax2/(12*nu+36))));
    }
    else if((l*l/ax)<0.5) {
      double beta=ax-0.5*M_PI*(l+1);
      jl=(cos(beta)*(1-(nu2-0.25)*(nu2-2.25)/(8*ax2)*(1-(nu2-6.25)*(nu2-12.25)/(48*ax2)))-
	  sin(beta)*(nu2-0.25)/(2*ax)*(1-(nu2-2.25)*(nu2-6.25)/(24*ax2)*
				       (1-(nu2-12.25)*(nu2-20.25)/(80*ax2))))/ax;
    }
    else {
      double l3=pow(nu,0.325);
      if(ax<nu-1.31*l3) {
	double cosb=nu/ax;
	double sx=sqrt(nu2-ax2);
	double cotb=nu/sx;
	double secb=ax/nu;
	double beta=log(cosb+sx/ax);
	double cot3b=cotb*cotb*cotb;
	double cot6b=cot3b*cot3b;
	double sec2b=secb*secb;
	double expterm=((2+3*sec2b)*cot3b/24
			-((4+sec2b)*sec2b*cot6b/16
			  +((16-(1512+(3654+375*sec2b)*sec2b)*sec2b)*cot3b/5760
			    +(32+(288+(232+13*sec2b)*sec2b)*sec2b)*sec2b*cot6b/(128*nu))*
			  cot6b/nu)/nu)/nu;
	jl=sqrt(cotb*cosb)/(2*nu)*exp(-nu*beta+nu/cotb-expterm);
      }
      else if(ax>nu+1.48*l3) {
	double cosb=nu/ax;
	double sx=sqrt(ax2-nu2);
	double cotb=nu/sx;
	double secb=ax/nu;
	double beta=acos(cosb);
	double cot3b=cotb*cotb*cotb;
	double cot6b=cot3b*cot3b;
	double sec2b=secb*secb;
	double trigarg=nu/cotb-nu*beta-0.25*M_PI-
	  ((2+3*sec2b)*cot3b/24+(16-(1512+(3654+375*sec2b)*sec2b)*sec2b)*
	   cot3b*cot6b/(5760*nu2))/nu;
	double expterm=((4+sec2b)*sec2b*cot6b/16-
			(32+(288+(232+13*sec2b)*sec2b)*sec2b)*
			sec2b*cot6b*cot6b/(128*nu2))/nu2;
	jl=sqrt(cotb*cosb)/nu*exp(-expterm)*cos(trigarg);
      }
      else {
	double beta=ax-nu;
	double beta2=beta*beta;
	double sx=6/ax;
	double sx2=sx*sx;
	double secb=pow(sx,0.3333333333333333333333);
	double sec2b=secb*secb;

	jl=(CLUCU_GAMMA1*secb+beta*CLUCU_GAMMA2*sec2b
	    -(beta2/18-1.0/45.0)*beta*sx*secb*CLUCU_GAMMA1
	    -((beta2-1)*beta2/36+1.0/420.0)*sx*sec2b*CLUCU_GAMMA2
	    +(((beta2/1620-7.0/3240.0)*beta2+1.0/648.0)*beta2-1.0/8100.0)*sx2*secb*CLUCU_GAMMA1
	    +(((beta2/4536-1.0/810.0)*beta2+19.0/11340.0)*beta2-13.0/28350.0)*beta*sx2*sec2b*CLUCU_GAMMA2
	    -((((beta2/349920-1.0/29160.0)*beta2+71.0/583200.0)*beta2-121.0/874800.0)*
	      beta2+7939.0/224532000.0)*beta*sx2*sx*secb*CLUCU_GAMMA1)*sqrt(sx)/CLUCU_ROOTPI12;
      }
    }
  }
  if((x<0)&&(l%2!=0)) jl=-jl;

  return jl;
}

void save_matrix(FILE *fp,double **matrix,double *x,double *y,int nx,int ny)
{
    fprintf(fp,"0.0        ");
	for(int ix=0;ix<nx;ix++)
		fprintf(fp,"%f               ",x[ix]);
	fprintf(fp,"\n");
	for(int iy=0;iy<ny;iy++)
	{
		fprintf(fp,"%f   ",y[iy]);
		for(int ix=0;ix<nx;ix++)
			fprintf(fp,"%.15le   ",matrix[ix][iy]);
		fprintf(fp,"\n");
	}
}


/* ------- ROUTINE: 数组求积分 ------
 * INPUTS: 
 *          n: 数组大小
 *          *x: 指向一维数组x，对x积分
 *          *y: y(x)
 *          a,b: 积分区间，如果b<a，则对所有x积分
 *          T：积分方法，一般为gsl_integration_cquad
 * TASK: int y(x) dx
 * OUTPUT: double 积分结果result
 * */
double clucu_integ_spline(int n,double *x,double *y,
                      double a, double b,
                      const gsl_interp_type *T)
{
    double result;
    if(b==a)
    {
        return 0.;
    }
    if(b<a)
    {
        b=x[n-1];
        a=x[0];
    }

    if((b>x[n-1]) || (a<x[0]))
    {
        CLUCU_RAISE_WARNING(CLUCU_ERROR_SPLINE,
                            "clucu_utils.c: clucu_integ_spline(): "
                            "integration limits beyond interpolated range (a,b=[%e,%e],x=[%e,%e])\n",a,b,x[0],x[n-1]);
    }

    gsl_interp_accel *ia = NULL;
    gsl_spline *s = NULL;
    s = gsl_spline_alloc(T, n);
    if(s == NULL)
        CLUCU_RAISE_WARNING(CLUCU_ERROR_SPLINE,
                        "clucu_utils.c: clucu_integ_spline(): "
                        "");

    ia = gsl_interp_accel_alloc();
    if(ia == NULL)
        CLUCU_RAISE_WARNING(CLUCU_ERROR_MEMORY,
                        "clucu_utils.c: clucu_integ_spline(): "
                        "");

    if(gsl_spline_init(s, x, y, n)) 
        CLUCU_RAISE_WARNING(CLUCU_ERROR_SPLINE,
                    "clucu_utils.c: clucu_integ_spline(): "
                    "");

    int sstat = gsl_spline_eval_integ_e(s, a, b, ia, &result);
    if(sstat)
        CLUCU_RAISE_WARNING(CLUCU_ERROR_SPLINE_EV,
                "clucu_utils.c: clucu_integ_spline(): "
                "");

    gsl_spline_free(s);
    gsl_interp_accel_free(ia);

    return result;
}