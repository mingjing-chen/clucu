#include "clucu_power.h"


/*-----------------------------------------------------------------------------*/
/*                                       EH98                                  */
/*-----------------------------------------------------------------------------*/
/*------ ROUTINE: tsqr_BBKS -----
INPUT: clucu_parameters and k wavenumber in Mpc^-1
TASK: provide the square of the BBKS transfer function with baryonic correction
NOTE: Bardeen et al. (1986) as implemented in Sugiyama (1995)
*/
//k in 1/Mpc
double tsqr_BBKS(clucu_cosmology *cosmo, double k) 
{
  double tfac = cosmo->param.T_CMB / 2.7;
  double q = tfac * tfac * k / (
    cosmo->param.Omega_m * cosmo->param.h * cosmo->param.h *
    exp(-cosmo->param.Omega_b * (1.0 + pow(2. * cosmo->param.h, .5) / cosmo->param.Omega_m)));
  return (
    pow(log(1. + 2.34*q) / (2.34*q), 2.0) /
    pow(1. + 3.89*q + pow(16.1*q, 2.0) + pow(5.46*q, 3.0) + pow(6.71*q, 4.0), 0.5));
}
/*------ ROUTINE: bbks_power -----
INPUT: clucu_parameters and k wavenumber in 1/Mpc
TASK: compute the unnormalized BBKS power spectrum
*/
double clucu_bbks_power(clucu_cosmology *cosmo, double k) {

    double kp = cosmo->param.k_pivot;//1/Mpc
    double n_1 = (cosmo->param.n_s-1.)+0.5*cosmo->param.alpha_s*log(k/kp);
    double As;
    if(!(isnan(cosmo->param.sigma8)) && isnan(cosmo->param.delta_zeta))
        As=1.;
    if(isnan(cosmo->param.sigma8) && !(isnan(cosmo->param.delta_zeta)))
        As = cosmo->param.delta_zeta * cosmo->param.delta_zeta;
int status=0;
    double H0 = clucu_hubble_parameter(cosmo,0.,&status);
    double Om = cosmo->param.Omega_m;
//    double pre = pow( 2.*k*k/5./H0/H0/Om,2.);
    double pre = pow( 0.3*k*k/H0/H0/Om,2.);

    double Delta2_k = As * pre * pow(k/kp,n_1) * tsqr_BBKS(cosmo, k);
    double result = 2.*M_PI*M_PI/k/k/k *Delta2_k;
    return result;
}
double clucu_bbks_power_sigma8(clucu_cosmology *cosmo, double k) {
  return pow(k, cosmo->param.n_s) * tsqr_BBKS(cosmo, k);
}
/*
 * Allocate a new struct for storing EH98 data
 * @param params Cosmological parameters
 * @param int, include BAO wiggles if not 0, smooth otherwise
 */
eh_struct* clucu_eh_struct_new(clucu_cosmology *cosmo, int wiggled) {
  //////
  // Computes Eisenstein & Hu parameters for
  // P_k and r_sound
  // see astro-ph/9709112 for the relevant equations
  double OMh2,OBh2;
  double th2p7;
  eh_struct *eh=malloc(sizeof(eh_struct));
  if(eh==NULL)
    return NULL;

  OMh2=cosmo->param.Omega_m*cosmo->param.h*cosmo->param.h; //Cosmo params scaled by h^2
  OBh2=cosmo->param.Omega_b*cosmo->param.h*cosmo->param.h;
  th2p7=cosmo->param.T_CMB/2.7; //This is Theta_{2.7} in E&Hu notation
  eh->th2p7=th2p7; //This is Theta_{2.7} in E&Hu notation
  eh->zeq=2.5E4*OMh2/pow(th2p7,4); //Eq. 2
  eh->keq=0.0746*OMh2/(th2p7*th2p7); //Eq. 3

  //This group corresponds to Eq. 4
  double b1,b2;
  b1=0.313*pow(OMh2,-0.419)*(1+0.607*pow(OMh2,0.674));
  b2=0.238*pow(OMh2,0.223);
  eh->zdrag=1291*pow(OMh2,0.251)*(1+b1*pow(OBh2,b2))/(1+0.659*pow(OMh2,0.828));

  //These are the baryon-to-photon ratios
  //at equality (Req) and drag (Rd) epochs
  //Eq. 5
  double Req,Rd;
  Req=31.5*OBh2*1000./(eh->zeq*pow(th2p7,4));
  Rd=31.5*OBh2*1000./((1+eh->zdrag)*pow(th2p7,4));
  eh->rsound=2/(3*eh->keq)*sqrt(6/Req)*
    log((sqrt(1+Rd)+sqrt(Rd+Req))/(1+sqrt(Req)));

  //This is Eq. 7 (but in 1/Mpc)
  eh->kSilk=1.6*pow(OBh2,0.52)*pow(OMh2,0.73)*(1+pow(10.4*OMh2,-0.95));

  //These are Eqs. 11
  double a1,a2,b_frac;
  a1=pow(46.9*OMh2,0.670)*(1+pow(32.1*OMh2,-0.532));
  a2=pow(12.0*OMh2,0.424)*(1+pow(45.0*OMh2,-0.582));
  b_frac=OBh2/OMh2;
  eh->alphac=pow(a1,-b_frac)*pow(a2,-b_frac*b_frac*b_frac);

  //These are Eqs. 12
  double bb1,bb2;
  bb1=0.944/(1+pow(458*OMh2,-0.708));
  bb2=pow(0.395*OMh2,-0.0266);
  eh->betac=1/(1+bb1*(pow(1-b_frac,bb2)-1));

  double y=eh->zeq/(1+eh->zdrag);
  double sqy=sqrt(1+y);
  double gy=y*(-6*sqy+(2+3*y)*log((sqy+1)/(sqy-1))); //Eq 15
  //Baryon suppression Eq. 14
  eh->alphab=2.07*eh->keq*eh->rsound*pow(1+Rd,-0.75)*gy;

  //Baryon envelope shift Eq. 24
  eh->betab=0.5+b_frac+(3-2*b_frac)*sqrt(pow(17.2*OMh2,2)+1);

  //Node shift parameter Eq. 23
  eh->bnode=8.41*pow(OMh2,0.435);

  //Approx for the sound horizon, Eq. 26
  eh->rsound_approx=cosmo->param.h*44.5*log(9.83/OMh2)/
    sqrt(1+10*pow(OBh2,0.75));

  eh->wiggled=wiggled;

  (eh_struct *)cosmo->eh;
  cosmo->eh = eh;
  cosmo->computed_transfer;

  return eh;
}

static double tkEH_0(double keq,double k,double a,double b)
{
  //////
  // Eisentstein & Hu's Tk_0
  // see astro-ph/9709112 for the relevant equations
  double q=k/(13.41*keq); //Eq 10
  double c=14.2/a+386./(1+69.9*pow(q,1.08)); //Eq 20
  double l=log(M_E+1.8*b*q); //Change of var for Eq 19
  return l/(l+c*q*q); //Returns Eq 19
}

static double tkEH_c(eh_struct *eh,double k)
{
  //////
  // Eisenstein & Hu's Tk_c
  // see astro-ph/9709112 for the relevant equations
  double f=1/(1+pow(k*eh->rsound/5.4,4)); //Eq 18
  double aa=tkEH_0(eh->keq,k,1,eh->betac);
  double bb=tkEH_0(eh->keq,k,eh->alphac,eh->betac);
  double cc=f*tkEH_0(eh->keq,k,1,eh->betac)+
    (1-f)*tkEH_0(eh->keq,k,eh->alphac,eh->betac);
  return f*tkEH_0(eh->keq,k,1,eh->betac)+
    (1-f)*tkEH_0(eh->keq,k,eh->alphac,eh->betac); //Returns Eq 17
}

static double jbes0(double x)
{
  double jl;
  double ax2=x*x;

  if(ax2<1e-4) jl=1-ax2*(1-ax2/20.)/6.;
  else jl=sin(x)/x;

  return jl;
}

static double tkEH_b(eh_struct *eh,double k)
{
  //////
  // Eisenstein & Hu's Tk_b (Eq 21)
  // see astro-ph/9709112 for the relevant equations
  double x_bessel,part1,part2;
  double x=k*eh->rsound;

  //First term of Eq. 21
  if(k==0) x_bessel=0;
  else {
    x_bessel=x*pow(1+eh->bnode*eh->bnode*eh->bnode/(x*x*x),-1./3.);
  }

  part1=tkEH_0(eh->keq,k,1,1)/(1+pow(x/5.2,2));

  //Second term of Eq. 21
  if(k==0)
    part2=0;
  else
    part2=eh->alphab/(1+pow(eh->betab/x,3))*exp(-pow(k/eh->kSilk,1.4));
    double aa=jbes0(x_bessel);
    double bb=jbes0(x_bessel)*(part1+part2);
  return jbes0(x_bessel)*(part1+part2);
}

//k in 1/Mpc
double tsqr_EH(clucu_cosmology *cosmo,eh_struct *eh,double k)
{
  //////
  // Eisenstein & Hu's Tk_c
  // see astro-ph/9709112 for the relevant equations
  // Notice the last parameter in eh_power controls
  // whether to introduce wiggles (BAO) in the power spectrum.
  // We do this by default when obtaining the power spectrum.
  double tk;
  double b_frac=cosmo->param.Omega_b/cosmo->param.Omega_m;
  if(eh->wiggled)
    //Case with baryons (Eq 8)
    tk=b_frac*tkEH_b(eh,k)+(1-b_frac)*tkEH_c(eh,k);
  else {
    //Zero baryon case (sec 4.2)
    double OMh2=cosmo->param.Omega_m*cosmo->param.h*cosmo->param.h;
    // Compute Eq. 31
    double alpha_gamma=1-0.328*log(431*OMh2)*b_frac+0.38*log(22.3*OMh2)*b_frac*b_frac;
    // Compute Eq. 30
    double gamma_eff=cosmo->param.Omega_m*cosmo->param.h*(alpha_gamma+(1-alpha_gamma)/
						(1+pow(0.43*k*eh->rsound_approx,4)));
    // Compute Eq. 28 (assume k in 1/Mpc)
    double q=k*eh->th2p7*eh->th2p7/gamma_eff;
    // Compute Eqs. 29
    double l0=log(2*M_E+1.8*q);
    double c0=14.2+731/(1+62.5*q);
    tk=l0/(l0+c0*q*q);  //T_0 of Eq. 29
  }

  return tk*tk; //Return T_0^2
}

/*
 * Compute the Eisenstein & Hu (1998) unnormalize power spectrum
 * @param params Cosmological parameters
 * @param p, an eh_struct instance
 * @param k, wavenumber in Mpc^-1
 */
double clucu_eh_power_sigma8(clucu_cosmology *cosmo, eh_struct* eh, double k) 
{
  //Wavenumber in units of Mpc^-1
  return pow(k, cosmo->param.n_s) * tsqr_EH(cosmo, eh, k);
}
double clucu_eh_power(clucu_cosmology *cosmo, eh_struct* eh, double k) {
  //Wavenumber in units of Mpc^-1
    double kp = cosmo->param.k_pivot;//1/Mpc
    double n_1 = cosmo->param.n_s-1.+0.5*cosmo->param.alpha_s*log(k/kp);
    double As;
    if(!(isnan(cosmo->param.sigma8)) && isnan(cosmo->param.delta_zeta))
        As=1.;
    if(isnan(cosmo->param.sigma8) && !(isnan(cosmo->param.delta_zeta)))
        As = cosmo->param.delta_zeta * cosmo->param.delta_zeta;

    double H0 = cosmo->param.h*100.* 1000. /clucu_constants.CLIGHT;
    double Om = cosmo->param.Omega_m;
    double pre = pow( 0.4*k*k/H0/H0/Om,2.);

    double Delta2_k = As * pre * pow(k/kp,n_1) * tsqr_EH(cosmo, eh, k);
    double result = 2.*M_PI*M_PI/k/k/k *Delta2_k;
    return result /(0.6*0.6)*(4.546324222012e-01*4.546324222012e-01);
    //return result;
}

//TODO:这里加D(z)好像有问题。体现在，用EH算功率谱和用CLASS算功率谱，得到的星系团数，红移越高差别越大
//TODO: 如果用EH的方法计算功率谱，还需要插值吗
//输入pk(不是log pk)，内部自动乘D(z)
static clucu_f2d_t *clucu_compute_linpower_analytic(clucu_cosmology *cosmo, void* par,
                                                double (*pk)(clucu_cosmology *cosmo,
                                                             void* p, double k),
                                                int* status) 
{
    if(isnan(cosmo->param.sigma8) && isnan(cosmo->param.delta_zeta))
    {
        *status = CLUCU_ERROR_INCONSISTENT;
        CLUCU_RAISE_WARNING(*status, "sigma8 not set, required for analytic power spectra");
        clucu_cosmology_set_status_message(cosmo,
                "clucu_power.c: clucu_compute_linpower_analytic(): "
                "sigma8 not set, required for analytic power spectra\n");
        return NULL;
    }
    clucu_f2d_t *psp_out = NULL;
    double sigma8,log_sigma8;
    //These are the limits of the splining range
    double kmin = cosmo->spline_param.K_MIN;
    double kmax = cosmo->spline_param.K_MAX;
    //Compute nk from number of decades and N_K = # k per decade
    double ndecades = log10(kmax) - log10(kmin);
    int nk = (int)ceil(ndecades*cosmo->spline_param.N_K);
    // Compute na using predefined spline spacing
    int nz = cosmo->data.spline_Nz;
    double *z = cosmo->data.spline_z;


    // The x array is initially k, but will later
    // be overwritten with log(k)
    double *x=NULL, *y=NULL, *y2d=NULL;
    x=clucu_log_spacing(kmin, kmax, nk);
    y=malloc(sizeof(double)*nk);
    y2d = malloc(nk * nz * sizeof(double));

    for (int i=0; i<nk; i++) 
    {
        y[i] = log((*pk)(cosmo, par, x[i]));//计算ln(pk), x=k
        x[i] = log(x[i]);//x=ln k
    }  

    for (int j = 0; j < nz; j++) 
    {
        double gfac = clucu_growth_factor(cosmo,z[j],status);
        double g2 = 2.*log(gfac);
        for (int i=0; i<nk; i++) {
        y2d[j*nk+i] = y[i]+g2;// ln(pk)+2ln(Dz)=ln(pk*Dz^2)
        }
    }

    if(*status==0) 
    {
    //输入lnk,z,ln(pkz)
    //is_factorizable=0，因为我们不是分开输入pk,pz
    //大k外推法：1，logk 一阶导数
    //小k外推法：2，logk 二阶导数
    //大z外推法：用clucugrowth算增长因子D, D ** growth_exponent
    //is_fka_log=1，输入的ln(pkz)是log形式，在外推时会判断是相乘还是相加。建议都输入log形式
    //大z外推法：指定增长因子为常数growth_factor_0=0，这个优先级低于clucu_f2d_clucugrowth
    //growth_exponent=2
    //二维插值方法：Bicubic
    psp_out=clucu_f2d_t_new(nz,z,nk,x,y2d,NULL,NULL,0,
                            1,2,clucu_f2d_clucugrowth,1,0,2,
                            clucu_f2d_3,status);
    }
    if(!(isnan(cosmo->param.sigma8)) && (isnan(cosmo->param.delta_zeta)))
    {
        //如果cosmo参数给的是sigma8，需要用sigma8，对log(pkz)做归一化
        if(*status==0) {
            sigma8 = clucu_sigma8(cosmo, psp_out, status);
        }
        if(*status==0) {
            // Calculate normalization factor using computed value of sigma8, then
            // recompute P(k, a) using this normalization
            log_sigma8 = 2*(log(cosmo->param.sigma8) - log(sigma8));
            for(int j=0;j<nz*nk;j++)
            y2d[j] += log_sigma8;
        }
        if(*status==0) {
            // Free the previous P(k,a) spline, and allocate a new one to store the
            // properly-normalized P(k,a)
            clucu_f2d_t_free(psp_out);
            psp_out = clucu_f2d_t_new(nz,z,nk,x,y2d,NULL,NULL,0,
                                    1,2,clucu_f2d_clucugrowth,1,0,2,
                                    clucu_f2d_3,status);
        }
    }
    else
        cosmo->param.sigma8 = clucu_sigma8(cosmo, psp_out, status);

    free(x);
    free(y);
    free(y2d);
    return psp_out;
}
// helper functions for BBKS and EH98
static double bbks_power(clucu_cosmology *cosmo, void *p, double k) 
{
    return clucu_bbks_power(cosmo, k);
}

static double eh_power(clucu_cosmology *cosmo, void *p, double k) 
{
    return clucu_eh_power(cosmo, (eh_struct*)p, k);
}
void clucu_compute_linpower_eh(clucu_cosmology *cosmo, int wiggled, int *status)
{
    clucu_f2d_t *psp = NULL;
    eh_struct *eh = NULL;
    eh = clucu_eh_struct_new(cosmo,wiggled);
    psp=clucu_compute_linpower_analytic(cosmo, eh,
                                    eh_power,
                                    status);
    cosmo->data.log_pkz = psp;
    cosmo->computed_power = true;
}

void clucu_compute_linpower_bbks(clucu_cosmology *cosmo, int *status)
{
    clucu_f2d_t *psp=clucu_compute_linpower_analytic(cosmo, NULL, bbks_power, status);
    cosmo->data.log_pkz = psp;
    cosmo->computed_power = true;
}
/*-----------------------------------------------------------------------------*/
/*                                       CLASS                                  */
/*-----------------------------------------------------------------------------*/

void clucu_compute_power_class_Pk(clucu_cosmology *cosmo, int *status)
{   
    if(cosmo->runned_class==false)
    {
        RunClass(cosmo);
    }
	//**pkz
    int nk=10000.;//先给个够大的Nk
    int nz=cosmo->data.spline_Nz;
    double *lnk=Create1Grid(nk);
    double **lnpkz=Create2Grid(nk,nz);
    double *y2d = malloc(nk * nz * sizeof(double));
    double kk,pp;
	for(int iz=0;iz<nz;iz++)
	{
		char nameclassout[150];
		if(cosmo->param.N_nu_mass > 0) 
			sprintf(nameclassout,"../data/%s/p%d_z%d_pk_cb.dat", cosmo->classname, cosmo->class_id,iz+1);
		else 
			sprintf(nameclassout,"../data/%s/p%d_z%d_pk.dat", cosmo->classname, cosmo->class_id,iz+1);
		FILE* fp=fopen(nameclassout,"r");
    	if(! fp)
        {
            cosmo->computed_power = false;
            return;
        }
		int ik=0;
		while(!feof(fp))
		{
			char a[1000];
			if(ik<4)
				fgets(a,1000,fp);
			else
			{
				fscanf(fp,"%lf  %lf\n",&kk,&pp);
                lnk[ik-4] =  log( kk *cosmo->param.h );// class输出格式带h，h/Mpc，换成不带h
                lnpkz[ik-4][iz] = log( pp/pow(cosmo->param.h,3.) );

			}
			ik++;
		}
		if(nk > ik-4) nk = ik-4;
		fclose(fp);
	}

    clucu_f2d_t *psp_out = NULL;

    for (int jz = 0; jz < nz; jz++) 
    {
        for (int ik = 0; ik < nk; ik++) 
        {
            y2d[jz*nk+ik] = lnpkz[ik][jz];// ln(pk)+2ln(Dz)=ln(pk*Dz^2)
        }
    }
  
  if(*status==0) {
    //输入lnk,z,ln(pkz)
    //is_factorizable=0，因为我们不是分开输入pk,pz
    //大k外推法：1，logk 一阶导数
    //小k外推法：2，logk 二阶导数
    //大z外推法：用clucu_f2d_clucugrowth算增长因子D, D ** growth_exponent
    //is_fka_log=1，输入的ln(pkz)是log形式，在外推时会判断是相乘还是相加。建议都输入log形式
    //大z外推法：指定增长因子为常数growth_factor_0=0，需要设置clucu_f2d_clucugrowth
    //growth_exponent=2
    //二维插值方法：Bicubic
    psp_out=clucu_f2d_t_new(nz,cosmo->data.spline_z,nk,lnk,y2d,NULL,NULL,0,
                          1,2,clucu_f2d_clucugrowth,1,0,2,
                          clucu_f2d_3,status);
    }
    if((isnan(cosmo->param.sigma8)) && !(isnan(cosmo->param.delta_zeta)))
        cosmo->param.sigma8 = clucu_sigma8(cosmo, psp_out, status);
        
    free(lnk);
    Free2Grid(lnpkz,nk,nz);
    free(y2d);
    cosmo->data.log_pkz = psp_out;

    cosmo->computed_power = true;
}

void clucu_compute_power_class_Tk(clucu_cosmology *cosmo, int *status)
{   
    if(cosmo->runned_class==false)
    {
        RunClass(cosmo);
    }
	//**pkz
    int nk=10000.;//先给个够大的Nk
    int nz=cosmo->data.spline_Nz;
    double *lnk=Create1Grid(nk);
    //&tg,&tb,&tc,&tur,   &tnu1,&tnu2,&tnu3
    double **Tg=Create2Grid(nk,nz);
    double **Tcb=Create2Grid(nk,nz);
    double **Tur=Create2Grid(nk,nz);
    double **Tnu1=Create2Grid(nk,nz);
    double **Tnu2=Create2Grid(nk,nz);
    double **Tnu3=Create2Grid(nk,nz);
    double **Pcb=Create2Grid(nk,nz);

    

    double *y2d_g = malloc(nk * nz * sizeof(double));
    double *y2d_ur = malloc(nk * nz * sizeof(double));
    double *y2d_cb = malloc(nk * nz * sizeof(double));
    double *y2d_nu1 = malloc(nk * nz * sizeof(double));
    double *y2d_nu2 = malloc(nk * nz * sizeof(double));
    double *y2d_nu3 = malloc(nk * nz * sizeof(double));
    double *y2d = malloc(nk * nz * sizeof(double));
    double kk,tg,tb,tc,tur,tnu1,tnu2,tnu3;
    if(!(isnan(cosmo->param.sigma8)) && (isnan(cosmo->param.delta_zeta)))
        cosmo->param.A_s=1.;

	for(int iz=0;iz<nz;iz++)
	{
		char nameclassout[150];
			sprintf(nameclassout,"../data/%s/p%d_z%d_tk.dat", cosmo->classname, cosmo->class_id,iz+1);
		FILE* fp=fopen(nameclassout,"r");
    	if(! fp)
        {
            cosmo->computed_power = false;
            return;
        }
		int ik=0;
        int skiprows=9;
		while(!feof(fp))
		{
			char a[1000];
			if(ik<skiprows)
				fgets(a,1000,fp);
			else
			{
                //有 无动态暗能量，fld都在第5列
                if(cosmo->param.N_nu_mass==0)
                    fscanf(fp,"%lf  %lf  %lf  %lf  %*lf       %lf  %*lf  %*lf  %*lf\n",&kk,&tg,&tb,&tc,&tur   );
                else if(cosmo->param.N_nu_mass==1)
                    fscanf(fp,"%lf  %lf  %lf  %lf  %*lf       %lf  %lf  %*lf  %*lf  %*lf\n",&kk,&tg,&tb,&tc,&tur,   &tnu1);
                else if(cosmo->param.N_nu_mass==2)
                    fscanf(fp,"%lf  %lf  %lf  %lf  %*lf       %lf  %lf  %lf  %*lf  %*lf  %*lf\n",&kk,&tg,&tb,&tc,&tur,   &tnu1,&tnu2);
                else if(cosmo->param.N_nu_mass==3)
                    fscanf(fp,"%lf  %lf  %lf  %lf  %*lf       %lf  %lf  %lf  %lf  %*lf  %*lf  %*lf\n",&kk,&tg,&tb,&tc,&tur,   &tnu1,&tnu2,&tnu3);
                kk = kk*cosmo->param.h;// class输出格式带h，h/Mpc，换成不带h
                lnk[ik-skiprows] =  log( kk  );
                Tg[ik-skiprows][iz] = tg;
                Tcb[ik-skiprows][iz] = cosmo->param.Omega_c/(cosmo->param.Omega_c+cosmo->param.Omega_b) *tc
                                        +cosmo->param.Omega_b/(cosmo->param.Omega_c+cosmo->param.Omega_b) *tb;
                Tur[ik-skiprows][iz] = tur;
                Tnu1[ik-skiprows][iz] = tnu1;
                Tnu2[ik-skiprows][iz] = tnu2;
                Tnu3[ik-skiprows][iz] = tnu3;
                Pcb[ik-skiprows][iz]=(2.*M_PI*M_PI)*cosmo->param.A_s*pow(kk/cosmo->param.k_pivot,cosmo->param.n_s-1.)/kk/kk/kk*Tcb[ik-skiprows][iz]*Tcb[ik-skiprows][iz];
                //Tg,Tur，Tnu,有>0有<0。对log(k)～T插值，一阶导数
                //Tcd,都是<0，换成>0，对log(k)~log(-T)，插值，一阶导数
			}
			ik++;
		}
		if(nk > ik-skiprows) nk = ik-skiprows;
		fclose(fp);
	}

    clucu_f2d_t *tsp_g = NULL;
    clucu_f2d_t *tsp_cb = NULL;
    clucu_f2d_t *tsp_ur = NULL;
    clucu_f2d_t *tsp_nu1 = NULL;
    clucu_f2d_t *tsp_nu2 = NULL;
    clucu_f2d_t *tsp_nu3 = NULL;
    clucu_f2d_t *psp = NULL;

    for (int jz = 0; jz < nz; jz++) 
    {
        for (int ik = 0; ik < nk; ik++) 
        {
            y2d_g[jz*nk+ik] = Tg[ik][jz];// ln(pk)+2ln(Dz)=ln(pk*Dz^2)
            y2d_ur[jz*nk+ik] = Tur[ik][jz];
            y2d_cb[jz*nk+ik] = log(-Tcb[ik][jz]);
            y2d_nu1[jz*nk+ik] = Tnu1[ik][jz];
            y2d_nu2[jz*nk+ik] = Tnu2[ik][jz];
            y2d_nu3[jz*nk+ik] = Tnu3[ik][jz];
            y2d[jz*nk+ik] = log(Pcb[ik][jz]);
        }
    }
  
  if(*status==0) {
    //输入lnk,z,ln(pkz)
    //is_factorizable=0，因为我们不是分开输入pk,pz
    //大k外推法：1，logk 一阶导数
    //小k外推法：2，logk 二阶导数
    //大z外推法：用clucu_f2d_clucugrowth算增长因子D, D ** growth_exponent
    //is_fka_log=1，输入的ln(pkz)是log形式，在外推时会判断是相乘还是相加。建议都输入log形式
    //大z外推法：指定增长因子为常数growth_factor_0=0，需要设置clucu_f2d_clucugrowth
    //growth_exponent=2
    //二维插值方法：Bicubic
    tsp_g=clucu_f2d_t_new(nz,cosmo->data.spline_z,nk,lnk,y2d_g,NULL,NULL,0,
                          1,2,clucu_f2d_clucugrowth,0,0,2,
                          clucu_f2d_3,status);
    tsp_ur=clucu_f2d_t_new(nz,cosmo->data.spline_z,nk,lnk,y2d_ur,NULL,NULL,0,
                          1,2,clucu_f2d_clucugrowth,0,0,2,
                          clucu_f2d_3,status);
    tsp_cb=clucu_f2d_t_new(nz,cosmo->data.spline_z,nk,lnk,y2d_cb,NULL,NULL,0,
                          1,2,clucu_f2d_clucugrowth,1,0,2,
                          clucu_f2d_3,status);
    tsp_nu1=clucu_f2d_t_new(nz,cosmo->data.spline_z,nk,lnk,y2d_nu1,NULL,NULL,0,
                          1,2,clucu_f2d_clucugrowth,0,0,2,
                          clucu_f2d_3,status);
    tsp_nu2=clucu_f2d_t_new(nz,cosmo->data.spline_z,nk,lnk,y2d_nu2,NULL,NULL,0,
                          1,2,clucu_f2d_clucugrowth,0,0,2,
                          clucu_f2d_3,status);
    tsp_nu3=clucu_f2d_t_new(nz,cosmo->data.spline_z,nk,lnk,y2d_nu3,NULL,NULL,0,
                          1,2,clucu_f2d_clucugrowth,0,0,2,
                          clucu_f2d_3,status);
    psp=clucu_f2d_t_new(nz,cosmo->data.spline_z,nk,lnk,y2d,NULL,NULL,0,
                          1,2,clucu_f2d_clucugrowth,1,0,2,
                          clucu_f2d_3,status);                                        
    }
    double sigma8,log_sigma8;
    if(!(isnan(cosmo->param.sigma8)) && (isnan(cosmo->param.delta_zeta)))
    {
        //如果cosmo参数给的是sigma8，需要用sigma8，对log(pkz)做归一化
        if(*status==0) {
            sigma8 = clucu_sigma8(cosmo, psp, status);
        }
        if(*status==0) {
            // Calculate normalization factor using computed value of sigma8, then
            // recompute P(k, a) using this normalization
            log_sigma8 = 2*(log(cosmo->param.sigma8) - log(sigma8));
            for(int j=0;j<nz*nk;j++)
            y2d[j] += log_sigma8;
        }
        if(*status==0) {
            // Free the previous P(k,a) spline, and allocate a new one to store the
            // properly-normalized P(k,a)
            clucu_f2d_t_free(psp);
            psp=clucu_f2d_t_new(nz,cosmo->data.spline_z,nk,lnk,y2d,NULL,NULL,0,
                          1,2,clucu_f2d_clucugrowth,1,0,2,
                          clucu_f2d_3,status);
        }
    }
    else
        cosmo->param.sigma8 = clucu_sigma8(cosmo, psp, status);
        
    free(lnk);
    Free2Grid(Pcb,nk,nz);
    Free2Grid(Tg,nk,nz);
    Free2Grid(Tur,nk,nz);
    Free2Grid(Tcb,nk,nz);
    Free2Grid(Tnu1,nk,nz);
    Free2Grid(Tnu2,nk,nz);
    Free2Grid(Tnu3,nk,nz);
    free(y2d);
    free(y2d_g);
    free(y2d_ur);
    free(y2d_cb);
    free(y2d_nu1);
    free(y2d_nu2);
    free(y2d_nu3);
    cosmo->data.log_pkz = psp;
    cosmo->data.transfer_g = tsp_g;
    cosmo->data.transfer_ur = tsp_ur;
    cosmo->data.transfer_cb = tsp_cb;
    cosmo->data.transfer_nu1 = tsp_nu1;
    cosmo->data.transfer_nu2 = tsp_nu2;
    cosmo->data.transfer_nu3 = tsp_nu3;
    
    //这里不能clucu_f2d_t_free，因为cosmo->data.log_pkz是指针，指向我们创建的psp

    cosmo->computed_power = true;
}

/*-----------------------------------------------------------------------------*/
/*                                       total                                 */
/*-----------------------------------------------------------------------------*/
void clucu_compute_power(clucu_cosmology *cosmo, int *status)
{
    switch(cosmo->config.matter_power_spectrum_method)
    {
        case clucu_boltzmann_class_Pk:clucu_compute_power_class_Pk(cosmo,status);break;
        case clucu_boltzmann_class_Tk:clucu_compute_power_class_Tk(cosmo,status);break;
        case clucu_bbks:clucu_compute_linpower_bbks(cosmo,status);break;
        case clucu_eisenstein_hu:clucu_compute_linpower_eh(cosmo,1,status);break;
    }
}

double clucu_power(clucu_cosmology *cosmo,double k, double z,int *status)
{
    double pkz = clucu_f2d_t_eval(cosmo->data.log_pkz, log(k), z,
                           cosmo, status);
    return pkz;
}
double clucu_transfer_label(clucu_cosmology *cosmo,double k, double z,clucu_species_x_label label)
{
    int status=0;
    double tk;
    switch(label)
    {
        case clucu_species_g_label: tk = clucu_f2d_t_eval(cosmo->data.transfer_g, log(k), z,cosmo, &status);break;
        case clucu_species_ur_label: tk = clucu_f2d_t_eval(cosmo->data.transfer_ur, log(k), z,cosmo, &status);break;
        case clucu_species_cb_label: tk = -clucu_f2d_t_eval(cosmo->data.transfer_cb, log(k), z,cosmo, &status);break;
        case clucu_species_nu1_label: tk = clucu_f2d_t_eval(cosmo->data.transfer_nu1, log(k), z,cosmo, &status);break;
        case clucu_species_nu2_label: tk = clucu_f2d_t_eval(cosmo->data.transfer_nu2, log(k), z,cosmo, &status);break;
        case clucu_species_nu3_label: tk = clucu_f2d_t_eval(cosmo->data.transfer_nu3, log(k), z,cosmo, &status);break;
    }
    return tk;
}
double clucu_transfer(clucu_cosmology *cosmo,double k,int *status)
{
    if(!cosmo->computed_transfer){
        clucu_eh_struct_new(cosmo, 1);
        cosmo->computed_transfer = true;
    }
    double Tk2 = tsqr_EH(cosmo,cosmo->eh,k);

    return sqrt(Tk2);

}


/*-----------------------------------------------------------------------------*/
/*                                       sigma                                  */
/*-----------------------------------------------------------------------------*/

//质量的单位是：M_sun；半径的单位：Mpc，不涉及质量定义
double M_to_R(clucu_cosmology *cosmo,double halomass,int *status)
{ 
    double rho, smooth_radius;
    // is not Comoving, matter density
    rho = clucu_rho_x(cosmo, 0., clucu_species_cb_label, 0, status);
    smooth_radius = pow((3.0*halomass) / (4*M_PI*rho), (1.0/3.0));
	return smooth_radius;
}
/*----------------------*/
/*        窗口函数       */
/*----------------------*/
//窗口函数,tophat
static double W_kR(double k, double R)
{
    double kR=k*R;
    double kR2=kR*kR;
	double w;
    // This is the Maclaurin expansion of W(x)=[sin(x)-xcos(x)]*3/x**3 to O(x^10), with x=kR.
    // Necessary numerically because at low x W(x) relies on the fine cancellation of two terms
    if(kR<0.1) {
        w= 1. + kR2*(-1.0/10.0 + kR2*(1.0/280.0 +
            kR2*(-1.0/15120.0 + kR2*(1.0/1330560.0 +
            kR2* (-1.0/172972800.0)))));
        }
    else
        w = 3.*(sin(kR) - kR*cos(kR))/(kR2*kR);

	return w;
}
/*----------------------*/
/*         sigma        */
/*----------------------*/
// Params for sigma(R) integrand
typedef struct {
  clucu_cosmology *cosmo;
  double R;
  double z;
  clucu_f2d_t *psp;
  int *status;
} SigmaR_pars;

//对log10k积分，k,R单位都 不 带h
// Integrand for sigmaR integral
static double sigmaR_integrand(double log10k, void *params)
{
    double k=pow(10.,log10k);
    SigmaR_pars *par=(SigmaR_pars *)params;
    double pk=clucu_f2d_t_eval(par->psp, log10k * M_LN10, par->z,
                           par->cosmo, par->status);
    double wk=W_kR(k, par->R);
    return pk * wk*wk *k*k *k*M_LN10;//对log10k积分
}
//sigma，R单位 不 带h
double clucu_sigmaR(clucu_cosmology *cosmo,double R, double z,clucu_f2d_t *psp, int *status)
{
    SigmaR_pars par;
    par.cosmo=cosmo;
    par.status=status;
    par.psp=psp;
    par.R=R;
    par.z=z;

    gsl_integration_cquad_workspace *workspace =  NULL;
    gsl_function F;
    F.function=&sigmaR_integrand;
    F.params=&par;
    double result;

    workspace = gsl_integration_cquad_workspace_alloc(cosmo->gsl_param.N_ITERATION);
    if (workspace == NULL) {
    *status = CLUCU_ERROR_MEMORY;
    }

    if (*status == 0) 
    {
        int gslstatus = gsl_integration_cquad(&F,
                                                log10(cosmo->spline_param.K_MIN),
                                                log10(cosmo->spline_param.K_MAX),
                                                0.0, cosmo->gsl_param.INTEGRATION_SIGMAR_EPSREL,
                                                workspace,&result,NULL,NULL);
        if(gslstatus != GSL_SUCCESS) {
            CLUCU_RAISE_GSL_WARNING(gslstatus,NULL);
            *status |= gslstatus;
        }
    }
    gsl_integration_cquad_workspace_free(workspace);
    return sqrt(result/(2*M_PI*M_PI));
}
double clucu_sigma8(clucu_cosmology *cosmo, clucu_f2d_t *psp, int *status)
{
  return clucu_sigmaR(cosmo, 8/cosmo->param.h, 0., psp, status);
}

//算一系列sigma(M,A)，存起来插值
void clucu_compute_logsigma(clucu_cosmology *cosmo,int *status)
{
  if(cosmo->computed_sigma)
    return;

    int nz = cosmo->data.spline_Nz;
    double *z = cosmo->data.spline_z;
    int nm = cosmo->spline_param.LOG10M_SPLINE_NM;
    double *m = clucu_linear_spacing(cosmo->spline_param.LOG10M_SPLINE_MIN,
                          cosmo->spline_param.LOG10M_SPLINE_MAX, nm);

  double *y = NULL;

  // create space for y, to be filled with sigma
  if (*status == 0) 
  {
    y = malloc(sizeof(double)*nm*nz);
    if (y == NULL) 
    {
      *status = CLUCU_ERROR_MEMORY;
      CLUCU_RAISE_WARNING(*status,NULL);
      clucu_cosmology_set_status_message(cosmo,
                                        "clucu_massfunc.c: clucu_compute_logsigma(): "
                                        "memory allocation\n");
    }
  }

  // fill in sigma, if no errors have been triggered at this time.
  if (*status == 0) 
  {
    int thread_num = (int)(omp_get_num_procs()*THREAD_RATE/cosmo->thread_num);
	if(CHECK) thread_num=1;
    omp_set_num_threads(thread_num);
    #pragma omp parallel default(none) shared(nz, z, nm, m, y, cosmo,thread_num,status) 
    {
      double redshift, smooth_radius,sigma,log10m;
      int iz,im;
      int thread_id = omp_get_thread_num();
      for(int number = thread_id;number<nz*nm;number+=thread_num)
      {
        iz = number/nm;
        im = number % nm;
        redshift = z[iz];
        log10m = m[im];
        smooth_radius = M_to_R( cosmo,pow(10,m[im]),status);
        sigma=clucu_sigmaR(cosmo,smooth_radius, redshift,cosmo->data.log_pkz,status);
        y[iz*nm + im] = log(sigma);
        printf("iz=%d,z=%.3f,sigma=%.3f\n",iz,z,sigma);

      }//end omp for
    }//end omp parallel
  }
  for(int i=0;i<nz;i++)
        printf("[%d]%.3f\n",i,z[i]);

  gsl_spline2d *lsM = NULL;
  if (*status == 0) 
  {
    lsM = gsl_spline2d_alloc(gsl_interp2d_bicubic, nm, nz);
    if (lsM == NULL) 
    {
      *status = CLUCU_ERROR_MEMORY;
      CLUCU_RAISE_WARNING(*status,"error allocating 2D spline");
      clucu_cosmology_set_status_message(cosmo,
                                        "clucu_massfunc.c: clucu_compute_logsigma(): "
                                        "error allocating 2D spline\n");
    }
  }

  if(*status== 0) 
  {
    int s2dstatus=gsl_spline2d_init(lsM, m, z, y, nm, nz);
    if (s2dstatus) 
    {
      *status = CLUCU_ERROR_SPLINE;
      CLUCU_RAISE_WARNING(*status,"error initializing spline");
      clucu_cosmology_set_status_message(cosmo,
                                      "clucu_massfunc.c: clucu_compute_logsigma(): "
                                      "error initializing spline\n");
    }
  }

  if (*status == 0) 
  {
    cosmo->computed_sigma = true;
    cosmo->data.spline_sigma = lsM;
  }
  else
    gsl_spline2d_free(lsM);

  free(m);
  free(y);
}

/*----- ROUTINE: clucu_sigma_MZ -----
INPUT: clucu_cosmology * cosmo, double halo mass in units of Msun, double scale factor
TASK: returns sigma from the sigmaM interpolation. Also computes the sigma interpolation if
necessary.
*/
double clucu_sigma(clucu_cosmology *cosmo, double halo_mass, double redshift,int *status)
{
    double log10_halomass=log10(halo_mass);
    // Check if sigma has already been calculated
    if (!cosmo->computed_sigma) {
        *status = CLUCU_ERROR_SIGMA_INIT;
        CLUCU_RAISE_WARNING(*status, "sigma(M) spline has not been computed!");
        clucu_cosmology_set_status_message(cosmo,
                                        "clucu_massfunc.c: clucu_sigma(): "
                                        "sigma(M) spline has not been computed!");
        return NAN;
    }

    double lnsigma;
    int gslstatus = gsl_spline2d_eval_e(cosmo->data.spline_sigma, log10_halomass,
                                        redshift, NULL, NULL, &lnsigma);

    if(gslstatus != GSL_SUCCESS) {
        CLUCU_RAISE_GSL_WARNING(gslstatus,"mass=%.2e,z=%.2e is outside interpolation range [%.2e,%.2e][%.2e,%.2e]",
                                        halo_mass,redshift,
                                        pow(10.,cosmo->data.spline_sigma->interp_object.xmin),
                                        pow(10.,cosmo->data.spline_sigma->interp_object.xmax),
                                        cosmo->data.spline_sigma->interp_object.ymin,
                                        cosmo->data.spline_sigma->interp_object.ymax);
        *status|= gslstatus;
    }

    return exp(lnsigma);
}

/*----- ROUTINE: clucu_dlnsigM_dLOG10M -----
INPUT: clucu_cosmology *cosmo, double halo mass in units of Msun
TASK: returns the value of the derivative of ln(sigma^-1) with respect to log10 in halo mass.
*/
double clucu_dlnsigma_dlnm(clucu_cosmology *cosmo, double halo_mass, double redshift,int *status)
{
    double log10_halomass=log10(halo_mass);
    // Check if sigma has already been calculated
    if (!cosmo->computed_sigma) {
    *status = CLUCU_ERROR_SIGMA_INIT;
    CLUCU_RAISE_WARNING(*status,"sigma(M) spline has not been computed!");
    clucu_cosmology_set_status_message(cosmo,
                                        "clucu_massfunc.c: clucu_dlnsigma_dlnm(): "
                                        "sigma(M) spline has not been computed!");
    return NAN;
    }

    double dlnsigma_dlog10M;
    int gslstatus = gsl_spline2d_eval_deriv_x_e(cosmo->data.spline_sigma,
                                                log10_halomass, redshift,
                                                NULL, NULL, &dlnsigma_dlog10M);
    if(gslstatus){ 
        CLUCU_RAISE_GSL_WARNING(gslstatus,"mass=%.2e,z=%.2e is outside interpolation range [%.2e,%.2e][%.2e,%.2e]",
                                        halo_mass,redshift,
                                        pow(10.,cosmo->data.spline_sigma->interp_object.xmin),
                                        pow(10.,cosmo->data.spline_sigma->interp_object.xmax),
                                        cosmo->data.spline_sigma->interp_object.ymin,
                                        cosmo->data.spline_sigma->interp_object.ymax);
        *status |= gslstatus;
    }
    return -dlnsigma_dlog10M/M_LN10;
}