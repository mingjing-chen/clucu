#include "clucu_covariance.h"

typedef struct {
    clucu_cosmology *cosmo;
    double z;
    double chi;
    int *status;

    double ell;
    double ln_mass1;
    double ln_mass2;
    double ln_mass11;
    double ln_mass22;
    double ln_mass_ob;
  
} cl_pars;
//Sibb'
double SurveyWindow(clucu_cosmology *cosmo,double ell)
{
    double Theta=sqrt(cosmo->survey.Delta_Omega / M_PI);
    double x = Theta*ell;
    return 2. * gsl_sf_bessel_J1(x)/x;
}
double SamplingVariance_integrand_log10ell(double log10ell,void *params)
{
    //cosmo,z,chi
    cl_pars *par=(cl_pars *)params;
    double ell = pow(10.,log10ell);
    double ws = SurveyWindow(par->cosmo,ell);
    double k = ell/par->chi;//单位 1/Mpc
    double pk =  clucu_power(par->cosmo,k,par->z,par->status);//单位Mpc^3
    if(*(par->status)){CLUCU_RAISE_WARNING(*(par->status),NULL);abort();}
    
    return ws*ws * pk * ell * ell*M_LN10 / (2.*M_PI); //Mpc^3
}
double SamplingVariance_integrand_chi(double chi,void *params)
{
    //cosmo,ln_mass1,ln_mass2,ln_mass11,ln_mass22,
    cl_pars *par=(cl_pars *)params;
    if(chi < 1.0e-15) 
        return 0;
    double z = clucu_chi_to_z(par->cosmo,chi,par->status);

    double wh = weight_h( par->cosmo,par->ln_mass1,par->ln_mass2, z,par->status);//还没除barn,除之后单位应该为单位1/Mpc
    double whh;
    if( (fabs(par->ln_mass11-par->ln_mass1)<1.0e-10) && (fabs(par->ln_mass22-par->ln_mass2)<1.0e-10) )
        whh = wh;
    else
        whh=weight_h( par->cosmo,par->ln_mass11,par->ln_mass22,z,par->status);//还没除barn,除之后单位应该为单位1/Mpc

    cl_pars par_in;
    par_in.cosmo = par->cosmo;
    par_in.z = z;
    par_in.chi = chi;
    par_in.status = par->status;

    gsl_integration_cquad_workspace *workspace =  NULL;
    workspace = gsl_integration_cquad_workspace_alloc(par->cosmo->gsl_param.N_ITERATION);
    gsl_function F;
    F.function=&SamplingVariance_integrand_log10ell;
    F.params=&par_in;
    //确定积分的上下限
	double x_min,x_max,y;
    x_min=-3.;
    x_max=3.;
    int gslstatus = gsl_integration_cquad(&F,
											x_min,
											x_max,
											0.0, par->cosmo->gsl_param.INTEGRATION_MASS_TRUE_EPSREL,
											workspace,&y,NULL,NULL);
	gsl_integration_cquad_workspace_free(workspace);

    return wh*whh * y /(chi*chi); //还没除barn,除之后单位应该为单位1/Mpc
}
//先对ell积分，再对chi积分
double SamplingVariance(clucu_cosmology *cosmo,double ln_mass1,double ln_mass2,double ln_mass11,double ln_mass22,double chi1,double chi2,int *status)
{
    cl_pars par;
    par.cosmo = cosmo;
    par.ln_mass1 = ln_mass1;
    par.ln_mass2 = ln_mass2;
    par.ln_mass11 = ln_mass11;
    par.ln_mass22 = ln_mass22;
    par.status = status;


    gsl_integration_cquad_workspace *workspace =  NULL;
    workspace = gsl_integration_cquad_workspace_alloc(par.cosmo->gsl_param.N_ITERATION);
    gsl_function F;
    F.function=&SamplingVariance_integrand_chi;
    F.params=&par;
    //确定积分的上下限
	double x_min,x_max,y;
    x_min=chi1;
    x_max=chi2;
    int gslstatus = gsl_integration_cquad(&F,
											x_min,
											x_max,
											0.0, par.cosmo->gsl_param.INTEGRATION_MASS_TRUE_EPSREL,
											workspace,&y,NULL,NULL);
	gsl_integration_cquad_workspace_free(workspace);
    return y * pow(cosmo->survey.Delta_Omega,2.);
    //最终的结果除一个n_bar又乘一个n_bar = 不用除n_bar = y
    //最终结果的单位 = n_bar的单位（红移bin内、单位立体角的星系团数目）* 立体角^2 = y * Delta_Omega^2
}
//Sibb'
void clucu_compute_samplingvariance(clucu_cosmology *cosmo,int *status)
{
    cosmo->data.CovNN_SV = Create3Grid(cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.Nzob);
    int Nm=cosmo->data.Nmob-1;
    int Nz=cosmo->data.Nzob-1;
    int thread_num = (int)( omp_get_num_procs() * THREAD_RATE/cosmo->thread_num);
	if(CHECK) thread_num=1;
	omp_set_num_threads(thread_num);
    #pragma omp parallel default(none) shared(Nm,Nz,cosmo,thread_num,status)
	{
        int local_status=*status;
		double ln_mass_ob2,ln_mass_ob1,ln_mass_ob22,ln_mass_ob11,chi1,chi2;
		int iz;
		int thread_id = omp_get_thread_num();
		for(int i = thread_id;i<Nz;i+=thread_num)
		{
            iz = i;
            double z1=cosmo->data.zob[iz];
            double z2=cosmo->data.zob[iz+1];
            chi1 = clucu_comoving_angular_distance(cosmo,cosmo->data.zob[iz],&local_status);
            chi2 = clucu_comoving_angular_distance(cosmo,cosmo->data.zob[iz+1],&local_status);
            for(int im=0;im<Nm;im++)
            {
                ln_mass_ob1 = cosmo->data.ln_mass_ob[im];
                ln_mass_ob2 = cosmo->data.ln_mass_ob[im+1];
                for(int imm=0;imm<=im;imm++)
                {
                    ln_mass_ob11 = cosmo->data.ln_mass_ob[imm];
                    ln_mass_ob22 = cosmo->data.ln_mass_ob[imm+1];
                    cosmo->data.CovNN_SV[im][imm][iz]=SamplingVariance(cosmo,ln_mass_ob1,ln_mass_ob2,ln_mass_ob11,ln_mass_ob22, chi1,chi2,&local_status);
                }
            }
            printf("cosmo_id[%d]:Covhh_SV[%d][%d][%d]=%le\n",cosmo->cosmo_id,Nm-1,Nm-1,iz,cosmo->data.CovNN_SV[Nm-1][Nm-1][iz]);
            //补全
            for(int im=0;im<Nm;im++)
            {
                for(int imm=im+1;imm<Nm;imm++)
                {
                    cosmo->data.CovNN_SV[im][imm][iz]=cosmo->data.CovNN_SV[imm][im][iz];
                }
            }
        }
        if(local_status) {
        #pragma omp atomic write
        *status=local_status;
        }
	}
    cosmo->computed_CovNN_SV=true;
}
void clucu_compute_cov_NN(clucu_cosmology *cosmo,int *status)
{
    cosmo->data.CovNN = Create4Grid(cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.Nzob,cosmo->data.Nzob);
    //Nbin_mzob
    //cosmo->data.CovNN_SV[im][imm][iz][izz]
    int Nm=cosmo->data.Nmob-1;
    int Nz=cosmo->data.Nzob-1;
    for(int im=0;im<Nm;im++)
    {
        for(int imm=0;imm<Nm;imm++)
        {
            for(int iz=0;iz<Nz;iz++)
            {
                for(int izz=0;izz<Nz;izz++)
                {
                    cosmo->data.CovNN[im][imm][iz][izz] = cosmo->data.Nbin_mzob[im][iz] * kronecker(im,imm) * kronecker(iz,izz)
                                                            + cosmo->data.CovNN_SV[im][imm][iz] * kronecker(iz,izz);   
                }
            }
        }
    }
    cosmo->computed_CovNN=true;
}
void SaveCov_NN(clucu_cosmology *cosmo)
{
	char namedir[50];
    char nameNlzob[100];
	sprintf(namedir,"../output/%s", cosmo->name);
	sprintf(nameNlzob,"../output/%s/Cov%d_NN.dat", cosmo->name,cosmo->cosmo_id);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameNlzob,"w");
	for(int iz=0;iz<cosmo->data.Nzob;iz++)
		fprintf(fp,"%f      %le\n",cosmo->data.zob[iz],cosmo->data.CovNN[0][0][iz][iz]);
	fclose(fp);
}
void clucu_compute_cov_hhhh(clucu_cosmology *cosmo,int *status)
{
    cosmo->data.Covhhhh = Create7Grid(cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.NzobP,cosmo->data.Nell);
    int Nm=cosmo->data.Nmob-1;
    int Nz=cosmo->data.NzobP-1;
    int Nell=cosmo->data.Nell-1;
    double delta_ell,ell,prefactor;
    for(int iell=0;iell<Nell;iell++)
    {
        ell = cosmo->data.ell[iell];
        delta_ell = cosmo->data.ell[iell+1]-ell;
        prefactor = (4.*M_PI/cosmo->survey.Delta_Omega) / ((2.*ell+1)*delta_ell);
        
        for(int iz=0;iz<Nz;iz++)
        {
            for(int jz=0;jz<Nz;jz++)
            {
                for(int bm=0;bm<Nm;bm++)
                {
                    for(int bam=0;bam<Nm;bam++)
                    {
                        for(int btm=0;btm<Nm;btm++)
                        {
                            for(int btam=0;btam<Nm;btam++)
                            {
                                cosmo->data.Covhhhh[bm][bam][btm][btam][iz][jz][iell] = prefactor * kronecker(iz,jz)
                                                                                    *( cosmo->data.cl_hh_hat[bm][btm][iz][iell] 
                                                                                    * cosmo->data.cl_hh_hat[bam][btam][iz][iell] 
                                                                                    + cosmo->data.cl_hh_hat[bm][btam][iz][iell] 
                                                                                    * cosmo->data.cl_hh_hat[bam][btm][iz][iell]);
                            }
                        }
                    }
                }
            }
        }
    }                   
    cosmo->computed_Covhhhh=true;
}
void SaveCov_hhhh(clucu_cosmology *cosmo)
{
	char namedir[50];
    char nameNlzob[100];
	sprintf(namedir,"../output/%s", cosmo->name);
	sprintf(nameNlzob,"../output/%s/Cov%d_hhhh.dat", cosmo->name,cosmo->cosmo_id);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameNlzob,"w");
    save_matrix(fp,cosmo->data.Covhhhh[0][0][0][0][0],cosmo->data.zobP,cosmo->data.ell,cosmo->data.NzobP,cosmo->data.Nell);
	fclose(fp);
}
void clucu_compute_cov_hkhh(clucu_cosmology *cosmo,int *status)
{
    cosmo->data.Covhkhh = Create6Grid(cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.NzobP,cosmo->data.Nell);
    int Nm=cosmo->data.Nmob-1;
    int Nz=cosmo->data.NzobP-1;
    int Nell=cosmo->data.Nell-1;
    double delta_ell,ell,prefactor;
    for(int iell=0;iell<Nell;iell++)
    {
        ell = cosmo->data.ell[iell];
        delta_ell = cosmo->data.ell[iell+1]-ell;
        prefactor = (4.*M_PI/cosmo->survey.Delta_Omega) / ((2.*ell+1)*delta_ell);
        
        for(int iz=0;iz<Nz;iz++)
        {
            for(int jz=0;jz<Nz;jz++)
            {
                for(int bm=0;bm<Nm;bm++)
                {
                    for(int bam=0;bam<Nm;bam++)
                    {
                        for(int btm=0;btm<Nm;btm++)
                        {
                            cosmo->data.Covhkhh[btm][bm][bam][iz][jz][iell] = prefactor * kronecker(iz,jz)
                                                                            *( cosmo->data.cl_hh_hat[btm][bm][iz][iell] 
                                                                            * cosmo->data.cl_hk[bam][iz][iell] 
                                                                            + cosmo->data.cl_hh_hat[btm][bam][iz][iell] 
                                                                            * cosmo->data.cl_hk[bm][iz][iell]);
                        }
                    }
                }
            }
        }
    }                   
    cosmo->computed_Covhkhh=true;
}
void clucu_compute_cov_hkhk(clucu_cosmology *cosmo,int *status)
{
    cosmo->data.Covhkhk = Create5Grid(cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.NzobP,cosmo->data.Nell);
    int Nm=cosmo->data.Nmob-1;
    int Nz=cosmo->data.NzobP-1;
    int Nell=cosmo->data.Nell-1;
    double delta_ell,ell,prefactor;
    double a,b,c,d,e,f,p;
    for(int iell=0;iell<Nell;iell++)
    {
        ell = cosmo->data.ell[iell];
        delta_ell = cosmo->data.ell[iell+1]-ell;
        prefactor = (4.*M_PI/cosmo->survey.Delta_Omega) / ((2.*ell+1)*delta_ell);
        
        for(int iz=0;iz<Nz;iz++)
        {
            for(int jz=0;jz<Nz;jz++)
            {
                for(int bm=0;bm<Nm;bm++)
                {
                    for(int bam=0;bam<Nm;bam++)
                    {
                        p=prefactor;
                        a=cosmo->data.cl_hh_hat[bm][bam][iz][iell] ;
                        b=cosmo->data.cl_kk_hat[iell] ;
                        c=cosmo->data.cl_hk[bm][iz][iell] ;
                        d=a*b;
                        e=c*c;
                        f=p*(d+e);
                        cosmo->data.Covhkhk[bm][bam][iz][jz][iell] = prefactor 
                                                                    *( cosmo->data.cl_hh_hat[bm][bam][iz][iell] 
                                                                    * cosmo->data.cl_kk_hat[iell] 
                                                                    * kronecker(iz,jz)
                                                                    + cosmo->data.cl_hk[bm][iz][iell] 
                                                                    * cosmo->data.cl_hk[bam][jz][iell]);

                    }
                }
            }
        }
    }                   
    cosmo->computed_Covhkhk=true;
}
void SaveCov_hkhk(clucu_cosmology *cosmo)
{
	char namedir[50];
    char nameNlzob[100];
	sprintf(namedir,"../output/%s", cosmo->name);
	sprintf(nameNlzob,"../output/%s/Cov%d_hkhk.dat", cosmo->name,cosmo->cosmo_id);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameNlzob,"w");
	for(int iell=0;iell<cosmo->data.Nell;iell++)
		fprintf(fp,"%f      %le\n",cosmo->data.ell[iell],cosmo->data.Covhkhk[0][0][0][0][iell]);
	fclose(fp);
}