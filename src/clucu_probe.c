/**
 * @file clucu_probe.c
 * @author Mingjing Chen (mingjing@mail.ustc.edu.cn)
 * @brief 
 * @version 0.1
 * @date 2022-09-07
 * 
 * @copyright Copyright (c) 2022
 */
#include "clucu_probe.h"


// Integrand for numbercounts integral: d ln_mass_true
double numbercounts_integrand_ln_mass_true(double ln_mass_true, void *params)
{
	int_pars *par=(int_pars *)params;
	double dn_dlnm = clucu_halo_dn_dlnm(par->cosmo,exp(ln_mass_true),par->z_true,par->status);
    double pm = mass_ob_probability(par->cosmo,par->ln_mass_ob,ln_mass_true,par->z_true);
    double v = clucu_volume_element(par->cosmo,par->z_true,par->status);
    double Omega = par->cosmo->survey.Delta_Omega;
    double pz = redshift_probability(par->cosmo,par->z_ob,par->z_true);
    return dn_dlnm * pm * v * Omega * pz;
}
void clucu_compute_cluster_numbercounts(clucu_cosmology *cosmo,int CEN_massob,int CEN_zob,int *status)
{
    cosmo->data.Nbin_mzob = Create2Grid(cosmo->data.Nmob,cosmo->data.Nzob);
	int thread_num = (int)( omp_get_num_procs() * THREAD_RATE/cosmo->thread_num);
	if(CHECK) thread_num=1;
	omp_set_num_threads(thread_num);
    #pragma omp parallel default(none) shared(cosmo,thread_num,CEN_massob,CEN_zob,status)
    {
        int local_status=*status;
        double ln_mass_ob2,ln_mass_ob1, z_ob2,z_ob1;
        int number_mob,number_zob;
        int thread_id = omp_get_thread_num();
        for(int number = thread_id;number<(cosmo->data.Nmob-1)*(cosmo->data.Nzob-1);number+=thread_num)
        {
        number_mob = number / (cosmo->data.Nzob-1);
        number_zob = number % (cosmo->data.Nzob-1);
        ln_mass_ob1 = cosmo->data.ln_mass_ob[number_mob];
        ln_mass_ob2 = cosmo->data.ln_mass_ob[number_mob+1];
        z_ob1 = cosmo->data.zob[number_zob];
        z_ob2 = cosmo->data.zob[number_zob+1];
        cosmo->data.Nbin_mzob[number_mob][number_zob] = func_integrand_mass_z(cosmo,numbercounts_integrand_ln_mass_true,1,1,ln_mass_ob1,ln_mass_ob2,CEN_massob,z_ob1,z_ob2,CEN_zob,&local_status);
        printf("cosmo_id[%d]:Nbin_mzob[%d][%d]=%f\n",cosmo->cosmo_id,number_mob,number_zob,cosmo->data.Nbin_mzob[number_mob][number_zob]);
        }//end omp for
        if(local_status) {
        #pragma omp atomic write
        *status=local_status;
        }
    }//end omp parallel
	cosmo->computed_cluster_numbercounts=true;
}


// Integrand for averagebias integral: d ln_mass_true
double average_hmf_integrand_ln_mass_true(double ln_mass_true, void *params)
{
	int_pars *par=(int_pars *)params;
	double dn_dlnm=clucu_halo_dn_dlnm(par->cosmo,exp(ln_mass_true),par->z_true,par->status);
    double pm = mass_ob_probability(par->cosmo,par->ln_mass_ob,ln_mass_true,par->z_true);
    return dn_dlnm * pm;
}
double average_hmfb_integrand_ln_mass_true(double ln_mass_true, void *params)
{
	int_pars *par=(int_pars *)params;	
	double dndlnm_bias=clucu_halo_dndlnm_bias(par->cosmo,exp(ln_mass_true),par->z_true,par->status);
    double pm = mass_ob_probability(par->cosmo,par->ln_mass_ob,ln_mass_true,par->z_true);
    return dndlnm_bias * pm;
}
double average_hmfbk_integrand_ln_mass_true(double ln_mass_true, void *params)
{
	int_pars *par=(int_pars *)params;	
	double dndlnm_bias=clucu_halo_dndlnm_bias_k(par->cosmo,exp(ln_mass_true),par->z_true,par->k,par->status);
    double pm = mass_ob_probability(par->cosmo,par->ln_mass_ob,ln_mass_true,par->z_true);
    return dndlnm_bias * pm;
}
//积分mass_true，不积分z_true，积分mass_ob，不积分z_ob
void clucu_compute_cluster_averagebias(clucu_cosmology *cosmo,int *status)
{
    cosmo->data.sumHMF = Create1Grid(cosmo->data.NzobP);
    cosmo->data.averageBias = Create1Grid(cosmo->data.NzobP);
    cosmo->data.Nbin_zobP = Create1Grid(cosmo->data.NzobP);
	double **HMF = Create2Grid(cosmo->data.Nmob-1,cosmo->data.NzobP);
	double **HMFb = Create2Grid(cosmo->data.Nmob-1,cosmo->data.NzobP);
	double **HMFN = Create2Grid(cosmo->data.Nmob-1,cosmo->data.NzobP);
	int thread_num = (int)( omp_get_num_procs() * THREAD_RATE/cosmo->thread_num);
	if(CHECK) thread_num=1;
	omp_set_num_threads(thread_num);
    #pragma omp parallel default(none) shared(cosmo,thread_num,HMF,HMFb,HMFN,status)
	{
        int local_status=*status;
		double ln_mass_ob2,ln_mass_ob1, zobP;
		int number_mob,number_zobP;
		int thread_id = omp_get_thread_num();
		for(int number = thread_id;number<(cosmo->data.Nmob-1)*cosmo->data.NzobP;number+=thread_num)
		{
			number_mob = number / cosmo->data.NzobP;
			number_zobP = number % cosmo->data.NzobP;
            ln_mass_ob1 = cosmo->data.ln_mass_ob[number_mob];
            ln_mass_ob2 = cosmo->data.ln_mass_ob[number_mob+1];
			zobP = cosmo->data.zobP[number_zobP];
			HMF[number_mob][number_zobP]=func_integrand_mass_z(cosmo,average_hmf_integrand_ln_mass_true,1,0,ln_mass_ob1,ln_mass_ob2,1,zobP,zobP,-1,&local_status);
			HMFb[number_mob][number_zobP]=func_integrand_mass_z(cosmo,average_hmfb_integrand_ln_mass_true,1,0,ln_mass_ob1,ln_mass_ob2,1,zobP,zobP,-1,&local_status);
			HMFN[number_mob][number_zobP]=HMF[number_mob][number_zobP]*clucu_volume_element(cosmo,zobP,&local_status) * cosmo->survey.Delta_Omega;
		}
        if(local_status) {
        #pragma omp atomic write
        *status=local_status;
        }
	}
	for(int number_zobP=0; number_zobP<cosmo->data.NzobP; number_zobP++)
	{
		double sumHMF=0.,sumHMFb=0.,sumHMFN=0.;
		for(int number_mob=0;number_mob<cosmo->data.Nmob-1;number_mob++)
		{
			sumHMF+=HMF[number_mob][number_zobP];
			sumHMFb+=HMFb[number_mob][number_zobP];
			sumHMFN+=HMFN[number_mob][number_zobP];
		}
        
		cosmo->data.averageBias[number_zobP]=sumHMFb/sumHMF;
		cosmo->data.sumHMF[number_zobP] = sumHMF;
        //严格来说，这里少了对V(z)dz的积分
		cosmo->data.Nbin_zobP[number_zobP] = sumHMFN * cosmo->survey.zobP_bin;
        printf("cosmo_id[%d]:cosmo->sumHMF[%d]=%le,Nbin_zobP[%d]=%f\n",cosmo->cosmo_id,number_zobP,cosmo->data.sumHMF[number_zobP],number_zobP,cosmo->data.Nbin_zobP[number_zobP]);
	}
	Free2Grid(HMF,cosmo->data.Nmob-1,cosmo->data.NzobP);
	Free2Grid(HMFb,cosmo->data.Nmob-1,cosmo->data.NzobP);
	Free2Grid(HMFN,cosmo->data.Nmob-1,cosmo->data.NzobP);


	cosmo->computed_cluster_averagebias=true;

}
//如果有GISDB，则beff(k,z)，还有k
void clucu_compute_cluster_averagebias_k(clucu_cosmology *cosmo,int *status)
{
    cosmo->data.sumHMF = Create1Grid(cosmo->data.NzobP);
    cosmo->data.averageBias = Create1Grid(cosmo->data.NzobP);
    cosmo->data.Nbin_zobP = Create1Grid(cosmo->data.NzobP);
    //k, 0-0.3差值
    int Nk = 121;
    double *k = clucu_linear_spacing(0.,0.3,Nk);
    int Nz = cosmo->data.NzobP;
	double **HMF = Create2Grid(cosmo->data.Nmob-1,Nz);
	double ***HMFb = Create3Grid(cosmo->data.Nmob-1,Nz,Nk);
	double **HMFN = Create2Grid(cosmo->data.Nmob-1,Nz);
    double *bias = Create1Grid(Nk*Nz);
	int thread_num = (int)( omp_get_num_procs() * THREAD_RATE/cosmo->thread_num);
	if(CHECK) thread_num=1;
	omp_set_num_threads(thread_num);
    #pragma omp parallel default(none) shared(cosmo,thread_num,HMF,HMFb,HMFN,status,Nk,k)
	{
        int local_status=*status;
		double ln_mass_ob2,ln_mass_ob1, zobP;
		int number_mob,number_zobP;
		int thread_id = omp_get_thread_num();
		for(int number = thread_id;number<(cosmo->data.Nmob-1)*cosmo->data.NzobP;number+=thread_num)
		{
			number_mob = number / cosmo->data.NzobP;
			number_zobP = number % cosmo->data.NzobP;
            ln_mass_ob1 = cosmo->data.ln_mass_ob[number_mob];
            ln_mass_ob2 = cosmo->data.ln_mass_ob[number_mob+1];
			zobP = cosmo->data.zobP[number_zobP];
			HMF[number_mob][number_zobP]=func_integrand_mass_z(cosmo,average_hmf_integrand_ln_mass_true,1,0,ln_mass_ob1,ln_mass_ob2,1,zobP,zobP,-1,&local_status);
			HMFN[number_mob][number_zobP]=HMF[number_mob][number_zobP]*clucu_volume_element(cosmo,zobP,&local_status) * cosmo->survey.Delta_Omega;
            for(int ik=0;ik<Nk;ik++){
                HMFb[number_mob][number_zobP][ik]=func_integrand_mass_z_k(cosmo,average_hmfbk_integrand_ln_mass_true,1,0,k[ik],ln_mass_ob1,ln_mass_ob2,1,zobP,zobP,-1,&local_status);
            }
        }
        if(local_status) {
        #pragma omp atomic write
        *status=local_status;
        }
	}
    for(int ik=0;ik<Nk;ik++)
    {
        for(int number_zobP=0; number_zobP<cosmo->data.NzobP; number_zobP++)
        {
            double sumHMF=0.,sumHMFb=0.,sumHMFN=0.;
            for(int number_mob=0;number_mob<cosmo->data.Nmob-1;number_mob++)
            {
                sumHMF+=HMF[number_mob][number_zobP];
                sumHMFb+=HMFb[number_mob][number_zobP][ik];
                sumHMFN+=HMFN[number_mob][number_zobP];
            }
            bias[number_zobP*Nk + ik] = sumHMFb/sumHMF;            
            cosmo->data.sumHMF[number_zobP] = sumHMF;
            //严格来说，这里少了对V(z)dz的积分
            cosmo->data.Nbin_zobP[number_zobP] = sumHMFN * cosmo->survey.zobP_bin;
            printf("cosmo_id[%d]:cosmo->sumHMF[%d]=%le,Nbin_zobP[%d]=%f\n",cosmo->cosmo_id,number_zobP,cosmo->data.sumHMF[number_zobP],number_zobP,cosmo->data.Nbin_zobP[number_zobP]);
        }
    }
    //对beff(k,z)二维插值
    gsl_spline2d *bkz_eff = NULL;
    if (*status == 0) 
    {
    bkz_eff = gsl_spline2d_alloc(gsl_interp2d_bicubic, Nk, Nz);
    if (bkz_eff == NULL) 
    {
        *status = CLUCU_ERROR_MEMORY;
        CLUCU_RAISE_WARNING(*status,"error allocating 2D spline");
    }
    }

    if(*status== 0) 
    {
    int s2dstatus=gsl_spline2d_init(bkz_eff, k, cosmo->data.zobP, bias, Nk, Nz);
    if (s2dstatus) 
    {
        *status = CLUCU_ERROR_SPLINE;
        CLUCU_RAISE_WARNING(*status,"error initializing spline");
    }
    }

    if (*status == 0) 
    {
    cosmo->computed_cluster_averagebias = true;
    cosmo->data.spline_averageBias_k = bkz_eff;
    }
    else
    gsl_spline2d_free(bkz_eff);

    free(k);
    free(bias);

	Free2Grid(HMF,cosmo->data.Nmob-1,cosmo->data.NzobP);
	Free3Grid(HMFb,cosmo->data.Nmob-1,cosmo->data.NzobP,Nk);
	Free2Grid(HMFN,cosmo->data.Nmob-1,cosmo->data.NzobP);
}
double clucu_averageBias_k(clucu_cosmology *cosmo, double k, double z,int *status)
{
    // Check if sigma has already been calculated
    if (!cosmo->computed_cluster_averagebias) {
        *status = CLUCU_ERROR_SIGMA_INIT;
        CLUCU_RAISE_WARNING(*status, "averagebias spline has not been computed!");
        return NAN;
    }

    double result;
    int gslstatus = gsl_spline2d_eval_e(cosmo->data.spline_averageBias_k, k,
                                        z, NULL, NULL, &result);

    if(gslstatus != GSL_SUCCESS) {
        CLUCU_RAISE_GSL_WARNING(gslstatus,"k=%.2e,z=%.2e is outside interpolation range [%.2e,%.2e][%.2e,%.2e]",
                                        k,z,
                                        cosmo->data.spline_averageBias_k->interp_object.xmin,
                                        cosmo->data.spline_averageBias_k->interp_object.xmax,
                                        cosmo->data.spline_averageBias_k->interp_object.ymin,
                                        cosmo->data.spline_averageBias_k->interp_object.ymax);
        *status|= gslstatus;
    }

    return result;
}



void clucu_compute_cluster_powerspectrum(clucu_cosmology *cosmo,int *status)
{
    cosmo->data.beta = Create1Grid(cosmo->data.Nzob);
    cosmo->data.bias = Create3Grid(cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);
    cosmo->data.pc = Create3Grid(cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);
    cosmo->data.pc0 = Create3Grid(cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);
    cosmo->data.lnkPc = Create3Grid(cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);
    if (!cosmo->computed_cluster_averagebias)
    {
        *status = CLUCU_ERROR_COMPUTE;
        CLUCU_RAISE_WARNING(*status,"averagebias has not been computed!");
        cosmo->status = CLUCU_ERROR_SIGMA_INIT;
        clucu_cosmology_set_status_message(cosmo,
                                        "clucu_probe.c: clucu_compute_cluster_volumeeffect(): "
                                        "averagebias has not been computed!");
    }
	for (int number_zobP=0;number_zobP<cosmo->data.NzobP;number_zobP++)
    {
        double zobP = cosmo->data.zobP[number_zobP];

        double s_perp = cosmo->data.shiftperp[number_zobP];
        double s_para = cosmo->data.shiftpara[number_zobP];
        double f=clucu_growth_rate(cosmo,zobP,status);
        double Hz = clucu_hubble_parameter(cosmo,zobP,status);//插值表的单位是1/Mpc，实际上是Hz/c
        double sigma = cosmo->survey.redshift_sigma0/Hz*(1+zobP);//单位Mpc
        for (int number_k1=0;number_k1<cosmo->data.Nk1;number_k1++)//k1循环
        {
            for (int number_k2=0;number_k2<cosmo->data.Nk2;number_k2++)//k2循环
            {
                //k_fid单位1/Mpc
                //给定的是k_fiducial，单位1/Mpc
                double k1_fid = cosmo->data.k1[number_k1];
                double k2_fid = cosmo->data.k2[number_k2];
                double k_fid = sqrt( k1_fid*k1_fid + k2_fid*k2_fid);
                double k1_true = k1_fid * s_perp ;
                double k2_true = k2_fid * s_para ;
                double k_true = sqrt( k1_true*k1_true + k2_true*k2_true);// 1/Mpc
                double bz;
                if(cosmo->GISDB==0)
                    bz=cosmo->data.averageBias[number_zobP];
                else if(cosmo->GISDB==1)
                    bz=clucu_averageBias_k(cosmo,k_true,zobP,status);
//TODO：这里的k要分方向吗？即，需要考虑GISDB效应对RSD的影响吗
//                bz=(bz-1.)*GISDB_factor(cosmo,k_true,zobP)+1.;
//这里得先对b^L+1积分！不能直接乘G(k,z)
                double beta=f/bz;
                cosmo->data.beta[number_zobP]=beta;
                double bias = ( 1.0+beta*pow(k2_true/k_true,2) ) * bz;
                cosmo->data.bias[number_k1][number_k2][number_zobP] = bias;

                double pkz = clucu_power(cosmo,k_true,zobP,status);//单位均为Mpc
                double pc0 = pow(bias,2) * pkz;
                double Pc_e = exp(   -pow(k2_true*sigma,2.0)    );
                double pc = pc0 * Pc_e;
                cosmo->data.pc0[number_k1][number_k2][number_zobP] = pc0;
                cosmo->data.pc[number_k1][number_k2][number_zobP] = pc;
                //double lnkpc = log( pow(k1_fid,2) * k2_fid * pc);
                //double lnkpc = log( pc);
                double lnkpc = log( pow(k1_true,2) * k2_true * pc);
                cosmo->data.lnkPc[number_k1][number_k2][number_zobP] = lnkpc;
                //printf("lnkpc[%d][%d][%d]=%le\n",number_zobP,number_k1,number_k2,lnkpc);
            }
        }
    }
    cosmo->computed_cluster_clusterpowerspectrum=true;
}
void clucu_compute_cluster_volumeeffect(clucu_cosmology *cosmo,int *status)
{
    cosmo->data.V_eff = Create3Grid(cosmo->data.Nk1,cosmo->data.Nk2,cosmo->data.NzobP);
    cosmo->data.V_k = Create2Grid(cosmo->data.Nk1,cosmo->data.Nk2);
    if (!cosmo->computed_cluster_averagebias)
    {
         *status = CLUCU_ERROR_COMPUTE;
        CLUCU_RAISE_WARNING(*status,"averagebias has not been computed!");
        cosmo->status = CLUCU_ERROR_SIGMA_INIT;
        clucu_cosmology_set_status_message(cosmo,
                                        "clucu_probe.c: clucu_compute_cluster_volumeeffect(): "
                                        "averagebias has not been computed!");
    }
    if (!cosmo->computed_cluster_clusterpowerspectrum)
    {
         *status = CLUCU_ERROR_COMPUTE;
        CLUCU_RAISE_WARNING(*status,"cluster power spectrum has not been computed!");
        cosmo->status = CLUCU_ERROR_SIGMA_INIT;
        clucu_cosmology_set_status_message(cosmo,
                                        "clucu_probe.c: clucu_compute_cluster_volumeeffect(): "
                                        "cluster power spectrum has not been computed!");
    }
	double hmf,pc,v;
    for (int number_k1=0;number_k1<cosmo->data.Nk1;number_k1++)//k1循环
    {
        for (int number_k2=0;number_k2<cosmo->data.Nk2;number_k2++)//k2循环
        {
            
            for(int number_zobP=0;number_zobP<cosmo->data.NzobP;number_zobP++)
            {
                hmf = cosmo->data.sumHMF[number_zobP];
                v = clucu_volume_element(cosmo,cosmo->data.zobP[number_zobP],status);
                pc = cosmo->data.pc[number_k1][number_k2][number_zobP];
                cosmo->data.V_eff[number_k1][number_k2][number_zobP] = pow(  hmf*pc / (1.0+hmf*pc) ,2.0  ) * v * cosmo->survey.Delta_Omega * cosmo->survey.zobP_bin;
            }
        }
    }

    double k1,k11;
    for (int number_k1=0;number_k1<cosmo->data.Nk1;number_k1++)//k1循环
    {
        k1 = cosmo->data.k1[number_k1] - 0.5*cosmo->survey.k1_bin;
        k11 = k1 + cosmo->survey.k1_bin;
        for (int number_k2=0;number_k2<cosmo->data.Nk2;number_k2++)//k2循环
        {
            cosmo->data.V_k[number_k1][number_k2]=(k11*k11-k1*k1) * cosmo->survey.k2_bin / pow(2*M_PI,2);
        }
    }
	cosmo->computed_cluster_volumeeffect=true;
}