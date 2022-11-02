#include "clucu_cls.h"

/*----------------------*/
/*        two kernel       */
/*----------------------*/
//调用probe.c中的函数算kernel
//double average_hmf_integrand(clucu_cosmology *cosmo,double ln_mass_ob,double z_true)
//double average_hmfb_integrand(clucu_cosmology *cosmo,double ln_mass_ob,double z_true)
//只对Mtrue,ztrue积分，得到结果后对Mob,zob插值，得到kernel(Mob,zob)
void clucu_compute_kernel_n(clucu_cosmology *cosmo,int *status)
{
    //按zob插值，最终结果是zobP
    int Nz=100;
    double *z=clucu_linear_spacing(cosmo->survey.zob_min,cosmo->survey.zob_max,Nz);
    int Nm=100;
    double *lnm=clucu_linear_spacing(log(cosmo->survey.mass_ob_min/10.),log(cosmo->survey.mass_ob_max*5.),Nm);
    double *y=malloc(sizeof(double)*Nm*Nz);

    int thread_num = (int)(omp_get_num_procs()*THREAD_RATE/cosmo->thread_num);
	if(CHECK) thread_num=1;
    omp_set_num_threads(thread_num);
    #pragma omp parallel default(none) shared(Nz, z, Nm, lnm, y, cosmo,thread_num,status)
    {
        int local_status=*status;
        double redshift, ln_mass_ob;
        int number_z,number_m;
        int thread_id = omp_get_thread_num();
        for(int number = thread_id;number<Nz*Nm;number+=thread_num)
        {
            number_z = number/Nm;
            number_m = number % Nm;
            redshift = z[number_z];
            ln_mass_ob = lnm[number_m];
            y[number_z*Nm + number_m] = func_integrand_mass_z(cosmo,average_hmf_integrand_ln_mass_true,1,1,ln_mass_ob,ln_mass_ob,-1,redshift,redshift,-1,&local_status);


        }//end omp for
        if(local_status) {
            #pragma omp atomic write
            *status=local_status;
            }
    }//end omp parallel
    gsl_spline2d *kernel = NULL;
    kernel = gsl_spline2d_alloc(gsl_interp2d_bicubic, Nm, Nz);
    gsl_spline2d_init(kernel, lnm, z, y, Nm, Nz);
    cosmo->computed_kernel_n = true;
    cosmo->data.spline_kernel_n = kernel;

    free(z);
    free(lnm);
    free(y);
}

void clucu_compute_kernel_nb(clucu_cosmology *cosmo,int *status)
{
    //按zob插值，最终结果是zobP
    int Nz=100;
    double *z=clucu_linear_spacing(cosmo->survey.zob_min,cosmo->survey.zob_max,Nz);
    int Nm=100;
    double *lnm=clucu_linear_spacing(log(cosmo->survey.mass_ob_min/10.),log(cosmo->survey.mass_ob_max*5.),Nm);
    double *y=malloc(sizeof(double)*Nm*Nz);

    int thread_num = (int)(omp_get_num_procs()*THREAD_RATE/cosmo->thread_num);
	if(CHECK) thread_num=1;
    omp_set_num_threads(thread_num);
    #pragma omp parallel default(none) shared(Nz, z, Nm, lnm, y, cosmo,thread_num,status)
    {
        int local_status=*status;
        double redshift, ln_mass_ob;
        int number_z,number_m;
        int thread_id = omp_get_thread_num();
        for(int number = thread_id;number<Nz*Nm;number+=thread_num)
        {
        number_z = number/Nm;
        number_m = number % Nm;
        redshift = z[number_z];
        ln_mass_ob = lnm[number_m];
        y[number_z*Nm + number_m] = func_integrand_mass_z(cosmo,average_hmfb_integrand_ln_mass_true,1,1,ln_mass_ob,ln_mass_ob,-1,redshift,redshift,-1,&local_status);

        }//end omp for
        if(local_status) {
        #pragma omp atomic write
        *status=local_status;
        }
    }//end omp parallel
    gsl_spline2d *kernel = NULL;
    kernel = gsl_spline2d_alloc(gsl_interp2d_bicubic, Nm, Nz);
    gsl_spline2d_init(kernel, lnm, z, y, Nm, Nz);
    cosmo->computed_kernel_nb = true;
    cosmo->data.spline_kernel_nb = kernel;

    free(z);
    free(lnm);
    free(y);
}
double kernel_n(clucu_cosmology *cosmo,double ln_mass_ob,double z_ob,int *status)
{
    double result;
    if(cosmo->data.spline_kernel_nb==NULL){
        *status = 243234;
        CLUCU_RAISE_WARNING(*status,"you shoule set kernel_n");
    }
    if(fabs(z_ob - cosmo->data.spline_kernel_n->interp_object.ymax)<=1.0e-15)
        z_ob = cosmo->data.spline_kernel_n->interp_object.ymax;
    if(ln_mass_ob < cosmo->data.spline_kernel_n->interp_object.xmin
    || ln_mass_ob > cosmo->data.spline_kernel_n->interp_object.xmax
    || z_ob < cosmo->data.spline_kernel_n->interp_object.ymin
    || z_ob > cosmo->data.spline_kernel_n->interp_object.ymax){
        *status = CLUCU_ERROR_SPLINE_EV;
        CLUCU_RAISE_WARNING(*status,"mass=%.2e,z=%.2e is outside interpolation range [%.2e,%.2e][%.2e,%.2e]. [cen-min,max-cen]: ln_mass:[%.2e,%.2e] z:[%.2e,%.2e]",
                                        exp(ln_mass_ob),z_ob,
                                        exp(cosmo->data.spline_kernel_n->interp_object.xmin),
                                        exp(cosmo->data.spline_kernel_n->interp_object.xmax),
                                        cosmo->data.spline_kernel_n->interp_object.ymin,
                                        cosmo->data.spline_kernel_n->interp_object.ymax,

                                        ln_mass_ob-cosmo->data.spline_kernel_n->interp_object.xmin,
                                        cosmo->data.spline_kernel_n->interp_object.xmax-ln_mass_ob,
                                        z_ob-cosmo->data.spline_kernel_n->interp_object.ymin,
                                        cosmo->data.spline_kernel_n->interp_object.ymax-z_ob);
        }
    int gslstatus = gsl_spline2d_eval_e(cosmo->data.spline_kernel_n, ln_mass_ob,z_ob, NULL, NULL,&result);
    if(gslstatus != GSL_SUCCESS)
    {
        CLUCU_RAISE_GSL_WARNING(gslstatus,NULL);
        *status|= gslstatus;
    }
    return result;
}
double kernel_nb(clucu_cosmology *cosmo,double ln_mass_ob,double z_ob,int *status)
{
    double result;
    if(cosmo->data.spline_kernel_nb==NULL){
        *status = 243234;
        CLUCU_RAISE_WARNING(*status,"you shoule set kernel_nb");
    }
    if(fabs(z_ob - cosmo->data.spline_kernel_nb->interp_object.ymax)<=1.0e-15)
        z_ob = cosmo->data.spline_kernel_nb->interp_object.ymax;
    if(ln_mass_ob < cosmo->data.spline_kernel_nb->interp_object.xmin
    || ln_mass_ob > cosmo->data.spline_kernel_nb->interp_object.xmax
    || z_ob < cosmo->data.spline_kernel_nb->interp_object.ymin
    || z_ob > cosmo->data.spline_kernel_nb->interp_object.ymax){
        *status = CLUCU_ERROR_SPLINE_EV;
        CLUCU_RAISE_WARNING(*status,"mass=%.2e,z=%.2e is outside interpolation range [%.2e,%.2e][%.2e,%.2e]. [cen-min,max-cen]: ln_mass:[%.2e,%.2e] z:[%.2e,%.2e]",
                                        exp(ln_mass_ob),z_ob,
                                        exp(cosmo->data.spline_kernel_nb->interp_object.xmin),
                                        exp(cosmo->data.spline_kernel_nb->interp_object.xmax),
                                        cosmo->data.spline_kernel_nb->interp_object.ymin,
                                        cosmo->data.spline_kernel_nb->interp_object.ymax,

                                        ln_mass_ob-cosmo->data.spline_kernel_nb->interp_object.xmin,
                                        cosmo->data.spline_kernel_nb->interp_object.xmax-ln_mass_ob,
                                        z_ob-cosmo->data.spline_kernel_nb->interp_object.ymin,
                                        cosmo->data.spline_kernel_nb->interp_object.ymax-z_ob);
        }
    int gslstatus = gsl_spline2d_eval_e(cosmo->data.spline_kernel_nb, ln_mass_ob,z_ob, NULL, NULL,&result);
    if(gslstatus != GSL_SUCCESS)
    {
        CLUCU_RAISE_GSL_WARNING(gslstatus,NULL);
        *status|= gslstatus;
    }
    return result;
}

double barn_integrand_ln_mass_ob(double ln_mass_ob,void *params)
{
    int_pars *par=(int_pars *)params;
    double v = clucu_volume_element(par->cosmo,par->z_ob,par->status);
    double n = kernel_n(par->cosmo,ln_mass_ob,par->z_ob,par->status);//单位 1//Mpc^3
    return v*n;//单位1
}
void clucu_compute_barn(clucu_cosmology *cosmo,int *status)
{
    cosmo->data.barn_mzobP=Create2Grid(cosmo->data.Nmob,cosmo->data.NzobP);
    int thread_num = (int)( omp_get_num_procs() * THREAD_RATE/cosmo->thread_num);
    if(CHECK) thread_num=1;
    omp_set_num_threads(thread_num);
    #pragma omp parallel default(none) shared(cosmo,thread_num,status)
    {
        int local_status=*status;
      double ln_mass_ob2,ln_mass_ob1, z_obP2,z_obP1;
      int number_mob,number_zobP;
      int thread_id = omp_get_thread_num();
      for(int number = thread_id;number<(cosmo->data.Nmob-1)*(cosmo->data.NzobP-1);number+=thread_num)
      {
        number_mob = number / (cosmo->data.NzobP-1);
        number_zobP = number % (cosmo->data.NzobP-1);
        ln_mass_ob1 = cosmo->data.ln_mass_ob[number_mob];
        ln_mass_ob2 = cosmo->data.ln_mass_ob[number_mob+1];
        z_obP1 = cosmo->data.zobP[number_zobP];
        z_obP2 = cosmo->data.zobP[number_zobP+1];
        //不用对mass_true,z_true积分，因为kernel_n已经积过了
        cosmo->data.barn_mzobP[number_mob][number_zobP] = func_integrand_mass_z(cosmo,barn_integrand_ln_mass_ob,0,0,ln_mass_ob1,ln_mass_ob2,0,z_obP1,z_obP2,0,&local_status);
        printf("cosmo_id[%d]:barn_mzobP[%d][%d]=%f\n",cosmo->cosmo_id,number_mob,number_zobP,cosmo->data.barn_mzobP[number_mob][number_zobP]);
      }//end omp for
      if(local_status) {
      #pragma omp atomic write
      *status=local_status;
      }
    }//end omp parallel
    cosmo->computed_cluster_barn=true;
}

/*----------------------*/
/*          Chh         */
/*----------------------*/
double weight_h_integrand(double ln_mass_ob,void *params)
{
    int_pars *par=(int_pars *)params;
    double v = clucu_volume_element(par->cosmo,par->z_ob,par->status);//单位 Mpc^3
    double nb = kernel_nb(par->cosmo,ln_mass_ob,par->z_ob,par->status);//单位 1/Mpc^3
    double Hz = clucu_hubble_parameter(par->cosmo,par->z_ob,par->status);//单位1/Mpc

    return Hz * v *nb;//单位1/Mpc
}
double weight_h(clucu_cosmology *cosmo,double ln_mass_ob1,double ln_mass_ob2,double z_ob,int *status)
{
    //不用对mass_true,z_true积分，因为kernel_nb已经积过了
    return func_integrand_mass_z(cosmo,weight_h_integrand,0,0,ln_mass_ob1,ln_mass_ob2,0,z_ob,z_ob,-1,status);
}
double cl_hh_integrand(double chi,void *params)
{
    int_pars *par=(int_pars *)params;
    if(chi < 1.0e-15) 
        return 0;
    double z = clucu_chi_to_z(par->cosmo,chi,par->status);
    double wh = weight_h( par->cosmo,par->ln_mass_ob1,par->ln_mass_ob2,z,par->status);//还没除barn,单位1/Mpc
    double whh = wh;//还没除barn,单位1/Mpc
    if(par->ln_mass_ob1!=par->ln_mass_ob11 || par->ln_mass_ob2!=par->ln_mass_ob22)
        whh = weight_h( par->cosmo,par->ln_mass_ob11,par->ln_mass_ob22,z,par->status);
    double k = par->ell/chi;//单位 1/Mpc
    double pk = clucu_power(par->cosmo,k,z,par->status);//单位Mpc^3
    return wh*whh * pk /(chi*chi); //单位1/Mpc。还没除barn,除后单位1/Mpc
}
double clucu_cl_hh(clucu_cosmology *cosmo,double ell,double ln_mass_ob1,double ln_mass_ob2,double ln_mass_ob11,double ln_mass_ob22,double chi1,double chi2,int *status)
{
    int_pars par;
    par.cosmo = cosmo;
    par.ell = ell;
    par.ln_mass_ob1 = ln_mass_ob1;
    par.ln_mass_ob2 = ln_mass_ob2;
    par.ln_mass_ob11 = ln_mass_ob11;
    par.ln_mass_ob22 = ln_mass_ob22;
    par.status = status;

    gsl_integration_cquad_workspace *workspace =  NULL;
    workspace = gsl_integration_cquad_workspace_alloc(cosmo->gsl_param.N_ITERATION);
    gsl_function F;
    F.function=&cl_hh_integrand;
    F.params=&par;
    //确定积分的上下限
	double x_min,x_max,y;
    x_min=chi1;
    x_max=chi2;
    int gslstatus = gsl_integration_cquad(&F,
											x_min,
											x_max,
											0.0, cosmo->gsl_param.INTEGRATION_MASS_TRUE_EPSREL,
											workspace,&y,NULL,NULL);
	gsl_integration_cquad_workspace_free(workspace);
    return y;//单位1/Mpc。还没除barn,除后单位1/Mpc
}
//单位1/Mpc.两层积分，mass_ob,chi
void clucu_compute_cl_hh(clucu_cosmology *cosmo,int *status)
{
    cosmo->data.cl_hh = Create4Grid(cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.Nell);
    cosmo->data.cl_hh_hat = Create4Grid(cosmo->data.Nmob,cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.Nell);
    int Nm=cosmo->data.Nmob-1;
    int Nz=cosmo->data.NzobP-1;
    int Nell=cosmo->data.Nell;
    int thread_num = (int)( omp_get_num_procs() * THREAD_RATE/cosmo->thread_num);
	if(CHECK) thread_num=1;
	omp_set_num_threads(thread_num);
    #pragma omp parallel default(none) shared(Nm,Nz,Nell,cosmo,thread_num,status)
	{
        int local_status=*status;
		double ln_mass_ob2,ln_mass_ob1,ln_mass_ob22,ln_mass_ob11,chi1,chi2,ell,bar_n2;
		int iz,iell;
		int thread_id = omp_get_thread_num();
		for(int i = thread_id;i<Nz*Nell;i+=thread_num)
		{
            iz = i / Nell;
            iell = i % Nell;
            chi1 = clucu_comoving_angular_distance(cosmo,cosmo->data.zobP[iz],&local_status);
            chi2 = clucu_comoving_angular_distance(cosmo,cosmo->data.zobP[iz+1],&local_status);
            ell = cosmo->data.ell[iell];
            for(int im=0;im<Nm;im++)
            {
                ln_mass_ob1 = cosmo->data.ln_mass_ob[im];
                ln_mass_ob2 = cosmo->data.ln_mass_ob[im+1];
                for(int imm=0;imm<=im;imm++)
                {
                    ln_mass_ob11 = cosmo->data.ln_mass_ob[imm];
                    ln_mass_ob22 = cosmo->data.ln_mass_ob[imm+1];
                    bar_n2 = cosmo->data.barn_mzobP[im][iz]*cosmo->data.barn_mzobP[imm][iz];
                    cosmo->data.cl_hh[im][imm][iz][iell]=clucu_cl_hh(cosmo, ell,ln_mass_ob1,ln_mass_ob2,ln_mass_ob11,ln_mass_ob22, chi1,chi2,&local_status)/bar_n2;
                    cosmo->data.cl_hh_hat[im][imm][iz][iell]=cosmo->data.cl_hh[im][imm][iz][iell]+kronecker(im,imm)/cosmo->data.barn_mzobP[im][iz];
                }
            }
            //补全
            for(int im=0;im<Nm;im++){
                for(int imm=im+1;imm<Nm;imm++){
                    cosmo->data.cl_hh[im][imm][iz][iell] = cosmo->data.cl_hh[imm][im][iz][iell];
                    cosmo->data.cl_hh_hat[im][imm][iz][iell] = cosmo->data.cl_hh_hat[imm][im][iz][iell];
                }
            }
            printf("cosmo_id[%d]:cl_hh[%d][%d][%d][%d]=%le\n",cosmo->cosmo_id,Nm-1,Nm-1,iz,iell,cosmo->data.cl_hh[Nm-1][Nm-1][iz][iell]);    
        }
        if(local_status) {
        #pragma omp atomic write
        *status=local_status;
        }
	}
    cosmo->computed_cl_hh=true;
}

/*----------------------*/
/*        Ckk       */
/*----------------------*/

double cl_kk_integrand(double chi,void *params)
{
    //cosmo,ell
    int_pars *par=(int_pars *)params;
    //chi不能等于0
    if(chi<=1.0e-15) return 0;
    double z = clucu_chi_to_z(par->cosmo,chi,par->status);
    double H0 = clucu_hubble_parameter(par->cosmo,0.,par->status); //1/Mpc
    double Om = par->cosmo->param.Omega_m;
    double wk = 1.5*Om*H0*H0 * (1.+z) * chi * (1.-chi*par->cosmo->data.average_one_over_chis); //1/Mpc
    double k = par->ell/chi;//单位 1/Mpc
    //TODO：用非线性功率谱
    double pk = clucu_power(par->cosmo,k,z,par->status);//单位Mpc^3
    return wk*wk * pk /(chi*chi); //单位1/Mpc
}
double clucu_cl_kk(clucu_cosmology *cosmo,double ell,double chi1,double chi2,int *status)
{
    int_pars par;
    par.cosmo = cosmo;
    par.ell = ell;
    par.status = status;

    gsl_integration_cquad_workspace *workspace =  NULL;
    workspace = gsl_integration_cquad_workspace_alloc(cosmo->gsl_param.N_ITERATION);
    gsl_function F;
    F.function=&cl_kk_integrand;
    F.params=&par;
    //确定积分的上下限
	double x_min,x_max,y;
    x_min=chi1;
    x_max=chi2;
    int gslstatus = gsl_integration_cquad(&F,
											x_min,
											x_max,
											0.0, cosmo->gsl_param.INTEGRATION_MASS_TRUE_EPSREL,
											workspace,&y,NULL,NULL);
	gsl_integration_cquad_workspace_free(workspace);
    return y;//单位1/Mpc
}
//一层积分chi
void clucu_compute_cl_kk(clucu_cosmology *cosmo,int *status)
{
    cosmo->data.cl_kk = Create1Grid(cosmo->data.Nell);
    cosmo->data.cl_kk_hat = Create1Grid(cosmo->data.Nell);
    double sigma_e=0.35;
    double n_gal = 5.;
    n_gal = n_gal *pow(60*180/M_PI,2.);
    int thread_num = (int)( omp_get_num_procs() * THREAD_RATE/cosmo->thread_num);
	if(CHECK) thread_num=1;
	omp_set_num_threads(thread_num);
    #pragma omp parallel default(none) shared(cosmo,thread_num,sigma_e,n_gal,status)
    {
        int local_status=*status;
        int izobP,iell;
        double chi1,chi2,ell;
        int thread_id = omp_get_thread_num();
        for(int i = thread_id; i<cosmo->data.Nell;i+=thread_num)
        {
            iell = i ;
            ell = cosmo->data.ell[iell];
            //TODO: 确定积分范围。
            //根据画的error bar的图，可能是0-1.5，或者0-1.0
            //即，积分上限是z_smin=1.5，还是z_mean=1.0
            //不可能是1.5-5，0-5
            chi1=clucu_comoving_angular_distance(cosmo,0.,&local_status);
            chi2=clucu_comoving_angular_distance(cosmo,1.5,&local_status);
            cosmo->data.cl_kk[iell] = clucu_cl_kk(cosmo,ell,chi1,chi2,&local_status);
            cosmo->data.cl_kk_hat[iell] = cosmo->data.cl_kk[iell] + 0.5*sigma_e*sigma_e/n_gal;
        }
        if(local_status) {
        #pragma omp atomic write
        *status=local_status;
        }
    }
    cosmo->computed_cl_kk = true;
}
/*----------------------*/
/*        Chk       */
/*----------------------*/
double cl_hk2_integrand_ln_mass_ob(double ln_mass_ob,void *params)
{
    int_pars *par=(int_pars *)params;
    double chi = clucu_comoving_angular_distance(par->cosmo,par->z_ob,par->status);//单位Mpc
    //chi不能等于0
    if(chi<=1.0e-15) chi=clucu_comoving_angular_distance(par->cosmo,0.001,par->status);

    double v = clucu_volume_element(par->cosmo,par->z_ob,par->status);//单位 Mpc^3
    double nb = kernel_nb(par->cosmo,ln_mass_ob,par->z_ob,par->status);//单位 1/Mpc^3
    double Hz = clucu_hubble_parameter(par->cosmo,par->z_ob,par->status);//单位1/Mpc
    double wh = Hz * v *nb;//还没除barn,单位1/Mpc

    double H0=clucu_hubble_parameter(par->cosmo,0.,par->status); //1/Mpc
    double Om=par->cosmo->param.Omega_m;
    double wk = 1.5*Om*H0*H0 * (1.+par->z_ob) * chi * (1.-chi*par->cosmo->data.average_one_over_chis); //1/Mpc

    double k = par->ell/chi;//单位 1/Mpc
    double pk = clucu_power(par->cosmo,k,par->z_ob,par->status);//单位Mpc^3
    double result = wh*wk * pk /(chi*chi); //还没除barn,单位1/Mpc
    return result /Hz; //对z积分
}
//两层积分mass_ob,chi
void clucu_compute_cl_hk2(clucu_cosmology *cosmo,int CEN_massob,int CEN_zob,int *status)
{
    cosmo->data.cl_hk2 = Create3Grid(cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.Nell);
    int Nm=cosmo->data.Nmob-1;
    int Nz=cosmo->data.NzobP-1;
    int Nell=cosmo->data.Nell;
    int thread_num = (int)( omp_get_num_procs() * THREAD_RATE/cosmo->thread_num);
	if(CHECK) thread_num=1;
	omp_set_num_threads(thread_num);
    #pragma omp parallel default(none) shared(Nm,Nz,Nell,cosmo,thread_num,status,CEN_massob,CEN_zob)
	{
        int local_status=*status;
		double ln_mass_ob2,ln_mass_ob1,chi1,chi2,ell,bar_n2;
        double z_ob1,z_ob2;
		int im,iz,iell;
		int thread_id = omp_get_thread_num();
		for(int i = thread_id;i<Nm*Nz*Nell;i+=thread_num)
		{
            im = i / (Nz*Nell);
            iz = (  i % (Nz*Nell)  ) / Nell;
            iell = (  i % (Nz*Nell)  ) % Nell;

            ln_mass_ob1 = cosmo->data.ln_mass_ob[im];
            ln_mass_ob2 = cosmo->data.ln_mass_ob[im+1];
            z_ob1 = cosmo->data.zobP[iz];
            z_ob2 = cosmo->data.zobP[iz+1];
            ell = cosmo->data.ell[iell];
            cosmo->data.cl_hk2[im][iz][iell]=func_integrand_mass_z_ell(cosmo,cl_hk2_integrand_ln_mass_ob,0,0, ell,ln_mass_ob1,ln_mass_ob2,CEN_massob, z_ob1,z_ob2,CEN_zob,&local_status)/cosmo->data.barn_mzobP[im][iz];

            printf("cosmo_id[%d]:hk2[%d][%d][%d]=%le\n",cosmo->cosmo_id,im,iz,iell,cosmo->data.cl_hk2[im][iz][iell]); 
        }
        if(local_status) {
        #pragma omp atomic write
        *status=local_status;
        }
	}    
    cosmo->computed_cl_hk2=true;
}
/*----------------------*/
/*        Chk       */
/*----------------------*/
double one_over_sigma_critical(clucu_cosmology *cosmo,double z,int *status)
{
    double chi=clucu_comoving_angular_distance(cosmo,z,status);//Mpc
    double one_over_chis = cosmo->data.average_one_over_chis;//Mpc
    double G = clucu_constants.GNEWT; //单位：m^3 / kg s^2
    double M_SUN = clucu_constants.SOLAR_MASS;//kg
    double MPC = clucu_constants.MPC_TO_METER; //m 
    double c = clucu_constants.CLIGHT;// m/s
    double G0 = G * M_SUN / pow(MPC, 3.0) /pow(c/MPC,2.);   //单位：Mpc^3/M_sun s^2
    return 4.*M_PI*G0 * chi *(1.-chi*one_over_chis) / (1.+z); // Mpc^2/Msun^2
}
//这里的质量是mass true
double kappa_lMz(clucu_cosmology *cosmo,double ell,double mass,double z,int *status)
{
    double chi = clucu_comoving_angular_distance(cosmo,z,status);//Mpc
    double uk = nfw_uk( cosmo, ell/chi, mass, z,status);//1,k也不带h
    double m = mass ; //Msun
    double sigma_inverse = one_over_sigma_critical(cosmo,z,status); //Mpc^2/Msun
    double kappa = m * uk * pow(1.+z,2.) *sigma_inverse/ (chi*chi);//1

    //off center
    double Mp=3.0e14/cosmo->param.h;//Msun
//    double fcen=0.75+0.05*log(mass/Mp);
    double fcen = cosmo->survey.f_cen0 + cosmo->survey.p_cenM*log(mass/Mp) + cosmo->survey.p_cenz*log(1.+z);
    double dA=clucu_angular_diameter_distance(cosmo,z,status);//Mpc
//    double sigma_s=0.42/ cosmo->param.h;//Mpc/h转换成Mpc
//    sigma_s /= dA;//单位1
    double sigma_s0 = cosmo->survey.sigma_s0/cosmo->param.h/dA;//Mpc/h转换成Mpc,转换成单位1
    double sigma_s = sigma_s0 + cosmo->survey.p_sigmaM*log(mass/Mp) + cosmo->survey.p_sigmaz*log(1.+z);
    double A = fcen + (1.-fcen)*exp(-0.5*sigma_s*ell*sigma_s*ell);

    kappa = A*kappa;
    return kappa;
}

double cl_hk1_integrand_ln_mass_true(double ln_mass_true,void *params)
{
    int_pars *par=(int_pars *)params;
    //如果z=0,V(z)=0，直接return 0
    if(par->z_ob<=1.0e-15) 
        return 0; 
    
    double chi = clucu_comoving_angular_distance(par->cosmo,par->z_ob,par->status);//单位Mpc

    double dndlnm=clucu_halo_dn_dlnm(par->cosmo,exp(ln_mass_true),par->z_ob,par->status);
    double pl = mass_ob_probability(par->cosmo,par->ln_mass_ob,ln_mass_true,par->z_ob);
    double n = dndlnm * pl;
    double v = clucu_volume_element(par->cosmo,par->z_ob,par->status);
    double Hz = clucu_hubble_parameter(par->cosmo,par->z_ob,par->status);//单位1/Mpc
    double whk1 = Hz * v *n;//还没除barn,单位1/Mpc

    double k = par->ell/chi;//单位 1/Mpc
    double kappa = kappa_lMz(par->cosmo,par->ell,exp(ln_mass_true),par->z_ob,par->status);//1
    
    double result = whk1 * kappa ;//还没除barn,单位1/Mpc
    return result / Hz;//对z积分
}
//3层积分，mass_true,mass_ob,chi。另外还有参数ell
void clucu_compute_cl_hk1(clucu_cosmology *cosmo,int CEN_massob,int CEN_zob,int *status)
{
    cosmo->data.cl_hk1 = Create3Grid(cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.Nell);
    cosmo->data.cl_hk = Create3Grid(cosmo->data.Nmob,cosmo->data.NzobP,cosmo->data.Nell);
    int Nm=cosmo->data.Nmob-1;
    int Nz=cosmo->data.NzobP-1;
    int Nell=cosmo->data.Nell;
    int thread_num = (int)( omp_get_num_procs() * THREAD_RATE/cosmo->thread_num);
	if(CHECK) thread_num=1;
	omp_set_num_threads(thread_num);
    #pragma omp parallel default(none) shared(Nm,Nz,Nell,cosmo,thread_num,status,CEN_massob,CEN_zob)
	{
        int local_status=*status;
		double ln_mass_ob2,ln_mass_ob1,chi1,chi2,ell,bar_n2;
        double z_ob1,z_ob2;
		int im,iz,iell;
		int thread_id = omp_get_thread_num();
		for(int i = thread_id;i<Nm*Nz*Nell;i+=thread_num)
		{
            im = i / (Nz*Nell);
            iz = (  i % (Nz*Nell)  ) / Nell;
            iell = (  i % (Nz*Nell)  ) % Nell;

            ln_mass_ob1 = cosmo->data.ln_mass_ob[im];
            ln_mass_ob2 = cosmo->data.ln_mass_ob[im+1];
            z_ob1 = cosmo->data.zobP[iz];
            z_ob2 = cosmo->data.zobP[iz+1];
            ell = cosmo->data.ell[iell];
            cosmo->data.cl_hk1[im][iz][iell]=func_integrand_mass_z_ell(cosmo,cl_hk1_integrand_ln_mass_true,1,0, ell,ln_mass_ob1,ln_mass_ob2,CEN_massob, z_ob1,z_ob2,CEN_zob,&local_status)/cosmo->data.barn_mzobP[im][iz];
            printf("cosmo_id[%d]:cl_hk1[%d][%d][%d]=%le\n",cosmo->cosmo_id,im,iz,iell,cosmo->data.cl_hk1[im][iz][iell]); 
        }
        if(local_status) {
        #pragma omp atomic write
        *status=local_status;
        }
	}

    for(int iz=0;iz<Nz;iz++)
    {
        for(int iell=0;iell<Nell;iell++)
        {
            for(int im=0;im<Nm;im++)
                cosmo->data.cl_hk[im][iz][iell]=cosmo->data.cl_hk1[im][iz][iell]+cosmo->data.cl_hk2[im][iz][iell];
        }
    }
    cosmo->computed_cl_hk1=true;
    cosmo->computed_cl_hk=true;
}

/*----------------------*/
/*        save       */
/*----------------------*/

void SaveClz_hh(clucu_cosmology *cosmo)
{
	char namedir[50];
    char nameNlzob[100];
	sprintf(namedir,"../output/%s", cosmo->name);
	sprintf(nameNlzob,"../output/%s/Clz%d_hh.dat", cosmo->name,cosmo->cosmo_id);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameNlzob,"w");
    save_matrix(fp,cosmo->data.cl_hh[0][0],cosmo->data.zobP,cosmo->data.ell,cosmo->data.NzobP,cosmo->data.Nell);
	fclose(fp);
}
void SaveClz_hh_hat(clucu_cosmology *cosmo)
{
	char namedir[50];
    char nameNlzob[100];
	sprintf(namedir,"../output/%s", cosmo->name);
	sprintf(nameNlzob,"../output/%s/Clz%d_hh_hat.dat", cosmo->name,cosmo->cosmo_id);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameNlzob,"w");
    save_matrix(fp,cosmo->data.cl_hh_hat[0][0],cosmo->data.zobP,cosmo->data.ell,cosmo->data.NzobP,cosmo->data.Nell);
	fclose(fp);
}

void SaveClz_kk(clucu_cosmology *cosmo)
{
	char namedir[50];
    char nameNlzob[100];
	sprintf(namedir,"../output/%s", cosmo->name);
	sprintf(nameNlzob,"../output/%s/Clz%d_kk.dat", cosmo->name,cosmo->cosmo_id);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameNlzob,"w");
    for(int iell=0;iell<cosmo->data.Nell;iell++)
		fprintf(fp,"%f      %le\n",cosmo->data.ell[iell],cosmo->data.cl_kk[iell]);
	fclose(fp);
}
void SaveClz_kk_hat(clucu_cosmology *cosmo)
{
	char namedir[50];
    char nameNlzob[100];
	sprintf(namedir,"../output/%s", cosmo->name);
	sprintf(nameNlzob,"../output/%s/Clz%d_kk_hat.dat", cosmo->name,cosmo->cosmo_id);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameNlzob,"w");
    for(int iell=0;iell<cosmo->data.Nell;iell++)
		fprintf(fp,"%f      %le\n",cosmo->data.ell[iell],cosmo->data.cl_kk_hat[iell]);
	fclose(fp);
}

void SaveClz_hk1(clucu_cosmology *cosmo)
{
	char namedir[50];
    char nameNlzob[100];
	sprintf(namedir,"../output/%s", cosmo->name);
	sprintf(nameNlzob,"../output/%s/Clz%d_hk1.dat", cosmo->name,cosmo->cosmo_id);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameNlzob,"w");
    save_matrix(fp,cosmo->data.cl_hk1[0],cosmo->data.zobP,cosmo->data.ell,cosmo->data.NzobP,cosmo->data.Nell);
	fclose(fp);
}
void SaveClz_hk2(clucu_cosmology *cosmo)
{
	char namedir[50];
    char nameNlzob[100];
	sprintf(namedir,"../output/%s", cosmo->name);
	sprintf(nameNlzob,"../output/%s/Clz%d_hk2.dat", cosmo->name,cosmo->cosmo_id);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameNlzob,"w");
    save_matrix(fp,cosmo->data.cl_hk2[0],cosmo->data.zobP,cosmo->data.ell,cosmo->data.NzobP,cosmo->data.Nell);
	fclose(fp);
}
void SaveClz_hk(clucu_cosmology *cosmo)
{
	char namedir[50];
    char nameNlzob[100];
	sprintf(namedir,"../output/%s", cosmo->name);
	sprintf(nameNlzob,"../output/%s/Clz%d_hk.dat", cosmo->name,cosmo->cosmo_id);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameNlzob,"w");
    save_matrix(fp,cosmo->data.cl_hk[0],cosmo->data.zobP,cosmo->data.ell,cosmo->data.NzobP,cosmo->data.Nell);
	fclose(fp);
}