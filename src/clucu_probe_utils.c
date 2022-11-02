#include "clucu_probe_utils.h"


/*-----------------------------------------------------------------------------*/
/*                                       积分算法                                  */
/*-----------------------------------------------------------------------------*/
// Integrand for numbercounts integral: d z_true，对dlnm积分，得到dz的被积函数
static double func_z_true(double z_true, void *params)
{
	int_pars *par=(int_pars *)params;
	int_pars par_in;
	par_in.cosmo=par->cosmo;
	par_in.ln_mass_ob=par->ln_mass_ob;
	par_in.z_ob=par->z_ob;
	par_in.z_true = z_true;
    par_in.status=par->status;
    par_in.ell=par->ell;
    par_in.k=par->k;

    if ( par_in.cosmo->survey.mass_ob_sigma0<=1.0e-15 || !par->if_int_masstrue)
        return par->func( par_in.ln_mass_ob, &par_in);

	gsl_integration_cquad_workspace *workspace =  NULL;
    gsl_function F;
    F.function=par->func;
    F.params=&par_in;

    double result;
	workspace = gsl_integration_cquad_workspace_alloc(par->cosmo->gsl_param.N_ITERATION);
    if (workspace == NULL) {
    par->cosmo->status = CLUCU_ERROR_MEMORY;
    }


    int gslstatus = gsl_integration_cquad(&F,
                                        log(pow(10.,par->cosmo->spline_param.LOG10M_SPLINE_MIN)),
                                        log(pow(10.,par->cosmo->spline_param.LOG10M_SPLINE_MAX)),
                                        0.0, par->cosmo->gsl_param.INTEGRATION_MASS_TRUE_EPSREL,
                                        workspace,&result,NULL,NULL);
    gsl_integration_cquad_workspace_free(workspace);

    return result;
}
// Integrand for numbercounts，对dz_true积分
static double func_ln_mass_ob(double ln_mass_ob,void *params)
{
    int_pars *par=(int_pars *)params;
	int_pars par_in;
	par_in.cosmo=par->cosmo;
	par_in.ln_mass_ob=ln_mass_ob;
	par_in.z_ob=par->z_ob;
    par_in.status=par->status;
    par_in.if_int_masstrue=par->if_int_masstrue;
    par_in.func=par->func;
    par_in.ell=par->ell;
    par_in.k=par->k;


    if ( par->cosmo->survey.redshift_sigma0<=1.0e-15 || !par->if_int_ztrue )
        return func_z_true( par->z_ob, &par_in);

	gsl_integration_cquad_workspace *workspace =  NULL;
    gsl_function F;
    F.function=func_z_true;
    F.params=&par_in;
    double result=0.;

	workspace = gsl_integration_cquad_workspace_alloc(par->cosmo->gsl_param.N_ITERATION);
    if (workspace == NULL) {
    par->cosmo->status = CLUCU_ERROR_MEMORY;
    }

	//确定积分的上下限
    double out[5];
	double x_min,x_max,y;
    zob_to_ztrue(par->cosmo,par->z_ob,out);
	for(int i=0;i<4;i++)
	{
		double x_min=out[i];
		double x_max=out[i+1];   
		if(fabs(x_max-x_min)<1.0e-10) continue;
		//开始对红移积分
		int gslstatus = gsl_integration_cquad(&F,
											x_min,
											x_max,
											0.0, par->cosmo->gsl_param.INTEGRATION_Z_TRUE_EPSREL,
											workspace,&y,NULL,NULL);
		result += y;
	}
	gsl_integration_cquad_workspace_free(workspace);

    return result;  
}
static double func_z_ob(double z_ob,void *params)
{
    int_pars *par=(int_pars *)params;
	int_pars par_in;
	par_in.cosmo=par->cosmo;
	par_in.z_ob=z_ob;
    par_in.status=par->status;
    par_in.if_int_ztrue=par->if_int_ztrue;
    par_in.if_int_masstrue=par->if_int_masstrue;
    par_in.func=par->func;
    par_in.ell=par->ell;
    par_in.k=par->k;
    
    //如果不对mass ob积分，ln_mass_ob1=ln_mass_ob2
    if (par->ln_mass_ob1==par->ln_mass_ob2)
        return func_ln_mass_ob(par->ln_mass_ob1,&par_in);

    if (par->CEN_massob && (log10(exp(par->ln_mass_ob2))-log10(exp(par->ln_mass_ob1)))< 0.3)
        return func_ln_mass_ob( 0.5*(par->ln_mass_ob1+par->ln_mass_ob2), &par_in)*(par->ln_mass_ob2-par->ln_mass_ob1);

	gsl_integration_cquad_workspace *workspace =  NULL;
    gsl_function F;
    F.function=func_ln_mass_ob;
    F.params=&par_in;

	workspace = gsl_integration_cquad_workspace_alloc(par->cosmo->gsl_param.N_ITERATION);
    if (workspace == NULL) {
    par->cosmo->status = CLUCU_ERROR_MEMORY;
    }

	double x_min,x_max,y;
    x_min=par->ln_mass_ob1;
    x_max=par->ln_mass_ob2;
    int gslstatus = gsl_integration_cquad(&F,
                                        x_min,
                                        x_max,
                                        0.0, par->cosmo->gsl_param.INTEGRATION_MASS_OB_EPSREL,
                                        workspace,&y,NULL,NULL);
	gsl_integration_cquad_workspace_free(workspace);

    return y;

}

/**
 * @brief 依次对mass_true,mass_ob,z_true,z_ob积分。
 * 如果ln_mass_ob1=ln_mass_ob2：不对mass_ob积分；
 * 如果z_ob1=z_ob2：不对z_ob积分；
 * 
 * @param cosmo         cosmo结构体
 * @param func          被积函数，要求x是mass_true， p是参数
 * @param if_int_masstrue   是否对mass_true积分，如果没有pm，默认不积分
 * @param if_int_ztrue      是否对z_true积分，如果没有pz，默认不积分
 * @param ln_mass_ob1   积分下限
 * @param ln_mass_ob2   积分上限
 * @param CEN_massob    对mass_ob积分时，是否直接用中间值的矩形积分法
 * @param z_ob1         积分下限
 * @param z_ob2         积分上限 
 * @param CEN_zob       对z_ob积分时，是否直接用中间值的矩形积分法 
 * @param status 
 * @return double 
 */
double func_integrand_mass_z(clucu_cosmology *cosmo,double (*func)(double x,void *p),int if_int_masstrue,int if_int_ztrue,double ln_mass_ob1,double ln_mass_ob2,int CEN_massob,double z_ob1,double z_ob2,int CEN_zob,int *status)
{
	int_pars par_in;
 	par_in.cosmo=cosmo;
	par_in.ln_mass_ob1=ln_mass_ob1;
    par_in.ln_mass_ob2=ln_mass_ob2;
    par_in.status=status;
    par_in.CEN_massob=CEN_massob;
    par_in.if_int_ztrue=if_int_ztrue;
    par_in.if_int_masstrue=if_int_masstrue;
    par_in.func=func;

    //如果不对z z_ob1=z_ob2
    if (z_ob1==z_ob2)
        return func_z_ob(z_ob1,&par_in);
        
    if (CEN_zob)
        return func_z_ob( 0.5*(z_ob1+z_ob2), &par_in)*(z_ob2-z_ob1);

	gsl_integration_cquad_workspace *workspace =  NULL;
    gsl_function F;
    F.function=func_z_ob;
    F.params=&par_in;

	workspace = gsl_integration_cquad_workspace_alloc(cosmo->gsl_param.N_ITERATION);
    if (workspace == NULL) {
    cosmo->status = CLUCU_ERROR_MEMORY;
    }
	double x_min,x_max,y;
    x_min=z_ob1;
    x_max=z_ob2;
    int gslstatus = gsl_integration_cquad(&F,
                                        x_min,
                                        x_max,
                                        0.0, cosmo->gsl_param.INTEGRATION_Z_OB_EPSREL,
                                        workspace,&y,NULL,NULL);
	gsl_integration_cquad_workspace_free(workspace);

    return y;
}
double func_integrand_mass_z_ell(clucu_cosmology *cosmo,double (*func)(double x,void *p),int if_int_masstrue,int if_int_ztrue,double ell,double ln_mass_ob1,double ln_mass_ob2,int CEN_massob,double z_ob1,double z_ob2,int CEN_zob,int *status)
{
	int_pars par_in;
 	par_in.cosmo=cosmo;
	par_in.ln_mass_ob1=ln_mass_ob1;
    par_in.ln_mass_ob2=ln_mass_ob2;
    par_in.status=status;
    par_in.CEN_massob=CEN_massob;
    par_in.if_int_ztrue=if_int_ztrue;
    par_in.if_int_masstrue=if_int_masstrue;
    par_in.func=func;
    par_in.ell=ell;

    //如果不对z z_ob1=z_ob2
    if (z_ob1==z_ob2)
        return func_z_ob(z_ob1,&par_in);
        
    if (CEN_zob)
        return func_z_ob( 0.5*(z_ob1+z_ob2), &par_in)*(z_ob2-z_ob1);

	gsl_integration_cquad_workspace *workspace =  NULL;
    gsl_function F;
    F.function=func_z_ob;
    F.params=&par_in;

	workspace = gsl_integration_cquad_workspace_alloc(cosmo->gsl_param.N_ITERATION);
    if (workspace == NULL) {
    cosmo->status = CLUCU_ERROR_MEMORY;
    }
	double x_min,x_max,y;
    x_min=z_ob1;
    x_max=z_ob2;
    int gslstatus = gsl_integration_cquad(&F,
                                        x_min,
                                        x_max,
                                        0.0, cosmo->gsl_param.INTEGRATION_Z_OB_EPSREL,
                                        workspace,&y,NULL,NULL);
	gsl_integration_cquad_workspace_free(workspace);

    return y;
}
double func_integrand_mass_z_k(clucu_cosmology *cosmo,double (*func)(double x,void *p),int if_int_masstrue,int if_int_ztrue,double k,double ln_mass_ob1,double ln_mass_ob2,int CEN_massob,double z_ob1,double z_ob2,int CEN_zob,int *status)
{
	int_pars par_in;
 	par_in.cosmo=cosmo;
	par_in.ln_mass_ob1=ln_mass_ob1;
    par_in.ln_mass_ob2=ln_mass_ob2;
    par_in.status=status;
    par_in.CEN_massob=CEN_massob;
    par_in.if_int_ztrue=if_int_ztrue;
    par_in.if_int_masstrue=if_int_masstrue;
    par_in.func=func;
    par_in.k=k;

    //如果不对z z_ob1=z_ob2
    if (z_ob1==z_ob2)
        return func_z_ob(z_ob1,&par_in);
        
    if (CEN_zob)
        return func_z_ob( 0.5*(z_ob1+z_ob2), &par_in)*(z_ob2-z_ob1);

	gsl_integration_cquad_workspace *workspace =  NULL;
    gsl_function F;
    F.function=func_z_ob;
    F.params=&par_in;

	workspace = gsl_integration_cquad_workspace_alloc(cosmo->gsl_param.N_ITERATION);
    if (workspace == NULL) {
    cosmo->status = CLUCU_ERROR_MEMORY;
    }
	double x_min,x_max,y;
    x_min=z_ob1;
    x_max=z_ob2;
    int gslstatus = gsl_integration_cquad(&F,
                                        x_min,
                                        x_max,
                                        0.0, cosmo->gsl_param.INTEGRATION_Z_OB_EPSREL,
                                        workspace,&y,NULL,NULL);
	gsl_integration_cquad_workspace_free(workspace);

    return y;
}