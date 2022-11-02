#include "clucu.h"

#include <time.h>
#include <string.h>

typedef struct {
    clucu_cosmology cosmo;
    double z;
    double chi;
    double ell;
    double ln_mass1;
    double ln_mass2;
    double ln_mass11;
    double ln_mass22;
    double ln_mass_ob;
} cl_pars;
int main()
{
    int status=0;
    omp_set_nested(1);
    //用于设置是否允许OpenMP进行嵌套并行，默认的设置为false。
	time_t start,end;
    start =time(NULL);//or time(&start);
    //创建 fiducial model
    clucu_cosmology fiducial=clucu_cosmology_create(config_oguri,param_oguri,survey_oguri);
//    fiducial.config.matter_power_spectrum_method=clucu_boltzmann_class;
     //fiducial.param.sigma8=0.79583;
     fiducial.param.sigma8=0.7959;
     fiducial.param.delta_zeta=NAN;
    fiducial.runned_class = true;
    strcpy(fiducial.classname,"oguri");
	strcpy(fiducial.name,"oguri1");
    
    //对 fiducial model 稍作修改
    fiducial.survey.mass_ob_min_h=pow(10.,15.);
    fiducial.survey.mass_ob_max_h=pow(10.,16.);
    fiducial.survey.mass_ob_maxc_h=NAN;


    fiducial.survey.dlog10_mass_ob=1.;
    fiducial.survey.zobP_min=0.7;
    fiducial.survey.zobP_max=0.8;
    fiducial.survey.zobP_bin=0.1;
    fiducial.survey.zob_min=0.0;
    fiducial.survey.zob_max=1.5;
    fiducial.survey.zob_bin=0.05;
	
	clucu_cosmology cosmo = fiducial;
    
	cosmo.cosmo_id = 0;
	cosmo.class_id = 0;
	InitializeCosmo(&cosmo);
    clucu_compute_background(&cosmo,clucu_background_class_label,&status);
    clucu_compute_power(&cosmo,&status);

    int Nm=cosmo.data.Nmob;
    int Nz=cosmo.data.NzobP;
    int Nell=cosmo.data.Nell;
    clucu_compute_logsigma(&cosmo,&status);//功率谱相同即可
    printf("log sigma --over\n");
    clucu_compute_kernel_n( &cosmo,&status);//需要mass-ob 
    printf("kernel_n --over\n");
    clucu_compute_kernel_nb( &cosmo,&status);
    printf("kernel_nb --over\n");
    clucu_compute_barn( &cosmo,&status);
    printf("barn --over\n");
    clucu_compute_cluster_averagebias( &cosmo,&status);
    printf("averagebias --over\n");

    //测试NN
    clucu_compute_cluster_numbercounts(&cosmo,0,0,&status);SaveNmzob(cosmo);
    clucu_compute_samplingvariance(&cosmo,&status);
    clucu_compute_cov_NN(&cosmo,&status);SaveCov_NN(&cosmo);
    
    //测试hh
    clucu_compute_cl_hh( &cosmo,&status);SaveClz_hh(&cosmo);SaveClz_hh_hat(&cosmo);
    clucu_compute_cov_hhhh(&cosmo,&status);SaveCov_hhhh(&cosmo);
    
    clucu_compute_average_one_over_chis(&cosmo);
    clucu_compute_cl_kk(&cosmo,&status);SaveClz_kk(&cosmo);SaveClz_kk_hat(&cosmo);

    clucu_compute_cl_hk2(&cosmo,0,1,&status);SaveClz_hk2(&cosmo);
    clucu_compute_cl_hk1(&cosmo,0,1,&status);SaveClz_hk1(&cosmo);SaveClz_hk(&cosmo);
    clucu_compute_cov_hkhk(&cosmo,&status);SaveCov_hkhk(&cosmo);


    double ln_mass_ob1 = cosmo.data.ln_mass_ob[1];
    double ln_mass_ob2 = cosmo.data.ln_mass_ob[1+1];
    double z_obP1 = cosmo.data.zobP[1];
    double z_obP2 = cosmo.data.zobP[1+1];

/*
    测试一些被积函数
    cl_pars par;
    par.cosmo=cosmo;
    par.ell=100.;
    par.ln_mass1 = ln_mass_ob1;
    par.ln_mass2 = ln_mass_ob2;
    par.ln_mass11 = ln_mass_ob1;
    par.ln_mass22 = ln_mass_ob2;

    char namedir[50];
    char nameNlzob[100];
	sprintf(namedir,"../output/%s", cosmo.name);
	sprintf(nameNlzob,"../output/%s/integrand1.dat", cosmo.name);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameNlzob,"w");
    double Nzlist=1000;
    double *ell=clucu_log10_spacing(1.,1.0e4,100);
    //double *z_list=clucu_log10_spacing(1.0e-5,100.,Nzlist);
    //for(int iz=0;iz<Nzlist;iz++)
    for(int iell=0;iell<100;iell++)
    {
        par.ell=ell[iell];
        double y=clucu_cl_kk(cosmo, ell[iell],0.,ComovingDistance(cosmo,1.5));
        fprintf(fp,"%e %le\n",ell[iell],y);
    }
    fclose(fp);
*/
    //测试质量转换
    //double M200m=M500cz_to_M200mz(cosmo,0.321533e14,0.);



/*
    printf("%f,%f\n",log(cosmo.survey.mass_ob_min),log(cosmo.survey.mass_ob_max));
    printf("%f,%f\n",cosmo.survey.zob_min,cosmo.survey.zob_max);

    //--check kernel_n-
    printf("kernel_n=%le,prob=%le\n",kernel_n(cosmo,ln_mass_ob1,z_obP1),average_hmf_integrand( cosmo, ln_mass_ob1, z_obP1));
    //--check kernel_nb-
    printf("kernel_nb=%le,prob=%le\n",kernel_nb(cosmo,ln_mass_ob1,z_obP1),average_hmfb_integrand( cosmo, ln_mass_ob1, z_obP1));
    //--check barn and sumHMF[iz]-
    double *barn=Create1Grid(Nz);
    for(int iz=0;iz<Nz-1;iz++)
    {
        double sum=0.;
        for(int im=0;im<Nm;im++)
            sum += cosmo.data.barn_mzobP[im][iz] ;
        barn[iz]=sum;
        double z=(cosmo.data.zobP[iz]+cosmo.data.zobP[iz+1])/2.;
        printf("%d,barn=%le,aven=%le\n",iz,barn[iz],cosmo.data.sumHMF[iz]*cosmo.survey.zobP_bin*CalculateVolumeElement(cosmo,z));
    }
*/






    

    


    end =time(NULL);
    printf("time=%fs\n",difftime(end,start));
	printf("time=%fm\n",difftime(end,start)/60.);
	printf("time=%fh\n",difftime(end,start)/3600.);
    return 1;
}