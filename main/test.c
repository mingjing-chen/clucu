#include <math.h>
#include <stdio.h>
#include "clucu.h"

int main()
{
    int status=0;
    omp_set_nested(1);
    //用于设置是否允许OpenMP进行嵌套并行，默认的设置为false。

    //创建 fiducial model
    clucu_cosmology fiducial=clucu_cosmology_create(config_default,param_WMAP5,survey_oguri);
    fiducial.config.matter_power_spectrum_method=clucu_boltzmann_class_Pk;
    fiducial.runned_class = true;
    strcpy(fiducial.classname,"Des");
	strcpy(fiducial.name,"Des");
	
	clucu_cosmology cosmo = fiducial;
    cosmo.param.f_nl=300.;
    
	cosmo.cosmo_id = 0;
	cosmo.class_id = 0;
	InitializeCosmo(&cosmo);
    clucu_compute_background(&cosmo,clucu_background_cla_label,&status);
    clucu_compute_power(&cosmo,&status);

    clucu_compute_logsigma(&cosmo,&status);
    printf("log sigma --over\n");
    clucu_compute_sigmaS3(&cosmo,&status);
    printf("sigmaS3 --over\n");

    char namedir[50];
    char namefile[100];
	sprintf(namedir,"../output/%s", cosmo.name);
	sprintf(namefile,"../output/%s/NGhmf_factor.dat", cosmo.name);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(namefile,"w");

    
    double z=1.;
    int Nm = 100;
    double *M = clucu_log10_spacing(1.0e13,1.0e16,Nm);
    double sigma,nu,factor,dndlnm;
    for(int i=0;i<Nm;i++){
        sigma = clucu_sigma(&cosmo,M[i],z,&status);
        nu = 1.686/sigma;
        factor = NGhmf_factor(&cosmo,sigma,M[i],z,&status);
        dndlnm = clucu_halo_dn_dlnm(&cosmo,M[i],z,&status);
        fprintf(fp,"%e %e %le %le\n",M[i]*cosmo.param.h,nu,factor,dndlnm/pow(cosmo.param.h,3.));
    }

    fclose(fp);
    return 1;

}