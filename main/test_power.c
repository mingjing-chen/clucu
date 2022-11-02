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
//    clucu_cosmology fiducial=clucu_cosmology_create(config_default,param_oguri,survey_oguri);
/*    fiducial.config.mass_observable_method =   clucu_oguri2011;//完全一样
    fiducial.config.bias_function_method = clucu_bias_bhattacharya2011;//完全一样
    fiducial.config.halo_concentration_method = clucu_concentration_duffy2008;
    fiducial.config.halo_define_method =  clucu_delta_vir;//完全一样
    fiducial.config.mass_function_method =         clucu_massfunc_bhattacharya2011;//变大一点点
*/
    fiducial.config.matter_power_spectrum_method=clucu_eisenstein_hu;//红移=0.6开始，尾部变大一点
//fiducial.config.matter_power_spectrum_method=clucu_bbks;
//fiducial.config.matter_power_spectrum_method=clucu_boltzmann_class;
    fiducial.runned_class = false;
	strcpy(fiducial.classname,"oguri");
	strcpy(fiducial.name,"test_power");
	clucu_cosmology cosmo = fiducial;
    
	cosmo.cosmo_id = 0;
	cosmo.class_id = 0;
	InitializeCosmo(&cosmo);
    clucu_compute_background(&cosmo,clucu_background_cla_label,&status);
    clucu_compute_power(&cosmo,&status);

    int nk=1000;
    double *k=clucu_log10_spacing(1.0e-5,1.0e3,nk);
    char namedir[50];
    char namefile[100];
	sprintf(namedir,"../output/%s", cosmo.name);
	sprintf(namefile,"../output/%s/bbks.dat", cosmo.name);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(namefile,"w");
    for(int ik=0;ik<nk;ik++)
    {
        fprintf(fp,"%e %le\n",k[ik],sqrt(tsqr_BBKS(&cosmo,k[ik])));
    }
    fclose(fp);

    sprintf(namefile,"../output/%s/eh.dat", cosmo.name);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	fp=fopen(namefile,"w");
    eh_struct *eh = NULL;
    eh = clucu_eh_struct_new(&cosmo,1);
    for(int ik=0;ik<nk;ik++)
    {
        fprintf(fp,"%e %le\n",k[ik],sqrt(tsqr_EH(&cosmo,eh,k[ik])));
    }
    fclose(fp);



    return 1;
}