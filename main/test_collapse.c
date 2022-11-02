#include "clucu.h"

#include <time.h>
#include <string.h>
// gcc main_csst.c  -I/home/mjchen/Program/cluster/newcode/include -I/home/mjchen/Program/gsl-2.7/gsl  -L/home/mjchen/Program/gsl-2.7/.libs -lgsl -lgslcblas -lm -fopenmp  -o test


const clucu_parameter param_collapse={
    .Neff = 3.046,
    //.N_nu_mass =3,
    //.Onuh2_mass = 0.,
    .m_nu = {0.0,0.0,0.0},
    //.m_nu_sum = 0.6,
    .neutrino_mass_split_label = clucu_nu_splited,

    .h = 0.7,
    .Obh2 = 0.02,
    .Och2 = 0.12,

    .w0 = NAN,


    .n_s = 0.96,
    .alpha_s = 0.,
    .k_pivot = 0.05,
    .sigma8 = NAN,
    .delta_zeta = 0.,
    .A_s = 2.0e-09,
    //对应的CLASS算出来的.sigma8 = 0.844482,

    .tau_reio = 0.089,
};


int main()
{
    int status=0;
	omp_set_nested(1);
    //用于设置是否允许OpenMP进行嵌套并行，默认的设置为false。
	time_t start,end;
    start =time(NULL);//or time(&start);
    //创建 fiducial model
    clucu_cosmology fiducial=clucu_cosmology_create(config_default,param_collapse,survey_csst2);
    //对 fiducial model 稍作修改
    fiducial.config.matter_power_spectrum_method = clucu_boltzmann_class_Tk,
    fiducial.config.mass_observable_method=clucu_my_csst;
    fiducial.survey.mass_ob_min_h=0.836*pow(10.,14.) ;


    //fiducial.param.Onuh2_mass = 0.6/93.14;

    //.m_nu = {0.02,0.02,0.02},
    fiducial.GISDB = 1;
	strcpy(fiducial.classname,"collapse");
	strcpy(fiducial.name,"nu00");
	fiducial.survey.redshift_sigma0 = 0.001;
	fiducial.runned_class = false;

    double mnu=0.06;
    double Onuh2_1 = density_WDM(mnu,pow(4./11.,1./3.)*2.7255)/clucu_constants.RHO_CRITICAL /mnu*94.;
    double Onuh2_2 = density_WDM(mnu,0.71611*2.7255)/clucu_constants.RHO_CRITICAL /mnu*93.14;
    //double Onuh2_2 = density_WDM(0.059,0.71611*2.725)/clucu_constants.RHO_CRITICAL;
    //0.71611
    //0.71599

	//创建一系列model
	int Np=15;
	int Nid = Np*2+1;
	double *partial_param=Create1Grid(Np);
	clucu_cosmology *cosmo =(clucu_cosmology *) malloc(sizeof(clucu_cosmology)*Nid);
    
	cosmo[0]=fiducial;
	cosmo[0].cosmo_id = 0;
	cosmo[0].class_id = 0;
	InitializeCosmo(&cosmo[0]);
    //末尾多加几个200
    int nztemp = 100;
    cosmo->data.spline_Nz = nztemp+4;
    cosmo->data.spline_z = Create1Grid(cosmo->data.spline_Nz);
    double *ztemp = clucu_linlog_spacing(0.,2.,198.,21,80);
    for(int iz=0;iz<cosmo->data.spline_Nz;iz++)
        cosmo->data.spline_z[iz]=ztemp[iz];
    cosmo->data.spline_z[nztemp]=199.;
    cosmo->data.spline_z[nztemp+1]=200.;
    cosmo->data.spline_z[nztemp+2]=201;
    cosmo->data.spline_z[nztemp+3]=202;

    clucu_compute_background(&cosmo[0],clucu_background_class_label,&status);
    clucu_compute_power(&cosmo[0],&status);

    //clucu_compute_logsigma(&cosmo[0],&status);
    printf("#%d:log sigma --over\n",0);  
    //clucu_sigma(&cosmo[0],1.0e15,1.,&status); 
    clucu_collapse(&cosmo[0],&status);

	
}


