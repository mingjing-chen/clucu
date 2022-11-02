#pragma once

#include "clucu_survey.h"

// Params for hmf() integrand
typedef struct {
    //double *func;
    double (*func)(double x,void *p);
    clucu_cosmology *cosmo;
    double ln_mass_ob;
    double ln_mass_ob1;
    double ln_mass_ob2;
    double ln_mass_ob11;
    double ln_mass_ob22;
    double z_ob;
    double z_true;
    double chi;
    double ell;
    double k;
    int CEN_massob;
    int if_int_masstrue;
    int if_int_ztrue;
    int *status;
} int_pars;
double func_integrand_mass_z(clucu_cosmology *cosmo,double (*func)(double x,void *p),int if_int_masstrue,int if_int_ztrue,double ln_mass_ob1,double ln_mass_ob2,int CEN_massob,double z_ob1,double z_ob2,int CEN_zob,int *status);
double func_integrand_mass_z_ell(clucu_cosmology *cosmo,double (*func)(double x,void *p),int if_int_masstrue,int if_int_ztrue,double ell,double ln_mass_ob1,double ln_mass_ob2,int CEN_massob,double z_ob1,double z_ob2,int CEN_zob,int *status);
double func_integrand_mass_z_k(clucu_cosmology *cosmo,double (*func)(double x,void *p),int if_int_masstrue,int if_int_ztrue,double k,double ln_mass_ob1,double ln_mass_ob2,int CEN_massob,double z_ob1,double z_ob2,int CEN_zob,int *status);
