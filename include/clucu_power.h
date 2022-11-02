#pragma once

#include "clucu_core.h"
#include "clucu_background.h"
#include "clucu_power_utils.h"



typedef struct eh_struct {
  double rsound;
  double zeq;
  double keq;
  double zdrag;
  double kSilk;
  double rsound_approx;
  double th2p7;
  double alphac;
  double alphab;
  double betac;
  double betab;
  double bnode;
  int wiggled;
} eh_struct;


double tsqr_BBKS(clucu_cosmology *cosmo, double k) ;
eh_struct* clucu_eh_struct_new(clucu_cosmology *cosmo, int wiggled);
double tsqr_EH(clucu_cosmology *cosmo,eh_struct *eh,double k);
double clucu_transfer(clucu_cosmology *cosmo,double k,int *status);

double clucu_transfer_label(clucu_cosmology *cosmo,double k, double z,clucu_species_x_label label);


void clucu_compute_linpower_eh(clucu_cosmology *cosmo, int wiggled, int *status);
void clucu_compute_linpower_bbks(clucu_cosmology *cosmo, int *status);
void clucu_compute_power_class_Pk(clucu_cosmology *cosmo, int *status);
void clucu_compute_power_class_Tk(clucu_cosmology *cosmo, int *status);
void clucu_compute_power(clucu_cosmology *cosmo, int *status);

double clucu_power(clucu_cosmology *cosmo,double k, double z,int *status);

double clucu_sigmaR(clucu_cosmology *cosmo,double R, double z,clucu_f2d_t *psp, int *status);
double clucu_sigma8(clucu_cosmology *cosmo, clucu_f2d_t *psp, int *status) ;

void clucu_compute_logsigma(clucu_cosmology *cosmo,int *status);

double clucu_sigma(clucu_cosmology *cosmo, double halo_mass, double redshift,int *status);
double clucu_dlnsigma_dlnm(clucu_cosmology *cosmo, double halo_mass, double redshift,int *status);

double M_to_R(clucu_cosmology *cosmo,double halomass,int *status);
