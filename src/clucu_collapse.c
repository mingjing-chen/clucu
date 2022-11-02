#define print_headers 1
#define MsuntoMpc 4.80e-20
#include "clucu_collapse.h"
//TODO:想一下，这里的sigma，需不需要先插值。插值需要时间，但如果sigma算的次数很多，还是插值更好。
int clucu_collapse(clucu_cosmology *cosmo,int *status){
    double H0 = clucu_hubble_parameter(cosmo,0,status);//1/Mpc
    double h = cosmo->param.h;
    double Omega_cb = cosmo->param.Omega_b + cosmo->param.Omega_c;
    double Omega_r = cosmo->param.Omega_nu_rel+cosmo->param.Omega_g;
//create k list, survey.k
    int Nk=200;
    double *log10k = clucu_linear_spacing(-5.,1.,Nk);
    //创建M
//    int Nm = cosmo->spline_param.LOG10M_SPLINE_NM;
//    double *log10M = clucu_linear_spacing(14.,cosmo->spline_param.LOG10M_SPLINE_MAX,Nm);
    int Nm = 2;
    double *log10M = clucu_linear_spacing(12.,18.,Nm);
    double **short_over_long = Create2Grid(Nm,Nk);
    double **b_L = Create2Grid(Nm,Nk);
    double *delta_crit = Create1Grid(Nm);


    double zip1 = z_initial - 1.;
    double zip2 = z_initial + 1.;

    

//这里可以改成其他light relic, 需要改温度
    double T0_nu = clucu_constants.TNCDM * clucu_constants.T_CMB;
	double rho_nu1_zi=density_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_initial));
	double rho_nu1_z0=density_WDM(cosmo->param.m_nu[0],T0_nu * (1.+0.));
    double rho_nu2_zi=density_WDM(cosmo->param.m_nu[1],T0_nu * (1.+z_initial));
	double rho_nu2_z0=density_WDM(cosmo->param.m_nu[1],T0_nu * (1.+0.));
    double rho_nu3_zi=density_WDM(cosmo->param.m_nu[2],T0_nu * (1.+z_initial));
	double rho_nu3_z0=density_WDM(cosmo->param.m_nu[2],T0_nu * (1.+0.));

	double rho_nu1_ratio_zi=0.;
    double rho_nu2_ratio_zi=0.;
    double rho_nu3_ratio_zi=0.;
	if(cosmo->param.m_nu[0]>0)
		rho_nu1_ratio_zi=rho_nu1_zi/rho_nu1_z0;
    if(cosmo->param.m_nu[1]>0)
		rho_nu2_ratio_zi=rho_nu2_zi/rho_nu2_z0;
	if(cosmo->param.m_nu[2]>0){
		rho_nu2_ratio_zi=rho_nu3_zi/rho_nu3_z0;
	}

    

///////////////////////////////////////////////////////////////////////////////////
////    we set the R(z) array that we will use to find extra species clustering////
///////////////////////////////////////////////////////////////////////////////////
//如果需要考虑中微子的clustring，则需要计算一个中微子团的质量
    const int Nz_solution=400;//number of z in solution of R(t). LOGSPACED.
	const double logz_solution_min=log(z_initial);
	const double logz_solution_max=log(z_collapse/2);//as in collapse code.
	const double dlogz_solution=(logz_solution_max-logz_solution_min)/(Nz_solution-1);

	double *Rhalo_solution=Create1Grid(Nz_solution);//radius of halo at different zs.
	double *Mnu1_solution=Create1Grid(Nz_solution); //neutrino 1 mass within R_halo, needed to find collapse.
    double *Mnu2_solution=Create1Grid(Nz_solution);; //neutrino 1 mass within R_halo, needed to find collapse.
	double *Mnu3_solution=Create1Grid(Nz_solution); //neutrino 1 mass within R_halo, needed to find collapse.
    //TODO:很奇怪，M会有整齐的、非零的初值。平时用创建动态数组的函数，初值都是1e-316
    //所以我们手动赋初值0
    for(int i=0;i<Nz_solution;i++){
        Mnu1_solution[i]=0.;
        Mnu2_solution[i]=0.;
        Mnu3_solution[i]=0.;
    }

	//z array logspaced between zi and zcollapse/2, to match collapse code.

//////////////////////////////////////////////////////
////    we set initial conditions									////
/////////////////////////////////////////////////////
	double OmL_i=cosmo->param.Omega_l;
	double OmM_i=Omega_cb*pow(1.+z_initial,3.);
	double OmR_i=Omega_r*pow(1.+z_initial,4.);
    double Omnu1_i=cosmo->param.Omega_nu_mass[0]*rho_nu1_ratio_zi;
    double Omnu2_i=cosmo->param.Omega_nu_mass[1]*rho_nu2_ratio_zi;
    double Omnu3_i=cosmo->param.Omega_nu_mass[2]*rho_nu3_ratio_zi;
    const double Hi=H0*sqrt(OmL_i + OmM_i + OmR_i + Omnu1_i + Omnu2_i + Omnu3_i);//H(zi) in Mpc-1
    //或者直接调用函数
    //const double Hi=clucu_hubble_parameter(cosmo,zi,status);


    //3个红移处的sigma8
    //注意，检查log_pkz是不是P_cb
	const double sigma8_0 = clucu_sigmaR(cosmo,8/cosmo->param.h,0.,cosmo->data.log_pkz,status);
	const double sigma8_collapse = clucu_sigmaR(cosmo,8/cosmo->param.h,z_collapse,cosmo->data.log_pkz,status);
	const double sigma8_i = clucu_sigmaR(cosmo,8/cosmo->param.h,z_initial,cosmo->data.log_pkz,status);


//we set delta_short.
 	const double delta_short_min_constant=delta_short_constant_centroid/(2.0+boost_initial_conditions); //initial conditions, to do bipartition. these are the benchmark.
 	const double delta_short_max_constant=delta_short_constant_centroid*(2.0+boost_initial_conditions);

	const double tolerance=0.0002/max(min(precision_scale,10.),1.); //relative error before we call it converged. 3e-4 should guarantee 0.1% in bias.
	const double tolerance_ini= tolerance + (do_clustering>0)*tolerance*4.0; //first try does not have to be as precise (if we do clustering), since it is used to find Mnu collapse only.
	const double tolerance_z = z_collapse * 0.1; //this is just to make sure that we are not artificially converging to a "bisection" \delta_crit if our initial conditions are bad.


///////////////////////////////////////////////////////////////////////////////////
////    we iterate over klong array defined above, as well as masss		      ////
/////////////////////////////////////////////////////////////////////////////////
    const int Nk_calc[N_delta_long] = {1, Nk}; //what k_longs we calculate over. For delta_long=0 only one is needed, otherwise one per kmode.
    //N_delta_long is defined in common.h

    //创建delta_L
    double *delta_long_list = Create1Grid(N_delta_long);
    const double delta_long_min=0.; //the maximum value is specified in common.h for convenience to change it.
    const double delta_long_step=(delta_long_max-delta_long_min)/(N_delta_long-1.);
    for(int i=0;i<N_delta_long;i++){
        delta_long_list[i]=delta_long_min+delta_long_step*i;
    }

    //创建delta_S_coll
    ////these variables are to save the delta_L and delta crit to files
    double **delta_short_collapse = Create2Grid(N_delta_long,Nk); //delta_short (initial) to collapse for each delta_long and k_long.
    double **delta_short_crit = Create2Grid(N_delta_long,Nk); //delta_short_critical, extrapolated to z_collapse.
	double **delta_long_collapse = Create2Grid(N_delta_long,Nk); //delta_long, extrapolated to z_collapse.

    
    //开始循环
    for(int im=0;im<Nm;im++){

        double mass = pow(10.,log10M[im]);//mass in M_sun
        double mass_Mpc = mass * MsuntoMpc;//mass in Mpc (*G)
        if(debug_mode>=0)  printf("M_halo = %.1le Msun \n", mass);

        //average R_i
       // double Ribar =  M_to_R(cosmo,mass,status);//R in Mpc.
        double Ribar = pow(Omega_cb*pow(1.+z_initial,3.)/2.*H0*H0/mass_Mpc,-1./3); //R in Mpc.
        //TODO: 这里是zi时的半径！
        double R0 = pow(Omega_cb/2.*H0*H0/mass_Mpc,-1./3);//算sigma时用这个半径。

        //sigma(M) at zi, z)collapse, and derivative (and z=0 just in case), for initial conditions of ODE. (option_initial_velocity=0)
        // double sigmaM_0 = clucu_sigma(cosmo,mass,0.,status);
        // double sigmaM_collapse = clucu_sigma(cosmo,mass,z_collapse,status);
        // double sigmaM_i = clucu_sigma(cosmo,mass,zi,status);
        // double sigmaM_zip1 = clucu_sigma(cosmo,mass,zip1,status);
        // double sigmaM_zip2 = clucu_sigma(cosmo,mass,zip2,status);
        // double d_sigmaM_i = (sigmaM_zip2 - sigmaM_zip1)/(zip2-zip1);
        double sigmaM_0 = clucu_sigmaR(cosmo,R0,0.,cosmo->data.log_pkz,status);
        double sigmaM_collapse = clucu_sigmaR(cosmo,R0,z_collapse,cosmo->data.log_pkz,status);
        double sigmaM_i = clucu_sigmaR(cosmo,R0,z_initial,cosmo->data.log_pkz,status);
        double sigmaM_zip1 = clucu_sigmaR(cosmo,R0,zip1,cosmo->data.log_pkz,status);
        double sigmaM_zip2 = clucu_sigmaR(cosmo,R0,zip2,cosmo->data.log_pkz,status);
        double d_sigmaM_i = (sigmaM_zip2 - sigmaM_zip1)/(zip2-zip1);


        if(debug_mode>0){
            printf("IC_1:%.3le\n", d_sigmaM_i/sigmaM_i);
            printf("Ribar=%.6le \n", Ribar);
            printf("R0=%.6le \n", R0);

            printf("We will run some checks: \n");
            printf("@z=%.1le \t sigma_M= %.6le \n",0., sigmaM_0);
            printf("@z=%.1le \t sigma_8= %.6le \n",0., sigma8_0);
            printf("@z=%.1le \t sigma_M= %.6le \n", z_collapse, sigmaM_collapse);
            printf("@z=%.1le \t sigma_8= %.6le \n", z_collapse, sigma8_collapse);
            printf("@z=%.1le \t sigma_M= %.6le \n",z_initial, sigmaM_i);
            printf("@z=%.1le  d sigma_M/dz= %.6le \n",z_initial, d_sigmaM_i);
            printf("@z=%.1le \t sigma_8= %.6le \n",z_initial, sigma8_i);
        }

        #ifdef _OPENMP
        omp_set_nested(1); //we need to nest the two parallel loops since the limits on the inner for() depend on the outer one, so the scheduler does not know how to do collapse(2).
        #endif
        
        // int thread_num = (int)(omp_get_num_procs()*THREAD_RATE/cosmo->thread_num);
        // if(CHECK) thread_num=1;
        // thread_num=1;
        // omp_set_num_threads(thread_num);
        #pragma omp parallel for
        for(int i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
            #pragma omp parallel for
            for(int ik=0;ik<Nk_calc[i_delta_long];ik++){
                //i_delta_long时，delta_long=0，k_long只用算1个, Nk_calc[0]=1
                double delta_long=delta_long_list[i_delta_long]; //long-wavelength perturbation for CDM+b , for photons , and for massless neutrinos.
                double k_long=pow(10.,log10k[ik]);
                double Ti_cb = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_cb_label);
                double dTi_cb =(clucu_transfer_label(cosmo,k_long,zip2,clucu_species_cb_label)-clucu_transfer_label(cosmo,k_long,zip1,clucu_species_cb_label))/(zip2-zip1);

                double ddelta_long_dt=delta_long*(-Hi*(1.+z_initial))*dTi_cb/Ti_cb; //d delta_long/dt, calculated as dT/dt*1/T*delta_long.

                double delta_short_min=delta_short_min_constant; //we restart the initial conditions, and redo the procedure.
                double delta_short_max=delta_short_max_constant;
                double delta_short=0.0,ddelta_short_dt; //private copies for each run
                double Ri, Ridot, Rpi, z_coll_iteration=0.0; //initial conditions and z of collapse of each iteration.

                if(debug_mode>1){
                    printf("deltashort_min=%.2le, deltashort_max=%.2le \n", delta_short_min, delta_short_max);
                }

                for(;delta_short_max/delta_short_min-1.>tolerance_ini;){
                    delta_short=(delta_short_max+delta_short_min)/2.0; //bisection, we update the min or max value.
                    ddelta_short_dt=delta_short*(-Hi*(1.+z_initial))*d_sigmaM_i/sigmaM_i; //since: 		ddelta_dt/delta = dsigma_dt/sigma = -H*(1+z) dsigma_dz/sigma

                    //these are the i.c. we will feed the function find_z_collapse_XXX.
                    Ri= Ribar*(1.0-1./3.*(delta_short+delta_long)); //R_i in Mpc.
                    Ridot= Ri*Hi*(1.0-1./3.*(ddelta_short_dt+ddelta_long_dt)/Hi); //dR/dt in Mpc/Mpc
                    Rpi= -Ridot/((1.0+z_initial)*Hi); //dR/dz in Mpc.

                    if(debug_mode > 0){
                        printf("Ri=%.6le, dRi/dz=%.6le, delta_short=%.6le \n", Ri, Rpi, delta_short);
    //					printf("i=%d, T_cdm_i=%.2le, T_g=%.2le, T_nu=%.2le \n",ik, Ti_cb[ik], transfer_gamma_klong[ik][33], transfer_nu_massless_klong[ik][33]);
                    }

                    //解方程，找塌缩红移。目前只考虑两种情况。
                    //we now call the collapse routine relevant for each case:
                    if(cosmo->param.N_nu_mass == 0){//nothing extra
                        if (debug_mode>0) printf("Collapse with nothing extra \n");
                        z_coll_iteration=find_z_collapse_nothing(cosmo, Ri, Rpi, delta_long, k_long, mass);
                        printf("z_coll=%.6e\n",z_coll_iteration);
                    }
                    else if(cosmo->param.N_nu_mass == 1){//only mnu1
                        if (debug_mode>0) printf("Collapse with nu1 \n");
                        z_coll_iteration=find_z_collapse_1nu(cosmo, Ri, Rpi, delta_long, k_long, mass,
                                        Nz_solution, Rhalo_solution, Mnu1_solution);
                    }
                    else if(cosmo->param.N_nu_mass == 2){//only mnu1
                        if (debug_mode>0) printf("Collapse with nu2 \n");
                        z_coll_iteration=find_z_collapse_2nu(cosmo, Ri, Rpi, delta_long, k_long, mass,
                                        Nz_solution, Rhalo_solution, Mnu1_solution,Mnu2_solution);
                    }
                    else if(cosmo->param.N_nu_mass == 3){//only mnu1
                        if (debug_mode>0) printf("Collapse with nu3 \n");
                        z_coll_iteration=find_z_collapse_3nu(cosmo, Ri, Rpi, delta_long, k_long, mass,
                                        Nz_solution, Rhalo_solution, Mnu1_solution,Mnu2_solution,Mnu3_solution);
                    }
                    

                    if(z_coll_iteration>z_collapse){ //collapses too quickly
                        delta_short_max=delta_short;
                    }
                    else{ 							//collapses too slowly
                        delta_short_min=delta_short;
                    }

                }	//end of delta_short loop

                if(debug_mode > 0){
                        printf("delta_short=%le, delta_long=%le (with k_long=%le) \n",delta_short, delta_long, k_long);
                        printf("z_coll_iteration= %le \n",z_coll_iteration);
                    }


                if(abs(z_coll_iteration-z_collapse)>tolerance_z){ //collapses too quickly
                    printf("Not converged to z_collapse. Initial conditions too narrow, make boost_initial_conditions bigger. \n");
                    printf("z_coll_it=%.3le and z_collapse=%.3le \n",z_coll_iteration,z_collapse);
                }

                delta_short_collapse[i_delta_long][ik]=delta_short;


                //			double percentage=1.*(ik+1.+cosmo->Nk*i_delta_long)/N_delta_long/cosmo->Nk;//not really a percentage, a one-centage.
                //			printProgress(percentage);
                //	we do not use percentage bar if parallel. It looks ugly.
            }
        }
        

        // since the delta_long=0 does not depend on k we can just copy it for all ks.
        for(int i_delta_long=0;i_delta_long<1;i_delta_long++){
            for(int ik=1;ik<Nk;ik++){
                delta_short_collapse[i_delta_long][ik]=delta_short_collapse[i_delta_long][0];
            }
        }

        if (debug_mode > 1){
            for (int i=0;i<Nz_solution;i++){
                printf("i=%d, z=%.1le, R(z)/Mpc=%.1le \n",i, logz_solution_min*exp(dlogz_solution*i),Rhalo_solution[i]);
            }
        }
        //TODO: relicfast还考虑了neutrino clustering，我们暂不考虑

        //////////////////////////////////////////////////////////////////////////////
        //// we now extrapolate  to find\delta_crit to the redshift of collapse	/////
        ////////////////////////////////////////////////////////////////////////////
//        #pragma omp parallel for collapse(2)
        for(int i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
            for(int ik=0;ik<Nk;ik++){
                double delta_long=delta_long_list[i_delta_long];
                double k_long=pow(10,log10k[ik]);

                double T_z_collapse_klong = clucu_transfer_label(cosmo,k_long,z_collapse,clucu_species_cb_label);
                double Ti_cb = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_cb_label);

                delta_long_collapse[i_delta_long][ik]=delta_long*T_z_collapse_klong/Ti_cb;	//this is right as well //			printf("klong=%.1le, delta_long_collapse/delta_long= %.3le \n", k_long, delta_long_collapse[i_delta_long][ik]/delta_long);

                delta_short_crit[i_delta_long][ik]=delta_short_collapse[i_delta_long][ik]*sigmaM_collapse/sigmaM_i;

            }
        }
        for(int ik=0;ik<Nk;ik++){
            short_over_long[im][ik] = (delta_short_crit[1][ik]-delta_short_crit[0][ik])/(delta_long_collapse[1][ik]-delta_long_collapse[0][ik]);
        }
        delta_crit[im]= delta_short_crit[0][0];
        ////////////////////////////////////////////
        ////and we save the results to files.	/////
        //////////////////////////////////////////
        int lengthname=200;
        char *filename; //To open files.
        filename=(char *)malloc((lengthname+1)*sizeof(char));
        FILE *fp;

        lengthname=sprintf(filename,"../output/%s",cosmo->name);
        mkdir("../output", 0755);
        mkdir(filename, 0755);

        lengthname=sprintf(filename,"../output/%s/delta_initial_z%.2f_M%.2f_Nk%d.dat",cosmo->name,z_collapse, log10(mass),Nk); //we save the delta_crit extrapolated to z=0.5
		fp=fopen(filename,"w");
		if(print_headers!=0){
			fprintf(fp,"delta_long_initial   k_long[1/Mpc]   delta_short_initial \n");
		}

		for(int i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
			double delta_long=delta_long_list[i_delta_long];
			for(int ik=0;ik<Nk;ik++){
				double k_long=pow(10.,log10k[ik]);
				fprintf(fp,"%le   %le   %le \n", delta_long, k_long, delta_short_collapse[i_delta_long][ik]);
			}
		}
		fclose(fp);


		lengthname=sprintf(filename,"../output/%s/delta_crit_z%.2f_M%.2f_Nk%d.dat",cosmo->name,z_collapse, log10(mass),Nk); //we save the delta_crit extrapolated to z=0.5
		fp=fopen(filename,"w");
		if(print_headers!=0){
			fprintf(fp,"delta_long   k_long[1/Mpc]   delta_crit \n");
		}
		for(int i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
			for(int ik=0;ik<Nk;ik++){
				double k_long=pow(10.,log10k[ik]);
				fprintf(fp,"%le   %le   %le \n", delta_long_collapse[i_delta_long][ik], k_long, delta_short_crit[i_delta_long][ik]);
			}
		}
		fclose(fp);

        //算一下bias
        double delta_crit_sim=0.;
		int counter=0;
		for(int i_delta_long=0;i_delta_long<N_delta_long;i_delta_long++){
			for(int i_klong=0;i_klong<Nk;i_klong++){
				counter++;
				delta_crit_sim+=(delta_long_collapse[i_delta_long][i_klong]+delta_short_crit[i_delta_long][i_klong]);
			}
		}
		delta_crit_sim/=(1.*counter);

        double delta_ref=1.686;
		double a_MICE=1.37*pow(1+z_collapse,-0.15);
		double b_MICE=0.3*pow(1+z_collapse,-0.084);
		double c_MICE=1.036*pow(1+z_collapse,-0.024);

		double HMF_MICE= -2. * c_MICE * delta_crit_sim/(delta_ref*delta_ref *sigmaM_collapse*sigmaM_collapse) +
						a_MICE/delta_crit_sim/(1. + b_MICE*pow(sigmaM_collapse*delta_ref/delta_crit_sim,a_MICE));
        
        double a_ST=0.707; //updated from original 0.707, from Ref.1005.2239 Table 3
		double p_ST=0.3;
		double HMF_ST= (1-a_ST*pow(delta_crit_sim/sigmaM_collapse,2.))/delta_crit_sim - (2*p_ST/delta_crit_sim)/(1.+pow(a_ST*pow(delta_crit_sim/sigmaM_collapse,2.),p_ST));    
        
        double HMF = HMF_ST;
        if(debug_mode>0){
			printf("delta_crit_sim=%le, sigmaM_coll=%le ,nu=%le  ,HMF=%le\n",delta_crit_sim,sigmaM_collapse,delta_crit_sim/sigmaM_collapse,HMF_ST);
		}
        for(int ik=0;ik<Nk;ik++){
			b_L[im][ik] = HMF * short_over_long[im][ik];
            if(debug_mode>0){
			printf("k/h=%.1le, HMF=%le ,derivative=%le  ,bL=%le\n",pow(10.,log10k[ik])/h,HMF,short_over_long[im][ik],b_L[im][ik]);
		    }
		}
        lengthname=sprintf(filename,"../output/%s/bL_z%.2f_M%.2f_Nk%d.dat",cosmo->name,z_collapse, log10(mass),Nk); //we save the delta_crit extrapolated to z=0.5
		fp=fopen(filename,"w");
		if(print_headers!=0){
			fprintf(fp,"k_long[1/Mpc]   short_over_long  bL  HMF=%.2e\n",HMF);
            fprintf(fp,"delta_crit_sim=%le, sigmaM_coll=%le ,nu=%le  ,HMF=%le\n",delta_crit_sim,sigmaM_collapse,delta_crit_sim/sigmaM_collapse,HMF);
		}
        for(int ik=0;ik<Nk;ik++){
            double k_long=pow(10.,log10k[ik]);
            fprintf(fp,"%le   %le   %le \n", k_long, short_over_long[im][ik],b_L[im][ik]);
        }
		fclose(fp);

        lengthname=sprintf(filename,"../output/%s/SoverL_z%.2f_M%.2f_Nk%d.dat",cosmo->name,z_collapse, log10(mass),Nk); //we save the delta_crit extrapolated to z=0.5
		fp=fopen(filename,"w");
		if(print_headers!=0){
			fprintf(fp,"k_long[1/Mpc]   SoverL \n");
		}
        for(int ik=0;ik<Nk;ik++){
            double k_long=pow(10.,log10k[ik]);
            fprintf(fp,"%le   %le \n", k_long, short_over_long[im][ik]);
        }
		fclose(fp);

        
    
    }//M循环

    int lengthname=200;
    char *filename; //To open files.
    filename=(char *)malloc((lengthname+1)*sizeof(char));
    FILE *fp;

    lengthname=sprintf(filename,"../output/%s",cosmo->name);
	mkdir("../output", 0755);
	mkdir(filename, 0755);

    lengthname=sprintf(filename,"../output/%s/SoverL.dat",cosmo->name); //we save the delta_crit extrapolated to z=0.5
    fp=fopen(filename,"w");
    if(print_headers!=0){
        fprintf(fp,"k_long[1/Mpc]\\Msun   SoverL \n");
    }
    save_matrix(fp,short_over_long,log10M,log10k,Nm,Nk);
	fclose(fp);

    lengthname=sprintf(filename,"../output/%s/delta.dat",cosmo->name); //we save the delta_crit extrapolated to z=0.5
    fp=fopen(filename,"w");
    if(print_headers!=0){
        fprintf(fp,"Msun   delta \n");
    }
    for(int im=0;im<Nm;im++){
        fprintf(fp,"%le   %le\n",log10M[im], delta_crit[im]);
    }
	fclose(fp);
    


    //插值 Ddelta_s_coll/Ddelta_l_coll(M,k)
    //插值(delta_l=0时)delta_crit(M)
    double *y1 = malloc(sizeof(double)*Nm*Nk);
    for(int im=0;im<Nm;im++){
        for(int ik=0;ik<Nk;ik++){
            y1[ik*Nm + im] = short_over_long[im][ik];
        }
    }
    gsl_spline2d *spl1 = NULL;
    spl1 = gsl_spline2d_alloc(gsl_interp2d_bicubic, Nm, Nk);
    if (spl1 == NULL) 
    {
        *status = CLUCU_ERROR_MEMORY;
        CLUCU_RAISE_WARNING(*status,"error allocating 2D spline");
    }
    int s2dstatus=gsl_spline2d_init(spl1, log10M, log10k, y1, Nm, Nk);
    if (s2dstatus) 
    {
        *status = CLUCU_ERROR_SPLINE;
        CLUCU_RAISE_WARNING(*status,"error initializing spline");
    }

    gsl_spline *spl2 = gsl_spline_alloc(gsl_interp_cspline, Nm);
    if (gsl_spline_init(spl2,log10M,delta_crit,Nm))
    {
        *status = CLUCU_ERROR_SPLINE;
        CLUCU_RAISE_WARNING(*status,"Error creating delta_crit spline");
    }

    if (*status == 0) {
        cosmo->GISDB = true;
        cosmo->data.spline_SoverL = spl1;
        cosmo->data.spline_delta_crit = spl2;
    }
    else{
        gsl_spline2d_free(spl1);
        gsl_spline_free(spl2);
    }
}
double clucu_SoverL(clucu_cosmology *cosmo,double mass,double k,int *status)
{
    if(!cosmo->GISDB){
        *status = CLUCU_ERROR_DISTANCES_INIT;
        CLUCU_RAISE_WARNING(*status,"SoverL splines have not been precomputed!");
        return NAN;
    }

    double result;
    int gslstatus = gsl_spline2d_eval_e(cosmo->data.spline_SoverL, log10(mass),log10(k), NULL,NULL, &result);
    if(gslstatus != GSL_SUCCESS){
        CLUCU_RAISE_GSL_WARNING(gslstatus,"mass=%.2e,k=%.2e is outside interpolation range [%.2e,%.2e][%.2e,%.2e]",
                                    mass,k,
                                    cosmo->data.spline_SoverL->interp_object.xmin,cosmo->data.spline_SoverL->interp_object.xmax,
                                    cosmo->data.spline_SoverL->interp_object.ymin,cosmo->data.spline_SoverL->interp_object.ymax);
        *status = gslstatus;
        return NAN;
    }
     return result;
}
double clucu_delta_crit(clucu_cosmology *cosmo,double mass,double k,int *status)
{
    if(!cosmo->GISDB){
        *status = CLUCU_ERROR_DISTANCES_INIT;
        CLUCU_RAISE_WARNING(*status,"delta_crit splines have not been precomputed!");
        return NAN;
    }

    double result;
    int gslstatus = gsl_spline_eval_e(cosmo->data.spline_delta_crit, log10(mass),NULL, &result);
    if(gslstatus != GSL_SUCCESS){
        CLUCU_RAISE_GSL_WARNING(gslstatus,"mass=%.2e is outside interpolation range [%.2e,%.2e]",
                                    mass,cosmo->data.spline_delta_crit->interp->xmin,cosmo->data.spline_delta_crit->interp->xmax);
        *status = gslstatus;
        return NAN;
    }
     return result;
}
double find_z_collapse_nothing
(clucu_cosmology * const cosmo, const double Ri, const double Rpi, const double delta_long,const double k_long,const double mass){
// This function returns the z at which an overdensity collapses. global cosmological parameters assumed.
// In the presence of a long-wavelength PHOTON, MASSLESS and 1 MASSIVE NU (OR OTHER EXTRA) perturbation.
//we read the instantaneous transfer function from CAMB at each redshift.
//and we save R(z) to perform BKT calculation of neutrino clustering, or we input Mnu(z) if we already calculated it.
//we only modify R_solution and Mnu_solution, the rest are read-only.

    int status=0;
    double H0 = clucu_hubble_parameter(cosmo,0.,&status);
    double h = cosmo->param.h;
    double T0_nu = clucu_constants.TNCDM * clucu_constants.T_CMB;
    double mass_Mpc = mass * MsuntoMpc;
    const double T_matter = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_cb_label);

	const int precision = (int) fmax(fmin(precision_normalization,10),1);//1<=precision<=20
	const long npoints = precision*30000; //number of points for the ODE solution to converge, ~10^4 logspaced in z.

	double zf_code=z_collapse/2.0;
	if(zf_code<z_collapse_min){
		printf("Warning, for z_collapse = 0 you need to change z binning to linear in collapse code \n");
    zf_code=z_collapse_min;
	}
	const double zstep_log=(log(zf_code/z_initial))/(npoints-1.); //log-spaced in z.
	double zstep_lin; //for integration of ODE, linear step.


	double R1, R2; //R and R_next.
	double Rp1, Rp2; //derivative of R and next value.
    double Rpp1, Rpp2; //d^2R(z)/dz^2
	double R1tilde, Rp1tilde; //the "guesses" of R and Rprime for  Heun's methood
	R1=R2=Ri; //initial conditions
	Rp1=Rp2=Rpi;

	//the nu_massive can be any new species, just make sure that the Omegas are well defined here
    // if((cosmo->Omega_extra>0) && (cosmo->mnu1>0)){
    // printf("Using wrong collapse function, XXX_1nu only accepts one species. \n");
    // }

//set the initial
//other,(we assme c_s^2=w=1/3 for photons and massless nu)
	double OmGbar= cosmo->param.Omega_g * pow(1.+z_initial,4.); //photon
	double Omnu_masslessbar= cosmo->param.Omega_nu_rel * pow(1.+z_initial,4.); //massless nu
	double OmRbar = OmGbar + Omnu_masslessbar; //all radiation (massless nu + photon)
	double OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_initial,3.); //matter
	double OmL = cosmo->param.Omega_l; //Omega_Lambda(z), we take it as z-independent.
    double T_gamma = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_g_label);
	double T_nu = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_ur_label);
	double OmG = OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
	double Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
	double OmR, OmR2 = OmG + Omnu_massless;
//H
    double H = H0 * sqrt(OmL + OmM + OmRbar); //H(zi), 1/Mpc,TODO，这里求H可以直接用函数代替
	


//for the next values, to do Heun's method. We assume Omega_L does not change. Symbol 2 denotes next step.
//we also get initial dE, for which we need H2 (next z)
	double z_next = z_initial * exp(zstep_log);//first step
//and the w and c_s^2 of nus at first step we calculate.
//other
    OmGbar = cosmo->param.Omega_g * pow(1.+z_next,4.);
	Omnu_masslessbar = cosmo->param.Omega_nu_rel * pow(1.+z_next,4.);
	OmRbar = OmGbar + Omnu_masslessbar;
	OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_next,3.);
	OmL = cosmo->param.Omega_l;
//H2
    double H2 = H0 * sqrt(OmL + OmM + OmRbar); //H(next), 1/Mpc
	double dE, dE2 = (H2-H)/H2/(z_initial*zstep_log);

//printf("H=%.6e,  H2=%.6e,  dE2=%.2e\n",H,H2,dE2);

///////////////////////////////////////////////////////////
////    here we solve for the collapse				  ////
/////////////////////////////////////////////////////////
    double z;
	for(long i=1; i<npoints-1 && R2>0.; i++){
        //we update the values from the previous step
        R1 = R2;
        Rp1 = Rp2;
        //including all z-dependent quantities, for Heun's method. Remember 2 means next step.
        z=z_next;
        H=H2;
        dE=dE2;
        OmR=OmR2;



		//current redshift is updated at the end
		z_next *= exp(zstep_log); //next redshift
		zstep_lin = z_next - z; //linear step, more precise than z*zstep_log and cheap.



	    //we find R''(z). Only depends on z and other variables calculated at the previous z.
		Rpp1 = - Rp1*(1./(1+z) + dE)
		- mass_Mpc/R1/R1/pow(H*(1.+z),2.)
		- R1/2.*(2*OmR - 2*OmL)*H0*H0/pow(H*(1.+z),2.); //R''(R1,t1)
//printf("Rp1=%.2e,  dE=%.2e  z=%.2e,  R1=%.2e,  mass_Mpc=%.2e,  OmR=%.2e,  OmL=%.2e,  H0=%.2e  H=%.2e\n",Rp1,dE,z,R1,mass_Mpc,OmR,OmL,H0,H);
        
        //we update all z-dependent stuff to the next z
        //other
        OmGbar= cosmo->param.Omega_g * pow(1.+z_next,4.);
		Omnu_masslessbar= cosmo->param.Omega_nu_rel * pow(1.+z_next,4.);
        OmRbar = OmGbar + Omnu_masslessbar;
        OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_next,3.);
		T_gamma = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_g_label);
		T_nu = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_ur_label);
        OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
		Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
		OmR2 = OmG + Omnu_massless;
		//H
        H2 = H0 * sqrt(OmL + OmM + OmRbar);
		dE2 = (H2-H)/H2/zstep_lin; //dE/dz*1/E (z_next). We update the stored value, although it's unnecessary.


		if(debug_mode > 1){ //check that derivatives are well defined, it should be for any reasonable value of npoints, but just in case:
			printf("H(z1)-H(z2)=%.3le , H(z2)=%.3le , dH/dz=%.3le \n", H-H2, H2, dE2*H2);
	    }


        //these are too annoying to keep even with debug_mode, activate manually if you want:
        //		printf("z=%le, H=%le , H_LCDM=%le , dE=%le , dE_LCDM=%le\n", z, H, H0*E_LCDM(cosmo, z), dE, dlogE_dz_LCDM(cosmo, z));
        //		printf("z=%.1le, OmL=%.1le , OmR=%.1le , R=%.1le \n", z, OmL, OmR, R1);

        //now we evolve R and R'.
		R1tilde = R1 + zstep_lin * Rp1;
		Rp1tilde = Rp1 + zstep_lin * Rpp1;

        //we need R''(z_next,R_tilde):
		Rpp2 = - Rp1tilde*(1./(1+z_next) + dE2)
		- mass_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- R1tilde/2.*(2*OmR - 2*OmL)*H0*H0/pow(H2*(1.+z_next),2.); //R''(R2,t2)

        //and these will be the next-z solutions.
		R2 = R1 + zstep_lin/2.0 * (Rp1 + Rp1tilde);
		Rp2 = Rp1 + zstep_lin/2.0 * (Rpp1 + Rpp2);

        //printf("R1=%.2e,  Rp1=%.2e  Rpp1=%.2e,  R2=%.2e,  Rp2=%.2e\n",R1,Rp1,Rpp1,R2,Rp2);
	}
    return z;

}
// double find_z_collapse_nothing
// (clucu_cosmology * const cosmo, const double Ri, const double Rpi, const double delta_long,const double k_long,const double mass)
// {
//     int precision = (int) fmax(fmin(precision_normalization,10),1);//1<=precision<=20


// 	const long npoints = precision*30000; //number of points for the ODE solution to converge, ~10^4 logspaced in z.



// 	double zf_code=z_collapse/2.0;
// 	if(zf_code<z_collapse_min){
// 		printf("Warning, for z_collapse = 0 you need to change z binning to linear in collapse code \n");
//     zf_code=z_collapse_min;
// 	}

// 	const double zstep_log=(log(zf_code/z_initial))/(npoints-1.); //log-spaced in z.
// 	double zstep_lin; //for integration of ODE, linear step.


// 	double R1, R2; //R and R_next.
// 	double Rp1, Rp2; //derivative of R and next value.
// 	double R1tilde, Rp1tilde; //the "guesses" of R and Rprime for  Heun's methood

// 	R1=R2=Ri; //initial conditions
// 	Rp1=Rp2=Rpi;


// 	long i;
// 	double z=z_initial;
// 	double dE, OmL, OmR, H, OmM; //z-dependent.
// 	double H2; //H for previous redshift, for derivatives.

// 	double OmG, Omnu_massless; //photon, massless and massive neutrino energy density.
// 	double OmRbar, OmGbar, Omnu_masslessbar; //average for all of them



// 	double T_gamma, T_nu;//transfer functions at each k and z.

// 	const double T_matter = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_cb_label);
// //TODO,check 
// int status=0;
// double H0_Mpc = clucu_hubble_parameter(cosmo,0.,&status);
// double Mhalo_Mpc = mass * MsuntoMpc;
// //we set the initial H
// 	OmGbar= cosmo->param.Omega_g * pow(1.+z_initial,4.);
// 	Omnu_masslessbar= cosmo->param.Omega_nu_rel * pow(1.+z_initial,4.);
// 	OmRbar = OmGbar + Omnu_masslessbar;
// 	OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_initial,3.);
// 	OmL = cosmo->param.Omega_l; //Omega_L(z), we take it as z-independent.
// 	H = H0_Mpc* sqrt(OmL + OmM + OmRbar); //H(z_initial)


// //and the initial OmR
// 	T_gamma = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_g_label);
// 	T_nu = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_ur_label);
// 	OmG= OmGbar * (1.+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
// 	Omnu_massless= Omnu_masslessbar * (1.+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
// 	double OmR2 = OmG + Omnu_massless;

// // printf("H=%.6e,  H0=%.6e,  M=%.6e,  g=%.6e,   ur=%.6e,  cb=%.6e,  l=%.6e\n",H,H0_Mpc,Mhalo_Mpc,cosmo->param.Omega_g,cosmo->param.Omega_nu_rel ,(cosmo->param.Omega_c+cosmo->param.Omega_b),cosmo->param.Omega_l);
// // printf("Tm=%.6e,   Tg=%.6e,  Tnu=%.6e\n",T_matter,T_gamma,T_nu);

// 	double z_next, dE2; //for the next values, to do Heun's method. We assume Omega_L does not change.

// //symbol 2 denotes next step.
// //we also get initial dE, for which we need H2 (next z)
// 	z_next = z_initial * exp(zstep_log);//first step
// 	OmGbar= cosmo->param.Omega_g * pow(1.+z_next,4.);
// 	Omnu_masslessbar= cosmo->param.Omega_nu_rel * pow(1.+z_next,4.);
// 	OmRbar = OmGbar + Omnu_masslessbar;
// 	OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_next,3.);
// 	OmL = cosmo->param.Omega_l;
// 	H2  = H0_Mpc * sqrt(OmL + OmM + OmRbar);// H(z_next)
// 	dE2 = (H2-H)/H2/(z_initial*zstep_log);

// 	double Rpp1, Rpp2; //d^2R(z)/dz^2
    
// //printf("H=%.6e,  H2=%.6e,  dE2=%.2e\n",H,H2,dE2);

// ///////////////////////////////////////////////////////////
// ////    here we solve for the collapse								////
// /////////////////////////////////////////////////////////
// 	for(int i=1; i<npoints-1 && R2>0.; i++){
// 		//we update the values from the previous step
// 			R1 = R2;
// 			Rp1 = Rp2;
// 		//including all z-dependent quantities, for Heun's method. Remember 2 means next step.
// 			z=z_next;
// 			H=H2;
// 			OmR=OmR2;
// 			dE=dE2;
        
//         // if(i==1)
//         //     printf("[%d]R1=%.6e,  Rp1=%.6e,  z=%.6e,  H=%.6e,  OmR=%.6e,  dE=%.6e\n",i,R1,Rp1,z,H,OmR,dE);

// 		//current redshift is updated at the end
// 		z_next *= exp(zstep_log); //next redshift
// 		zstep_lin = z_next - z; //linear step, more precise than z*zstep_log and cheap.

// 	//we find R''(z). Only depends on z and other variables calculated at the previous z.
// 		Rpp1 = - Rp1*(1./(1+z) + dE)
// 		- Mhalo_Mpc/R1/R1/pow(H*(1.+z),2.)
// 		- R1/2.*(2*OmR - 2*OmL)*H0_Mpc*H0_Mpc/pow(H*(1.+z),2.); //R''(R1,t1)
// // if((i-1)%1000==0 && i>70000)
// //     printf("[%d]z=%.6e,  R1=%.6e,  Rp1=%.6e,  Rpp1=%.6e,  dE=%.6e  OmR=%.6e,  OmL=%.6e,  H=%.6e\n",i,z,R1,Rp1,Rpp1,dE,OmR,OmL,H);


// //we update all z-dependent stuff to the next z
// 		OmGbar= cosmo->param.Omega_g * pow(1.+z_next,4.);
// 		Omnu_masslessbar= cosmo->param.Omega_nu_rel * pow(1.+z_next,4.);
// 		OmRbar = OmGbar + Omnu_masslessbar;
// 		OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_next,3.);


// 		T_gamma = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_g_label);
// 		T_nu = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_ur_label);

// 		OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
// 		Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
// 		OmR2 = OmG + Omnu_massless;


// 		H2  = H0_Mpc * sqrt(OmL + OmM + OmRbar);// H(z_next)

// 		dE2 = (H2-H)/H2/zstep_lin; //dE/dz*1/E (z_next). We update the stored value at the end.



// //these are too annoying to keep even with debug_mode, activate manually if you want to debug:
// //		printf("z=%le, H=%le , H_LCDM=%le , dE=%le , dE_LCDM=%le\n", z, H, H0_Mpc*E_LCDM(cosmo, z), dE, dlogE_dz_LCDM(cosmo, z));
// //		printf("z=%.1le, OmL=%.1le , OmR=%.1le , R=%.1le \n", z, OmL, OmR, R1);


// //now we evolve R and R'.

// 		R1tilde = R1 + zstep_lin * Rp1;
// 		Rp1tilde = Rp1 + zstep_lin * Rpp1;


// //we need R''(z_next,R_tilde):
// 		Rpp2 = - Rp1tilde*(1./(1+z_next) + dE2)
// 		- Mhalo_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
// 		- R1tilde/2.*(2*OmR2 - 2*OmL)*H0_Mpc*H0_Mpc/pow(H2*(1.+z_next),2.); //R''(R2,t2)

// //and these will be the next-z solutions.
// 		R2 = R1 + zstep_lin/2.0 * (Rp1 + Rp1tilde);
// 		Rp2 = Rp1 + zstep_lin/2.0 * (Rpp1 + Rpp2);


// //printf("R1=%.2e,  Rp1=%.2e  Rpp1=%.2e,  R2=%.2e,  Rp2=%.2e\n",R1,Rp1,Rpp1,R2,Rp2);


// 	}

// 	return z;
// }

double find_z_collapse_1nu
(clucu_cosmology * const cosmo, const double Ri, const double Rpi, const double delta_long,const double k_long,const double mass,
const long Nz_solution, double *R_solution, double *Mnu_solution){
// This function returns the z at which an overdensity collapses. global cosmological parameters assumed.
// In the presence of a long-wavelength PHOTON, MASSLESS and 1 MASSIVE NU (OR OTHER EXTRA) perturbation.
//we read the instantaneous transfer function from CAMB at each redshift.
//and we save R(z) to perform BKT calculation of neutrino clustering, or we input Mnu(z) if we already calculated it.
//we only modify R_solution and Mnu_solution, the rest are read-only.
for(int i=0;i<Nz_solution;i++){
    if(Mnu_solution[i]>0)
        printf("[[%d]]M=%.2e\n",i,Mnu_solution[i]);
}
    int status=0;
    double H0 = clucu_hubble_parameter(cosmo,0.,&status);
    double h = cosmo->param.h;
    double T0_nu = clucu_constants.TNCDM * clucu_constants.T_CMB;
    double mass_Mpc = mass * MsuntoMpc;
    const double T_matter = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_cb_label);

	const int precision = (int) fmax(fmin(precision_normalization,10),1);//1<=precision<=20
	const long npoints = precision*30000; //number of points for the ODE solution to converge, ~10^4 logspaced in z.

	double zf_code=z_collapse/2.0;
	if(zf_code<z_collapse_min){
		printf("Warning, for z_collapse = 0 you need to change z binning to linear in collapse code \n");
    zf_code=z_collapse_min;
	}
	const double zstep_log=(log(zf_code/z_initial))/(npoints-1.); //log-spaced in z.
	double zstep_lin; //for integration of ODE, linear step.


	double R1, R2; //R and R_next.
	double Rp1, Rp2; //derivative of R and next value.
    double Rpp1, Rpp2; //d^2R(z)/dz^2
	double R1tilde, Rp1tilde; //the "guesses" of R and Rprime for  Heun's methood
	R1=R2=Ri; //initial conditions
	Rp1=Rp2=Rpi;

	//the nu_massive can be any new species, just make sure that the Omegas are well defined here
    // if((cosmo->Omega_extra>0) && (cosmo->mnu1>0)){
    // printf("Using wrong collapse function, XXX_1nu only accepts one species. \n");
    // }

//set the initial
//nu1
	double rho_nu1_z0=density_WDM(cosmo->param.m_nu[0],T0_nu * (1.+0.));
    double rho_nu1_z=density_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_initial));
    double rho_nu1_ratio_z=rho_nu1_ratio_z=rho_nu1_z/rho_nu1_z0;
    double Omnu1bar_z=cosmo->param.Omega_nu_mass[0]*rho_nu1_ratio_z;
    double p_nu1_z = pressure_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_initial)); //pressure of nu1 at z
	double w_nu1_z = p_nu1_z/rho_nu1_z; //equation of state
	double T_nu1_z = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_nu1_label);
	double delta_nu1_z = delta_long*T_nu1_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next
//other,(we assme c_s^2=w=1/3 for photons and massless nu)
	double OmGbar= cosmo->param.Omega_g * pow(1.+z_initial,4.); //photon
	double Omnu_masslessbar= cosmo->param.Omega_nu_rel * pow(1.+z_initial,4.); //massless nu
	double OmRbar = OmGbar + Omnu_masslessbar; //all radiation (massless nu + photon)
	double OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_initial,3.); //matter
	double OmL = cosmo->param.Omega_l; //Omega_Lambda(z), we take it as z-independent.
    double T_gamma = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_g_label);
	double T_nu = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_ur_label);
	double OmG = OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
	double Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
	double OmR, OmR2 = OmG + Omnu_massless;
//H
    double H = H0 * sqrt(OmL + OmM + OmRbar + Omnu1bar_z); //H(z_initial), 1/Mpc,TODO，这里求H可以直接用函数代替
	


//for the next values, to do Heun's method. We assume Omega_L does not change. Symbol 2 denotes next step.
//we also get initial dE, for which we need H2 (next z)
	double z_next = z_initial * exp(zstep_log);//first step
//and the w and c_s^2 of nus at first step we calculate.
//for nu1 
    rho_nu1_z = density_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_next));
    rho_nu1_ratio_z=rho_nu1_z/rho_nu1_z0;
    Omnu1bar_z=cosmo->param.Omega_nu_mass[0]*rho_nu1_ratio_z;
    p_nu1_z = pressure_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_next));
	double w_nu1_z2=p_nu1_z/rho_nu1_z;
	double dw_nu1_z = (w_nu1_z2-w_nu1_z)/(z_next-z_initial);
	double csq_ad_nu1_z = w_nu1_z2 + dw_nu1_z/(3.0*(1+w_nu1_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
//other
    OmGbar = cosmo->param.Omega_g * pow(1.+z_next,4.);
	Omnu_masslessbar = cosmo->param.Omega_nu_rel * pow(1.+z_next,4.);
	OmRbar = OmGbar + Omnu_masslessbar;
	OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_next,3.);
	OmL = cosmo->param.Omega_l;
//H2
    double H2 = H0 * sqrt(OmL + OmM + OmRbar + Omnu1bar_z); //H(next), 1/Mpc
	double dE, dE2 = (H2-H)/H2/(z_initial*zstep_log);


	int i_solution=0; //we only check for the solution of neutrino clustering every few steps
	double mass_nu1_Mpc=0; //M of the nu1 halo
	int collapse_steps= (int) npoints/Nz_solution;	//how many steps we take before saving R(z) and/or updating Mnu
	//we save manually the first value of R and mass to avoid running the loop at i=0.
    R_solution[i_solution]=R1;
    mass_nu1_Mpc = Mnu_solution[i_solution];
    // printf("collapse_steps= (int) npoints/Nz_solution = %d/%d=%d\n",npoints,Nz_solution,collapse_steps);
    // printf("[[%d]]Rs=%.2e,Ms=%.2e\n",i_solution,R_solution[i_solution],Mnu_solution[i_solution]);
    i_solution++;
    

///////////////////////////////////////////////////////////
////    here we solve for the collapse				  ////
/////////////////////////////////////////////////////////
    double z;
    long i;
	for(i=1; i<npoints-1 && R2>0.; i++){
        //we update the values from the previous step
        R1 = R2;
        Rp1 = Rp2;
        //including all z-dependent quantities, for Heun's method. Remember 2 means next step.
        z=z_next;
        H=H2;
        dE=dE2;
        OmR=OmR2;
        w_nu1_z=w_nu1_z2;

        //we save R and read Mnu, only every few steps since it doesn't vary much.
        if (i % collapse_steps == 0){
            R_solution[i_solution]=R1;
            mass_nu1_Mpc = Mnu_solution[i_solution];
            // printf("i%%collapse_steps=%d%%%d\n",i,collapse_steps);
            // printf("[[%d]]Rs=%.2e,Ms=%.2e\n",i_solution,R_solution[i_solution],Mnu_solution[i_solution]);
            i_solution++;
        if(debug_mode > 1){
            printf("z=%.1le, R=%.1le, Mnu/mass=%1le \n\n", z, R_solution[i_solution], mass_nu1_Mpc/mass);
            }
        }


		//current redshift is updated at the end
		z_next *= exp(zstep_log); //next redshift
		zstep_lin = z_next - z; //linear step, more precise than z*zstep_log and cheap.

        //these are the terms that go inside the sum, define them outside for clarity.
		double Onu1 = (1.0+3.0*w_nu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z;


	    //we find R''(z). Only depends on z and other variables calculated at the previous z.
		Rpp1 = - Rp1*(1./(1+z) + dE)
		- mass_Mpc/R1/R1/pow(H*(1.+z),2.)
		- mass_nu1_Mpc/R1/R1/pow(H*(1.+z),2.)
		- R1/2.*(2*OmR + Onu1 - 2*OmL)*H0*H0/pow(H*(1.+z),2.); //R''(R1,t1)
//if((i-1)%1000==0)
    //printf("[%d]Rpp1=%.2e,  Rp1=%.2e,  mass_nu1_Mpc=%.2e\n",i,Rpp1,Rp1,mass_nu1_Mpc);
//if(Rp1<0)
    //printf("[%d]Rpp1=%.2e,  Rp1=%.2e,  dE=%.2e  z=%.2e,  R1=%.2e,  mass_Mpc=%.2e,  OmR=%.2e,  Onu1=%.2e,  OmL=%.2e,  H0=%.2e  H=%.2e\n",i,Rpp1,Rp1,dE,z,R1,mass_Mpc,OmR,Onu1,OmL,H0,H);

        //we update all z-dependent stuff to the next z
        //nu1
        rho_nu1_z = density_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_next));
        rho_nu1_ratio_z=rho_nu1_z/rho_nu1_z0;
        Omnu1bar_z=cosmo->param.Omega_nu_mass[0]*rho_nu1_ratio_z;
        T_nu1_z = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_nu1_label);
		delta_nu1_z = delta_long*T_nu1_z/T_matter; //long-wavlength neutrino perturbation
        p_nu1_z = pressure_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_next));
	    w_nu1_z2=p_nu1_z/rho_nu1_z;
	    dw_nu1_z = (w_nu1_z2-w_nu1_z)/zstep_lin;
	    csq_ad_nu1_z = w_nu1_z2 + dw_nu1_z/(3.0*(1+w_nu1_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
        Onu1 = (1.0+3.0*w_nu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z;
        //other
        OmGbar= cosmo->param.Omega_g * pow(1.+z_next,4.);
		Omnu_masslessbar= cosmo->param.Omega_nu_rel * pow(1.+z_next,4.);
        OmRbar = OmGbar + Omnu_masslessbar;
        OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_next,3.);
		T_gamma = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_g_label);
		T_nu = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_ur_label);
        OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
		Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
		OmR2 = OmG + Omnu_massless;
		//H
        H2 = H0 * sqrt(OmL + OmM + OmRbar + Omnu1bar_z);
		dE2 = (H2-H)/H2/zstep_lin; //dE/dz*1/E (z_next). We update the stored value, although it's unnecessary.


		if(debug_mode > 1){ //check that derivatives are well defined, it should be for any reasonable value of npoints, but just in case:
			printf("w(z1)-w(z2)=%.3le , w(z2)=%.3le , dw/dz=%.3le \n", w_nu1_z-w_nu1_z2, w_nu1_z2, dw_nu1_z);
			printf("H(z1)-H(z2)=%.3le , H(z2)=%.3le , dH/dz=%.3le \n", H-H2, H2, dE2*H2);
			printf("z=%.3le, cs^2/w(z2)=%.3le \n", z_next, csq_ad_nu1_z/w_nu1_z2); //cs^2 and w are within 10% of each other.
	    }


        //these are too annoying to keep even with debug_mode, activate manually if you want:
        //		printf("z=%le, H=%le , H_LCDM=%le , dE=%le , dE_LCDM=%le\n", z, H, H0*E_LCDM(cosmo, z), dE, dlogE_dz_LCDM(cosmo, z));
        //		printf("z=%.1le, OmL=%.1le , OmR=%.1le , R=%.1le \n", z, OmL, OmR, R1);

        //now we evolve R and R'.
		R1tilde = R1 + zstep_lin * Rp1;
		Rp1tilde = Rp1 + zstep_lin * Rpp1;

        //we need R''(z_next,R_tilde):
		Rpp2 = - Rp1tilde*(1./(1+z_next) + dE2)
		- mass_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- mass_nu1_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- R1tilde/2.*(2*OmR + Onu1 - 2*OmL)*H0*H0/pow(H2*(1.+z_next),2.); //R''(R2,t2)

        //and these will be the next-z solutions.
		R2 = R1 + zstep_lin/2.0 * (Rp1 + Rp1tilde);
		Rp2 = Rp1 + zstep_lin/2.0 * (Rpp1 + Rpp2);
        // if((i-1)%1000==0)
        //     printf("[%d]z=%.2e,  R2=%.2e  R1=%.2e,  Rp1=%.2e,  Rp1tilde=%.2e,  Rpp1=%.2e\n",i,z,R2,R1,Rp1,Rp1tilde,Rpp1);
        // if(i==1208 || i==1209)
        //     printf("[%d]z=%.2e,  R2=%.2e  R1=%.2e,  Rp1=%.2e,  Rp1tilde=%.2e,  Rpp1=%.2e\n",i,z,R2,R1,Rp1,Rp1tilde,Rpp1);
	}
    //printf("[%d]Rp1=%.2e,  dE=%.2e  z=%.2e,  R1=%.2e,  mass_Mpc=%.2e,  OmR=%.2e,  OmL=%.2e,  H0=%.2e  H=%.2e\n",i,Rp1,dE,z,R1,mass_Mpc,OmR,OmL,H0,H);
    return z;

}
double find_z_collapse_2nu
(clucu_cosmology * const cosmo, const double Ri, const double Rpi, const double delta_long,const double k_long,const double mass,
const long Nz_solution, double *R_solution, double *Mnu1_solution, double *Mnu2_solution){
// This function returns the z at which an overdensity collapses. global cosmological parameters assumed.
// In the presence of a long-wavelength PHOTON, MASSLESS and 1 MASSIVE NU (OR OTHER EXTRA) perturbation.
//we read the instantaneous transfer function from CAMB at each redshift.
//and we save R(z) to perform BKT calculation of neutrino clustering, or we input Mnu(z) if we already calculated it.
//we only modify R_solution and Mnu_solution, the rest are read-only.

    int status=0;
    double H0 = clucu_hubble_parameter(cosmo,0.,&status);
    double h = cosmo->param.h;
    double T0_nu = clucu_constants.TNCDM * clucu_constants.T_CMB;
    double mass_Mpc = mass * MsuntoMpc;
    const double T_matter = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_cb_label);

	const int precision = (int) fmax(fmin(precision_normalization,10),1);//1<=precision<=20
	const long npoints = precision*30000; //number of points for the ODE solution to converge, ~10^4 logspaced in z.

	double zf_code=z_collapse/2.0;
	if(zf_code<z_collapse_min){
		printf("Warning, for z_collapse = 0 you need to change z binning to linear in collapse code \n");
    zf_code=z_collapse_min;
	}
	const double zstep_log=(log(zf_code/z_initial))/(npoints-1.); //log-spaced in z.
	double zstep_lin; //for integration of ODE, linear step.


	double R1, R2; //R and R_next.
	double Rp1, Rp2; //derivative of R and next value.
    double Rpp1, Rpp2; //d^2R(z)/dz^2
	double R1tilde, Rp1tilde; //the "guesses" of R and Rprime for  Heun's methood
	R1=R2=Ri; //initial conditions
	Rp1=Rp2=Rpi;

	//the nu_massive can be any new species, just make sure that the Omegas are well defined here
    // if((cosmo->Omega_extra>0) && (cosmo->mnu1>0)){
    // printf("Using wrong collapse function, XXX_1nu only accepts one species. \n");
    // }

//set the initial
//nu1
	double rho_nu1_z0=density_WDM(cosmo->param.m_nu[0],T0_nu * (1.+0.));
    double rho_nu1_z=density_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_initial));
    double rho_nu1_ratio_z=rho_nu1_z/rho_nu1_z0;
    double Omnu1bar_z=cosmo->param.Omega_nu_mass[0]*rho_nu1_ratio_z;
    double p_nu1_z = pressure_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_initial)); //pressure of nu1 at z
	double w_nu1_z = p_nu1_z/rho_nu1_z; //equation of state
	double T_nu1_z = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_nu1_label);
	double delta_nu1_z = delta_long*T_nu1_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next
//nu2
	double rho_nu2_z0=density_WDM(cosmo->param.m_nu[1],T0_nu * (1.+0.));
    double rho_nu2_z=density_WDM(cosmo->param.m_nu[1],T0_nu * (1.+z_initial));
    double rho_nu2_ratio_z=rho_nu2_z/rho_nu2_z0;
    double Omnu2bar_z=cosmo->param.Omega_nu_mass[1]*rho_nu2_ratio_z;
    double p_nu2_z = pressure_WDM(cosmo->param.m_nu[1],T0_nu * (1.+z_initial)); //pressure of nu1 at z
	double w_nu2_z = p_nu2_z/rho_nu2_z; //equation of state
	double T_nu2_z = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_nu2_label);
	double delta_nu2_z = delta_long*T_nu2_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next
//other,(we assme c_s^2=w=1/3 for photons and massless nu)
	double OmGbar= cosmo->param.Omega_g * pow(1.+z_initial,4.); //photon
	double Omnu_masslessbar= cosmo->param.Omega_nu_rel * pow(1.+z_initial,4.); //massless nu
	double OmRbar = OmGbar + Omnu_masslessbar; //all radiation (massless nu + photon)
	double OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_initial,3.); //matter
	double OmL = cosmo->param.Omega_l; //Omega_Lambda(z), we take it as z-independent.
    double T_gamma = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_g_label);
	double T_nu = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_ur_label);
	double OmG = OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
	double Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
	double OmR, OmR2 = OmG + Omnu_massless;
//H
    double H = H0 * sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z); //H(z_initial), 1/Mpc,TODO，这里求H可以直接用函数代替
	


//for the next values, to do Heun's method. We assume Omega_L does not change. Symbol 2 denotes next step.
//we also get initial dE, for which we need H2 (next z)
	double z_next = z_initial * exp(zstep_log);//first step
//and the w and c_s^2 of nus at first step we calculate.
//for nu1 
    rho_nu1_z = density_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_next));
    rho_nu1_ratio_z=rho_nu1_z/rho_nu1_z0;
    Omnu1bar_z=cosmo->param.Omega_nu_mass[0]*rho_nu1_ratio_z;
    p_nu1_z = pressure_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_next));
	double w_nu1_z2=p_nu1_z/rho_nu1_z;
	double dw_nu1_z = (w_nu1_z2-w_nu1_z)/(z_next-z_initial);
	double csq_ad_nu1_z = w_nu1_z2 + dw_nu1_z/(3.0*(1+w_nu1_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
//for nu2
    rho_nu2_z = density_WDM(cosmo->param.m_nu[1],T0_nu * (1.+z_next));
    rho_nu2_ratio_z=rho_nu2_z/rho_nu2_z0;
    Omnu2bar_z=cosmo->param.Omega_nu_mass[1]*rho_nu2_ratio_z;
	p_nu2_z = pressure_WDM(cosmo->param.m_nu[1],T0_nu * (1.+z_next));
	double w_nu2_z2=p_nu2_z/rho_nu2_z;
	double dw_nu2_z = (w_nu2_z2-w_nu2_z)/(z_next-z_initial);
	double csq_ad_nu2_z = w_nu2_z2 + dw_nu2_z/(3.0*(1+w_nu2_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
//other
    OmGbar = cosmo->param.Omega_g * pow(1.+z_next,4.);
	Omnu_masslessbar = cosmo->param.Omega_nu_rel * pow(1.+z_next,4.);
	OmRbar = OmGbar + Omnu_masslessbar;
	OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_next,3.);
	OmL = cosmo->param.Omega_l;
//H2
    double H2 = H0 * sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z); //H(next), 1/Mpc
	double dE, dE2 = (H2-H)/H2/(z_initial*zstep_log);


	int i_solution=0; //we only check for the solution of neutrino clustering every few steps
	double mass_nu1_Mpc=0; //M of the nu1 halo
    double mass_nu2_Mpc=0; //M of the nu2 halo
	int collapse_steps= (int) npoints/Nz_solution;	//how many steps we take before saving R(z) and/or updating Mnu
	//we save manually the first value of R and mass to avoid running the loop at i=0.
    R_solution[i_solution]=R1;
    mass_nu1_Mpc = Mnu1_solution[i_solution];
    mass_nu2_Mpc = Mnu2_solution[i_solution];
    i_solution++;

///////////////////////////////////////////////////////////
////    here we solve for the collapse				  ////
/////////////////////////////////////////////////////////
    double z;
	for(long i=1; i<npoints-1 && R2>0.; i++){
        //we update the values from the previous step
        R1 = R2;
        Rp1 = Rp2;
        //including all z-dependent quantities, for Heun's method. Remember 2 means next step.
        z=z_next;
        H=H2;
        dE=dE2;
        OmR=OmR2;
        w_nu1_z=w_nu1_z2;
        w_nu2_z=w_nu2_z2;

        //we save R and read Mnu, only every few steps since it doesn't vary much.
        if (i % collapse_steps == 0){
            R_solution[i_solution]=R1;
            mass_nu1_Mpc = Mnu1_solution[i_solution];
            mass_nu2_Mpc = Mnu2_solution[i_solution];
            i_solution++;
        if(debug_mode > 1){
            printf("z=%.1le, R=%.1le, Mnu/mass=%1le \n\n", z, R_solution[i_solution], mass_nu1_Mpc/mass);
            }
        }


		//current redshift is updated at the end
		z_next *= exp(zstep_log); //next redshift
		zstep_lin = z_next - z; //linear step, more precise than z*zstep_log and cheap.

        //these are the terms that go inside the sum, define them outside for clarity.
		double Onu1 = (1.0+3.0*w_nu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z;
		double Onu2 = (1.0+3.0*w_nu2_z + delta_nu2_z*(1.0+3.0*csq_ad_nu2_z))*Omnu2bar_z;

	    //we find R''(z). Only depends on z and other variables calculated at the previous z.
		Rpp1 = - Rp1*(1./(1+z) + dE)
		- mass_Mpc/R1/R1/pow(H*(1.+z),2.)
		- mass_nu1_Mpc/R1/R1/pow(H*(1.+z),2.)
        - mass_nu2_Mpc/R1/R1/pow(H*(1.+z),2.)
		- R1/2.*(2*OmR + Onu1 + Onu2 - 2*OmL)*H0*H0/pow(H*(1.+z),2.); //R''(R1,t1)
// if((i-1)%1000==0)
//     printf("[%d]Rpp1=%.2e,  Rp1=%.2e,  dE=%.2e  z=%.2e,  R1=%.2e,  mass_Mpc=%.2e,  OmR=%.2e,  Onu1=%.2e,  OmL=%.2e,  H0=%.2e  H=%.2e\n",i,Rpp1,Rp1,dE,z,R1,mass_Mpc,OmR,Onu1,OmL,H0,H);

        //we update all z-dependent stuff to the next z
        //nu1
        rho_nu1_z = density_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_next));
        rho_nu1_ratio_z=rho_nu1_z/rho_nu1_z0;
        Omnu1bar_z=cosmo->param.Omega_nu_mass[0]*rho_nu1_ratio_z;
        T_nu1_z = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_nu1_label);
		delta_nu1_z = delta_long*T_nu1_z/T_matter; //long-wavlength neutrino perturbation
        p_nu1_z = pressure_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_next));
	    w_nu1_z2=p_nu1_z/rho_nu1_z;
	    dw_nu1_z = (w_nu1_z2-w_nu1_z)/zstep_lin;
	    csq_ad_nu1_z = w_nu1_z2 + dw_nu1_z/(3.0*(1+w_nu1_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
        Onu1 = (1.0+3.0*w_nu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z;
        //nu2
        rho_nu2_z = density_WDM(cosmo->param.m_nu[1],T0_nu * (1.+z_next));
        rho_nu2_ratio_z=rho_nu2_z/rho_nu2_z0;
        Omnu2bar_z=cosmo->param.Omega_nu_mass[1]*rho_nu2_ratio_z;
        T_nu2_z = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_nu2_label);
        delta_nu2_z = delta_long*T_nu2_z/T_matter; //long-wavlength neutrino perturbation
        p_nu2_z = pressure_WDM(cosmo->param.m_nu[1],T0_nu * (1.+z_next));
        w_nu2_z2=p_nu2_z/rho_nu2_z;
        dw_nu2_z = (w_nu2_z2-w_nu2_z)/zstep_lin;
	    csq_ad_nu2_z = w_nu2_z2 + dw_nu2_z/(3.0*(1+w_nu2_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
        Onu2 = (1.0+3.0*w_nu2_z + delta_nu2_z*(1.0+3.0*csq_ad_nu2_z))*Omnu2bar_z;
        //other
        OmGbar= cosmo->param.Omega_g * pow(1.+z_next,4.);
		Omnu_masslessbar= cosmo->param.Omega_nu_rel * pow(1.+z_next,4.);
        OmRbar = OmGbar + Omnu_masslessbar;
        OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_next,3.);
		T_gamma = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_g_label);
		T_nu = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_ur_label);
        OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
		Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
		OmR2 = OmG + Omnu_massless;
		//H
        H2 = H0 * sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z);
		dE2 = (H2-H)/H2/zstep_lin; //dE/dz*1/E (z_next). We update the stored value, although it's unnecessary.


		if(debug_mode > 1){ //check that derivatives are well defined, it should be for any reasonable value of npoints, but just in case:
			printf("w(z1)-w(z2)=%.3le , w(z2)=%.3le , dw/dz=%.3le \n", w_nu1_z-w_nu1_z2, w_nu1_z2, dw_nu1_z);
			printf("H(z1)-H(z2)=%.3le , H(z2)=%.3le , dH/dz=%.3le \n", H-H2, H2, dE2*H2);
			printf("z=%.3le, cs^2/w(z2)=%.3le \n", z_next, csq_ad_nu1_z/w_nu1_z2); //cs^2 and w are within 10% of each other.
	    }


        //these are too annoying to keep even with debug_mode, activate manually if you want:
        //		printf("z=%le, H=%le , H_LCDM=%le , dE=%le , dE_LCDM=%le\n", z, H, H0*E_LCDM(cosmo, z), dE, dlogE_dz_LCDM(cosmo, z));
        //		printf("z=%.1le, OmL=%.1le , OmR=%.1le , R=%.1le \n", z, OmL, OmR, R1);

        //now we evolve R and R'.
		R1tilde = R1 + zstep_lin * Rp1;
		Rp1tilde = Rp1 + zstep_lin * Rpp1;

        //we need R''(z_next,R_tilde):
		Rpp2 = - Rp1tilde*(1./(1+z_next) + dE2)
		- mass_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- mass_nu1_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
        - mass_nu2_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- R1tilde/2.*(2*OmR + Onu1 + Onu2 - 2*OmL)*H0*H0/pow(H2*(1.+z_next),2.); //R''(R2,t2)

        //and these will be the next-z solutions.
		R2 = R1 + zstep_lin/2.0 * (Rp1 + Rp1tilde);
		Rp2 = Rp1 + zstep_lin/2.0 * (Rpp1 + Rpp2);
	}
    return z;

}
double find_z_collapse_3nu
(clucu_cosmology * const cosmo, const double Ri, const double Rpi, const double delta_long,const double k_long,const double mass,
const long Nz_solution, double *R_solution, double *Mnu1_solution, double *Mnu2_solution, double *Mnu3_solution){
// This function returns the z at which an overdensity collapses. global cosmological parameters assumed.
// In the presence of a long-wavelength PHOTON, MASSLESS and 1 MASSIVE NU (OR OTHER EXTRA) perturbation.
//we read the instantaneous transfer function from CAMB at each redshift.
//and we save R(z) to perform BKT calculation of neutrino clustering, or we input Mnu(z) if we already calculated it.
//we only modify R_solution and Mnu_solution, the rest are read-only.

    int status=0;
    double H0 = clucu_hubble_parameter(cosmo,0.,&status);
    double h = cosmo->param.h;
    double T0_nu = clucu_constants.TNCDM * clucu_constants.T_CMB;
    double mass_Mpc = mass * MsuntoMpc;
    const double T_matter = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_cb_label);

	const int precision = (int) fmax(fmin(precision_normalization,10),1);//1<=precision<=20
	const long npoints = precision*30000; //number of points for the ODE solution to converge, ~10^4 logspaced in z.

	double zf_code=z_collapse/2.0;
	if(zf_code<z_collapse_min){
		printf("Warning, for z_collapse = 0 you need to change z binning to linear in collapse code \n");
    zf_code=z_collapse_min;
	}
	const double zstep_log=(log(zf_code/z_initial))/(npoints-1.); //log-spaced in z.
	double zstep_lin; //for integration of ODE, linear step.


	double R1, R2; //R and R_next.
	double Rp1, Rp2; //derivative of R and next value.
    double Rpp1, Rpp2; //d^2R(z)/dz^2
	double R1tilde, Rp1tilde; //the "guesses" of R and Rprime for  Heun's methood
	R1=R2=Ri; //initial conditions
	Rp1=Rp2=Rpi;

	//the nu_massive can be any new species, just make sure that the Omegas are well defined here
    // if((cosmo->Omega_extra>0) && (cosmo->mnu1>0)){
    // printf("Using wrong collapse function, XXX_1nu only accepts one species. \n");
    // }

//set the initial
//nu1
	double rho_nu1_z0=density_WDM(cosmo->param.m_nu[0],T0_nu * (1.+0.));
    double rho_nu1_z=density_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_initial));
    double rho_nu1_ratio_z=rho_nu1_z/rho_nu1_z0;
    double Omnu1bar_z=cosmo->param.Omega_nu_mass[0]*rho_nu1_ratio_z;
    double p_nu1_z = pressure_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_initial)); //pressure of nu1 at z
	double w_nu1_z = p_nu1_z/rho_nu1_z; //equation of state
	double T_nu1_z = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_nu1_label);
	double delta_nu1_z = delta_long*T_nu1_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next
//nu2
	double rho_nu2_z0=density_WDM(cosmo->param.m_nu[1],T0_nu * (1.+0.));
    double rho_nu2_z=density_WDM(cosmo->param.m_nu[1],T0_nu * (1.+z_initial));
    double rho_nu2_ratio_z=rho_nu2_z/rho_nu2_z0;
    double Omnu2bar_z=cosmo->param.Omega_nu_mass[1]*rho_nu2_ratio_z;
    double p_nu2_z = pressure_WDM(cosmo->param.m_nu[1],T0_nu * (1.+z_initial)); //pressure of nu1 at z
	double w_nu2_z = p_nu2_z/rho_nu2_z; //equation of state
	double T_nu2_z = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_nu2_label);
	double delta_nu2_z = delta_long*T_nu2_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next
//nu3
	double rho_nu3_z0=density_WDM(cosmo->param.m_nu[2],T0_nu * (1.+0.));
    double rho_nu3_z=density_WDM(cosmo->param.m_nu[2],T0_nu * (1.+z_initial));
    double rho_nu3_ratio_z=rho_nu3_z/rho_nu3_z0;
    double Omnu3bar_z=cosmo->param.Omega_nu_mass[2]*rho_nu3_ratio_z;
    double p_nu3_z = pressure_WDM(cosmo->param.m_nu[2],T0_nu * (1.+z_initial)); //pressure of nu1 at z
	double w_nu3_z = p_nu3_z/rho_nu3_z; //equation of state
	double T_nu3_z = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_nu3_label);
	double delta_nu3_z = delta_long*T_nu3_z/T_matter; //the long-wavelength nu1 perturbation, at z and z_next
//other,(we assme c_s^2=w=1/3 for photons and massless nu)
	double OmGbar= cosmo->param.Omega_g * pow(1.+z_initial,4.); //photon
	double Omnu_masslessbar= cosmo->param.Omega_nu_rel * pow(1.+z_initial,4.); //massless nu
	double OmRbar = OmGbar + Omnu_masslessbar; //all radiation (massless nu + photon)
	double OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_initial,3.); //matter
	double OmL = cosmo->param.Omega_l; //Omega_Lambda(z), we take it as z-independent.
    double T_gamma = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_g_label);
	double T_nu = clucu_transfer_label(cosmo,k_long,z_initial,clucu_species_ur_label);
	double OmG = OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
	double Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
	double OmR, OmR2 = OmG + Omnu_massless;
//H
    double H = H0 * sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z + Omnu3bar_z); //H(z_initial), 1/Mpc,TODO，这里求H可以直接用函数代替
	


//for the next values, to do Heun's method. We assume Omega_L does not change. Symbol 2 denotes next step.
//we also get initial dE, for which we need H2 (next z)
	double z_next = z_initial * exp(zstep_log);//first step
//and the w and c_s^2 of nus at first step we calculate.
//for nu1 
    rho_nu1_z = density_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_next));
    rho_nu1_ratio_z=rho_nu1_z/rho_nu1_z0;
    Omnu1bar_z=cosmo->param.Omega_nu_mass[0]*rho_nu1_ratio_z;
    p_nu1_z = pressure_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_next));
	double w_nu1_z2=p_nu1_z/rho_nu1_z;
	double dw_nu1_z = (w_nu1_z2-w_nu1_z)/(z_next-z_initial);
	double csq_ad_nu1_z = w_nu1_z2 + dw_nu1_z/(3.0*(1+w_nu1_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
//for nu2
    rho_nu2_z = density_WDM(cosmo->param.m_nu[1],T0_nu * (1.+z_next));
    rho_nu2_ratio_z=rho_nu2_z/rho_nu2_z0;
    Omnu2bar_z=cosmo->param.Omega_nu_mass[1]*rho_nu2_ratio_z;
	p_nu2_z = pressure_WDM(cosmo->param.m_nu[1],T0_nu * (1.+z_next));
	double w_nu2_z2=p_nu2_z/rho_nu2_z;
	double dw_nu2_z = (w_nu2_z2-w_nu2_z)/(z_next-z_initial);
	double csq_ad_nu2_z = w_nu2_z2 + dw_nu2_z/(3.0*(1+w_nu2_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
//for nu3
    rho_nu3_z = density_WDM(cosmo->param.m_nu[2],T0_nu * (1.+z_next));
    rho_nu3_ratio_z=rho_nu3_z/rho_nu3_z0;
    Omnu3bar_z=cosmo->param.Omega_nu_mass[2]*rho_nu3_ratio_z;
	p_nu3_z = pressure_WDM(cosmo->param.m_nu[2],T0_nu * (1.+z_next));
	double w_nu3_z2=p_nu3_z/rho_nu3_z;
	double dw_nu3_z = (w_nu3_z2-w_nu3_z)/(z_next-z_initial);
	double csq_ad_nu3_z = w_nu3_z2 + dw_nu3_z/(3.0*(1+w_nu3_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
//other
    OmGbar = cosmo->param.Omega_g * pow(1.+z_next,4.);
	Omnu_masslessbar = cosmo->param.Omega_nu_rel * pow(1.+z_next,4.);
	OmRbar = OmGbar + Omnu_masslessbar;
	OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_next,3.);
	OmL = cosmo->param.Omega_l;
//H2
    double H2 = H0 * sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z + Omnu3bar_z); //H(next), 1/Mpc
	double dE, dE2 = (H2-H)/H2/(z_initial*zstep_log);


	int i_solution=0; //we only check for the solution of neutrino clustering every few steps
	double mass_nu1_Mpc=0; //M of the nu1 halo
    double mass_nu2_Mpc=0; //M of the nu2 halo
    double mass_nu3_Mpc=0; //M of the nu3 halo
	int collapse_steps= (int) npoints/Nz_solution;	//how many steps we take before saving R(z) and/or updating Mnu
	//we save manually the first value of R and mass to avoid running the loop at i=0.
    R_solution[i_solution]=R1;
    mass_nu1_Mpc = Mnu1_solution[i_solution];
    mass_nu2_Mpc = Mnu2_solution[i_solution];
    mass_nu3_Mpc = Mnu3_solution[i_solution];
    i_solution++;

///////////////////////////////////////////////////////////
////    here we solve for the collapse				  ////
/////////////////////////////////////////////////////////
    double z;
	for(long i=1; i<npoints-1 && R2>0.; i++){
        //we update the values from the previous step
        R1 = R2;
        Rp1 = Rp2;
        //including all z-dependent quantities, for Heun's method. Remember 2 means next step.
        z=z_next;
        H=H2;
        dE=dE2;
        OmR=OmR2;
        w_nu1_z=w_nu1_z2;
        w_nu2_z=w_nu2_z2;
        w_nu3_z=w_nu3_z2;

        //we save R and read Mnu, only every few steps since it doesn't vary much.
        if (i % collapse_steps == 0){
            R_solution[i_solution]=R1;
            mass_nu1_Mpc = Mnu1_solution[i_solution];
            mass_nu2_Mpc = Mnu2_solution[i_solution];
            mass_nu3_Mpc = Mnu3_solution[i_solution];
            i_solution++;
        if(debug_mode > 1){
            printf("z=%.1le, R=%.1le, Mnu/mass=%1le \n\n", z, R_solution[i_solution], mass_nu1_Mpc/mass);
            }
        }


		//current redshift is updated at the end
		z_next *= exp(zstep_log); //next redshift
		zstep_lin = z_next - z; //linear step, more precise than z*zstep_log and cheap.

        //these are the terms that go inside the sum, define them outside for clarity.
		double Onu1 = (1.0+3.0*w_nu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z;
		double Onu2 = (1.0+3.0*w_nu2_z + delta_nu2_z*(1.0+3.0*csq_ad_nu2_z))*Omnu2bar_z;
        double Onu3 = (1.0+3.0*w_nu3_z + delta_nu3_z*(1.0+3.0*csq_ad_nu3_z))*Omnu3bar_z;

	    //we find R''(z). Only depends on z and other variables calculated at the previous z.
		Rpp1 = - Rp1*(1./(1+z) + dE)
		- mass_Mpc/R1/R1/pow(H*(1.+z),2.)
		- mass_nu1_Mpc/R1/R1/pow(H*(1.+z),2.)
        - mass_nu2_Mpc/R1/R1/pow(H*(1.+z),2.)
        - mass_nu3_Mpc/R1/R1/pow(H*(1.+z),2.)
		- R1/2.*(2*OmR + Onu1 + Onu2 + Onu3 - 2*OmL)*H0*H0/pow(H*(1.+z),2.); //R''(R1,t1)
// if((i-1)%1000==0)
//     printf("[%d]Rpp1=%.2e,  Rp1=%.2e,  dE=%.2e  z=%.2e,  R1=%.2e,  mass_Mpc=%.2e,  OmR=%.2e,  Onu1=%.2e,  OmL=%.2e,  H0=%.2e  H=%.2e\n",i,Rpp1,Rp1,dE,z,R1,mass_Mpc,OmR,Onu1,OmL,H0,H);

        //we update all z-dependent stuff to the next z
        //nu1
        rho_nu1_z = density_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_next));
        rho_nu1_ratio_z=rho_nu1_z/rho_nu1_z0;
        Omnu1bar_z=cosmo->param.Omega_nu_mass[0]*rho_nu1_ratio_z;
        T_nu1_z = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_nu1_label);
		delta_nu1_z = delta_long*T_nu1_z/T_matter; //long-wavlength neutrino perturbation
        p_nu1_z = pressure_WDM(cosmo->param.m_nu[0],T0_nu * (1.+z_next));
	    w_nu1_z2=p_nu1_z/rho_nu1_z;
	    dw_nu1_z = (w_nu1_z2-w_nu1_z)/zstep_lin;
	    csq_ad_nu1_z = w_nu1_z2 + dw_nu1_z/(3.0*(1+w_nu1_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
        Onu1 = (1.0+3.0*w_nu1_z + delta_nu1_z*(1.0+3.0*csq_ad_nu1_z))*Omnu1bar_z;
        //nu2
        rho_nu2_z = density_WDM(cosmo->param.m_nu[1],T0_nu * (1.+z_next));
        rho_nu2_ratio_z=rho_nu2_z/rho_nu2_z0;
        Omnu2bar_z=cosmo->param.Omega_nu_mass[1]*rho_nu2_ratio_z;
        T_nu2_z = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_nu2_label);
        delta_nu2_z = delta_long*T_nu2_z/T_matter; //long-wavlength neutrino perturbation
        p_nu2_z = pressure_WDM(cosmo->param.m_nu[1],T0_nu * (1.+z_next));
        w_nu2_z2=p_nu2_z/rho_nu2_z;
        dw_nu2_z = (w_nu2_z2-w_nu2_z)/zstep_lin;
	    csq_ad_nu2_z = w_nu2_z2 + dw_nu2_z/(3.0*(1+w_nu2_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
        Onu2 = (1.0+3.0*w_nu2_z + delta_nu2_z*(1.0+3.0*csq_ad_nu2_z))*Omnu2bar_z;
        //nu3
        rho_nu3_z = density_WDM(cosmo->param.m_nu[2],T0_nu * (1.+z_next));
        rho_nu3_ratio_z=rho_nu3_z/rho_nu3_z0;
        Omnu3bar_z=cosmo->param.Omega_nu_mass[2]*rho_nu3_ratio_z;
        T_nu3_z = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_nu3_label);
        delta_nu3_z = delta_long*T_nu3_z/T_matter; //long-wavlength neutrino perturbation
        p_nu3_z = pressure_WDM(cosmo->param.m_nu[2],T0_nu * (1.+z_next));
        w_nu3_z2=p_nu3_z/rho_nu3_z;
        dw_nu3_z = (w_nu3_z2-w_nu3_z)/zstep_lin;
	    csq_ad_nu3_z = w_nu3_z2 + dw_nu3_z/(3.0*(1+w_nu3_z2))*(1.0+z_next); //sound speed squared, calculated as w - w'/(3*(1+w)*(1+z)).
        Onu3 = (1.0+3.0*w_nu3_z + delta_nu3_z*(1.0+3.0*csq_ad_nu3_z))*Omnu3bar_z;        
        //other
        OmGbar= cosmo->param.Omega_g * pow(1.+z_next,4.);
		Omnu_masslessbar= cosmo->param.Omega_nu_rel * pow(1.+z_next,4.);
        OmRbar = OmGbar + Omnu_masslessbar;
        OmM = (cosmo->param.Omega_c+cosmo->param.Omega_b) * pow(1.+z_next,3.);
		T_gamma = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_g_label);
		T_nu = clucu_transfer_label(cosmo,k_long,z_next,clucu_species_ur_label);
        OmG= OmGbar * (1.0+delta_long*T_gamma/T_matter); //Omega_gamma(z) (1+\delta_long_gamma)
		Omnu_massless= Omnu_masslessbar * (1.0+delta_long*T_nu/T_matter); //Omega_massless neutrino(z) (1+\delta_long_nu)
		OmR2 = OmG + Omnu_massless;
		//H
        H2 = H0 * sqrt(OmL + OmM + OmRbar + Omnu1bar_z + Omnu2bar_z + Omnu3bar_z);
		dE2 = (H2-H)/H2/zstep_lin; //dE/dz*1/E (z_next). We update the stored value, although it's unnecessary.


		if(debug_mode > 1){ //check that derivatives are well defined, it should be for any reasonable value of npoints, but just in case:
			printf("w(z1)-w(z2)=%.3le , w(z2)=%.3le , dw/dz=%.3le \n", w_nu1_z-w_nu1_z2, w_nu1_z2, dw_nu1_z);
			printf("H(z1)-H(z2)=%.3le , H(z2)=%.3le , dH/dz=%.3le \n", H-H2, H2, dE2*H2);
			printf("z=%.3le, cs^2/w(z2)=%.3le \n", z_next, csq_ad_nu1_z/w_nu1_z2); //cs^2 and w are within 10% of each other.
	    }


        //these are too annoying to keep even with debug_mode, activate manually if you want:
        //		printf("z=%le, H=%le , H_LCDM=%le , dE=%le , dE_LCDM=%le\n", z, H, H0*E_LCDM(cosmo, z), dE, dlogE_dz_LCDM(cosmo, z));
        //		printf("z=%.1le, OmL=%.1le , OmR=%.1le , R=%.1le \n", z, OmL, OmR, R1);

        //now we evolve R and R'.
		R1tilde = R1 + zstep_lin * Rp1;
		Rp1tilde = Rp1 + zstep_lin * Rpp1;

        //we need R''(z_next,R_tilde):
		Rpp2 = - Rp1tilde*(1./(1+z_next) + dE2)
		- mass_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- mass_nu1_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
        - mass_nu2_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
        - mass_nu3_Mpc/R1tilde/R1tilde/pow(H2*(1.+z_next),2.)
		- R1tilde/2.*(2*OmR + Onu1 + Onu2 + Onu3 - 2*OmL)*H0*H0/pow(H2*(1.+z_next),2.); //R''(R2,t2)

        //and these will be the next-z solutions.
		R2 = R1 + zstep_lin/2.0 * (Rp1 + Rp1tilde);
		Rp2 = Rp1 + zstep_lin/2.0 * (Rpp1 + Rpp2);
	}
    return z;

}