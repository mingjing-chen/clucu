/**
 * @file main_oguri.c
 * @author Mingjing Chen (mingjing@mail.ustc.edu.cn)
 * @brief 计算oguri2011的table2
 * @version 0.1
 * @date 2022-09-09
 * 
 * @copyright Copyright (c) 2022
 * 
 */
#include "clucu.h"

#include <time.h>
#include <string.h>
//gcc -std=c99  main_oguri.c -fopenmp /home/mjchen/Program/cluster/newcode/src/clucu_*.c -I/home/mjchen/Program/cluster/newcode/include -I/home/mjchen/Program/gsl-2.7/gsl -L/home/mjchen/Program/gsl-2.7/.libs -lgsl -lgslcblas -lm -o test

#define CASEEND cosmo[id1].class_id=0;cosmo[id2].class_id=0;cosmo[id1].runned_class=true;cosmo[id2].runned_class=true;break
#define CASEEND1(x) cosmo[id1].survey.x *= 1-delta;cosmo[id2].survey.x *= 1+delta;partial=2.0*fiducial.survey.x*delta;CASEEND
#define CASEEND0(x,dx) cosmo[id1].survey.x -= dx;cosmo[id2].survey.x += dx;partial=2.0*dx;CASEEND
//#define CASEEND0(x,dx) cosmo[id1].survey.x -= 0.1;cosmo[id2].survey.x += 0.1;partial=2.0*0.1;CASEEND
int main()
{
    int status=0;
    omp_set_nested(1);
    //用于设置是否允许OpenMP进行嵌套并行，默认的设置为false。
	time_t start,end;
    start =time(NULL);//or time(&start);
    //创建 fiducial model
    clucu_cosmology fiducial=clucu_cosmology_create(config_oguri,param_oguri,survey_oguri);
    fiducial.config.matter_power_spectrum_method=clucu_boltzmann_class_Pk;
    fiducial.runned_class = true;
    strcpy(fiducial.classname,"oguri");
	strcpy(fiducial.name,"oguri_cen");
    fiducial.survey.mass_ob_maxc_h = pow(10.,14.8);
    //测试期间改少bin数
//    fiducial.survey.mass_ob_max_h=pow(10.,14.6);
//    fiducial.survey.zobP_max=1.1;
	
	//创建一系列model
	int Np=33;
    //int Np=0;
	int Nid = Np*2+1;
	double *partial_param=Create1Grid(Np);
	clucu_cosmology *cosmo =(clucu_cosmology *) malloc(sizeof(clucu_cosmology)*Nid);
    
	cosmo[0]=fiducial;
	cosmo[0].cosmo_id = 0;
	cosmo[0].class_id = 0;
	InitializeCosmo(&cosmo[0]);
    clucu_compute_background(&cosmo[0],clucu_background_class_label,&status);
    clucu_compute_power(&cosmo[0],&status);

    
	double delta,partial;
	int id1,id2;
    double lnx,lnx1,lnx2;
	for(int number_p=0;number_p<Np;number_p++)
	{
		delta=1.0e-3;
		id2=number_p*2+2;
		id1=number_p*2+1;
		cosmo[id1]=fiducial;
		cosmo[id2]=fiducial;
		cosmo[id1].cosmo_id = id1;
		cosmo[id2].cosmo_id = id2;
		cosmo[id1].class_id = id1;
		cosmo[id2].class_id = id2;
		switch(number_p)
		{
            case 0:cosmo[id1].param.Omega_l *= 1-delta;cosmo[id2].param.Omega_l *= 1+delta;partial=2.0*fiducial.param.Omega_l*delta;break;
            case 1:cosmo[id1].param.w0 *= 1-delta;cosmo[id2].param.w0 *= 1+delta;partial=2.0*fiducial.param.w0*delta;break;
            case 2:cosmo[id1].param.wa -= delta;cosmo[id2].param.wa += delta;partial=2.0*delta;break;           
            case 3:
            {
                lnx = log(cosmo[0].param.Omh2);
                lnx1 =lnx*(1-delta);
                lnx2 =lnx*(1+delta);
                cosmo[id1].param.Omh2 = exp(lnx1);
                cosmo[id2].param.Omh2 = exp(lnx2);
                partial=2.0*lnx*delta;
                break;
            }
			case 4:
            {
                lnx = log(cosmo[0].param.Obh2);
                lnx1 =lnx*(1-delta);
                lnx2 =lnx*(1+delta);
                cosmo[id1].param.Obh2 = exp(lnx1);
                cosmo[id2].param.Obh2 = exp(lnx2);
                partial=2.0*lnx*delta;
                break;
            }
            case 5:cosmo[id1].param.n_s *= 1-delta;cosmo[id2].param.n_s *= 1+delta;partial=2.0*fiducial.param.n_s*delta;break;
            case 6:cosmo[id1].param.alpha_s -= delta;cosmo[id2].param.alpha_s += delta;partial=2.0*delta;break;
            case 7:
            {
                lnx = log(cosmo[0].param.delta_zeta);
                lnx1 =lnx*(1-delta);
                lnx2 =lnx*(1+delta);
                cosmo[id1].param.delta_zeta = exp(lnx1);
                cosmo[id2].param.delta_zeta = exp(lnx2);
                partial=2.0*lnx*delta;
                break;
            }

            case 8:cosmo[id1].param.f_nl -= 10.;cosmo[id2].param.f_nl += 10.;partial=2.0*10.;CASEEND;//hmf 开始用到。算sigma时还不用

            case 9:CASEEND0(lnM_b0,1.0e-2);//使星系团数，1e-1
            case 10:CASEEND0(q_b[0],1.0e-2);//使星系团数，1e-2
            case 11:CASEEND0(q_b[1],1.0e-2);//使星系团数，1e-2
            case 12:CASEEND0(q_b[2],1.0e-3);//使星系团数，1e-3
            case 13:CASEEND0(s_b[0],0.1);//使星系团数，1e-1
            case 14:CASEEND0(s_b[1],0.1);//使星系团数，1e-1
            case 15:CASEEND0(s_b[2],0.1);//使星系团数，1e-1
            case 16:CASEEND1(sigma_lnM0);
            case 17:CASEEND0(q_sigma_lnm[0],1.0e-2);//使sigma>0，1e-2；星系团数，1e-2
            case 18:CASEEND0(q_sigma_lnm[1],1.0e-3);//使sigma>0，1e-4；星系团数，1e-3
            case 19:CASEEND0(q_sigma_lnm[2],1.0e-4);//使sigma>0，1e-4；星系团数，1e-4
            case 20:CASEEND0(s_sigma_lnm[0],0.1);//使星系团数，1e-1
            case 21:CASEEND0(s_sigma_lnm[1],0.1);//使星系团数，1e-1
            case 22:CASEEND0(s_sigma_lnm[2],0.1);//使星系团数，1e-1
            //算number counts，Chh, 到这， 

            case 23:CASEEND1(source_redshift);// kk只跟这个参数有关，只能限制宇宙学参数+这1个。hk2到这。

            case 24:CASEEND1(A_vir);
            case 25:CASEEND1(B_vir);
            case 26:CASEEND1(C_vir);//

            case 27:CASEEND1(f_cen0);
            case 28:CASEEND1(p_cenM);
            case 29:CASEEND0(p_cenz,0.01);
            case 30:CASEEND1(sigma_s0);
            case 31:CASEEND0(p_sigmaM,0.01);
            case 32:CASEEND0(p_sigmaz,0.01);//hk1到这。


		}
		partial_param[number_p] = partial;
		InitializeCosmo(&cosmo[id1]);
		InitializeCosmo(&cosmo[id2]);
	}

	int Nz=cosmo[0].data.NzobP-1;
	int Nk1=cosmo[0].data.Nk1;
	int Nk2=cosmo[0].data.Nk2;
	int Nm=cosmo[0].data.Nmob-1;
    int Nell=cosmo[0].data.Nell;

//	int thread_num=Nid-1;
	int thread_num=5;
	//if(CHECK) thread_num=1;
	omp_set_num_threads(thread_num);
	#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
		for(int id = thread_id+1;id<Nid;id+=thread_num)
		{
			clucu_compute_background(&cosmo[id],clucu_background_class_label,&status);
            clucu_compute_power(&cosmo[id],&status);
		}
	}


	for(int id=0;id<Nid;id++)
	{
        cosmo[id].data.shiftperp[0]=clucu_angular_diameter_distance(&cosmo[0],0.0001,&status)/clucu_angular_diameter_distance(&cosmo[id],0.0001,&status);
		cosmo[id].data.shiftpara[0]=clucu_hubble_parameter(&cosmo[id],0,&status)/clucu_hubble_parameter(&cosmo[0],0,&status);
        for(int number_zobP=1;number_zobP<Nz;number_zobP++)
        {
            double zobP=cosmo[0].data.zobP[number_zobP];
            cosmo[id].data.shiftperp[number_zobP]=clucu_angular_diameter_distance(&cosmo[0],zobP,&status)/clucu_angular_diameter_distance(&cosmo[id],zobP,&status);
            cosmo[id].data.shiftpara[number_zobP]=clucu_hubble_parameter(&cosmo[id],zobP,&status)/clucu_hubble_parameter(&cosmo[0],zobP,&status);
		}		
	}
    // int iiid=31;
    // clucu_compute_logsigma(&cosmo[iiid],&cosmo[iiid].status);//功率谱相同即可
    // clucu_compute_cluster_numbercounts(&cosmo[iiid],1,1,&cosmo[iiid].status);
    // double sum=0.;
    // for(int iz=0;iz<Nz;iz++)
    // {
    //     for(int im=0;im<Nm;im++)
    //         sum+=cosmo[iiid].data.Nbin_mzob[im][iz];
    // }
    // printf("2.5e4[%d]=%e\n",iiid,sum);

	thread_num=Nid;
	if(CHECK) thread_num=1;
	omp_set_num_threads(thread_num);
	#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
		for(int id = thread_id;id<Nid;id+=thread_num)
		{
			cosmo[id].thread_num=thread_num;
            clucu_compute_sigmaS3(&cosmo[id],&cosmo[id].status);//拟合式
            clucu_compute_logsigma(&cosmo[id],&cosmo[id].status);//功率谱相同即可
            printf("log sigma --over\n");
            clucu_compute_kernel_n(&cosmo[id],&cosmo[id].status);//需要mass-ob 
            printf("kernel_n --over\n");
            clucu_compute_kernel_nb(&cosmo[id],&cosmo[id].status);
            printf("kernel_nb --over\n");
            clucu_compute_barn(&cosmo[id],&cosmo[id].status);
            printf("barn --over\n");
            // clucu_compute_cluster_averagebias(&cosmo[id],&cosmo[id].status);
            // printf("averagebias --over\n");
		}
		printf("exit#%d\n", thread_id);
	}

    // //测试NN
    // if(CHECK) thread_num=1;
	// omp_set_num_threads(thread_num);
	// #pragma omp parallel
	// {
	// 	int thread_id = omp_get_thread_num();
	// 	for(int id = thread_id;id<Nid;id+=thread_num)
	// 	{
	// 		cosmo[id].thread_num=thread_num;
    //         clucu_compute_cluster_numbercounts(&cosmo[id],1,1,&cosmo[id].status);
    //         printf("NN --over\n");
	// 	}
	// 	printf("exit#%d\n", thread_id);
	// }
    // for(int id = 0;id<Nid;id++){
    //     double sum=0.;
    //     for(int iz=0;iz<Nz;iz++)
    //     {
    //         for(int im=0;im<Nm;im++)
    //             sum+=cosmo[id].data.Nbin_mzob[im][iz];
    //     }
    //     printf("[%d]=%e\n",id,sum);
    // }

    // clucu_compute_samplingvariance(&cosmo[0],&cosmo[0].status);

    
    //测试hh
    if(CHECK) thread_num=1;
	omp_set_num_threads(thread_num);
	#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
		for(int id = thread_id;id<Nid;id+=thread_num)
		{
			cosmo[id].thread_num=thread_num;
            clucu_compute_cl_hh(&cosmo[id],&cosmo[id].status);//40s
            clucu_compute_average_one_over_chis(&cosmo[id]);
            clucu_compute_cl_kk(&cosmo[id],&cosmo[id].status);
            clucu_compute_cl_hk2(&cosmo[id],1,1,&cosmo[id].status);
            clucu_compute_cl_hk1(&cosmo[id],1,1,&cosmo[id].status);
            //clucu_compute_cov_hhhh(&cosmo[id],&cosmo[id].status);//3min。太久了，为啥，是不是循环太多了
            printf("hh --over\n");
		}
		printf("exit#%d\n", thread_id);
	}
    
    // //set data. N
    // double ***partial_N=Create3Grid(Np,Nz,Nm);
    // for(int number_p=0;number_p<Np;number_p++){
	// 	id2=number_p*2+2;
	// 	id1=number_p*2+1;
    //     for(int iz=0;iz<Nz;iz++){
    //         for(int im=0;im<Nm;im++){
    //             partial_N[number_p][iz][im] =  ( cosmo[id2].data.Nbin_mzob[im][iz] - cosmo[id1].data.Nbin_mzob[im][iz] )/partial_param[number_p];
    //         }
    //     }
    // }
    // //每个z，都有一个cov,
    // double ***CovN = Create3Grid(Nz,Nm,Nm);
    // double ***invCovN = Create3Grid(Nz,Nm,Nm);
    // for(int iz=0;iz<Nz;iz++){
    //     for(int im=0;im<Nm;im++){
    //         for(int jm=0;jm<Nm;jm++){
    //             CovN[iz][im][jm] = cosmo[0].data.Nbin_mzob[im][iz] * kronecker(im,jm) + cosmo[0].data.CovNN_SV[im][jm][iz];
    //         }
    //     }
    //     clucu_matrix_inv(Nm,CovN[iz],invCovN[iz]);
    // }
    // //对z求sum
    // double **fisherN = Create2Grid(Np,Np);
    // for(int p1=0;p1<Np;p1++)
    // {
    //     for(int p2=0;p2<Np;p2++)
    //     {
    //         double sum=0.;
    //         for(int iz=0;iz<Nz;iz++){
    //             for(int i=0;i<Nm;i++){
    //                 for(int j=0;j<Nm;j++){
    //                     sum +=  partial_N[p1][iz][i] * invCovN[iz][i][j] * partial_N[p2][iz][j];
    //                 }
    //             }
    //         }
    //         fisherN[p1][p2] = sum;
    //     }
    // }
	
    

	//set data
    //cosmo->data.cl_hh[im][imm][iz][iell]
    //cosmo->data.cl_hk[im][iz][iell]
	double *****partial_C_hh=Create5Grid(Np,Nm,Nm,Nz,Nell);
    double ****partial_C_hk=Create4Grid(Np,Nm,Nz,Nell);
    int Nhh = Nz*Nm*Nm;
    int Nhk = Nz*Nm;
    int Ndata = Nhh + Nhk;
    double ***partial_C_data=Create3Grid(Np,Nell,Ndata);
	for(int number_p=0;number_p<Np;number_p++)
	{
		id2=number_p*2+2;
		id1=number_p*2+1;
		for(int iell=0;iell<Nell;iell++)
		{
            int idata=0;
            for(int iz=0;iz<Nz;iz++){
                for(int im=0;im<Nm;im++){
                    for(int imm=0;imm<Nm;imm++){
                        partial_C_hh[number_p][im][imm][iz][iell] = ( cosmo[id2].data.cl_hh[im][imm][iz][iell] - cosmo[id1].data.cl_hh[im][imm][iz][iell] )/partial_param[number_p];
                        partial_C_data[number_p][iell][idata] = partial_C_hh[number_p][im][imm][iz][iell];
                        idata++;
                        if(isnan(partial_C_data[number_p][iell][idata]))
                            CLUCU_RAISE_WARNING(11,"[%d][%d][%d]]:%le",number_p,iell,idata,partial_C_data[number_p][iell][idata]);
                    }
                }
            }
            for(int iz=0;iz<Nz;iz++){
                for(int im=0;im<Nm;im++){
                    partial_C_hk[number_p][im][iz][iell] = ( cosmo[id2].data.cl_hk[im][iz][iell] - cosmo[id1].data.cl_hk[im][iz][iell] )/partial_param[number_p];
                    partial_C_data[number_p][iell][idata] = partial_C_hk[number_p][im][iz][iell];
                    idata++;
                    if(isnan(partial_C_data[number_p][iell][idata]))
                            CLUCU_RAISE_WARNING(11,"[%d][%d][%d]]:%le",number_p,iell,idata,partial_C_data[number_p][iell][idata]);
                }
            }
		}
	}

    //每个ell，都有一个cov,
    double ***Cov = Create3Grid(Nell,Ndata,Ndata);
    double ***invCov = Create3Grid(Nell,Ndata,Ndata);
    //单独的hh
    double ***Covhh = Create3Grid(Nell,Nhh,Nhh);
    double ***invCovhh = Create3Grid(Nell,Nhh,Nhh);
    //单独的hk
    double ***Covhk = Create3Grid(Nell,Nhk,Nhk);
    double ***invCovhk = Create3Grid(Nell,Nhk,Nhk);

    if(CHECK) thread_num=1;
	omp_set_num_threads(thread_num);
	#pragma omp parallel
	{
        for(int iell=0;iell<Nell;iell++){
            double ell = cosmo[0].data.ell[iell];
            double delta_ell = cosmo[0].data.ell[iell+1]-ell;
            double prefactor = (4.*M_PI/cosmo[0].survey.Delta_Omega) / ((2.*ell+1)*delta_ell);
            //set hh-hh
            int iz,im,imm,jz,jm,jmm;
            for(int i=0;i<Nhh;i++){
                iz = i / (Nm*Nm);
                im = (i % (Nm*Nm)) / Nm;
                imm = (i % (Nm*Nm)) % Nm;
                for(int j=0;j<=i;j++){
                    jz = j / (Nm*Nm);
                    jm = (j % (Nm*Nm)) / Nm;
                    jmm = (j % (Nm*Nm)) % Nm;
                    Cov[iell][i][j]= prefactor * kronecker(iz,jz)
                                        *( cosmo[0].data.cl_hh_hat[im][imm][iz][iell] 
                                        * cosmo[0].data.cl_hh_hat[jm][jmm][iz][iell] 
                                        + cosmo[0].data.cl_hh_hat[im][jmm][iz][iell] 
                                        * cosmo[0].data.cl_hh_hat[jm][imm][iz][iell]);
                    Covhh[iell][i][j]=Cov[iell][i][j];
                    if(isnan(Cov[iell][i][j]))
                        CLUCU_RAISE_WARNING(11,"[%d][%d][%d][%d]]:%le",iell,i,j,Cov[iell][i][j]);
                }
            }
            //set hk-hh
            for(int i=Nhh;i<Ndata;i++){
                iz = (i-Nhh) / Nm;
                im = (i-Nhh) % Nm;
                for(int j=0;j<Nhh;j++){
                    jz = j / (Nm*Nm);
                    jm = (j % (Nm*Nm)) / Nm;
                    jmm = (j % (Nm*Nm)) % Nm;
                    Cov[iell][i][j]= prefactor * kronecker(iz,jz)
                                        *( cosmo[0].data.cl_hh_hat[im][jm][iz][iell] 
                                        * cosmo[0].data.cl_hk[jmm][iz][iell] 
                                        + cosmo[0].data.cl_hh_hat[im][jmm][iz][iell] 
                                        * cosmo[0].data.cl_hk[jm][iz][iell]);
                    if(isnan(Cov[iell][i][j]))
                        CLUCU_RAISE_WARNING(11,"[%d][%d][%d][%d]]:%le",iell,i,j,Cov[iell][i][j]);
                                                                                
                }
            }
            //set hk-hk
            for(int i=Nhh;i<Ndata;i++){
                iz = (i-Nhh) / Nm;
                im = (i-Nhh) % Nm;
                for(int j=Nhh;j<=i;j++){
                    jz = (j-Nhh) / Nm;
                    jm = (j-Nhh) % Nm;
                    Cov[iell][i][j] = prefactor 
                                    *( cosmo[0].data.cl_hh_hat[im][jm][iz][iell] 
                                    * cosmo[0].data.cl_kk_hat[iell] 
                                    * kronecker(iz,jz)
                                    + cosmo[0].data.cl_hk[im][iz][iell] 
                                    * cosmo[0].data.cl_hk[jm][jz][iell]);
                    Covhk[iell][i-Nhh][j-Nhh]=Cov[iell][i][j];
                    if(isnan(Cov[iell][i][j]))
                        CLUCU_RAISE_WARNING(11,"[%d][%d][%d][%d]]:%le",iell,i,j,Cov[iell][i][j]);
                    
                }
            }
            //补全
            for(int i=0;i<Ndata;i++){
                for(int j=i+1;j<Ndata;j++){
                    Cov[iell][i][j] = Cov[iell][j][i];
                }
            }
            for(int i=0;i<Nhh;i++){
                for(int j=i+1;j<Nhh;j++){
                    //printf("hh:[%d][%d][%d]=%e\n",iell,i,j,Covhh[iell][j][i]);
                    Covhh[iell][i][j] = Covhh[iell][j][i];
                }
            }
            for(int i=0;i<Nhk;i++){
                for(int j=i+1;j<Nhk;j++){
                    //printf("hk:[%d][%d][%d]=%e\n",iell,i,j,Covhk[iell][j][i]);
                    Covhk[iell][i][j] = Covhk[iell][j][i];
                }
            }
            clucu_matrix_inv(Ndata,Cov[iell],invCov[iell]);
            clucu_matrix_inv(Nhh,Covhh[iell],invCovhh[iell]);
            clucu_matrix_inv(Nhk,Covhk[iell],invCovhk[iell]);
            
            for(int i=Nhh;i<Ndata;i++){
                for(int j=Nhh;j<i;j++){
                    if(isnan(invCov[iell][i][j]))
                        CLUCU_RAISE_WARNING(11,"[%d][%d][%d][%d]]:%le",iell,i,j,invCov[iell][i][j]);
                }
            }
            
        }
    }
    /*
    FILE *fp2;
    char name2[100];
    double *ydata = clucu_linear_spacing(0,Ndata-1.,Ndata);
    double *yhh = clucu_linear_spacing(0,Nhh-1.,Nhh);
    double *yhk = clucu_linear_spacing(0,Nhk-1.,Nhk);
    sprintf(name2,"../output/%s/Cov.dat", cosmo[0].name);
    fp2=fopen(name2,"w");
    save_matrix(fp2,Cov[0],ydata,ydata,Ndata,Ndata);
    fclose(fp2);
    sprintf(name2,"../output/%s/Covhh.dat", cosmo[0].name);
    fp2=fopen(name2,"w");
    save_matrix(fp2,Covhh[0],yhh,yhh,Nhh,Nhh);
    fclose(fp2);
    sprintf(name2,"../output/%s/Covhk.dat", cosmo[0].name);
    fp2=fopen(name2,"w");
    save_matrix(fp2,Covhk[0],yhk,yhk,Nhk,Nhk);
    fclose(fp2);
    */

    double **fisherC = Create2Grid(Np,Np);
    double **fisherChh = Create2Grid(Np,Np);
    double **fisherChk = Create2Grid(Np,Np);
    //计算fisher matrix
    for(int p1=0;p1<Np;p1++)
    {
        for(int p2=0;p2<Np;p2++)
        {
            double sum=0.;
            for(int iell=0;iell<Nell;iell++){
                for(int i=0;i<Ndata;i++){
                    for(int j=0;j<Ndata;j++){
                        sum +=  partial_C_data[p1][iell][i] * invCov[iell][i][j] * partial_C_data[p2][iell][j];
                    }
                }
            }
            double sumhh=0.;
            for(int iell=0;iell<Nell;iell++){
                for(int i=0;i<Nhh;i++){
                    for(int j=0;j<Nhh;j++){
                        //这里不可以用invCov的子矩阵哦，取逆之后就不一样了哦
                        sumhh +=  partial_C_data[p1][iell][i] * invCovhh[iell][i][j] * partial_C_data[p2][iell][j];
                    }
                }
            }
            double sumhk=0.;
            for(int iell=0;iell<Nell;iell++){
                for(int i=0;i<Nhk;i++){
                    for(int j=0;j<Nhk;j++){
                        sumhk +=  partial_C_data[p1][iell][i+Nhh] * invCovhk[iell][i][j] * partial_C_data[p2][iell][j+Nhh];
                    }
                }
            }
            fisherC[p1][p2]=sum;
            fisherChh[p1][p2]=sumhh;
            fisherChk[p1][p2]=sumhk;
        }
    }
/*
    char name1[100];
    char afdsf[10]="8";
    double *ydata = clucu_linear_spacing(0,Ndata-1.,Ndata);
    FILE *fp1;
    for(int i=0;i<Np;i++){
        sprintf(name1,"../output/%s/partial%dp_%s.dat", cosmo[0].name,i,afdsf);
        fp1=fopen(name1,"w");
        save_matrix(fp1,partial_C_data[i],cosmo[0].data.ell,ydata,Nell,Ndata);
        fclose(fp1);
    }
    
    FILE *fp2;
    char name2[100];
    for(int i=0;i<21;i=i+10){
        sprintf(name2,"../output/%s/Cov%dl_%s.dat", cosmo[0].name,i,afdsf);
        fp2=fopen(name2,"w");
        save_matrix(fp2,Cov[i],ydata,ydata,Ndata,Ndata);
	    fclose(fp2);

        sprintf(name2,"../output/%s/invCov%dl_%s.dat", cosmo[0].name,i,afdsf);
        fp2=fopen(name2,"w");
        save_matrix(fp2,invCov[i],ydata,ydata,Ndata,Ndata);
	    fclose(fp2);
    }
*/


	char namedir[50];
    char nameInfo[100];
	sprintf(namedir,"../output/%s", cosmo[0].name);
	sprintf(nameInfo,"../output/%s/info.dat", cosmo[0].name);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameInfo,"a");
    end =time(NULL);
    fprintf(fp,"\n\ntime=%fs\n",difftime(end,start));
	fprintf(fp,"time=%fm\n",difftime(end,start)/60.);
	fprintf(fp,"time=%fh\n",difftime(end,start)/3600.);
	fclose(fp);

    char nameMatrix[100];
	sprintf(namedir,"../output/%s", cosmo[0].name);
	sprintf(nameMatrix,"../output/%s/matrix.dat", cosmo[0].name);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	//FILE* fp=fopen("./matrix.dat","w");
	freopen(nameMatrix,"w",stdout);

    //PrintMatrix(fisherN,Np,Np);
	PrintMatrix(fisherC,Np,Np);
    PrintMatrix(fisherChh,Np,Np);
    PrintMatrix(fisherChk,Np,Np);

    //double **covC = Create2Grid(Np,Np);
    //clucu_matrix_inv(Np,fisherC,covC);
    //PrintMatrix(covC,Np,Np);
    

}