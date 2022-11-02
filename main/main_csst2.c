#include "clucu.h"

#include <time.h>
#include <string.h>
// gcc main_csst.c  -I/home/mjchen/Program/cluster/newcode/include -I/home/mjchen/Program/gsl-2.7/gsl  -L/home/mjchen/Program/gsl-2.7/.libs -lgsl -lgslcblas -lm -fopenmp  -o test
/**
 * @brief 
 * 20220716: mass_ob = mass_true. sigma_lnM^2 = sigma_lnM0^2 -1 + (1+z)^2beta
 * 参数共有：Al,Bl,Cl,Dl,beta
 * 分别对应：mass_ob_A，mass_ob_B，mass_ob_Bz,mass_ob_sigma0,mass_ob_qz
 * 参数共有：B_M0,alpha,kappa
 * 分别对应：lnM_b0,s_b[0],mass_ob_qz
 * log10 mass_ob =[14,16], delta = 0.2
 * @return int 
 */
#define END1(x) cosmo[id1].x *= 1-delta;cosmo[id2].x *= 1+delta;partial=2.0*fiducial.x*delta;break
#define END0(x,dx) cosmo[id1].x -= dx;cosmo[id2].x += dx;partial=2.0*dx;break
#define CASEEND cosmo[id1].class_id=0;cosmo[id2].class_id=0;cosmo[id1].runned_class=true;cosmo[id2].runned_class=true;break
#define CASEEND1(x) cosmo[id1].x *= 1-delta;cosmo[id2].x *= 1+delta;partial=2.0*fiducial.x*delta;CASEEND
#define CASEEND0(x,dx) cosmo[id1].x -= dx;cosmo[id2].x += dx;partial=2.0*dx;CASEEND
int main()
{
    int status=0;
	omp_set_nested(1);
    //用于设置是否允许OpenMP进行嵌套并行，默认的设置为false。
	time_t start,end;
    start =time(NULL);//or time(&start);
    //创建 fiducial model
    clucu_cosmology fiducial=clucu_cosmology_create(config_default,param_default,survey_csst2);
    //对 fiducial model 稍作修改
    fiducial.config.mass_observable_method=clucu_my_csst;
    fiducial.survey.mass_ob_min_h=0.836*pow(10.,14.);
    fiducial.GISDB = 1;
	strcpy(fiducial.classname,"csst2");
	strcpy(fiducial.name,"csst2_1");
	fiducial.survey.redshift_sigma0 = 0.001;
	fiducial.runned_class = true;

	//创建一系列model
	int Np=15;
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
            case 0:END1(param.Obh2);
            case 1:END1(param.Omh2);
            case 2:END1(param.Onuh2_mass_sum);
            case 3:END1(param.Omega_l);
            case 4:END1(param.w0);
            case 5:END0(param.wa,delta);
            case 6:END1(param.n_s);
            case 7:END1(param.sigma8);
            case 8:CASEEND1(survey.mass_ob_A);
            case 9:CASEEND1(survey.mass_ob_B);
            case 10:CASEEND1(survey.mass_ob_Bz);
            case 11:CASEEND1(survey.mass_ob_sigma0);
            case 12:CASEEND1(survey.mass_ob_qz);
            case 13:CASEEND0(survey.lnM_b0,delta);
            case 14:CASEEND0(survey.s_b[0],delta);
		}
		partial_param[number_p] = partial;
		InitializeCosmo(&cosmo[id1]);
		InitializeCosmo(&cosmo[id2]);
	}

	int Nzob=cosmo[0].data.Nzob-1;
	int NzobP=cosmo[0].data.NzobP;
	int Nk1=cosmo[0].data.Nk1;
	int Nk2=cosmo[0].data.Nk2;
	int Nmob=cosmo[0].data.Nmob-1;

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
        for(int number_zobP=1;number_zobP<NzobP;number_zobP++)
        {
            double zobP=cosmo[0].data.zobP[number_zobP];
            cosmo[id].data.shiftperp[number_zobP]=clucu_angular_diameter_distance(&cosmo[0],zobP,&status)/clucu_angular_diameter_distance(&cosmo[id],zobP,&status);
            cosmo[id].data.shiftpara[number_zobP]=clucu_hubble_parameter(&cosmo[id],zobP,&status)/clucu_hubble_parameter(&cosmo[0],zobP,&status);
		}		
	}
    

	thread_num=Nid;
	if(CHECK) thread_num=1;
	omp_set_num_threads(thread_num);
	#pragma omp parallel
	{
		int thread_id = omp_get_thread_num();
		for(int id = thread_id;id<Nid;id+=thread_num)
		{
			cosmo[id].thread_num=thread_num;
            clucu_compute_logsigma(&cosmo[id],&status);
            printf("#%d:log sigma --over\n",id);
			clucu_compute_cluster_numbercounts(&cosmo[id],1,1,&status);
			printf("#%d:number counts --over\n",id);
			if(cosmo->GISDB==0)
                clucu_compute_cluster_averagebias(&cosmo[id],&status);
            else if(cosmo->GISDB==1)
                clucu_compute_cluster_averagebias_k(&cosmo[id],&status);
			printf("#%d:average bias --over\n",id);
            clucu_compute_cluster_powerspectrum(&cosmo[id],&status);
			printf("#%d:power spectrum--over\n",id);
			clucu_compute_cluster_volumeeffect(&cosmo[id],&status);
			printf("#%d:volume effect --over\n",id);
		}
		printf("exit#%d\n", thread_id);
	}
//    SaveNmzob(cosmo[0]);
	


	printf("\n\ncalculate fisherN, fisherPc, fisherNPc\n");
    double **fisherN = Create2Grid(Np,Np);//定义星系团计数的fisher matrix，一个8x8的矩阵
    double **fisherPc = Create2Grid(Np,Np);//定义星系团功率谱的fisher matrix，一个8x8的矩阵

	//计算N
	double ***partial_N=Create3Grid(Np,Nmob,Nzob);
	for(int number_p=0;number_p<Np;number_p++)
	{
		id2=number_p*2+2;
		id1=number_p*2+1;
		for(int number_zob=0;number_zob<Nzob;number_zob++)
		{
			for(int number_mob=0;number_mob<Nmob;number_mob++)
			{
				partial_N[number_p][number_mob][number_zob]= (  cosmo[id2].data.Nbin_mzob[number_mob][number_zob] - cosmo[id1].data.Nbin_mzob[number_mob][number_zob]  ) /partial_param[number_p];
			}
		}
	}
	for(int p1=0;p1<Np;p1++)
	{
		for(int p2=0;p2<Np;p2++)
		{
			double sumfisherN = 0;
			int zmin=1;
			if(cosmo[0].survey.redshift_sigma0>0.0) zmin=0;
			for(int number_zob=zmin;number_zob<Nzob;number_zob++)//这里不能从z=0开始，但LSST的最小红移不是零。
			{
				for(int number_mob=0;number_mob<Nmob;number_mob++)
				{
					sumfisherN=sumfisherN+partial_N[p1][number_mob][number_zob]*partial_N[p2][number_mob][number_zob]/cosmo[0].data.Nbin_mzob[number_mob][number_zob]; 
				}
			}
			fisherN[p1][p2]=sumfisherN;
		}
	}


	//计算P
	double aa,bb,cc,dd;
	double ****partial_kPc=Create4Grid(Np,Nk1,Nk2,NzobP);
	for(int number_p=0;number_p<Np;number_p++)
	{
		id2=number_p*2+2;
		id1=number_p*2+1;
		for(int number_zobP=0;number_zobP<NzobP;number_zobP++)
		{
			// 设置partial_kPc0
			for(int number_k1=0;number_k1<Nk1;number_k1++)//k1循环
			{
				for(int number_k2=0;number_k2<Nk2;number_k2++)//k2循环
				{
					aa=cosmo[id2].data.lnkPc[number_k1][number_k2][number_zobP];
					bb=cosmo[id1].data.lnkPc[number_k1][number_k2][number_zobP];
					cc=aa-bb;
					partial_kPc[number_p][number_k1][number_k2][number_zobP] = (cosmo[id2].data.lnkPc[number_k1][number_k2][number_zobP]-cosmo[id1].data.lnkPc[number_k1][number_k2][number_zobP] )/ partial_param[number_p];
				}
			}
		}
	}
	for(int p1=0;p1<Np;p1++)
    {
        for(int p2=0;p2<Np;p2++)
        {
            double sumfisherPc=0;
            for (int number_k1=0;number_k1<Nk1;number_k1++)//k1循环
            {
                for (int number_k2=0;number_k2<Nk2;number_k2++)//k2循环
                {
                    for (int number_zobP=0;number_zobP<NzobP;number_zobP++)
                    {
						aa=partial_kPc[p1][number_k1][number_k2][number_zobP];
						bb=partial_kPc[p2][number_k1][number_k2][number_zobP];
						cc=cosmo[0].data.V_k[number_k1][number_k2];
						dd=cosmo[0].data.V_eff[number_k1][number_k2][number_zobP];
						sumfisherPc=sumfisherPc+partial_kPc[p1][number_k1][number_k2][number_zobP]*partial_kPc[p2][number_k1][number_k2][number_zobP]*cosmo[0].data.V_k[number_k1][number_k2]*cosmo[0].data.V_eff[number_k1][number_k2][number_zobP]/2.;
					}
                }
            }
            fisherPc[p1][p2]=sumfisherPc;
        }
    }


	end =time(NULL);

	double sumN1=0.,sumN2=0.,sumN3=0.;
	for(int number_zob=0;number_zob<Nzob;number_zob++)
	{
		sumN1+=cosmo[0].data.Nbin_zob[number_zob];
	}
	for(int number_zob=0;number_zob<Nzob;number_zob++)
	{
		for(int number_mob=0;number_mob<Nmob;number_mob++)
			sumN2+=cosmo[0].data.Nbin_mzob[number_mob][number_zob];
	}
	for(int number_zobP=0;number_zobP<cosmo[0].data.NzobP;number_zobP++)
	{
		sumN3+=cosmo[0].data.Nbin_zobP[number_zobP];
	}
	printf("time=%fs\n",difftime(end,start));
	printf("time=%fm\n",difftime(end,start)/60.);
	printf("time=%fh\n",difftime(end,start)/3600.);
	printf("sumNbinzob=%f\nsumNlz=%f\nsumNzP=%f\n",sumN1,sumN2,sumN3);

	char namedir[50];
    char nameInfo[100];
	sprintf(namedir,"../output/%s", cosmo[0].name);
	sprintf(nameInfo,"../output/%s/info.dat", cosmo[0].name);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameInfo,"a");
    fprintf(fp,"time=%fs\n",difftime(end,start));
	fprintf(fp,"time=%fm\n",difftime(end,start)/60.);
	fprintf(fp,"time=%fh\n",difftime(end,start)/3600.);
	fprintf(fp,"sumNbinzob=%f\nsumNlz=%f\nsumNzP=%f\n",sumN1,sumN2,sumN3);
	fclose(fp);

    char nameMatrix[100];
	sprintf(namedir,"../output/%s", cosmo[0].name);
	sprintf(nameMatrix,"../output/%s/matrix.dat", cosmo[0].name);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	//FILE* fp=fopen("./matrix.dat","w");
	freopen(nameMatrix,"w",stdout);

	PrintMatrix(fisherN,Np,Np);
	PrintMatrix(fisherPc,Np,Np);
}


