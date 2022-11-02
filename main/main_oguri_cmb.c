#include "clucu.h"

#include <time.h>
#include <string.h>

// gcc main_csst.c  -I/home/mjchen/Program/cluster/newcode/include -I/home/mjchen/Program/gsl-2.7/gsl  -L/home/mjchen/Program/gsl-2.7/.libs -lgsl -lgslcblas -lm -fopenmp  -o test
/**
 * @brief 
 * 20220716: mass_ob = mass_true. sigma_lnM^2 = sigma_lnM0^2 -1 + (1+z)^2beta
 * 参数共有：Al,Bl,Cl,Dl,beta
 * 分别对应：mass_ob_A，mass_ob_B，mass_ob_Bz,mass_ob_sigma0,mass_ob_qz
 * log10 mass_ob =[14,16], delta = 0.2
 * @return int 
 */
void RunClassCMB(clucu_cosmology *cosmo)
{
    if(cosmo->runned_class == true) return;  
	char namedir[50];
    char nameclassini[50];
    char cmd[100];
	sprintf(namedir,"../data/%s", cosmo->classname);
	mkdir("../data", 0755);
	mkdir(namedir, 0755);
    sprintf(nameclassini,"../data/%s/p%d.ini", cosmo->classname, cosmo->class_id);
    sprintf(cmd," ../src/class  ../data/%s/p%d.ini", cosmo->classname, cosmo->class_id);


    //创建ini文件
    FILE* fp1=fopen(nameclassini,"w");
    
    //基本参数
//	fprintf(fp1,"h = %.15lf  \nT_cmb =%lf  \nomega_b = %lf   \nomega_cdm = %lf    \n",
//                    cosmo->param.h,cosmo->param.T_CMB,  cosmo->param.Obh2,  cosmo->param.Och2);
    fprintf(fp1,"h = %.15lf  \nOmega_b = %lf   \nOmega_cdm = %lf    \n",
                    cosmo->param.h,  cosmo->param.Omega_b,  cosmo->param.Omega_c);
    fprintf(fp1,"Omega_k = %lf    \nT_cmb =%lf  \n",cosmo->param.Omega_k,cosmo->param.T_CMB);
    
    //暗能量参数
    if(!isnan(cosmo->param.w0))
        fprintf(fp1,"Omega_Lambda = 0   \nfluid_equation_of_state = CLP     \nw0_fld = %lf    \nwa_fld = %lf    \n",cosmo->param.w0,cosmo->param.wa);
    else
        fprintf(fp1,"Omega_Lambda = %lf     \n",cosmo->param.Omega_l);
    
    //功率谱参数
    fprintf(fp1,"k_pivot=%lf  \nn_s = %lf    \nalpha_s = %lf    \n",cosmo->param.k_pivot,cosmo->param.n_s,cosmo->param.alpha_s);
    if(isnan(cosmo->param.sigma8) && !isnan(cosmo->param.A_s))
        fprintf(fp1,"A_s = %le    \n",cosmo->param.A_s);
    else if(!isnan(cosmo->param.sigma8) && isnan(cosmo->param.A_s))
        fprintf(fp1,"sigma8 = %lf    \n",cosmo->param.sigma8);
    
    //中微子参数
//    if(cosmo->param.N_nu_mass > 0)
//        fprintf(fp1,"N_ur = 2.0328       \nN_ncdm = 1        \nm_ncdm = %.15lf\n",cosmo->param.sum_nu_masses);
    if(cosmo->param.N_nu_mass > 0)
    {
        fprintf(fp1,"N_ur = %lf    \nN_ncdm = %d        \n\n",cosmo->param.N_nu_rel,cosmo->param.N_nu_mass);
        fprintf(fp1,"m_ncdm = ");
        for(int i=0;i<cosmo->param.N_nu_mass-1;i++)
            fprintf(fp1,"%lf,",cosmo->param.m_nu[i]);
        fprintf(fp1,"%lf\n",cosmo->param.m_nu[cosmo->param.N_nu_mass-1]);
    }
    
    fprintf(fp1,"tau_reio = %lf    \n",cosmo->param.tau_reio);
    fprintf(fp1,"\noutput = tCl,pCl  \nroot = ../data/%s/p%d_    \n",cosmo->classname, cosmo->class_id);
    fprintf(fp1,"format = camb\n");

    fclose(fp1);

    //调用class，得到输出文件。这里只读取H(z)
    system(cmd);
    cosmo->runned_class = true;
    return;
}
void LoadCMB(clucu_cosmology *cosmo,int *status)
//H,dC,dA,dL,
{   
    if(cosmo->runned_class==false)
    {
        RunClassCMB(cosmo);
    }
    char nameclassbg[150];
    sprintf(nameclassbg,"../data/%s/p%d_cl.dat", cosmo->classname, cosmo->class_id);

    Free1Grid(cosmo->data.ell,cosmo->data.Nell);
    int Nell=1991;
    cosmo->data.Nell = Nell;
    cosmo->data.ell = NULL;
    cosmo->data.ell = Create1Grid(cosmo->data.Nell);
    cosmo->data.CTT =Create1Grid(cosmo->data.Nell);
    cosmo->data.CTT_hat =Create1Grid(cosmo->data.Nell);
    cosmo->data.CEE =Create1Grid(cosmo->data.Nell);
    cosmo->data.CEE_hat =Create1Grid(cosmo->data.Nell);
    cosmo->data.CTE =Create1Grid(cosmo->data.Nell);
    cosmo->data.CovTT =Create1Grid(cosmo->data.Nell);
    cosmo->data.CovEE =Create1Grid(cosmo->data.Nell);
    cosmo->data.CovTE =Create1Grid(cosmo->data.Nell);

    double f_sky = 0.65;
    double theta_nu1=8.0;//单位arcmin
    double theta_nu2=5.5;
    double sigmaT_nu1=5.2;//单位uK
    double sigmaT_nu2=11.7;
    double sigmaP_nu1=10.8;
    double sigmaP_nu2=24.3;
    theta_nu1 *= M_PI/(60.*180.);//换成单位1
    theta_nu2 *= M_PI/(60.*180.);
    double wT_nu1 = pow(theta_nu1*sigmaT_nu1,-2.);//单位1/uK^2
    double wT_nu2 = pow(theta_nu2*sigmaT_nu2,-2.);
    double wP_nu1 = pow(theta_nu1*sigmaP_nu1,-2.);
    double wP_nu2 = pow(theta_nu2*sigmaP_nu2,-2.);
    double B2_nu1,B2_nu2,noiseT,noiseP;

    FILE* fp=fopen(nameclassbg,"r");

    double ctt,cee,cte,ctth,ceeh;
    //int N_Abandon = 15;
    int N_Abandon = 14;
    for(int table_i=0;table_i<N_Abandon+Nell;table_i++)
    {
        char abandon[1000];
        if(table_i<N_Abandon)
            fgets(abandon,1000,fp);
        else
        {
            int i =table_i-N_Abandon;
            fscanf(fp,"%lf  %le  %le  %*le  %le\n",&cosmo->data.ell[i],&cosmo->data.CTT[i],&cosmo->data.CEE[i],&cosmo->data.CTE[i]);
            //CLASS输出的功率谱带l，化成不带l
            //CLASS输出的功率谱单位为1，化成uK^2
            //cosmo->data.CTT[i] *= 2.*M_PI/(cosmo->data.ell[i]*(cosmo->data.ell[i]+1) ) *pow(cosmo->param.T_CMB * 1.0e6,2.);
            cosmo->data.CTT[i] *= 2.*M_PI/(cosmo->data.ell[i]*(cosmo->data.ell[i]+1) );
            cosmo->data.CEE[i] *= 2.*M_PI/(cosmo->data.ell[i]*(cosmo->data.ell[i]+1) );
            cosmo->data.CTE[i] *= 2.*M_PI/(cosmo->data.ell[i]*(cosmo->data.ell[i]+1) );
            B2_nu1 = exp(  -cosmo->data.ell[i]*(cosmo->data.ell[i]+1)  *theta_nu1*theta_nu1 /8./log(2.)  ); 
            B2_nu2 = exp(  -cosmo->data.ell[i]*(cosmo->data.ell[i]+1)  *theta_nu2*theta_nu2 /8./log(2.)  ); 
            noiseT = 1./(wT_nu1 * B2_nu1 + wT_nu2 * B2_nu2);//单位1/uK^2
            noiseP = 1./(wP_nu1 * B2_nu1 + wP_nu2 * B2_nu2);
            cosmo->data.CTT_hat[i] = cosmo->data.CTT[i] + noiseT;
            cosmo->data.CEE_hat[i] = cosmo->data.CEE[i] + noiseP;

            cosmo->data.CovTT[i] = pow(cosmo->data.CTT_hat[i],2.)   *2./(2.*cosmo->data.ell[i]+1.)/f_sky ;
            cosmo->data.CovEE[i] = pow(cosmo->data.CEE_hat[i],2.)   *2./(2.*cosmo->data.ell[i]+1.)/f_sky ;
            cosmo->data.CovTE[i] = pow(cosmo->data.CTE[i],2.)   *2./(2.*cosmo->data.ell[i]+1.)/f_sky ;


            ctt=cosmo->data.CTT[i];
            cee=cosmo->data.CEE[i];
            cte=cosmo->data.CTE[i];
            ctth=cosmo->data.CTT_hat[i];
            ceeh=cosmo->data.CEE_hat[i];
            }
        }
    fclose(fp);
}
void SaveC_cmb(clucu_cosmology *cosmo)
{
	char namedir[50];
    char nameNlzob[100];
	sprintf(namedir,"../output/%s", cosmo->name);
	sprintf(nameNlzob,"../output/%s/C%d_cmb.dat", cosmo->name,cosmo->cosmo_id);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameNlzob,"w");
	for(int i=0;i<cosmo->data.Nell;i++)
		fprintf(fp,"%f      %le      %le      %le\n",cosmo->data.ell[i],cosmo->data.CTT[i],cosmo->data.CEE[i],cosmo->data.CTE[i]);
	fclose(fp);
}
void SaveCov_cmb(clucu_cosmology *cosmo)
{
	char namedir[50];
    char nameNlzob[100];
	sprintf(namedir,"../output/%s", cosmo->name);
	sprintf(nameNlzob,"../output/%s/Cov%d_cmb.dat", cosmo->name,cosmo->cosmo_id);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	FILE *fp=fopen(nameNlzob,"w");
	for(int i=0;i<cosmo->data.Nell;i++)
		fprintf(fp,"%f      %le      %le      %le\n",cosmo->data.ell[i],cosmo->data.CovTT[i],cosmo->data.CovEE[i],cosmo->data.CovTE[i]);
	fclose(fp);
}

int main()
{
    int status=0;
	omp_set_nested(1);
    //用于设置是否允许OpenMP进行嵌套并行，默认的设置为false。
	time_t start,end;
    start =time(NULL);//or time(&start);
    //创建 fiducial model
    clucu_cosmology fiducial=clucu_cosmology_create(config_oguri,param_oguri,survey_oguri);
	strcpy(fiducial.classname,"oguri_cmb");
	strcpy(fiducial.name,"oguri_cmb");
    fiducial.runned_class = true;

	//创建一系列model
	int Np=9;
	int Nid = Np*2+1;
	double *partial_param=Create1Grid(Np);
	clucu_cosmology *cosmo =(clucu_cosmology *) malloc(sizeof(clucu_cosmology)*Nid);
    
	cosmo[0]=fiducial;
	cosmo[0].cosmo_id = 0;
	cosmo[0].class_id = 0;
	InitializeCosmo(&cosmo[0]);
    RunClassCMB(&cosmo[0]);
    LoadCMB(&cosmo[0],&status);
    SaveC_cmb(&cosmo[0]);
    SaveCov_cmb(&cosmo[0]);
    


    
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
            case 8:cosmo[id1].param.tau_reio -= delta;cosmo[id2].param.tau_reio += delta;partial=2.0*delta;break;
		}
		partial_param[number_p] = partial;
		InitializeCosmo(&cosmo[id1]);
		InitializeCosmo(&cosmo[id2]);
        RunClassCMB(&cosmo[id1]);
        RunClassCMB(&cosmo[id2]);
        LoadCMB(&cosmo[id1],&status);
        LoadCMB(&cosmo[id2],&status);
	}

    int Nell=cosmo[0].data.Nell;

    double **fisherC = Create2Grid(Np,Np);//定义星系团计数的fisher matrix，一个8x8的矩阵

	//计算N
	double **partial_C_TT=Create2Grid(Np,Nell);
    double **partial_C_EE=Create2Grid(Np,Nell);
    double **partial_C_TE=Create2Grid(Np,Nell);
	for(int number_p=0;number_p<Np;number_p++)
	{
		id2=number_p*2+2;
		id1=number_p*2+1;
		for(int number_ell=0;number_ell<Nell;number_ell++)
		{
			partial_C_TT[number_p][number_ell]= (  cosmo[id2].data.CTT[number_ell] - cosmo[id1].data.CTT[number_ell]  ) /partial_param[number_p];
            partial_C_EE[number_p][number_ell]= (  cosmo[id2].data.CEE[number_ell] - cosmo[id1].data.CEE[number_ell]  ) /partial_param[number_p];
            partial_C_TE[number_p][number_ell]= (  cosmo[id2].data.CTE[number_ell] - cosmo[id1].data.CTE[number_ell]  ) /partial_param[number_p];
		}
	}

    double f_sky = 0.65;
    for(int p1=0;p1<Np;p1++)
    {
        for(int p2=0;p2<Np;p2++)
        {
            double sum=0.;
            for(int l=0;l<Nell;l++)
            {
                gsl_matrix *mat_cl= gsl_matrix_alloc (2, 2);
                gsl_matrix *matder1= gsl_matrix_alloc (2, 2);
                gsl_matrix *matder2= gsl_matrix_alloc (2, 2);

                gsl_matrix_set (mat_cl, 0, 0, cosmo[0].data.CTT_hat[l]);
                gsl_matrix_set (mat_cl, 0, 1, cosmo[0].data.CTE[l]);
                gsl_matrix_set (mat_cl, 1, 0, cosmo[0].data.CTE[l]);
                gsl_matrix_set (mat_cl, 1, 1, cosmo[0].data.CEE_hat[l]);

                gsl_matrix_set (matder1, 0, 0, partial_C_TT[p1][l]);
                gsl_matrix_set (matder1, 0, 1, partial_C_TE[p1][l]);
                gsl_matrix_set (matder1, 1, 0, partial_C_TE[p1][l]);
                gsl_matrix_set (matder1, 1, 1, partial_C_EE[p1][l]);

                gsl_matrix_set (matder2, 0, 0, partial_C_TT[p2][l]);
                gsl_matrix_set (matder2, 0, 1, partial_C_TE[p2][l]);
                gsl_matrix_set (matder2, 1, 0, partial_C_TE[p2][l]);
                gsl_matrix_set (matder2, 1, 1, partial_C_EE[p2][l]);

                gsl_matrix *inv_cl = gsl_matrix_inv(mat_cl);//取逆
                gsl_matrix *A = gsl_matrix_dot(inv_cl,matder1);
                gsl_matrix *B = gsl_matrix_dot(inv_cl,matder2);
                gsl_matrix *C = gsl_matrix_dot(A,B);
                sum += gsl_matrix_trace(C) *(2*cosmo[0].data.ell[l]+1)*f_sky/2.;
                
                gsl_matrix_free (mat_cl);
                gsl_matrix_free (matder1);
                gsl_matrix_free (matder2);
                gsl_matrix_free (inv_cl);
                gsl_matrix_free (A);
                gsl_matrix_free (B);
                gsl_matrix_free (C);

            }
            fisherC[p1][p2]=sum;
        }
    }

    char namedir[50];
    char nameMatrix[100];
	sprintf(namedir,"../output/%s", cosmo[0].name);
	sprintf(nameMatrix,"../output/%s/matrix.dat", cosmo[0].name);
	mkdir("../output", 0755);
	mkdir(namedir, 0755);
	//FILE* fp=fopen("./matrix.dat","w");
	freopen(nameMatrix,"w",stdout);

	PrintMatrix(fisherC,Np,Np);

}


