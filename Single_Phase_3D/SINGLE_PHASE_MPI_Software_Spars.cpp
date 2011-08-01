#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<math.h> 
# include "mpi.h"



#define MAXN 100
#define eps 1e-12
#define zero(x) (fabs(x)<eps)

struct mat{
    int n,m;
    double data[MAXN][MAXN];
};



//D3Q19 STANDARD LATTICE MRT LATTICE BOLTZMANN METHOD
//UNIFORM MESH, LATTICE VELOCITY 1


using namespace std;        
const int Q=19;          


double u_max,u_ave,gx,gy,gz,porosity;

//----------
double s_e;
double s_eps;
double s_q;
//----------

double s_v;
double q_p;

int cl,cr;

int EI;


int Count,NX,NY,NZ;

int mirX=0;
int mirY=0;
int mirZ=0;
int mir=0;


double M[19][19]=
{{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
{-30,-11,-11,-11,-11,-11,-11,8,8,8,8,8,8,8,8,8,8,8,8},
{12,-4,-4,-4,-4,-4,-4,1,1,1,1,1,1,1,1,1,1,1,1},
{0,1,-1,0,0,0,0,1,-1,1,-1,0,0,0,0,1,-1,1,-1},
{0,-4,4,0,0,0,0,1,-1,1,-1,0,0,0,0,1,-1,1,-1},
{0,0,0,1,-1,0,0,1,1,-1,-1,1,-1,1,-1,0,0,0,0},
{0,0,0,-4,4,0,0,1,1,-1,-1,1,-1,1,-1,0,0,0,0},
{0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1},
{0,0,0,0,0,-4,4,0,0,0,0,1,1,-1,-1,1,1,-1,-1},
{0,2,2,-1,-1,-1,-1,1,1,1,1,-2,-2,-2,-2,1,1,1,1},
{0,-4,-4,2,2,2,2,1,1,1,1,-2,-2,-2,-2,1,1,1,1},
{0,0,0,1,1,-1,-1,1,1,1,1,0,0,0,0,-1,-1,-1,-1},
{0,0,0,-2,-2,2,2,1,1,1,1,0,0,0,0,-1,-1,-1,-1},
{0,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1},
{0,0,0,0,0,0,0,1,-1,1,-1,0,0,0,0,-1,1,-1,1},
{0,0,0,0,0,0,0,-1,-1,1,1,1,-1,1,-1,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,1,1,1,-1,-1}};


double MI[19][19];







double m[19];
double meq[19];





double uMax,c,Re,dx,dy,Lx,Ly,dt,rho0,P0,tau_f,niu,error,SFx,SFy,reso;

void Read_Rock(int***,double*,char[128]);

void tests();

void init_Sparse(int***,int***,int*, int*);

void init(double*, double**, double**,double*, double*, double*);

void periodic_streaming(double** ,double** ,int* ,int***,int*, int*,double*, double**);

void periodic_streaming_MR(double** ,double** ,int* ,int*** ,int* ,int* ,double* ,double** );

void standard_bounceback_boundary(int,double**);

void collision(double*,double** ,double** ,double*, double* ,double* , int* ,int***);

void comput_macro_variables( double* ,double**,double** ,double** ,double** ,double* , double* , double* ,int* ,int***);

void comput_macro_variables_IMR( double* ,double** ,double** ,double**,double** ,int*,int*** , double* , double*, double*, int ,double ,double* ,double* , double* );

double Error(double** ,double** ,double*, double*);

void boundary_velocity(int,double,int , double ,int ,double ,int ,double ,int ,double ,int ,double ,double** ,double* ,double** ,int*** );

void boundary_pressure(int ,double ,int , double ,int ,double ,int ,double ,int ,double ,int ,double ,double** ,double** ,double* ,int*** );

void output_velocity(int ,double* ,double** ,int ,int ,int ,int ,int*** );

void output_density(int ,double* ,int ,int ,int ,int ,int*** );	

void Geometry(int*** );	

double Comput_Perm(double** u,double*,int);

double S[19];

void Comput_MI(double[19][19], double[19][19]);
int inverse(mat &a);
double feq(int,double, double[3]);

void Suppliment(int*,int***);

int e[19][3]=
{{0,0,0},{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};
double w[19]={1.0/3.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};

int LR[19]={0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15};

int n,nx_l,n_max,in_IMR,PerDir,freRe,freDe,freVe,Par_Geo,Par_nx,Par_ny,Par_nz;
int Zoom;


int wr_per,pre_xp,pre_xn,pre_yp,pre_yn,pre_zp,pre_zn;
int vel_xp,vel_xn,vel_yp,vel_yn,vel_zp,vel_zn;
double in_vis,p_xp,p_xn,p_yp,p_yn,p_zp,p_zn;
double inivx,inivy,inivz,v_xp,v_xn,v_yp,v_yn,v_zp,v_zn;
double error_perm;

char outputfile[128]="./";

int main(int argc , char *argv [])
{	

MPI :: Init (argc , argv );
MPI_Status status ;

double start , finish;

int rank = MPI :: COMM_WORLD . Get_rank ();
int para_size=MPI :: COMM_WORLD . Get_size ();

int dif;
 

	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	int NCHAR=128;
	char     filename[128], dummy[128+1];
	int      dummyInt;

	if (rank==0)
	{
	ifstream fin(argv[1]);
	
	fin >> filename;				fin.getline(dummy, NCHAR);
	fin >> NX >> NY >> NZ;				fin.getline(dummy, NCHAR);
	fin >> n_max;					fin.getline(dummy, NCHAR);
	fin >> reso;					fin.getline(dummy, NCHAR);
	fin >> in_IMR;					fin.getline(dummy, NCHAR);
	fin >> mirX >> mirY >> mirZ;			fin.getline(dummy, NCHAR);
	fin >> gx >> gy >> gz;				fin.getline(dummy, NCHAR);
	fin >> pre_xp >> p_xp >> pre_xn >> p_xn;	fin.getline(dummy, NCHAR);
	fin >> pre_yp >> p_yp >> pre_yn >> p_yn;	fin.getline(dummy, NCHAR);
	fin >> pre_zp >> p_zp >> pre_zn >> p_zn;	fin.getline(dummy, NCHAR);
	fin >> vel_xp >> v_xp >> vel_xn >> v_xn;	fin.getline(dummy, NCHAR);
	fin >> vel_yp >> v_yp >> vel_yn >> v_yn;	fin.getline(dummy, NCHAR);
	fin >> vel_zp >> v_zp >> vel_zn >> v_zn;	fin.getline(dummy, NCHAR);
	fin >> in_vis;					fin.getline(dummy, NCHAR);
	fin >> inivx >> inivy >> inivz;			fin.getline(dummy, NCHAR);
							fin.getline(dummy, NCHAR);
	fin >> wr_per;					fin.getline(dummy, NCHAR);
	fin >> PerDir;					fin.getline(dummy, NCHAR);
	fin >> freRe;					fin.getline(dummy, NCHAR);
	fin >> freVe;					fin.getline(dummy, NCHAR);
	fin >> freDe;					fin.getline(dummy, NCHAR);
	fin >> mir;					fin.getline(dummy, NCHAR);
							fin.getline(dummy, NCHAR);
	fin >> Par_Geo >> Par_nx >> Par_ny >> Par_nz;	fin.getline(dummy, NCHAR);
	fin >> Zoom;					fin.getline(dummy, NCHAR);
	fin >> outputfile;				fin.getline(dummy, NCHAR);
	fin >> EI;					fin.getline(dummy, NCHAR);
	fin >> q_p;					fin.getline(dummy, NCHAR);
	fin.close();
	
	//cout<<NX<<"    asdfa "<<endl;
	NX=NX-1;NY=NY-1;NZ=NZ-1;
	}

	MPI_Bcast(&filename,128,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&mirX,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&NX,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&mirY,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&NY,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&mirZ,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&NZ,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&n_max,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&reso,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&in_IMR,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&gx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&gy,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&gz,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_xp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&p_xp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_xn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&p_xn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(&pre_yp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&p_yp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_yn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&p_yn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_zp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&p_zp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_zn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&p_zn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&vel_xp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&v_xp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&vel_xn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&v_xn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&vel_yp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&v_yp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&vel_yn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&v_yn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&vel_zp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&v_zp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&vel_zn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&v_zn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(&inivx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&inivy,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&inivz,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&wr_per,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&PerDir,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&freRe,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&freVe,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&freDe,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&mir,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&in_vis,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Zoom,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&outputfile,128,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&EI,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&q_p,1,MPI_DOUBLE,0,MPI_COMM_WORLD);


      
int U_max_ref=0;


if (mirX==1)
	NX=NX*2+1;
if (mirY==1)
	NY=NY*2+1;
if (mirZ==1)
	NZ=NZ*2+1;

/*
//*************THIN SOLID BOUNDARY MESH REFINEDMENT**************
if (Zoom>1)
	{	
	NX=(NX)*Zoom;
	NY=(NY)*Zoom;
	NZ=(NZ)*Zoom;
	}
//***************************************************************
*/

if (Zoom>1)
	{	
	NX=(NX+1)*Zoom-1;
	NY=(NY+1)*Zoom-1;
	NZ=(NZ+1)*Zoom-1;
	}

//if (Zoom>1)	
//	reso=reso/Zoom;



	nx_l=(int)((NX+1)/para_size);
	dif=(NX+1)-nx_l*para_size;
	
	if (rank>para_size-1-dif)
		nx_l+=1;

//	if (rank==para_size-1)
//		nx_l+=(NX+1)%para_size;

	double* Permia;
	double* rho;
	double** u;
	double**f;
	double**F;
	double**u0;
	int* SupInv;
	double* forcex;
	double* forcey;
	double* forcez;

	int*** Solids;
	int*** Solid;
	int*  Sl;
	int*  Sr;

		
	
	Solid = new int**[nx_l];
	Sl = new int[(NY+1)*(NZ+1)];
	Sr = new int[(NY+1)*(NZ+1)];

	Solids = new int**[(NX+1)/Zoom];

	for (int i=0;i<(NX+1)/Zoom;i++)
		{		
		Solids[i] = new int*[(NY+1)/Zoom];
	
			for (int j=0;j<(NY+1)/Zoom;j++)
			{
			Solids[i][j]= new int[(NZ+1)/Zoom];
			

			}
		}


	for (int i=0;i<nx_l;i++)
		{
		Solid[i] = new int*[NY+1];
			for (int j=0;j<=NY;j++)
			Solid[i][j]= new int[NZ+1];
		}

	

	
	
	
//	Count = new int[rank];
//	if (!(filename=="NOSOLID"))
		Read_Rock(Solids,&porosity,filename);
//	else
//		{
//		for(int i=0;i<=NX;i++)	
//			for(int j=0;j<=NY;j++)
//				for(int k=0;k<=NZ;k++)
//				Solids[i][j][k]=0;
//		}

	

	init_Sparse(Solids,Solid,Sl,Sr);

	for (int i=0;i<(NX+1)/Zoom;i++)
		{
		for (int j=0;j<(NY+1)/Zoom;j++)
			delete [] Solids[i][j];
		delete [] Solids[i];
		}
	delete [] Solids;


	//***************************************************
	//WARRING: SPARSE MATRIX STARTS FROM INDEX 1 NOT 0!!!
	//***************************************************

	Permia = new double[3];
	rho = new double[Count+1];
	forcex = new double[Count+1];
	forcey = new double[Count+1];
	forcez = new double[Count+1];
	u = new double*[Count+1];
	f = new double*[Count+1];
	F = new double*[Count+1];
	u0 = new double*[Count+1];
	SupInv = new int[Count+1];

	for (int i=0;i<=Count;i++)
		{
		u[i] = new double[3];
		f[i] = new double[19];
		u0[i] = new double[3];
		F[i] = new double[19];
		}

	Comput_MI(M,MI);
	
	Suppliment(SupInv,Solid);

	MPI_Barrier(MPI_COMM_WORLD);
	
	
	Geometry(Solid);

	init(rho,u,f,forcex,forcey,forcez);

if (rank==0)
		cout<<"Porosity= "<<porosity<<endl;

/*
char FileName[128]="Results.txt";
char FileName2[128]="Permeability.txt";
char FileName3[128]="bodyforce.txt";
*/




//========================================================
char FileName[128];strcpy(FileName,outputfile);
char FileName2[128];strcpy(FileName2,outputfile);
char FileName3[128];strcpy(FileName3,outputfile);

strcat(FileName,"Results.txt");
strcat(FileName2,"Permeability.txt");
strcat(FileName3,"bodyforce.txt");
//========================================================



ofstream fins;	
	fins.open(FileName,ios::trunc);
	fins.close();

if (wr_per==1)
	{
	fins.open(FileName2,ios::trunc);
	fins.close();
	}


	fins.open(FileName3,ios::out);
	fins.close();




	
	for(n=0;n<=n_max;n++)
	{
	
	//cout<<"@@@@@@@@@@@   "<<n<<endl;
	collision(rho,u,f,forcex,forcey,forcez,SupInv,Solid);//cout<<"markcollision"<<endl;

	
	if (EI==0)
		{periodic_streaming(f,F,SupInv,Solid,Sl,Sr,rho,u);}//cout<<"mark"<<endl;}
	else
		periodic_streaming_MR(f,F,SupInv,Solid,Sl,Sr,rho,u);
	
//	if ((1-pre_xp)*(1-pre_xn)*(1-pre_yp)*(1-pre_yn)*(1-pre_zp)*(1-pre_zn)==0)
//		boundary_pressure(pre_xp,p_xp,pre_xn,p_xn,pre_yp,p_yp,pre_yn,p_yn,pre_zp,p_zp,pre_zn,p_zn,f,u,rho,Solid);
	
//	if ((1-vel_xp)*(1-vel_xn)*(1-vel_yp)*(1-vel_yn)*(1-vel_zp)*(1-vel_zn)==0)
//		boundary_velocity(vel_xp,v_xp,vel_xn,v_xn,vel_yp,v_yp,vel_yn,v_yn,vel_zp,v_zp,vel_zn,v_zn,f,rho,u,Solid);

	if (in_IMR==1)
		comput_macro_variables_IMR(rho,u,u0,f,F,SupInv,Solid,forcex,forcey,forcez,n,porosity,&gx,&gy,&gz);
	else
  		comput_macro_variables(rho,u,u0,f,F,forcex,forcey,forcez,SupInv,Solid); 

	if ((1-pre_xp)*(1-pre_xn)*(1-pre_yp)*(1-pre_yn)*(1-pre_zp)*(1-pre_zn)==0)
		boundary_pressure(pre_xp,p_xp,pre_xn,p_xn,pre_yp,p_yp,pre_yn,p_yn,pre_zp,p_zp,pre_zn,p_zn,f,u,rho,Solid);
	
	if ((1-vel_xp)*(1-vel_xn)*(1-vel_yp)*(1-vel_yn)*(1-vel_zp)*(1-vel_zn)==0)
		boundary_velocity(vel_xp,v_xp,vel_xn,v_xn,vel_yp,v_yp,vel_yn,v_yn,vel_zp,v_zp,vel_zn,v_zn,f,rho,u,Solid);
	
	if(n%freRe==0)
		{       
			
			error=Error(u,u0,&u_max,&u_ave);if (u_max>=10.0)	U_max_ref+=1;
			error_perm=Comput_Perm(u,Permia,PerDir);
			if (rank==0)
			{
			finish = MPI_Wtime();
			ofstream fin(FileName,ios::app);
			fin<<"The"<<n<<"th computation result:"<<endl;
		//=============================================================================================
			fin<<"The permiability is: "<<Permia[0]*reso*reso*1000<<", "<<Permia[1]*reso*reso*1000<<", "<<Permia[2]*reso*reso*1000<<endl;
			fin<<"The relative error of permiability computing is: "<<error_perm<<endl;
		//==============================================================================================

		//==============================================================================================
			Re=u_ave*(NY+1)/(1.0/3.0*(1/s_v-0.5));
			fin<<"The Maximum velocity is: "<<setprecision(6)<<u_max<<"   Re="<<Re<<endl;
			
		//===============================================================================================
			fin<<"The max relative error of velocity is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			fin<<"Elapsed time is "<< finish-start <<" seconds"<<endl;
			fin<<endl;
			fin.close();


		if (wr_per==1)
			{
			ofstream finfs(FileName2,ios::app);
			switch(PerDir)
				{
				case 1:
				finfs<<Permia[0]*reso*reso*1000<<endl;break;
				case 2:
				finfs<<Permia[1]*reso*reso*1000<<endl;break;
				case 3:
				finfs<<Permia[2]*reso*reso*1000<<endl;break;
				default:
				finfs<<Permia[0]*reso*reso*1000<<endl;break;
				}
			finfs.close();
			}

		
			ofstream finf3(FileName3,ios::app);
			finf3<<gx<<endl;
			finf3.close();
		

			cout<<"The"<<n<<"th computation result:"<<endl;
		//=============================================================================================
			cout<<"The permiability is: "<<Permia[0]*reso*reso*1000<<", "<<Permia[1]*reso*reso*1000<<", "<<Permia[2]*reso*reso*1000<<endl;
			cout<<"The relative error of permiability computing is: "<<error_perm<<endl;
		//==============================================================================================

		//==============================================================================================
			
			Re=u_ave*(NY+1)/(1.0/3.0*(1/s_v-0.5));
			cout<<"The Maximum velocity is: "<<setprecision(6)<<u_max<<"   Re="<<Re<<endl;
			
		//===============================================================================================
			cout<<"The max relative error of uv is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			cout<<"Elapsed time is "<< finish-start <<" seconds"<<endl;
			cout<<endl;
			}
			
			if ((freDe>=0) and (n%freDe==0))
			output_density(n,rho,mirX,mirY,mirZ,mir,Solid);
			if ((freVe>=0) and (n%freVe==0))
			output_velocity(n,rho,u,mirX,mirY,mirZ,mir,Solid);

			if(error!=error) {cout<<"PROGRAM STOP"<<endl;break;};
			if(U_max_ref>=5) {cout<<"PROGRAM STOP DUE TO HIGH VELOCITY"<<endl;break;}
		}	
	}


	

	for (int i=0;i<=Count;i++)
		{
		delete [] u[i];
		delete [] u0[i];
		delete [] f[i];
		delete [] F[i];
		}
	
	delete [] f;
	delete [] u;
	delete [] F;
	delete [] u0;
	delete [] rho;
	delete [] forcex;
	delete [] forcey;
	delete [] forcez;
	delete [] SupInv;

	delete [] Sl;
	delete [] Sr;

	delete [] Permia;
//	delete [] Count;
	finish = MPI_Wtime();

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0)
			{
    			cout<<"Elapsed time is "<< finish-start <<" seconds"<<endl;
    			cout<<"Accuracy: "<<MPI_Wtick()<<" Second"<<endl;
			ofstream fin(FileName,ios::app);
			fin<<"Elapsed time is "<< finish-start <<" seconds"<<endl;
			fin<<"Accuracy: "<<MPI_Wtick()<<" Second"<<endl;
			}

	MPI :: Finalize ();

	
}


int inverse(mat &a){
    double t;
    int i,j,k,is[MAXN],js[MAXN];
    if(a.n!=a.m) return 0;
    for(k=0;k<a.n;k++){
        for(t=0,i=k;i<a.n;i++)
            for(j=k;j<a.n;j++)
                if(fabs(a.data[i][j])>t)
                    t=fabs(a.data[is[k]=i][js[k]=j]);
        if(zero(t)) return 0;
        if(is[k]!=k)
            for(j=0;j<a.n;j++)
                t=a.data[k][j],a.data[k][j]=a.data[is[k]][j],a.data[is[k]][j]=t;
        if(js[k]!=k)
            for(i=0;i<a.n;i++)
                t=a.data[i][k],a.data[i][k]=a.data[i][js[k]],a.data[i][js[k]]=t;
        a.data[k][k]=1/a.data[k][k];
        for(j=0;j<a.n;j++)
            if(j!=k)
                a.data[k][j]*=a.data[k][k];
        for(i=0;i<a.n;i++)
            if(i!=k)
                for(j=0;j<a.n;j++)
                    if(j!=k)
                        a.data[i][j]-=a.data[i][k]*a.data[k][j];
        for(i=0;i<a.n;i++)
            if(i!=k)
                a.data[i][k]*=-a.data[k][k];
    }
    for(k=a.n-1;k>=0;k--){
        for(j=0;j<a.n;j++)
            if(js[k]!=k)
                t=a.data[k][j],a.data[k][j]=a.data[js[k]][j],a.data[js[k]][j]=t;
        for(i=0;i<a.n;i++)
            if(is[k]!=k)
                t=a.data[i][k],a.data[i][k]=a.data[i][is[k]],a.data[i][is[k]]=t;
    }
    return 1;
}


void Comput_MI(double M[19][19], double MI[19][19])
{

double mim[19][19];

mat a;
    int i,j;
    int n_s=19;
        for(int i=0;i<n_s;i++)
            for(int j=0;j<n_s;j++)
                a.data[i][j]=M[i][j];
	a.m=a.n=n_s;
        if(inverse(a))
            for(int i=0;i<n_s;i++)
                for(int j=0;j<n_s;j++)
                    MI[i][j]=a.data[i][j];
               
            
        else
            puts("NO");


}



void init_Sparse(int*** Solids, int*** Solid,int* Sl,int* Sr)
{	
	MPI_Status status[4] ;
	MPI_Request request[4];

	
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

bool mark;
int kk,ip,jp,kp,mean_l,mpi_test,s_c;
	
	mean_l=(int)((NX+1)/mpi_size);

	s_c=0;
	for (int i=0;i<rank;i++)
		if (i<mpi_size-((NX+1)-mean_l*mpi_size))
			s_c+=mean_l;
		else
			s_c+=mean_l+1;

		


	int* Sl_send;
	int* Sr_send;
	Sl_send = new int[(NY+1)*(NZ+1)];
	Sr_send = new int[(NY+1)*(NZ+1)];

	Count=1;

//cout<<Solids[3][4][5]<<"      @@@@@@@@@@@@@@@@@@"<<endl;


//for(int i=0;i<nx_l;i++)	
//	for(int j=0;j<=NY;j++)
//		for(int k=0;k<=NZ;k++)
//		cout<<Solids[int((s_c+i-(s_c+i)%Zoom)/Zoom)][int((j-j%Zoom)/Zoom)][int((k-k%Zoom)/Zoom)]<<endl;


for(int i=0;i<nx_l;i++)	
	for(int j=0;j<=NY;j++)
		for(int k=0;k<=NZ;k++)
		{
			//cout<<(s_c+i-(s_c+i)%Zoom)/Zoom<<"   "<<(j-j%Zoom)/Zoom<<"  "<<(k-k%Zoom)/Zoom<<endl;

			if (Solids[int((s_c+i-(s_c+i)%Zoom)/Zoom)][int((j-j%Zoom)/Zoom)][int((k-k%Zoom)/Zoom)]==0)
				{
			
				Solid[i][j][k]=Count;
				Count++;
				}
			else
				Solid[i][j][k]=0;

		}

	Count-=1;
	
	cl=0;cr=0;	
	
	for (int j=0;j<=NY;j++)
		for (int k=0;k<=NZ;k++)
		{
			
			if (Solid[0][j][k]>0)
				{
				cl++;
				Sl_send[j*(NZ+1)+k]=cl;
				}
			else
				Sl_send[j*(NZ+1)+k]=0;
			

		
			if (Solid[nx_l-1][j][k]>0)
				{
				cr++;
				Sr_send[j*(NZ+1)+k]=cr;
				}
			else
				Sr_send[j*(NZ+1)+k]=0;
			
		}




if (rank==0)
		{
		MPI_Isend(Sr_send, (NY+1)*(NZ+1), MPI_INT, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(Sl_send, (NY+1)*(NZ+1), MPI_INT, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(Sr , (NY+1)*(NZ+1), MPI_INT, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(Sl, (NY+1)*(NZ+1), MPI_INT, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );
		
		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(Sl_send, (NY+1)*(NZ+1), MPI_INT, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(Sr_send, (NY+1)*(NZ+1), MPI_INT, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(Sl, (NY+1)*(NZ+1), MPI_INT, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(Sr, (NY+1)*(NZ+1), MPI_INT, 0, 0, MPI_COMM_WORLD,&request[3]);
			
			}
			else
			{
			MPI_Isend(Sl_send, (NY+1)*(NZ+1), MPI_INT, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(Sr_send, (NY+1)*(NZ+1), MPI_INT, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(Sl, (NY+1)*(NZ+1), MPI_INT, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(Sr, (NY+1)*(NZ+1), MPI_INT, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);
			
			};

	
	MPI_Waitall(4,request, status);
	MPI_Testall(4,request,&mpi_test,status);

	delete [] Sl_send;
	delete [] Sr_send;	
}

void Suppliment(int* SupInv,int*** Solid)
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	for(int i=0;i<nx_l;i++)	
		for(int j=0;j<=NY;j++)
			for(int k=0;k<=NZ;k++)	
			if (Solid[i][j][k]>0)
				SupInv[abs(Solid[i][j][k])]=i*(NY+1)*(NZ+1)+j*(NZ+1)+k;


}

void Read_Rock(int*** Solids,double* porosity,char poreFileName[128])
{


int rank = MPI :: COMM_WORLD . Get_rank ();
int mpi_size=MPI :: COMM_WORLD . Get_size ();



int nx0=NX+1;
int ny0=NY+1;
int nz0=NZ+1;


int nx=NX+1;
int ny=NY+1;
int nz=NZ+1;

int* Solid_Int;
int nx_a,ny_a,nz_a;

if (Zoom>1)
	{
	nx0=(nx0)/Zoom;
	ny0=(ny0)/Zoom;
	nz0=(nz0)/Zoom;

	nx=(nx)/Zoom;
	ny=(ny)/Zoom;
	nz=(nz)/Zoom;	
	}



if (mirX==1)
	nx0=(nx0)/2;
if (mirY==1)
	ny0=(ny0)/2;
if (mirZ==1)
	nz0=(nz0)/2;


double pore;
int i, j, k,ir,jr,kr;

Solid_Int = new int[nx*ny*nz];



	

if (rank==0)
{


if (Par_Geo==0)
	{
	nx_a=nx0;
	ny_a=ny0;
	nz_a=nz0;
	}
else
	{
	nx_a=Par_nx;
	ny_a=Par_ny;
	nz_a=Par_nz;
	}

FILE *ftest;
	ifstream fin;
	
	ftest = fopen(poreFileName, "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file (" << poreFileName <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest);

	fin.open(poreFileName);


	
	// Reading pore geometry
	for(k=0 ; k<nz_a ; k++)
	for(j=0 ; j<ny_a ; j++)
	for(i=0 ; i<nx_a ; i++)
	
	{
		while(true)
		{	
			fin >> pore;
			if( pore == 0.0 || pore == 1.0) break;
		}
		if ((pore == 0.0) && (i<nx0) && (j<ny0) && (k<nz0))	Solid_Int[i*ny*nz+j*nz+k] = 0;
		//else			Solid_Int[i][j][k] = 1;
		if ((pore == 1.0) && (i<nx0) && (j<ny0) && (k<nz0))	Solid_Int[i*ny*nz+j*nz+k] = 1;
	}
	fin.close();

	// Mirroring the rock
	if(mirX==1){
		for(i=nx0 ; i<nx ; i++)
		for(j=0   ; j<ny ; j++)
		for(k=0   ; k<nz ; k++)
				Solid_Int[i*ny*nz+j*nz+k] = Solid_Int[(nx-i-1)*ny*nz+j*nz+k];
	}

	if(mirY==1){
		for(j=ny0 ; j<ny ; j++)
		for(i=0   ; i<nx ; i++)
		for(k=0   ; k<nz ; k++)
				Solid_Int[i*ny*nz+j*nz+k] = Solid_Int[i*ny*nz+(ny-j-1)*nz+k];
	}

	if(mirZ==1){
		for(k=nz0 ; k<nz ; k++)
		for(j=0   ; j<ny ; j++)
		for(i=0   ; i<nx ; i++)
				Solid_Int[i*ny*nz+j*nz+k] = Solid_Int[i*ny*nz+j*nz+nz-k-1];
	}


	//MESH REFINEMENT
	


	//double porosity;

	// Calculate Porosity
	int nNodes = 0;
	for(i=0 ; i<nx*ny*nz ; i++)
		if(Solid_Int[i] == 0) nNodes++;

	*porosity = (double)nNodes / (nx*ny*nz);

}


	
	MPI_Bcast(Solid_Int,nx*ny*nz,MPI_INT,0,MPI_COMM_WORLD);

	
	cout<<"INPUT FILE READING COMPLETE.  THE POROSITY IS: "<<*porosity<<endl;

	//cout<<nx<<"  "<<ny<<"  "<<nz<<"  zoom "<<Zoom<<endl;

	
	for (i=0;i<nx;i++)
		for (j=0;j<ny;j++)
			for (k=0;k<nz;k++)
			Solids[i][j][k]=Solid_Int[i*(ny)*(nz)+j*(nz)+k];

//cout<<"asdfasdfasdfasdfa"<<endl;
MPI_Barrier(MPI_COMM_WORLD);



	
/*
//*************THIN SOLID BOUNDARY MESH REFINEDMENT**************
	
	for (i=0;i<nx;i++)
		for (j=0;j<ny;j++)
			for (k=0;k<nz;k++)
			Solids[i*2][j*2][k*2]=Solid_Int[i*(ny)*(nz)+j*(nz)+k];



	cout<<nx<<" "<<ny<<" "<<nz<<endl;

	for (k=0;k<nz;k++)   //k==2*k
		{
		for (j=0;j<ny;j++)	//j=2*j
			for(i=0;i<nx-1;i++)
				{
				ir=2*i+1;jr=2*j;kr=2*k;
				if ((Solids[ir-1][jr][kr]==1) and (Solids[ir+1][jr][kr]==1))
					Solids[ir][jr][kr]=1;
				else
					Solids[ir][jr][kr]=0;
				}

		for (i=0;i<nx;i++)
			for (j=0;j<ny-1;j++)
				{
				ir=2*i;jr=2*j+1;kr=k*2;
				if ((Solids[ir][jr-1][kr]==1) and (Solids[ir][jr+1][kr]==1))
					Solids[ir][jr][kr]=1;
				else
					Solids[ir][jr][kr]=0;
				}

		}
	//cout<<"@@@@@@@@@@@@@"<<endl;
	for (k=0;k<nz-1;k++)
		{
		for (i=0;i<=NX;i++)
			for (j=0;j<=NY;j++)
			{
			kr=2*k+1;
			if ((Solids[i][j][kr-1]==1) and (Solids[i][j][kr+1]==1))
				Solids[i][j][kr]=1;
			else
					Solids[i][j][kr]=0;
			}
		}
//******************************************************************************
*/
/*
	for (i=0;i<nx;i++)
		for (j=0;j<(ny);j++)
			for (k=0;k<(nz);k++)
				for (int iz=0;iz<=Zoom-1;iz++)
					for (int jz=0;jz<=Zoom-1;jz++)
						for (int kz=0;kz<=Zoom-1;kz++)
						Solids[Zoom*i+iz][Zoom*j+jz][Zoom*k+kz]=Solid_Int[i*(ny)*(nz)+j*(nz)+k];
				


	}
	
*/
	delete [] Solid_Int;


}




void init(double* rho, double** u, double** f,double* forcex,double* forcey, double* forcez)
{	
      


	double Cylinder_r=7;
	double usqr,vsqr;

	
	rho0=1.0;dt=1.0/Zoom;dx=1.0/Zoom;
 	uMax=0.0;
	Re=50.0;
	//forcex=gx;forcey=gy;
	//niu=U*Lx/Re;
	niu=uMax*Cylinder_r*2/Re;
	niu=in_vis;
	tau_f=3.0*niu/dt+0.5;
	//gx=3*uMax*niu/((NY/2)*(NY/2));
	//tau_f=1.0;
	s_v=1/tau_f;
        //tau_f=1.0;
	//cout<<"tau_f= "<<tau_f<<endl;
	double pr; //raduis of the obstacles
        //if (Zoom>1)
	//s_v=1.0/((1.0/s_v-0.5)/Zoom+0.5);

	//srand(time(0));                        //RANDOM NUMBER GENERATION SEED
	//cylinder_creation(50,103,Cylinder_r);
	double s_other=8*(2-s_v)/(8-s_v);
	double u_tmp[3];
/*
	S[0]=s_v;
	S[1]=s_v;
	S[2]=s_v;
	S[3]=s_v;
	S[4]=s_v;
	S[5]=s_v;
	S[6]=s_v;
	S[7]=s_v;
	S[8]=s_v;
	S[9]=s_v;
	S[10]=s_v;
	S[11]=s_v;
	S[12]=s_v;
	S[13]=s_v;
	S[14]=s_v;
	S[15]=s_v;

	S[16]=1;
	S[17]=1;
	S[18]=1;

*/
	S[0]=0;
	S[1]=s_v;
	S[2]=s_v;
	S[3]=0;
	S[4]=s_other;
	S[5]=0;
	S[6]=s_other;
	S[7]=0;
	S[8]=s_other;
	S[9]=s_v;
	S[10]=s_v;
	S[11]=s_v;
	S[12]=s_v;
	S[13]=s_v;
	S[14]=s_v;
	S[15]=s_v;
	S[16]=s_other;
	S[17]=s_other;
	S[18]=s_other;




	// for(int i=0;i<=NX;i++)	
	//	for(int j=0;j<=NY;j++)
	//		for(int k=0;k<=NZ;k++)
	//
	//		Solid[i][j][k]=0;
	
	for (int i=1;i<=Count;i++)	
			
		{
			u[i][0]=inivx;
			u[i][1]=inivy;
			u[i][2]=inivz;
			u_tmp[0]=u[i][0];
			u_tmp[1]=u[i][1];
			u_tmp[2]=u[i][2];

			rho[i]=1.0;
			
			
			

			//***********************************************************************

			forcex[i]=gx;
			forcey[i]=gy;
			forcez[i]=gz;







			//***********************************************************************


			//INITIALIZATION OF m and f

			for (int lm=0;lm<19;lm++)
				//if (Solid[(int)(SupInv[i]/((NY+1)*(NZ+1)))][(int)((SupInv[i]%((NY+1)*(NZ+1)))/(NZ+1))][SupInv[i]%(NZ+1)]<0)
				//	f[i][lm]=0.0;
				//else
					f[i][lm]=feq(lm,rho[i],u_tmp);
				

		
		

	}

	

	 	
}


double feq(int k,double rho, double u[3])
{
	double eu,uv,feq;
        double c2,c4;
	double c=1;
	c2=c*c;c4=c2*c2;
	eu=(e[k][0]*u[0]+e[k][1]*u[1]+e[k][2]*u[2]);
	uv=(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);// WITH FORCE TERM:GRAVITY IN X DIRECTION
	feq=w[k]*rho*(1.0+3.0*eu/c2+4.5*eu*eu/c4-1.5*uv/c2);
	return feq;
}



void periodic_streaming_MR(double** f,double** F,int* SupInv,int*** Solid,int* Sl,int* Sr,double* rho,double** u)
{
	MPI_Status status[4] ;
	MPI_Request request[4];

	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	int ip,jp,kp,in,jn,kn,i,j,k,mpi_test;
	double rho_ls,v_ls[3],v_ls2[3];
	double qprim=(1-2*q_p)*1.2;

	int* Gcl = new int[mpi_size];
	int* Gcr = new int[mpi_size];

	
	MPI_Gather(&cl,1,MPI_INT,Gcl,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(Gcl,mpi_size,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Gather(&cr,1,MPI_INT,Gcr,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(Gcr,mpi_size,MPI_INT,0,MPI_COMM_WORLD);


	double* recvl,*recvl_rho,*recvl_u,*sendl_rho,*sendl_u;
	double* recvr,*recvr_rho,*recvr_u,*sendr_rho,*sendr_u;

	double* sendl = new double[Gcl[rank]*19];
	double* sendr = new double[Gcr[rank]*19];


		for(k=0;k<19;k++)
			{
			for(i=1;i<=Gcl[rank];i++)
				sendl[(i-1)*19+k]=f[i][k];
			for(j=Count-Gcr[rank]+1;j<=Count;j++)
				sendr[(j-(Count-Gcr[rank]+1))*19+k]=f[j][k];
			}


	
	if (rank==0)
		{
		recvl = new double[Gcr[mpi_size-1]*19];
		recvr = new double[Gcl[rank+1]*19];
		}	
		else
		if (rank==mpi_size-1)
			{
			recvl = new double[Gcr[rank-1]*19];
			recvr = new double[Gcl[0]*19];
			}
			else
			{
			recvl = new double[Gcr[rank-1]*19];

			recvr = new double[Gcl[rank+1]*19];
			}

if (q_p<0.5)
	{

	sendl_rho = new double[Gcl[rank]];
	sendr_rho = new double[Gcr[rank]];
		
	sendl_u = new double[Gcl[rank]*3];
	sendr_u = new double[Gcr[rank]*3];


	for(k=0;k<3;k++)
			{
			for(i=1;i<=Gcl[rank];i++)
				sendl_u[(i-1)*3+k]=u[i][k];
			for(j=Count-Gcr[rank]+1;j<=Count;j++)
				sendr_u[(j-(Count-Gcr[rank]+1))*3+k]=u[j][k];
			}

	for(i=1;i<=Gcl[rank];i++)
		sendl_rho[(i-1)]=rho[i];
	for(j=Count-Gcr[rank]+1;j<=Count;j++)
		sendr_rho[(j-(Count-Gcr[rank]+1))]=rho[j];


	if (rank==0)
		{
		recvl_rho = new double[Gcr[mpi_size-1]];
		recvr_rho = new double[Gcl[rank+1]];
		recvl_u = new double[Gcr[mpi_size-1]*3];
		recvr_u = new double[Gcl[rank+1]*3];
		}	
		else
		if (rank==mpi_size-1)
			{
			recvl_rho = new double[Gcr[rank-1]];
			recvr_rho = new double[Gcl[0]];
			recvl_u = new double[Gcr[rank-1]*3];
			recvr_u = new double[Gcl[0]*3];
			}
			else
			{
			recvl_rho = new double[Gcr[rank-1]];
			recvr_rho = new double[Gcl[rank+1]];
			recvl_u = new double[Gcr[rank-1]*3];
			recvr_u = new double[Gcl[rank+1]*3];
			}

	}

MPI_Barrier(MPI_COMM_WORLD);


if (rank==0)
		{
		
		MPI_Isend(sendr, Gcr[0]*19, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(sendl, Gcl[0]*19, MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(recvr , Gcl[1]*19, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(recvl, Gcr[mpi_size-1]*19, MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );
		
		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(sendl, Gcl[rank]*19, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr, Gcr[rank]*19, MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl, Gcr[rank-1]*19, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(recvr, Gcl[0]*19, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request[3]);
			
			}
			else
			{

			MPI_Isend(sendl, Gcl[rank]*19, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr, Gcr[rank]*19, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl, Gcr[rank-1]*19, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(recvr, Gcl[rank+1]*19, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);

			
			};

	
	MPI_Waitall(4,request, status);

	MPI_Testall(4,request,&mpi_test,status);

	delete [] sendl;
	delete [] sendr;


	if (q_p<0.5)
	{
		
	if (rank==0)
		{
		
		MPI_Isend(sendr_rho, Gcr[0], MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(sendl_rho, Gcl[0], MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(recvr_rho , Gcl[1], MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(recvl_rho, Gcr[mpi_size-1], MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );
		
		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(sendl_rho, Gcl[rank], MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr_rho, Gcr[rank], MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl_rho, Gcr[rank-1], MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(recvr_rho, Gcl[0], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request[3]);
			
			}
			else
			{

			MPI_Isend(sendl_rho, Gcl[rank], MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr_rho, Gcr[rank], MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl_rho, Gcr[rank-1], MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(recvr_rho, Gcl[rank+1], MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);

			
			};

	
	MPI_Waitall(4,request, status);

	MPI_Testall(4,request,&mpi_test,status);

	delete [] sendl_rho;
	delete [] sendr_rho;

	if (rank==0)
		{
		
		MPI_Isend(sendr_u, Gcr[0]*3, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(sendl_u, Gcl[0]*3, MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(recvr_u , Gcl[1]*3, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(recvl_u, Gcr[mpi_size-1]*3, MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );
		
		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(sendl_u, Gcl[rank]*3, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr_u, Gcr[rank]*3, MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl_u, Gcr[rank-1]*3, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(recvr_u, Gcl[0]*3, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request[3]);
			
			}
			else
			{

			MPI_Isend(sendl_u, Gcl[rank]*3, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr_u, Gcr[rank]*3, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl_u, Gcr[rank-1]*3, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(recvr_u, Gcl[rank+1]*3, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);

			
			};

	
	MPI_Waitall(4,request, status);

	MPI_Testall(4,request,&mpi_test,status);

	delete [] sendl_u;
	delete [] sendr_u;

	}



if (q_p>=0.5)


	{
		for(int ci=1;ci<=Count;ci++)	
		{
				i=(int)(SupInv[ci]/((NY+1)*(NZ+1)));
				j=(int)((SupInv[ci]%((NY+1)*(NZ+1)))/(NZ+1));
				k=(int)(SupInv[ci]%(NZ+1));			
			for(int lm=0;lm<19;lm++)
			{       

				ip=i-e[lm][0];
				jp=j-e[lm][1];if (jp<0) {jp=NY;}; if (jp>NY) {jp=0;};
				kp=k-e[lm][2];if (kp<0) {kp=NZ;}; if (kp>NZ) {kp=0;};

			
				if (ip<0)
					{
					if (Sl[jp*(NZ+1)+kp]>0)
						F[ci][lm]=recvl[(Sl[jp*(NZ+1)+kp]-1)*19+lm];
					else
						{
						v_ls[0]=u[ci][0];
						v_ls[1]=u[ci][1];
						v_ls[2]=u[ci][2];
						F[ci][lm]=f[ci][LR[lm]]-feq(LR[lm],rho[ci],v_ls)+(2*q_p-1)/q_p*w[LR[lm]]+(1-q_p)/q_p*feq(LR[lm],rho[ci],v_ls);
						//F[ci][lm]=f[ci][LR[lm]];
						//cout<<i<<" "<<j<<" "<<k<<endl;
						}
					}

				if (ip>=nx_l)
					{
					if (Sr[jp*(NZ+1)+kp]>0)
						F[ci][lm]=recvr[(Sr[jp*(NZ+1)+kp]-1)*19+lm];
					else
						{
						v_ls[0]=u[ci][0];
						v_ls[1]=u[ci][1];
						v_ls[2]=u[ci][2];
						F[ci][lm]=f[ci][LR[lm]]-feq(LR[lm],rho[ci],v_ls)+(2*q_p-1)/q_p*w[LR[lm]]+(1-q_p)/q_p*feq(LR[lm],rho[ci],v_ls);
						//F[ci][lm]=f[ci][LR[lm]];
						}
					}

				if ((ip>=0) and (ip<nx_l))
					{
					if (Solid[ip][jp][kp]>0)
						F[ci][lm]=f[abs(Solid[ip][jp][kp])][lm];
					else
						{
						v_ls[0]=u[ci][0];
						v_ls[1]=u[ci][1];
						v_ls[2]=u[ci][2];
						F[ci][lm]=f[ci][LR[lm]]-feq(LR[lm],rho[ci],v_ls)+(2*q_p-1)/q_p*w[LR[lm]]+(1-q_p)/q_p*feq(LR[lm],rho[ci],v_ls);
						//F[ci][lm]=f[ci][LR[lm]];
						}
					}
				
			}	
				
		}

	}
else

	
for(int ci=1;ci<=Count;ci++)	
{
	i=(int)(SupInv[ci]/((NY+1)*(NZ+1)));
	j=(int)((SupInv[ci]%((NY+1)*(NZ+1)))/(NZ+1));
	k=(int)(SupInv[ci]%(NZ+1));			
	
	if ((i>=1) and (i<nx_l-1))
	{

		for(int lm=0;lm<19;lm++)
		{       
		ip=i-e[lm][0];
		jp=j-e[lm][1];if (jp<0) {jp=NY;}; if (jp>NY) {jp=0;};
		kp=k-e[lm][2];if (kp<0) {kp=NZ;}; if (kp>NZ) {kp=0;};
	
		in=i+e[lm][0];
		jn=j+e[lm][1];if (jn<0) {jn=NY;}; if (jn>NY) {jn=0;};
		kn=k+e[lm][2];if (kn<0) {kn=NZ;}; if (kn>NZ) {kn=0;};

		if (Solid[ip][jp][kp]>0)
			F[ci][lm]=f[abs(Solid[ip][jp][kp])][lm];
		else
			{
			v_ls[0]=u[ci][0];
			v_ls[1]=u[ci][1];
			v_ls[2]=u[ci][2];
			if (Solid[in][jn][kn]>0)
				F[ci][lm]=f[ci][LR[lm]]-feq(LR[lm],rho[ci],v_ls)+(1-2*q_p)*feq(LR[lm],rho[Solid[in][jn][kn]],u[Solid[in][jn][kn]])+2*q_p*feq(LR[lm],rho[ci],v_ls);
			else
				F[ci][lm]=(qprim+2*q_p-1)/qprim*(feq(LR[lm],rho[ci],v_ls))+(1-2*q_p)/qprim*w[LR[lm]];
						
						
			} 
		}
	}
	else
		for(int lm=0;lm<19;lm++)
		{       
		ip=i-e[lm][0];
		jp=j-e[lm][1];if (jp<0) {jp=NY;}; if (jp>NY) {jp=0;};
		kp=k-e[lm][2];if (kp<0) {kp=NZ;}; if (kp>NZ) {kp=0;};

		in=i+e[lm][0];
		jn=j+e[lm][1];if (jn<0) {jn=NY;}; if (jn>NY) {jn=0;};
		kn=k+e[lm][2];if (kn<0) {kn=NZ;}; if (kn>NZ) {kn=0;};

		if (i==0)
		{
			if (ip<0)
			{
			if (Sl[jp*(NZ+1)+kp]>0)
				F[ci][lm]=recvl[(Sl[jp*(NZ+1)+kp]-1)*19+lm];
			else
				{
				v_ls[0]=u[ci][0];v_ls2[0]=u[Solid[in][jn][kn]][0];
				v_ls[1]=u[ci][1];v_ls2[1]=u[Solid[in][jn][kn]][1];
				v_ls[2]=u[ci][2];v_ls2[2]=u[Solid[in][jn][kn]][2];
		
				if (Solid[in][jn][kn]>0)
					F[ci][lm]=f[ci][LR[lm]]-feq(LR[lm],rho[ci],v_ls)+(1-2*q_p)*feq(LR[lm],rho[Solid[in][jn][kn]],v_ls2)+2*q_p*feq(LR[lm],rho[ci],v_ls);
				else
					F[ci][lm]=(qprim+2*q_p-1)/qprim*(feq(LR[lm],rho[ci],v_ls))+(1-2*q_p)/qprim*w[LR[lm]];
								
						
				}
			}

			if (in<0)
			{
			if (Solid[ip][jp][kp]>0)
				F[ci][lm]=f[abs(Solid[ip][jp][kp])][lm];
			else
				{
				v_ls[0]=u[ci][0];v_ls2[0]=recvl_u[(Sl[jn*(NZ+1)+kn]-1)*3];
				v_ls[1]=u[ci][1];v_ls2[1]=recvl_u[(Sl[jn*(NZ+1)+kn]-1)*3+1];
				v_ls[2]=u[ci][2];v_ls2[2]=recvl_u[(Sl[jn*(NZ+1)+kn]-1)*3+2];
	
				if (Sl[jn*(NZ+1)+kn]>0)
				F[ci][lm]=f[ci][LR[lm]]-feq(LR[lm],rho[ci],v_ls)+(1-2*q_p)*feq(LR[lm],recvl_rho[Sl[jn*(NZ+1)+kn]-1],v_ls2)+2*q_p*feq(LR[lm],rho[ci],v_ls);
				else
				F[ci][lm]=(qprim+2*q_p-1)/qprim*(feq(LR[lm],rho[ci],v_ls))+(1-2*q_p)/qprim*w[LR[lm]];
						

				}
			}
					
			if ((ip>=0) and (in>=0))
			{
					
			if (Solid[ip][jp][kp]>0)
				F[ci][lm]=f[abs(Solid[ip][jp][kp])][lm];
			else
				{
				v_ls[0]=u[ci][0];v_ls2[0]=u[Solid[in][jn][kn]][0];
				v_ls[1]=u[ci][1];v_ls2[1]=u[Solid[in][jn][kn]][1];
				v_ls[2]=u[ci][2];v_ls2[2]=u[Solid[in][jn][kn]][2];
	
				if (Solid[in][jn][kn]>0)
					F[ci][lm]=f[ci][LR[lm]]-feq(LR[lm],rho[ci],v_ls)+(1-2*q_p)*feq(LR[lm],rho[Solid[in][jn][kn]],v_ls2)+2*q_p*feq(LR[lm],rho[ci],v_ls);
				else
					F[ci][lm]=(qprim+2*q_p-1)/qprim*(feq(LR[lm],rho[ci],v_ls))+(1-2*q_p)/qprim*w[LR[lm]];
						
						
				} 	


			}
			


		}



		if (i==nx_l-1)
		{
			if (ip>=nx_l)
			{
				if (Sr[jp*(NZ+1)+kp]>0)
					F[ci][lm]=recvr[(Sr[jp*(NZ+1)+kp]-1)*19+lm];
				else
					{
					v_ls[0]=u[ci][0];v_ls2[0]=u[Solid[in][jn][kn]][0];
					v_ls[1]=u[ci][1];v_ls2[1]=u[Solid[in][jn][kn]][1];
					v_ls[2]=u[ci][2];v_ls2[2]=u[Solid[in][jn][kn]][2];
		
					if (Solid[in][jn][kn]>0)
						F[ci][lm]=f[ci][LR[lm]]-feq(LR[lm],rho[ci],v_ls)+(1-2*q_p)*feq(LR[lm],rho[Solid[in][jn][kn]],v_ls2)+2*q_p*feq(LR[lm],rho[ci],v_ls);
					else
						F[ci][lm]=(qprim+2*q_p-1)/qprim*(feq(LR[lm],rho[ci],v_ls))+(1-2*q_p)/qprim*w[LR[lm]];
								
					
					}
			}

			if (in>=nx_l)
			{
				if (Solid[ip][jp][kp]>0)
					F[ci][lm]=f[abs(Solid[ip][jp][kp])][lm];
				else
					{
					v_ls[0]=u[ci][0];v_ls2[0]=recvr_u[(Sr[jn*(NZ+1)+kn]-1)*3];
					v_ls[1]=u[ci][1];v_ls2[1]=recvr_u[(Sr[jn*(NZ+1)+kn]-1)*3+1];
					v_ls[2]=u[ci][2];v_ls2[2]=recvr_u[(Sr[jn*(NZ+1)+kn]-1)*3+2];
					if (Sr[jn*(NZ+1)+kn]>0)
						F[ci][lm]=f[ci][LR[lm]]-feq(LR[lm],rho[ci],v_ls)+(1-2*q_p)*feq(LR[lm],recvr_rho[Sr[jn*(NZ+1)+kn]-1],v_ls2)+2*q_p*feq(LR[lm],rho[ci],v_ls);
					else
						F[ci][lm]=(qprim+2*q_p-1)/qprim*(feq(LR[lm],rho[ci],v_ls))+(1-2*q_p)/qprim*w[LR[lm]];
						

					}
			}
					
			if ((ip<nx_l) and (in<nx_l))
			{
				
				if (Solid[ip][jp][kp]>0)
					F[ci][lm]=f[abs(Solid[ip][jp][kp])][lm];
				else
					{
					v_ls[0]=u[ci][0];v_ls2[0]=u[Solid[in][jn][kn]][0];
					v_ls[1]=u[ci][1];v_ls2[1]=u[Solid[in][jn][kn]][1];
					v_ls[2]=u[ci][2];v_ls2[2]=u[Solid[in][jn][kn]][2];
	
					if (Solid[in][jn][kn]>0)
						F[ci][lm]=f[ci][LR[lm]]-feq(LR[lm],rho[ci],v_ls)+(1-2*q_p)*feq(LR[lm],rho[Solid[in][jn][kn]],v_ls2)+2*q_p*feq(LR[lm],rho[ci],v_ls);
					else
						F[ci][lm]=(qprim+2*q_p-1)/qprim*(feq(LR[lm],rho[ci],v_ls))+(1-2*q_p)/qprim*w[LR[lm]];
						
						
					} 	


			}
			
		}
		
		}

	}






	//for(int ci=1;ci<=Count;ci++)	
	//	for(int lm=0;lm<19;lm++)
        //        	f[ci][lm]=F[ci][lm];


	delete [] Gcl;
	delete [] Gcr;

	
	delete [] recvl;
	delete [] recvr;
			
	if (q_p<0.5)
		{
		delete [] recvl_rho;
		delete [] recvl_u;
		delete [] recvr_rho;
		delete [] recvr_u;
		
		}						
						
	

}

void periodic_streaming(double** f,double** F,int* SupInv,int*** Solid,int* Sl,int* Sr,double* rho,double** u)
{
	MPI_Status status[4] ;
	MPI_Request request[4];

	//cout<<"@@@@@@@@@@@   "<<n<<endl;
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	int ip,jp,kp,i,j,k,mpi_test;
	double rho_ls,v_ls[3];
	
	int* Gcl = new int[mpi_size];
	int* Gcr = new int[mpi_size];

	
	MPI_Gather(&cl,1,MPI_INT,Gcl,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(Gcl,mpi_size,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Gather(&cr,1,MPI_INT,Gcr,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(Gcr,mpi_size,MPI_INT,0,MPI_COMM_WORLD);


	double* recvl;
	double* recvr;

	double* sendl = new double[Gcl[rank]*19];
	double* sendr = new double[Gcr[rank]*19];


		for(k=0;k<19;k++)
			{
			for(i=1;i<=Gcl[rank];i++)
				sendl[(i-1)*19+k]=f[i][k];
			for(j=Count-Gcr[rank]+1;j<=Count;j++)
				sendr[(j-(Count-Gcr[rank]+1))*19+k]=f[j][k];
			}


	
	if (rank==0)
		{
		recvl = new double[Gcr[mpi_size-1]*19];
		recvr = new double[Gcl[rank+1]*19];
		}	
		else
		if (rank==mpi_size-1)
			{
			recvl = new double[Gcr[rank-1]*19];
			recvr = new double[Gcl[0]*19];
			}
			else
			{
			recvl = new double[Gcr[rank-1]*19];

			recvr = new double[Gcl[rank+1]*19];
			}

MPI_Barrier(MPI_COMM_WORLD);

//cout<<"@@@@@@@@@@@   "<<n<<endl;
if (rank==0)
		{
		
		MPI_Isend(sendr, Gcr[0]*19, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(sendl, Gcl[0]*19, MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(recvr , Gcl[1]*19, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(recvl, Gcr[mpi_size-1]*19, MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );
		
		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(sendl, Gcl[rank]*19, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr, Gcr[rank]*19, MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl, Gcr[rank-1]*19, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(recvr, Gcl[0]*19, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request[3]);
			
			}
			else
			{

			MPI_Isend(sendl, Gcl[rank]*19, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr, Gcr[rank]*19, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl, Gcr[rank-1]*19, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(recvr, Gcl[rank+1]*19, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);

			
			};

	
	MPI_Waitall(4,request, status);

	MPI_Testall(4,request,&mpi_test,status);


//cout<<"@@@@@@@@@@@   "<<n<<endl;
		for(int ci=1;ci<=Count;ci++)	
		{
				i=(int)(SupInv[ci]/((NY+1)*(NZ+1)));
				j=(int)((SupInv[ci]%((NY+1)*(NZ+1)))/(NZ+1));
				k=(int)(SupInv[ci]%(NZ+1));                 
		       for(int lm=0;lm<19;lm++)


			{       
				

				ip=i-e[lm][0];
				jp=j-e[lm][1];if (jp<0) {jp=NY;}; if (jp>NY) {jp=0;};
				kp=k-e[lm][2];if (kp<0) {kp=NZ;}; if (kp>NZ) {kp=0;};

			
				if (ip<0) 
					if (Sl[jp*(NZ+1)+kp]>0)
						F[ci][lm]=recvl[(Sl[jp*(NZ+1)+kp]-1)*19+lm];
					else
						F[ci][lm]=f[ci][LR[lm]];
					
						
					
					
				if (ip>=nx_l)
					if (Sr[jp*(NZ+1)+kp]>0)
						F[ci][lm]=recvr[(Sr[jp*(NZ+1)+kp]-1)*19+lm];
					else
						F[ci][lm]=f[ci][LR[lm]];

				if ((ip>=0) and (ip<nx_l)) 
					if (Solid[ip][jp][kp]>0)
						F[ci][lm]=f[Solid[ip][jp][kp]][lm];
					else
						F[ci][lm]=f[ci][LR[lm]];
	
					
			}
		}
			

/*
	double sum=0;
	for(int ci=1;ci<=Count;ci++)	
		{
		for(int lm=0;lm<19;lm++)
                	{f[ci][lm]=F[ci][lm];sum+=F[ci][lm];}
		cout<<sum<<endl;sum=0;
		}
*/	

	delete [] Gcl;
	delete [] Gcr;

	delete [] sendl;
	delete [] sendr;
	delete [] recvl;
	delete [] recvr;
			
							
						

}


void collision(double* rho,double** u,double** f,double* forcex, double* forcey, double* forcez, int* SupInv,int*** Solid)
{

double lm0,lm1;
double usqr,vsqr;
double F_hat[19],GuoF[19],f_eq[19],u_tmp[3];
double m_l[19];

	
	

	for(int i=1;i<=Count;i++)	
	

		{	
				
			
			

			//=================FORCE TERM_GUO=========================================
			for (int k=0;k<19;k++)
			{	
			lm0=((e[k][0]-u[i][0])*forcex[i]+(e[k][1]-u[i][1])*forcey[i]+(e[k][2]-u[i][2])*forcez[i])*3;
			lm1=(e[k][0]*u[i][0]+e[k][1]*u[i][1]+e[k][2]*u[i][2])*(e[k][0]*forcex[i]+e[k][1]*forcey[i]+e[k][2]*forcez[i])*9;
			GuoF[k]=w[k]*(lm0+lm1);
			//GuoF[k]=0.0;
			}


			
			//=====================equilibrium of moment=================================
			u_tmp[0]=u[i][0];
			u_tmp[1]=u[i][1];
			u_tmp[2]=u[i][2];

			for(int k=0;k<19;k++)
				{
				f_eq[k]=feq(k,rho[i],u_tmp);
				}
			for (int l=0;l<19;l++)
				{
				meq[l]=0;
				for(int lm=0;lm<19;lm++)
				meq[l]+=M[l][lm]*f_eq[lm];				
				}

			//============================================================================

			
			// ==================   m=Mf matrix calculation  =============================
			// ==================   F_hat=(I-.5*S)MGuoF =====================================
				for (int mi=0; mi<19; mi++)
					{m_l[mi]=0;F_hat[mi]=0;
					for (int mj=0; mj<19; mj++)
						{
						m_l[mi]+=M[mi][mj]*f[i][mj];
						F_hat[mi]+=M[mi][mj]*GuoF[mj];
						}
					F_hat[mi]*=(1-0.5*S[mi]);
					}
			//============================================================================


//			if (i==20)
//				{
//				cout<<"Collison"<<endl;							
//							for(int lm=0;lm<19;lm++)
//									cout<<f[20][lm]<<", ";
//							cout<<endl;
//				}
				
			
				
				for (int sk=0;sk<19;sk++)
				//m[sk]=m[sk]-S[sk]*(m[sk]-meq[sk]);
				//m[sk]=m[sk]-S[sk]*(m[sk]-meq[sk])+(1-S[sk]*F_hat[sk]/2);
				m_l[sk]=m_l[sk]-S[sk]*(m_l[sk]-meq[sk])+dt*F_hat[sk];

			//}

			// ==================   f=M_-1m matrix calculation  =============================
				for (int mi=0; mi<19; mi++)
					{f[i][mi]=0;
					for (int mj=0; mj<19; mj++)
						f[i][mi]+=MI[mi][mj]*m_l[mj];
					}
			//============================================================================
			
			
			}
                        
/*	
	double sum=0;
	for(int ci=1;ci<=Count;ci++)	
		{
		for(int lm=0;lm<19;lm++)
                	{sum+=f[ci][lm];}
		cout<<sum<<endl;sum=0;
		}
			
*/			
		

}


void standard_bounceback_boundary(int it,double** f)
{

	double tmp;
			tmp = f[it][1];f[it][1] = f[it][2];f[it][2] = tmp;
			tmp = f[it][3];f[it][3] = f[it][4];f[it][4] = tmp;
                        tmp = f[it][5];f[it][5] = f[it][6];f[it][6] = tmp;
			tmp = f[it][7];f[it][7] = f[it][10];f[it][10] = tmp;
			tmp = f[it][8];f[it][8] = f[it][9];f[it][9] = tmp;
			tmp = f[it][11];f[it][11] = f[it][14];f[it][14] = tmp;
                        tmp = f[it][12];f[it][12] = f[it][13];f[it][13] = tmp;
			tmp = f[it][15];f[it][15] = f[it][18];f[it][18] = tmp;
			tmp = f[it][16];f[it][16] = f[it][17];f[it][17] = tmp;


			

}



void comput_macro_variables_IMR( double* rho,double** u,double** u0,double** f,double** F,int* SupInv,int*** Solid, double* forcex, double* forcey, double* forcez, int n,double porosity,double* gx,double* gy, double* gz)
{
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	double rho0=1.0;
	double dp[3],rhok;
	double *rbuf;
	rbuf=new double[mpi_size*3];

	for(int i=1;i<=Count;i++)	
                  
			{ 

				
				rhok=rho[i];
				u0[i][0]=u[i][0];
				u0[i][1]=u[i][1];
				u0[i][2]=u[i][2];
				rho[i]=0;
				u[i][0]=0;
				u[i][1]=0;
				u[i][2]=0;
	
				for(int k=0;k<19;k++)
					{
					f[i][k]=F[i][k];
					rho[i]+=f[i][k];
					u[i][0]+=e[k][0]*f[i][k];
					u[i][1]+=e[k][1]*f[i][k];
					u[i][2]+=e[k][2]*f[i][k];
					}
			

				u[i][0]=(u[i][0]+dt*forcex[i]/2)/rho[i];
				u[i][1]=(u[i][1]+dt*forcey[i]/2)/rho[i];
				u[i][2]=(u[i][2]+dt*forcez[i]/2)/rho[i];
					

			

				dp[0]+=forcex[i]-(u[i][0]-u0[i][0])*rhok;
				dp[1]+=forcey[i]-(u[i][1]-u0[i][1])*rhok;
				dp[2]+=forcez[i]-(u[i][2]-u0[i][2])*rhok;
				
			}
			

                     

	if ((n%1000==0) and (n>100))
	{
	
	MPI_Barrier(MPI_COMM_WORLD);   
	MPI_Gather(dp,3,MPI_DOUBLE,rbuf,3,MPI_DOUBLE,0,MPI_COMM_WORLD);

	if (rank==0)
		{
		dp[0]=0;dp[1]=0;dp[2]=0;
		for (int i=0;i<mpi_size;i++)
			{
			dp[0]+=rbuf[i*3+0];
			dp[1]+=rbuf[i*3+1];
			dp[2]+=rbuf[i*3+2];
			}
		dp[0]/=(NX+1)*(NY+1)*(NZ+1)*porosity;
		dp[1]/=(NX+1)*(NY+1)*(NZ+1)*porosity;
		dp[2]/=(NX+1)*(NY+1)*(NZ+1)*porosity;
		}
	

		
	MPI_Bcast(dp,3,MPI_DOUBLE,0,MPI_COMM_WORLD);	
		for (int i=1;i<=Count;i++)
			{
		switch(PerDir)
		{
		case 1:
			forcex[i]=dp[0];break;
		case 2:
			forcey[i]=dp[1];break;
		case 3:
			forcez[i]=dp[2];break;
		default:
			forcex[i]=dp[0];
		}

			}
	*gx=dp[0];
	*gy=dp[1];
	*gz=dp[2];
	
	}
	delete [] rbuf;

}


void comput_macro_variables( double* rho,double** u,double** u0,double** f,double** F,double* forcex, double* forcey, double* forcez,int* SupInv,int*** Solid)
{
//	int rank = MPI :: COMM_WORLD . Get_rank ();
//	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	
	

	for(int i=1;i<=Count;i++)	
                   
			{
			
				u0[i][0]=u[i][0];
				u0[i][1]=u[i][1];
				u0[i][2]=u[i][2];
				rho[i]=0;
				u[i][0]=0;
				u[i][1]=0;
				u[i][2]=0;
	
				for(int k=0;k<19;k++)
					{
					
					f[i][k]=F[i][k];
					rho[i]+=f[i][k];
					u[i][0]+=e[k][0]*f[i][k];
					u[i][1]+=e[k][1]*f[i][k];
					u[i][2]+=e[k][2]*f[i][k];
					}
				
				//cout<<f[i][2]<<endl;

				u[i][0]=(u[i][0]+dt*forcex[i]/2)/rho[i];
				u[i][1]=(u[i][1]+dt*forcey[i]/2)/rho[i];
				u[i][2]=(u[i][2]+dt*forcez[i]/2)/rho[i];
				
		//=============================DEBUG=======================================
	/*
				if ((n%10==0) and (n<=1000))
					{
					forcex[i]=gx*0.01*n/10;
					forcey[i]=gy*0.01*n/10;
					forcez[i]=gz*0.01*n/10;
					
					}
*/

		//===========================================================================


		
			}
			
	
                        
			
	//MPI_Barrier(MPI_COMM_WORLD); 

}



void boundary_velocity(int xp,double v_xp,int xn, double v_xn,int yp,double v_yp,int yn,double v_yn,int zp,double v_zp,int zn,double v_zn,double** f,double* rho,double** u,int*** Solid)

// PLEASE SPECIFY THE BOUNDARY LOCATION BEFORE USE
// FOR EXAMPLE: IF SOUTH BOUNDARY IS INCLUDE, THEN THE BOOL VARIABLE SOUTH=1, AND ALSO SPECIFY THE LOCAL RHO
// IF THIS BOUNDARY CONDITION IS NOT GOING TO BE USED, PLEASE SET THE CORRESPONDING BOOL VARIABLE AS 0
// FORMAT:
// (NORTH BOUNDARY MARK,LOCAL RHO VALUE,SOUTH BOUNDARY MARK, LOCAL RHO VALUE,EAST BOUNDARY MARK, LOCAL RHO VALUE
// WEST BOUNDARY MARK, LOCAL RHO VALUE)

{

int Q=19;
//double ux0,uy0,uz0;

double u_xp[3]={v_xp,0,0};
double u_xn[3]={v_xn,0,0};
double u_yp[3]={0,v_yp,0};
double u_yn[3]={0,v_yn,0};
double u_zp[3]={0,0,v_zp};
double u_zn[3]={0,0,v_zn};

int rank = MPI :: COMM_WORLD . Get_rank ();
int mpi_size=MPI :: COMM_WORLD . Get_size ();


if (yp==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
			{
			if (Solid[i][NY][k]>0)
				if (Solid[i][NY-1][k]>0) 
					f[Solid[i][NY][k]][ks]=feq(ks,1.0,u_yp)+f[Solid[i][NY-1][k]][ks]-feq(ks,1.0,u[Solid[i][NY-1][k]]);
				else
					f[Solid[i][NY][k]][ks]=feq(ks,1.0,u_yp);

			f[Solid[i][NY][k]][ks]=feq(ks,1.0,u_yp);
			}


if (yn==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
			{
			if (Solid[i][0][k]>0)
				if (Solid[i][1][k]>0) 
					f[Solid[i][0][k]][ks]=feq(ks,1.0,u_yn)+f[Solid[i][1][k]][ks]-feq(ks,1.0,u[Solid[i][1][k]]);
				else
					f[Solid[i][0][k]][ks]=feq(ks,1.0,u_yn);
			f[Solid[i][0][k]][ks]=feq(ks,1.0,u_yn);
			}



if (zp==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
		for (int ks=0;ks<Q;ks++)
			{
			if (Solid[i][j][NZ]>0)
				if (Solid[i][j][NZ-1]>0) 
					f[Solid[i][j][NZ]][ks]=feq(ks,1.0,u_zp)+f[Solid[i][j][NZ-1]][ks]-feq(ks,1.0,u[Solid[i][j][NZ-1]]);
				else
					f[Solid[i][j][NZ]][ks]=feq(ks,1.0,u_zp);
			f[Solid[i][j][NZ]][ks]=feq(ks,rho[abs(Solid[i][j][NZ-1])],u_zp);
			}



if (zn==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
		for (int ks=0;ks<Q;ks++)
			{
			if (Solid[i][j][0]>0)
				if (Solid[i][j][1]>0) 
					f[Solid[i][j][0]][ks]=feq(ks,1.0,u_zn)+f[Solid[i][j][1]][ks]-feq(ks,1.0,u[Solid[i][j][1]]);
				else
					f[Solid[i][j][0]][ks]=feq(ks,1.0,u_zn);
			f[Solid[i][j][0]][ks]=feq(ks,rho[abs(Solid[i][j][1])],u_zn);
			}




if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
			{
			if (Solid[nx_l-1][k][j]>0)
				if (Solid[nx_l-2][j][k]>0) 
					f[Solid[nx_l-1][j][k]][ks]=feq(ks,1.0,u_xp)+f[Solid[nx_l-2][j][k]][ks]-feq(ks,1.0,u[Solid[nx_l-2][j][k]]);
				else
					f[Solid[nx_l-1][j][k]][ks]=feq(ks,1.0,u_xp);
			f[Solid[nx_l-1][j][k]][ks]=feq(ks,rho[abs(Solid[nx_l-2][j][k])],u_xp);
			}



if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
			{
			if (Solid[0][j][k]>0)
				if (Solid[1][j][k]>0) 
					f[Solid[0][j][k]][ks]=feq(ks,1.0,u_xn)+f[Solid[1][j][k]][ks]-feq(ks,1.0,u[Solid[1][j][k]]);
				else
					f[Solid[0][j][k]][ks]=feq(ks,1.0,u_xn);
			f[Solid[0][j][k]][ks]=feq(ks,rho[abs(Solid[1][j][k])],u_xn);
			}
	


}



void boundary_pressure(int xp,double rho_xp,int xn, double rho_xn,int yp,double rho_yp,int yn,double rho_yn,int zp,double rho_zp,int zn,double rho_zn,double** f,double** u,double* rho,int*** Solid)

// PLEASE SPECIFY THE BOUNDARY LOCATION BEFORE USE
// FOR EXAMPLE: IF SOUTH BOUNDARY IS INCLUDE, THEN THE BOOL VARIABLE SOUTH=1, AND ALSO SPECIFY THE LOCAL RHO
// IF THIS BOUNDARY CONDITION IS NOT GOING TO BE USED, PLEASE SET THE CORRESPONDING BOOL VARIABLE AS 0
// FORMAT:
// (NORTH BOUNDARY MARK,LOCAL RHO VALUE,SOUTH BOUNDARY MARK, LOCAL RHO VALUE,EAST BOUNDARY MARK, LOCAL RHO VALUE
// WEST BOUNDARY MARK, LOCAL RHO VALUE)

{

int Q=19;
//double ux0,uy0,uz0;
double u_ls[3]={0,0,0};
int rank = MPI :: COMM_WORLD . Get_rank ();
int mpi_size=MPI :: COMM_WORLD . Get_size ();


if (yp==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
	{ 	

		for (int ks=0;ks<Q;ks++)
			
			if (Solid[i][NY][k]>0)
			{
				if (Solid[i][NY-1][k]>0)
					{
					u_ls[0]=u[Solid[i][NY-1][k]][0];
					u_ls[1]=u[Solid[i][NY-1][k]][1];
					u_ls[2]=u[Solid[i][NY-1][k]][2];
					//f[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls)+f[Solid[i][NY-1][k]][ks]-feq(ks,rho[Solid[i][NY-1][k]],u[Solid[i][NY-1][k]]);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					//f[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls);
					}
			f[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls);
			}
	}


if (yn==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
	{ 	
	
		for (int ks=0;ks<Q;ks++)
			if (Solid[i][0][k]>0)
			{
				if (Solid[i][1][k]>0)
					{
					u_ls[0]=u[Solid[i][1][k]][0];
					u_ls[1]=u[Solid[i][1][k]][1];
					u_ls[2]=u[Solid[i][1][k]][2];
					//f[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls)+f[Solid[i][1][k]][ks]-feq(ks,rho[Solid[i][1][k]],u[Solid[i][1][k]]);
					}
			else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					//f[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls);
					}
			f[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls);
			}
	}


if (zp==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	{ 	
		
		for (int ks=0;ks<Q;ks++)
			
			if (Solid[i][j][NZ]>0)
			{			
				if (Solid[i][j][NZ-1]>0)
					{
					u_ls[0]=u[Solid[i][j][NZ-1]][0];
					u_ls[1]=u[Solid[i][j][NZ-1]][1];
					u_ls[2]=u[Solid[i][j][NZ-1]][2];
					//f[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls)+f[Solid[i][j][NZ-1]][ks]-feq(ks,rho[Solid[i][j][NZ-1]],u[Solid[i][j][NZ-1]]);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					//f[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls);
					}

			f[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls);
			}
	}



if (zn==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	{ 	
		
		for (int ks=0;ks<Q;ks++)
			
			if (Solid[i][j][0]>0)
			{
				if (Solid[i][j][1]>0)
					{
					u_ls[0]=u[Solid[i][j][1]][0];
					u_ls[1]=u[Solid[i][j][1]][1];
					u_ls[2]=u[Solid[i][j][1]][2];
					//f[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls)+f[Solid[i][j][1]][ks]-feq(ks,rho[Solid[i][j][1]],u[Solid[i][j][1]]);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;//f[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls);
					}
			f[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls);
			}
	}



if ((xp==1) && (rank==mpi_size-1))
      for (int j=0;j<=NY;j++)
		for (int k=0;k<=NZ;k++)
	{	
		
		for (int ks=0;ks<Q;ks++)
			
			if (Solid[nx_l-1][j][k]>0)
			{
				if (Solid[nx_l-2][j][k]>0)
					{
					u_ls[0]=u[Solid[nx_l-2][j][k]][0];
					u_ls[1]=u[Solid[nx_l-2][j][k]][1];
					u_ls[2]=u[Solid[nx_l-2][j][k]][2];
					//f[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls)+f[Solid[nx_l-2][j][k]][ks]-feq(ks,rho[Solid[nx_l-2][j][k]],u[Solid[nx_l-2][j][k]]);
					//f[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls);
					}
				else	
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;
					u_ls[2]=0.0;
					//f[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls);
					}
			f[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls);
			}
			
			

	}


if ((xn==1) && (rank==0))
      for (int j=0;j<=NY;j++)
		for(int k=0;k<=NZ;k++)
	{	
		//cout<<"@@@@@@@@@@@"<<endl;
		for (int ks=0;ks<Q;ks++)
			
			if (Solid[0][j][k]>0)
			{
				if(Solid[1][j][k]>0)
					{
					u_ls[0]=u[Solid[1][j][k]][0];
					u_ls[1]=u[Solid[1][j][k]][1];
					u_ls[2]=u[Solid[1][j][k]][1];
					//f[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls)+f[Solid[1][j][k]][ks]-feq(ks,rho[Solid[1][j][k]],u[Solid[1][j][k]]);
					//f[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;
					u_ls[2]=0.0;
					//f[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls);
					}

			f[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls);
			}
			//f[0][j][k][ks]=feq(ks,rho[1][j][k],u_ls,w,e);
			
	

	}

}


double Error(double** u,double** u0,double *v_max,double* u_average)
{	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	double *rbuf,*um,*uave,u_compt;
	double temp1,temp2,temp3;
	temp1=0;
	temp2=0;
	temp3=0;
	double error_in;
	double u_max;
	

	rbuf=new double[mpi_size];
	um = new double[mpi_size];
	uave = new double[mpi_size];
	u_max=0;
	*v_max=0;




//MPI_Barrier(MPI_COMM_WORLD);



for(int i=1; i<Count; i++)
			{	

			temp1+=(u[i][0]-u0[i][0])*(u[i][0]-u0[i][0])+(u[i][1]-u0[i][1])*(u[i][1]-u0[i][1])+(u[i][2]-u0[i][2])*(u[i][2]-u0[i][2]);
			temp2 += u[i][0]*u[i][0]+u[i][1]*u[i][1]+u[i][2]*u[i][2];	
			temp3+=sqrt(u[i][0]*u[i][0]+u[i][1]*u[i][1]+u[i][2]*u[i][2]);
			if (u[i][0]*u[i][0]+u[i][1]*u[i][1]+u[i][2]*u[i][2]>u_max)
				u_max=u[i][0]*u[i][0]+u[i][1]*u[i][1]+u[i][2]*u[i][2];
			
			}
		
		temp1=sqrt(temp1);
		temp2=sqrt(temp2);
		error_in=temp1/(temp2+1e-30);
		u_max=sqrt(u_max);

		MPI_Barrier(MPI_COMM_WORLD);
		
	MPI_Gather(&error_in,1,MPI_DOUBLE,rbuf,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
	MPI_Gather(&u_max,1,MPI_DOUBLE,um,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Gather(&temp3,1,MPI_DOUBLE,uave,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	u_compt=0;
	if (rank==0)
	    for (int i=0;i<mpi_size;i++)
		{
		if (rbuf[i]>error_in)
			error_in=rbuf[i];
		if (um[i]>*v_max)
			*v_max=um[i];
		u_compt+=uave[i];

		}

	u_compt/=(NX+1)*(NY+1)*(NZ+1);
	*u_average=u_compt;
	
	delete [] rbuf;
	delete [] um;
	delete [] uave;

	//MPI_Barrier(MPI_COMM_WORLD);

	return(error_in);

}


/*
void Geometry(int*** Solid)	
{	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	const int root_rank=0;

	
	int* nx_g = new int[mpi_size];
	int* disp = new int[mpi_size];

	
	MPI_Gather(&nx_l,1,MPI_INT,nx_g,1,MPI_INT,root_rank,MPI_COMM_WORLD);
	
	
	if (rank==root_rank)
		{
		disp[0]=0;
		for (int i=0;i<mpi_size;i++)
			nx_g[i]*=(NY+1)*(NZ+1);

		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
		}


		
	
	int* Solid_storage= new int[nx_l*(NY+1)*(NZ+1)];
	int* rbuf;

	for(int i=0;i<nx_l;i++)
		for(int j=0;j<=NY;j++)
			for(int k=0;k<=NZ;k++)
			if (Solid[i][j][k]<=0)
				Solid_storage[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=1;
			else
				Solid_storage[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=0;

	if (rank==root_rank)
		rbuf= new int[(NX+1)*(NY+1)*(NZ+1)];
	
	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gatherv(Solid_storage,nx_l*(NY+1)*(NZ+1),MPI_INT,rbuf,nx_g,disp,MPI_INT,root_rank,MPI_COMM_WORLD);
	
	int NX0=NX+1;
	int NY0=NY+1;
	int NZ0=NZ+1;

if (mir==0)
	{	
	if (mirX==1)
		NX0=NX0/2;
	if (mirY==1)
		NY0=NY0/2;
	if (mirZ==1)
		NZ0=NZ0/2;
	}



	if (rank==root_rank)
	{
	ostringstream name;
	name<<outputfile<<"LBM_Geometry"<<".vtk";
	ofstream out;
	out.open(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Geometry"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<NX0<<"         "<<NY0<<"         "<<NZ0<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<NX0*NY0*NZ0<<endl;
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
for(int k=0;k<NZ0;k++)
        for(int j=0; j<NY0; j++)
		for(int i=0;i<NX0;i++)
			//for(int k=0;k<=NZ;k++)
			out<<"		"<<rbuf[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]<<endl;
	out.close();
	
	}
		
	delete [] Solid_storage;
	if (rank==root_rank)
		delete [] rbuf;

	delete [] nx_g;
	delete [] disp;

		
}

*/



void Geometry(int*** Solid)	
{	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	const int root_rank=0;

	
	int nx_g[mpi_size];
	int disp[mpi_size];

	
	MPI_Gather(&nx_l,1,MPI_INT,nx_g,1,MPI_INT,root_rank,MPI_COMM_WORLD);
	
	if (rank==root_rank)
		{
		disp[0]=0;
	
		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
		}

	
	MPI_Bcast(&disp,mpi_size,MPI_INT,0,MPI_COMM_WORLD);
	
	
	
	int NX0=NX+1;
	int NY0=NY+1;
	int NZ0=NZ+1;

if (mir==0)
	{	
	if (mirX==1)
		NX0=NX0/2;
	if (mirY==1)
		NY0=NY0/2;
	if (mirZ==1)
		NZ0=NZ0/2;
	}


	ostringstream name;
	name<<outputfile<<"LBM_Geometry"<<".vtk";
	if (rank==root_rank)
	{

	ofstream out;
	out.open(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Geometry"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<NZ0<<"         "<<NY0<<"         "<<NX0<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<NX0*NY0*NZ0<<endl;
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	out.close();
	}
	
	


	for (int processor=0;processor<mpi_size;processor++)
	{
		if (rank==processor) 
		{
		ofstream out(name.str().c_str(),ios::app);
		for(int i=0;i<nx_l;i++)
			if (i+disp[rank]<NX0)
        		for(int j=0; j<NY0; j++)
				for(int k=0;k<NZ0;k++)
				//	out<<1<<endl;
				if (Solid[i][j][k]<=0)
					out<<1<<endl;
				else
					out<<0<<endl;
		out.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	
		
}




void output_velocity(int m,double* rho,double** u,int MirX,int MirY,int MirZ,int mir,int*** Solid)	
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	const int mpi_size=MPI :: COMM_WORLD . Get_size ();
	const int root_rank=0;
	
	int nx_g[mpi_size];
	int disp[mpi_size];

	
	MPI_Gather(&nx_l,1,MPI_INT,nx_g,1,MPI_INT,root_rank,MPI_COMM_WORLD);
	
	
	if (rank==root_rank)
		{
		disp[0]=0;

		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
		}
	

	

	int NX0=NX+1;
	int NY0=NY+1;
	int NZ0=NZ+1;

if (mir==0)
	{	
	if (MirX==1)
		NX0=NX0/2;
	if (MirY==1)
		NY0=NY0/2;
	if (MirZ==1)
		NZ0=NZ0/2;
	}


	ostringstream name;
	name<<outputfile<<"LBM_velocity_Vector_"<<m<<".vtk";
	if (rank==root_rank)
	{

	ofstream out;
	out.open(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Velocity"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<NZ0<<"         "<<NY0<<"         "<<NX0<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;

	out<<"POINT_DATA     "<<NX0*NY0*NZ0<<endl;
	out<<"VECTORS sample_vectors double"<<endl;
	out<<endl;

	out.close();
	}


	
			
	
	for (int processor=0;processor<mpi_size;processor++)
	{
		if (rank==processor) 
		{
		ofstream out(name.str().c_str(),ios::app);
		for(int i=0;i<nx_l;i++)
			if (i+disp[rank]<NX0)
        		for(int j=0; j<NY0; j++)
				for(int k=0;k<NZ0;k++)
				if (Solid[i][j][k]>0)
				out<<u[Solid[i][j][k]][2]<<" "<<u[Solid[i][j][k]][1]<<" "<<u[Solid[i][j][k]][0]<<endl;
				else
				out<<0.0<<" "<<0.0<<" "<<0.0<<endl;


		out.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}



	
/*	ostringstream name2;
	name2<<"LBM_velocity_"<<m<<".out";
	ofstream out2(name2.str().c_str());
	for (int j=0;j<=NY;j++)
		{
		if (Solid[1][j][1]>0)
			out2<<u[Solid[2][j][1]][0]<<endl;
		else
			out2<<0.0<<endl;
		}
*/
	
	

	

		
}


void output_density(int m,double* rho,int MirX,int MirY,int MirZ,int mir,int*** Solid)	
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	const int mpi_size=MPI :: COMM_WORLD . Get_size ();
	const int root_rank=0;
	

	int nx_g[mpi_size];
	int disp[mpi_size];

	
	MPI_Gather(&nx_l,1,MPI_INT,nx_g,1,MPI_INT,root_rank,MPI_COMM_WORLD);
	
	
	if (rank==root_rank)
		{
		disp[0]=0;
	
		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
		}
	

	
	int NX0=NX+1;
	int NY0=NY+1;
	int NZ0=NZ+1;

	if (mir==0)
	{	
	if (MirX==1)
		NX0=NX0/2;
	if (MirY==1)
		NY0=NY0/2;
	if (MirZ==1)
		NZ0=NZ0/2;
	}

	ostringstream name;
	name<<outputfile<<"LBM_Density_"<<m<<".vtk";
	if (rank==root_rank)
	{
	ofstream out;
	out.open(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<NZ0<<"         "<<NY0<<"         "<<NX0<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<NX0*NY0*NZ0<<endl;
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	out.close();

	}
        

				
	for (int processor=0;processor<mpi_size;processor++)
	{
		if (rank==processor) 
		{
		ofstream out(name.str().c_str(),ios::app);
		for(int i=0;i<nx_l;i++)
			if (i+disp[rank]<NX0)
        		for(int j=0; j<NY0; j++)
				for(int k=0;k<NZ0;k++)
				if (Solid[i][j][k]>0)
				out<<rho[Solid[i][j][k]]<<endl;
				else
				out<<1.0<<endl;
		out.close();
		}
		MPI_Barrier(MPI_COMM_WORLD);
	}
	
	

		
}



double Comput_Perm(double** u,double* Permia,int PerDIr)
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	double *rbuf;
	rbuf=new double[mpi_size*3];
	double Perm[3];
	double error;
	double Q[3]={0.0,0.0,0.0};

	
	for (int i=1;i<=Count;i++)
		{
		Q[0]+=u[i][0];
		Q[1]+=u[i][1];
		Q[2]+=u[i][2];
		}

	MPI_Barrier(MPI_COMM_WORLD);

	//Qx/=(NX+1)/mpi_size*(NY+1)*(NZ+1);
	//Qy/=(NX+1)/mpi_size*(NY+1)*(NZ+1);
	//Qz/=(NX+1)/mpi_size*(NY+1)*(NZ+1);

	MPI_Gather(&Q,3,MPI_DOUBLE,rbuf,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	if (rank==0)
		{
		Q[0]=0;Q[1]=0;Q[2]=0;
		for (int i=0;i<mpi_size;i++)
			{
			Q[0]+=rbuf[i*3+0];
			Q[1]+=rbuf[i*3+1];
			Q[2]+=rbuf[i*3+2];
			}

		//Perm[0]=Q[0]*(1.0/Zoom)*(1.0/Zoom)/((NX+1)/Zoom*(NY+1)/Zoom*(NZ+1)/Zoom)*(in_vis)/gx;
		//Perm[1]=Q[1]*(1.0/Zoom)*(1.0/Zoom)/((NX+1)/Zoom*(NY+1)/Zoom*(NZ+1)/Zoom)*(in_vis)/gy;
		//Perm[2]=Q[2]*(1.0/Zoom)*(1.0/Zoom)/((NX+1)/Zoom*(NY+1)/Zoom*(NZ+1)/Zoom)*(in_vis)/gz;

		Perm[0]=Q[0]/((NX+1)*(NY+1)*(NZ+1))*(in_vis)/gx;
		Perm[1]=Q[1]/((NX+1)*(NY+1)*(NZ+1))*(in_vis)/gy;
		Perm[2]=Q[2]/((NX+1)*(NY+1)*(NZ+1))*(in_vis)/gz;
		
		switch(PerDIr)
		{
		case 1:
			error=(Perm[0]-Permia[0])/Permia[0];break;
		case 2:
			error=(Perm[1]-Permia[1])/Permia[1];break;
		case 3:
			error=(Perm[2]-Permia[2])/Permia[2];break;
		default:
			error=(Perm[0]-Permia[0])/Permia[0];
		}


		Permia[0]=Perm[0];
		Permia[1]=Perm[1];
		Permia[2]=Perm[2];

		}
	
	delete [] rbuf;
	
	return (error);
	

}

