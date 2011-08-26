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

double s_v=1.0;
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
double S[19];






//double m[19];
//double meq[19];





double uMax,c,Re,dx,dy,Lx,Ly,dt,rho0,P0,tau_f,niu,error,SFx,SFy,reso;

void Read_Rock(int***,double***,double*,char[128],char[128]);

void tests();

void init_Sparse(int***,int***,double***, double***, int*, int*);

void init(double*, double**, double**, double**,double*,double*, double*, double*,double***,int*);

void periodic_streaming(double** ,double** ,double**, double**, int* ,int***,int*, int*,double*, double**);

void periodic_streaming_Speed(double** ,double** ,double**, double**, int* ,int***,int*, int*,double*, double**);

void standard_bounceback_boundary(int,double**);

void collision(double* ,double* ,double** ,double** , double** ,double* , double* , double* , int*,int* ,int*,int***);

void comput_macro_variables(double* , double* ,double** ,double**,double** ,double** , double**,double** ,double* , double*, double*,int*,int***,double***);

double Error(double** ,double** ,double*, double*);

void output_velocity(int ,double* ,double** ,int ,int ,int ,int ,int*** );

void output_density(int ,double* ,int ,int ,int ,int ,int*** );	

void output_psi(int ,double* ,int ,int ,int ,int ,int*** );

void Geometry(int*** );	

void output_velocity_b(int ,double* ,double** ,int ,int ,int ,int ,int*** );

void output_density_b(int ,double* ,int ,int ,int ,int ,int*** );	

void output_psi_b(int ,double* ,int ,int ,int ,int ,int*** );

void Geometry_b(int*** );

double Comput_Perm(double* ,double** ,double* ,double* ,int );


double Comput_Saturation(double* ,int***);


void Comput_MI(double[19][19], double[19][19]);
int inverse(mat &a);
double feq(int,double, double[3]);

void Suppliment(int*,int***);

int e[19][3]=
{{0,0,0},{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};
double w[19]={1.0/3.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};

int LR[19]={0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15};

int n,nx_l,n_max,in_BC,PerDir,freRe,freDe,freVe,frePsi,Par_Geo,Par_nx,Par_ny,Par_nz;
int Zoom;

double ca=0.04;
double kappa,CM;

int wr_per,pre_xp,pre_xn,pre_yp,pre_yn,pre_zp,pre_zn;
int vel_xp,vel_xn,vel_yp,vel_yn,vel_zp,vel_zn,Sub_BC,Out_Mode;
double niu_l,niu_g,p_xp,p_xn,p_yp,p_yn,p_zp,p_zn,Re_l,Re_g;
double inivx,inivy,inivz,v_xp,v_xn,v_yp,v_yn,v_zp,v_zn,S_l,S_g,ContactAngle_parameter;
double error_perm,Permeability;

char outputfile[128]="./";


//=============================
double in_vis;
//=============================

int main(int argc , char *argv [])
{	

MPI :: Init (argc , argv );
double start , finish;

int rank = MPI :: COMM_WORLD . Get_rank ();
int para_size=MPI :: COMM_WORLD . Get_size ();

int dif;
double Per_l[3],Per_g[3];
double v_max,error_Per;


	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	int NCHAR=128;
	char     filename[128], dummy[128+1],filenamepsi[128];
	int      dummyInt;

	if (rank==0)
	{
	ifstream fin(argv[1]);
	                                                        fin.getline(dummy, NCHAR);
	fin >> filename;				fin.getline(dummy, NCHAR);
	fin >> filenamepsi;				fin.getline(dummy, NCHAR);
	fin >> NX >> NY >> NZ;			fin.getline(dummy, NCHAR);
	fin >> n_max;					fin.getline(dummy, NCHAR);
	fin >> reso;					fin.getline(dummy, NCHAR);
	fin >> in_BC;					fin.getline(dummy, NCHAR);
	fin >> mirX >> mirY >> mirZ;			fin.getline(dummy, NCHAR);
	fin >> gx >> gy >> gz;				fin.getline(dummy, NCHAR);
	fin >> pre_xp >> p_xp >> pre_xn >> p_xn;	fin.getline(dummy, NCHAR);
	fin >> pre_yp >> p_yp >> pre_yn >> p_yn;	fin.getline(dummy, NCHAR);
	fin >> pre_zp >> p_zp >> pre_zn >> p_zn;	fin.getline(dummy, NCHAR);
	fin >> vel_xp >> v_xp >> vel_xn >> v_xn;	fin.getline(dummy, NCHAR);
	fin >> vel_yp >> v_yp >> vel_yn >> v_yn;	fin.getline(dummy, NCHAR);
	fin >> vel_zp >> v_zp >> vel_zn >> v_zn;	fin.getline(dummy, NCHAR);
	fin >> niu_l;					fin.getline(dummy, NCHAR);
	fin >> niu_g;					fin.getline(dummy, NCHAR);
	fin >> ContactAngle_parameter;	fin.getline(dummy, NCHAR);
	fin >> kappa;					fin.getline(dummy, NCHAR);
	fin >> CM;					fin.getline(dummy, NCHAR);
	fin >> inivx >> inivy >> inivz;		fin.getline(dummy, NCHAR);
	fin >> Permeability;				fin.getline(dummy, NCHAR);
						        	fin.getline(dummy, NCHAR);
	fin >> wr_per;					fin.getline(dummy, NCHAR);
	fin >> PerDir;					fin.getline(dummy, NCHAR);
	fin >> Out_Mode;				fin.getline(dummy, NCHAR);
	fin >> freRe;					fin.getline(dummy, NCHAR);
	fin >> freVe;					fin.getline(dummy, NCHAR);
	fin >> freDe;					fin.getline(dummy, NCHAR);
	fin >> frePsi;					fin.getline(dummy, NCHAR);
	fin >> mir;					fin.getline(dummy, NCHAR);
						        	fin.getline(dummy, NCHAR);
	fin >> Par_Geo >> Par_nx >> Par_ny >> Par_nz;	fin.getline(dummy, NCHAR);
	fin >> Zoom;					fin.getline(dummy, NCHAR);
	fin >> outputfile;				fin.getline(dummy, NCHAR);

	fin.close();
		
	//cout<<CM<<"    @@@@@@@@@@asdfa "<<endl;
	NX=NX-1;NY=NY-1;NZ=NZ-1;
	}

	
	
	
	MPI_Bcast(&filename,128,MPI_CHAR,0,MPI_COMM_WORLD);MPI_Bcast(&filenamepsi,128,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&mirX,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&NX,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&mirY,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&NY,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&mirZ,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&NZ,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&n_max,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&reso,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&in_BC,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&gx,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
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
	MPI_Bcast(&mir,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&niu_l,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&niu_g,1,MPI_DOUBLE,0,MPI_COMM_WORLD); MPI_Bcast(&ContactAngle_parameter,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&frePsi,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&CM,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&Permeability,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Zoom,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&outputfile,128,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&Out_Mode,1,MPI_INT,0,MPI_COMM_WORLD);

      
int U_max_ref=0;


if (mirX==1)
	NX=NX*2+1;
if (mirY==1)
	NY=NY*2+1;
if (mirZ==1)
	NZ=NZ*2+1;


//cout<<NX<<"    @@@@@@@@@@asdfa  "<<Zoom<<endl;
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
//cout<<"@@@@@@@@@@@@@@@@@     "<<(NX+1)/Zoom<<"  "<<para_size<<endl;


	nx_l=(int)((NX+1)/para_size);
	dif=(NX+1)-nx_l*para_size;
	
	if (rank>para_size-1-dif)
		nx_l+=1;

//	if (rank==para_size-1)
//		nx_l+=(NX+1)%para_size;

	double* Permia;
	double* rho;
	double* psi;
	double** u;
	double** f;
	double** g;
	double** Fg;
	double** F;
	double** u0;
	int* SupInv;
	double* forcex;
	double* forcey;
	double* forcez;

	int*** Solids;
	double*** Psis;
	int*** Solid;
	double*** Psi_local;
	int*  Sl;
	int*  Sr;

	//f,g,Fg,Psis,psi
	//cout<<"asdfaaaafff      "<<(NX+1)/Zoom<<endl;
	
	Solid = new int**[nx_l];
	Psi_local = new double**[nx_l];
	Sl = new int[(NY+1)*(NZ+1)];
	Sr = new int[(NY+1)*(NZ+1)];

	Solids = new int**[(NX+1)/Zoom];
	Psis = new double**[(NX+1)/Zoom];

	for (int i=0;i<(NX+1)/Zoom;i++)
		{		
		Solids[i] = new int*[(NY+1)/Zoom];
		Psis[i] = new double*[(NY+1)/Zoom];
		
			for (int j=0;j<(NY+1)/Zoom;j++)
			{
			Solids[i][j]= new int[(NZ+1)/Zoom];
			Psis[i][j]= new double[(NZ+1)/Zoom];

			}
		}
		//cout<<"asdfaaaafff"<<endl;

	for (int i=0;i<nx_l;i++)
		{
		Solid[i] = new int*[NY+1];
		Psi_local[i] = new double*[NY+1];
			for (int j=0;j<=NY;j++)
			{
			Solid[i][j]= new int[NZ+1];
			Psi_local[i][j] = new double[NZ+1];
			}
		}

	

	
	
	
//	Count = new int[rank];
//	if (!(filename=="NOSOLID"))
		Read_Rock(Solids,Psis,&porosity,filename,filenamepsi);
//	else
//		{
//		for(int i=0;i<=NX;i++)	
//			for(int j=0;j<=NY;j++)
//				for(int k=0;k<=NZ;k++)
//				Solids[i][j][k]=0;
//		}

	

	init_Sparse(Solids,Solid,Psis,Psi_local,Sl,Sr);

	for (int i=0;i<(NX+1)/Zoom;i++)
		{
		for (int j=0;j<(NY+1)/Zoom;j++)
		        {
			delete [] Solids[i][j];
			delete [] Psis[i][j];
			}
		delete [] Solids[i];
		delete [] Psis[i];
		}
	delete [] Solids;
	delete [] Psis;
	

	//***************************************************
	//WARRING: SPARSE MATRIX STARTS FROM INDEX 1 NOT 0!!!
	//***************************************************

	Permia = new double[3];
	rho = new double[Count+1];
	psi = new double[Count+1];
	forcex = new double[Count+1];
	forcey = new double[Count+1];
	forcez = new double[Count+1];
	u = new double*[Count+1];
	f = new double*[Count+1];
	g =  new double*[Count+1];
	Fg = new double*[Count+1];
	F = new double*[Count+1];
	u0 = new double*[Count+1];
	SupInv = new int[Count+1];

	for (int i=0;i<=Count;i++)
		{
		u[i] = new double[3];
		f[i] = new double[19];
		u0[i] = new double[3];
		F[i] = new double[19];
		g[i] = new double[19];
		Fg[i] = new double[19];
		}

	Comput_MI(M,MI);
	
	Suppliment(SupInv,Solid);

	MPI_Barrier(MPI_COMM_WORLD);
	
	if (Out_Mode==1)
		Geometry(Solid);
	else
		Geometry_b(Solid);

	init(rho,u,f,g,psi,forcex,forcey,forcez,Psi_local,SupInv);

/*	
	for (int i=0;i<nx_l;i++)
		{
		for (int j=0;j<(NY+1);j++)
			delete [] Psi_local[i][j];
		delete [] Psi_local[i];
		}
	delete [] Psi_local;
*/	

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
char FileName4[128];strcpy(FileName4,outputfile);

strcat(FileName,"Results.txt");
strcat(FileName2,"Relative_Permeability_Component1.txt");
strcat(FileName4,"Relative_Permeability_Component2.txt");
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


//periodic_streamings(f,F,SupInv,Solid,Sl,Sr,rho,u);

	
	for(n=0;n<=n_max;n++)
	{
	
	//cout<<"@@@@@@@@@@@   "<<n<<endl;
	collision(rho,psi,u,f,g,forcex,forcey,forcez,Sl,Sr,SupInv,Solid);//cout<<"markcollision"<<endl;

	periodic_streaming(f,F,g,Fg,SupInv,Solid,Sl,Sr,rho,u);

/*	
	if ((1-pre_xp)*(1-pre_xn)*(1-pre_yp)*(1-pre_yn)*(1-pre_zp)*(1-pre_zn)==0)
		boundary_pressure(pre_xp,p_xp,pre_xn,p_xn,pre_yp,p_yp,pre_yn,p_yn,pre_zp,p_zp,pre_zn,p_zn,f,F,u,rho,Solid);

	
	if ((1-vel_xp)*(1-vel_xn)*(1-vel_yp)*(1-vel_yn)*(1-vel_zp)*(1-vel_zn)==0)
		boundary_velocity(vel_xp,v_xp,vel_xn,v_xn,vel_yp,v_yp,vel_yn,v_yn,vel_zp,v_zp,vel_zn,v_zn,f,F,rho,u,Solid);
*/

//	if (in_IMR==1)
//		comput_macro_variables_IMR(rho,u,u0,f,F,SupInv,Solid,forcex,forcey,forcez,n,porosity,&gx,&gy,&gz);
//	else

        comput_macro_variables(psi,rho,u,u0,f,g,F,Fg,forcex, forcey,forcez,SupInv,Solid,Psi_local);


	//if ((1-pre_xp)*(1-pre_xn)*(1-pre_yp)*(1-pre_yn)*(1-pre_zp)*(1-pre_zn)==0)
	//	boundary_pressure(pre_xp,p_xp,pre_xn,p_xn,pre_yp,p_yp,pre_yn,p_yn,pre_zp,p_zp,pre_zn,p_zn,f,u,rho,Solid);
		
	//if ((1-vel_xp)*(1-vel_xn)*(1-vel_yp)*(1-vel_yn)*(1-vel_zp)*(1-vel_zn)==0)
	//        boundary_velocity(vel_xp,v_xp,vel_xn,v_xn,vel_yp,v_yp,vel_yn,v_yn,vel_zp,v_zp,vel_zn,v_zn,f,rho,u,Solid);
	
	
	
	
	if(n%freRe==0)
		{       
		        error=Error(u,u0,&v_max,&u_ave);
			error_Per=Comput_Perm(psi,u,Per_l,Per_g,PerDir);
			S_l=Comput_Saturation(psi,Solid);
			
			
			if (rank==0)
			{

			ofstream fin(FileName,ios::app);
			fin<<"The"<<n<<"th computation result:"<<endl;
		//=============================================================================================
			//fin<<"The permiability is: "<<Permia[0]*reso*reso*100<<", "<<Permia[1]*reso*reso*100<<", "<<Permia[2]*reso*reso*100<<endl;
			//fin<<"The relative error of permiability computing is: "<<error_perm<<endl;
		//==============================================================================================

		//==============================================================================================
			//cout<<"The Density of point(NX/2,NY/2,NZ/2) is: "<<setprecision(6)
			//	<<rho[int((NX+1)/para_size/2)][NY/2][NZ/2]<<endl;
			Re_l=u_ave*(NY+1)/niu_l;Re_g=u_ave*(NY+1)/niu_g;
			fin<<"The Maximum velocity is: "<<setprecision(6)<<v_max<<"   Re_l="<<Re_l<<"   Re_g="<<Re_g<<endl;
			
		//===============================================================================================
			fin<<"The max relative error of velocity is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			fin<<"The relative permeability of component 1 is "<<Per_l[0]*reso*reso*1000/Permeability<<", "<<Per_l[1]*reso*reso*1000/Permeability<<", "<<Per_l[2]*reso*reso*1000/Permeability<<endl;
			fin<<"The relative permeability of component 2 is "<<Per_g[0]*reso*reso*1000/Permeability<<", "<<Per_g[1]*reso*reso*1000/Permeability<<", "<<Per_g[2]*reso*reso*1000/Permeability<<endl;
			fin<<"Satuation of Component 1: "<<S_l<<", "<<"The satuation of Component 2: "<<1-S_l<<endl;
			fin<<"The relative error of permiability computing is: "<<error_Per<<endl;
			fin<<endl;
			fin.close();


		if (wr_per==1)
			{
			ofstream finfs(FileName2,ios::app);
			switch(PerDir)
				{
				case 1:
				finfs<<Per_l[0]*reso*reso*1000/Permeability<<endl;break;
				case 2:
				finfs<<Per_l[1]*reso*reso*1000/Permeability<<endl;break;
				case 3:
				finfs<<Per_l[2]*reso*reso*1000/Permeability<<endl;break;
				default:
				finfs<<Per_l[0]*reso*reso*1000/Permeability<<endl;break;
				}
			finfs.close();
			
			ofstream finfs2(FileName4,ios::app);
			switch(PerDir)
				{
				case 1:
				finfs2<<Per_g[0]*reso*reso*1000/Permeability<<endl;break;
				case 2:
				finfs2<<Per_g[1]*reso*reso*1000/Permeability<<endl;break;
				case 3:
				finfs2<<Per_g[2]*reso*reso*1000/Permeability<<endl;break;
				default:
				finfs2<<Per_g[0]*reso*reso*1000/Permeability<<endl;break;
				}
			finfs2.close();
			}

			
//============================================================================================================

			cout<<"The"<<n<<"th computation result:"<<endl;
			cout<<"The Density of point(NX/2,NY/2,NZ/2) is: "<<setprecision(6)
				<<rho[Solid[((NX+1)/para_size/2)][NY/2][NZ/2]]<<endl;
			
			cout<<"The Maximum velocity is: "<<setprecision(6)<<v_max<<"   Re_l="<<Re_l<<"   Re_g="<<Re_g<<endl;
			
			cout<<"The max relative error of uv is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			cout<<"The relative permeability of component 1 is "<<Per_l[0]*reso*reso*1000/Permeability<<", "<<Per_l[1]*reso*reso*1000/Permeability<<", "<<Per_l[2]*reso*reso*1000/Permeability<<endl;
			cout<<"The relative permeability of component 2 is "<<Per_g[0]*reso*reso*1000/Permeability<<", "<<Per_g[1]*reso*reso*1000/Permeability<<", "<<Per_g[2]*reso*reso*1000/Permeability<<endl;
			cout<<"Satuation of Component 1: "<<S_l<<", "<<"The satuation of Component 2: "<<1-S_l<<endl;
			cout<<"The relative error of permiability computing is: "<<error_Per<<endl;
			cout<<endl;
                        }


			
			if ((freDe>=0) and (n%freDe==0))
				if (Out_Mode==1)
					output_density(n,rho,mirX,mirY,mirZ,mir,Solid);
				else
					output_density_b(n,rho,mirX,mirY,mirZ,mir,Solid);

			if ((freVe>=0) and (n%freVe==0))
				if (Out_Mode==1)
					output_velocity(n,rho,u,mirX,mirY,mirZ,mir,Solid);
				else
					output_velocity_b(n,rho,u,mirX,mirY,mirZ,mir,Solid);
				
			if ((freVe>=0) and (n%freVe==0))
			        if (Out_Mode==1)
			                output_psi(n,psi,mirX,mirY,mirZ,mir,Solid);
			        else
			                output_psi_b(n,psi,mirX,mirY,mirZ,mir,Solid);

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
		delete [] g[i];
		delete [] Fg[i];
		delete [] psi;
		}
		
	delete [] g;
	delete [] Fg;
	delete [] psi;
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



void init_Sparse(int*** Solids, int*** Solid, double***Psis, double*** Psi_local,int* Sl,int* Sr)
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




for(int i=0;i<nx_l;i++)	
	for(int j=0;j<=NY;j++)
		for(int k=0;k<=NZ;k++)
		{
			
		        Psi_local[i][j][k]=Psis[int((s_c+i-(s_c+i)%Zoom)/Zoom)][int((j-j%Zoom)/Zoom)][int((k-k%Zoom)/Zoom)];
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

//****************************
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
//*******************************









void Read_Rock(int*** Solids,double*** Psis,double* porosity,char poreFileName[128], char psiFileName[128])
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
double* Psi_Int;

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
Psi_Int = new double[nx*ny*nz];


	

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
		
if (rank==0)
        
{
FILE *ftest2;
	ifstream fin2;
	
	ftest2 = fopen(psiFileName, "r");

	if(ftest2 == NULL)
	{
		cout << "\n The concentration file (" << psiFileName <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest2);

	fin2.open(psiFileName);


	
	// Reading pore geometry
	for(k=0 ; k<nz_a ; k++)
	for(j=0 ; j<ny_a ; j++)
	for(i=0 ; i<nx_a ; i++)
	
	{
		while(true)
		{	
			fin2 >> pore;
			if( pore == -1.0 || pore == 1.0) break;
		}
		if ((pore == -1.0) && (i<nx0) && (j<ny0) && (k<nz0))	Psi_Int[i*ny*nz+j*nz+k] = -1.0;
		//else			Solid_Int[i][j][k] = 1;
		if ((pore == 1.0) && (i<nx0) && (j<ny0) && (k<nz0))	Psi_Int[i*ny*nz+j*nz+k] = 1.0;
	}
	fin2.close();

	// Mirroring the concentration
	if(mirX==1){
		for(i=nx0 ; i<nx ; i++)
		for(j=0   ; j<ny ; j++)
		for(k=0   ; k<nz ; k++)
				Psi_Int[i*ny*nz+j*nz+k] = Psi_Int[(nx-i-1)*ny*nz+j*nz+k];
	}

	if(mirY==1){
		for(j=ny0 ; j<ny ; j++)
		for(i=0   ; i<nx ; i++)
		for(k=0   ; k<nz ; k++)
				Psi_Int[i*ny*nz+j*nz+k] = Psi_Int[i*ny*nz+(ny-j-1)*nz+k];
	}

	if(mirZ==1){
		for(k=nz0 ; k<nz ; k++)
		for(j=0   ; j<ny ; j++)
		for(i=0   ; i<nx ; i++)
				Psi_Int[i*ny*nz+j*nz+k] = Psi_Int[i*ny*nz+j*nz+nz-k-1];
	}


}


	
	MPI_Bcast(Psi_Int,nx*ny*nz,MPI_DOUBLE,0,MPI_COMM_WORLD);

	
	cout<<"Concentration FILE READING COMPLETE. "<<endl;

	
	for (i=0;i<nx;i++)
		for (j=0;j<ny;j++)
			for (k=0;k<nz;k++)
			{
			   //     if (Psi_Int[i*(ny)*(nz)+j*(nz)+k]>0)
			      //          cout<<i<<" "<<j<<" "<<k<<endl;
			Psis[i][j][k]=Psi_Int[i*(ny)*(nz)+j*(nz)+k];
			}	
		
		
		

//cout<<"asdfasdfasdfasdfa"<<endl;
MPI_Barrier(MPI_COMM_WORLD);

	delete [] Solid_Int;
	delete [] Psi_Int;


}




void init(double* rho, double** u, double** f,double** g, double* psi,  double* forcex,double* forcey, double* forcez, double*** Psi_local, int* SupInv)
{	
    
	double usqr,vsqr;

	
	rho0=1.0;dt=1.0/Zoom;dx=1.0/Zoom;
 	uMax=0.0;

	double pr; //raduis of the obstacles
	double s_other=8*(2-s_v)/(8-s_v);
	double u_tmp[3];
	
	S[0]=0;
	S[1]=1;
	S[2]=1;
	S[3]=0;
	S[4]=1;
	S[5]=0;
	S[6]=1;
	S[7]=0;
	S[8]=1;
	S[9]=s_v;
	S[10]=1;
	S[11]=s_v;
	S[12]=1;
	S[13]=s_v;
	S[14]=s_v;
	S[15]=s_v;
	S[16]=1;
	S[17]=1;
	S[18]=1;
	
	
/*	
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
*/


	
	for (int i=1;i<=Count;i++)	
			
		{
			u[i][0]=inivx;
			u[i][1]=inivy;
			u[i][2]=inivz;
			u_tmp[0]=u[i][0];
			u_tmp[1]=u[i][1];
			u_tmp[2]=u[i][2];
			psi[i]=Psi_local[(int)(SupInv[i]/((NY+1)*(NZ+1)))][(int)((SupInv[i]%((NY+1)*(NZ+1)))/(NZ+1))][SupInv[i]%(NZ+1)];
			rho[i]=1.0;
		
			

			//***********************************************************************

			forcex[i]=gx;
			forcey[i]=gy;
			forcez[i]=gz;


			s_v=niu_g+(psi[i]+1.0)/2.0*(niu_l-niu_g);
			s_v=1.0/(3*s_v/dt+0.5);
			
			S[9]=s_v;S[11]=s_v;S[13]=s_v;S[14]=s_v;S[15]=s_v;
						

	}

	

	 	
}



void periodic_streaming(double** f,double** F,double** g, double** Fg, int* SupInv,int*** Solid,int* Sl,int* Sr,double* rho,double** u)
{
	MPI_Status status[4] ;
	MPI_Request request[4];

	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	int ip,jp,kp,i,j,k,mpi_test;
	
	
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
			
			
			
for(k=0;k<19;k++)
			{
			for(i=1;i<=Gcl[rank];i++)
				sendl[(i-1)*19+k]=f[i][k];
			for(j=Count-Gcr[rank]+1;j<=Count;j++)
				sendr[(j-(Count-Gcr[rank]+1))*19+k]=f[j][k];
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
		
		
for(k=0;k<19;k++)
			{
			for(i=1;i<=Gcl[rank];i++)
				sendl[(i-1)*19+k]=g[i][k];
			for(j=Count-Gcr[rank]+1;j<=Count;j++)
				sendr[(j-(Count-Gcr[rank]+1))*19+k]=g[j][k];
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
						Fg[ci][lm]=recvl[(Sl[jp*(NZ+1)+kp]-1)*19+lm];
					else
						Fg[ci][lm]=g[ci][LR[lm]];
					
						
					
					
				if (ip>=nx_l)
					if (Sr[jp*(NZ+1)+kp]>0)
						Fg[ci][lm]=recvr[(Sr[jp*(NZ+1)+kp]-1)*19+lm];
					else
						Fg[ci][lm]=g[ci][LR[lm]];

				if ((ip>=0) and (ip<nx_l)) 
					if (Solid[ip][jp][kp]>0)
						Fg[ci][lm]=g[Solid[ip][jp][kp]][lm];
					else
						Fg[ci][lm]=g[ci][LR[lm]];
	
					
			}
		}
		
	delete [] Gcl;
	delete [] Gcr;

	delete [] sendl;
	delete [] sendr;
	delete [] recvl;
	delete [] recvr;
			
							
						

}


void periodic_streaming_Speed(double** f,double** F,double** g, double** Fg, int* SupInv,int*** Solid,int* Sl,int* Sr,double* rho,double** u)
{
	MPI_Status status[8] ;
	MPI_Request request[8];


	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	int ip,jp,kp,i,j,k,mpi_test;
	
	
	int* Gcl = new int[mpi_size];
	int* Gcr = new int[mpi_size];

	
	MPI_Gather(&cl,1,MPI_INT,Gcl,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(Gcl,mpi_size,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Gather(&cr,1,MPI_INT,Gcr,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(Gcr,mpi_size,MPI_INT,0,MPI_COMM_WORLD);


	double* recvl;
	double* recvr;
	double* recvl_g;
	double* recvr_g;

	double* sendl = new double[Gcl[rank]*19];
	double* sendr = new double[Gcr[rank]*19];
	double* sendl_g = new double[Gcl[rank]*19];
	double* sendr_g = new double[Gcr[rank]*19];
	
	
	
	if (rank==0)
		{
		recvl = new double[Gcr[mpi_size-1]*19];
		recvr = new double[Gcl[rank+1]*19];
		recvl_g = new double[Gcr[mpi_size-1]*19];
		recvr_g = new double[Gcl[rank+1]*19];
		}	
		else
		if (rank==mpi_size-1)
			{
			recvl = new double[Gcr[rank-1]*19];
			recvr = new double[Gcl[0]*19];
			recvl_g = new double[Gcr[rank-1]*19];
			recvr_g = new double[Gcl[0]*19];
			}
			else
			{
			recvl = new double[Gcr[rank-1]*19];
			recvr = new double[Gcl[rank+1]*19];
			recvl_g = new double[Gcr[rank-1]*19];
			recvr_g = new double[Gcl[rank+1]*19];
			}
			
			
			
for(k=0;k<19;k++)
			{
			for(i=1;i<=Gcl[rank];i++)
			        {
			        sendl[(i-1)*19+k]=f[i][k];
			        sendl_g[(i-1)*19+k]=f[i][k];
			        }
			for(j=Count-Gcr[rank]+1;j<=Count;j++)
			        {
				sendr[(j-(Count-Gcr[rank]+1))*19+k]=f[j][k];
				sendr_g[(j-(Count-Gcr[rank]+1))*19+k]=f[j][k];
				}
			}
MPI_Barrier(MPI_COMM_WORLD);


if (rank==0)
		{
		
		MPI_Isend(sendr, Gcr[0]*19, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(sendl, Gcl[0]*19, MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(recvr , Gcl[1]*19, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(recvl, Gcr[mpi_size-1]*19, MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );
		MPI_Isend(sendr_g, Gcr[0]*19, MPI_DOUBLE, rank+1, rank*2+1+10000, MPI_COMM_WORLD,&request[4]);
      		MPI_Isend(sendl_g, Gcl[0]*19, MPI_DOUBLE, mpi_size-1, rank*2+10000, MPI_COMM_WORLD,&request[5]);
		MPI_Irecv(recvr_g , Gcl[1]*19, MPI_DOUBLE, rank+1, (rank+1)*2+10000, MPI_COMM_WORLD,&request[6]);		
      		MPI_Irecv(recvl_g, Gcr[mpi_size-1]*19, MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1+10000, MPI_COMM_WORLD,&request[7] );
		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(sendl, Gcl[rank]*19, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr, Gcr[rank]*19, MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl, Gcr[rank-1]*19, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(recvr, Gcl[0]*19, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request[3]);
			MPI_Isend(sendl_g, Gcl[rank]*19, MPI_DOUBLE, rank-1, rank*2+10000, MPI_COMM_WORLD,&request[4]);
      			MPI_Isend(sendr_g, Gcr[rank]*19, MPI_DOUBLE, 0, rank*2+1+10000, MPI_COMM_WORLD,&request[5]);
			MPI_Irecv(recvl_g, Gcr[rank-1]*19, MPI_DOUBLE, rank-1, (rank-1)*2+1+10000, MPI_COMM_WORLD,&request[6] );
      			MPI_Irecv(recvr_g, Gcl[0]*19, MPI_DOUBLE, 0, 0+10000, MPI_COMM_WORLD,&request[7]);
			}
			else
			{

			MPI_Isend(sendl, Gcl[rank]*19, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr, Gcr[rank]*19, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl, Gcr[rank-1]*19, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(recvr, Gcl[rank+1]*19, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);
      			MPI_Isend(sendl_g, Gcl[rank]*19, MPI_DOUBLE, rank-1, rank*2+10000, MPI_COMM_WORLD,&request[4]);
      			MPI_Isend(sendr_g, Gcr[rank]*19, MPI_DOUBLE, rank+1, rank*2+1+10000, MPI_COMM_WORLD,&request[5]);
			MPI_Irecv(recvl_g, Gcr[rank-1]*19, MPI_DOUBLE, rank-1, (rank-1)*2+1+10000, MPI_COMM_WORLD,&request[6]);
      			MPI_Irecv(recvr_g, Gcl[rank+1]*19, MPI_DOUBLE, rank+1, (rank+1)*2+10000, MPI_COMM_WORLD,&request[7]);
			
			};

	
	MPI_Waitall(8,request, status);

	MPI_Testall(8,request,&mpi_test,status);



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
					{
						F[ci][lm]=recvl[(Sl[jp*(NZ+1)+kp]-1)*19+lm];
						Fg[ci][lm]=recvl_g[(Sl[jp*(NZ+1)+kp]-1)*19+lm];
					}
					else
					{
						F[ci][lm]=f[ci][LR[lm]];
						Fg[ci][lm]=g[ci][LR[lm]];
					}
						
					
					
				if (ip>=nx_l)
					if (Sr[jp*(NZ+1)+kp]>0)
					{
						F[ci][lm]=recvr[(Sr[jp*(NZ+1)+kp]-1)*19+lm];
						Fg[ci][lm]=recvr_g[(Sr[jp*(NZ+1)+kp]-1)*19+lm];
					}
					else
					{
						F[ci][lm]=f[ci][LR[lm]];
						Fg[ci][lm]=g[ci][LR[lm]];
					}

				if ((ip>=0) and (ip<nx_l)) 
					if (Solid[ip][jp][kp]>0)
					{
						F[ci][lm]=f[Solid[ip][jp][kp]][lm];
						Fg[ci][lm]=g[Solid[ip][jp][kp]][lm];
					}
					else
					{
						F[ci][lm]=f[ci][LR[lm]];
						Fg[ci][lm]=g[ci][LR[lm]];
					}
	
					
			}
		}
		

	delete [] Gcl;
	delete [] Gcr;

	delete [] sendl;
	delete [] sendr;
	delete [] recvl;
	delete [] recvr;
	delete [] sendl_g;
	delete [] sendr_g;
	delete [] recvl_g;
	delete [] recvr_g;
							
						

}










void collision(double* rho,double* psi,double** u,double** f, double** g,double* forcex, double* forcey, double* forcez, int* Sl,int* Sr,int* SupInv,int*** Solid)
{
        MPI_Status status[8] ;
	MPI_Request request[8];
	int mpi_test;

	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

double lm0,lm1,s_v;
double A,B,C,D,E;
double F_hat[19],GuoF[19],f_eq[19];
double m_l[19],meq[19];;
double feq[19],feq_g[19],Fi[19];
int i,j,m;

A=1.0/6.0;
B=1.0/12.0;
C=1.0/3.0;
D=1.0/6.0;
E=-6*C-12*D;

//=========================CONSTANT DEFINITION FOR EQUILIBRIUM COMPUTING=====================================
double px[19]={0.0,A,-A,0.0,0.0,0.0,0.0,B,-B,B,-B,0.0,0.0,0.0,0.0,B,-B,B,-B};
double pz[19]={0.0,0.0,0.0,0.0,0.0,A,-A,0.0,0.0,0.0,0.0,B,B,-B,-B,B,B,-B,-B};
double py[19]={0.0,0.0,0.0,A,-A,0.0,0.0,B,B,-B,-B,B,-B,B,-B,0.0,0.0,0.0,0.0};
double lap[19]={E,C,C,C,C,C,C,D,D,D,D,D,D,D,D,D,D,D,D};
double delta[3][3]={{1.0,0.0,0.0},{0.0,1.0,0.0},{0.0,0.0,1.0}};
double ome[19][6]={{0.0,0.0,0.0,0.0,0.0,0.0},
		{5.0/12.0,-1.0/3.0,-1.0/3.0,0.0,0.0,0.0},{5.0/12.0,-1.0/3.0,-1.0/3.0,0.0,0.0,0.0},
		{-1.0/3.0,5.0/12.0,-1.0/3.0,0.0,0.0,0.0},{-1.0/3.0,5.0/12.0,-1.0/3.0,0.0,0.0,0.0},
		{-1.0/3.0,-1.0/3.0,5.0/12.0,0.0,0.0,0.0},{-1.0/3.0,-1.0/3.0,5.0/12.0,0.0,0.0,0.0},
		{-1.0/24.0,-1.0/24.0,1.0/12.0,1.0/4.0,0.0,0.0},{-1.0/24.0,-1.0/24.0,1.0/12.0,-1.0/4.0,0.0,0.0},
		{-1.0/24.0,-1.0/24.0,1.0/12.0,-1.0/4.0,0.0,0.0},{-1.0/24.0,-1.0/24.0,1.0/12.0,1.0/4.0,0.0,0.0},
		{1.0/12.0,-1.0/24.0,-1.0/24.0,0.0,1.0/4.0,0.0},{1.0/12.0,-1.0/24.0,-1.0/24.0,0.0,-1.0/4.0,0.0},
		{1.0/12.0,-1.0/24.0,-1.0/24.0,0.0,-1.0/4.0,0.0},{1.0/12.0,-1.0/24.0,-1.0/24.0,0.0,1.0/4.0,0.0},
		{-1.0/24.0,1.0/12.0,-1.0/24.0,0.0,0.0,1.0/4.0},{-1.0/24.0,1.0/12.0,-1.0/24.0,0.0,0.0,-1.0/4.0},
		{-1.0/24.0,1.0/12.0,-1.0/24.0,0.0,0.0,-1.0/4.0},{-1.0/24.0,1.0/12.0,-1.0/24.0,0.0,0.0,1.0/4.0}};


double p0;
double T=0.56; 		//TEMPRETURE FOR NONIDEAL SYSTEM
double a=9.0/49.0;  	//PARAMETER FOR DEFINITION OF EQUATION OF STATE
double b=2.0/21.0;	//PARAMETER FOR DEFINITION OF EQUATION OF STATE
double par[3],par_r[3],delxy,delxy_r;
//================================================================





	int* Gcl = new int[mpi_size];
	int* Gcr = new int[mpi_size];

	
	MPI_Gather(&cl,1,MPI_INT,Gcl,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(Gcl,mpi_size,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Gather(&cr,1,MPI_INT,Gcr,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(Gcr,mpi_size,MPI_INT,0,MPI_COMM_WORLD);

	double* recvl_r;
	double* recvr_r;
	double* recvl_psi;
	double* recvr_psi;

	double* sendl_r = new double[Gcl[rank]];
	double* sendr_r = new double[Gcr[rank]];
	double* sendl_psi = new double[Gcl[rank]];
	double* sendr_psi = new double[Gcr[rank]];
	
	
	
	if (rank==0)
		{
		recvl_r= new double[Gcr[mpi_size-1]];
		recvr_r= new double[Gcl[rank+1]];
		recvl_psi = new double[Gcr[mpi_size-1]];
		recvr_psi = new double[Gcl[rank+1]];
		}	
		else
		if (rank==mpi_size-1)
			{
			recvl_r = new double[Gcr[rank-1]];
			recvr_r = new double[Gcl[0]];
			recvl_psi = new double[Gcr[rank-1]];
			recvr_psi = new double[Gcl[0]];
			}
			else
			{
			recvl_r = new double[Gcr[rank-1]];
			recvr_r = new double[Gcl[rank+1]];
			recvl_psi = new double[Gcr[rank-1]];
			recvr_psi = new double[Gcl[rank+1]];
			}
			
			

			for(i=1;i<=Gcl[rank];i++)
			        {
			        sendl_r[i-1]=rho[i];
			        sendl_psi[i-1]=psi[i];
			        }
			for(j=Count-Gcr[rank]+1;j<=Count;j++)
			        {
				sendr_r[j-(Count-Gcr[rank]+1)]=rho[j];
				sendr_psi[j-(Count-Gcr[rank]+1)]=psi[j];
				}
			
MPI_Barrier(MPI_COMM_WORLD);

//cout<<"@@@@@@@@@@@   "<<n<<endl;
if (rank==0)
		{
		
		MPI_Isend(sendr_r, Gcr[0], MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(sendl_r, Gcl[0], MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(recvr_r, Gcl[1], MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(recvl_r, Gcr[mpi_size-1], MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );
		MPI_Isend(sendr_psi, Gcr[0], MPI_DOUBLE, rank+1, rank*2+1+10000, MPI_COMM_WORLD,&request[4]);
      		MPI_Isend(sendl_psi, Gcl[0], MPI_DOUBLE, mpi_size-1, rank*2+10000, MPI_COMM_WORLD,&request[5]);
		MPI_Irecv(recvr_psi, Gcl[1], MPI_DOUBLE, rank+1, (rank+1)*2+10000, MPI_COMM_WORLD,&request[6]);		
      		MPI_Irecv(recvl_psi, Gcr[mpi_size-1], MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1+10000, MPI_COMM_WORLD,&request[7] );
		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(sendl_r, Gcl[rank], MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr_r, Gcr[rank], MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl_r, Gcr[rank-1], MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(recvr_r, Gcl[0], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request[3]);
			MPI_Isend(sendl_psi, Gcl[rank], MPI_DOUBLE, rank-1, rank*2+10000, MPI_COMM_WORLD,&request[4]);
      			MPI_Isend(sendr_psi, Gcr[rank], MPI_DOUBLE, 0, rank*2+1+10000, MPI_COMM_WORLD,&request[5]);
			MPI_Irecv(recvl_psi, Gcr[rank-1], MPI_DOUBLE, rank-1, (rank-1)*2+1+10000, MPI_COMM_WORLD,&request[6] );
      			MPI_Irecv(recvr_psi, Gcl[0], MPI_DOUBLE, 0, 0+10000, MPI_COMM_WORLD,&request[7]);
			}
			else
			{

			MPI_Isend(sendl_r, Gcl[rank], MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr_r, Gcr[rank], MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl_r, Gcr[rank-1], MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(recvr_r, Gcl[rank+1], MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);
      			MPI_Isend(sendl_psi, Gcl[rank], MPI_DOUBLE, rank-1, rank*2+10000, MPI_COMM_WORLD,&request[4]);
      			MPI_Isend(sendr_psi, Gcr[rank], MPI_DOUBLE, rank+1, rank*2+1+10000, MPI_COMM_WORLD,&request[5]);
			MPI_Irecv(recvl_psi, Gcr[rank-1], MPI_DOUBLE, rank-1, (rank-1)*2+1+10000, MPI_COMM_WORLD,&request[6]);
      			MPI_Irecv(recvr_psi, Gcl[rank+1], MPI_DOUBLE, rank+1, (rank+1)*2+10000, MPI_COMM_WORLD,&request[7]);
			
			};

	
	MPI_Waitall(8,request, status);

	MPI_Testall(8,request,&mpi_test,status);
			
			
double term1,term2,term3,tempvar;
double lambda,miu,niu_temp;  	//NEED TO BE DETERMINED  *******
int interi,interj,interk;			
			
			

	for(int ci=1;ci<=Count;ci++)	
	

		{	
				i=(int)(SupInv[ci]/((NY+1)*(NZ+1)));
				j=(int)((SupInv[ci]%((NY+1)*(NZ+1)))/(NZ+1));
				m=(int)(SupInv[ci]%(NZ+1));   
				

//=====================================BULK PRESSURE DEFINITION=====================================================

//==============================NUMERICAL DERIVATIVES COMPUTING=(BINARY SYSTEM)============================
		par[1]=0;par[0]=0;par[2]=0;delxy=0;par_r[1]=0;par_r[0]=0;par_r[2]=0;delxy_r=0;
		for (int tmpi=0;tmpi<19;tmpi++)
		{

			//-------------------PERIODIC BOUNDARY CONDITION---------------------------
		if (in_BC>0)
		{	     
			interi=i+e[tmpi][0];
			if (((pre_xn-1)* (vel_xn-1)==0) and (rank==0) and (interi<0))
			        interi=1;
			if (((pre_xp-1)*(vel_xp-1)==0) and (rank==mpi_size-1) and (interi>=nx_l))
			        interi=nx_l-2;
			
			interj=j+e[tmpi][1];
			if ((pre_yn-1)*(vel_yn-1)==0)
			        {if (interj<0) {interj=1;}}
			else
			        {if (interj<0) {interj=NY;}}
			 
			if ((pre_yp-1)*(vel_yp-1)==0)
			        {if (interj>NY) {interj=NY-1;}}
			else
			        {if (interj>NY) {interj=0;}}
			
			
			interk=m+e[tmpi][2];
			if ((pre_zn-1)*(vel_zn-1)==0)
			        {if (interk<0) {interk=1;}}
			else
			        {if (interk<0) {interk=NZ;}}
			
			if ((pre_zp-1)*(vel_zp-1)==0)
			        {if (interk>NZ) {interk=NZ-1;}}
			else
			        {if (interk>NZ) {interk=0;}}
		}
		else
		        {
		        interi=i+e[tmpi][0];
			interj=j+e[tmpi][1];if (interj<0) {interj=NY;}; if (interj>NY) {interj=0;};
			interk=m+e[tmpi][2];if (interk<0) {interk=NZ;}; if (interk>NZ) {interk=0;};          
		                
		        }
		
			//-------------------------------------------------------------------------
			
			
			
			
	/*
		if (ip>=nx_l)
					if (Sr[jp*(NZ+1)+kp]>0)
					{
						F[ci][lm]=recvr[(Sr[jp*(NZ+1)+kp]-1)*19+lm];
						Fg[ci][lm]=recvr_g[(Sr[jp*(NZ+1)+kp]-1)*19+lm];
					}
					else
					{
						F[ci][lm]=f[ci][LR[lm]];
						Fg[ci][lm]=g[ci][LR[lm]];
					}

						
	*/
	
	
			
			if (interi<0)
			{
			        if (Sl[interj*(NZ+1)+interk]>0)
			        {
			                par_r[0]+=px[tmpi]*recvl_r[Sl[interj*(NZ+1)+interk]-1];
			                par_r[1]+=py[tmpi]*recvl_r[Sl[interj*(NZ+1)+interk]-1];
			                par_r[2]+=pz[tmpi]*recvl_r[Sl[interj*(NZ+1)+interk]-1];

			                par[0]+=px[tmpi]*recvl_psi[Sl[interj*(NZ+1)+interk]-1];
			                par[1]+=py[tmpi]*recvl_psi[Sl[interj*(NZ+1)+interk]-1];
			                par[2]+=pz[tmpi]*recvl_psi[Sl[interj*(NZ+1)+interk]-1];

			                delxy+=lap[tmpi]*recvl_psi[Sl[interj*(NZ+1)+interk]-1];
			                delxy_r+=lap[tmpi]*recvl_r[Sl[interj*(NZ+1)+interk]-1];
			        }
			else
			        {
			                par_r[0]+=px[tmpi]*(1.0);
			                par_r[1]+=py[tmpi]*(1.0);
			                par_r[2]+=pz[tmpi]*(1.0);

			                par[0]+=px[tmpi]*(psi[ci]-2*ContactAngle_parameter);
			                par[1]+=py[tmpi]*(psi[ci]-2*ContactAngle_parameter);
			                par[2]+=pz[tmpi]*(psi[ci]-2*ContactAngle_parameter);

			                delxy+=lap[tmpi]*(psi[ci]-2*ContactAngle_parameter);
			                delxy_r+=lap[tmpi]*1.0;     
			        }
			}
			
			
			if (interi>=nx_l)
			{
			        if (Sr[interj*(NZ+1)+interk]>0)
			        {
			                par_r[0]+=px[tmpi]*recvr_r[Sr[interj*(NZ+1)+interk]-1];
			                par_r[1]+=py[tmpi]*recvr_r[Sr[interj*(NZ+1)+interk]-1];
			                par_r[2]+=pz[tmpi]*recvr_r[Sr[interj*(NZ+1)+interk]-1];

			                par[0]+=px[tmpi]*recvr_psi[Sr[interj*(NZ+1)+interk]-1];
			                par[1]+=py[tmpi]*recvr_psi[Sr[interj*(NZ+1)+interk]-1];
			                par[2]+=pz[tmpi]*recvr_psi[Sr[interj*(NZ+1)+interk]-1];

			                delxy+=lap[tmpi]*recvr_psi[Sr[interj*(NZ+1)+interk]-1];
			                delxy_r+=lap[tmpi]*recvr_r[Sr[interj*(NZ+1)+interk]-1];
			        }
			        else
			                
			                {
			                par_r[0]+=px[tmpi]*1.0;
			                par_r[1]+=py[tmpi]*1.0;
			                par_r[2]+=pz[tmpi]*1.0;

			                par[0]+=px[tmpi]*(psi[ci]-2*ContactAngle_parameter);
			                par[1]+=py[tmpi]*(psi[ci]-2*ContactAngle_parameter);
			                par[2]+=pz[tmpi]*(psi[ci]-2*ContactAngle_parameter);

			                delxy+=lap[tmpi]*(psi[ci]-2*ContactAngle_parameter);
			                delxy_r+=lap[tmpi]*1.0;
      			                }
      			}
      			
      			
      			
			if ((interi>=0) and (interi<nx_l))
			{
			        if (Solid[interi][interj][interk]>0)
			        {
			                par_r[0]+=px[tmpi]*rho[Solid[interi][interj][interk]];
			                par_r[1]+=py[tmpi]*rho[Solid[interi][interj][interk]];
			                par_r[2]+=pz[tmpi]*rho[Solid[interi][interj][interk]];

			                par[0]+=px[tmpi]*psi[Solid[interi][interj][interk]];        
			                par[1]+=py[tmpi]*psi[Solid[interi][interj][interk]];
			                par[2]+=pz[tmpi]*psi[Solid[interi][interj][interk]];

			                delxy+=lap[tmpi]*psi[Solid[interi][interj][interk]];
			                delxy_r+=lap[tmpi]*rho[Solid[interi][interj][interk]];
			        }        
			        else
			        {
			                
			                par_r[0]+=px[tmpi]*1.0;
			                par_r[1]+=py[tmpi]*1.0;
			                par_r[2]+=pz[tmpi]*1.0;

			                par[0]+=px[tmpi]*(psi[ci]-2*ContactAngle_parameter);        
			                par[1]+=py[tmpi]*(psi[ci]-2*ContactAngle_parameter);
			                par[2]+=pz[tmpi]*(psi[ci]-2*ContactAngle_parameter);

			                delxy+=lap[tmpi]*(psi[ci]-2*ContactAngle_parameter);
			                delxy_r+=lap[tmpi]*1.0;  
			                
			        }
			
			
			}
			}
	tempvar=psi[ci]*psi[ci];
	p0=rho[ci]/3.0+ca*(-0.5*tempvar+0.75*tempvar*tempvar);
	lambda=0.0;
	miu=ca*(-psi[ci]+tempvar*psi[ci])-kappa*delxy;
			
	
	for (int k=1;k<19;k++)
		{
		feq[k]=0;feq_g[k]=0;
		
		feq[k]+=p0-kappa*psi[ci]*delxy+(e[k][0]*u[ci][0]+e[k][1]*u[ci][1]+e[k][2]*u[ci][2])*rho[ci];
		
		term1=0;term2=0;term3=0;
		
		for(int ii=0;ii<=2;ii++)
			for (int jj=0;jj<=2;jj++)
			{
		term1+=(e[k][ii]*e[k][jj]-(1.0/3.0)*(delta[ii][jj]))*(rho[ci]*u[ci][ii]*u[ci][jj]+lambda*(u[ci][ii]*par_r[jj]+u[ci][jj]*par_r[ii]+delta[ii][jj]*(u[ci][0]*par_r[0]+u[ci][1]*par_r[1]+u[ci][2]*par_r[2])));

			term2+=1.5*(e[k][ii]*e[k][jj]-1.0/3.0*delta[ii][jj])*psi[ci]*u[ci][ii]*u[ci][jj];
		
			}
		feq_g[k]=w[k]*3.0*(2.0*CM*miu+e[k][0]*psi[ci]*u[ci][0]+e[k][1]*psi[ci]*u[ci][1]+e[k][2]*psi[ci]*u[ci][2]+term2);

		term1*=3.0/2.0;
		term2=kappa*(ome[k][0]*par[0]*par[0]+ome[k][1]*par[1]*par[1]+ome[k][2]*par[2]*par[2]+ome[k][3]*par[0]*par[1]+ome[k][4]*par[1]*par[2]+ome[k][5]*par[2]*par[0]);

		
		feq[k]=(term1+feq[k])*w[k]*3.0+term2;	
		
		}

	feq[0]=rho[ci]-(feq[1]+feq[2]+feq[3]+feq[4]+feq[5]+feq[6]+feq[7]+feq[8]+feq[9]+feq[10]+feq[11]+feq[12]+feq[13]+feq[14]+feq[15]+feq[16]+feq[17]+feq[18]);

	feq_g[0]=psi[ci]-(feq_g[1]+feq_g[2]+feq_g[3]+feq_g[4]+feq_g[5]+feq_g[6]+feq_g[7]+feq_g[8]+feq_g[9]+feq_g[10]+feq_g[11]+feq_g[12]+feq_g[13]+feq_g[14]+feq_g[15]+feq_g[16]+feq_g[17]+feq_g[18]);

	
	
//==========================================EQILIBRIUM FUNCTION COMPUTING FOR G FUNCTION===========================
//***THIS PART IS FOR G FUNCTION RELAXATION TIME =1 , FOR OTHER RELAXATION TIME,MODIFICATION SHOULD BE INVOLVED****





			//=================FORCE TERM_GUO=========================================
			for (int k=0;k<19;k++)
			{	
			lm0=((e[k][0]-u[ci][0])*forcex[ci]+(e[k][1]-u[ci][1])*forcey[ci]+(e[k][2]-u[ci][2])*forcez[ci])*3;
			lm1=(e[k][0]*u[ci][0]+e[k][1]*u[ci][1]+e[k][2]*u[ci][2])*(e[k][0]*forcex[ci]+e[k][1]*forcey[ci]+e[k][2]*forcez[ci])*9;
			GuoF[k]=w[k]*(lm0+lm1);
			//GuoF[k]=0.0;
			}

			//=====================equilibrium of moment=================================
			

			
			for (int l=0;l<19;l++)
				{
				meq[l]=0;
				for(int lm=0;lm<19;lm++)
				meq[l]+=M[l][lm]*feq[lm];				
				}

			//============================================================================

			s_v=niu_g+(psi[ci]+1.0)/2.0*(niu_l-niu_g);
			s_v=1.0/(3*s_v/dt+0.5);

			S[9]=s_v;S[11]=s_v;S[13]=s_v;S[14]=s_v;S[15]=s_v;
			
		// ==================   m=Mf matrix calculation  =============================
			// ==================   F_hat=(I-.5*S)MGuoF =====================================
				for (int mi=0; mi<19; mi++)
					{m_l[mi]=0;F_hat[mi]=0;
					for (int mj=0; mj<19; mj++)
						{
						m_l[mi]+=M[mi][mj]*f[ci][mj];
						F_hat[mi]+=M[mi][mj]*GuoF[mj];
						}
					F_hat[mi]*=(1-0.5*S[mi]);
					}	
					
					
			for (int sk=0;sk<19;sk++)
				//m[sk]=m[sk]-S[sk]*(m[sk]-meq[sk]);
				//m[sk]=m[sk]-S[sk]*(m[sk]-meq[sk])+(1-S[sk]*F_hat[sk]/2);
				m_l[sk]=m_l[sk]-S[sk]*(m_l[sk]-meq[sk])+F_hat[sk];

			//}

			// ==================   f=M_-1m matrix calculation  =============================
				for (int mi=0; mi<19; mi++)
					{f[ci][mi]=0;
					for (int mj=0; mj<19; mj++)
						f[ci][mi]+=MI[mi][mj]*m_l[mj];
					}
			//============================================================================
	
			g[ci][0]=feq_g[0];g[ci][1]=feq_g[1];g[ci][2]=feq_g[2];g[ci][3]=feq_g[3];
			g[ci][4]=feq_g[4];g[ci][5]=feq_g[5];g[ci][6]=feq_g[6];g[ci][7]=feq_g[7];
			g[ci][8]=feq_g[8];g[ci][9]=feq_g[9];g[ci][10]=feq_g[10];g[ci][11]=feq_g[11];
			g[ci][12]=feq_g[12];g[ci][13]=feq_g[13];g[ci][14]=feq_g[14];g[ci][15]=feq_g[15];
			g[ci][16]=feq_g[16];g[ci][17]=feq_g[17];g[ci][18]=feq_g[18];

			if (n==0)
			{
			f[ci][0]=feq[0];f[ci][1]=feq[1];f[ci][2]=feq[2];f[ci][3]=feq[3];
			f[ci][4]=feq[4];f[ci][5]=feq[5];f[ci][6]=feq[6];f[ci][7]=feq[7];
			f[ci][8]=feq[8];f[ci][9]=feq[9];f[ci][10]=feq[10];f[ci][11]=feq[11];
			f[ci][12]=feq[12];f[ci][13]=feq[13];f[ci][14]=feq[14];f[ci][15]=feq[15];
			f[ci][16]=feq[16];f[ci][17]=feq[17];f[ci][18]=feq[18];
			
			}
			
			if (in_BC>0)
			        {
			        if ((((pre_xn-1)*(vel_xn-1)==0) and (rank==0) and (i==0)) or
			                (((pre_xp-1)*(vel_xp-1)==0) and (rank==mpi_size-1) and (i==nx_l)) or
			                (((pre_yn-1)*(vel_yn-1)==0) and (j==0)) or
			                (((pre_yp-1)*(vel_yp-1)==0) and (j==NY)) or
			                (((pre_zn-1)*(vel_zn-1)==0) and (m==0)) or
			                (((pre_zp-1)*(vel_zp-1)==0) and (m==NZ))) 
			                {
			                 f[ci][0]=feq[0];f[ci][1]=feq[1];f[ci][2]=feq[2];f[ci][3]=feq[3];
			                 f[ci][4]=feq[4];f[ci][5]=feq[5];f[ci][6]=feq[6];f[ci][7]=feq[7];
			                 f[ci][8]=feq[8];f[ci][9]=feq[9];f[ci][10]=feq[10];f[ci][11]=feq[11];
			                 f[ci][12]=feq[12];f[ci][13]=feq[13];f[ci][14]=feq[14];f[ci][15]=feq[15];
			                 f[ci][16]=feq[16];f[ci][17]=feq[17];f[ci][18]=feq[18];      
			                }
			                
			  		        
			        
			        }
			
			
			
			
			
			
			
		}	
		
	delete [] sendl_r;
	delete [] sendr_r;
	delete [] recvl_r;
	delete [] recvr_r;
	delete [] sendl_psi;
	delete [] sendr_psi;
	delete [] recvl_psi;
	delete [] recvr_psi;	

	//delete pointes*******	
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


void comput_macro_variables(double* psi, double* rho,double** u,double** u0,double** f,double** g, double** F,double** Fg,double* forcex, double* forcey, double* forcez,int* SupInv,int*** Solid,double*** Psi_local)
{
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	int lsi;

	for(int i=1;i<=Count;i++)	
                   
			{
			
				u0[i][0]=u[i][0];
				u0[i][1]=u[i][1];
				u0[i][2]=u[i][2];
				rho[i]=0;psi[i]=0;
				u[i][0]=0;
				u[i][1]=0;
				u[i][2]=0;
	
				for(int k=0;k<19;k++)
					{
					
					f[i][k]=F[i][k];
					g[i][k]=Fg[i][k];
					rho[i]+=f[i][k];
					psi[i]+=g[i][k];
					u[i][0]+=e[k][0]*f[i][k];
					u[i][1]+=e[k][1]*f[i][k];
					u[i][2]+=e[k][2]*f[i][k];
					}
				

				u[i][0]=(u[i][0]+dt*forcex[i]/2)/rho[i];
				u[i][1]=(u[i][1]+dt*forcey[i]/2)/rho[i];
				u[i][2]=(u[i][2]+dt*forcez[i]/2)/rho[i];
				
				
		
			}
if (in_BC>0)
{
if (((pre_xn-1)*(vel_xn-1)==0) and (rank==0))
for (int j=0;j<=NY;j++)
        for (int k=0;k<=NZ;k++)
        { 
                lsi=Solid[0][j][k];
                if (lsi>0)
                        if (pre_xn==1)
                        {
                                psi[lsi]=Psi_local[0][j][k];
                                rho[lsi]=p_xn;
                                if (Solid[1][j][k]>0)
                                {
                                        u[lsi][0]=u[Solid[1][j][k]][0];
                                        u[lsi][1]=u[Solid[1][j][k]][1];
                                        u[lsi][2]=u[Solid[1][j][k]][2];
                                }
                                else
                                {u[lsi][0]=0.0;u[lsi][1]=0.0;u[lsi][2]=0.0;}
                        }
                        else
                        if (vel_xn==1)
                        {        
                                psi[lsi]=Psi_local[0][j][k];
                                u[lsi][0]=v_xn;u[lsi][1]=0.0;u[lsi][2]=0.0;
                                //if (Solid[1][j][k]>0)
                                //        rho[lsi]=rho[Solid[1][j][k]];
                                //else
                                        rho[lsi]=1.0;
                        }
                }
   
                
 if (((pre_xp-1)*(vel_xp-1)==0) and (rank==mpi_size-1))
 for (int j=0;j<=NY;j++)
         for (int k=0;k<=NZ;k++)
        { 
                lsi=Solid[nx_l-1][j][k];
                if (lsi>0)
                        if (pre_xp==1)
                        {
                                psi[lsi]=Psi_local[nx_l-1][j][k];
                                rho[lsi]=p_xp;
                                if (Solid[nx_l-2][j][k]>0)
                                {
                                        u[lsi][0]=u[Solid[nx_l-1][j][k]][0];
                                        u[lsi][1]=u[Solid[nx_l-1][j][k]][1];
                                        u[lsi][2]=u[Solid[nx_l-1][j][k]][2];
                                }
                                else
                                {u[lsi][0]=0.0;u[lsi][1]=0.0;u[lsi][2]=0.0;}
                        }
                        else
                        if (vel_xp==1)
                        {        
                                psi[lsi]=Psi_local[nx_l-1][j][k];
                                u[lsi][0]=v_xp;u[lsi][1]=0.0;u[lsi][2]=0.0;
                               // if (Solid[nx_l-2][j][k]>0)
                               //         rho[lsi]=rho[Solid[nx_l-1][j][k]];
                               // else
                                        rho[lsi]=1.0;
                        }
                }                                       
 // cout<<"@@@@@@@@@  "<<pre_yn<<"   "<<pre_yp<<endl;                                      
if ((pre_yn-1)*(vel_yn-1)==0)
for (int i=0;i<nx_l;i++)
        for (int k=0;k<=NZ;k++)
        {  //cout<<"asdfasdf   "<<p_yn<<"   "<<p_yp<<endl;
                lsi=Solid[i][0][k];
                if (lsi>0)
                        if (pre_yn==1)
                        {
                                psi[lsi]=Psi_local[i][0][k];
                                rho[lsi]=p_yn;
                                if (Solid[i][1][k]>0)
                                {
                                        u[lsi][0]=u[Solid[i][1][k]][0];
                                        u[lsi][1]=u[Solid[i][1][k]][1];
                                        u[lsi][2]=u[Solid[i][1][k]][2];
                                }
                                else
                                {u[lsi][0]=0.0;u[lsi][1]=0.0;u[lsi][2]=0.0;}
                        }
                        else
                        if (vel_yn==1)
                        {        
                                psi[lsi]=Psi_local[i][0][k];
                                u[lsi][0]=0.0;u[lsi][1]=v_yn;u[lsi][2]=0.0;
                                //if (Solid[i][1][k]>0)
                                //        rho[lsi]=rho[Solid[i][1][k]];
                                //else
                                        rho[lsi]=1.0;
                        }
                }
   
                
 if ((pre_yp-1)*(vel_yp-1)==0)
 for (int i=0;i<nx_l;i++)
         for (int k=0;k<=NZ;k++)
        { 
                lsi=Solid[i][NY][k];
                if (lsi>0)
                        if (pre_yp==1)
                        {
                                psi[lsi]=Psi_local[i][NY][k];
                                rho[lsi]=p_yp;
                                if (Solid[i][NY-1][k]>0)
                                {
                                        u[lsi][0]=u[Solid[i][NY-1][k]][0];
                                        u[lsi][1]=u[Solid[i][NY-1][k]][1];
                                        u[lsi][2]=u[Solid[i][NY-1][k]][2];
                                }
                                else
                                {u[lsi][0]=0.0;u[lsi][1]=0.0;u[lsi][2]=0.0;}
                        }
                        else
                        if (vel_yp==1)
                        {        
                                psi[lsi]=Psi_local[i][NY][k];
                                u[lsi][0]=0.0;u[lsi][1]=v_yp;u[lsi][2]=0.0;
                                //if (Solid[i][NY-1][k]>0)
                                //        rho[lsi]=rho[Solid[i][NY-1][k]];
                               // else
                                        rho[lsi]=1.0;
                        }
                }                                                 
                        

if ((pre_zn-1)*(vel_zn-1)==0)
for (int i=0;i<nx_l;i++)
        for (int j=0;j<=NY;j++)
        { 
                lsi=Solid[i][j][0];
                if (lsi>0)
                        if (pre_zn==1)
                        {
                                psi[lsi]=Psi_local[i][j][0];
                                rho[lsi]=p_zn;
                                if (Solid[i][j][1]>0)
                                {
                                        u[lsi][0]=u[Solid[i][j][1]][0];
                                        u[lsi][1]=u[Solid[i][j][1]][1];
                                        u[lsi][2]=u[Solid[i][j][1]][2];
                                }
                                else
                                {u[lsi][0]=0.0;u[lsi][1]=0.0;u[lsi][2]=0.0;}
                        }
                        else
                        if (vel_zn==1)
                        {        
                                psi[lsi]=Psi_local[i][j][0];
                                u[lsi][0]=0.0;u[lsi][1]=0.0;u[lsi][2]=v_zn;
                                //if (Solid[i][j][1]>0)
                                //        rho[lsi]=rho[Solid[i][j][1]];
                                //else
                                        rho[lsi]=1.0;
                        }
                }
   
                
 if ((pre_zp-1)*(vel_zp-1)==0)
 for (int i=0;i<nx_l;i++)
         for (int j=0;j<=NY;j++)
        { 
                lsi=Solid[i][j][NZ];
                if (lsi>0)
                        if (pre_zp==1)
                        {
                                psi[lsi]=Psi_local[i][j][NZ];
                                rho[lsi]=p_zp;
                                if (Solid[i][j][NZ-1]>0)
                                {
                                        u[lsi][0]=u[Solid[i][j][NZ-1]][0];
                                        u[lsi][1]=u[Solid[i][j][NZ-1]][1];
                                        u[lsi][2]=u[Solid[i][j][NZ-1]][2];
                                }
                                else
                                {u[lsi][0]=0.0;u[lsi][1]=0.0;u[lsi][2]=0.0;}
                        }
                        else
                        if (vel_zp==1)
                        {        
                                psi[lsi]=Psi_local[i][j][NZ];
                                u[lsi][0]=0.0;u[lsi][1]=0.0;u[lsi][2]=v_zp;
                                //if (Solid[i][j][NZ-1]>0)
                                //       rho[lsi]=rho[Solid[i][j][NZ-1]];
                                //else
                                        rho[lsi]=1.0;
                        }
                }                                                 
                           
                                        
			             
}
	//if (rank==0)				
	//cout<<rho[Solid[0][3][2]]<<endl;
                     
			
	//MPI_Barrier(MPI_COMM_WORLD); 

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




void Geometry(int*** Solid)	
{	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	MPI_Status status;
	MPI_Request request;
	
	const int root_rank=0;

	int* send;
	int* rece;
	
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
	
	
	send = new int[nx_l*NY0*NZ0];
	for (int i=0;i<nx_l;i++)
		for (int j=0;j<NY0;j++)
			for (int k=0;k<NZ0;k++)
			send[i*NY0*NZ0+j*NZ0+k]=Solid[i][j][k];
		
	
	if (rank==0)
	{
	ofstream out(name.str().c_str(),ios::app);
	for(int i=0;i<nx_l;i++)
        	for(int j=0; j<NY0; j++)
			for(int k=0;k<NZ0;k++)
			if (Solid[i][j][k]<=0)
				out<<1<<endl;
			else
				out<<0<<endl;
	out.close();
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	for (int processor=1;processor<mpi_size;processor++)
	{	
		
		if (rank==0)
			rece = new int[nx_g[processor]*NY0*NZ0];
		
		
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank==processor)
			MPI_Isend(send,nx_l*NY0*NZ0,MPI_INT,0,processor,MPI_COMM_WORLD,&request);
			
		
		
		if (rank==0)
			MPI_Irecv(rece,nx_g[processor]*NY0*NZ0,MPI_INT,processor,processor,MPI_COMM_WORLD,&request);
		
		
		if ((rank==0) or (rank==processor))
			MPI_Wait(&request,&status);
			
		
		if (rank==0)
		{
		ofstream out(name.str().c_str(),ios::app);
		
		for(int i=0;i<nx_g[processor];i++)
			if (i+disp[processor]<NX0)
        		for(int j=0; j<NY0; j++)
				for(int k=0;k<NZ0;k++)

				{
				if (rece[i*NY0*NZ0+j*NZ0+k]<=0)
					out<<1<<endl;
				else
					out<<0<<endl;
				}
	
		out.close();
		}
		
		
		if (rank==0)
			delete [] rece;
		
	}

	delete [] send;
	
	
		
}




void output_velocity(int m,double* rho,double** u,int MirX,int MirY,int MirZ,int mir,int*** Solid)	
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	const int mpi_size=MPI :: COMM_WORLD . Get_size ();
	const int root_rank=0;
	
	int nx_g[mpi_size];
	int disp[mpi_size];
	        
	MPI_Status status;
	MPI_Request request;

	double* send;
	double* rece;
	
	
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


	send = new double[nx_l*NY0*NZ0*3];
	for (int i=0;i<nx_l;i++)
		for (int j=0;j<NY0;j++)
			for (int k=0;k<NZ0;k++)
			if (Solid[i][j][k]>0)
			{
			send[i*NY0*NZ0*3+j*NZ0*3+k*3]=u[Solid[i][j][k]][0];
			send[i*NY0*NZ0*3+j*NZ0*3+k*3+1]=u[Solid[i][j][k]][1];
			send[i*NY0*NZ0*3+j*NZ0*3+k*3+2]=u[Solid[i][j][k]][2];
			}
			else
			        {
			        send[i*NY0*NZ0*3+j*NZ0*3+k*3]=0.0;
			        send[i*NY0*NZ0*3+j*NZ0*3+k*3+1]=0.0;
			        send[i*NY0*NZ0*3+j*NZ0*3+k*3+2]=0.0;        
			        }
			        
			        
	
			        
	
	if (rank==0)
	{
	ofstream out(name.str().c_str(),ios::app);
	for(int i=0;i<nx_l;i++)
        	for(int j=0; j<NY0; j++)
			for(int k=0;k<NZ0;k++)
			if (Solid[i][j][k]>0)
				out<<u[Solid[i][j][k]][2]<<" "<<u[Solid[i][j][k]][1]<<" "<<u[Solid[i][j][k]][0]<<endl;
			else
				out<<0.0<<" "<<0.0<<" "<<0.0<<endl;
			
	out.close();
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	for (int processor=1;processor<mpi_size;processor++)
	{	
		
		if (rank==0)
			rece = new double[nx_g[processor]*NY0*NZ0*3];
		
		
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank==processor)
			MPI_Isend(send,nx_l*NY0*NZ0*3,MPI_DOUBLE,0,processor,MPI_COMM_WORLD,&request);
		
		
		
		if (rank==0)
		{
		        
			MPI_Irecv(rece,nx_g[processor]*NY0*NZ0*3,MPI_DOUBLE,processor,processor,MPI_COMM_WORLD,&request);
		}
		
		if ((rank==0) or (rank==processor))
			MPI_Wait(&request,&status);
			
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (rank==0)
		{
		ofstream out(name.str().c_str(),ios::app);
		
		for(int i=0;i<nx_g[processor];i++)
			if (i+disp[processor]<NX0)
        		for(int j=0; j<NY0; j++)
				for(int k=0;k<NZ0;k++)
				out<<rece[i*NY0*NZ0*3+j*NZ0*3+k*3+2]<<" "<<rece[i*NY0*NZ0*3+j*NZ0*3+k*3+1]<<" "<<rece[i*NY0*NZ0*3+j*NZ0*3+k*3]<<endl;	

		
		out.close();
		}
		
		
		if (rank==0)
			delete [] rece;
		
	}

	delete [] send;
	

	
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
	
	double rho_0=1.0;
	
	
	MPI_Status status;
	MPI_Request request;

	double* send;
	double* rece;

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
        
	send = new double[nx_l*NY0*NZ0];
	for (int i=0;i<nx_l;i++)
		for (int j=0;j<NY0;j++)
			for (int k=0;k<NZ0;k++)
			if (Solid[i][j][k]>0)
			send[i*NY0*NZ0+j*NZ0+k]=rho[Solid[i][j][k]];
			else
			send[i*NY0*NZ0+j*NZ0+k]=rho_0;
		
		
		
	if (rank==0)
	{
	ofstream out(name.str().c_str(),ios::app);
	for(int i=0;i<nx_l;i++)
        	for(int j=0; j<NY0; j++)
			for(int k=0;k<NZ0;k++)
			if (Solid[i][j][k]>0)
				out<<rho[Solid[i][j][k]]<<endl;
			else
				out<<rho_0<<endl;
			
	out.close();
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	for (int processor=1;processor<mpi_size;processor++)
	{	
		
		if (rank==0)
			rece = new double[nx_g[processor]*NY0*NZ0];
		
		
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank==processor)
			MPI_Isend(send,nx_l*NY0*NZ0,MPI_DOUBLE,0,processor,MPI_COMM_WORLD,&request);
		
		
		
		if (rank==0)
		{
		        
			MPI_Irecv(rece,nx_g[processor]*NY0*NZ0,MPI_DOUBLE,processor,processor,MPI_COMM_WORLD,&request);
		}
		
		if ((rank==0) or (rank==processor))
			MPI_Wait(&request,&status);
			
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (rank==0)
		{
		ofstream out(name.str().c_str(),ios::app);
		
		for(int i=0;i<nx_g[processor];i++)
			if (i+disp[processor]<NX0)
        		for(int j=0; j<NY0; j++)
				for(int k=0;k<NZ0;k++)
				out<<rece[i*NY0*NZ0+j*NZ0+k]<<endl;	

		
		out.close();
		}
		
		
		if (rank==0)
			delete [] rece;
		
	}

	delete [] send;

}

void Geometry_b(int*** Solid)	
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


void output_velocity_b(int m,double* rho,double** u,int MirX,int MirY,int MirZ,int mir,int*** Solid)	
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	const int mpi_size=MPI :: COMM_WORLD . Get_size ();
	const int root_rank=0;
	
	int* nx_g = new int[mpi_size];
	int* disp = new int[mpi_size];

	
	MPI_Gather(&nx_l,1,MPI_INT,nx_g,1,MPI_INT,root_rank,MPI_COMM_WORLD);
	
	
	if (rank==root_rank)
		{
		disp[0]=0;
		for (int i=0;i<mpi_size;i++)
			nx_g[i]*=(NY+1)*(NZ+1)*3;

		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
		}
	

	double* rbuf_v;
	double* v_storage = new double[nx_l*(NY+1)*(NZ+1)*3];


	for (int i=0;i<nx_l;i++)
		for (int j=0;j<=NY;j++)
			for(int k=0;k<=NZ;k++)
			{
			if (Solid[i][j][k]>0)
				{
				v_storage[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3]=u[Solid[i][j][k]][0];
				v_storage[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+1]=u[Solid[i][j][k]][1];
				v_storage[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+2]=u[Solid[i][j][k]][2];
				}				
			else
				{
				v_storage[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3]=0;
				v_storage[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+1]=0;
				v_storage[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+2]=0;
				}
			}

	if (rank==root_rank)
		rbuf_v= new double[(NX+1)*(NY+1)*(NZ+1)*3];


	//MPI_Barrier(MPI_COMM_WORLD);
	MPI_Gatherv(v_storage,nx_l*(NY+1)*(NZ+1)*3,MPI_DOUBLE,rbuf_v,nx_g,disp,MPI_DOUBLE,root_rank,MPI_COMM_WORLD);


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



	if (rank==root_rank)
	{
	ostringstream name;
	name<<outputfile<<"LBM_velocity_Vector_"<<m<<".vtk";
	ofstream out;
	out.open(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Velocity"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<NX0<<"         "<<NY0<<"         "<<NZ0<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;

	out<<"POINT_DATA     "<<NX0*NY0*NZ0<<endl;
	out<<"VECTORS sample_vectors double"<<endl;
	out<<endl;
	//out<<"LOOKUP_TABLE default"<<endl;
	for(int k=0;k<NZ0;k++)
      		for(int j=0; j<NY0; j++)
			{
			for(int i=0;i<NX0;i++)
        		out<<rbuf_v[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3]<<" "<<rbuf_v[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+1]<<" "<<rbuf_v[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+2]<<" "<<endl;
			//out<<endl;
			}
			
	out.close();
	
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

	if (rank==root_rank)
		{		
		delete [] rbuf_v;
		}
	delete [] nx_g;
	delete [] disp;
	delete [] v_storage;
	

		
}


void output_density_b(int m,double* rho,int MirX,int MirY,int MirZ,int mir,int*** Solid)	
{
        
	int rank = MPI :: COMM_WORLD . Get_rank ();
	const int mpi_size=MPI :: COMM_WORLD . Get_size ();
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
	

	double* rbuf_rho;
	double* rho_storage = new double[nx_l*(NY+1)*(NZ+1)];


	for (int i=0;i<nx_l;i++)
		for (int j=0;j<=NY;j++)
			for(int k=0;k<=NZ;k++)
			{
			if (Solid[i][j][k]>0)
				rho_storage[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=rho[Solid[i][j][k]];
			else
				rho_storage[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=1.0;
			}

	if (rank==root_rank)
		rbuf_rho= new double[(NX+1)*(NY+1)*(NZ+1)];

	
	//MPI_Barrier(MPI_COMM_WORLD);	
	MPI_Gatherv(rho_storage,nx_l*(NY+1)*(NZ+1),MPI_DOUBLE,rbuf_rho,nx_g,disp,MPI_DOUBLE,root_rank,MPI_COMM_WORLD);

	
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

	if (rank==root_rank)
	{
	
	ostringstream name;
	name<<outputfile<<"LBM_Density_"<<m<<".vtk";
	ofstream out;
	out.open(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
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
				out<<"		"<<rbuf_rho[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]<<endl;

	out.close();
				
	}
	
	//int lss=0;
	//if (rank==0)	
	//	for (int i=1;i<=SumCount;i++)
//			lss+=rbuf_rho[i];
//		cout<<lss<<endl;
	//cout<<SumCount<<endl;
	
	if (rank==root_rank)
		{		
		delete [] rbuf_rho;
		}
	delete [] nx_g;
	delete [] disp;
	delete [] rho_storage;

		
}



void output_psi_b(int m,double* psi,int MirX,int MirY,int MirZ,int mir,int*** Solid)	
{
       
	int rank = MPI :: COMM_WORLD . Get_rank ();
	const int mpi_size=MPI :: COMM_WORLD . Get_size ();
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
	

	double* rbuf_psi;
	double* psi_storage = new double[nx_l*(NY+1)*(NZ+1)];


	for (int i=0;i<nx_l;i++)
		for (int j=0;j<=NY;j++)
			for(int k=0;k<=NZ;k++)
			{
			if (Solid[i][j][k]>0)
				psi_storage[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=psi[Solid[i][j][k]];
			else
				psi_storage[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=0.0;
			}

	if (rank==root_rank)
		rbuf_psi= new double[(NX+1)*(NY+1)*(NZ+1)];

	
	//MPI_Barrier(MPI_COMM_WORLD);	
	MPI_Gatherv(psi_storage,nx_l*(NY+1)*(NZ+1),MPI_DOUBLE,rbuf_psi,nx_g,disp,MPI_DOUBLE,root_rank,MPI_COMM_WORLD);

	
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

	if (rank==root_rank)
	{
	
	ostringstream name;
	name<<outputfile<<"LBM_psi_"<<m<<".vtk";
	ofstream out;
	out.open(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
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
				out<<"		"<<rbuf_psi[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]<<endl;

	out.close();
				
	}
	
	//int lss=0;
	//if (rank==0)	
	//	for (int i=1;i<=SumCount;i++)
//			lss+=rbuf_rho[i];
//		cout<<lss<<endl;
	//cout<<SumCount<<endl;
	
	if (rank==root_rank)
		{		
		delete [] rbuf_psi;
		}
	delete [] nx_g;
	delete [] disp;
	delete [] psi_storage;

		
}

void output_psi(int m,double* psi,int MirX,int MirY,int MirZ,int mir,int*** Solid)	
{
       
        
	int rank = MPI :: COMM_WORLD . Get_rank ();
	const int mpi_size=MPI :: COMM_WORLD . Get_size ();
	const int root_rank=0;
	
	double psi_0=0.0;
	
	
	MPI_Status status;
	MPI_Request request;

	double* send;
	double* rece;

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
	name<<outputfile<<"LBM_psi_"<<m<<".vtk";
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
        
	send = new double[nx_l*NY0*NZ0];
	for (int i=0;i<nx_l;i++)
		for (int j=0;j<NY0;j++)
			for (int k=0;k<NZ0;k++)
			if (Solid[i][j][k]>0)
			send[i*NY0*NZ0+j*NZ0+k]=psi[Solid[i][j][k]];
			else
			send[i*NY0*NZ0+j*NZ0+k]=psi_0;
		
		
		
	if (rank==0)
	{
	ofstream out(name.str().c_str(),ios::app);
	for(int i=0;i<nx_l;i++)
        	for(int j=0; j<NY0; j++)
			for(int k=0;k<NZ0;k++)
			if (Solid[i][j][k]>0)
				out<<psi[Solid[i][j][k]]<<endl;
			else
				out<<psi_0<<endl;
			
	out.close();
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	for (int processor=1;processor<mpi_size;processor++)
	{	
		
		if (rank==0)
			rece = new double[nx_g[processor]*NY0*NZ0];
		
		
		MPI_Barrier(MPI_COMM_WORLD);
		if (rank==processor)
			MPI_Isend(send,nx_l*NY0*NZ0,MPI_DOUBLE,0,processor,MPI_COMM_WORLD,&request);
		
		
		
		if (rank==0)
		{
		        
			MPI_Irecv(rece,nx_g[processor]*NY0*NZ0,MPI_DOUBLE,processor,processor,MPI_COMM_WORLD,&request);
		}
		
		if ((rank==0) or (rank==processor))
			MPI_Wait(&request,&status);
			
		MPI_Barrier(MPI_COMM_WORLD);
		
		if (rank==0)
		{
		ofstream out(name.str().c_str(),ios::app);
		
		for(int i=0;i<nx_g[processor];i++)
			if (i+disp[processor]<NX0)
        		for(int j=0; j<NY0; j++)
				for(int k=0;k<NZ0;k++)
				out<<rece[i*NY0*NZ0+j*NZ0+k]<<endl;	

		
		out.close();
		}
		
		
		if (rank==0)
			delete [] rece;
		
	}

	delete [] send;

}




double Comput_Perm(double* psi,double** u,double* Per_l,double* Per_g,int PerDIr)
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	double *rbuf_l, *rbuf_g;
	rbuf_l=new double[mpi_size*3];
	rbuf_g=new double[mpi_size*3];
	
	
	double Perm_l[3];
	double Perm_g[3];
	double error;
	double Q_l[3]={0.0,0.0,0.0};
	double Q_g[3]={0.0,0.0,0.0};

	
	for (int i=1;i<=Count;i++)
	        if (psi[i]>0)
		        {
		                Q_l[0]+=u[i][0];
		                Q_l[1]+=u[i][1];
		                Q_l[2]+=u[i][2];
		        }
		else
		        {
                                Q_g[0]+=u[i][0];
                                Q_g[1]+=u[i][1];
                                Q_g[2]+=u[i][2];
		        
		        
		        }

	MPI_Barrier(MPI_COMM_WORLD);

	

	MPI_Gather(&Q_l,3,MPI_DOUBLE,rbuf_l,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Gather(&Q_g,3,MPI_DOUBLE,rbuf_g,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	
	if (rank==0)
		{
		Q_l[0]=0;Q_l[1]=0;Q_l[2]=0;
		Q_g[0]=0;Q_g[1]=0;Q_g[2]=0;
		for (int i=0;i<mpi_size;i++)
			{
			Q_g[0]+=rbuf_g[i*3+0];
			Q_g[1]+=rbuf_g[i*3+1];
			Q_g[2]+=rbuf_g[i*3+2];
			Q_l[0]+=rbuf_l[i*3+0];
			Q_l[1]+=rbuf_l[i*3+1];
			Q_l[2]+=rbuf_l[i*3+2];
			}

		Perm_l[0]=Q_l[0]/((NX+1)*(NY+1)*(NZ+1))*(niu_l)/gx;
		Perm_l[1]=Q_l[1]/((NX+1)*(NY+1)*(NZ+1))*(niu_l)/gy;
		Perm_l[2]=Q_l[2]/((NX+1)*(NY+1)*(NZ+1))*(niu_l)/gz;

		Perm_g[0]=Q_g[0]/((NX+1)*(NY+1)*(NZ+1))*(niu_g)/gx;
		Perm_g[1]=Q_g[1]/((NX+1)*(NY+1)*(NZ+1))*(niu_g)/gy;
		Perm_g[2]=Q_g[2]/((NX+1)*(NY+1)*(NZ+1))*(niu_g)/gz;
		
		switch(PerDIr)
		{
		case 1:
			error=(Perm_l[0]-Per_l[0])/Per_l[0];break;
		case 2:
			error=(Perm_l[1]-Per_l[1])/Per_l[1];break;
		case 3:
			error=(Perm_l[2]-Per_l[2])/Per_l[2];break;
		default:
			error=(Perm_l[0]-Per_l[0])/Per_l[0];
		}


		Per_l[0]=Perm_l[0];Per_g[0]=Perm_g[0];
		Per_l[1]=Perm_l[1];Per_g[1]=Perm_g[1];
		Per_l[2]=Perm_l[2];Per_g[2]=Perm_g[2];

		}
	
	delete [] rbuf_l;
	delete [] rbuf_g;
	
	
	return (error);
	

}



double Comput_Saturation(double* psi,int*** Solid)
{
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	double S_l,S_g;
	

	
double *rbuf_l,*rbuf_g;

	rbuf_l=new double[mpi_size];
	rbuf_g=new double[mpi_size];

	S_l=0;S_g=0;
	for (int i=1;i<=Count;i++)
			{				
			if (psi[i]>=0) 
			S_l+=1;
			else
			S_g+=1;
			}

	
		
	

		MPI_Gather(&S_l,1,MPI_DOUBLE,rbuf_l,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
		MPI_Gather(&S_g,1,MPI_DOUBLE,rbuf_g,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	if (rank==0)
	{S_l=0;S_g=0;
	for (int i=0;i<mpi_size;i++)
			{
			S_l+=rbuf_l[i];S_g+=rbuf_g[i];
			}
	
	
	
	S_l=S_l/((NX+1)*(NY+1)*(NZ+1)*porosity);
	S_g=S_g/((NX+1)*(NY+1)*(NZ+1)*porosity);
	}

	
	delete [] rbuf_l;
	delete [] rbuf_g;
	
	
	return (S_l);
			

}

