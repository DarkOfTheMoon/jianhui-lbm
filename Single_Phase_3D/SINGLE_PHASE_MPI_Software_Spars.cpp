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


double u_max,u_ave,u_ave2,gx,gy,gz,porosity;

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





double uMax,lat_c,c_s,c_s2,Re,dx,dy,Lx,Ly,dt,rho0,P0,tau_f,niu,error,SFx,SFy,reso;


void tests();

void init_Sparse_read_rock_parallel(int***,int*, int*);

void init(double*, double**, double**, int***);

void periodic_streaming(double** ,double** ,int* ,int***,int*, int*,double*, double**);

void periodic_streaming_MR(double** ,double** ,int* ,int*** ,int* ,int* ,double* ,double** );

void standard_bounceback_boundary(int,double**);

void collision(double*,double** ,double** ,double** , int* ,int***,int*, int*);

void comput_macro_variables( double* ,double**,double** ,double** ,double**  ,int* ,int***);

double Error(double** ,double** ,double*, double*);

void boundary_velocity(int,double,int , double ,int ,double ,int ,double ,int ,double ,int ,double ,double** ,double**,double* ,double** ,int*** );

void boundary_pressure(int ,double ,int , double ,int ,double ,int ,double ,int ,double ,int ,double ,double** ,double**,double** ,double* ,int*** );

void output_velocity(int ,double* ,double** ,int ,int ,int ,int ,int*** );

void output_density(int ,double* ,int ,int ,int ,int ,int*** );	

void Geometry(int*** );	

void output_velocity_b(int ,double* ,double** ,int ,int ,int ,int ,int*** );

void output_density_b(int ,double* ,int ,int ,int ,int ,int*** );	

void Geometry_b(int*** );

void Backup(int ,double* ,double**, double**);

double Comput_Perm(double** u,double*,int,int*);

double S[19];

void Comput_MI(double[19][19], double[19][19]);

int inverse(mat &a);

double feq(int,double, double[3]);

void Suppliment(int*,int***);

void Backup_init(double* rho, double** u, double** f, char[128], char[128],char[128]);

void Parallelize_Geometry();

void Comput_Grop_Perm(double** ,double* ,int ,int* );


const int e[19][3]=
{{0,0,0},{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};

double elat[19][3]=
{{0,0,0},{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};

const double w[19]={1.0/3.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};

const int LR[19]={0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15};
const int FRP[19]={0,0,0,0,0,0,0,1,0,2,0,0,0,0,0,3,0,4,0};
const int FLN[19]={0,0,0,0,0,0,0,0,1,0,2,0,0,0,0,0,3,0,4};
const int RP[5]={1,7,9,15,17};
const int LN[5]={2,8,10,16,18};


int n,nx_l,n_max,in_BC,PerDir,freRe,freDe,freVe,Par_Geo,Par_nx,Par_ny,Par_nz;
int Zoom;


int wr_per,pre_xp,pre_xn,pre_yp,pre_yn,pre_zp,pre_zn,fre_backup,lattice_v;
int vel_xp,vel_xn,vel_yp,vel_yn,vel_zp,vel_zn,Sub_BC,Out_Mode,mode_backup_ini;
double in_vis,p_xp,p_xn,p_yp,p_yn,p_zp,p_zn,dx_input,dt_input;
double inivx,inivy,inivz,v_xp,v_xn,v_yp,v_yn,v_zp,v_zn;
double error_perm;
int par_per_x,par_per_y,par_per_z,per_xp,per_xn,per_yp,per_yn,per_zp,per_zn;

char outputfile[128]="./";
int NCHAR=128;
	char     filename[128], dummy[128+1], backup_rho[128], backup_velocity[128],backup_f[128];
	int      dummyInt;
	
int*** Solid;	
	
int main(int argc , char *argv [])
{	

MPI :: Init (argc , argv );
MPI_Status status ;

double start , finish,remain,elaps;

int rank = MPI :: COMM_WORLD . Get_rank ();
int para_size=MPI :: COMM_WORLD . Get_size ();

int dif,ts,th,tm;
int tse,the,tme;
 

	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

//	int NCHAR=128;
//	char     filename[128], dummy[128+1], backup_rho[128], backup_velocity[128],backup_f[128];
//	int      dummyInt;

	if (rank==0)
	{
	ifstream fin(argv[1]);


							fin.getline(dummy, NCHAR);
	fin >> filename;				fin.getline(dummy, NCHAR);
	fin >> NX >> NY >> NZ;				fin.getline(dummy, NCHAR);
	fin >> n_max;					fin.getline(dummy, NCHAR);
	fin >> reso;					fin.getline(dummy, NCHAR);
	fin >> in_BC;					fin.getline(dummy, NCHAR);
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
	fin >> Out_Mode;				fin.getline(dummy, NCHAR);
	fin >> freVe;					fin.getline(dummy, NCHAR);
	fin >> freDe;					fin.getline(dummy, NCHAR);
							fin.getline(dummy, NCHAR);
	fin >> lattice_v >> dx_input >> dt_input;	fin.getline(dummy, NCHAR);
	fin >> outputfile;				fin.getline(dummy, NCHAR);
	fin >> Sub_BC;					fin.getline(dummy, NCHAR);
	fin >> par_per_x >> par_per_y >>par_per_z;	fin.getline(dummy, NCHAR);
	fin >> per_xp >> per_xn;			fin.getline(dummy, NCHAR);
	fin >> per_yp >> per_yn;			fin.getline(dummy, NCHAR);
	fin >> per_zp >> per_zn;			fin.getline(dummy, NCHAR);
							fin.getline(dummy, NCHAR);
	fin >> fre_backup;                        	fin.getline(dummy, NCHAR);
	fin >>mode_backup_ini;                		fin.getline(dummy, NCHAR);
	fin >> backup_rho;                   		fin.getline(dummy, NCHAR);
	fin >> backup_velocity;                		fin.getline(dummy, NCHAR);
	fin >> backup_f;                        	fin.getline(dummy, NCHAR);
	
	//fin >> EI;					fin.getline(dummy, NCHAR);
	//fin >> q_p;					fin.getline(dummy, NCHAR);
	fin.close();
	
	//cout<<NX<<"    asdfa "<<endl;
	NX=NX-1;NY=NY-1;NZ=NZ-1;
	}

	MPI_Bcast(&filename,128,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&NX,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&NY,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&NZ,1,MPI_INT,0,MPI_COMM_WORLD);
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
	MPI_Bcast(&in_vis,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&outputfile,128,MPI_CHAR,0,MPI_COMM_WORLD);
	//MPI_Bcast(&EI,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&q_p,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&fre_backup,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&mode_backup_ini,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Sub_BC,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&Out_Mode,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&backup_rho,128,MPI_CHAR,0,MPI_COMM_WORLD);MPI_Bcast(&backup_velocity,128,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&backup_f,128,MPI_CHAR,0,MPI_COMM_WORLD);

	MPI_Bcast(&lattice_v,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&dx_input,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&dt_input,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	MPI_Bcast(&par_per_x,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&par_per_y,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&par_per_z,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&per_yp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&per_xp,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&per_yn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&per_zp,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&per_zn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&per_xn,1,MPI_INT,0,MPI_COMM_WORLD);

	
      
int U_max_ref=0;

mirX=0;
mirY=0;
mirZ=0;
mir=1;

Par_Geo=0;

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
Zoom=1;
if (Zoom>1)
	{	
	NX=(NX+1)*Zoom-1;
	NY=(NY+1)*Zoom-1;
	NZ=(NZ+1)*Zoom-1;
	}

//if (Zoom>1)	
//	reso=reso/Zoom;



//	nx_l=(int)((NX+1)/para_size);
//	dif=(NX+1)-nx_l*para_size;
	
//	if (rank>para_size-1-dif)
//		nx_l+=1;
	
	
	

//	if (rank==para_size-1)
//		nx_l+=(NX+1)%para_size;

	double* Permia;
	double* rho;
	double** u;
	double**f;
	double**F;
	double**u0;
	int* SupInv;
	//double* forcex;
	//double* forcey;
	//double* forcez;


	
	int*  Sl;
	int*  Sr;

	
	
//	Solid = new int**[nx_l];
	Sl = new int[(NY+1)*(NZ+1)];
	Sr = new int[(NY+1)*(NZ+1)];

//	Solids = new int**[(NX+1)/Zoom];
//
//	for (int i=0;i<(NX+1)/Zoom;i++)
//		{		
//		Solids[i] = new int*[(NY+1)/Zoom];
//	
//			for (int j=0;j<(NY+1)/Zoom;j++)
//			{
//			Solids[i][j]= new int[(NZ+1)/Zoom];
//			
//
//			}
//		}


//	for (int i=0;i<nx_l;i++)
//		{
//		Solid[i] = new int*[NY+1];
//			for (int j=0;j<=NY;j++)
//			Solid[i][j]= new int[NZ+1];
//		}

	

	
	
	
//	Count = new int[rank];
//	if (!(filename=="NOSOLID"))
//		Read_Rock(Solids,&porosity,filename);
//	else
//		{
//		for(int i=0;i<=NX;i++)	
//			for(int j=0;j<=NY;j++)
//				for(int k=0;k<=NZ;k++)
//				Solids[i][j][k]=0;
//		}

	
        Parallelize_Geometry();
        
	
        
        MPI_Barrier(MPI_COMM_WORLD);
        
	init_Sparse_read_rock_parallel(Solid,Sl,Sr);
	
	

	/*
	for (int i=0;i<(NX+1)/Zoom;i++)
		{
		for (int j=0;j<(NY+1)/Zoom;j++)
			delete [] Solids[i][j];
		delete [] Solids[i];
		}
	delete [] Solids;
*/

	//***************************************************
	//WARRING: SPARSE MATRIX STARTS FROM INDEX 1 NOT 0!!!
	//***************************************************

	Permia = new double[3];
	rho = new double[Count+1];
	//forcex = new double[Count+1];
	//forcey = new double[Count+1];
	//forcez = new double[Count+1];
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
	
	if ((freVe>0) or (freDe>0))
	{
	if (Out_Mode==1)
		Geometry(Solid);
	else
		Geometry_b(Solid);
	}

	
	
	
if (mode_backup_ini==0)
	init(rho,u,f,Solid);
else
        Backup_init(rho,u,f,backup_rho,backup_velocity,backup_f);


if (rank==0)
{
        cout<<endl;
        cout<<"INITIALIZATION COMPLETED"<<endl;
        cout<<endl;
}



//========================================================
char FileName[128];strcpy(FileName,outputfile);
char FileName2[128];strcpy(FileName2,outputfile);
char FileName3[128];strcpy(FileName3,outputfile);
char FileName4[128];strcpy(FileName4,outputfile);


strcat(FileName,"Results.txt");
strcat(FileName2,"Permeability.txt");
strcat(FileName3,"bodyforce.txt");
strcat(FileName4,"Velocity_ave_max.txt");
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

	
	fins.open(FileName4,ios::out);
	fins.close();


//cout<<Sr[1*(NZ+1)+0]<<"  RRRR   "<<rank<<endl;
	
	for(n=0;n<=n_max;n++)
	{
	
	//cout<<"@@@@@@@@@@@   "<<n<<endl;
	collision(rho,u,f,F,SupInv,Solid,Sl,Sr);//cout<<"markcollision"<<endl;

	
	//if (EI==0)
	//	{periodic_streaming(f,F,SupInv,Solid,Sl,Sr,rho,u);}//cout<<"mark"<<endl;}
	//else
	//	periodic_streaming_MR(f,F,SupInv,Solid,Sl,Sr,rho,u);
	
	if ((1-pre_xp)*(1-pre_xn)*(1-pre_yp)*(1-pre_yn)*(1-pre_zp)*(1-pre_zn)==0)
		boundary_pressure(pre_xp,p_xp,pre_xn,p_xn,pre_yp,p_yp,pre_yn,p_yn,pre_zp,p_zp,pre_zn,p_zn,f,F,u,rho,Solid);

	
	if ((1-vel_xp)*(1-vel_xn)*(1-vel_yp)*(1-vel_yn)*(1-vel_zp)*(1-vel_zn)==0)
		boundary_velocity(vel_xp,v_xp,vel_xn,v_xn,vel_yp,v_yp,vel_yn,v_yn,vel_zp,v_zp,vel_zn,v_zn,f,F,rho,u,Solid);


  		comput_macro_variables(rho,u,u0,f,F,SupInv,Solid); 


	
	if(n%freRe==0)
		{       
			
			
			if (rank==0)
			{
			        
			 ofstream fin(FileName,ios::out);       
			 fin<<"The"<<n-freRe<<"th computation result:"<<endl;
			fin<<"The permiability is: "<<Permia[0]*reso*reso*1000<<", "<<Permia[1]*reso*reso*1000<<", "<<Permia[2]*reso*reso*1000<<endl;
			fin<<"The relative error of permiability computing is: "<<error_perm<<endl;	
			fin<<"The Maximum velocity is: "<<setprecision(6)<<u_max<<"   Re="<<Re<<"     Courant Number="<<u_max*dt/dx<<endl;
			fin<<"The max relative error of velocity is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			fin<<"Elapsed time is "<< the<<"h"<<tme<<"m"<<tse<<"s"<<endl;
			fin<<"The expected completion time is "<<th<<"h"<<tm<<"m"<<ts<<"s"<<endl;
			fin<<endl; 
			 fin.close();     
			}

			 error=Error(u,u0,&u_max,&u_ave);if (u_max>=10.0)	U_max_ref+=1;
			error_perm=Comput_Perm(u,Permia,PerDir,SupInv); 
			
			
			//Comput_Grop_Perm(u,Permia,PerDir,SupInv);
			
			
			 if (rank==0)
			{ 
			    
			ofstream fin(FileName,ios::app);          
			finish = MPI_Wtime();
			
			remain=(n_max-n)*((finish-start)/n);
			
			th=int(remain/3600);
			tm=int((remain-th*3600)/60);
			ts=int(remain-(th*3600+tm*60));

			elaps=finish-start;
			the=int(elaps/3600);
			tme=int((elaps-the*3600)/60);
			tse=int(elaps-(the*3600+tme*60));


			fin<<"The"<<n<<"th computation result:"<<endl;
		//=============================================================================================
			fin<<"The permiability is: "<<Permia[0]*reso*reso*1000<<", "<<Permia[1]*reso*reso*1000<<", "<<Permia[2]*reso*reso*1000<<endl;
			fin<<"The relative error of permiability computing is: "<<error_perm<<endl;
		//==============================================================================================

		//==============================================================================================
			Re=u_ave*(NY+1)/(1.0/3.0*(1/s_v-0.5));
			fin<<"The Maximum velocity is: "<<setprecision(6)<<u_max<<"   Re="<<Re<<"     Courant Number="<<u_max*dt/dx<<endl;
			
		//===============================================================================================
			fin<<"The max relative error of velocity is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			fin<<"Elapsed time is "<< the<<"h"<<tme<<"m"<<tse<<"s"<<endl;
			fin<<"The expected completion time is "<<th<<"h"<<tm<<"m"<<ts<<"s"<<endl;
			fin<<endl;
			fin.close();


		if (wr_per==1)
			{
			ofstream finfs(FileName2,ios::app);
			switch(PerDir)
				{
				case 1:
				finfs<<Permia[0]*reso*reso*1000<<" "<<error_perm<<endl;break;
				case 2:
				finfs<<Permia[1]*reso*reso*1000<<" "<<error_perm<<endl;break;
				case 3:
				finfs<<Permia[2]*reso*reso*1000<<" "<<error_perm<<endl;break;
				default:
				finfs<<Permia[0]*reso*reso*1000<<" "<<error_perm<<endl;break;
				}
			finfs.close();
			}

		
			ofstream finf3(FileName3,ios::app);
			finf3<<gx<<endl;
			finf3.close();
			ofstream finf4(FileName4,ios::app);
			finf4<<u_ave<<"  "<<u_max<<" "<<error<<endl;
			finf4.close();
		

			cout<<"The"<<n<<"th computation result:"<<endl;
			cout<<"The Density of point(NX/2,NY/2,NZ/2) is: "<<setprecision(6)
				<<rho[(int)Count/2]<<endl;
		//=============================================================================================
			cout<<"The permiability is: "<<Permia[0]*reso*reso*1000<<", "<<Permia[1]*reso*reso*1000<<", "<<Permia[2]*reso*reso*1000<<endl;
			cout<<"The relative error of permiability computing is: "<<error_perm<<endl;
		//==============================================================================================

		//==============================================================================================
			
			Re=u_ave*(NY+1)/(1.0/3.0*(1/s_v-0.5));
			cout<<"The Maximum velocity is: "<<setprecision(6)<<u_max<<"   Re="<<Re<<"     Courant Number="<<u_max*dt/dx<<endl;
			
		//===============================================================================================
			cout<<"The max relative error of uv is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			cout<<"Elapsed time is "<<the<<"h"<<tme<<"m"<<tse<<"s"<<endl;
			cout<<"The expected completion time is "<<th<<"h"<<tm<<"m"<<ts<<"s"<<endl;
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
				
			if ((fre_backup>=0) and (n%fre_backup==0) and (n>0))
			        Backup(n,rho,u,f);
			        
			if(error!=error) {cout<<"PROGRAM STOP"<<endl;break;};
			if(U_max_ref>=5) {cout<<"PROGRAM STOP DUE TO HIGH VELOCITY"<<endl;break;}
		}	
	}

	if (fre_backup>=0)
			        Backup(n_max,rho,u,f);
	

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
	//delete [] forcex;
	//delete [] forcey;
	//delete [] forcez;
	delete [] SupInv;

	delete [] Sl;
	delete [] Sr;

//	delete [] Permia;
//	delete [] Count;
	finish = MPI_Wtime();
	
	

	MPI_Barrier(MPI_COMM_WORLD);
	if(rank==0)
			{
    			cout<<"Elapsed time is "<< the<<"h"<<tme<<"m"<<tse<<"s"<<endl;
    			cout<<"Accuracy: "<<MPI_Wtick()<<" Second"<<endl;
			ofstream fin(FileName,ios::app);
			fin<<"Elapsed time is "<< the<<"h"<<tme<<"m"<<tse<<"s"<<endl;
			fin<<"Accuracy: "<<MPI_Wtick()<<" Second"<<endl;
			}


	//==========FOR===SCRIPT==RUN=====================
	char FileNameG[128];strcpy(FileNameG,outputfile);
	strcat(FileNameG,"General_Results.txt");
	if (rank==0)
	{
	ofstream fin(FileNameG,ios::app);
	fin<<filename<<"  "<<Permia[0]*reso*reso*1000<<" "<<Permia[1]*reso*reso*1000<<" "<<Permia[2]*reso*reso*1000<<endl;
	fin<<endl;
	fin.close();
	}
	//==============================================

	delete [] Permia;
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

void Parallelize_Geometry()
{
        int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();


int nx_g[mpi_size];
int disp[mpi_size];
int* Solid_rank0;
int pore;
int loc_por[NX+1];
int sum=0;
int sum2=0;
double ave_nx;
int nx_pre,nx_aft,n_i,sum_nx;
int* recv_solid;

int bufsize[mpi_size];
int bufloc[mpi_size];

	for (int i=0;i<=NX;i++)
	        loc_por[i]=0;
	
		if (par_per_x==0)
		{per_xn=0;per_xp=NX;}
	if (par_per_y==0)
		{per_yn=0;per_yp=NY;}
	if (par_per_z==0)
		{per_zn=0;per_zp=NZ;}
	
	
if (rank==0)
{
	FILE *ftest;
	ifstream fin;
	
	ftest = fopen(filename, "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file (" << filename <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest);

	Solid_rank0 = new int[(NX+1)*(NY+1)*(NZ+1)];
	
	fin.open(filename);
	for(int k=0 ; k<=NZ ; k++)
	for(int j=0 ; j<=NY ; j++)
	for(int i=0 ; i<=NX ; i++)
	
	{
		while(true)
		{	
			fin >> pore;
			if( pore == 0.0 || pore == 1.0) break;
		}
		
		Solid_rank0[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=pore;
		if (pore==0)
		        {
		                sum+=1;
		                loc_por[i]+=1;
		                if ((i>=per_xn) and (i<=per_xp) and (j>=per_yn) and (j<=per_yp) and (k>=per_zn) and (k<=per_zp))
		                        sum2+=1;
		        }
	}
	fin.close();
}

        MPI_Bcast(loc_por,NX+1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&sum,1,MPI_INT,0,MPI_COMM_WORLD);
        
         
	nx_pre=0;nx_aft=0;sum_nx=0;
	ave_nx=(double)sum/(mpi_size);
	disp[0]=0;bufloc[0]=0;
	porosity=(double)sum2/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1));
	
	for (int i=0;i<mpi_size-1;i++)
	        {
	        nx_pre=0;nx_aft=0;n_i=0;
	        while (nx_aft<ave_nx*(i+1))
	                {
	                nx_pre=nx_aft;
	                nx_aft+=loc_por[n_i];
	                n_i+=1;
	                }
	        if ((double)(nx_aft-ave_nx*(i+1))>(double)(ave_nx*(i+1)-nx_pre))
	                disp[i+1]=n_i-1;
	        else
	                disp[i+1]=n_i;
	        
	        nx_g[i]=disp[i+1]-disp[i];
	        
	          bufsize[i]=nx_g[i]*(NY+1)*(NZ+1);
	        bufloc[i+1]=bufloc[i]+bufsize[i];
	        
	        }
	  
	  nx_g[mpi_size-1]=(NX+1)-disp[mpi_size-1];
	  bufsize[mpi_size-1]=nx_g[mpi_size-1]*(NY+1)*(NZ+1);
	  nx_l=nx_g[rank];


        Solid = new int**[nx_l];
	for (int i=0;i<nx_l;i++)
	{
		Solid[i] = new int*[NY+1];
			for (int j=0;j<=NY;j++)
			Solid[i][j]= new int[NZ+1];
			
		
	}
	        recv_solid = new int[nx_l*(NY+1)*(NZ+1)];
		
	  
MPI_Barrier(MPI_COMM_WORLD);

        
MPI_Scatterv(Solid_rank0,bufsize,bufloc,MPI_INT,recv_solid,nx_l*(NY+1)*(NZ+1),MPI_INT,0,MPI_COMM_WORLD);

cout<<"GEOMETRY INPUT FILE PARTITIONING FOR PARALLEL READING DONE   Processor No."<<rank<<endl;	  
	  
for (int i=0;i<nx_l;i++)
	                for (int j=0;j<=NY;j++)
	                for (int k=0;k<=NZ;k++)
	                Solid[i][j][k]=recv_solid[i*(NY+1)*(NZ+1)+j*(NZ+1)+k];
	               
	           
	 delete [] recv_solid;	  
	if (rank==0)
	       delete [] Solid_rank0;	

}

void init_Sparse_read_rock_parallel(int*** Solid,int* Sl,int* Sr)
{	
	MPI_Status status[4] ;
	MPI_Request request[4];

	
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

bool mark;
int kk,ip,jp,kp,mean_l,mpi_test;



	int* Sl_send;
	int* Sr_send;
	Sl_send = new int[(NY+1)*(NZ+1)];
	Sr_send = new int[(NY+1)*(NZ+1)];

	Count=1;




for(int i=0;i<nx_l;i++)	
	for(int j=0;j<=NY;j++)
		for(int k=0;k<=NZ;k++)
		{
			if (Solid[i][j][k]==0)
				{
			
				Solid[i][j][k]=Count;
				Count++;
				}
			else
				Solid[i][j][k]=0;

		}
		

	
//	cout<<"GEOMETRY DATA READING DONE-----PROCESSOR No."<<rank<<endl;
	
	
	
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




void init(double* rho, double** u, double** f,int*** Solid)
{	
      

//	int rank = MPI :: COMM_WORLD . Get_rank ();
//	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	

	
	rho0=1.0;dt=1.0/Zoom;dx=1.0/Zoom;
	
	if (lattice_v==1)
		{dx=dx_input;dt=dt_input;}

	lat_c=dx/dt;
	c_s=lat_c/sqrt(3);
	c_s2=lat_c*lat_c/3;

 	
	niu=in_vis;
	tau_f=niu/(c_s2*dt)+0.5;
	//tau_f=3.0*niu/dt+0.5;
	
	s_v=1/tau_f;
      
	double s_other=8*(2-s_v)/(8-s_v);
	double u_tmp[3];

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

	if (lattice_v==1)
	for (int i=0;i<19;i++)
		for (int j=0;j<3;j++)
		elat[i][j]=e[i][j]*lat_c;



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

			//forcex[i]=gx;
			//forcey[i]=gy;
			//forcez[i]=gz;
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

	c2=lat_c*lat_c;c4=c2*c2;
	eu=(elat[k][0]*u[0]+elat[k][1]*u[1]+elat[k][2]*u[2]);
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


void collision(double* rho,double** u,double** f,double** F,int* SupInv,int*** Solid, int* Sl, int* Sr)
{

	MPI_Status status[4] ;
	MPI_Request request[4];
	int mpi_test;
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	int* Gcl = new int[mpi_size];
	int* Gcr = new int[mpi_size];

	
	MPI_Gather(&cl,1,MPI_INT,Gcl,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(Gcl,mpi_size,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Gather(&cr,1,MPI_INT,Gcr,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(Gcr,mpi_size,MPI_INT,0,MPI_COMM_WORLD);


double lm0,lm1,sum;
double usqr,vsqr;
double F_hat[19],GuoF[19],f_eq[19],u_tmp[3];
double m_l[19],m_inv_l[19];
int i,j,m,ip,jp,kp;


const double c_l=lat_c;	


	double* sendl;
	double* sendr;

	double* recvl= new double[Gcl[rank]*5];
	double* recvr = new double[Gcr[rank]*5];

	if (rank==0)
		{
		sendl = new double[Gcr[mpi_size-1]*5];
		sendr = new double[Gcl[rank+1]*5];
		
		
		for (int ka=0;ka<Gcr[mpi_size-1]*5;ka++)
		                sendl[ka]=0;
		                
		 for (int ka=0;ka<Gcl[rank+1]*5;ka++)
		                sendr[ka]=0;
		                
		      
		
		
		}	
		else
		if (rank==mpi_size-1)
			{
			
			sendl = new double[Gcr[rank-1]*5];
			sendr = new double[Gcl[0]*5];
			for (int ka=0;ka<Gcr[rank-1]*5;ka++)
		                sendl[ka]=0;
		                
		        for (int ka=0;ka<Gcl[0]*5;ka++)
		                sendr[ka]=0;

			
			}
			else
			{
	
			sendl = new double[Gcr[rank-1]*5];
			sendr = new double[Gcl[rank+1]*5];
			
			for (int ka=0;ka<Gcr[rank-1]*5;ka++)
		                sendl[ka]=0;
	
		        for (int ka=0;ka<Gcl[rank+1]*5;ka++)
		                sendr[ka]=0;
		

			}
			
	

	for(int ci=1;ci<=Count;ci++)	
	

		{	
				
			i=(int)(SupInv[ci]/((NY+1)*(NZ+1)));
			j=(int)((SupInv[ci]%((NY+1)*(NZ+1)))/(NZ+1));
			m=(int)(SupInv[ci]%(NZ+1));  
			

			//=================FORCE TERM_GUO=========================================
/*
			for (int k=0;k<19;k++)
			{	
			lm0=((elat[k][0]-u[ci][0])*gx+(elat[k][1]-u[ci][1])*gy+(elat[k][2]-u[ci][2])*gz)/c_s2;
			lm1=(elat[k][0]*u[ci][0]+elat[k][1]*u[ci][1]+elat[k][2]*u[ci][2])*(elat[k][0]*gx+elat[k][1]*gy+elat[k][2]*gz)/(c_s2*c_s2);
			GuoF[k]=w[k]*(lm0+lm1);
			//GuoF[k]=0.0;
			}

*/



//==========================
lm0=((+0.000*c_l-u[ci][0])*gx+(+0.000*c_l-u[ci][1])*gy+(+0.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(+0.000*c_l*gx+(+0.000*c_l)*gy+(+0.000*c_l)*gz)/(c_s2*c_s2);
GuoF[0]=w[0]*(lm0+lm1);

lm0=((+1.000*c_l-u[ci][0])*gx+(+0.000*c_l-u[ci][1])*gy+(+0.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(+1.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(+1.000*c_l*gx+(+0.000*c_l)*gy+(+0.000*c_l)*gz)/(c_s2*c_s2);
GuoF[1]=w[1]*(lm0+lm1);

lm0=((-1.000*c_l-u[ci][0])*gx+(+0.000*c_l-u[ci][1])*gy+(+0.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(-1.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(-1.000*c_l*gx+(+0.000*c_l)*gy+(+0.000*c_l)*gz)/(c_s2*c_s2);
GuoF[2]=w[2]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gx+(+1.000*c_l-u[ci][1])*gy+(+0.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(+1.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(+0.000*c_l*gx+(+1.000*c_l)*gy+(+0.000*c_l)*gz)/(c_s2*c_s2);
GuoF[3]=w[3]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gx+(-1.000*c_l-u[ci][1])*gy+(+0.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(-1.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(+0.000*c_l*gx+(-1.000*c_l)*gy+(+0.000*c_l)*gz)/(c_s2*c_s2);
GuoF[4]=w[4]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gx+(+0.000*c_l-u[ci][1])*gy+(+1.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(+1.000*c_l)*u[ci][2])*(+0.000*c_l*gx+(+0.000*c_l)*gy+(+1.000*c_l)*gz)/(c_s2*c_s2);
GuoF[5]=w[5]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gx+(+0.000*c_l-u[ci][1])*gy+(-1.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(-1.000*c_l)*u[ci][2])*(+0.000*c_l*gx+(+0.000*c_l)*gy+(-1.000*c_l)*gz)/(c_s2*c_s2);
GuoF[6]=w[6]*(lm0+lm1);

lm0=((+1.000*c_l-u[ci][0])*gx+(+1.000*c_l-u[ci][1])*gy+(+0.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(+1.000*c_l*u[ci][0]+(+1.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(+1.000*c_l*gx+(+1.000*c_l)*gy+(+0.000*c_l)*gz)/(c_s2*c_s2);
GuoF[7]=w[7]*(lm0+lm1);

lm0=((-1.000*c_l-u[ci][0])*gx+(+1.000*c_l-u[ci][1])*gy+(+0.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(-1.000*c_l*u[ci][0]+(+1.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(-1.000*c_l*gx+(+1.000*c_l)*gy+(+0.000*c_l)*gz)/(c_s2*c_s2);
GuoF[8]=w[8]*(lm0+lm1);

lm0=((+1.000*c_l-u[ci][0])*gx+(-1.000*c_l-u[ci][1])*gy+(+0.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(+1.000*c_l*u[ci][0]+(-1.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(+1.000*c_l*gx+(-1.000*c_l)*gy+(+0.000*c_l)*gz)/(c_s2*c_s2);
GuoF[9]=w[9]*(lm0+lm1);

lm0=((-1.000*c_l-u[ci][0])*gx+(-1.000*c_l-u[ci][1])*gy+(+0.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(-1.000*c_l*u[ci][0]+(-1.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(-1.000*c_l*gx+(-1.000*c_l)*gy+(+0.000*c_l)*gz)/(c_s2*c_s2);
GuoF[10]=w[10]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gx+(+1.000*c_l-u[ci][1])*gy+(+1.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(+1.000*c_l)*u[ci][1]+(+1.000*c_l)*u[ci][2])*(+0.000*c_l*gx+(+1.000*c_l)*gy+(+1.000*c_l)*gz)/(c_s2*c_s2);
GuoF[11]=w[11]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gx+(-1.000*c_l-u[ci][1])*gy+(+1.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(-1.000*c_l)*u[ci][1]+(+1.000*c_l)*u[ci][2])*(+0.000*c_l*gx+(-1.000*c_l)*gy+(+1.000*c_l)*gz)/(c_s2*c_s2);
GuoF[12]=w[12]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gx+(+1.000*c_l-u[ci][1])*gy+(-1.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(+1.000*c_l)*u[ci][1]+(-1.000*c_l)*u[ci][2])*(+0.000*c_l*gx+(+1.000*c_l)*gy+(-1.000*c_l)*gz)/(c_s2*c_s2);
GuoF[13]=w[13]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gx+(-1.000*c_l-u[ci][1])*gy+(-1.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(-1.000*c_l)*u[ci][1]+(-1.000*c_l)*u[ci][2])*(+0.000*c_l*gx+(-1.000*c_l)*gy+(-1.000*c_l)*gz)/(c_s2*c_s2);
GuoF[14]=w[14]*(lm0+lm1);

lm0=((+1.000*c_l-u[ci][0])*gx+(+0.000*c_l-u[ci][1])*gy+(+1.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(+1.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(+1.000*c_l)*u[ci][2])*(+1.000*c_l*gx+(+0.000*c_l)*gy+(+1.000*c_l)*gz)/(c_s2*c_s2);
GuoF[15]=w[15]*(lm0+lm1);

lm0=((-1.000*c_l-u[ci][0])*gx+(+0.000*c_l-u[ci][1])*gy+(+1.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(-1.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(+1.000*c_l)*u[ci][2])*(-1.000*c_l*gx+(+0.000*c_l)*gy+(+1.000*c_l)*gz)/(c_s2*c_s2);
GuoF[16]=w[16]*(lm0+lm1);

lm0=((+1.000*c_l-u[ci][0])*gx+(+0.000*c_l-u[ci][1])*gy+(-1.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(+1.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(-1.000*c_l)*u[ci][2])*(+1.000*c_l*gx+(+0.000*c_l)*gy+(-1.000*c_l)*gz)/(c_s2*c_s2);
GuoF[17]=w[17]*(lm0+lm1);

lm0=((-1.000*c_l-u[ci][0])*gx+(+0.000*c_l-u[ci][1])*gy+(-1.000*c_l-u[ci][2])*gz)/c_s2;
lm1=(-1.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(-1.000*c_l)*u[ci][2])*(-1.000*c_l*gx+(+0.000*c_l)*gy+(-1.000*c_l)*gz)/(c_s2*c_s2);
GuoF[18]=w[18]*(lm0+lm1);

//====================




			//============================================================================


			//=====================equilibrium of moment=================================
			u_tmp[0]=u[ci][0];
			u_tmp[1]=u[ci][1];
			u_tmp[2]=u[ci][2];

			for(int k=0;k<19;k++)
				{
				f_eq[k]=feq(k,rho[ci],u_tmp);
				}
			
			//for (int l=0;l<19;l++)
			//	{
			//	meq[l]=0;
			//	for(int lm=0;lm<19;lm++)
			//	meq[l]+=M[l][lm]*f_eq[lm];				
			//	}

			//============================================================================

			/*
			// ==================   m=Mf matrix calculation  =============================
			// ==================   F_hat=(I-.5*S)MGuoF =====================================
				for (int mi=0; mi<19; mi++)
					{
					m_l[mi]=0;F_hat[mi]=0;meq[mi]=0;
					for (int mj=0; mj<19; mj++)
						{
						m_l[mi]+=M[mi][mj]*f[ci][mj];
						F_hat[mi]+=M[mi][mj]*GuoF[mj];
						meq[mi]+=M[mi][mj]*f_eq[mj];
						}
					F_hat[mi]*=(1-0.5*S[mi]);
					m_l[mi]=m_l[mi]-S[mi]*(m_l[mi]-meq[mi])+dt*F_hat[mi];
					}
			//============================================================================
			*/
			

//==========================
m_l[0]=+1.000*f[ci][0]+1.000*f[ci][1]+1.000*f[ci][2]+1.000*f[ci][3]+1.000*f[ci][4]+1.000*f[ci][5]+1.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]+1.000*f[ci][9]+1.000*f[ci][10]+1.000*f[ci][11]+1.000*f[ci][12]+1.000*f[ci][13]+1.000*f[ci][14]+1.000*f[ci][15]+1.000*f[ci][16]+1.000*f[ci][17]+1.000*f[ci][18];

F_hat[0]=+1.000*GuoF[0]+1.000*GuoF[1]+1.000*GuoF[2]+1.000*GuoF[3]+1.000*GuoF[4]+1.000*GuoF[5]+1.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]+1.000*GuoF[9]+1.000*GuoF[10]+1.000*GuoF[11]+1.000*GuoF[12]+1.000*GuoF[13]+1.000*GuoF[14]+1.000*GuoF[15]+1.000*GuoF[16]+1.000*GuoF[17]+1.000*GuoF[18];

meq[0]=+1.000*f_eq[0]+1.000*f_eq[1]+1.000*f_eq[2]+1.000*f_eq[3]+1.000*f_eq[4]+1.000*f_eq[5]+1.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]+1.000*f_eq[9]+1.000*f_eq[10]+1.000*f_eq[11]+1.000*f_eq[12]+1.000*f_eq[13]+1.000*f_eq[14]+1.000*f_eq[15]+1.000*f_eq[16]+1.000*f_eq[17]+1.000*f_eq[18];

F_hat[0]*=(1-0.5*S[0]);
m_l[0]=m_l[0]-S[0]*(m_l[0]-meq[0])+dt*F_hat[0];
//=======================================

m_l[1]=-30.000*f[ci][0]-11.000*f[ci][1]-11.000*f[ci][2]-11.000*f[ci][3]-11.000*f[ci][4]-11.000*f[ci][5]-11.000*f[ci][6]+8.000*f[ci][7]+8.000*f[ci][8]+8.000*f[ci][9]+8.000*f[ci][10]+8.000*f[ci][11]+8.000*f[ci][12]+8.000*f[ci][13]+8.000*f[ci][14]+8.000*f[ci][15]+8.000*f[ci][16]+8.000*f[ci][17]+8.000*f[ci][18];

F_hat[1]=-30.000*GuoF[0]-11.000*GuoF[1]-11.000*GuoF[2]-11.000*GuoF[3]-11.000*GuoF[4]-11.000*GuoF[5]-11.000*GuoF[6]+8.000*GuoF[7]+8.000*GuoF[8]+8.000*GuoF[9]+8.000*GuoF[10]+8.000*GuoF[11]+8.000*GuoF[12]+8.000*GuoF[13]+8.000*GuoF[14]+8.000*GuoF[15]+8.000*GuoF[16]+8.000*GuoF[17]+8.000*GuoF[18];

meq[1]=-30.000*f_eq[0]-11.000*f_eq[1]-11.000*f_eq[2]-11.000*f_eq[3]-11.000*f_eq[4]-11.000*f_eq[5]-11.000*f_eq[6]+8.000*f_eq[7]+8.000*f_eq[8]+8.000*f_eq[9]+8.000*f_eq[10]+8.000*f_eq[11]+8.000*f_eq[12]+8.000*f_eq[13]+8.000*f_eq[14]+8.000*f_eq[15]+8.000*f_eq[16]+8.000*f_eq[17]+8.000*f_eq[18];

F_hat[1]*=(1-0.5*S[1]);
m_l[1]=m_l[1]-S[1]*(m_l[1]-meq[1])+dt*F_hat[1];
//=======================================

m_l[2]=+12.000*f[ci][0]-4.000*f[ci][1]-4.000*f[ci][2]-4.000*f[ci][3]-4.000*f[ci][4]-4.000*f[ci][5]-4.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]+1.000*f[ci][9]+1.000*f[ci][10]+1.000*f[ci][11]+1.000*f[ci][12]+1.000*f[ci][13]+1.000*f[ci][14]+1.000*f[ci][15]+1.000*f[ci][16]+1.000*f[ci][17]+1.000*f[ci][18];

F_hat[2]=+12.000*GuoF[0]-4.000*GuoF[1]-4.000*GuoF[2]-4.000*GuoF[3]-4.000*GuoF[4]-4.000*GuoF[5]-4.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]+1.000*GuoF[9]+1.000*GuoF[10]+1.000*GuoF[11]+1.000*GuoF[12]+1.000*GuoF[13]+1.000*GuoF[14]+1.000*GuoF[15]+1.000*GuoF[16]+1.000*GuoF[17]+1.000*GuoF[18];

meq[2]=+12.000*f_eq[0]-4.000*f_eq[1]-4.000*f_eq[2]-4.000*f_eq[3]-4.000*f_eq[4]-4.000*f_eq[5]-4.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]+1.000*f_eq[9]+1.000*f_eq[10]+1.000*f_eq[11]+1.000*f_eq[12]+1.000*f_eq[13]+1.000*f_eq[14]+1.000*f_eq[15]+1.000*f_eq[16]+1.000*f_eq[17]+1.000*f_eq[18];

F_hat[2]*=(1-0.5*S[2]);
m_l[2]=m_l[2]-S[2]*(m_l[2]-meq[2])+dt*F_hat[2];
//=======================================

m_l[3]=+0.000*f[ci][0]+1.000*f[ci][1]-1.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+1.000*f[ci][7]-1.000*f[ci][8]+1.000*f[ci][9]-1.000*f[ci][10]+0.000*f[ci][11]+0.000*f[ci][12]+0.000*f[ci][13]+0.000*f[ci][14]+1.000*f[ci][15]-1.000*f[ci][16]+1.000*f[ci][17]-1.000*f[ci][18];

F_hat[3]=+0.000*GuoF[0]+1.000*GuoF[1]-1.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+1.000*GuoF[7]-1.000*GuoF[8]+1.000*GuoF[9]-1.000*GuoF[10]+0.000*GuoF[11]+0.000*GuoF[12]+0.000*GuoF[13]+0.000*GuoF[14]+1.000*GuoF[15]-1.000*GuoF[16]+1.000*GuoF[17]-1.000*GuoF[18];

meq[3]=+0.000*f_eq[0]+1.000*f_eq[1]-1.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+1.000*f_eq[7]-1.000*f_eq[8]+1.000*f_eq[9]-1.000*f_eq[10]+0.000*f_eq[11]+0.000*f_eq[12]+0.000*f_eq[13]+0.000*f_eq[14]+1.000*f_eq[15]-1.000*f_eq[16]+1.000*f_eq[17]-1.000*f_eq[18];

F_hat[3]*=(1-0.5*S[3]);
m_l[3]=m_l[3]-S[3]*(m_l[3]-meq[3])+dt*F_hat[3];
//=======================================

m_l[4]=+0.000*f[ci][0]-4.000*f[ci][1]+4.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+1.000*f[ci][7]-1.000*f[ci][8]+1.000*f[ci][9]-1.000*f[ci][10]+0.000*f[ci][11]+0.000*f[ci][12]+0.000*f[ci][13]+0.000*f[ci][14]+1.000*f[ci][15]-1.000*f[ci][16]+1.000*f[ci][17]-1.000*f[ci][18];

F_hat[4]=+0.000*GuoF[0]-4.000*GuoF[1]+4.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+1.000*GuoF[7]-1.000*GuoF[8]+1.000*GuoF[9]-1.000*GuoF[10]+0.000*GuoF[11]+0.000*GuoF[12]+0.000*GuoF[13]+0.000*GuoF[14]+1.000*GuoF[15]-1.000*GuoF[16]+1.000*GuoF[17]-1.000*GuoF[18];

meq[4]=+0.000*f_eq[0]-4.000*f_eq[1]+4.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+1.000*f_eq[7]-1.000*f_eq[8]+1.000*f_eq[9]-1.000*f_eq[10]+0.000*f_eq[11]+0.000*f_eq[12]+0.000*f_eq[13]+0.000*f_eq[14]+1.000*f_eq[15]-1.000*f_eq[16]+1.000*f_eq[17]-1.000*f_eq[18];

F_hat[4]*=(1-0.5*S[4]);
m_l[4]=m_l[4]-S[4]*(m_l[4]-meq[4])+dt*F_hat[4];
//=======================================

m_l[5]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+1.000*f[ci][3]-1.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]-1.000*f[ci][9]-1.000*f[ci][10]+1.000*f[ci][11]-1.000*f[ci][12]+1.000*f[ci][13]-1.000*f[ci][14]+0.000*f[ci][15]+0.000*f[ci][16]+0.000*f[ci][17]+0.000*f[ci][18];

F_hat[5]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+1.000*GuoF[3]-1.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]-1.000*GuoF[9]-1.000*GuoF[10]+1.000*GuoF[11]-1.000*GuoF[12]+1.000*GuoF[13]-1.000*GuoF[14]+0.000*GuoF[15]+0.000*GuoF[16]+0.000*GuoF[17]+0.000*GuoF[18];

meq[5]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+1.000*f_eq[3]-1.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]-1.000*f_eq[9]-1.000*f_eq[10]+1.000*f_eq[11]-1.000*f_eq[12]+1.000*f_eq[13]-1.000*f_eq[14]+0.000*f_eq[15]+0.000*f_eq[16]+0.000*f_eq[17]+0.000*f_eq[18];

F_hat[5]*=(1-0.5*S[5]);
m_l[5]=m_l[5]-S[5]*(m_l[5]-meq[5])+dt*F_hat[5];
//=======================================

m_l[6]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]-4.000*f[ci][3]+4.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]-1.000*f[ci][9]-1.000*f[ci][10]+1.000*f[ci][11]-1.000*f[ci][12]+1.000*f[ci][13]-1.000*f[ci][14]+0.000*f[ci][15]+0.000*f[ci][16]+0.000*f[ci][17]+0.000*f[ci][18];

F_hat[6]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]-4.000*GuoF[3]+4.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]-1.000*GuoF[9]-1.000*GuoF[10]+1.000*GuoF[11]-1.000*GuoF[12]+1.000*GuoF[13]-1.000*GuoF[14]+0.000*GuoF[15]+0.000*GuoF[16]+0.000*GuoF[17]+0.000*GuoF[18];

meq[6]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]-4.000*f_eq[3]+4.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]-1.000*f_eq[9]-1.000*f_eq[10]+1.000*f_eq[11]-1.000*f_eq[12]+1.000*f_eq[13]-1.000*f_eq[14]+0.000*f_eq[15]+0.000*f_eq[16]+0.000*f_eq[17]+0.000*f_eq[18];

F_hat[6]*=(1-0.5*S[6]);
m_l[6]=m_l[6]-S[6]*(m_l[6]-meq[6])+dt*F_hat[6];
//=======================================

m_l[7]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+1.000*f[ci][5]-1.000*f[ci][6]+0.000*f[ci][7]+0.000*f[ci][8]+0.000*f[ci][9]+0.000*f[ci][10]+1.000*f[ci][11]+1.000*f[ci][12]-1.000*f[ci][13]-1.000*f[ci][14]+1.000*f[ci][15]+1.000*f[ci][16]-1.000*f[ci][17]-1.000*f[ci][18];

F_hat[7]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+1.000*GuoF[5]-1.000*GuoF[6]+0.000*GuoF[7]+0.000*GuoF[8]+0.000*GuoF[9]+0.000*GuoF[10]+1.000*GuoF[11]+1.000*GuoF[12]-1.000*GuoF[13]-1.000*GuoF[14]+1.000*GuoF[15]+1.000*GuoF[16]-1.000*GuoF[17]-1.000*GuoF[18];

meq[7]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+1.000*f_eq[5]-1.000*f_eq[6]+0.000*f_eq[7]+0.000*f_eq[8]+0.000*f_eq[9]+0.000*f_eq[10]+1.000*f_eq[11]+1.000*f_eq[12]-1.000*f_eq[13]-1.000*f_eq[14]+1.000*f_eq[15]+1.000*f_eq[16]-1.000*f_eq[17]-1.000*f_eq[18];

F_hat[7]*=(1-0.5*S[7]);
m_l[7]=m_l[7]-S[7]*(m_l[7]-meq[7])+dt*F_hat[7];
//=======================================

m_l[8]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]-4.000*f[ci][5]+4.000*f[ci][6]+0.000*f[ci][7]+0.000*f[ci][8]+0.000*f[ci][9]+0.000*f[ci][10]+1.000*f[ci][11]+1.000*f[ci][12]-1.000*f[ci][13]-1.000*f[ci][14]+1.000*f[ci][15]+1.000*f[ci][16]-1.000*f[ci][17]-1.000*f[ci][18];

F_hat[8]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]-4.000*GuoF[5]+4.000*GuoF[6]+0.000*GuoF[7]+0.000*GuoF[8]+0.000*GuoF[9]+0.000*GuoF[10]+1.000*GuoF[11]+1.000*GuoF[12]-1.000*GuoF[13]-1.000*GuoF[14]+1.000*GuoF[15]+1.000*GuoF[16]-1.000*GuoF[17]-1.000*GuoF[18];

meq[8]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]-4.000*f_eq[5]+4.000*f_eq[6]+0.000*f_eq[7]+0.000*f_eq[8]+0.000*f_eq[9]+0.000*f_eq[10]+1.000*f_eq[11]+1.000*f_eq[12]-1.000*f_eq[13]-1.000*f_eq[14]+1.000*f_eq[15]+1.000*f_eq[16]-1.000*f_eq[17]-1.000*f_eq[18];

F_hat[8]*=(1-0.5*S[8]);
m_l[8]=m_l[8]-S[8]*(m_l[8]-meq[8])+dt*F_hat[8];
//=======================================

m_l[9]=+0.000*f[ci][0]+2.000*f[ci][1]+2.000*f[ci][2]-1.000*f[ci][3]-1.000*f[ci][4]-1.000*f[ci][5]-1.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]+1.000*f[ci][9]+1.000*f[ci][10]-2.000*f[ci][11]-2.000*f[ci][12]-2.000*f[ci][13]-2.000*f[ci][14]+1.000*f[ci][15]+1.000*f[ci][16]+1.000*f[ci][17]+1.000*f[ci][18];

F_hat[9]=+0.000*GuoF[0]+2.000*GuoF[1]+2.000*GuoF[2]-1.000*GuoF[3]-1.000*GuoF[4]-1.000*GuoF[5]-1.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]+1.000*GuoF[9]+1.000*GuoF[10]-2.000*GuoF[11]-2.000*GuoF[12]-2.000*GuoF[13]-2.000*GuoF[14]+1.000*GuoF[15]+1.000*GuoF[16]+1.000*GuoF[17]+1.000*GuoF[18];

meq[9]=+0.000*f_eq[0]+2.000*f_eq[1]+2.000*f_eq[2]-1.000*f_eq[3]-1.000*f_eq[4]-1.000*f_eq[5]-1.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]+1.000*f_eq[9]+1.000*f_eq[10]-2.000*f_eq[11]-2.000*f_eq[12]-2.000*f_eq[13]-2.000*f_eq[14]+1.000*f_eq[15]+1.000*f_eq[16]+1.000*f_eq[17]+1.000*f_eq[18];

F_hat[9]*=(1-0.5*S[9]);
m_l[9]=m_l[9]-S[9]*(m_l[9]-meq[9])+dt*F_hat[9];
//=======================================

m_l[10]=+0.000*f[ci][0]-4.000*f[ci][1]-4.000*f[ci][2]+2.000*f[ci][3]+2.000*f[ci][4]+2.000*f[ci][5]+2.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]+1.000*f[ci][9]+1.000*f[ci][10]-2.000*f[ci][11]-2.000*f[ci][12]-2.000*f[ci][13]-2.000*f[ci][14]+1.000*f[ci][15]+1.000*f[ci][16]+1.000*f[ci][17]+1.000*f[ci][18];

F_hat[10]=+0.000*GuoF[0]-4.000*GuoF[1]-4.000*GuoF[2]+2.000*GuoF[3]+2.000*GuoF[4]+2.000*GuoF[5]+2.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]+1.000*GuoF[9]+1.000*GuoF[10]-2.000*GuoF[11]-2.000*GuoF[12]-2.000*GuoF[13]-2.000*GuoF[14]+1.000*GuoF[15]+1.000*GuoF[16]+1.000*GuoF[17]+1.000*GuoF[18];

meq[10]=+0.000*f_eq[0]-4.000*f_eq[1]-4.000*f_eq[2]+2.000*f_eq[3]+2.000*f_eq[4]+2.000*f_eq[5]+2.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]+1.000*f_eq[9]+1.000*f_eq[10]-2.000*f_eq[11]-2.000*f_eq[12]-2.000*f_eq[13]-2.000*f_eq[14]+1.000*f_eq[15]+1.000*f_eq[16]+1.000*f_eq[17]+1.000*f_eq[18];

F_hat[10]*=(1-0.5*S[10]);
m_l[10]=m_l[10]-S[10]*(m_l[10]-meq[10])+dt*F_hat[10];
//=======================================

m_l[11]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+1.000*f[ci][3]+1.000*f[ci][4]-1.000*f[ci][5]-1.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]+1.000*f[ci][9]+1.000*f[ci][10]+0.000*f[ci][11]+0.000*f[ci][12]+0.000*f[ci][13]+0.000*f[ci][14]-1.000*f[ci][15]-1.000*f[ci][16]-1.000*f[ci][17]-1.000*f[ci][18];

F_hat[11]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+1.000*GuoF[3]+1.000*GuoF[4]-1.000*GuoF[5]-1.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]+1.000*GuoF[9]+1.000*GuoF[10]+0.000*GuoF[11]+0.000*GuoF[12]+0.000*GuoF[13]+0.000*GuoF[14]-1.000*GuoF[15]-1.000*GuoF[16]-1.000*GuoF[17]-1.000*GuoF[18];

meq[11]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+1.000*f_eq[3]+1.000*f_eq[4]-1.000*f_eq[5]-1.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]+1.000*f_eq[9]+1.000*f_eq[10]+0.000*f_eq[11]+0.000*f_eq[12]+0.000*f_eq[13]+0.000*f_eq[14]-1.000*f_eq[15]-1.000*f_eq[16]-1.000*f_eq[17]-1.000*f_eq[18];

F_hat[11]*=(1-0.5*S[11]);
m_l[11]=m_l[11]-S[11]*(m_l[11]-meq[11])+dt*F_hat[11];
//=======================================

m_l[12]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]-2.000*f[ci][3]-2.000*f[ci][4]+2.000*f[ci][5]+2.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]+1.000*f[ci][9]+1.000*f[ci][10]+0.000*f[ci][11]+0.000*f[ci][12]+0.000*f[ci][13]+0.000*f[ci][14]-1.000*f[ci][15]-1.000*f[ci][16]-1.000*f[ci][17]-1.000*f[ci][18];

F_hat[12]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]-2.000*GuoF[3]-2.000*GuoF[4]+2.000*GuoF[5]+2.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]+1.000*GuoF[9]+1.000*GuoF[10]+0.000*GuoF[11]+0.000*GuoF[12]+0.000*GuoF[13]+0.000*GuoF[14]-1.000*GuoF[15]-1.000*GuoF[16]-1.000*GuoF[17]-1.000*GuoF[18];

meq[12]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]-2.000*f_eq[3]-2.000*f_eq[4]+2.000*f_eq[5]+2.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]+1.000*f_eq[9]+1.000*f_eq[10]+0.000*f_eq[11]+0.000*f_eq[12]+0.000*f_eq[13]+0.000*f_eq[14]-1.000*f_eq[15]-1.000*f_eq[16]-1.000*f_eq[17]-1.000*f_eq[18];

F_hat[12]*=(1-0.5*S[12]);
m_l[12]=m_l[12]-S[12]*(m_l[12]-meq[12])+dt*F_hat[12];
//=======================================

m_l[13]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+1.000*f[ci][7]-1.000*f[ci][8]-1.000*f[ci][9]+1.000*f[ci][10]+0.000*f[ci][11]+0.000*f[ci][12]+0.000*f[ci][13]+0.000*f[ci][14]+0.000*f[ci][15]+0.000*f[ci][16]+0.000*f[ci][17]+0.000*f[ci][18];

F_hat[13]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+1.000*GuoF[7]-1.000*GuoF[8]-1.000*GuoF[9]+1.000*GuoF[10]+0.000*GuoF[11]+0.000*GuoF[12]+0.000*GuoF[13]+0.000*GuoF[14]+0.000*GuoF[15]+0.000*GuoF[16]+0.000*GuoF[17]+0.000*GuoF[18];

meq[13]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+1.000*f_eq[7]-1.000*f_eq[8]-1.000*f_eq[9]+1.000*f_eq[10]+0.000*f_eq[11]+0.000*f_eq[12]+0.000*f_eq[13]+0.000*f_eq[14]+0.000*f_eq[15]+0.000*f_eq[16]+0.000*f_eq[17]+0.000*f_eq[18];

F_hat[13]*=(1-0.5*S[13]);
m_l[13]=m_l[13]-S[13]*(m_l[13]-meq[13])+dt*F_hat[13];
//=======================================

m_l[14]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+0.000*f[ci][7]+0.000*f[ci][8]+0.000*f[ci][9]+0.000*f[ci][10]+1.000*f[ci][11]-1.000*f[ci][12]-1.000*f[ci][13]+1.000*f[ci][14]+0.000*f[ci][15]+0.000*f[ci][16]+0.000*f[ci][17]+0.000*f[ci][18];

F_hat[14]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+0.000*GuoF[7]+0.000*GuoF[8]+0.000*GuoF[9]+0.000*GuoF[10]+1.000*GuoF[11]-1.000*GuoF[12]-1.000*GuoF[13]+1.000*GuoF[14]+0.000*GuoF[15]+0.000*GuoF[16]+0.000*GuoF[17]+0.000*GuoF[18];

meq[14]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+0.000*f_eq[7]+0.000*f_eq[8]+0.000*f_eq[9]+0.000*f_eq[10]+1.000*f_eq[11]-1.000*f_eq[12]-1.000*f_eq[13]+1.000*f_eq[14]+0.000*f_eq[15]+0.000*f_eq[16]+0.000*f_eq[17]+0.000*f_eq[18];

F_hat[14]*=(1-0.5*S[14]);
m_l[14]=m_l[14]-S[14]*(m_l[14]-meq[14])+dt*F_hat[14];
//=======================================

m_l[15]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+0.000*f[ci][7]+0.000*f[ci][8]+0.000*f[ci][9]+0.000*f[ci][10]+0.000*f[ci][11]+0.000*f[ci][12]+0.000*f[ci][13]+0.000*f[ci][14]+1.000*f[ci][15]-1.000*f[ci][16]-1.000*f[ci][17]+1.000*f[ci][18];

F_hat[15]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+0.000*GuoF[7]+0.000*GuoF[8]+0.000*GuoF[9]+0.000*GuoF[10]+0.000*GuoF[11]+0.000*GuoF[12]+0.000*GuoF[13]+0.000*GuoF[14]+1.000*GuoF[15]-1.000*GuoF[16]-1.000*GuoF[17]+1.000*GuoF[18];

meq[15]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+0.000*f_eq[7]+0.000*f_eq[8]+0.000*f_eq[9]+0.000*f_eq[10]+0.000*f_eq[11]+0.000*f_eq[12]+0.000*f_eq[13]+0.000*f_eq[14]+1.000*f_eq[15]-1.000*f_eq[16]-1.000*f_eq[17]+1.000*f_eq[18];

F_hat[15]*=(1-0.5*S[15]);
m_l[15]=m_l[15]-S[15]*(m_l[15]-meq[15])+dt*F_hat[15];
//=======================================

m_l[16]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+1.000*f[ci][7]-1.000*f[ci][8]+1.000*f[ci][9]-1.000*f[ci][10]+0.000*f[ci][11]+0.000*f[ci][12]+0.000*f[ci][13]+0.000*f[ci][14]-1.000*f[ci][15]+1.000*f[ci][16]-1.000*f[ci][17]+1.000*f[ci][18];

F_hat[16]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+1.000*GuoF[7]-1.000*GuoF[8]+1.000*GuoF[9]-1.000*GuoF[10]+0.000*GuoF[11]+0.000*GuoF[12]+0.000*GuoF[13]+0.000*GuoF[14]-1.000*GuoF[15]+1.000*GuoF[16]-1.000*GuoF[17]+1.000*GuoF[18];

meq[16]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+1.000*f_eq[7]-1.000*f_eq[8]+1.000*f_eq[9]-1.000*f_eq[10]+0.000*f_eq[11]+0.000*f_eq[12]+0.000*f_eq[13]+0.000*f_eq[14]-1.000*f_eq[15]+1.000*f_eq[16]-1.000*f_eq[17]+1.000*f_eq[18];

F_hat[16]*=(1-0.5*S[16]);
m_l[16]=m_l[16]-S[16]*(m_l[16]-meq[16])+dt*F_hat[16];
//=======================================

m_l[17]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]-1.000*f[ci][7]-1.000*f[ci][8]+1.000*f[ci][9]+1.000*f[ci][10]+1.000*f[ci][11]-1.000*f[ci][12]+1.000*f[ci][13]-1.000*f[ci][14]+0.000*f[ci][15]+0.000*f[ci][16]+0.000*f[ci][17]+0.000*f[ci][18];

F_hat[17]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]-1.000*GuoF[7]-1.000*GuoF[8]+1.000*GuoF[9]+1.000*GuoF[10]+1.000*GuoF[11]-1.000*GuoF[12]+1.000*GuoF[13]-1.000*GuoF[14]+0.000*GuoF[15]+0.000*GuoF[16]+0.000*GuoF[17]+0.000*GuoF[18];

meq[17]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]-1.000*f_eq[7]-1.000*f_eq[8]+1.000*f_eq[9]+1.000*f_eq[10]+1.000*f_eq[11]-1.000*f_eq[12]+1.000*f_eq[13]-1.000*f_eq[14]+0.000*f_eq[15]+0.000*f_eq[16]+0.000*f_eq[17]+0.000*f_eq[18];

F_hat[17]*=(1-0.5*S[17]);
m_l[17]=m_l[17]-S[17]*(m_l[17]-meq[17])+dt*F_hat[17];
//=======================================

m_l[18]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+0.000*f[ci][7]+0.000*f[ci][8]+0.000*f[ci][9]+0.000*f[ci][10]-1.000*f[ci][11]-1.000*f[ci][12]+1.000*f[ci][13]+1.000*f[ci][14]+1.000*f[ci][15]+1.000*f[ci][16]-1.000*f[ci][17]-1.000*f[ci][18];

F_hat[18]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+0.000*GuoF[7]+0.000*GuoF[8]+0.000*GuoF[9]+0.000*GuoF[10]-1.000*GuoF[11]-1.000*GuoF[12]+1.000*GuoF[13]+1.000*GuoF[14]+1.000*GuoF[15]+1.000*GuoF[16]-1.000*GuoF[17]-1.000*GuoF[18];

meq[18]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+0.000*f_eq[7]+0.000*f_eq[8]+0.000*f_eq[9]+0.000*f_eq[10]-1.000*f_eq[11]-1.000*f_eq[12]+1.000*f_eq[13]+1.000*f_eq[14]+1.000*f_eq[15]+1.000*f_eq[16]-1.000*f_eq[17]-1.000*f_eq[18];

F_hat[18]*=(1-0.5*S[18]);
m_l[18]=m_l[18]-S[18]*(m_l[18]-meq[18])+dt*F_hat[18];



//==========================
m_inv_l[0]=+((double)0X1.AF286BCA1AF28P-5)*m_l[0]+((double)-0X1.9AA066A819A9EP-7)*m_l[1]+((double)0X1.8618618618619P-5)*m_l[2]+((double)0X0P+0)*m_l[3]+((double)0X0P+0)*m_l[4]+((double)0X0P+0)*m_l[5]+((double)0X0P+0)*m_l[6]+((double)0X0P+0)*m_l[7]+((double)0X0P+0)*m_l[8]+((double)0X1P-58)*m_l[9]+((double)0X1P-59)*m_l[10]+((double)0X0P+0)*m_l[11]+((double)0X0P+0)*m_l[12]+((double)0X0P+0)*m_l[13]+((double)0X0P+0)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)0X0P+0)*m_l[16]+((double)-0X0P+0)*m_l[17]+((double)0X0P+0)*m_l[18];

m_inv_l[1]=+((double)0X1.AF286BCA1AF2EP-5)*m_l[0]+((double)-0X1.2D204B4812D21P-8)*m_l[1]+((double)-0X1.0410410410411P-6)*m_l[2]+((double)0X1.999999999999AP-4)*m_l[3]+((double)-0X1.999999999999AP-4)*m_l[4]+((double)0X0P+0)*m_l[5]+((double)0X0P+0)*m_l[6]+((double)0X0P+0)*m_l[7]+((double)0X0P+0)*m_l[8]+((double)0X1.C71C71C71C71CP-5)*m_l[9]+((double)-0X1.C71C71C71C71BP-5)*m_l[10]+((double)0X0P+0)*m_l[11]+((double)0X0P+0)*m_l[12]+((double)0X0P+0)*m_l[13]+((double)0X0P+0)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)0X0P+0)*m_l[16]+((double)-0X0P+0)*m_l[17]+((double)0X0P+0)*m_l[18];

m_inv_l[2]=+((double)0X1.AF286BCA1AF2EP-5)*m_l[0]+((double)-0X1.2D204B4812D2P-8)*m_l[1]+((double)-0X1.041041041041P-6)*m_l[2]+((double)-0X1.999999999999AP-4)*m_l[3]+((double)0X1.999999999999AP-4)*m_l[4]+((double)0X0P+0)*m_l[5]+((double)0X0P+0)*m_l[6]+((double)0X0P+0)*m_l[7]+((double)0X0P+0)*m_l[8]+((double)0X1.C71C71C71C71BP-5)*m_l[9]+((double)-0X1.C71C71C71C71DP-5)*m_l[10]+((double)0X0P+0)*m_l[11]+((double)0X0P+0)*m_l[12]+((double)0X0P+0)*m_l[13]+((double)0X0P+0)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)0X0P+0)*m_l[16]+((double)0X0P+0)*m_l[17]+((double)0X0P+0)*m_l[18];

m_inv_l[3]=+((double)0X1.AF286BCA1AF2AP-5)*m_l[0]+((double)-0X1.2D204B4812D2P-8)*m_l[1]+((double)-0X1.041041041041P-6)*m_l[2]+((double)0X0P+0)*m_l[3]+((double)0X0P+0)*m_l[4]+((double)0X1.9999999999999P-4)*m_l[5]+((double)-0X1.999999999999AP-4)*m_l[6]+((double)-0X1.999999999999AP-58)*m_l[7]+((double)-0X1.999999999999AP-60)*m_l[8]+((double)-0X1.C71C71C71C71CP-6)*m_l[9]+((double)0X1.C71C71C71C71CP-6)*m_l[10]+((double)0X1.5555555555555P-4)*m_l[11]+((double)-0X1.5555555555556P-4)*m_l[12]+((double)0X0P+0)*m_l[13]+((double)0X1P-57)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)0X0P+0)*m_l[16]+((double)-0X1P-57)*m_l[17]+((double)0X0P+0)*m_l[18];

m_inv_l[4]=+((double)0X1.AF286BCA1AF2AP-5)*m_l[0]+((double)-0X1.2D204B4812D2P-8)*m_l[1]+((double)-0X1.041041041041P-6)*m_l[2]+((double)0X0P+0)*m_l[3]+((double)0X0P+0)*m_l[4]+((double)-0X1.9999999999999P-4)*m_l[5]+((double)0X1.999999999999AP-4)*m_l[6]+((double)-0X1.999999999999AP-58)*m_l[7]+((double)-0X1.999999999999AP-60)*m_l[8]+((double)-0X1.C71C71C71C719P-6)*m_l[9]+((double)0X1.C71C71C71C71EP-6)*m_l[10]+((double)0X1.5555555555555P-4)*m_l[11]+((double)-0X1.5555555555556P-4)*m_l[12]+((double)0X0P+0)*m_l[13]+((double)-0X1P-57)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)0X0P+0)*m_l[16]+((double)-0X0P+0)*m_l[17]+((double)0X0P+0)*m_l[18];

m_inv_l[5]=+((double)0X1.AF286BCA1AF29P-5)*m_l[0]+((double)-0X1.2D204B4812D21P-8)*m_l[1]+((double)-0X1.041041041041P-6)*m_l[2]+((double)0X0P+0)*m_l[3]+((double)0X0P+0)*m_l[4]+((double)0X0P+0)*m_l[5]+((double)0X0P+0)*m_l[6]+((double)0X1.999999999999AP-4)*m_l[7]+((double)-0X1.999999999999AP-4)*m_l[8]+((double)-0X1.C71C71C71C71BP-6)*m_l[9]+((double)0X1.C71C71C71C71DP-6)*m_l[10]+((double)-0X1.5555555555555P-4)*m_l[11]+((double)0X1.5555555555556P-4)*m_l[12]+((double)0X0P+0)*m_l[13]+((double)0X0P+0)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)0X0P+0)*m_l[16]+((double)0X0P+0)*m_l[17]+((double)0X0P+0)*m_l[18];

m_inv_l[6]=+((double)0X1.AF286BCA1AF2AP-5)*m_l[0]+((double)-0X1.2D204B4812D21P-8)*m_l[1]+((double)-0X1.041041041041P-6)*m_l[2]+((double)0X0P+0)*m_l[3]+((double)0X0P+0)*m_l[4]+((double)0X0P+0)*m_l[5]+((double)0X0P+0)*m_l[6]+((double)-0X1.9999999999999P-4)*m_l[7]+((double)0X1.999999999999AP-4)*m_l[8]+((double)-0X1.C71C71C71C717P-6)*m_l[9]+((double)0X1.C71C71C71C71EP-6)*m_l[10]+((double)-0X1.5555555555555P-4)*m_l[11]+((double)0X1.5555555555556P-4)*m_l[12]+((double)0X0P+0)*m_l[13]+((double)0X0P+0)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)0X0P+0)*m_l[16]+((double)0X0P+0)*m_l[17]+((double)0X0P+0)*m_l[18];

m_inv_l[7]=+((double)0X1.AF286BCA1AF2AP-5)*m_l[0]+((double)0X1.B6006D801B603P-9)*m_l[1]+((double)0X1.0410410410412P-8)*m_l[2]+((double)0X1.999999999999AP-4)*m_l[3]+((double)0X1.999999999999AP-6)*m_l[4]+((double)0X1.999999999999AP-4)*m_l[5]+((double)0X1.999999999999AP-6)*m_l[6]+((double)0X0P+0)*m_l[7]+((double)0X0P+0)*m_l[8]+((double)0X1.C71C71C71C71DP-6)*m_l[9]+((double)0X1.C71C71C71C71DP-7)*m_l[10]+((double)0X1.5555555555555P-4)*m_l[11]+((double)0X1.5555555555555P-5)*m_l[12]+((double)0X1P-2)*m_l[13]+((double)0X0P+0)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)0X1P-3)*m_l[16]+((double)-0X1P-3)*m_l[17]+((double)0X0P+0)*m_l[18];

m_inv_l[8]=+((double)0X1.AF286BCA1AF2AP-5)*m_l[0]+((double)0X1.B6006D801B605P-9)*m_l[1]+((double)0X1.0410410410412P-8)*m_l[2]+((double)-0X1.999999999999AP-4)*m_l[3]+((double)-0X1.999999999999AP-6)*m_l[4]+((double)0X1.999999999999AP-4)*m_l[5]+((double)0X1.999999999999AP-6)*m_l[6]+((double)0X0P+0)*m_l[7]+((double)0X0P+0)*m_l[8]+((double)0X1.C71C71C71C71DP-6)*m_l[9]+((double)0X1.C71C71C71C71DP-7)*m_l[10]+((double)0X1.5555555555555P-4)*m_l[11]+((double)0X1.5555555555555P-5)*m_l[12]+((double)-0X1P-2)*m_l[13]+((double)-0X0P+0)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)-0X1P-3)*m_l[16]+((double)-0X1P-3)*m_l[17]+((double)0X0P+0)*m_l[18];

m_inv_l[9]=+((double)0X1.AF286BCA1AF2AP-5)*m_l[0]+((double)0X1.B6006D801B603P-9)*m_l[1]+((double)0X1.0410410410412P-8)*m_l[2]+((double)0X1.999999999999AP-4)*m_l[3]+((double)0X1.999999999999AP-6)*m_l[4]+((double)-0X1.999999999999AP-4)*m_l[5]+((double)-0X1.999999999999AP-6)*m_l[6]+((double)0X0P+0)*m_l[7]+((double)0X0P+0)*m_l[8]+((double)0X1.C71C71C71C71DP-6)*m_l[9]+((double)0X1.C71C71C71C71DP-7)*m_l[10]+((double)0X1.5555555555555P-4)*m_l[11]+((double)0X1.5555555555555P-5)*m_l[12]+((double)-0X1P-2)*m_l[13]+((double)0X0P+0)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)0X1P-3)*m_l[16]+((double)0X1P-3)*m_l[17]+((double)0X0P+0)*m_l[18];

m_inv_l[10]=+((double)0X1.AF286BCA1AF2AP-5)*m_l[0]+((double)0X1.B6006D801B605P-9)*m_l[1]+((double)0X1.0410410410412P-8)*m_l[2]+((double)-0X1.999999999999AP-4)*m_l[3]+((double)-0X1.999999999999AP-6)*m_l[4]+((double)-0X1.999999999999AP-4)*m_l[5]+((double)-0X1.999999999999AP-6)*m_l[6]+((double)0X0P+0)*m_l[7]+((double)0X0P+0)*m_l[8]+((double)0X1.C71C71C71C71DP-6)*m_l[9]+((double)0X1.C71C71C71C71DP-7)*m_l[10]+((double)0X1.5555555555555P-4)*m_l[11]+((double)0X1.5555555555555P-5)*m_l[12]+((double)0X1P-2)*m_l[13]+((double)0X0P+0)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)-0X1P-3)*m_l[16]+((double)0X1P-3)*m_l[17]+((double)0X0P+0)*m_l[18];

m_inv_l[11]=+((double)0X1.AF286BCA1AF27P-5)*m_l[0]+((double)0X1.B6006D801B6P-9)*m_l[1]+((double)0X1.041041041041P-8)*m_l[2]+((double)0X0P+0)*m_l[3]+((double)0X0P+0)*m_l[4]+((double)0X1.9999999999998P-4)*m_l[5]+((double)0X1.9999999999998P-6)*m_l[6]+((double)0X1.9999999999998P-4)*m_l[7]+((double)0X1.9999999999998P-6)*m_l[8]+((double)-0X1.C71C71C71C71CP-5)*m_l[9]+((double)-0X1.C71C71C71C71CP-6)*m_l[10]+((double)0X0P+0)*m_l[11]+((double)0X0P+0)*m_l[12]+((double)0X0P+0)*m_l[13]+((double)0X1P-2)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)0X0P+0)*m_l[16]+((double)0X1.0000000000001P-3)*m_l[17]+((double)-0X1.0000000000001P-3)*m_l[18];

m_inv_l[12]=+((double)0X1.AF286BCA1AF2CP-5)*m_l[0]+((double)0X1.B6006D801B602P-9)*m_l[1]+((double)0X1.0410410410412P-8)*m_l[2]+((double)0X0P+0)*m_l[3]+((double)0X0P+0)*m_l[4]+((double)-0X1.999999999999AP-4)*m_l[5]+((double)-0X1.999999999999AP-6)*m_l[6]+((double)0X1.999999999999AP-4)*m_l[7]+((double)0X1.999999999999AP-6)*m_l[8]+((double)-0X1.C71C71C71C71CP-5)*m_l[9]+((double)-0X1.C71C71C71C71CP-6)*m_l[10]+((double)0X0P+0)*m_l[11]+((double)0X0P+0)*m_l[12]+((double)0X0P+0)*m_l[13]+((double)-0X1P-2)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)0X0P+0)*m_l[16]+((double)-0X1P-3)*m_l[17]+((double)-0X1P-3)*m_l[18];

m_inv_l[13]=+((double)0X1.AF286BCA1AF28P-5)*m_l[0]+((double)0X1.B6006D801B602P-9)*m_l[1]+((double)0X1.0410410410412P-8)*m_l[2]+((double)0X0P+0)*m_l[3]+((double)0X0P+0)*m_l[4]+((double)0X1.999999999999AP-4)*m_l[5]+((double)0X1.999999999999AP-6)*m_l[6]+((double)-0X1.999999999999AP-4)*m_l[7]+((double)-0X1.999999999999AP-6)*m_l[8]+((double)-0X1.C71C71C71C71DP-5)*m_l[9]+((double)-0X1.C71C71C71C71DP-6)*m_l[10]+((double)0X0P+0)*m_l[11]+((double)0X0P+0)*m_l[12]+((double)0X0P+0)*m_l[13]+((double)-0X1P-2)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)0X0P+0)*m_l[16]+((double)0X1P-3)*m_l[17]+((double)0X1P-3)*m_l[18];

m_inv_l[14]=+((double)0X1.AF286BCA1AF2CP-5)*m_l[0]+((double)0X1.B6006D801B602P-9)*m_l[1]+((double)0X1.0410410410412P-8)*m_l[2]+((double)0X0P+0)*m_l[3]+((double)0X0P+0)*m_l[4]+((double)-0X1.999999999999AP-4)*m_l[5]+((double)-0X1.999999999999AP-6)*m_l[6]+((double)-0X1.999999999999AP-4)*m_l[7]+((double)-0X1.999999999999AP-6)*m_l[8]+((double)-0X1.C71C71C71C71BP-5)*m_l[9]+((double)-0X1.C71C71C71C71BP-6)*m_l[10]+((double)0X0P+0)*m_l[11]+((double)0X0P+0)*m_l[12]+((double)0X0P+0)*m_l[13]+((double)0X1P-2)*m_l[14]+((double)-0X0P+0)*m_l[15]+((double)0X0P+0)*m_l[16]+((double)-0X1P-3)*m_l[17]+((double)0X1P-3)*m_l[18];

m_inv_l[15]=+((double)0X1.AF286BCA1AF2AP-5)*m_l[0]+((double)0X1.B6006D801B601P-9)*m_l[1]+((double)0X1.0410410410412P-8)*m_l[2]+((double)0X1.999999999999AP-4)*m_l[3]+((double)0X1.999999999999AP-6)*m_l[4]+((double)0X0P+0)*m_l[5]+((double)0X0P+0)*m_l[6]+((double)0X1.999999999999AP-4)*m_l[7]+((double)0X1.999999999999AP-6)*m_l[8]+((double)0X1.C71C71C71C71DP-6)*m_l[9]+((double)0X1.C71C71C71C71DP-7)*m_l[10]+((double)-0X1.5555555555555P-4)*m_l[11]+((double)-0X1.5555555555555P-5)*m_l[12]+((double)-0X0P+0)*m_l[13]+((double)0X0P+0)*m_l[14]+((double)0X1P-2)*m_l[15]+((double)-0X1P-3)*m_l[16]+((double)0X0P+0)*m_l[17]+((double)0X1P-3)*m_l[18];

m_inv_l[16]=+((double)0X1.AF286BCA1AF2AP-5)*m_l[0]+((double)0X1.B6006D801B603P-9)*m_l[1]+((double)0X1.0410410410412P-8)*m_l[2]+((double)-0X1.999999999999AP-4)*m_l[3]+((double)-0X1.999999999999AP-6)*m_l[4]+((double)0X0P+0)*m_l[5]+((double)0X0P+0)*m_l[6]+((double)0X1.999999999999AP-4)*m_l[7]+((double)0X1.999999999999AP-6)*m_l[8]+((double)0X1.C71C71C71C71DP-6)*m_l[9]+((double)0X1.C71C71C71C71DP-7)*m_l[10]+((double)-0X1.5555555555555P-4)*m_l[11]+((double)-0X1.5555555555555P-5)*m_l[12]+((double)0X0P+0)*m_l[13]+((double)0X0P+0)*m_l[14]+((double)-0X1P-2)*m_l[15]+((double)0X1P-3)*m_l[16]+((double)0X0P+0)*m_l[17]+((double)0X1P-3)*m_l[18];

m_inv_l[17]=+((double)0X1.AF286BCA1AF2AP-5)*m_l[0]+((double)0X1.B6006D801B601P-9)*m_l[1]+((double)0X1.0410410410412P-8)*m_l[2]+((double)0X1.999999999999AP-4)*m_l[3]+((double)0X1.999999999999AP-6)*m_l[4]+((double)0X0P+0)*m_l[5]+((double)0X0P+0)*m_l[6]+((double)-0X1.999999999999AP-4)*m_l[7]+((double)-0X1.999999999999AP-6)*m_l[8]+((double)0X1.C71C71C71C71DP-6)*m_l[9]+((double)0X1.C71C71C71C71DP-7)*m_l[10]+((double)-0X1.5555555555555P-4)*m_l[11]+((double)-0X1.5555555555555P-5)*m_l[12]+((double)0X0P+0)*m_l[13]+((double)-0X0P+0)*m_l[14]+((double)-0X1P-2)*m_l[15]+((double)-0X1P-3)*m_l[16]+((double)0X0P+0)*m_l[17]+((double)-0X1P-3)*m_l[18];

m_inv_l[18]=+((double)0X1.AF286BCA1AF2AP-5)*m_l[0]+((double)0X1.B6006D801B603P-9)*m_l[1]+((double)0X1.0410410410412P-8)*m_l[2]+((double)-0X1.999999999999AP-4)*m_l[3]+((double)-0X1.999999999999AP-6)*m_l[4]+((double)0X0P+0)*m_l[5]+((double)0X0P+0)*m_l[6]+((double)-0X1.999999999999AP-4)*m_l[7]+((double)-0X1.999999999999AP-6)*m_l[8]+((double)0X1.C71C71C71C71DP-6)*m_l[9]+((double)0X1.C71C71C71C71DP-7)*m_l[10]+((double)-0X1.5555555555555P-4)*m_l[11]+((double)-0X1.5555555555555P-5)*m_l[12]+((double)-0X0P+0)*m_l[13]+((double)0X0P+0)*m_l[14]+((double)0X1P-2)*m_l[15]+((double)0X1P-3)*m_l[16]+((double)0X0P+0)*m_l[17]+((double)-0X1P-3)*m_l[18];

//====================

			// ==================   f=M_-1m matrix calculation and streaming =============================
		for (int mi=0; mi<19; mi++)
			{
			//sum=0;	
			//for (int mj=0; mj<19; mj++)
			//	sum+=MI[mi][mj]*m_l[mj];
			sum=m_inv_l[mi];
			
			ip=i+e[mi][0];
			jp=j+e[mi][1];if (jp<0) {jp=NY;}; if (jp>NY) {jp=0;};
			kp=m+e[mi][2];if (kp<0) {kp=NZ;}; if (kp>NZ) {kp=0;};

				

			if (ip<0) 
				if (Sl[jp*(NZ+1)+kp]>0)
				{
				sendl[(Sl[jp*(NZ+1)+kp]-1)*5+FLN[mi]]=sum;
				
				}
				else
				{
				F[ci][LR[mi]]=sum;
				
				
				}
					
			
		
			if (ip>=nx_l)
				if (Sr[jp*(NZ+1)+kp]>0)
				{
				sendr[(Sr[jp*(NZ+1)+kp]-1)*5+FRP[mi]]=sum;
				
				}
				else
				{
				F[ci][LR[mi]]=sum;
				
				}

			if ((ip>=0) and (ip<nx_l)) 
				if (Solid[ip][jp][kp]>0)
				{
				F[Solid[ip][jp][kp]][mi]=sum;
				
				}
				else
				{
				F[ci][LR[mi]]=sum;
				}



			
			}
			//============================================================================
			
			
			}

/*
	
	int dest_l,dest_r;
	if (rank+1>mpi_size-1)
		dest_r=0;
	else
		dest_r=rank+1;

	if (rank-1<0)
		dest_l=mpi_size-1;
	else
		dest_l=rank-1;

MPI_Barrier(MPI_COMM_WORLD);
MPI_Sendrecv(sendr,Gcl[dest_r]*5,MPI_DOUBLE,dest_r,rank*2+1,recvl,Gcl[rank]*5,MPI_DOUBLE,dest_l,(dest_l)*2+1, MPI_COMM_WORLD,&status[0]);
MPI_Sendrecv(sendl,Gcr[dest_l]*5,MPI_DOUBLE,dest_l,rank*2,recvr,Gcr[rank]*5,MPI_DOUBLE,dest_r,(dest_r)*2, MPI_COMM_WORLD,&status[1]);
*/



	if (rank==0)
		{
		
		
		MPI_Isend(sendr, Gcl[1]*5, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(sendl, Gcr[mpi_size-1]*5, MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(recvr, Gcr[0]*5, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(recvl, Gcl[0]*5, MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );
		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(sendl, Gcr[rank-1]*5, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr, Gcl[0]*5,  MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl, Gcl[rank]*5, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(recvr, Gcr[rank]*5, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request[3]);
			
			}
			else
			{

      			MPI_Isend(sendl, Gcr[rank-1]*5, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr, Gcl[rank+1]*5, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl, Gcl[rank]*5, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(recvr, Gcr[rank]*5, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);
			
			
			};

	
	MPI_Waitall(4,request, status);

	MPI_Testall(4,request,&mpi_test,status);	



		
			for(i=1;i<=Gcl[rank];i++)
				for (int lm=0;lm<5;lm++)
				if (recvl[(i-1)*5+lm]>0)
			        F[i][RP[lm]]=recvl[(i-1)*5+lm];
			        
			for(j=Count-Gcr[rank]+1;j<=Count;j++)
				for (int lm=0;lm<5;lm++)
				if (recvr[(j-(Count-Gcr[rank]+1))*5+lm]>0)
			        F[j][LN[lm]]=recvr[(j-(Count-Gcr[rank]+1))*5+lm];
			        
			

                        

	delete [] Gcl;
	delete [] Gcr;
	delete [] sendl;
	delete [] sendr;
	delete [] recvl;
	delete [] recvr;
		

}



void comput_macro_variables( double* rho,double** u,double** u0,double** f,double** F,int* SupInv,int*** Solid)
{
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();






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
					u[i][0]+=elat[k][0]*f[i][k];
					u[i][1]+=elat[k][1]*f[i][k];
					u[i][2]+=elat[k][2]*f[i][k];
					}
				

				u[i][0]=(u[i][0]+dt*gx/2)/rho[i];
				u[i][1]=(u[i][1]+dt*gy/2)/rho[i];
				u[i][2]=(u[i][2]+dt*gz/2)/rho[i];
				
				
		
			}




}



void boundary_velocity(int xp,double v_xp,int xn, double v_xn,int yp,double v_yp,int yn,double v_yn,int zp,double v_zp,int zn,double v_zn,double** f,double** F,double* rho,double** u,int*** Solid)

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


//Equilibrium boundary condition (Use equilibrium distribution to update the distributions of particles on the boundaries)
if (Sub_BC==0)
{
if ((yp-1)*(yn-1)==0)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
		{
		if ((yp==1)  && (Solid[i][NY][k]>0))   
		        F[Solid[i][NY][k]][ks]=feq(ks,1.0,u_yp);
		if ((yn==1) && (Solid[i][0][k]>0))
		       F[Solid[i][0][k]][ks]=feq(ks,1.0,u_yn);      
		}

if ((zp-1)*(zn-1)==0)		
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
		for (int ks=0;ks<Q;ks++)
		{
		if ((zp==1) && (Solid[i][j][NZ]>0))    
		        F[Solid[i][j][NZ]][ks]=feq(ks,1.0,u_zp); 
		if ((zn==1) && (Solid[i][j][0]>0))
		        F[Solid[i][j][0]][ks]=feq(ks,1.0,u_zn);
		}


if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
			F[Solid[nx_l-1][j][k]][ks]=feq(ks,1.0,u_xp);
			



if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
			F[Solid[0][j][k]][ks]=feq(ks,1.0,u_xn);
			
	

}

//Equilibrium boundary the value of the pressure of the boundary particles are set by using the value of neighbouring points
if (Sub_BC==1)
{
if ((yp-1)*(yn-1)==0)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
		{
		if ((yp==1)  && (Solid[i][NY][k]>0))
		        //if (Solid[i][NY-1][k]>0)
		         //        F[Solid[i][NY][k]][ks]=feq(ks,rho[Solid[i][NY-1][k]],u_yp);
		         //else
		                 F[Solid[i][NY][k]][ks]=feq(ks,1.0,u_yp);
		if ((yn==1) && (Solid[i][0][k]>0))
		        //if (Solid[i][1][k]>0)
		        //         F[Solid[i][0][k]][ks]=feq(ks,rho[Solid[i][1][k]],u_yn);
		        // else
		                 F[Solid[i][0][k]][ks]=feq(ks,1.0,u_yn);      
		}

if ((zp-1)*(zn-1)==0)		
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
		for (int ks=0;ks<Q;ks++)
		{
		if ((zp==1) && (Solid[i][j][NZ]>0)) 
		       // if (Solid[i][j][NZ-1]>0)
		        //        F[Solid[i][j][NZ]][ks]=feq(ks,rho[Solid[i][j][NZ-1]],u_zp);
		        //else
		                F[Solid[i][j][NZ]][ks]=feq(ks,1.0,u_zp); 
		if ((zn==1) && (Solid[i][j][0]>0))
		        //if (Solid[i][j][1]>0)
		       //         F[Solid[i][j][0]][ks]=feq(ks,rho[Solid[i][j][1]],u_zn);
		       // else
		                F[Solid[i][j][0]][ks]=feq(ks,1.0,u_zn);
		}


if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
		        if (Solid[nx_l-1][j][k]>0)
			//        F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho[Solid[nx_l-2][j][k]],u_xp);
			//else
			         F[Solid[nx_l-1][j][k]][ks]=feq(ks,1.0,u_xp);
			



if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
		        if (Solid[0][j][k]>0)
		        //        F[Solid[0][j][k]][ks]=feq(ks,rho[Solid[1][j][k]],u_xn);
		        //else
		                F[Solid[0][j][k]][ks]=feq(ks,1.0,u_xn);
			
	

}


//Non-Equilibrium boundary condition the distribution of BC points is updated using static pressure: rho=1.0
if (Sub_BC==2)
{
if (yp==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
		if (Solid[i][NY][k]>0)
		for (int ks=0;ks<Q;ks++)
		{
		if (NY+e[ks][1]<NY)
			//if (Solid[i][NY-1][k]>0) 
			//	F[Solid[i][NY][k]][ks]=feq(LR[ks],rho[Solid[i][NY-1][k]],u_yp)-F[Solid[i][NY][k]][LR[ks]]+feq(ks,rho[Solid[i][NY-1][k]],u_yp);
			//else
				F[Solid[i][NY][k]][ks]=feq(LR[ks],1.0,u_yp)-F[Solid[i][NY][k]][LR[ks]]-feq(ks,1.0,u_yp);

			}


if (yn==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
		if (Solid[i][0][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
			if (e[ks][1]>0)
				//if (Solid[i][1][k]>0) 
				//	F[Solid[i][0][k]][ks]=feq(LR[ks],rho[Solid[i][1][k]],u_yn)-F[Solid[i][0][k]][LR[ks]]+feq(ks,rho[Solid[i][1][k]],u_yn);
				//else
					F[Solid[i][0][k]][ks]=feq(LR[ks],1.0,u_yn)-F[Solid[i][0][k]][LR[ks]]+feq(ks,1.0,u_yn);

			}



if (zp==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][NZ]>0)
		for (int ks=0;ks<Q;ks++)
			{
			if (e[ks][2]<0)
				//if (Solid[i][j][NZ-1]>0) 
				//	F[Solid[i][j][NZ]][ks]=feq(LR[ks],rho[Solid[i][j][NZ-1]],u_zp)-F[Solid[i][j][NZ]][LR[ks]]+feq(ks,rho[Solid[i][j][NZ-1]],u_zp);
				//else
					F[Solid[i][j][NZ]][ks]=feq(LR[ks],1.0,u_zp)-F[Solid[i][j][NZ]][LR[ks]]+feq(ks,1.0,u_zp);
		
			}



if (zn==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][0]>0)
		for (int ks=0;ks<Q;ks++)
			{
			if (e[ks][2]>0)
				//if (Solid[i][j][1]>0) 
				//	F[Solid[i][j][0]][ks]=feq(LR[ks],rho[Solid[i][j][1]],u_zn)-F[Solid[i][j][0]][LR[ks]]+feq(ks,rho[Solid[i][j][1]],u_zn);
				//else
					F[Solid[i][j][0]][ks]=feq(LR[ks],1.0,u_zn)-F[Solid[i][j][0]][LR[ks]]+feq(ks,1.0,u_zn);
	
			}




if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
	if (Solid[nx_l-1][k][j]>0)
		for (int ks=0;ks<Q;ks++)
			{
			if (e[ks][0]<0)
				//if (Solid[nx_l-2][j][k]>0) 
				//	F[Solid[nx_l-1][j][k]][ks]=feq(LR[ks],rho[Solid[nx_l-2][j][k]],u_xp)-F[Solid[nx_l-1][j][k]][LR[ks]]+feq(ks,rho[Solid[nx_l-2][j][k]],u_xp);
				//else
					F[Solid[nx_l-1][j][k]][ks]=feq(LR[ks],1.0,u_xp)-F[Solid[nx_l-1][j][k]][LR[ks]]+feq(ks,1.0,u_xp);
		
			}



if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
	if (Solid[0][j][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
			if (e[ks][0]>0)
				//if (Solid[1][j][k]>0) 
				//	F[Solid[0][j][k]][ks]=feq(LR[ks],rho[Solid[1][j][k]],u_xn)-F[Solid[0][j][k]][LR[ks]]+feq(ks,rho[Solid[0][j][k]],u_xn);
				//else
					F[Solid[0][j][k]][ks]=feq(LR[ks],1.0,u_xn)-F[Solid[0][j][k]][LR[ks]]+feq(ks,1.0,u_xn);
		
			}
	
}

//Non-Equilibrium boundary condition use neighbouring points' value to update the BC points
if (Sub_BC==3)
{
if (yp==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
	if (Solid[i][NY][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
				//if (Solid[i][NY-1][k]>0) 
				//	F[Solid[i][NY][k]][ks]=feq(ks,rho[Solid[i][NY-1][k]],u_yp)+f[Solid[i][NY-1][k]][ks]-feq(ks,rho[Solid[i][NY-1][k]],u[Solid[i][NY-1][k]]);
				//else
					F[Solid[i][NY][k]][ks]=feq(ks,rho[Solid[i][NY-1][k]],u_yp);

			}


if (yn==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
	if (Solid[i][0][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
				//if (Solid[i][1][k]>0) 
				//	F[Solid[i][0][k]][ks]=feq(ks,rho[Solid[i][1][k]],u_yn)+f[Solid[i][1][k]][ks]-feq(ks,rho[Solid[i][1][k]],u[Solid[i][1][k]]);
				//else
					F[Solid[i][0][k]][ks]=feq(ks,rho[Solid[i][1][k]],u_yn);

			}



if (zp==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][NZ]>0)
		for (int ks=0;ks<Q;ks++)
			{
				//if (Solid[i][j][NZ-1]>0) 
				//	F[Solid[i][j][NZ]][ks]=feq(ks,rho[Solid[i][j][NZ-1]],u_zp)+f[Solid[i][j][NZ-1]][ks]-feq(ks,rho[Solid[i][j][NZ-1]],u[Solid[i][j][NZ-1]]);
				//else
					F[Solid[i][j][NZ]][ks]=feq(ks,rho[Solid[i][j][NZ-1]],u_zp);
		
			}



if (zn==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][0]>0)
		for (int ks=0;ks<Q;ks++)
			{
				//if (Solid[i][j][1]>0) 
				//	F[Solid[i][j][0]][ks]=feq(ks,rho[Solid[i][j][1]],u_zn)+f[Solid[i][j][1]][ks]-feq(ks,rho[Solid[i][j][1]],u[Solid[i][j][1]]);
				//else
					F[Solid[i][j][0]][ks]=feq(ks,rho[Solid[i][j][1]],u_zn);
	
			}




if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
	if (Solid[nx_l-1][k][j]>0)
		for (int ks=0;ks<Q;ks++)
			{
				//if (Solid[nx_l-2][j][k]>0) 
				//	F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho[Solid[nx_l-2][j][k]],u_xp)+f[Solid[nx_l-2][j][k]][ks]-feq(ks,rho[Solid[nx_l-2][j][k]],u[Solid[nx_l-2][j][k]]);
				//else
					F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho[Solid[nx_l-2][j][k]],u_xp);
		
			}



if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
	if (Solid[0][j][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
				//if (Solid[1][j][k]>0) 
				//	F[Solid[0][j][k]][ks]=feq(ks,rho[Solid[1][j][k]],u_xn)+f[Solid[1][j][k]][ks]-feq(ks,rho[Solid[1][j][k]],u[Solid[1][j][k]]);
				//else
					F[Solid[0][j][k]][ks]=feq(ks,rho[Solid[1][j][k]],u_xn);
		
			}
	
}


}



void boundary_pressure(int xp,double rho_xp,int xn, double rho_xn,int yp,double rho_yp,int yn,double rho_yn,int zp,double rho_zp,int zn,double rho_zn,double** f,double** F,double** u,double* rho,int*** Solid)


{

int Q=19;
double u_ls[3]={0,0,0};
int rank = MPI :: COMM_WORLD . Get_rank ();
int mpi_size=MPI :: COMM_WORLD . Get_size ();
double m_l[19];

//Equilibriun boundary condition. velocities of boundary particles are set as 0.0
if (Sub_BC==0)
{
if ((yp-1)*(yn-1)==0)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)	
	        for (int ks=0;ks<Q;ks++)
		{
		if ((yp==1) && (Solid[i][NY][k]>0))
		        F[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls);
		if ((yn==1) && (Solid[i][0][k]>0))
		        F[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls);
		
		}
	        
if ((zp-1)*(zn-1)==0)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
		for (int ks=0;ks<Q;ks++)
		{
		if ((zp==1) && (Solid[i][j][NZ]>0))
		        F[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls);
		if ((zn==1) && (Solid[i][j][0]>0))
		        F[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls);
		}

		
if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
        for (int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
		if (Solid[nx_l-1][j][k]>0)
			F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls);
			

if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)		
		if (Solid[0][j][k]>0)
			F[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls);
}
       


//Equilibruim boundary condition: Use the value of neighbouring points values to update the distributions of BC points
 if (Sub_BC==1)
{
if (yp==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
	if (Solid[i][NY][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
				if (Solid[i][NY-1][k]>0)
					{
					u_ls[0]=u[Solid[i][NY-1][k]][0];
					u_ls[1]=u[Solid[i][NY-1][k]][1];
					u_ls[2]=u[Solid[i][NY-1][k]][2];
					F[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls);
					}
			}
	


if (yn==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
	if (Solid[i][0][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
				if (Solid[i][1][k]>0)
					{
					u_ls[0]=u[Solid[i][1][k]][0];
					u_ls[1]=u[Solid[i][1][k]][1];
					u_ls[2]=u[Solid[i][1][k]][2];
					F[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls);
					}
			else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls);
					}
			}



if (zp==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][NZ]>0)
		for (int ks=0;ks<Q;ks++)
			{			
				if (Solid[i][j][NZ-1]>0)
					{
					u_ls[0]=u[Solid[i][j][NZ-1]][0];
					u_ls[1]=u[Solid[i][j][NZ-1]][1];
					u_ls[2]=u[Solid[i][j][NZ-1]][2];
					F[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls);
					}

			}



if (zn==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][0]>0)
		for (int ks=0;ks<Q;ks++)
			{
				if (Solid[i][j][1]>0)
					{
					u_ls[0]=u[Solid[i][j][1]][0];
					u_ls[1]=u[Solid[i][j][1]][1];
					u_ls[2]=u[Solid[i][j][1]][2];
					F[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls);
					}
			
			}
	



if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
	if (Solid[nx_l-1][j][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
				if (Solid[nx_l-2][j][k]>0)
					{
					u_ls[0]=u[Solid[nx_l-2][j][k]][0];
					u_ls[1]=u[Solid[nx_l-2][j][k]][1];
					u_ls[2]=u[Solid[nx_l-2][j][k]][2];
					//f[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls)+f[Solid[nx_l-2][j][k]][ks]-feq(ks,rho[Solid[nx_l-2][j][k]],u[Solid[nx_l-2][j][k]]);
					F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls);
					}
				else	
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;
					u_ls[2]=0.0;
					F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls);
					}

			}
			
			



if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
	if (Solid[0][j][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
				if(Solid[1][j][k]>0)
					{
					u_ls[0]=u[Solid[1][j][k]][0];
					u_ls[1]=u[Solid[1][j][k]][1];
					u_ls[2]=u[Solid[1][j][k]][1];
					//f[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls)+f[Solid[1][j][k]][ks]-feq(ks,rho[Solid[1][j][k]],u[Solid[1][j][k]]);
					F[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;
					u_ls[2]=0.0;
					F[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls);
					}
			}
	

}	  


//Non-equilibrium boundary condition  
if (Sub_BC==2)
{
if (yp==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)	
	        if (Solid[i][NY][k]>0)
			if (Solid[i][NY-1][k]>0)
			{
			for (int mi=0; mi<19; mi++)
					{
					m_l[mi]=0;
					for (int mj=0; mj<19; mj++)
						m_l[mi]+=M[mi][mj]*f[Solid[i][NY-1][k]][mj];
					}

			m_l[0]=rho_yp;
			for (int mi=0; mi<19; mi++)
				{
				F[Solid[i][NY][k]][mi]=0;
				for (int mj=0; mj<19; mj++)
					F[Solid[i][NY][k]][mi]+=MI[mi][mj]*m_l[mj];
				}
			}
			else
			{
			u_ls[0]=0.0;u_ls[1]=0.0;u_ls[2]=0.0;
			for (int ks=0;ks<Q;ks++)
				F[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls);
			
			}

		

if (yn==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
	if (Solid[i][0][k]>0)	
		if (Solid[i][1][k]>0)
			{
			for (int mi=0; mi<19; mi++)
					{
					m_l[mi]=0;
					for (int mj=0; mj<19; mj++)
						m_l[mi]+=M[mi][mj]*f[Solid[i][1][k]][mj];
					}

			m_l[0]=rho_yn;
			for (int mi=0; mi<19; mi++)
				{
				F[Solid[i][0][k]][mi]=0;
				for (int mj=0; mj<19; mj++)
					F[Solid[i][0][k]][mi]+=MI[mi][mj]*m_l[mj];
				}
			}
			else
			{
			u_ls[0]=0.0;u_ls[1]=0.0;u_ls[2]=0.0;
			for (int ks=0;ks<Q;ks++)
				F[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls);
			
			}




if (zp==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][NZ]>0)
		if (Solid[i][j][NZ-1]>0)
			{
			for (int mi=0; mi<19; mi++)
					{
					m_l[mi]=0;
					for (int mj=0; mj<19; mj++)
						m_l[mi]+=M[mi][mj]*f[Solid[i][j][NZ-1]][mj];
					}

			m_l[0]=rho_zp;
			for (int mi=0; mi<19; mi++)
				{
				F[Solid[i][j][NZ]][mi]=0;
				for (int mj=0; mj<19; mj++)
					F[Solid[i][j][NZ]][mi]+=MI[mi][mj]*m_l[mj];
				}
			}
			else
			{
			u_ls[0]=0.0;u_ls[1]=0.0;u_ls[2]=0.0;
			for (int ks=0;ks<Q;ks++)
				F[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls);
			
			}




if (zn==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][0]>0)
		if (Solid[i][j][1]>0)
			{
			for (int mi=0; mi<19; mi++)
					{
					m_l[mi]=0;
					for (int mj=0; mj<19; mj++)
						m_l[mi]+=M[mi][mj]*f[Solid[i][j][1]][mj];
					}

			m_l[0]=rho_zn;
			for (int mi=0; mi<19; mi++)
				{
				F[Solid[i][j][0]][mi]=0;
				for (int mj=0; mj<19; mj++)
					F[Solid[i][j][0]][mi]+=MI[mi][mj]*m_l[mj];
				}
			}
			else
			{
			u_ls[0]=0.0;u_ls[1]=0.0;u_ls[2]=0.0;
			for (int ks=0;ks<Q;ks++)
				F[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls);
			
			}



if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
	if (Solid[nx_l-1][j][k]>0)
		if (Solid[nx_l-2][j][k]>0)
			{
			for (int mi=0; mi<19; mi++)
					{
					m_l[mi]=0;
					for (int mj=0; mj<19; mj++)
						m_l[mi]+=M[mi][mj]*f[Solid[nx_l-2][j][k]][mj];
					}

			m_l[0]=rho_xp;
			for (int mi=0; mi<19; mi++)
				{
				F[Solid[nx_l-1][j][k]][mi]=0;
				for (int mj=0; mj<19; mj++)
					F[Solid[nx_l-1][j][k]][mi]+=MI[mi][mj]*m_l[mj];
				}
			}
			else
			{
			u_ls[0]=0.0;u_ls[1]=0.0;u_ls[2]=0.0;
			for (int ks=0;ks<Q;ks++)
				F[Solid[nx_l-1][j][NZ]][ks]=feq(ks,rho_xp,u_ls);
			
			}
		
		


if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
	if (Solid[0][j][k]>0)
		if (Solid[1][j][k]>0)
			{
			for (int mi=0; mi<19; mi++)
					{
					m_l[mi]=0;
					for (int mj=0; mj<19; mj++)
						m_l[mi]+=M[mi][mj]*f[Solid[1][j][k]][mj];
					}

			m_l[0]=rho_xn;
			for (int mi=0; mi<19; mi++)
				{
				F[Solid[0][j][k]][mi]=0;
				for (int mj=0; mj<19; mj++)
					F[Solid[0][j][k]][mi]+=MI[mi][mj]*m_l[mj];
				}
			}
			else
			{
			u_ls[0]=0.0;u_ls[1]=0.0;u_ls[2]=0.0;
			for (int ks=0;ks<Q;ks++)
				F[Solid[0][j][NZ]][ks]=feq(ks,rho_xn,u_ls);
			
			}
	


}		

//Non-equilibrium boundary condition
if (Sub_BC==3)
{
if (yp==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)	
	if (Solid[i][NY][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
				if (Solid[i][NY-1][k]>0)
					{
					u_ls[0]=u[Solid[i][NY-1][k]][0];
					u_ls[1]=u[Solid[i][NY-1][k]][1];
					u_ls[2]=u[Solid[i][NY-1][k]][2];
					F[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls)+f[Solid[i][NY-1][k]][ks]-feq(ks,rho[Solid[i][NY-1][k]],u[Solid[i][NY-1][k]]);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls);
					}
			}


if (yn==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)	
	if (Solid[i][0][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
				if (Solid[i][1][k]>0)
					{
					u_ls[0]=u[Solid[i][1][k]][0];
					u_ls[1]=u[Solid[i][1][k]][1];
					u_ls[2]=u[Solid[i][1][k]][2];
					F[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls)+f[Solid[i][1][k]][ks]-feq(ks,rho[Solid[i][1][k]],u[Solid[i][1][k]]);
					}
			else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls);
					}
			}

if (zp==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][NZ]>0)
		for (int ks=0;ks<Q;ks++)
			{			
			if (Solid[i][j][NZ-1]>0)
				{
				u_ls[0]=u[Solid[i][j][NZ-1]][0];
				u_ls[1]=u[Solid[i][j][NZ-1]][1];
				u_ls[2]=u[Solid[i][j][NZ-1]][2];
				F[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls)+f[Solid[i][j][NZ-1]][ks]-feq(ks,rho[Solid[i][j][NZ-1]],u[Solid[i][j][NZ-1]]);
				}
				else
				{
				u_ls[0]=0.0;
				u_ls[1]=0.0;u_ls[2]=0.0;
				F[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls);
				}
			}




if (zn==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][0]>0)
		for (int ks=0;ks<Q;ks++)
		{
			if (Solid[i][j][1]>0)
			{
			u_ls[0]=u[Solid[i][j][1]][0];
			u_ls[1]=u[Solid[i][j][1]][1];
			u_ls[2]=u[Solid[i][j][1]][2];
			F[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls)+f[Solid[i][j][1]][ks]-feq(ks,rho[Solid[i][j][1]],u[Solid[i][j][1]]);
			}
			else
			{
			u_ls[0]=0.0;
			u_ls[1]=0.0;u_ls[2]=0.0;
			F[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls);
			}
		
		}

if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
	if (Solid[nx_l-1][j][k]>0)
		for (int ks=0;ks<Q;ks++)			
			{
				if (Solid[nx_l-2][j][k]>0)
					{
					u_ls[0]=u[Solid[nx_l-2][j][k]][0];
					u_ls[1]=u[Solid[nx_l-2][j][k]][1];
					u_ls[2]=u[Solid[nx_l-2][j][k]][2];
					F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls)+f[Solid[nx_l-2][j][k]][ks]-feq(ks,rho[Solid[nx_l-2][j][k]],u[Solid[nx_l-2][j][k]]);
					//F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls);
					}
				else	
					{
					u_ls[0]=0.0;u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls);
					}

			}
		


if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
	if (Solid[0][j][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
				if(Solid[1][j][k]>0)
					{
					u_ls[0]=u[Solid[1][j][k]][0];
					u_ls[1]=u[Solid[1][j][k]][1];
					u_ls[2]=u[Solid[1][j][k]][1];
					F[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls)+f[Solid[1][j][k]][ks]-feq(ks,rho[Solid[1][j][k]],u[Solid[1][j][k]]);
					//f[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;
					u_ls[2]=0.0;
					F[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls);
					}

			}

	
	
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

	//u_compt/=(NX+1)*(NY+1)*(NZ+1);
	u_compt/=(NX+1)*(NY+1)*(NZ+1)*porosity;
	*u_average=u_compt;
	
	delete [] rbuf;
	delete [] um;
	delete [] uave;

	//MPI_Barrier(MPI_COMM_WORLD);

	return(error_in);

}




//OUTPUT SUBROUTAINS:
//ALL THE OUTPUTS ARE TRANSFERED TO PROCESSOR 0, AND EXPORT TO DAT FILE BY PROCESSOR 0


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
			out<<rbuf[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]<<endl;
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

//===================================================================
/*	
	ostringstream name2;
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
//=================================================================	
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
				out<<rbuf_rho[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]<<endl;

	out.close();
				
	}
	
	
	if (rank==root_rank)
		{		
		delete [] rbuf_rho;
		}
	delete [] nx_g;
	delete [] disp;
	delete [] rho_storage;

		
}


void Backup(int m,double* rho,double** u, double** f)
{

        int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	
	ostringstream name;
	name<<outputfile<<"LBM_Backup_Velocity_"<<m<<"."<<rank<<".input";
	ofstream out;
	out.open(name.str().c_str());
	for (int i=1;i<=Count;i++)
        		out<<u[i][0]<<" "<<u[i][1]<<" "<<u[i][2]<<" "<<endl;
		
			
	out.close();
	

	
	
	ostringstream name2;
	name2<<outputfile<<"LBM_Backup_Density_"<<m<<"."<<rank<<".input";
	ofstream out2;
	out2.open(name2.str().c_str());
	
	for (int i=1;i<=Count;i++)
        		out2<<rho[i]<<endl;
		
			
	out2.close();
	
	
	
	ostringstream name4;
	name4<<outputfile<<"LBM_Backup_f_"<<m<<"."<<rank<<".input";
	//ofstream out;
	out.open(name4.str().c_str());
	
	for (int i=1;i<=Count;i++)
	{
	        for (int j=0;j<19;j++)
        		out<<f[i][j]<<" ";
        out<<endl;
        }
                
        out.close();
        
	
	
}

double Comput_Perm(double** u,double* Permia,int PerDIr,int* SupInv)
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	int nx_g[mpi_size];
	int disp[mpi_size];
	int si,sj,sm;
	
	double avex,avey,avez;
	
	MPI_Gather(&nx_l,1,MPI_INT,nx_g,1,MPI_INT,0,MPI_COMM_WORLD);
	
	
	if (rank==0)
		{
		disp[0]=0;
	
		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
		
		}

	MPI_Bcast(disp,mpi_size,MPI_INT,0,MPI_COMM_WORLD);




	double *rbuf;
	rbuf=new double[mpi_size*3];
	double Perm[3];
	double error;
	double Q[3]={0.0,0.0,0.0};

	double dp;
	if (in_BC==0)
	        dp=0;
	else
	switch(PerDIr)
		{
		case 1:
			dp=abs(p_xp-p_xn)*c_s2/(NX+1)/dx;break;
		case 2:
			dp=abs(p_yp-p_yn)*c_s2/(NY+1)/dx;break;
		case 3:
			dp=abs(p_zp-p_zn)*c_s2/(NZ+1)/dx;break;
		default:
			dp=abs(p_xp-p_xn)*c_s2/(NX+1)/dx;
		}

		
	if ((par_per_x-1)*(par_per_y-1)*(par_per_z-1)==0)	
	for (int i=1;i<=Count;i++)
	{
		si=(int)(SupInv[i]/((NY+1)*(NZ+1)));
		sj=(int)((SupInv[i]%((NY+1)*(NZ+1)))/(NZ+1));
		sm=(int)(SupInv[i]%(NZ+1)); 
		si+=disp[rank];
		//if (rank==1)
		//cout<<rank<<"  "<<si<<" "<<sj<<" "<<sm<<endl;
		//cout<<si<<"  "<<per_xn<<"  "<<per_xp<<endl;
		if ((si>=per_xn) and (si<=per_xp) and (sj>=per_yn) and (sj<=per_yp) and (sm>=per_zn) and (sm<=per_zp))
		{
	        Q[0]+=u[i][0];
		Q[1]+=u[i][1];
		Q[2]+=u[i][2];
		}


	}
	else
	for (int i=1;i<=Count;i++)
	       {
		Q[0]+=u[i][0];
		Q[1]+=u[i][1];
		Q[2]+=u[i][2];
		} 

	
		
//	for (int i=1;i<=Count;i++)
//		{
//		Q[0]+=u[i][0];
//		Q[1]+=u[i][1];
//		Q[2]+=u[i][2];
//		}

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

		Perm[0]=Q[0]/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1))*(in_vis)/(gx+dp);
		Perm[1]=Q[1]/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1))*(in_vis)/(gy+dp);
		Perm[2]=Q[2]/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1))*(in_vis)/(gz+dp);
		
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

		avex=Q[0]/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1)*porosity);
		avey=Q[1]/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1)*porosity);
		avez=Q[2]/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1)*porosity);

		u_ave=sqrt(avex*avex+avey*avey+avez*avez);


		}
	
	delete [] rbuf;
	
	return (error);
	

}



void Backup_init(double* rho, double** u, double** f, char backup_rho[128], char backup_velocity[128], char backup_f[128])
{	
      
        int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	
	
	
	double usqr,vsqr,ls_rho,ls_v0,ls_v1,ls_v2;

	
	rho0=1.0;dt=1.0/Zoom;dx=1.0/Zoom;
 	uMax=0.0;

	if (lattice_v==1)
		{dx=dx_input;dt=dt_input;}

	lat_c=dx/dt;
	c_s=lat_c/sqrt(3);
	c_s2=lat_c*lat_c/3;

	niu=in_vis;
	tau_f=niu/(c_s2*dt)+0.5;
	//tau_f=3.0*niu/dt+0.5;
	
	s_v=1/tau_f;
       
	double s_other=8*(2-s_v)/(8-s_v);
	
	if (lattice_v==1)
	for (int i=0;i<19;i++)
		for (int j=0;j<3;j++)
		elat[i][j]=e[i][j]*lat_c;

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

	ostringstream name4;
	name4<<backup_f<<"."<<rank<<".input";
	
 	ostringstream name2;
	name2<<backup_velocity<<"."<<rank<<".input";
	ostringstream name;
	name<<backup_rho<<"."<<rank<<".input";
	
	 
	ifstream fin;
	fin.open(name.str().c_str());
	
        	for(int i=1;i<=Count;i++)
        	        fin >> rho[i];
  
       fin.close();
       
   
	fin.open(name2.str().c_str());
	
        	for(int i=1;i<=Count;i++)
        	        fin >> u[i][0] >> u[i][1] >> u[i][2];
  
       fin.close();
       
      
       fin.open(name4.str().c_str());
	
        	for(int i=1;i<=Count;i++)
        	        fin >> f[i][0] >> f[i][1] >> f[i][2] >> f[i][3] >> f[i][4] >> f[i][5] >>f[i][6] >> f[i][7] >> f[i][8] >> f[i][9] >> f[i][10] >> f[i][11] >> f[i][12] >> f[i][13] >> f[i][14] >> f[i][15] >> f[i][16] >>f[i][17] >> f[i][18];
  
       fin.close();
       
      
	

	 	
}




void Comput_Grop_Perm(double** u,double* Permia,int PerDIr,int* SupInv)
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	int nx_g[mpi_size];
	int disp[mpi_size];
	int si,sj,sm;
	int per_xn1,per_xp1,per_yn1,per_yp1,per_zn1,per_zp1;
	int General_size=512;
	int sub_size,por_loc,vol_g;
	char File[128];
	strcpy(File,outputfile);

	//========================================
	int loop[4]={1,1,2,4};
	int loop_size[4]={512,256,128,64};
	int loop_s[4][3]={{0,1,1},{128,129,129},{128,129,129},{128,129,129}};
	//========================================


	MPI_Gather(&nx_l,1,MPI_INT,nx_g,1,MPI_INT,0,MPI_COMM_WORLD);

	
	
	if (rank==0)
		{
		disp[0]=0;
	
		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
		
		}

	MPI_Bcast(disp,mpi_size,MPI_INT,0,MPI_COMM_WORLD);



	double *rbuf;
	int *rbuf_p;
	rbuf=new double[mpi_size*3];
	rbuf_p=new int[mpi_size];
	double Perm[3];
	double error;
	double Q[3]={0.0,0.0,0.0};

	double dp;
	if (in_BC==0)
	        dp=0;
	else
	switch(PerDIr)
		{
		case 1:
			dp=abs(p_xp-p_xn)*c_s2/(NX+1)/dx;break;
		case 2:
			dp=abs(p_yp-p_yn)*c_s2/(NY+1)/dx;break;
		case 3:
			dp=abs(p_zp-p_zn)*c_s2/(NZ+1)/dx;break;
		default:
			dp=abs(p_xp-p_xn)*c_s2/(NX+1)/dx;
		}

ostringstream name;
for (int divi=0;divi<4;divi++)	
{

if (rank==0)
{
name.str(""); 
name<<File<<"Perm_"<<loop_size[divi]<<".output";
	
//-----------------------------------------	 
//	ifstream fin;
//	fin.open(name.str().c_str());
//	
  //      	for(int i=1;i<=Count;i++)
    //    	        fin >> rho[i];
  
      // fin.close();

//ofstream out;
//	out.open(name.str().c_str());

ofstream fin(name.str().c_str(),ios::out);
fin.close();
//ofstream fin(FileName,ios::app);
//----------------------------------------------------
}


vol_g=loop_size[divi]*loop_size[divi]*loop_size[divi];


for (int in_x=0;in_x<loop[divi];in_x++)
for (int in_y=0;in_y<loop[divi];in_y++)
for (int in_z=0;in_z<loop[divi];in_z++)

{
	
	per_xn1=loop_s[divi][0]+in_x*loop_size[divi];
	per_yn1=loop_s[divi][1]+in_y*loop_size[divi];
	per_zn1=loop_s[divi][2]+in_z*loop_size[divi];

	per_xp1=loop_s[divi][0]+in_x*loop_size[divi]+loop_size[divi]-1;
	per_yp1=loop_s[divi][1]+in_y*loop_size[divi]+loop_size[divi]-1;
	per_zp1=loop_s[divi][2]+in_z*loop_size[divi]+loop_size[divi]-1;


	por_loc=0;
	Q[0]=0;Q[1]=0;Q[2]=0;
	for (int i=1;i<=Count;i++)
	{
		si=(int)(SupInv[i]/((NY+1)*(NZ+1)));
		sj=(int)((SupInv[i]%((NY+1)*(NZ+1)))/(NZ+1));
		sm=(int)(SupInv[i]%(NZ+1)); 
		si+=disp[rank];
		

		if ((si>=per_xn1) and (si<=per_xp1) and (sj>=per_yn1) and (sj<=per_yp1) and (sm>=per_zn1) and (sm<=per_zp1))
		{
	        Q[0]+=u[i][0];
		Q[1]+=u[i][1];
		Q[2]+=u[i][2];
		por_loc+=1;
		}


	}
	
	
		


	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Gather(&por_loc,1,MPI_INT,rbuf_p,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Gather(&Q,3,MPI_DOUBLE,rbuf,3,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	if (rank==0)
		{
		Q[0]=0;Q[1]=0;Q[2]=0;por_loc=0;
		for (int i=0;i<mpi_size;i++)
			{
			Q[0]+=rbuf[i*3+0];
			Q[1]+=rbuf[i*3+1];
			Q[2]+=rbuf[i*3+2];
			por_loc+=rbuf_p[i];
			}

	
		Perm[0]=Q[0]/(vol_g)*(in_vis)/(gx+dp);
		Perm[1]=Q[1]/(vol_g)*(in_vis)/(gy+dp);
		Perm[2]=Q[2]/(vol_g)*(in_vis)/(gz+dp);
		
		

		}
	MPI_Barrier(MPI_COMM_WORLD);


	//================INDICATED PERMEABILITY CALCULATION DIRECTION HERE!!!!!====================
	if (rank==0)
	{
	ofstream fin(name.str().c_str(),ios::app);
	fin<<Perm[0]*reso*reso*1000<<" "<< (double)por_loc/vol_g <<endl;
	fin.close();
	}



	}
}
	
	delete [] rbuf;
	delete [] rbuf_p;
	
	
	

}

