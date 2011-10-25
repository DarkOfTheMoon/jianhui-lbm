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
{-1,0,0,0,0,0,0,1,1,1,1,1,1,1,1,1,1,1,1},
{1,-2,-2,-2,-2,-2,-2,1,1,1,1,1,1,1,1,1,1,1,1},
{0,1,-1,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0},
{0,-2,2,0,0,0,0,1,-1,1,-1,1,-1,1,-1,0,0,0,0},
{0,0,0,1,-1,0,0,1,-1,-1,1,0,0,0,0,1,-1,1,-1},
{0,0,0,-2,2,0,0,1,-1,-1,1,0,0,0,0,1,-1,1,-1},
{0,0,0,0,0,1,-1,0,0,0,0,1,-1,-1,1,1,-1,-1,1},
{0,0,0,0,0,-2,2,0,0,0,0,1,-1,-1,1,1,-1,-1,1},
{0,2,2,-1,-1,-1,-1,1,1,1,1,1,1,1,1,-2,-2,-2,-2},
{0,-2,-2,1,1,1,1,1,1,1,1,1,1,1,1,-2,-2,-2,-2},
{0,0,0,1,1,-1,-1,1,1,1,1,-1,-1,-1,-1,0,0,0,0},
{0,0,0,-1,-1,1,1,1,1,1,1,-1,-1,-1,-1,0,0,0,0},
{0,0,0,0,0,0,0,1,1,-1,-1,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,-1,-1},
{0,0,0,0,0,0,0,0,0,0,0,1,1,-1,-1,0,0,0,0},
{0,0,0,0,0,0,0,1,-1,1,-1,-1,1,-1,1,0,0,0,0},
{0,0,0,0,0,0,0,-1,1,1,-1,0,0,0,0,1,-1,1,-1},
{0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1,-1,1,1,-1}};


double MI[19][19];

double M_c[19];






double m[19];
double meq[19];





double uMax,c,Re,dx,dy,Lx,Ly,dt,rho0,P0,tau_f,niu,error,SFx,SFy,reso;

void Read_Rock(int***,double***, double***, double*,char[128],char[128],char[128]);

void tests();

void init_Sparse(int***,int***,double***, double***, double***, double***, int*, int*);

void init(double*, double**, double**,double*, double*,double*, double*, double*,double***,double***, int*);

void standard_bounceback_boundary(int,double**);

void collision(double*,double** ,double** ,double**, double*, double*, double*, double* ,double* , int* ,int***,int* ,int*);

void comput_macro_variables( double* ,double**,double** ,double** ,double** ,double*, double*, double*, double*, double* ,int* ,int***,double***);

double Error(double** ,double** ,double*, double*);

void boundary_velocity(int,double,int , double ,int ,double ,int ,double ,int ,double ,int ,double ,double** ,double**,double* ,double** ,int*** );

void boundary_pressure(int ,double ,int , double ,int ,double ,int ,double ,int ,double ,int ,double ,double** ,double**,double** ,double* ,int*** );

void output_velocity(int ,double* ,double** ,int ,int ,int ,int ,int*** );

void output_density(int ,double* ,int ,int ,int ,int ,int*** );	

void Geometry(int*** );	

void output_velocity_b(int ,double* ,double** ,int ,int ,int ,int ,int*** );

void output_psi(int ,double* ,int ,int ,int ,int ,int***,char[128] );

void output_psi_b(int ,double* ,int ,int ,int ,int ,int*** ,char[128]);

void output_density_b(int ,double* ,int ,int ,int ,int ,int*** );	

void Geometry_b(int*** );

void Solute_Constant_BC(int,double,int, double,int ,double ,int ,double ,int ,double ,int,double ,double**,int***,double*,double**,double**);

void Solute_ZeroFlux_BC(int,double, int,double,int,double,int,double,int,double,int,double,double**,int***,double*,double**,double**);


double Comput_Perm(double* ,double** ,double* ,double* ,int );

double Comput_Saturation(double* ,int***);

double S[19];

void Comput_MI(double[19][19], double[19][19]);

int inverse(mat &a);

double feq(int,double, double[3]);

double feq_psi(int,double, double[3]);

void Suppliment(int*,int***);

void Backup_init(double* , double** , double** ,double* ,double* , double* , double*, double* ,char[128],char[128],char[128], char[128], char[128], char[128]);

void Backup(int ,double*, double*, double**, double**);




int e[19][3]=
{{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},{1,0,1},
{-1,0,-1},{1,0,-1},{-1,0,1},{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1}};

double elat[19][3]=
{{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},{1,0,1},
{-1,0,-1},{1,0,-1},{-1,0,1},{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1}};


double w[19]={1.0/3.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};


//========*************************===============
//int LR[19]={0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15};
int LR[19]={0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17};
int FRP[19]={0,0,0,0,0,0,0,1,0,2,0,3,0,4,0,0,0,0,0};
int FLN[19]={0,0,0,0,0,0,0,0,1,0,2,0,3,0,4,0,0,0,0};
int RP[5]={1,7,9,11,13};
int LN[5]={2,8,10,12,14};
int SYN[5]={3,7,10,15,17};
int SYP[5]={4,8,9,16,18};
int SZN[5]={5,11,14,15,18};
int SZP[5]={6,12,13,16,17};
int SXN[5]={1,7,9,11,13};
int SXP[5]={2,8,10,12,14};
//=========================================


int n,nx_l,n_max,in_BC,PerDir,freRe,freDe,freVe,frePsi,Par_Geo,Par_nx,Par_ny,Par_nz;
int Zoom,lattice_v,freHS,Sub_BC_psi;


int wr_per,pre_xp,pre_xn,pre_yp,pre_yn,pre_zp,pre_zn,stab,stab_time,fre_backup;
int vel_xp,vel_xn,vel_yp,vel_yn,vel_zp,vel_zn,Sub_BC,Out_Mode,mode_backup_ini;
double in_vis,p_xp,p_xn,p_yp,p_yn,p_zp,p_zn,niu_l,niu_g,ContactAngle_parameter,CapA;
double inivx,inivy,inivz,v_xp,v_xn,v_yp,v_yn,v_zp,v_zn,Re_l,Re_g,Capillary;
double error_Per,Permeability,psi_solid,S_l,gxs,gys,gzs,c_s,c_s2,dx_input,dt_input,lat_c;
double niu_s1,niu_s2,tau_s,Gravity,ref_psi;
double Buoyancy_parameter=1.0;

int sol_c_xp,sol_c_xn,sol_c_yp,sol_c_yn,sol_c_zp,sol_c_zn;
double c_xp,c_xn,c_yp,c_yn,c_zp,c_zn;

int sol_zf_xp,sol_zf_xn,sol_zf_yp,sol_zf_yn,sol_zf_zp,sol_zf_zn;
double zf_xp,zf_xn,zf_yp,zf_yn,zf_zp,zf_zn;


double** fg;
double** Fg;
double* rhoh;


char outputfile[128]="./";

int main(int argc , char *argv [])
{	

MPI :: Init (argc , argv );
MPI_Status status ;

double start , finish,remain;

int rank = MPI :: COMM_WORLD . Get_rank ();
int para_size=MPI :: COMM_WORLD . Get_size ();

int dif,th,tm,ts;

int tse,the,tme;
double elaps;

double Per_l[3],Per_g[3];
double v_max,error_Per;




	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	int NCHAR=128;
	char     filename[128], dummy[128+1],filenamepsi[128], filenamepsi2[128], backup_rho[128], backup_velocity[128],backup_psi[128],backup_f[128],backup_fg[128], backup_hs[128];
	int      dummyInt;

	if (rank==0)
	{
	ifstream fin(argv[1]);
	                                                fin.getline(dummy, NCHAR);
	fin >> filename;				fin.getline(dummy, NCHAR);
	fin >> filenamepsi;				fin.getline(dummy, NCHAR);
	fin >> filenamepsi2;				fin.getline(dummy, NCHAR);
	fin >> NX >> NY >> NZ;				fin.getline(dummy, NCHAR);
	fin >> n_max;					fin.getline(dummy, NCHAR);
	fin >> reso;					fin.getline(dummy, NCHAR);
	fin >> in_BC;					fin.getline(dummy, NCHAR);
	fin >> Buoyancy_parameter >>ref_psi;		fin.getline(dummy, NCHAR);
	fin >> mirX >> mirY >> mirZ;			fin.getline(dummy, NCHAR);
	fin >> Gravity;					fin.getline(dummy, NCHAR);
	fin >> gx >> gy >> gz;				fin.getline(dummy, NCHAR);
	fin >> pre_xp >> p_xp >> pre_xn >> p_xn;	fin.getline(dummy, NCHAR);
	fin >> pre_yp >> p_yp >> pre_yn >> p_yn;	fin.getline(dummy, NCHAR);
	fin >> pre_zp >> p_zp >> pre_zn >> p_zn;	fin.getline(dummy, NCHAR);
	fin >> vel_xp >> v_xp >> vel_xn >> v_xn;	fin.getline(dummy, NCHAR);
	fin >> vel_yp >> v_yp >> vel_yn >> v_yn;	fin.getline(dummy, NCHAR);
	fin >> vel_zp >> v_zp >> vel_zn >> v_zn;	fin.getline(dummy, NCHAR);
	fin >> sol_c_xp >> c_xp >> sol_c_xn >> c_xn;	fin.getline(dummy, NCHAR);
	fin >> sol_c_yp >> c_yp >> sol_c_yn >> c_yn;	fin.getline(dummy, NCHAR);
	fin >> sol_c_zp >> c_zp >> sol_c_zn >> c_zn;	fin.getline(dummy, NCHAR);
	fin >> sol_zf_xp >> zf_xp >> sol_zf_xn >> zf_xn;	fin.getline(dummy, NCHAR);
	fin >> sol_zf_yp >> zf_yp >> sol_zf_yn >> zf_yn;	fin.getline(dummy, NCHAR);
	fin >> sol_zf_zp >> zf_zp >> sol_zf_zn >> zf_zn;	fin.getline(dummy, NCHAR);
	fin >> niu_l;					fin.getline(dummy, NCHAR);
	fin >> niu_g;					fin.getline(dummy, NCHAR);
	fin >> niu_s1;					fin.getline(dummy, NCHAR);
	fin >> niu_s2;					fin.getline(dummy, NCHAR);
	fin >> ContactAngle_parameter;			fin.getline(dummy, NCHAR);
	fin >> CapA;					fin.getline(dummy, NCHAR);
	fin >> inivx >> inivy >> inivz;			fin.getline(dummy, NCHAR);
	fin >> Permeability;				fin.getline(dummy, NCHAR);
							fin.getline(dummy, NCHAR);
	fin >> wr_per;					fin.getline(dummy, NCHAR);
	fin >> PerDir;					fin.getline(dummy, NCHAR);
	fin >> freRe;					fin.getline(dummy, NCHAR);
	fin >> Out_Mode;				fin.getline(dummy, NCHAR);
	fin >> freVe;					fin.getline(dummy, NCHAR);
	fin >> freDe;					fin.getline(dummy, NCHAR);
	fin >> frePsi;					fin.getline(dummy, NCHAR);
	fin >> freHS;					fin.getline(dummy, NCHAR);
	fin >> mir;					fin.getline(dummy, NCHAR);
							fin.getline(dummy, NCHAR);
	fin >> Par_Geo >> Par_nx >> Par_ny >> Par_nz;	fin.getline(dummy, NCHAR);
	fin >> Zoom;					fin.getline(dummy, NCHAR);
	fin >> lattice_v >> dx_input >> dt_input;	fin.getline(dummy, NCHAR);
	fin >> outputfile;				fin.getline(dummy, NCHAR);
	fin >> Sub_BC;					fin.getline(dummy, NCHAR);
	fin >> Sub_BC_psi;				fin.getline(dummy, NCHAR);
	fin >> stab >> stab_time;			fin.getline(dummy, NCHAR);
	fin >> fre_backup;                        fin.getline(dummy, NCHAR);
	fin >>mode_backup_ini;                fin.getline(dummy, NCHAR);
	fin >> backup_rho;                        fin.getline(dummy, NCHAR);
	fin >> backup_velocity;                fin.getline(dummy, NCHAR);
	fin >> backup_psi;                        fin.getline(dummy, NCHAR);
	fin >> backup_f;                        fin.getline(dummy, NCHAR);
	fin >> backup_fg;                        fin.getline(dummy, NCHAR);
	fin >> backup_hs;                        fin.getline(dummy, NCHAR);
	fin.close();
	
	//cout<<q_p<<"    asdfa "<<endl;
	NX=NX-1;NY=NY-1;NZ=NZ-1;
	}

	MPI_Bcast(&filename,128,MPI_CHAR,0,MPI_COMM_WORLD);
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
	MPI_Bcast(&mir,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&CapA,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Zoom,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&outputfile,128,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&Sub_BC,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&Out_Mode,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&niu_l,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&niu_g,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
      	MPI_Bcast(&ContactAngle_parameter,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Permeability,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&frePsi,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&stab,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&stab_time,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&fre_backup,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&backup_rho,128,MPI_CHAR,0,MPI_COMM_WORLD);MPI_Bcast(&backup_velocity,128,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&mode_backup_ini,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&backup_psi,128,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&backup_f,128,MPI_CHAR,0,MPI_COMM_WORLD);

	MPI_Bcast(&lattice_v,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&dx_input,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&dt_input,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&niu_s1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&ref_psi,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&niu_s2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Buoyancy_parameter,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Gravity,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&freHS,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&backup_fg,128,MPI_CHAR,0,MPI_COMM_WORLD);MPI_Bcast(&backup_hs,128,MPI_CHAR,0,MPI_COMM_WORLD);

	MPI_Bcast(&sol_c_xp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&c_xp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sol_c_xn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&c_xn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sol_c_yp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&c_yp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sol_c_yn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&c_yn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sol_c_zp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&c_zp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sol_c_zn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&c_zn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sol_zf_xp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&zf_xp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sol_zf_xn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&zf_xn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sol_zf_yp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&zf_yp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sol_zf_yn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&zf_yn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sol_zf_zp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&zf_zp,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&sol_zf_zn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&zf_zn,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(&Sub_BC_psi,1,MPI_INT,0,MPI_COMM_WORLD);
	
int U_max_ref=0;


if (mirX==1)
	NX=NX*2+1;
if (mirY==1)
	NY=NY*2+1;
if (mirZ==1)
	NZ=NZ*2+1;


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
	double* rho_r;
	double* rho_b;
	double* rhor;
	double* rhob;
	double** u;
	double**f;
	double**F;
	double**u0;
	int* SupInv;
	//double* forcex;
	//double* forcey;
	//double* forcez;
	
	double* psi;
	double*** Psi_local;
	double*** Psis;

	int*** Solids;
	int*** Solid;

	double*** Psi_local2;
	double*** Psis2;


	
	int*  Sl;
	int*  Sr;

	
	
	Solid = new int**[nx_l];
	Psi_local = new double**[nx_l];
	Sl = new int[(NY+1)*(NZ+1)];
	Sr = new int[(NY+1)*(NZ+1)];

	Solids = new int**[(NX+1)/Zoom];
	Psis = new double**[(NX+1)/Zoom];

	Psi_local2 = new double**[nx_l];
	Psis2 = new double**[(NX+1)/Zoom];

	for (int i=0;i<(NX+1)/Zoom;i++)
		{		
		Solids[i] = new int*[(NY+1)/Zoom];
		Psis[i] = new double*[(NY+1)/Zoom];
		Psis2[i] = new double*[(NY+1)/Zoom];
			for (int j=0;j<(NY+1)/Zoom;j++)
			{
			Solids[i][j]= new int[(NZ+1)/Zoom];
			Psis[i][j]= new double[(NZ+1)/Zoom];
			Psis2[i][j]= new double[(NZ+1)/Zoom];

			}
		}


	for (int i=0;i<nx_l;i++)
		{
		Solid[i] = new int*[NY+1];
		Psi_local[i] = new double*[NY+1];
		Psi_local2[i] = new double*[NY+1];
			for (int j=0;j<=NY;j++)
			{
			Solid[i][j]= new int[NZ+1];
			Psi_local[i][j] = new double[NZ+1];
			Psi_local2[i][j] = new double[NZ+1];
			}
		}


		Read_Rock(Solids,Psis,Psis2,&porosity,filename,filenamepsi,filenamepsi2);


	init_Sparse(Solids,Solid,Psis,Psi_local,Psis2,Psi_local2,Sl,Sr);

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
	rhoh = new double[Count+1];
	rho_r = new double[Count+1];
	rho_b = new double[Count+1];
	rhor = new double[Count+1];
	rhob = new double[Count+1];
	psi = new double[Count+1];
	//forcex = new double[Count+1];
	//forcey = new double[Count+1];
	//forcez = new double[Count+1];
	u = new double*[Count+1];
	f = new double*[Count+1];
	F = new double*[Count+1];
	fg = new double*[Count+1];
	Fg = new double*[Count+1];
	u0 = new double*[Count+1];
	SupInv = new int[Count+1];

	for (int i=0;i<=Count;i++)
		{
		u[i] = new double[3];
		f[i] = new double[19];
		fg[i] = new double[19];
		u0[i] = new double[3];
		F[i] = new double[19];
		Fg[i] = new double[19];
		}

	Comput_MI(M,MI);
	
	Suppliment(SupInv,Solid);

	MPI_Barrier(MPI_COMM_WORLD);
	
	if (Out_Mode==1)
		Geometry(Solid);
	else
		Geometry_b(Solid);

	if (mode_backup_ini==0)
	        init(rho,u,f,psi,rho_r,rho_b,rhor, rhob, Psi_local,Psi_local2,SupInv);
	else
	      Backup_init( rho, u, f,psi,rho_r, rho_b, rhor, rhob, backup_rho, backup_velocity, backup_psi,backup_f,backup_fg,backup_hs);  

if (rank==0)
		cout<<"Porosity= "<<porosity<<endl;




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

	
	
	for(n=0;n<=n_max;n++)
	{
	
	if ((stab==1) and (n==stab_time))
		{gxs=gx;gys=gy;gzs=gz;}
			

	collision(rho,u,f,F,psi,rho_r,rho_b,rhor,rhob,SupInv,Solid,Sl,Sr);

	
	if ((1-pre_xp)*(1-pre_xn)*(1-pre_yp)*(1-pre_yn)*(1-pre_zp)*(1-pre_zn)==0)
		boundary_pressure(pre_xp,p_xp,pre_xn,p_xn,pre_yp,p_yp,pre_yn,p_yn,pre_zp,p_zp,pre_zn,p_zn,f,F,u,rho,Solid);

	
	if ((1-vel_xp)*(1-vel_xn)*(1-vel_yp)*(1-vel_yn)*(1-vel_zp)*(1-vel_zn)==0)
		boundary_velocity(vel_xp,v_xp,vel_xn,v_xn,vel_yp,v_yp,vel_yn,v_yn,vel_zp,v_zp,vel_zn,v_zn,f,F,rho,u,Solid);

  	if ((1-sol_c_xp)*(1-sol_c_xn)*(1-sol_c_yp)*(1-sol_c_yn)*(1-sol_c_zp)*(1-sol_c_zn)==0)
	Solute_Constant_BC(sol_c_xp,c_xp,sol_c_xn,c_xn,sol_c_yp,c_yp,sol_c_yn,c_yn,sol_c_zp,c_zp,sol_c_zn,c_zn,Fg,Solid,rhoh,u,fg);
	
	if ((1-sol_zf_xp)*(1-sol_zf_xn)*(1-sol_zf_yp)*(1-sol_zf_yn)*(1-sol_zf_zp)*(1-sol_zf_zn)==0)
	Solute_ZeroFlux_BC(sol_zf_xp,zf_xp,sol_zf_xn,zf_xn,sol_zf_yp,zf_yp,sol_zf_yn,zf_yn,sol_zf_zp,zf_zp,sol_zf_zn,zf_zn,Fg,Solid,rhoh,u,fg);
	 
		comput_macro_variables(rho,u,u0,f,F,rho_r,rho_b,rhor,rhob,psi,SupInv,Solid,Psi_local);

	
	if(n%freRe==0)
		{       
			error=Error(u,u0,&u_max,&u_ave);
			if (u_max>=10.0)	U_max_ref+=1;
			error_Per=Comput_Perm(psi,u,Per_l,Per_g,PerDir);
			S_l=Comput_Saturation(psi,Solid);
			if (rank==0)
			{
			ofstream fin(FileName,ios::out);
			fin<<"The"<<n-freRe<<"th computation result:"<<endl;
			Re_l=u_ave*(NY+1)/niu_l;Re_g=u_ave*(NY+1)/niu_g;
		        fin<<"The Maximum velocity is: "<<setprecision(6)<<u_max<<"   Re_l="<<Re_l<<"   Re_g="<<Re_g<<endl;
			fin<<"Courant Number="<<u_max*dt/dx<<"	 Capillary Num="<<Capillary<<endl;
			fin<<"Peclet Number="<<u_max*dx/niu_s2<<endl;
			fin<<"The max relative error of velocity is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			fin<<"The relative permeability of component 1 is "<<Per_l[0]*reso*reso*1000/Permeability<<", "<<Per_l[1]*reso*reso*1000/Permeability<<", "<<Per_l[2]*reso*reso*1000/Permeability<<endl;
			fin<<"The relative permeability of component 2 is "<<Per_g[0]*reso*reso*1000/Permeability<<", "<<Per_g[1]*reso*reso*1000/Permeability<<", "<<Per_g[2]*reso*reso*1000/Permeability<<endl;
			fin<<"Satuation of Component 1: "<<S_l<<", "<<"The satuation of Component 2: "<<1-S_l<<endl;
			fin<<"The relative error of permiability computing is: "<<error_Per<<endl;
			fin<<"Elapsed time is "<< the<<"h"<<tme<<"m"<<tse<<"s"<<endl;
			fin<<"The expected completion time is "<<th<<"h"<<tm<<"m"<<ts<<"s"<<endl;
			fin<<endl;
			
			
			finish = MPI_Wtime();

			
			fin<<"The"<<n<<"th computation result:"<<endl;

			Re_l=u_ave*(NY+1)/niu_l;Re_g=u_ave*(NY+1)/niu_g;
		//=============================================================================================
		        remain=(n_max-n)*((finish-start)/n);
			th=int(remain/3600);
			tm=int((remain-th*3600)/60);
			ts=int(remain-(th*3600+tm*60));

			elaps=finish-start;
			the=int(elaps/3600);
			tme=int((elaps-the*3600)/60);
			tse=int(elaps-(the*3600+tme*60));

		//==============================================================================================
		        fin<<"The Maximum velocity is: "<<setprecision(6)<<u_max<<"   Re_l="<<Re_l<<"   Re_g="<<Re_g<<endl;
			fin<<"Courant Number="<<u_max*dt/dx<<"	 Capillary Num="<<Capillary<<endl;
			fin<<"Peclet Number="<<u_max*dx/niu_s2<<endl;
		//===============================================================================================
			fin<<"The max relative error of velocity is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			fin<<"The relative permeability of component 1 is "<<Per_l[0]*reso*reso*1000/Permeability<<", "<<Per_l[1]*reso*reso*1000/Permeability<<", "<<Per_l[2]*reso*reso*1000/Permeability<<endl;
			fin<<"The relative permeability of component 2 is "<<Per_g[0]*reso*reso*1000/Permeability<<", "<<Per_g[1]*reso*reso*1000/Permeability<<", "<<Per_g[2]*reso*reso*1000/Permeability<<endl;
			fin<<"Satuation of Component 1: "<<S_l<<", "<<"The satuation of Component 2: "<<1-S_l<<endl;
			fin<<"The relative error of permiability computing is: "<<error_Per<<endl;
			fin<<"Elapsed time is "<< the<<"h"<<tme<<"m"<<tse<<"s"<<endl;
			fin<<"The expected completion time is "<<th<<"h"<<tm<<"m"<<ts<<"s"<<endl;
			fin<<endl;
			fin.close();
		//==============================================================================================
			
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

			
			cout<<"The"<<n<<"th computation result:"<<endl;
			cout<<"The Density of point(NX/2,NY/2,NZ/2) is: "<<setprecision(6)
				<<rho[Solid[((NX+1)/para_size/2)][NY/2][NZ/2]]<<endl;
			
			cout<<"The Maximum velocity is: "<<setprecision(6)<<u_max<<"   Re_l="<<Re_l<<"   Re_g="<<Re_g<<endl;
			cout<<"Courant Number="<<u_max*dt/dx<<"	 Capillary Num="<<Capillary<<endl;
			cout<<"Peclet Number="<<u_max*dx/niu_s2<<endl;
			cout<<"The max relative error of uv is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			cout<<"The relative permeability of component 1 is "<<Per_l[0]*reso*reso*1000/Permeability<<", "<<Per_l[1]*reso*reso*1000/Permeability<<", "<<Per_l[2]*reso*reso*1000/Permeability<<endl;
			cout<<"The relative permeability of component 2 is "<<Per_g[0]*reso*reso*1000/Permeability<<", "<<Per_g[1]*reso*reso*1000/Permeability<<", "<<Per_g[2]*reso*reso*1000/Permeability<<endl;
			cout<<"Satuation of Component 1: "<<S_l<<", "<<"The satuation of Component 2: "<<1-S_l<<endl;
			cout<<"The relative error of permiability computing is: "<<error_Per<<endl;
			cout<<"Elapsed time is "<< the<<"h"<<tme<<"m"<<tse<<"s"<<endl;
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
			
			//===================================	
			if ((frePsi>=0) and (n%frePsi==0))
				if (Out_Mode==1)
					output_psi(n,psi,mirX,mirY,mirZ,mir,Solid,"psi");
				else
					output_psi_b(n,psi,mirX,mirY,mirZ,mir,Solid,"psi");


			if ((freHS>=0) and (n%freHS==0))
				if (Out_Mode==1)
					output_psi(n,rhoh,mirX,mirY,mirZ,mir,Solid,"HS");
				else
					output_psi_b(n,rhoh,mirX,mirY,mirZ,mir,Solid,"HS");
			//===================================
			
			if ((fre_backup>0) and (n%fre_backup==0)  and (n>0))
			        Backup(n,rho,psi,u,f);
			
			
			if(error!=error) {cout<<"PROGRAM STOP"<<endl;break;};
			if(U_max_ref>=5) {cout<<"PROGRAM STOP DUE TO HIGH VELOCITY"<<endl;break;}
		}	
	}

	
	
	if (fre_backup>=0)
			        Backup(n_max,rho,psi,u,f);

	MPI_Barrier(MPI_COMM_WORLD);

	for (int i=0;i<=Count;i++)
		{
		delete [] u[i];
		delete [] u0[i];
		delete [] f[i];
		delete [] F[i];
		delete [] fg[i];
		delete [] Fg[i];
		}
	
	delete [] f;
	delete [] psi;
	delete [] u;
	delete [] F;
	delete [] fg;
	delete [] Fg;
	delete [] rhoh;
	delete [] u0;
	delete [] rho;
	delete [] rho_r;
	delete [] rho_b;
	delete [] rhor;
	delete [] rhob;
	//delete [] forcex;
	//delete [] forcey;
	//delete [] forcez;
	delete [] SupInv;

	delete [] Sl;
	delete [] Sr;

	delete [] Permia;
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

void init_Sparse(int*** Solids, int*** Solid, double***Psis, double*** Psi_local, double*** Psis2, double*** Psi_local2, int* Sl,int* Sr)
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
			Psi_local2[i][j][k]=Psis2[int((s_c+i-(s_c+i)%Zoom)/Zoom)][int((j-j%Zoom)/Zoom)][int((k-k%Zoom)/Zoom)];
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






void Read_Rock(int*** Solids,double*** Psis,double*** Psis2, double* porosity,char poreFileName[128], char psiFileName[128], char psiFileName2[128])
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
double* Psi_Int2;

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
Psi_Int2 = new double[nx*ny*nz];

	

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

if (rank==0)	
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

if (rank==0)	
	cout<<"Concentration FILE READING COMPLETE. "<<endl;

	
	for (i=0;i<nx;i++)
		for (j=0;j<ny;j++)
			for (k=0;k<nz;k++)
			{
			   
			Psis[i][j][k]=Psi_Int[i*(ny)*(nz)+j*(nz)+k];
			}	
		
		
if (rank==0)
        
{
FILE *ftest3;
	ifstream fin3;
	
	ftest3 = fopen(psiFileName2, "r");

	if(ftest3 == NULL)
	{
		cout << "\n The concentration file (" << psiFileName <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest3);

	fin3.open(psiFileName2);


	
	// Reading pore geometry
	for(k=0 ; k<nz_a ; k++)
	for(j=0 ; j<ny_a ; j++)
	for(i=0 ; i<nx_a ; i++)
	
	{
		while(true)
		{	
			fin3 >> pore;
			if( (pore>-10000.0) and (pore<10000)) break;
		}
		if ((i<nx0) && (j<ny0) && (k<nz0))	Psi_Int2[i*ny*nz+j*nz+k] = pore;
		
	}
	fin3.close();

	// Mirroring the concentration
	if(mirX==1){
		for(i=nx0 ; i<nx ; i++)
		for(j=0   ; j<ny ; j++)
		for(k=0   ; k<nz ; k++)
				Psi_Int2[i*ny*nz+j*nz+k] = Psi_Int2[(nx-i-1)*ny*nz+j*nz+k];
	}

	if(mirY==1){
		for(j=ny0 ; j<ny ; j++)
		for(i=0   ; i<nx ; i++)
		for(k=0   ; k<nz ; k++)
				Psi_Int2[i*ny*nz+j*nz+k] = Psi_Int2[i*ny*nz+(ny-j-1)*nz+k];
	}

	if(mirZ==1){
		for(k=nz0 ; k<nz ; k++)
		for(j=0   ; j<ny ; j++)
		for(i=0   ; i<nx ; i++)
				Psi_Int2[i*ny*nz+j*nz+k] = Psi_Int2[i*ny*nz+j*nz+nz-k-1];
	}


}


	
	MPI_Bcast(Psi_Int2,nx*ny*nz,MPI_DOUBLE,0,MPI_COMM_WORLD);

if (rank==0)
	{	
	cout<<"Concentration FILE READING COMPLETE. "<<endl;
	cout<<endl;
	}

	
	for (i=0;i<nx;i++)
		for (j=0;j<ny;j++)
			for (k=0;k<nz;k++)
			{
			   
			Psis2[i][j][k]=Psi_Int2[i*(ny)*(nz)+j*(nz)+k];
			}	
		


MPI_Barrier(MPI_COMM_WORLD);

	delete [] Solid_Int;
	delete [] Psi_Int;
	delete [] Psi_Int2;

}



void init(double* rho, double** u, double** f,double* psi,double* rho_r, double* rho_b, double* rhor, double* rhob,double*** Psi_local, double*** Psi_local2, int* SupInv)
{	
      


	
	double usqr,vsqr;
	double c2,c4;
	
	rho0=1.0;dt=1.0/Zoom;dx=1.0/Zoom;
 
	if (lattice_v==1)
		{dx=dx_input;dt=dt_input;}

	lat_c=dx/dt;
	c_s=lat_c/sqrt(3);
	c_s2=lat_c*lat_c/3;

	c2=lat_c*lat_c;c4=c2*c2;
	
	
	niu=in_vis;
	
	
	tau_f=niu/(c_s2*dt)+0.5;
	s_v=1/tau_f;

        //tau_s1=niu_s1/1(c_s2*dt)+0.5;
	//tau_s2=niu_s2/1(c_s2*dt)+0.5;

	double eu; 
       
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
	{
	
	M_c[0]=1.0;
	M_c[1]=lat_c*lat_c;
	M_c[2]=lat_c*lat_c*lat_c*lat_c;
	M_c[3]=lat_c;
	M_c[4]=lat_c*lat_c*lat_c;
	M_c[5]=lat_c;
	M_c[6]=lat_c*lat_c*lat_c;
	M_c[7]=lat_c;	
	M_c[8]=lat_c*lat_c*lat_c;
	M_c[9]=lat_c*lat_c;
	M_c[10]=lat_c*lat_c*lat_c*lat_c;
	M_c[11]=lat_c*lat_c;
	M_c[12]=lat_c*lat_c*lat_c*lat_c;
	M_c[13]=lat_c*lat_c;
	M_c[14]=lat_c*lat_c;
	M_c[15]=lat_c*lat_c;
	M_c[16]=lat_c*lat_c*lat_c;
	M_c[17]=lat_c*lat_c*lat_c;
	M_c[18]=lat_c*lat_c*lat_c;



	for (int i=0;i<19;i++)
		for (int j=0;j<3;j++)
		elat[i][j]=e[i][j]*lat_c;

	for (int i=0;i<19;i++)
		for (int j=0;j<19;j++)
		M[i][j]*=M_c[i];

	Comput_MI(M,MI);

	}



	psi_solid=ContactAngle_parameter;

	

	for (int i=1;i<=Count;i++)	
			
		{
			u[i][0]=inivx;
			u[i][1]=inivy;
			u[i][2]=inivz;
			u_tmp[0]=u[i][0];
			u_tmp[1]=u[i][1];
			u_tmp[2]=u[i][2];
			psi[i]=Psi_local[(int)(SupInv[i]/((NY+1)*(NZ+1)))][(int)((SupInv[i]%((NY+1)*(NZ+1)))/(NZ+1))][SupInv[i]%(NZ+1)];
			rhoh[i]=Psi_local2[(int)(SupInv[i]/((NY+1)*(NZ+1)))][(int)((SupInv[i]%((NY+1)*(NZ+1)))/(NZ+1))][SupInv[i]%(NZ+1)];
			rho[i]=1.0;
			rho_r[i]=(psi[i]*rho[i]+rho[i])/2;
			rho_b[i]=rho[i]-rho_r[i];
			rhor[i]=0;
			rhob[i]=0;
			

			//forcex[i]=gx;
			//forcey[i]=gy;
			//forcez[i]=gz;
			if (stab==1)
				{gxs=0;gys=0;gzs=0;}
			else
				{gxs=gx;gys=gy;gzs=gz;}

			

			//INITIALIZATION OF m and f

			for (int lm=0;lm<19;lm++)
				{
					f[i][lm]=feq(lm,rho[i],u_tmp);
					eu=elat[lm][0]*u[i][0]+elat[lm][1]*u[i][1]+elat[lm][2]*u[i][2];
					fg[i][lm]=w[lm]*rhoh[i]*(1+3*eu/c2);
				}

		
		

	}

	

	 	
}



double feq_psi(int k,double rho, double u[3])
{

	double ux,uy,uz;
	double eu,uv,feq;
        double c2,c4,ls;
	double rho_0=1.0;
	
	c2=lat_c*lat_c;c4=c2*c2;
	eu=(elat[k][0]*u[0]+elat[k][1]*u[1]+elat[k][2]*u[2]);
	uv=(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
	feq=w[k]*rho*(1.0+3.0*eu/c2+4.5*eu*eu/c4-1.5*uv/c2);
	return feq;
}



double feq(int k,double rho, double u[3])
{

	double ux,uy,uz;
	double eu,uv,feq;
        double c2,c4,ls;
	
	double rho_0=1.0;
	
	c2=lat_c*lat_c;c4=c2*c2;
	eu=(elat[k][0]*u[0]+elat[k][1]*u[1]+elat[k][2]*u[2]);
	uv=(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);// WITH FORCE TERM:GRAVITY IN X DIRECTION
	feq=w[k]*rho*(1.0+3.0*eu/c2+4.5*eu*eu/c4-1.5*uv/c2);

			ux=u[0];
			uy=u[1];
			uz=u[2];

	for(int s=0;s<19;s++)
		meq[s]=0;
			meq[0]=rho;meq[3]=rho_0*ux;meq[5]=rho_0*uy;meq[7]=rho_0*uz;
			meq[1]=rho_0*(ux*ux+uy*uy+uz*uz);
			meq[9]=rho_0*(2*ux*ux-uy*uy-uz*uz);
			meq[11]=rho_0*(uy*uy-uz*uz);
			meq[13]=rho_0*ux*uy;meq[14]=rho_0*uy*uz;
			meq[15]=rho_0*ux*uz;

	feq=0;
	for (int j=0;j<19;j++)
		feq+=MI[k][j]*meq[j];		





	return feq;

}






void collision(double* rho,double** u,double** f,double** F,double* psi, double* rho_r, double* rho_b, double* rhor, double* rhob, int* SupInv,int*** Solid,int* Sl, int* Sr)
{

	MPI_Status status[4] ;
	MPI_Request request[4];
	MPI_Status status2[16] ;
	MPI_Request request2[16];
	int mpi_test;
	int mpi_test2;
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

double g_h[19];
double C[3];
double g_r[19],g_b[19];
double rho_0=1.0;
double lm0,lm1,cc,sum,uu;
double ux,uy,uz,nx,ny,nz;
double usqr,vsqr,eu,ef,cospsi,s_other;
double F_hat[19],GuoF[19],f_eq[19],u_tmp[3];
double m_l[19];
int i,j,m,ind_S;
int interi,interj,interk,ip,jp,kp;
double c2,c4;
double gxf,gyf,gzf;

	c2=lat_c*lat_c;c4=c2*c2;


	int* Gcl = new int[mpi_size];
	int* Gcr = new int[mpi_size];

	
	MPI_Gather(&cl,1,MPI_INT,Gcl,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(Gcl,mpi_size,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Gather(&cr,1,MPI_INT,Gcr,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(Gcr,mpi_size,MPI_INT,0,MPI_COMM_WORLD);

	double* recvl_psi;
	double* recvr_psi;
	double* sendl_rhor;
	double* sendr_rhor;
	double* sendl_rhob;
	double* sendr_rhob;
	
	double* sendl;
	double* sendr;

	double* sendl_s;
	double* sendr_s;

	double* recvl= new double[Gcl[rank]*5];
	double* recvr = new double[Gcr[rank]*5];



	double* sendl_psi = new double[Gcl[rank]];
	double* sendr_psi = new double[Gcr[rank]];
	double* recvl_rhor= new double[Gcl[rank]];
	double* recvr_rhor = new double[Gcr[rank]];
	double* recvl_rhob= new double[Gcl[rank]];
	double* recvr_rhob = new double[Gcr[rank]];

	double* recvl_s= new double[Gcl[rank]*5];
	double* recvr_s = new double[Gcr[rank]*5];

if (rank==0)
		{
		recvl_psi = new double[Gcr[mpi_size-1]];
		recvr_psi = new double[Gcl[rank+1]];
		sendl_rhor = new double[Gcr[mpi_size-1]];
		sendr_rhor = new double[Gcl[rank+1]];
		sendl_rhob = new double[Gcr[mpi_size-1]];
		sendr_rhob = new double[Gcl[rank+1]];
		sendl = new double[Gcr[mpi_size-1]*5];
		sendr = new double[Gcl[rank+1]*5];

		sendl_s = new double[Gcr[mpi_size-1]*5];
		sendr_s = new double[Gcl[rank+1]*5];


		for (int ka=0;ka<Gcr[mpi_size-1];ka++)
		        {
		                sendl_rhor[ka]=0;sendl_rhob[ka]=0;
				for (int kb=0;kb<5;kb++)
				{sendl[ka*5+kb]=0;sendl_s[ka*5+kb]=-1;}
		                
		        }
		 for (int ka=0;ka<Gcl[rank+1];ka++)
		        {
		                sendr_rhor[ka]=0;sendr_rhob[ka]=0;
		                for (int kb=0;kb<5;kb++)
				{sendr[ka*5+kb]=0;sendr_s[ka*5+kb]=-1;}
		        }       
		
		
		}	
		else
		if (rank==mpi_size-1)
			{
			recvl_psi = new double[Gcr[rank-1]];
			recvr_psi = new double[Gcl[0]];
			sendl_rhor = new double[Gcr[rank-1]];
			sendr_rhor = new double[Gcl[0]];
			sendl_rhob = new double[Gcr[rank-1]];
			sendr_rhob = new double[Gcl[0]];
			sendl = new double[Gcr[rank-1]*5];
			sendr = new double[Gcl[0]*5];

			sendl_s = new double[Gcr[rank-1]*5];
			sendr_s = new double[Gcl[0]*5];

			for (int ka=0;ka<Gcr[rank-1];ka++)
		        {
		                sendl_rhor[ka]=0;sendl_rhob[ka]=0;
				for (int kb=0;kb<5;kb++)
				{sendl[ka*5+kb]=0;sendl_s[ka*5+kb]=-1;}
		                
		        }
		        for (int ka=0;ka<Gcl[0];ka++)
		        {
		                sendr_rhor[ka]=0;sendr_rhob[ka]=0;
				for (int kb=0;kb<5;kb++)
				{sendr[ka*5+kb]=0;sendr_s[ka*5+kb]=-1;}
		                
		        }       
			
			
			}
			else
			{
			recvl_psi = new double[Gcr[rank-1]];
			recvr_psi = new double[Gcl[rank+1]];
			sendl_rhor = new double[Gcr[rank-1]];
			sendr_rhor = new double[Gcl[rank+1]];
			sendl_rhob = new double[Gcr[rank-1]];
			sendr_rhob = new double[Gcl[rank+1]];
			sendl = new double[Gcr[rank-1]*5];
			sendr = new double[Gcl[rank+1]*5];

			sendl_s = new double[Gcr[rank-1]*5];
			sendr_s = new double[Gcl[rank+1]*5];

			

			for (int ka=0;ka<Gcr[rank-1];ka++)
		        {
		                sendl_rhor[ka]=0;sendl_rhob[ka]=0;
				for (int kb=0;kb<5;kb++)
				{sendl[ka*5+kb]=0;sendl_s[ka*5+kb]=-1;}
		                
		        }
		        for (int ka=0;ka<Gcl[rank+1];ka++)
		        {
		                sendr_rhor[ka]=0;sendr_rhob[ka]=0;
				for (int kb=0;kb<5;kb++)
				{sendr[ka*5+kb]=0;sendr_s[ka*5+kb]=-1;}
			}
			
			}
			
			

			for(i=1;i<=Gcl[rank];i++)
			        {
			        sendl_psi[i-1]=psi[i];
			        }
			for(j=Count-Gcr[rank]+1;j<=Count;j++)
			        {
				sendr_psi[j-(Count-Gcr[rank]+1)]=psi[j];
				}
			
		
MPI_Barrier(MPI_COMM_WORLD);

//cout<<"@@@@@@@@@@@   "<<n<<endl;
if (rank==0)
		{
		
		
		MPI_Isend(sendr_psi, Gcr[0], MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(sendl_psi, Gcl[0], MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(recvr_psi, Gcl[1], MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(recvl_psi, Gcr[mpi_size-1], MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );
		
		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(sendl_psi, Gcl[rank], MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr_psi, Gcr[rank], MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl_psi, Gcr[rank-1], MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(recvr_psi, Gcl[0], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request[3]);
			
			}
			else
			{

      			MPI_Isend(sendl_psi, Gcl[rank], MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(sendr_psi, Gcr[rank], MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(recvl_psi, Gcr[rank-1], MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(recvr_psi, Gcl[rank+1], MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);
			
			
			};

	
	MPI_Waitall(4,request, status);

	MPI_Testall(4,request,&mpi_test,status);	



	for(int ci=1;ci<=Count;ci++)	
	

		{	
				i=(int)(SupInv[ci]/((NY+1)*(NZ+1)));
				j=(int)((SupInv[ci]%((NY+1)*(NZ+1)))/(NZ+1));
				m=(int)(SupInv[ci]%(NZ+1));   

		//cout<<i<<" "<<j<<" "<<m<<" /before "<<rhor[ci]<<" nth= "<<n<<"  the ci= "<<ci<<"  rank=  "<<rank<<endl;


			C[0]=0;C[1]=0;C[2]=0;ind_S=0;
	for (int tmpi=0;tmpi<19;tmpi++)
		{
		        //cout<<f[ci][tmpi]<<endl;
			//-------------------PERIODIC BOUNDARY CONDITION---------------------------
		if (in_BC>0)
		{	     
			interi=i+e[tmpi][0];
			if (((pre_xn-1)* (vel_xn-1)==0) and (rank==0) and (interi<0))
			        interi=0;
			if (((pre_xp-1)*(vel_xp-1)==0) and (rank==mpi_size-1) and (interi>=nx_l))
			        interi=nx_l-1;
			
			interj=j+e[tmpi][1];
			if ((pre_yn-1)*(vel_yn-1)==0)
			        {if (interj<0) {interj=0;}}
			else
			        {if (interj<0) {interj=NY;}}
			 
			if ((pre_yp-1)*(vel_yp-1)==0)
			        {if (interj>NY) {interj=NY;}}
			else
			        {if (interj>NY) {interj=0;}}
			
			
			interk=m+e[tmpi][2];
			if ((pre_zn-1)*(vel_zn-1)==0)
			        {if (interk<0) {interk=0;}}
			else
			        {if (interk<0) {interk=NZ;}}
			
			if ((pre_zp-1)*(vel_zp-1)==0)
			        {if (interk>NZ) {interk=NZ;}}
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
			
			if (interi<0)
			{
			        if (Sl[interj*(NZ+1)+interk]>0)
					{
					C[0]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][0]*recvl_psi[Sl[interj*(NZ+1)+interk]-1];
					C[1]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][1]*recvl_psi[Sl[interj*(NZ+1)+interk]-1];
					C[2]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][2]*recvl_psi[Sl[interj*(NZ+1)+interk]-1];

					}
				else
			        	{
			        	ind_S=1;
					C[0]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][0]*psi_solid;
					C[1]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][1]*psi_solid;
					C[2]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][2]*psi_solid;
			              
			        	}
			}
			
			
			if (interi>=nx_l)
			{
			        if (Sr[interj*(NZ+1)+interk]>0)
			        	{
					C[0]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][0]*recvr_psi[Sr[interj*(NZ+1)+interk]-1];
					C[1]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][1]*recvr_psi[Sr[interj*(NZ+1)+interk]-1];
					C[2]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][2]*recvr_psi[Sr[interj*(NZ+1)+interk]-1];

			        	}
			        else
			                
			                {
			                ind_S=1;
					C[0]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][0]*psi_solid;
					C[1]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][1]*psi_solid;
					C[2]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][2]*psi_solid;
      			                }
      			}
      			
      			
      			
			if ((interi>=0) and (interi<nx_l))
			{
			        if (Solid[interi][interj][interk]>0)
			        	{
					C[0]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][0]*psi[Solid[interi][interj][interk]];
					C[1]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][1]*psi[Solid[interi][interj][interk]];
					C[2]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][2]*psi[Solid[interi][interj][interk]];
			        	}        
			        else
			        	{
                                        ind_S=1;
			               	C[0]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][0]*psi_solid;
					C[1]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][1]*psi_solid;
					C[2]+=3.0/(lat_c*lat_c*dt)*w[tmpi]*elat[tmpi][2]*psi_solid;
			        	}
			
			
			}
		}

		uu=u[ci][0]*u[ci][0]+u[ci][1]*u[ci][1]+u[ci][2]*u[ci][2];

		if ((sqrt((rho_r[ci]-rho_b[ci])*(rho_r[ci]-rho_b[ci]))>=0.9) and (ind_S=1))
		{C[0]=0;C[1]=0;C[2]=0;}

		//C[0]=0;C[1]=0;C[2]=0;
			cc=sqrt(C[0]*C[0]+C[1]*C[1]+C[2]*C[2]);
			if (cc>0)
			        {nx=C[0]/cc;ny=C[1]/cc;nz=C[2]/cc;}
			else
			        {nx=0;ny=0;nz=0;}

			gxf=gxs;
			gyf=gys+Gravity*Buoyancy_parameter*(rhoh[ci]-ref_psi);
			gzf=gzs;
			
			//=================FORCE TERM_GUO=========================================
			for (int k=0;k<19;k++)
			{	
			lm0=((elat[k][0]-u[ci][0])*gxf+(elat[k][1]-u[ci][1])*gyf+(elat[k][2]-u[ci][2])*gzf)/c_s2;
			lm1=(elat[k][0]*u[ci][0]+elat[k][1]*u[ci][1]+elat[k][2]*u[ci][2])*(elat[k][0]*gxf+elat[k][1]*gyf+elat[k][2]*gzf)/(c_s2*c_s2);
			GuoF[k]=w[k]*(lm0+lm1);
			//GuoF[k]=0.0;
			}


			
			//=====================equilibrium of moment=================================
			ux=u[ci][0];
			uy=u[ci][1];
			uz=u[ci][2];
			
			for(int k=0;k<19;k++)
				meq[k]=0;

			
			
			meq[0]=rho[ci];meq[3]=rho_0*ux;meq[5]=rho_0*uy;meq[7]=rho_0*uz;
			meq[1]=rho_0*(ux*ux+uy*uy+uz*uz)+CapA*cc;
			meq[9]=rho_0*(2*ux*ux-uy*uy-uz*uz)+0.5*CapA*cc*(2*nx*nx-ny*ny-nz*nz);
			meq[11]=rho_0*(uy*uy-uz*uz)+0.5*CapA*cc*(ny*ny-nz*nz);
			meq[13]=rho_0*ux*uy+0.5*CapA*cc*(nx*ny);
			meq[14]=rho_0*uy*uz+0.5*CapA*cc*(ny*nz);
			meq[15]=rho_0*ux*uz+0.5*CapA*cc*(nx*nz);
			
			
			s_v=niu_g+(psi[ci]+1.0)/2.0*(niu_l-niu_g);
			s_v=1.0/(s_v/(c_s2*dt)+0.5);
			s_other=8*(2-s_v)/(8-s_v);

			tau_s=niu_s2+(psi[ci]+1.0)/2.0*(niu_s1-niu_s2);
			tau_s=tau_s/(c_s2*dt)+0.5;
			
	//cout<<"@@@@@@@@@   "<<s_v<<"  "<<C[0]<<"   "<<C[1]<<"  "<<C[2]<<endl;
	//cout<<"@@@@@@@@@   "<<s_v<<"  "<<ux<<"   "<<uy<<"  "<<uz<<"  "<<rho_r[ci]<<" "<<rho_b[ci]<<endl;
	
	
	S[1]=s_v;S[2]=s_v;S[4]=s_other;S[6]=s_other;S[8]=s_other;S[9]=s_v;
	S[10]=s_v;S[11]=s_v;S[12]=s_v;S[13]=s_v;S[14]=s_v;S[15]=s_v;S[16]=s_other;
	S[17]=s_other;S[18]=s_other;



			//============================================================================
			uu=u[ci][0]*u[ci][0]+u[ci][1]*u[ci][1]+u[ci][2]*u[ci][2];
			//============================================================================



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
					m_l[mi]=m_l[mi]-S[mi]*(m_l[mi]-meq[mi])+dt*F_hat[mi];
					eu=elat[mi][0]*u[ci][0]+elat[mi][1]*u[ci][1]+elat[mi][2]*u[ci][2];

			//=========================EVOLUTION OF G FOR SIMPLE BOUNDARY========================
               		//	g_h[mi]=fg[ci][mi]-(fg[ci][mi]-w[mi]*rhoh[ci]*(1+3*eu))/tau_s;
			//=========================EVOLUTION OF G FOR COMPLEX BOUNDARY=======================
			g_h[mi]=fg[ci][mi]-(fg[ci][mi]-w[mi]*rhoh[ci]*(1+3*eu/c2+4.5*eu*eu/c4-1.5*uu/c2))/tau_s;
			//============================================================================
					}
			//============================================================================


		for (int mi=0; mi<19; mi++)
			{
			sum=0;
			for (int mj=0; mj<19; mj++)
				sum+=MI[mi][mj]*m_l[mj];

			//F[ci][mi]=0;
			ip=i+e[mi][0];
			jp=j+e[mi][1];if (jp<0) {jp=NY;}; if (jp>NY) {jp=0;};
			kp=m+e[mi][2];if (kp<0) {kp=NZ;}; if (kp>NZ) {kp=0;};


			if (ip<0) 
				if (Sl[jp*(NZ+1)+kp]>0)
				{
				sendl[(Sl[jp*(NZ+1)+kp]-1)*5+FLN[mi]]=sum;
				sendl_s[(Sl[jp*(NZ+1)+kp]-1)*5+FLN[mi]]=g_h[mi];
				//sendl_rhob[Sl[jp*(NZ+1)+kp]-1]+=g_b[lm];
				//cout<<g_r[lm]<<"    1"<<endl;
				}
				else
				{
				F[ci][LR[mi]]=sum;Fg[ci][LR[mi]]=g_h[mi];
				}
					
						
					
					
			if (ip>=nx_l)
				if (Sr[jp*(NZ+1)+kp]>0)
				{
				sendr[(Sr[jp*(NZ+1)+kp]-1)*5+FRP[mi]]=sum;
				sendr_s[(Sr[jp*(NZ+1)+kp]-1)*5+FRP[mi]]=g_h[mi];
				//sendr_rhob[Sr[jp*(NZ+1)+kp]-1]+=g_b[lm];
				//cout<<g_r[lm]<<"    2"<<endl;
				}
				else
				{
				F[ci][LR[mi]]=sum;
				Fg[ci][LR[mi]]=g_h[mi];
				}

			if ((ip>=0) and (ip<nx_l)) 
				if (Solid[ip][jp][kp]>0)
				{
				F[Solid[ip][jp][kp]][mi]=sum;
				Fg[Solid[ip][jp][kp]][mi]=g_h[mi];
				
				}
				else
				{
				F[ci][LR[mi]]=sum;Fg[ci][LR[mi]]=g_h[mi];
				}
		//=======================G streaming=================================================
		//for(int lm=0;lm<19;lm++)
                //{
                 eu=elat[mi][0]*u[ci][0]+elat[mi][1]*u[ci][1]+elat[mi][2]*u[ci][2];
                 //g_r[mi]=w[mi]*rho_r[ci]*(1+3*eu/c2);
                 //g_b[mi]=w[mi]*rho_b[ci]*(1+3*eu/c2);

		 g_r[mi]=w[mi]*rho_r[ci]*(1+3*eu/c2+4.5*eu*eu/c4-1.5*uu/c2);
                 g_b[mi]=w[mi]*rho_b[ci]*(1+3*eu/c2+4.5*eu*eu/c4-1.5*uu/c2);
		
			
		//	cout<<" "<<g_r[lm]<<" "<<g_b[lm]<<"  the number "<<n<<"  vector "<<lm<<endl;
                 }

		
                 
           if (cc>0)
           for(int kk=1;kk<19;kk+=2)
                {
                ef=elat[kk][0]*C[0]+elat[kk][1]*C[1]+elat[kk][2]*C[2];
                cospsi=g_r[kk]<g_r[kk+1]?g_r[kk]:g_r[kk+1];
                cospsi=cospsi<g_b[kk]?cospsi:g_b[kk];
                cospsi=cospsi<g_b[kk+1]?cospsi:g_b[kk+1];
                cospsi*=ef/cc;
		
		//cout<<"@@@@@     "<<ef/cc<<endl;

                g_r[kk]+=cospsi;
                g_r[kk+1]-=cospsi;
                g_b[kk]-=cospsi;
                g_b[kk+1]+=cospsi;
		
                }      
			

			   
		       for(int lm=0;lm<19;lm++)


			{       
				//cout<<f[ci][lm]<<endl;

				ip=i+e[lm][0];
				jp=j+e[lm][1];if (jp<0) {jp=NY;}; if (jp>NY) {jp=0;};
				kp=m+e[lm][2];if (kp<0) {kp=NZ;}; if (kp>NZ) {kp=0;};

			
				if (ip<0) 
					if (Sl[jp*(NZ+1)+kp]>0)
						{
						sendl_rhor[Sl[jp*(NZ+1)+kp]-1]+=g_r[lm];
						sendl_rhob[Sl[jp*(NZ+1)+kp]-1]+=g_b[lm];
						//cout<<g_r[lm]<<"    1"<<endl;
						}
					else
						{
						rhor[ci]+=g_r[lm];
						rhob[ci]+=g_b[lm];
						}
					
						
					
					
				if (ip>=nx_l)
					if (Sr[jp*(NZ+1)+kp]>0)
						{
						sendr_rhor[Sr[jp*(NZ+1)+kp]-1]+=g_r[lm];
						sendr_rhob[Sr[jp*(NZ+1)+kp]-1]+=g_b[lm];
						//cout<<g_r[lm]<<"    2"<<endl;
						}
					else
						{
						rhor[ci]+=g_r[lm];
						rhob[ci]+=g_b[lm];
						}

				if ((ip>=0) and (ip<nx_l)) 
					if (Solid[ip][jp][kp]>0)
						{
						rhor[Solid[ip][jp][kp]]+=g_r[lm];
						rhob[Solid[ip][jp][kp]]+=g_b[lm];
						
						}
					else
						{
						rhor[ci]+=g_r[lm];
						rhob[ci]+=g_b[lm];
						}
					
					
			}
		
		
			
			}
                        
	
	if (rank==0)
		{
		MPI_Isend(sendr_s, Gcl[1]*5, MPI_DOUBLE, rank+1, rank*2+1+30000, MPI_COMM_WORLD,&request2[12]);
      		MPI_Isend(sendl_s, Gcr[mpi_size-1]*5, MPI_DOUBLE, mpi_size-1, rank*2+30000, MPI_COMM_WORLD,&request2[13]);
		MPI_Irecv(recvr_s, Gcr[0]*5, MPI_DOUBLE, rank+1, (rank+1)*2+30000, MPI_COMM_WORLD,&request2[14]);		
      		MPI_Irecv(recvl_s, Gcl[0]*5, MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1+30000, MPI_COMM_WORLD,&request2[15] );
		
		MPI_Isend(sendr_rhor, Gcl[1], MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request2[0]);
      		MPI_Isend(sendl_rhor, Gcr[mpi_size-1], MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request2[1]);
		MPI_Irecv(recvr_rhor, Gcr[0], MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request2[2]);		
      		MPI_Irecv(recvl_rhor, Gcl[0], MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request2[3] );
		MPI_Isend(sendr_rhob, Gcl[1], MPI_DOUBLE, rank+1, rank*2+1+10000, MPI_COMM_WORLD,&request2[4]);
      		MPI_Isend(sendl_rhob, Gcr[mpi_size-1], MPI_DOUBLE, mpi_size-1, rank*2+10000, MPI_COMM_WORLD,&request2[5]);
		MPI_Irecv(recvr_rhob, Gcr[0], MPI_DOUBLE, rank+1, (rank+1)*2+10000, MPI_COMM_WORLD,&request2[6]);		
      		MPI_Irecv(recvl_rhob, Gcl[0], MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1+10000, MPI_COMM_WORLD,&request2[7] );
		MPI_Isend(sendr, Gcl[1]*5, MPI_DOUBLE, rank+1, rank*2+1+20000, MPI_COMM_WORLD,&request2[8]);
      		MPI_Isend(sendl, Gcr[mpi_size-1]*5, MPI_DOUBLE, mpi_size-1, rank*2+20000, MPI_COMM_WORLD,&request2[9]);
		MPI_Irecv(recvr, Gcr[0]*5, MPI_DOUBLE, rank+1, (rank+1)*2+20000, MPI_COMM_WORLD,&request2[10]);		
      		MPI_Irecv(recvl, Gcl[0]*5, MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1+20000, MPI_COMM_WORLD,&request2[11] );
		}
		else
		if (rank==mpi_size-1)
			{
			
			MPI_Isend(sendl_s, Gcr[rank-1]*5, MPI_DOUBLE, rank-1, rank*2+30000, MPI_COMM_WORLD,&request2[12]);
      			MPI_Isend(sendr_s, Gcl[0]*5,  MPI_DOUBLE, 0, rank*2+1+30000, MPI_COMM_WORLD,&request2[13]);
			MPI_Irecv(recvl_s, Gcl[rank]*5, MPI_DOUBLE, rank-1, (rank-1)*2+1+30000, MPI_COMM_WORLD,&request2[14] );
      			MPI_Irecv(recvr_s, Gcr[rank]*5, MPI_DOUBLE, 0, 0+30000, MPI_COMM_WORLD,&request2[15]);

			MPI_Isend(sendl_rhor, Gcr[rank-1], MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request2[0]);
      			MPI_Isend(sendr_rhor, Gcl[0],  MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request2[1]);
			MPI_Irecv(recvl_rhor, Gcl[rank], MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request2[2] );
      			MPI_Irecv(recvr_rhor, Gcr[rank], MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request2[3]);
			MPI_Isend(sendl_rhob, Gcr[rank-1], MPI_DOUBLE, rank-1, rank*2+10000, MPI_COMM_WORLD,&request2[4]);
      			MPI_Isend(sendr_rhob, Gcl[0],  MPI_DOUBLE, 0, rank*2+1+10000, MPI_COMM_WORLD,&request2[5]);
			MPI_Irecv(recvl_rhob, Gcl[rank], MPI_DOUBLE, rank-1, (rank-1)*2+1+10000, MPI_COMM_WORLD,&request2[6] );
      			MPI_Irecv(recvr_rhob, Gcr[rank], MPI_DOUBLE, 0, 0+10000, MPI_COMM_WORLD,&request2[7]);
			MPI_Isend(sendl, Gcr[rank-1]*5, MPI_DOUBLE, rank-1, rank*2+20000, MPI_COMM_WORLD,&request2[8]);
      			MPI_Isend(sendr, Gcl[0]*5,  MPI_DOUBLE, 0, rank*2+1+20000, MPI_COMM_WORLD,&request2[9]);
			MPI_Irecv(recvl, Gcl[rank]*5, MPI_DOUBLE, rank-1, (rank-1)*2+1+20000, MPI_COMM_WORLD,&request2[10] );
      			MPI_Irecv(recvr, Gcr[rank]*5, MPI_DOUBLE, 0, 0+20000, MPI_COMM_WORLD,&request2[11]);
			}
			else
			{
			MPI_Isend(sendl_s, Gcr[rank-1]*5, MPI_DOUBLE, rank-1, rank*2+30000, MPI_COMM_WORLD,&request2[12]);
      			MPI_Isend(sendr_s, Gcl[rank+1]*5, MPI_DOUBLE, rank+1, rank*2+1+30000, MPI_COMM_WORLD,&request2[13]);
			MPI_Irecv(recvl_s, Gcl[rank]*5, MPI_DOUBLE, rank-1, (rank-1)*2+1+30000, MPI_COMM_WORLD,&request2[14]);
      			MPI_Irecv(recvr_s, Gcr[rank]*5, MPI_DOUBLE, rank+1, (rank+1)*2+30000, MPI_COMM_WORLD,&request2[15]);


      			MPI_Isend(sendl_rhor, Gcr[rank-1], MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request2[0]);
      			MPI_Isend(sendr_rhor, Gcl[rank+1], MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request2[1]);
			MPI_Irecv(recvl_rhor, Gcl[rank], MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request2[2]);
      			MPI_Irecv(recvr_rhor, Gcr[rank], MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request2[3]);
			MPI_Isend(sendl_rhob, Gcr[rank-1], MPI_DOUBLE, rank-1, rank*2+10000, MPI_COMM_WORLD,&request2[4]);
      			MPI_Isend(sendr_rhob, Gcl[rank+1], MPI_DOUBLE, rank+1, rank*2+1+10000, MPI_COMM_WORLD,&request2[5]);
			MPI_Irecv(recvl_rhob, Gcl[rank], MPI_DOUBLE, rank-1, (rank-1)*2+1+10000, MPI_COMM_WORLD,&request2[6]);
      			MPI_Irecv(recvr_rhob, Gcr[rank], MPI_DOUBLE, rank+1, (rank+1)*2+10000, MPI_COMM_WORLD,&request2[7]);
			MPI_Isend(sendl, Gcr[rank-1]*5, MPI_DOUBLE, rank-1, rank*2+20000, MPI_COMM_WORLD,&request2[8]);
      			MPI_Isend(sendr, Gcl[rank+1]*5, MPI_DOUBLE, rank+1, rank*2+1+20000, MPI_COMM_WORLD,&request2[9]);
			MPI_Irecv(recvl, Gcl[rank]*5, MPI_DOUBLE, rank-1, (rank-1)*2+1+20000, MPI_COMM_WORLD,&request2[10]);
      			MPI_Irecv(recvr, Gcr[rank]*5, MPI_DOUBLE, rank+1, (rank+1)*2+20000, MPI_COMM_WORLD,&request2[11]);
			};

	
	MPI_Waitall(16,request2, status2);

	MPI_Testall(16,request2,&mpi_test2,status2);	
		
			for(i=1;i<=Gcl[rank];i++)
			        {
			
			        rhor[i]+=recvl_rhor[i-1];
			        rhob[i]+=recvl_rhob[i-1];
				
				for (int lm=0;lm<5;lm++)
					if (recvl[(i-1)*5+lm]>0)
				        {F[i][RP[lm]]=recvl[(i-1)*5+lm];
					Fg[i][RP[lm]]=recvl_s[(i-1)*5+lm];}
			        }
			for(j=Count-Gcr[rank]+1;j<=Count;j++)
			        {
				rhor[j]+=recvr_rhor[j-(Count-Gcr[rank]+1)];
				rhob[j]+=recvr_rhob[j-(Count-Gcr[rank]+1)];
				for (int lm=0;lm<5;lm++)
					if (recvr[(j-(Count-Gcr[rank]+1))*5+lm]>0)
			        	{F[j][LN[lm]]=recvr[(j-(Count-Gcr[rank]+1))*5+lm];
					Fg[j][LN[lm]]=recvr_s[(j-(Count-Gcr[rank]+1))*5+lm];}
				}
			
			
			
			
			
	delete [] recvl_psi;
	delete [] recvr_psi;
	delete [] sendr_psi;
	delete [] sendl_psi;
	delete [] Gcl;
	delete [] Gcr;
	delete [] sendl_rhor;
	delete [] sendr_rhor;
	delete [] recvl_rhor;
	delete [] recvr_rhor;		
	delete [] sendl_rhob;
	delete [] sendr_rhob;
	delete [] recvl_rhob;
	delete [] recvr_rhob;
	delete [] sendl;
	delete [] sendr;
	delete [] recvl;
	delete [] recvr;

	delete [] sendl_s;
	delete [] sendr_s;
	delete [] recvl_s;
	delete [] recvr_s;
	

}





void comput_macro_variables( double* rho,double** u,double** u0,double** f,double** F,double* rho_r, double* rho_b, double* rhor, double* rhob, double* psi,int* SupInv,int*** Solid,double*** Psi_local)
{
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	for(int i=1;i<=Count;i++)	
                   
			{
			
				u0[i][0]=u[i][0];
				u0[i][1]=u[i][1];
				u0[i][2]=u[i][2];
				rho[i]=0;rhoh[i]=0;
				u[i][0]=0;
				u[i][1]=0;
				u[i][2]=0;
	
				for(int k=0;k<19;k++)
					{
					fg[i][k]=Fg[i][k];
					f[i][k]=F[i][k];
					rho[i]+=f[i][k];
					rhoh[i]+=fg[i][k];
					u[i][0]+=elat[k][0]*f[i][k];
					u[i][1]+=elat[k][1]*f[i][k];
					u[i][2]+=elat[k][2]*f[i][k];
					}
				
				rho_r[i]=rhor[i];
				rho_b[i]=rhob[i];
				rhor[i]=0;
				rhob[i]=0;
				u[i][0]=(u[i][0]+dt*gxs)/rho[i];
				u[i][1]=(u[i][1]+dt*gys)/rho[i];
				u[i][2]=(u[i][2]+dt*gzs)/rho[i];
				
				
				psi[i]=(rho_r[i]-rho_b[i])/(rho_r[i]+rho_b[i]);
				
			}
			
		
	if (in_BC==1)
                {
                   if (((pre_xn==1) or (vel_xn==1)) and (rank==0))    
                           {
                           for(int j=0;j<=NY;j++)
                                   for (int k=0;k<=NZ;k++)
                                   {
                                   //rho_r[Solid[0][j][k]]=(Psi_local[0][j][k]*1.0+1.0)/2;
                                   //rho_b[Solid[0][j][k]]=1.0-rho_r[Solid[0][j][k]];
                                   //psi[Solid[0][j][k]]=Psi_local[0][j][k];
                                   psi[Solid[0][j][k]]=psi[Solid[1][j][k]];
                                   rho_r[Solid[0][j][k]]=(psi[Solid[0][j][k]]*1.0+1.0)/2;
                                   rho_b[Solid[0][j][k]]=1.0-rho_r[Solid[0][j][k]];
                                   
                                   }
                           }
                        
                      if (((pre_xp==1) or (vel_xp==1)) and (rank==mpi_size-1))    
                           {
                           for(int j=0;j<=NY;j++)
                                   for (int k=0;k<=NZ;k++)
                                   {
                                   //rho_r[Solid[nx_l-1][j][k]]=(Psi_local[nx_l-1][j][k]*1.0+1.0)/2;
                                   //rho_b[Solid[nx_l-1][j][k]]=1.0-rho_r[Solid[nx_l-1][j][k]];
                                   //psi[Solid[nx_l-1][j][k]]=Psi_local[nx_l-1][j][k];
                                   psi[Solid[nx_l-1][j][k]]=psi[Solid[nx_l-2][j][k]];
                                   rho_r[Solid[nx_l-1][j][k]]=(psi[Solid[nx_l-1][j][k]]*1.0+1.0)/2;
                                   rho_b[Solid[nx_l-1][j][k]]=1.0-rho_r[Solid[nx_l-1][j][k]];
                                   }
                           }  
                        
                        
                         if ((pre_yn==1) or (vel_yn==1))    
                           {
                           for(int i=0;i<nx_l;i++)
                                   for (int k=0;k<=NZ;k++)
                                   {
                                   //rho_r[Solid[i][0][k]]=(Psi_local[i][0][k]*1.0+1.0)/2;
                                   //rho_b[Solid[i][0][k]]=1.0-rho_r[Solid[i][0][k]];
                                   //psi[Solid[i][0][k]]=Psi_local[i][0][k];
                                   psi[Solid[i][0][k]]=psi[Solid[i][1][k]];
                                   rho_r[Solid[i][0][k]]=(psi[Solid[i][0][k]]*1.0+1.0)/2;
                                   rho_b[Solid[i][0][k]]=1.0-rho_r[Solid[i][0][k]];
                                   }
                           }  
                        
                        if ((pre_yp==1) or (vel_yp==1))    
                           {
                           for(int i=0;i<nx_l;i++)
                                   for (int k=0;k<=NZ;k++)
                                   {
                                   //rho_r[Solid[i][NY][k]]=(Psi_local[i][NY][k]*1.0+1.0)/2;
                                   //rho_b[Solid[i][NY][k]]=1.0-rho_r[Solid[i][NY][k]];
                                   //psi[Solid[i][NY][k]]=Psi_local[i][NY][k];
                                   psi[Solid[i][NY][k]]=psi[Solid[i][NY-1][k]];
                                   rho_r[Solid[i][NY][k]]=(psi[Solid[i][NY][k]]*1.0+1.0)/2;
                                   rho_b[Solid[i][NY][k]]=1.0-rho_r[Solid[i][NY][k]];
                                   }
                           }  
                        
                        if ((pre_zn==1) or (vel_zn==1))    
                           {
                           for(int i=0;i<nx_l;i++)
                                   for (int j=0;j<=NY;j++)
                                   {
                                   //rho_r[Solid[i][j][0]]=(Psi_local[i][j][0]*1.0+1.0)/2;
                                   //rho_b[Solid[i][j][0]]=1.0-rho_r[Solid[i][j][0]];
                                   //psi[Solid[i][j][0]]=Psi_local[i][j][0];
                                   psi[Solid[i][j][0]]=psi[Solid[i][j][1]];
                                   rho_r[Solid[i][j][0]]=(psi[Solid[i][j][0]]*1.0+1.0)/2;
                                   rho_b[Solid[i][j][0]]=1.0-rho_r[Solid[i][j][0]];
                                   }
                           }  
                        
                    if ((pre_zp==1) or (vel_zp==1))    
                           {
                           for(int i=0;i<nx_l;i++)
                                   for (int j=0;j<=NY;j++)
                                   {
                                   //rho_r[Solid[i][j][NZ]]=(Psi_local[i][j][NZ]*1.0+1.0)/2;
                                   //rho_b[Solid[i][j][NZ]]=1.0-rho_r[Solid[i][j][NZ]];
                                   //psi[Solid[i][j][NZ]]=Psi_local[i][j][NZ];
                                   psi[Solid[i][j][NZ]]=psi[Solid[i][j][NZ-1]];
                                   rho_r[Solid[i][j][NZ]]=(psi[Solid[i][j][NZ]]*1.0+1.0)/2;
                                   rho_b[Solid[i][j][NZ]]=1.0-rho_r[Solid[i][j][NZ]];
                                   }
                           }      
                        
                        
                }
	
                     
			
	//MPI_Barrier(MPI_COMM_WORLD); 

}




void Solute_ZeroFlux_BC(int sxp,double df_sxp, int sxn,double df_sxn,int syp,double df_syp,int syn,double df_syn,int szp,double df_szp,int szn,double df_szn,double** Fg,int*** Solid,double* rho_r,double** u,double** fg)
{
int Q=19;
int rank = MPI :: COMM_WORLD . Get_rank ();
int mpi_size=MPI :: COMM_WORLD . Get_size ();
double rho_t,rho_p;
double srho_syp,srho_syn,srho_szp,srho_szn,srho_sxp,srho_sxn;
double u_ls[3];

if (Sub_BC_psi==0)
{
if ((syp-1)*(syn-1)==0)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)	
	  	{
		if ((syp==1) and (Solid[i][NY][k]>0))
			{
			rho_t=df_syp;
			for (int ks=0;ks<=4;ks++)
				rho_t+=Fg[Solid[i][NY][k]][SYN[ks]];

			rho_p=(rho_t)/(1.0/6.0);
			Fg[Solid[i][NY][k]][SYP[0]]=rho_p/18.0;			
			for (int ks=1;ks<=4;ks++)
				Fg[Solid[i][NY][k]][SYP[ks]]=rho_p/36.0;
			}
		        
		if ((syn==1) and (Solid[i][0][k]>0))
		        {
			rho_t=df_syn;
			for (int ks=0;ks<=4;ks++)
				rho_t+=Fg[Solid[i][0][k]][SYP[ks]];

			rho_p=(rho_t)/(1.0/6.0);
			Fg[Solid[i][0][k]][SYN[0]]=rho_p/18.0;			
			for (int ks=1;ks<=4;ks++)
				Fg[Solid[i][0][k]][SYN[ks]]=rho_p/36.0;
			}
		
		}
	        
if ((szp-1)*(szn-1)==0)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
		{
		if ((szp==1) and (Solid[i][j][NZ]>0))
			{
			rho_t=df_szp;
			for (int ks=0;ks<=4;ks++)
				rho_t+=Fg[Solid[i][j][NZ]][SZN[ks]];//######

			rho_p=(rho_t)/(1.0/6.0);//######
			Fg[Solid[i][j][NZ]][SZP[0]]=rho_p/18.0;//######			
			for (int ks=1;ks<=4;ks++)
				Fg[Solid[i][j][NZ]][SZP[ks]]=rho_p/36.0;//######
			}

		if ((szn==1) and (Solid[i][j][0]>0))
			{
			rho_t=df_szn;
			for (int ks=0;ks<=4;ks++)
				rho_t+=Fg[Solid[i][j][0]][SZP[ks]];//######

			rho_p=(rho_t)/(1.0/6.0);//######
			Fg[Solid[i][j][0]][SZN[0]]=rho_p/18.0;//######			
			for (int ks=1;ks<=4;ks++)
				Fg[Solid[i][j][0]][SZN[ks]]=rho_p/36.0;//######
			}
		       
		}

		
if ((sxp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
        for (int k=0;k<=NZ;k++)
			if (Solid[nx_l-1][j][k]>0)
			{
			rho_t=df_sxp;
			for (int ks=0;ks<=4;ks++)
				rho_t+=Fg[Solid[nx_l-1][j][k]][SXN[ks]];//######

			rho_p=(rho_t)/(1.0/6.0);//######
			Fg[Solid[nx_l-1][j][k]][SXP[0]]=rho_p/18.0;//######			
			for (int ks=1;ks<=4;ks++)
				Fg[Solid[nx_l-1][j][k]][SXP[ks]]=rho_p/36.0;//######
			}
		
			

if ((sxn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
			if (Solid[0][j][k]>0)
			{
			rho_t=df_sxn;
			for (int ks=0;ks<=4;ks++)
				rho_t+=Fg[Solid[0][j][k]][SXP[ks]];//######

			rho_p=(rho_t)/(1.0/6.0);//######
			Fg[Solid[0][j][k]][SXN[0]]=rho_p/18.0;//######			
			for (int ks=1;ks<=4;ks++)
				Fg[Solid[0][j][k]][SXN[ks]]=rho_p/36.0;//######

			}
       
}


if (Sub_BC_psi==1)
{

if ((syp-1)*(syn-1)==0)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)	
	  	{
		if ((syp==1) && (Solid[i][NY][k]>0))
			{
			if (Solid[i][NY-1][k]>0)
					{
					srho_syp=df_syp*dx+rho_r[Solid[i][NY-1][k]];
					for (int ks=0;ks<19;ks++)
					Fg[Solid[i][NY][k]][ks]=feq_psi(ks,srho_syp,u[Solid[i][NY][k]])+fg[Solid[i][NY-1][k]][ks]-feq_psi(ks,rho_r[Solid[i][NY-1][k]],u[Solid[i][NY-1][k]]);
					}
				else	
					for (int ks=0;ks<19;ks++)
					{
					srho_syp=df_syp*dx+ref_psi;
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					Fg[Solid[i][NY][k]][ks]=feq_psi(ks,srho_syp,u_ls);
					}

			}
		        
		if ((syn==1) && (Solid[i][0][k]>0))
		    {
			if (Solid[i][1][k]>0)
					{
					srho_syn=rho_r[Solid[i][1][k]]-df_syn*dx;
					for (int ks=0;ks<19;ks++)
					Fg[Solid[i][0][k]][ks]=feq_psi(ks,srho_syn,u[Solid[i][0][k]])+fg[Solid[i][1][k]][ks]-feq_psi(ks,rho_r[Solid[i][1][k]],u[Solid[i][1][k]]);
					}
				else
					for (int ks=0;ks<19;ks++)
					{
					srho_syn=ref_psi-df_syn*dx;
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					Fg[Solid[i][0][k]][ks]=feq_psi(ks,srho_syn,u_ls);
					}

			}   
		
		}

if ((szp-1)*(szn-1)==0)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
		{
		if ((szp==1) and (Solid[i][j][NZ]>0))
			{
			if (Solid[i][j][NZ-1]>0)
					{
					srho_szp=df_szp*dx+rho_r[Solid[i][j][NZ-1]];
					for (int ks=0;ks<19;ks++)
					Fg[Solid[i][j][NZ]][ks]=feq_psi(ks,srho_szp,u[Solid[i][j][NZ]])+fg[Solid[i][j][NZ-1]][ks]-feq_psi(ks,rho_r[Solid[i][j][NZ-1]],u[Solid[i][j][NZ-1]]);
					}
				else
					for (int ks=0;ks<19;ks++)
					{
					srho_szp=df_szp*dx+ref_psi;
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					Fg[Solid[i][j][NZ]][ks]=feq_psi(ks,srho_szp,u_ls);
					}

			}

		if ((szn==1) and (Solid[i][j][0]>0))
			{
			if (Solid[i][j][1]>0)
					{
					srho_szn=rho_r[Solid[i][j][1]]-df_szn*dx;
					for (int ks=0;ks<19;ks++)
					Fg[Solid[i][j][0]][ks]=feq_psi(ks,srho_szn,u[Solid[i][j][0]])+fg[Solid[i][j][1]][ks]-feq_psi(ks,rho_r[Solid[i][j][1]],u[Solid[i][j][1]]);
					}
				else
					for (int ks=0;ks<19;ks++)
					{
					srho_szn=ref_psi-df_szn*dx;
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					Fg[Solid[i][j][0]][ks]=feq_psi(ks,srho_szn,u_ls);
					}
			}
		       
		}

if ((sxp==1) && (rank==mpi_size-1)) 
for (int j=0;j<=NY;j++)
        for (int k=0;k<=NZ;k++) 
			if (Solid[nx_l-1][j][k]>0)
			{
			if (Solid[nx_l-2][j][k]>0)
					{
					srho_sxp=df_sxp*dx+rho_r[Solid[nx_l-2][j][k]];
					for (int ks=0;ks<19;ks++)
					Fg[Solid[nx_l-1][j][k]][ks]=feq_psi(ks,srho_sxp,u[Solid[nx_l-1][j][k]])+fg[Solid[nx_l-2][j][k]][ks]-feq_psi(ks,rho_r[Solid[nx_l-2][j][k]],u[Solid[nx_l-2][j][k]]);
					}
				else
					for (int ks=0;ks<19;ks++)
					{
					srho_sxp=df_sxp*dx+ref_psi;
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					Fg[Solid[nx_l-1][j][k]][ks]=feq_psi(ks,srho_sxp,u_ls);
					}
			}
		
			

if ((sxn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++) 
			if (Solid[0][j][k]>0)
			{
			if (Solid[1][j][k]>0)
					{
					srho_sxn=rho_r[Solid[1][j][k]]-df_sxn*dx;
					for (int ks=0;ks<19;ks++)
					Fg[Solid[0][j][k]][ks]=feq_psi(ks,srho_sxn,u[Solid[0][j][k]])+fg[Solid[1][j][k]][ks]-feq_psi(ks,rho_r[Solid[1][j][k]],u[Solid[1][j][k]]);
					}
				else
					for (int ks=0;ks<19;ks++)
					{
					srho_sxn=ref_psi-df_sxn*dx;
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					Fg[Solid[0][j][k]][ks]=feq_psi(ks,srho_sxn,u_ls);
					}
			

			}
      



}



}



void Solute_Constant_BC(int sxp,double srho_sxp,int sxn, double srho_sxn,int syp,double srho_syp,int syn,double srho_syn,int szp,double srho_szp,int szn,double srho_szn,double** Fg,int*** Solid,double* rho_r,double** u,double** fg)
{
int Q=19;
int rank = MPI :: COMM_WORLD . Get_rank ();
int mpi_size=MPI :: COMM_WORLD . Get_size ();
double rho_t,rho_p;
double u_ls[3];

if (Sub_BC_psi==0)
{
if ((syp-1)*(syn-1)==0)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)	
	  	{
		if ((syp==1) and (Solid[i][NY][k]>0))
			{
			rho_t=0;
			for (int ks=0;ks<19;ks++)
				rho_t+=Fg[Solid[i][NY][k]][ks];

			rho_t-=(Fg[Solid[i][NY][k]][SYP[0]]+Fg[Solid[i][NY][k]][SYP[1]]+Fg[Solid[i][NY][k]][SYP[2]]+Fg[Solid[i][NY][k]][SYP[3]]+Fg[Solid[i][NY][k]][SYP[4]]);
			rho_p=(srho_syp-rho_t)/(1.0/6.0);
			Fg[Solid[i][NY][k]][SYP[0]]=rho_p/18.0;			
			for (int ks=1;ks<=4;ks++)
				Fg[Solid[i][NY][k]][SYP[ks]]=rho_p/36.0;
			}
		        
		if ((syn==1) and (Solid[i][0][k]>0))
		        {
			rho_t=0;
			for (int ks=0;ks<19;ks++)
				rho_t+=Fg[Solid[i][0][k]][ks];

			rho_t-=(Fg[Solid[i][0][k]][SYN[0]]+Fg[Solid[i][0][k]][SYN[1]]+Fg[Solid[i][0][k]][SYN[2]]+Fg[Solid[i][0][k]][SYN[3]]+Fg[Solid[i][0][k]][SYN[4]]);
			rho_p=(srho_syn-rho_t)/(1.0/6.0);
			Fg[Solid[i][0][k]][SYN[0]]=rho_p/18.0;			
			for (int ks=1;ks<=4;ks++)
				Fg[Solid[i][0][k]][SYN[ks]]=rho_p/36.0;
			}
		
		}
	        
if ((szp-1)*(szn-1)==0)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
		{
		if ((szp==1) and (Solid[i][j][NZ]>0))
			{
			rho_t=0;
			for (int ks=0;ks<19;ks++)
				rho_t+=Fg[Solid[i][j][NZ]][ks];//######

			rho_t-=(Fg[Solid[i][j][NZ]][SZP[0]]+Fg[Solid[i][j][NZ]][SZP[1]]+Fg[Solid[i][j][NZ]][SZP[2]]+Fg[Solid[i][j][NZ]][SZP[3]]+Fg[Solid[i][j][NZ]][SZP[4]]);//######
			rho_p=(srho_szp-rho_t)/(1.0/6.0);//######
			Fg[Solid[i][j][NZ]][SZP[0]]=rho_p/18.0;//######			
			for (int ks=1;ks<=4;ks++)
				Fg[Solid[i][j][NZ]][SZP[ks]]=rho_p/36.0;//######
			}

		if ((szn==1) and (Solid[i][j][0]>0))
			{
			rho_t=0;
			for (int ks=0;ks<19;ks++)
				rho_t+=Fg[Solid[i][j][0]][ks];//######

			rho_t-=(Fg[Solid[i][j][0]][SZN[0]]+Fg[Solid[i][j][0]][SZN[1]]+Fg[Solid[i][j][0]][SZN[2]]+Fg[Solid[i][j][0]][SZN[3]]+Fg[Solid[i][j][0]][SZN[4]]);//######
			rho_p=(srho_szn-rho_t)/(1.0/6.0);//######
			Fg[Solid[i][j][0]][SZN[0]]=rho_p/18.0;//######			
			for (int ks=1;ks<=4;ks++)
				Fg[Solid[i][j][0]][SZN[ks]]=rho_p/36.0;//######
			}
		       
		}

		
if ((sxp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
        for (int k=0;k<=NZ;k++)
			if (Solid[nx_l-1][j][k]>0)
			{
			rho_t=0;
			for (int ks=0;ks<19;ks++)
				rho_t+=Fg[Solid[nx_l-1][j][k]][ks];//######

			rho_t-=(Fg[Solid[nx_l-1][j][k]][SXP[0]]+Fg[Solid[nx_l-1][j][k]][SXP[1]]+Fg[Solid[nx_l-1][j][k]][SXP[2]]+Fg[Solid[nx_l-1][j][k]][SXP[3]]+Fg[Solid[nx_l-1][j][k]][SXP[4]]);//######
			rho_p=(srho_sxp-rho_t)/(1.0/6.0);//######
			Fg[Solid[nx_l-1][j][k]][SXP[0]]=rho_p/18.0;//######			
			for (int ks=1;ks<=4;ks++)
				Fg[Solid[nx_l-1][j][k]][SXP[ks]]=rho_p/36.0;//######
			}
		
			

if ((sxn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
			if(Solid[0][j][k]>0)
			{
			rho_t=0;
			for (int ks=0;ks<19;ks++)
				rho_t+=Fg[Solid[0][j][k]][ks];//######

			rho_t-=(Fg[Solid[0][j][k]][SXN[0]]+Fg[Solid[0][j][k]][SXN[1]]+Fg[Solid[0][j][k]][SXN[2]]+Fg[Solid[0][j][k]][SXN[3]]+Fg[Solid[0][j][k]][SXN[4]]);//######
			rho_p=(srho_sxn-rho_t)/(1.0/6.0);//######
			Fg[Solid[0][j][k]][SXN[0]]=rho_p/18.0;//######			
			for (int ks=1;ks<=4;ks++)
				Fg[Solid[0][j][k]][SXN[ks]]=rho_p/36.0;//######

			}
       
}



if (Sub_BC_psi==1)
{
if ((syp-1)*(syn-1)==0)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)	
	  	{
		if ((syp==1) && (Solid[i][NY][k]>0))
			{
			if (Solid[i][NY-1][k]>0)
					{
					for (int ks=0;ks<19;ks++)
					Fg[Solid[i][NY][k]][ks]=feq_psi(ks,srho_syp,u[Solid[i][NY][k]])+fg[Solid[i][NY-1][k]][ks]-feq_psi(ks,rho_r[Solid[i][NY-1][k]],u[Solid[i][NY-1][k]]);
					}
				else	
					for (int ks=0;ks<19;ks++)
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					Fg[Solid[i][NY][k]][ks]=feq_psi(ks,srho_syp,u_ls);
					}

			}
		        
		if ((syn==1) && (Solid[i][0][k]>0))
		    {
			if (Solid[i][1][k]>0)
					{
					for (int ks=0;ks<19;ks++)
					Fg[Solid[i][0][k]][ks]=feq_psi(ks,srho_syn,u[Solid[i][0][k]])+fg[Solid[i][1][k]][ks]-feq_psi(ks,rho_r[Solid[i][1][k]],u[Solid[i][1][k]]);
					}
				else
					for (int ks=0;ks<19;ks++)
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					Fg[Solid[i][0][k]][ks]=feq_psi(ks,srho_syn,u_ls);
					}

			}   
		
		}

if ((szp-1)*(szn-1)==0)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
		{
		if ((szp==1) and (Solid[i][j][NZ]>0))
			{
			if (Solid[i][j][NZ-1]>0)
					{
					for (int ks=0;ks<19;ks++)
					Fg[Solid[i][j][NZ]][ks]=feq_psi(ks,srho_szp,u[Solid[i][j][NZ]])+fg[Solid[i][j][NZ-1]][ks]-feq_psi(ks,rho_r[Solid[i][j][NZ-1]],u[Solid[i][j][NZ-1]]);
					}
				else
					for (int ks=0;ks<19;ks++)
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					Fg[Solid[i][j][NZ]][ks]=feq_psi(ks,srho_szp,u_ls);
					}

			}

		if ((szn==1) and (Solid[i][j][0]>0))
			{
			if (Solid[i][j][1]>0)
					{
					for (int ks=0;ks<19;ks++)
					Fg[Solid[i][j][0]][ks]=feq_psi(ks,srho_szn,u[Solid[i][j][0]])+fg[Solid[i][j][1]][ks]-feq_psi(ks,rho_r[Solid[i][j][1]],u[Solid[i][j][1]]);
					}
				else
					for (int ks=0;ks<19;ks++)
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					Fg[Solid[i][j][0]][ks]=feq_psi(ks,srho_szn,u_ls);
					}
			}
		       
		}

if ((sxp==1) && (rank==mpi_size-1)) 
for (int j=0;j<=NY;j++)
        for (int k=0;k<=NZ;k++) 
			if (Solid[nx_l-1][j][k]>0)
			{
			if (Solid[nx_l-2][j][k]>0)
					{
					for (int ks=0;ks<19;ks++)
					Fg[Solid[nx_l-1][j][k]][ks]=feq_psi(ks,srho_sxp,u[Solid[nx_l-1][j][k]])+fg[Solid[nx_l-2][j][k]][ks]-feq_psi(ks,rho_r[Solid[nx_l-2][j][k]],u[Solid[nx_l-2][j][k]]);
					}
				else
					for (int ks=0;ks<19;ks++)
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					Fg[Solid[nx_l-1][j][k]][ks]=feq_psi(ks,srho_sxp,u_ls);
					}
			}
		
			

if ((sxn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++) 
			if (Solid[0][j][k]>0)
			{
			if (Solid[1][j][k]>0)
					{
					for (int ks=0;ks<19;ks++)
					Fg[Solid[0][j][k]][ks]=feq_psi(ks,srho_sxn,u[Solid[0][j][k]])+fg[Solid[1][j][k]][ks]-feq_psi(ks,rho_r[Solid[1][j][k]],u[Solid[1][j][k]]);
					}
				else
					for (int ks=0;ks<19;ks++)
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					Fg[Solid[0][j][k]][ks]=feq_psi(ks,srho_sxn,u_ls);
					}
			

			}
       

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
		        if (Solid[i][NY-1][k]>0)
		                 F[Solid[i][NY][k]][ks]=feq(ks,rho[Solid[i][NY-1][k]],u_yp);
		         else
		                 F[Solid[i][NY][k]][ks]=feq(ks,1.0,u_yp);
		if ((yn==1) && (Solid[i][0][k]>0))
		        if (Solid[i][1][k]>0)
		                 F[Solid[i][0][k]][ks]=feq(ks,rho[Solid[i][1][k]],u_yn);
		         else
		                 F[Solid[i][0][k]][ks]=feq(ks,1.0,u_yn);      
		}

if ((zp-1)*(zn-1)==0)		
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
		for (int ks=0;ks<Q;ks++)
		{
		if ((zp==1) && (Solid[i][j][NZ]>0)) 
		        if (Solid[i][j][NZ-1]>0)
		                F[Solid[i][j][NZ]][ks]=feq(ks,rho[Solid[i][j][NZ-1]],u_zp);
		        else
		                F[Solid[i][j][NZ]][ks]=feq(ks,1.0,u_zp); 
		if ((zn==1) && (Solid[i][j][0]>0))
		        if (Solid[i][j][1]>0)
		                F[Solid[i][j][0]][ks]=feq(ks,rho[Solid[i][j][1]],u_zn);
		        else
		                F[Solid[i][j][0]][ks]=feq(ks,1.0,u_zn);
		}


if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
		        if (Solid[nx_l-2][j][k]>0)
			        F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho[Solid[nx_l-2][j][k]],u_xp);
			else
			         F[Solid[nx_l-1][j][k]][ks]=feq(ks,1.0,u_xp);
			



if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
		        if (Solid[1][j][k]>0)
		                F[Solid[0][j][k]][ks]=feq(ks,rho[Solid[1][j][k]],u_xn);
		        else
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
			if (Solid[i][NY-1][k]>0) 
				F[Solid[i][NY][k]][ks]=feq(LR[ks],rho[Solid[i][NY-1][k]],u_yp)-F[Solid[i][NY][k]][LR[ks]]+feq(ks,rho[Solid[i][NY-1][k]],u_yp);
			else
				F[Solid[i][NY][k]][ks]=feq(LR[ks],1.0,u_yp)-F[Solid[i][NY][k]][LR[ks]]-feq(ks,1.0,u_yp);

			}


if (yn==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
		if (Solid[i][0][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
			if (e[ks][1]>0)
				if (Solid[i][1][k]>0) 
					F[Solid[i][0][k]][ks]=feq(LR[ks],rho[Solid[i][1][k]],u_yn)-F[Solid[i][0][k]][LR[ks]]+feq(ks,rho[Solid[i][1][k]],u_yn);
				else
					F[Solid[i][0][k]][ks]=feq(LR[ks],1.0,u_yn)-F[Solid[i][0][k]][LR[ks]]+feq(ks,1.0,u_yn);

			}



if (zp==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][NZ]>0)
		for (int ks=0;ks<Q;ks++)
			{
			if (e[ks][2]<0)
				if (Solid[i][j][NZ-1]>0) 
					F[Solid[i][j][NZ]][ks]=feq(LR[ks],rho[Solid[i][j][NZ-1]],u_zp)-F[Solid[i][j][NZ]][LR[ks]]+feq(ks,rho[Solid[i][j][NZ-1]],u_zp);
				else
					F[Solid[i][j][NZ]][ks]=feq(LR[ks],1.0,u_zp)-F[Solid[i][j][NZ]][LR[ks]]+feq(ks,1.0,u_zp);
		
			}



if (zn==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][0]>0)
		for (int ks=0;ks<Q;ks++)
			{
			if (e[ks][2]>0)
				if (Solid[i][j][1]>0) 
					F[Solid[i][j][0]][ks]=feq(LR[ks],rho[Solid[i][j][1]],u_zn)-F[Solid[i][j][0]][LR[ks]]+feq(ks,rho[Solid[i][j][1]],u_zn);
				else
					F[Solid[i][j][0]][ks]=feq(LR[ks],1.0,u_zn)-F[Solid[i][j][0]][LR[ks]]+feq(ks,1.0,u_zn);
	
			}




if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
	if (Solid[nx_l-1][k][j]>0)
		for (int ks=0;ks<Q;ks++)
			{
			if (e[ks][0]<0)
				if (Solid[nx_l-2][j][k]>0) 
					F[Solid[nx_l-1][j][k]][ks]=feq(LR[ks],rho[Solid[nx_l-2][j][k]],u_xp)-F[Solid[nx_l-1][j][k]][LR[ks]]+feq(ks,rho[Solid[nx_l-2][j][k]],u_xp);
				else
					F[Solid[nx_l-1][j][k]][ks]=feq(LR[ks],1.0,u_xp)-F[Solid[nx_l-1][j][k]][LR[ks]]+feq(ks,1.0,u_xp);
		
			}



if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
	if (Solid[0][j][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
			if (e[ks][0]>0)
				if (Solid[1][j][k]>0) 
					F[Solid[0][j][k]][ks]=feq(LR[ks],rho[Solid[1][j][k]],u_xn)-F[Solid[0][j][k]][LR[ks]]+feq(ks,rho[Solid[0][j][k]],u_xn);
				else
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
				if (Solid[i][NY-1][k]>0) 
					F[Solid[i][NY][k]][ks]=feq(ks,rho[Solid[i][NY-1][k]],u_yp)+f[Solid[i][NY-1][k]][ks]-feq(ks,rho[Solid[i][NY-1][k]],u[Solid[i][NY-1][k]]);
				else
					F[Solid[i][NY][k]][ks]=feq(ks,rho[Solid[i][NY-1][k]],u_yp);

			}


if (yn==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
	if (Solid[i][0][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
				if (Solid[i][1][k]>0) 
					F[Solid[i][0][k]][ks]=feq(ks,rho[Solid[i][1][k]],u_yn)+f[Solid[i][1][k]][ks]-feq(ks,rho[Solid[i][1][k]],u[Solid[i][1][k]]);
				else
					F[Solid[i][0][k]][ks]=feq(ks,rho[Solid[i][1][k]],u_yn);

			}



if (zp==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][NZ]>0)
		for (int ks=0;ks<Q;ks++)
			{
				if (Solid[i][j][NZ-1]>0) 
					F[Solid[i][j][NZ]][ks]=feq(ks,rho[Solid[i][j][NZ-1]],u_zp)+f[Solid[i][j][NZ-1]][ks]-feq(ks,rho[Solid[i][j][NZ-1]],u[Solid[i][j][NZ-1]]);
				else
					F[Solid[i][j][NZ]][ks]=feq(ks,rho[Solid[i][j][NZ-1]],u_zp);
		
			}



if (zn==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][0]>0)
		for (int ks=0;ks<Q;ks++)
			{
				if (Solid[i][j][1]>0) 
					F[Solid[i][j][0]][ks]=feq(ks,rho[Solid[i][j][1]],u_zn)+f[Solid[i][j][1]][ks]-feq(ks,rho[Solid[i][j][1]],u[Solid[i][j][1]]);
				else
					F[Solid[i][j][0]][ks]=feq(ks,rho[Solid[i][j][1]],u_zn);
	
			}




if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
	if (Solid[nx_l-1][k][j]>0)
		for (int ks=0;ks<Q;ks++)
			{
				if (Solid[nx_l-2][j][k]>0) 
					F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho[Solid[nx_l-2][j][k]],u_xp)+f[Solid[nx_l-2][j][k]][ks]-feq(ks,rho[Solid[nx_l-2][j][k]],u[Solid[nx_l-2][j][k]]);
				else
					F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho[Solid[nx_l-2][j][k]],u_xp);
		
			}



if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
	if (Solid[0][j][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
				if (Solid[1][j][k]>0) 
					F[Solid[0][j][k]][ks]=feq(ks,rho[Solid[1][j][k]],u_xn)+f[Solid[1][j][k]][ks]-feq(ks,rho[Solid[1][j][k]],u[Solid[1][j][k]]);
				else
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
		for (int ks=0;ks<Q;ks++)
			
			if (e[ks][1]<0)
			{
				if (Solid[i][NY-1][k]>0)
					{
					u_ls[0]=u[Solid[i][NY-1][k]][0];
					u_ls[1]=u[Solid[i][NY-1][k]][1];
					u_ls[2]=u[Solid[i][NY-1][k]][2];
					F[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls)-F[Solid[i][NY][k]][LR[ks]]+feq(LR[ks],rho_yp,u_ls);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls)-F[Solid[i][NY][k]][LR[ks]]+feq(LR[ks],rho_yp,u_ls);
					}
			}


if (yn==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
	if (Solid[i][0][k]>0)	
		for (int ks=0;ks<Q;ks++)
			if (e[ks][1]>0)
			{
				if (Solid[i][1][k]>0)
					{
					u_ls[0]=u[Solid[i][1][k]][0];
					u_ls[1]=u[Solid[i][1][k]][1];
					u_ls[2]=u[Solid[i][1][k]][2];
					F[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls)-F[Solid[i][0][k]][LR[ks]]+feq(LR[ks],rho_yn,u_ls);
					}
			else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls)-F[Solid[i][0][k]][LR[ks]]+feq(LR[ks],rho_yn,u_ls);
					}
			}

if (zp==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][NZ]>0)
		for (int ks=0;ks<Q;ks++)
		if (e[ks][2]<0)
			{			
			if (Solid[i][j][NZ-1]>0)
				{
				u_ls[0]=u[Solid[i][j][NZ-1]][0];
				u_ls[1]=u[Solid[i][j][NZ-1]][1];
				u_ls[2]=u[Solid[i][j][NZ-1]][2];
				F[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls)-F[Solid[i][j][NZ]][LR[ks]]+feq(LR[ks],rho_zp,u_ls);
				}
				else
				{
				u_ls[0]=0.0;
				u_ls[1]=0.0;u_ls[2]=0.0;
				F[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls)-F[Solid[i][j][NZ]][LR[ks]]+feq(LR[ks],rho_zp,u_ls);
				}
			}




if (zn==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][0]>0)
		for (int ks=0;ks<Q;ks++)
		if (e[ks][2]>0)
		{
			if (Solid[i][j][1]>0)
			{
			u_ls[0]=u[Solid[i][j][1]][0];
			u_ls[1]=u[Solid[i][j][1]][1];
			u_ls[2]=u[Solid[i][j][1]][2];
			F[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls)-F[Solid[i][j][0]][LR[ks]]+feq(LR[ks],rho_zn,u_ls);
			}
			else
			{
			u_ls[0]=0.0;
			u_ls[1]=0.0;u_ls[2]=0.0;
			F[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls)-F[Solid[i][j][0]][LR[ks]]+feq(LR[ks],rho_zn,u_ls);
			}
		
		}

if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
	if (Solid[nx_l-1][j][k]>0)
		for (int ks=0;ks<Q;ks++)			
			if (e[ks][0]<0)
			{
				if (Solid[nx_l-2][j][k]>0)
					{
					u_ls[0]=u[Solid[nx_l-2][j][k]][0];
					u_ls[1]=u[Solid[nx_l-2][j][k]][1];
					u_ls[2]=u[Solid[nx_l-2][j][k]][2];
					F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls)-F[Solid[nx_l-1][j][k]][LR[ks]]+feq(LR[ks],rho_xp,u_ls);
					}
				else	
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;
					u_ls[2]=0.0;
					F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls)-F[Solid[nx_l-1][j][k]][LR[ks]]+feq(LR[ks],rho_xp,u_ls);
					}

			}
		


if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
	if (Solid[0][j][k]>0)
		for (int ks=0;ks<Q;ks++)
			
			if (e[ks][0]>0)
			{
				if(Solid[1][j][k]>0)
					{
					u_ls[0]=u[Solid[1][j][k]][0];
					u_ls[1]=u[Solid[1][j][k]][1];
					u_ls[2]=u[Solid[1][j][k]][1];
					F[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls)-F[Solid[0][j][k]][LR[ks]]+feq(LR[ks],rho_xn,u_ls);
					//f[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;
					u_ls[2]=0.0;
					F[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls)-F[Solid[0][j][k]][LR[ks]]+feq(LR[ks],rho_xn,u_ls);
					}

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
/*
//=======================================	
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
//===================================
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
				out<<rbuf_rho[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]<<endl;

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

void output_psi_b(int m,double* psi,int MirX,int MirY,int MirZ,int mir,int*** Solid,char out_name[128])	
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
	name<<outputfile<<"LBM_"<<out_name<<"_"<<m<<".vtk";
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
				out<<rbuf_psi[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]<<endl;

	out.close();
				
	
	

/*
//========================================================================================
	ostringstream name2;
	name2<<"LBM_HS_"<<m<<".out";
	ofstream out2(name2.str().c_str());
	for (int j=0;j<=NY;j++)
		{
		if (Solid[2][j][1]>0)
			out2<<j<<" "<<rbuf_psi[2*(NY+1)*(NZ+1)+j*(NZ+1)+1]<<endl;
		else
			out2<<0.0<<endl;
		}
	
	out2.close();
//======================================================================================
*/
	}

	if (rank==root_rank)
		{		
		delete [] rbuf_psi;
		}
	delete [] nx_g;
	delete [] disp;
	delete [] psi_storage;

		
}

void output_psi(int m,double* psi,int MirX,int MirY,int MirZ,int mir,int*** Solid,char out_name[128])	
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
	name<<outputfile<<"LBM_"<<out_name<<"_"<<m<<".vtk";
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
	double error,vxl,vyl,vzl;
	double Q_l[3]={0.0,0.0,0.0};
	double Q_g[3]={0.0,0.0,0.0};
	
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

		Perm_l[0]=Q_l[0]/((NX+1)*(NY+1)*(NZ+1))*(niu_l)/(gx+dp);
		Perm_l[1]=Q_l[1]/((NX+1)*(NY+1)*(NZ+1))*(niu_l)/(gy+dp);
		Perm_l[2]=Q_l[2]/((NX+1)*(NY+1)*(NZ+1))*(niu_l)/(gz+dp);

		Perm_g[0]=Q_g[0]/((NX+1)*(NY+1)*(NZ+1))*(niu_g)/(gx+dp);
		Perm_g[1]=Q_g[1]/((NX+1)*(NY+1)*(NZ+1))*(niu_g)/(gy+dp);
		Perm_g[2]=Q_g[2]/((NX+1)*(NY+1)*(NZ+1))*(niu_g)/(gz+dp);


		vxl=(Q_l[0]+Q_g[0])/((NX+1)*(NY+1)*(NZ+1));
		vyl=(Q_l[1]+Q_g[1])/((NX+1)*(NY+1)*(NZ+1));
		vzl=(Q_l[2]+Q_g[2])/((NX+1)*(NY+1)*(NZ+1));
		Capillary=sqrt(vxl*vxl+vyl*vyl+vzl*vzl)*niu_l/CapA;

		
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




void Backup_init(double* rho, double** u, double** f,double* psi,double* rho_r, double* rho_b, double* rhor, double* rhob,char backup_rho[128], char backup_velocity[128], char backup_psi[128], char backup_f[128], char backup_fg[128], char backup_hs[128])
{	
      
        int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	
	
	double usqr,vsqr,eu;
	double c2,c4;
	
	rho0=1.0;dt=1.0/Zoom;dx=1.0/Zoom;
 	uMax=0.0;

	if (lattice_v==1)
		{dx=dx_input;dt=dt_input;}

	lat_c=dx/dt;
	c_s=lat_c/sqrt(3);
	c_s2=lat_c*lat_c/3;

	c2=lat_c*lat_c;c4=c2*c2;

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
	{
	
	M_c[0]=1.0;
	M_c[1]=lat_c*lat_c;
	M_c[2]=lat_c*lat_c*lat_c*lat_c;
	M_c[3]=lat_c;
	M_c[4]=lat_c*lat_c*lat_c;
	M_c[5]=lat_c;
	M_c[6]=lat_c*lat_c*lat_c;
	M_c[7]=lat_c;	
	M_c[8]=lat_c*lat_c*lat_c;
	M_c[9]=lat_c*lat_c;
	M_c[10]=lat_c*lat_c*lat_c*lat_c;
	M_c[11]=lat_c*lat_c;
	M_c[12]=lat_c*lat_c*lat_c*lat_c;
	M_c[13]=lat_c*lat_c;
	M_c[14]=lat_c*lat_c;
	M_c[15]=lat_c*lat_c;
	M_c[16]=lat_c*lat_c*lat_c;
	M_c[17]=lat_c*lat_c*lat_c;
	M_c[18]=lat_c*lat_c*lat_c;



	for (int i=0;i<19;i++)
		for (int j=0;j<3;j++)
		elat[i][j]=e[i][j]*lat_c;

	for (int i=0;i<19;i++)
		for (int j=0;j<19;j++)
		M[i][j]*=M_c[i];

	Comput_MI(M,MI);

	}


	psi_solid=ContactAngle_parameter;
	
	ostringstream name6;
	name6<<backup_hs<<"."<<rank<<".input";
	ostringstream name5;
	name5<<backup_fg<<"."<<rank<<".input";
 	ostringstream name4;
	name4<<backup_f<<"."<<rank<<".input";
	ostringstream name3;
	name3<<backup_psi<<"."<<rank<<".input";
 	ostringstream name2;
	name2<<backup_velocity<<"."<<rank<<".input";
	ostringstream name;
	name<<backup_rho<<"."<<rank<<".input";
	
	 
	ifstream fin;
	fin.open(name.str().c_str());
	
        	for(int i=1;i<=Count;i++)
        	        fin >> rho[i];
  
       fin.close();

       
   	fin.open(name6.str().c_str());
	
        	for(int i=1;i<=Count;i++)
        	        fin >> rhoh[i];
  
       fin.close();


	fin.open(name2.str().c_str());
	
        	for(int i=1;i<=Count;i++)
        	        fin >> u[i][0] >> u[i][1] >> u[i][2];
  
       fin.close();
       
       fin.open(name3.str().c_str());
	
        	for(int i=1;i<=Count;i++)
        	{
        	        fin >> psi[i];
        	        rho_r[i]=(psi[i]*rho[i]+rho[i])/2;
			rho_b[i]=rho[i]-rho_r[i];
			rhor[i]=0;
			rhob[i]=0;
		}
  
       fin.close();
       
       fin.open(name4.str().c_str());
	
        	for(int i=1;i<=Count;i++)
        	        fin >> f[i][0] >> f[i][1] >> f[i][2] >> f[i][3] >> f[i][4] >> f[i][5] >>f[i][6] >> f[i][7] >> f[i][8] >> f[i][9] >> f[i][10] >> f[i][11] >> f[i][12] >> f[i][13] >> f[i][14] >> f[i][15] >> f[i][16] >>f[i][17] >> f[i][18] ;
  
       fin.close();
       
      fin.open(name5.str().c_str());
	
        	for(int i=1;i<=Count;i++)
        	        fin >> fg[i][0] >> fg[i][1] >> fg[i][2] >> fg[i][3] >> fg[i][4] >> fg[i][5] >>fg[i][6] >> fg[i][7] >> fg[i][8] >> fg[i][9] >> fg[i][10] >> fg[i][11] >> fg[i][12] >> fg[i][13] >> fg[i][14] >> fg[i][15] >> fg[i][16] >>fg[i][17] >> fg[i][18] ;
  
       fin.close();
	
	
	
		       if (stab==1)
				{gxs=0;gys=0;gzs=0;}
			else
				{gxs=gx;gys=gy;gzs=gz;}
			
	
	 	
}




void Backup(int m,double* rho,double* psi, double** u, double** f)
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
	

	
	ostringstream name3;
	name3<<outputfile<<"LBM_Backup_Concentration_"<<m<<"."<<rank<<".input";
	//ofstream out;
	out.open(name3.str().c_str());
	
	for (int i=1;i<=Count;i++)
        		out<<psi[i]<<endl;
		
			
	out.close();
	
	
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
        
	
	ostringstream name5;
	name5<<outputfile<<"LBM_Backup_fg_"<<m<<"."<<rank<<".input";
	//ofstream out;
	out.open(name5.str().c_str());
	
	for (int i=1;i<=Count;i++)
	{
	        for (int j=0;j<19;j++)
        		out<<fg[i][j]<<" ";
        out<<endl;
        }
                
        out.close();


	ostringstream name6;
	name6<<outputfile<<"LBM_Backup_HS_"<<m<<"."<<rank<<".input";
	//ofstream out;
	out.open(name6.str().c_str());
	
	for (int i=1;i<=Count;i++)
        		out<<rhoh[i]<<endl;
		
			
	out.close();


}


