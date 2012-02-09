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

int preci=2;
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

/*
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
*/



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

void Read_Rock(int***,double***, double*,char[128],char[128]);

void tests();

void init_Sparse_read_rock_parallel(int*, int*);

void init(double*, double**, double**,double*, double*,double*, double*, double*,int*);

void periodic_streaming(double** ,double** ,int* ,int***,int*, int*,double*, double**);

void standard_bounceback_boundary(int,double**);

void collision(double*,double** ,double** ,double**, double*, double*, double*, double* ,double* , int* ,int***,int* ,int*);

void comput_macro_variables( double* ,double**,double** ,double** ,double** ,double*, double*, double*, double*, double* ,int* );

double Error(double** ,double** ,double*, double*);

void boundary_velocity(int,double,int , double ,int ,double ,int ,double ,int ,double ,int ,double ,double** ,double**,double* ,double** ,int*** );

void boundary_pressure(int ,double ,int , double ,int ,double ,int ,double ,int ,double ,int ,double ,double** ,double**,double** ,double* ,int*** );

void output_velocity(int ,double* ,double** ,int ,int ,int ,int ,int*** );

void output_density(int ,double* ,int ,int ,int ,int ,int*** );	

void Geometry(int*** );	

void output_velocity_b(int ,double* ,double** ,int ,int ,int ,int ,int*** );

void output_psi(int ,double* ,int ,int ,int ,int ,int*** );

void output_psi_b(int ,double* ,int ,int ,int ,int ,int*** );

void output_density_b(int ,double* ,int ,int ,int ,int ,int*** );	

void Geometry_b(int*** );

double Comput_Perm(double* ,double** ,double* ,double* ,int,int* );

double Comput_Saturation(double* ,int***,int*);

double S[19];

void Comput_MI(double[19][19], double[19][19]);

int inverse(mat &a);

double feq(int,double, double[3]);

void Suppliment(int*,int***);

void Backup_init(double* , double** , double** ,double* ,double* , double* , double*, double* );

void Backup(int ,double*, double*, double**, double**,double*,double*);

void Parallelize_Geometry();

void pressure_bodyforce_change();

void Comput_Perm_LOCAL(double* ,double** ,double* ,double* ,int );


const int e[19][3]=
{{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},{1,0,1},
{-1,0,-1},{1,0,-1},{-1,0,1},{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1}};

double elat[19][3]=
{{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},{1,0,1},
{-1,0,-1},{1,0,-1},{-1,0,1},{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1}};


const double w[19]={1.0/3.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};


//========*************************===============
//int LR[19]={0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15};
const int LR[19]={0,2,1,4,3,6,5,8,7,10,9,12,11,14,13,16,15,18,17};
const int FRP[19]={0,0,0,0,0,0,0,1,0,2,0,3,0,4,0,0,0,0,0};
const int FLN[19]={0,0,0,0,0,0,0,0,1,0,2,0,3,0,4,0,0,0,0};
const int RP[5]={1,7,9,11,13};
const int LN[5]={2,8,10,12,14};
//=========================================

//======FOR==SWAP==STREAMING===========
//SWAPE [LOCAL]-->LADD

const int SWAPE_INV[19]={0,2,4,6,8,10,12,14,16,18,1,3,5,7,9,11,13,15,17};



int n,nx_l,n_max,in_BC,PerDir,freRe,freDe,freVe,frePsi,Par_Geo,Par_nx,Par_ny,Par_nz;
int Zoom,lattice_v,in_psi_BC,par_per_x,par_per_y,par_per_z;


int wr_per,pre_xp,pre_xn,pre_yp,pre_yn,pre_zp,pre_zn,stab,stab_time,fre_backup,psi_xp,psi_xn,psi_yp;
int psi_yn,per_xn,per_yn,per_zn,ini_buf;
int vel_xp,vel_xn,vel_yp,vel_yn,vel_zp,vel_zn,Sub_BC,Out_Mode,mode_backup_ini,psi_zp,psi_zn,per_xp,per_yp,per_zp;
double in_vis,p_xp,p_xn,p_yp,p_yn,p_zp,p_zn,niu_l,niu_g,ContactAngle_parameter,CapA;
double inivx,inivy,inivz,v_xp,v_xn,v_yp,v_yn,v_zp,v_zn,Re_l,Re_g,Capillary,ini_Sat,var_rho;
double error_Per,Permeability,psi_solid,S_l,gxs,gys,gzs,c_s,c_s2,dx_input,dt_input,lat_c;


char outputfile[128]="./";
int NCHAR=128;
	char     filename[128], dummy[128+1],filenamepsi[128], backup_rho[128], backup_velocity[128],backup_psi[128],backup_f[128];
	int      dummyInt;
	
int*** Solid;
float*** Psi_local;	
char pfix[128];
int decbin;


int pressure_change,pre_chan_pb,pre_chan_1,pre_chan_2,pre_chan_3,pre_chan_4,pre_chan_5;
double pre_chan_pn1,pre_chan_pn2,pre_chan_pn3,pre_chan_pn4,pre_chan_pn5;
double pre_chan_pp1,pre_chan_pp2,pre_chan_pp3,pre_chan_pp4,pre_chan_pp5;
double pre_chan_f1,pre_chan_f2,pre_chan_f3,pre_chan_f4,pre_chan_f5;
double rel_perm_psi;

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
double Per_l_LOCAL[3],Per_g_LOCAL[3];


double v_max,error_Per;

        strcpy(pfix,"./");
        if (argc>2)
                strcpy(pfix,argv[2]);
        


	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	if (rank==0)
	{
	ifstream fin(argv[1]);
	                                                fin.getline(dummy, NCHAR);
	fin >> filename;				fin.getline(dummy, NCHAR);
	fin >> filenamepsi;				fin.getline(dummy, NCHAR);
	fin >> NX >> NY >> NZ;				fin.getline(dummy, NCHAR);
	fin >> n_max;					fin.getline(dummy, NCHAR);
	fin >> reso;					fin.getline(dummy, NCHAR);
	fin >> in_BC;					fin.getline(dummy, NCHAR);
	fin >> in_psi_BC;				fin.getline(dummy, NCHAR);
	fin >> gx >> gy >> gz;				fin.getline(dummy, NCHAR);
	fin >> pre_xp >> p_xp >> pre_xn >> p_xn;	fin.getline(dummy, NCHAR);
	fin >> pre_yp >> p_yp >> pre_yn >> p_yn;	fin.getline(dummy, NCHAR);
	fin >> pre_zp >> p_zp >> pre_zn >> p_zn;	fin.getline(dummy, NCHAR);
	fin >> vel_xp >> v_xp >> vel_xn >> v_xn;	fin.getline(dummy, NCHAR);
	fin >> vel_yp >> v_yp >> vel_yn >> v_yn;	fin.getline(dummy, NCHAR);
	fin >> vel_zp >> v_zp >> vel_zn >> v_zn;	fin.getline(dummy, NCHAR);
	fin >> psi_xp >> psi_xn;			fin.getline(dummy, NCHAR);
	fin >> psi_yp >> psi_yn;			fin.getline(dummy, NCHAR);
	fin >> psi_zp >> psi_zn;			fin.getline(dummy, NCHAR);
	fin >> niu_l;					fin.getline(dummy, NCHAR);
	fin >> niu_g;					fin.getline(dummy, NCHAR);
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
							fin.getline(dummy, NCHAR);
	fin >> lattice_v >> dx_input >> dt_input;	fin.getline(dummy, NCHAR);
	fin >> outputfile;				fin.getline(dummy, NCHAR);
	fin >> Sub_BC;					fin.getline(dummy, NCHAR);
	fin >> stab >> stab_time;			fin.getline(dummy, NCHAR);
	fin >> ini_Sat;                                 fin.getline(dummy, NCHAR);
	fin >> ini_buf;                                 fin.getline(dummy, NCHAR);
	fin >> rel_perm_psi;				fin.getline(dummy, NCHAR);
	fin >> par_per_x >> par_per_y >>par_per_z;	fin.getline(dummy, NCHAR);
	fin >> per_xp >> per_xn;			fin.getline(dummy, NCHAR);
	fin >> per_yp >> per_yn;			fin.getline(dummy, NCHAR);
	fin >> per_zp >> per_zn;			fin.getline(dummy, NCHAR);
	                                      		fin.getline(dummy, NCHAR);
	fin >> fre_backup;                        	fin.getline(dummy, NCHAR);
	fin >>mode_backup_ini;                		fin.getline(dummy, NCHAR);
							fin.getline(dummy, NCHAR);
	fin >> decbin;					fin.getline(dummy, NCHAR);
							fin.getline(dummy, NCHAR);
	fin >> pressure_change >>pre_chan_pb;		fin.getline(dummy, NCHAR);
	fin >> pre_chan_1>> pre_chan_pn1 >> pre_chan_pp1>>pre_chan_f1;	fin.getline(dummy, NCHAR);
	fin >> pre_chan_2>> pre_chan_pn2 >> pre_chan_pp2>>pre_chan_f2;	fin.getline(dummy, NCHAR);
	fin >> pre_chan_3>> pre_chan_pn3 >> pre_chan_pp3>>pre_chan_f3;	fin.getline(dummy, NCHAR);
	fin >> pre_chan_4>> pre_chan_pn4 >> pre_chan_pp4>>pre_chan_f4;	fin.getline(dummy, NCHAR);
	fin >> pre_chan_5>> pre_chan_pn5 >> pre_chan_pp5>>pre_chan_f5;	fin.getline(dummy, NCHAR);



//	fin >> backup_rho;                        	fin.getline(dummy, NCHAR);
//	fin >> backup_velocity;                		fin.getline(dummy, NCHAR);
//	fin >> backup_psi;                        	fin.getline(dummy, NCHAR);
//	fin >> backup_f;                        	fin.getline(dummy, NCHAR);
	fin.close();
	
	//cout<<par_per_x<<"    asdfa    "<<backup_f<<endl;
	NX=NX-1;NY=NY-1;NZ=NZ-1;
	}

	MPI_Bcast(&filename,128,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&NX,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&var_rho,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&NY,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&ini_buf,1,MPI_INT,0,MPI_COMM_WORLD);
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
	MPI_Bcast(&dt_input,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&ini_Sat,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Bcast(&in_psi_BC,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&psi_xp,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&psi_xn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&psi_yp,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&psi_yn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&psi_zp,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&psi_zn,1,MPI_INT,0,MPI_COMM_WORLD);

	
	MPI_Bcast(&par_per_x,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&par_per_y,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&par_per_z,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&per_yp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&per_xp,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&per_yn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&per_zp,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&per_zn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&per_xn,1,MPI_INT,0,MPI_COMM_WORLD);


	MPI_Bcast(&pressure_change,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&pre_chan_pb,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_chan_1,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&pre_chan_pn1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_chan_pp1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&pre_chan_f1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_chan_2,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&pre_chan_pn2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_chan_pp2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&pre_chan_f2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_chan_3,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&pre_chan_pn3,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_chan_pp3,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&pre_chan_f3,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_chan_4,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&pre_chan_pn4,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_chan_pp4,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&pre_chan_f4,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_chan_5,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&pre_chan_pn5,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&pre_chan_pp5,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&pre_chan_f5,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&rel_perm_psi,1,MPI_DOUBLE,0,MPI_COMM_WORLD);








mirX=0;mirY=0;mirZ=0;
mir=1;Zoom=1;


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


	
	int*  Sl;
	int*  Sr;

	
	Parallelize_Geometry();
	
	Sl = new int[(NY+1)*(NZ+1)];
	Sr = new int[(NY+1)*(NZ+1)];

	
	MPI_Barrier(MPI_COMM_WORLD);
		

	init_Sparse_read_rock_parallel(Sl,Sr);
	MPI_Barrier(MPI_COMM_WORLD);


	
	//***************************************************
	//WARRING: SPARSE MATRIX STARTS FROM INDEX 1 NOT 0!!!
	//***************************************************

	Permia = new double[3];
	rho = new double[Count+1];
	rho_r = new double[Count+1];
	rho_b = new double[Count+1];
	rhor = new double[Count+1];
	rhob = new double[Count+1];
	psi = new double[Count+1];
	//forcex = new double[Count+1];
	//forcey = new double[Count+1];
	//forcez = new double[Count+1];
	u = new double*[Count+1];
	u[0] = new double[(Count+1)*3];
	        for (int i=1;i<=Count;i++)
	              u[i] = u[i-1]+3;
	      
	      
	f = new double*[Count+1];
	f[0] =new double[(Count+1)*19];
	        for (int i=1;i<=Count;i++)
	                f[i] = f[i-1]+19;
	        
	        
	F = new double*[Count+1];
	F[0] =new double[(Count+1)*19];
	for (int i=1;i<=Count;i++)
		F[i] = F[i-1]+19;
	
	
	u0 = new double*[Count+1];
	u0[0] = new double[(Count+1)*3];
	        for (int i=1;i<=Count;i++)
		u0[i] = u0[i-1]+3;
	
	
	SupInv = new int[Count+1];

/*	for (int i=0;i<=Count;i++)
		{
		u[i] = new double[3];
		f[i] = new double[19];
		u0[i] = new double[3];
		F[i] = new double[19];
		}
*/
	Comput_MI(M,MI);
	
	Suppliment(SupInv,Solid);

	MPI_Barrier(MPI_COMM_WORLD);
	
	if (Out_Mode==1)
		Geometry(Solid);
	else
		Geometry_b(Solid);

	if (mode_backup_ini==0)
	        init(rho,u,f,psi,rho_r,rho_b,rhor, rhob,SupInv);
	else
	      Backup_init(rho, u, f,psi,rho_r, rho_b, rhor, rhob);  

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
char FileName5[128];strcpy(FileName5,outputfile);

strcat(FileName,"Results.txt");
strcat(FileName2,"Relative_Permeability_Component1.txt");
strcat(FileName4,"Relative_Permeability_Component2.txt");
strcat(FileName3,"bodyforce.txt");
strcat(FileName5,"Velocity_ave_max.txt");
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

	fins.open(FileName5,ios::out);
	fins.close();

	
	
	for(n=0;n<=n_max;n++)
	{
	
	if ((stab==1) and (n==stab_time))
		{gxs=gx;gys=gy;gzs=gz;}
			
	//cout<<"       qgragf       "<<endl;
	collision(rho,u,f,F,psi,rho_r,rho_b,rhor,rhob,SupInv,Solid,Sl,Sr);

	if (pressure_change>0)
		pressure_bodyforce_change();
	//==========================PRESSURE AND BODYFORCE CHANGE ALONG WITH TIME ==================
	if (pressure_change>0)
	{
		if (pressure_change==1)
		{
		if (n==pre_chan_1)
			if (pre_chan_pb==1)
			{
			p_xn=pre_chan_pn1;
			p_xp=pre_chan_pp1;
			}
			else
				gxs=pre_chan_f1;

		if (n==pre_chan_1)
			if (pre_chan_pb==1)
			{
			p_xn=pre_chan_pn1;
			p_xp=pre_chan_pp1;
			}
			else
				gxs=pre_chan_f1;

		if (n==pre_chan_1)
			if (pre_chan_pb==1)
			{
			p_xn=pre_chan_pn1;
			p_xp=pre_chan_pp1;
			}
			else
				gxs=pre_chan_f1;

		if (n==pre_chan_1)
			if (pre_chan_pb==1)
			{
			p_xn=pre_chan_pn1;
			p_xp=pre_chan_pp1;
			}
			else
				gxs=pre_chan_f1;

		if (n==pre_chan_1)
			if (pre_chan_pb==1)
			{
			p_xn=pre_chan_pn1;
			p_xp=pre_chan_pp1;
			}
			else
				gxs=pre_chan_f1;

	 
		}

	}
	//==========================================================================================





	//periodic_streaming(f,F,SupInv,Solid,Sl,Sr,rho,u);	
	
	if ((stab==0) or ((stab==1) and (n>=stab_time)))
	{
	if ((1-pre_xp)*(1-pre_xn)*(1-pre_yp)*(1-pre_yn)*(1-pre_zp)*(1-pre_zn)==0)
		boundary_pressure(pre_xp,p_xp,pre_xn,p_xn,pre_yp,p_yp,pre_yn,p_yn,pre_zp,p_zp,pre_zn,p_zn,f,F,u,rho,Solid);

	
	if ((1-vel_xp)*(1-vel_xn)*(1-vel_yp)*(1-vel_yn)*(1-vel_zp)*(1-vel_zn)==0)
		boundary_velocity(vel_xp,v_xp,vel_xn,v_xn,vel_yp,v_yp,vel_yn,v_yn,vel_zp,v_zp,vel_zn,v_zn,f,F,rho,u,Solid);
	}
  		 
		comput_macro_variables(rho,u,u0,f,F,rho_r,rho_b,rhor,rhob,psi,SupInv);


	//if ((1-pre_xp)*(1-pre_xn)*(1-pre_yp)*(1-pre_yn)*(1-pre_zp)*(1-pre_zn)==0)
	//	boundary_pressure(pre_xp,p_xp,pre_xn,p_xn,pre_yp,p_yp,pre_yn,p_yn,pre_zp,p_zp,pre_zn,p_zn,f,u,rho,Solid);
		
	//if ((1-vel_xp)*(1-vel_xn)*(1-vel_yp)*(1-vel_yn)*(1-vel_zp)*(1-vel_zn)==0)
	//        boundary_velocity(vel_xp,v_xp,vel_xn,v_xn,vel_yp,v_yp,vel_yn,v_yn,vel_zp,v_zp,vel_zn,v_zn,f,rho,u,Solid);
	
	if(n%freRe==0)
		{       
			
			if (rank==0)
			{
			ofstream fin(FileName,ios::out);
			fin<<"The"<<n-freRe<<"th computation result:"<<endl;
			Re_l=u_ave*(NY+1)/niu_l;Re_g=u_ave*(NY+1)/niu_g;
		        fin<<"The Maximum velocity is: "<<setprecision(6)<<u_max<<"   Re_l="<<Re_l<<"   Re_g="<<Re_g<<endl;
			fin<<"Courant Number="<<u_max*dt/dx<<"	 Capillary Num="<<Capillary<<endl;
			fin<<"The max relative error of velocity is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			fin<<"The relative permeability of component 1 is "<<Per_l[0]*reso*reso*1000/Permeability<<", "<<Per_l[1]*reso*reso*1000/Permeability<<", "<<Per_l[2]*reso*reso*1000/Permeability<<endl;
			fin<<"The relative permeability of component 2 is "<<Per_g[0]*reso*reso*1000/Permeability<<", "<<Per_g[1]*reso*reso*1000/Permeability<<", "<<Per_g[2]*reso*reso*1000/Permeability<<endl;
			fin<<"The LOCAL relative permeability of component 1 is "<<Per_l_LOCAL[0]*reso*reso*1000/Permeability<<", "<<Per_l_LOCAL[1]*reso*reso*1000/Permeability<<", "<<Per_l_LOCAL[2]*reso*reso*1000/Permeability<<endl;
			fin<<"The LOCAL relative permeability of component 2 is "<<Per_g_LOCAL[0]*reso*reso*1000/Permeability<<", "<<Per_g_LOCAL[1]*reso*reso*1000/Permeability<<", "<<Per_g_LOCAL[2]*reso*reso*1000/Permeability<<endl;
			fin<<"Satuation of Component 1: "<<S_l<<", "<<"The satuation of Component 2: "<<1-S_l<<endl;
			fin<<"The relative error of permiability computing is: "<<error_Per<<endl;
			fin<<"Elapsed time is "<< the<<"h"<<tme<<"m"<<tse<<"s"<<endl;
			fin<<"The expected completion time is "<<th<<"h"<<tm<<"m"<<ts<<"s"<<endl;
			fin<<endl;
			fin.close();
			}
			
			
			error=Error(u,u0,&u_max,&u_ave);
			if (u_max>=10.0)	U_max_ref+=1;
			error_Per=Comput_Perm(psi,u,Per_l,Per_g,PerDir,SupInv);
			Comput_Perm_LOCAL(psi,u,Per_l_LOCAL,Per_g_LOCAL,PerDir);
			S_l=Comput_Saturation(psi,Solid,SupInv);
			

			if (rank==0)
			{
			ofstream fin(FileName,ios::app);
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
		//===============================================================================================
			fin<<"The max relative error of velocity is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			fin<<"The relative permeability of component 1 is "<<Per_l[0]*reso*reso*1000/Permeability<<", "<<Per_l[1]*reso*reso*1000/Permeability<<", "<<Per_l[2]*reso*reso*1000/Permeability<<endl;
			fin<<"The relative permeability of component 2 is "<<Per_g[0]*reso*reso*1000/Permeability<<", "<<Per_g[1]*reso*reso*1000/Permeability<<", "<<Per_g[2]*reso*reso*1000/Permeability<<endl;
			fin<<"The LOCAL relative permeability of component 1 is "<<Per_l_LOCAL[0]*reso*reso*1000/Permeability<<", "<<Per_l_LOCAL[1]*reso*reso*1000/Permeability<<", "<<Per_l_LOCAL[2]*reso*reso*1000/Permeability<<endl;
			fin<<"The LOCAL relative permeability of component 2 is "<<Per_g_LOCAL[0]*reso*reso*1000/Permeability<<", "<<Per_g_LOCAL[1]*reso*reso*1000/Permeability<<", "<<Per_g_LOCAL[2]*reso*reso*1000/Permeability<<endl;
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
				finfs<<Per_l[0]*reso*reso*1000/Permeability<<" "<<S_l<<" "<<Per_l_LOCAL[0]*reso*reso*1000/Permeability<<endl;break;
				case 2:
				finfs<<Per_l[1]*reso*reso*1000/Permeability<<" "<<S_l<<" "<<Per_l_LOCAL[1]*reso*reso*1000/Permeability<<endl;break;
				case 3:
				finfs<<Per_l[2]*reso*reso*1000/Permeability<<" "<<S_l<<" "<<Per_l_LOCAL[2]*reso*reso*1000/Permeability<<endl;break;
				default:
				finfs<<Per_l[0]*reso*reso*1000/Permeability<<" "<<S_l<<" "<<Per_l_LOCAL[0]*reso*reso*1000/Permeability<<endl;break;
				}
			finfs.close();

		ofstream finfs2(FileName4,ios::app);
			switch(PerDir)
				{
				case 1:
				finfs2<<Per_g[0]*reso*reso*1000/Permeability<<" "<<1-S_l<<" "<<Per_g_LOCAL[0]*reso*reso*1000/Permeability<<endl;break;
				case 2:
				finfs2<<Per_g[1]*reso*reso*1000/Permeability<<" "<<1-S_l<<" "<<Per_g_LOCAL[1]*reso*reso*1000/Permeability<<endl;break;
				case 3:
				finfs2<<Per_g[2]*reso*reso*1000/Permeability<<" "<<1-S_l<<" "<<Per_g_LOCAL[2]*reso*reso*1000/Permeability<<endl;break;
				default:
				finfs2<<Per_g[0]*reso*reso*1000/Permeability<<" "<<1-S_l<<" "<<Per_g_LOCAL[0]*reso*reso*1000/Permeability<<endl;break;
				}
			finfs2.close();
			}

			//==========for bodyforce output===========
			ofstream finf3(FileName3,ios::app);
			finf3<<S_l<<" "<<1-S_l<<endl;
			finf3.close();
			ofstream finf5(FileName5,ios::app);
			finf5<<u_ave<<" "<<u_max<<"  "<<error<<endl;
			finf5.close();

			//cout<<"The"<<n<<"th computation result:"<<endl;
		//=============================================================================================
			//cout<<"The permiability is: "<<Permia[0]*reso*reso*1000<<", "<<Permia[1]*reso*reso*1000<<", "<<Permia[2]*reso*reso*1000<<endl;
			//cout<<"The relative error of permiability computing is: "<<error_perm<<endl;
		//==============================================================================================

		//==============================================================================================
			cout<<"The"<<n<<"th computation result:"<<endl;
			cout<<"The Density of point(NX/2,NY/2,NZ/2) is: "<<setprecision(6)
				<<rho[(int)(Count/2)]<<endl;
			
			cout<<"The Maximum velocity is: "<<setprecision(6)<<u_max<<"   Re_l="<<Re_l<<"   Re_g="<<Re_g<<endl;
			cout<<"Courant Number="<<u_max*dt/dx<<"	 Capillary Num="<<Capillary<<endl;
			cout<<"The max relative error of uv is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			cout<<"The relative permeability of component 1 is "<<Per_l[0]*reso*reso*1000/Permeability<<", "<<Per_l[1]*reso*reso*1000/Permeability<<", "<<Per_l[2]*reso*reso*1000/Permeability<<endl;
			cout<<"The relative permeability of component 2 is "<<Per_g[0]*reso*reso*1000/Permeability<<", "<<Per_g[1]*reso*reso*1000/Permeability<<", "<<Per_g[2]*reso*reso*1000/Permeability<<endl;
			cout<<"The LOCAL relative permeability of component 1 is "<<Per_l_LOCAL[0]*reso*reso*1000/Permeability<<", "<<Per_l_LOCAL[1]*reso*reso*1000/Permeability<<", "<<Per_l_LOCAL[2]*reso*reso*1000/Permeability<<endl;
			cout<<"The LOCAL relative permeability of component 2 is "<<Per_g_LOCAL[0]*reso*reso*1000/Permeability<<", "<<Per_g_LOCAL[1]*reso*reso*1000/Permeability<<", "<<Per_g_LOCAL[2]*reso*reso*1000/Permeability<<endl;
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
					output_psi(n,psi,mirX,mirY,mirZ,mir,Solid);
				else
					output_psi_b(n,psi,mirX,mirY,mirZ,mir,Solid);
			//===================================
			
			if ((fre_backup>0) and (n%fre_backup==0)  and (n>0))
			        Backup(n,rho,psi,u,f,rho_r,rho_b);
			
			
			if(error!=error) {cout<<"PROGRAM STOP"<<endl;break;};
			if(U_max_ref>=5) {cout<<"PROGRAM STOP DUE TO HIGH VELOCITY"<<endl;break;}
		}	
	}

	
	
	if (fre_backup>=0)
			        Backup(n_max,rho,psi,u,f,rho_r,rho_b);

	MPI_Barrier(MPI_COMM_WORLD);

	/*
	for (int i=0;i<=Count;i++)
		{
		delete [] u[i];
		delete [] u0[i];
		delete [] f[i];
		delete [] F[i];
		}
	
	delete [] f;
	delete [] psi;
	delete [] u;
	delete [] F;
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

*/



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

void pressure_bodyforce_change()
{



if (pressure_change==1)
		{
		if (n==pre_chan_1)
			if (pre_chan_pb==1)
			{
			p_xn=pre_chan_pn1;
			p_xp=pre_chan_pp1;
			}
			else
				gxs=pre_chan_f1;

		if (n==pre_chan_2)
			if (pre_chan_pb==1)
			{
			p_xn=pre_chan_pn2;
			p_xp=pre_chan_pp2;
			}
			else
				gxs=pre_chan_f2;

		if (n==pre_chan_3)
			if (pre_chan_pb==1)
			{
			p_xn=pre_chan_pn3;
			p_xp=pre_chan_pp3;
			}
			else
				gxs=pre_chan_f3;

		if (n==pre_chan_4)
			if (pre_chan_pb==1)
			{
			p_xn=pre_chan_pn4;
			p_xp=pre_chan_pp4;
			}
			else
				gxs=pre_chan_f4;

		if (n==pre_chan_5)
			if (pre_chan_pb==1)
			{
			p_xn=pre_chan_pn5;
			p_xp=pre_chan_pp5;
			}
			else
				gxs=pre_chan_f5;


	 
		}

if (pressure_change==2)
		{
		if (n==pre_chan_1)
			if (pre_chan_pb==1)
			{
			p_yn=pre_chan_pn1;
			p_yp=pre_chan_pp1;
			}
			else
				gys=pre_chan_f1;

		if (n==pre_chan_2)
			if (pre_chan_pb==1)
			{
			p_yn=pre_chan_pn2;
			p_yp=pre_chan_pp2;
			}
			else
				gys=pre_chan_f2;

		if (n==pre_chan_3)
			if (pre_chan_pb==1)
			{
			p_yn=pre_chan_pn3;
			p_yp=pre_chan_pp3;
			}
			else
				gys=pre_chan_f3;

		if (n==pre_chan_4)
			if (pre_chan_pb==1)
			{
			p_yn=pre_chan_pn4;
			p_yp=pre_chan_pp4;
			}
			else
				gys=pre_chan_f4;

		if (n==pre_chan_5)
			if (pre_chan_pb==1)
			{
			p_yn=pre_chan_pn5;
			p_yp=pre_chan_pp5;
			}
			else
				gys=pre_chan_f5;


	 
		}


if (pressure_change==3)
		{
		if (n==pre_chan_1)
			if (pre_chan_pb==1)
			{
			p_zn=pre_chan_pn1;
			p_zp=pre_chan_pp1;
			}
			else
				gzs=pre_chan_f1;

		if (n==pre_chan_2)
			if (pre_chan_pb==1)
			{
			p_zn=pre_chan_pn2;
			p_zp=pre_chan_pp2;
			}
			else
				gzs=pre_chan_f2;

		if (n==pre_chan_3)
			if (pre_chan_pb==1)
			{
			p_zn=pre_chan_pn3;
			p_zp=pre_chan_pp3;
			}
			else
				gzs=pre_chan_f3;

		if (n==pre_chan_4)
			if (pre_chan_pb==1)
			{
			p_zn=pre_chan_pn4;
			p_zp=pre_chan_pp4;
			}
			else
				gzs=pre_chan_f4;

		if (n==pre_chan_5)
			if (pre_chan_pb==1)
			{
			p_zn=pre_chan_pn5;
			p_zp=pre_chan_pp5;
			}
			else
				gzs=pre_chan_f5;


	 
		}





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
	MPI_Status status;
	MPI_Request request;

int nx_g[mpi_size];
int disp[mpi_size];
int pore;
float pore2;
int loc_por[NX+1];
int sum=0;
int sum2=0;
double ave_nx;
int nx_pre,nx_aft,n_i,sum_nx;

//====================
int* Solid_rank0;
float* Psi_rank0;
int* recv_solid;
float* recv_psi;

int bufsize[mpi_size];
int bufloc[mpi_size];
//=====================


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
	fstream fin;
	if (decbin==0)
	{
	FILE *ftest;
	//ifstream fin;
	
	ftest = fopen(filename, "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file (" << filename <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(-1);
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
	else
	{
	//cout<<"aaaaaaaaaaaafffffffffffsssssssssssss"<<endl;
	fin.open(filename,ios::in);
	if (fin.fail())
	        {
	        cout<<"\n file open error on " << filename<<endl;
	        exit(-1);
	        }
	Solid_rank0 = new int[(NX+1)*(NY+1)*(NZ+1)];
	fin.read((char *)(&Solid_rank0[0]), sizeof(int)*(NX+1)*(NY+1)*(NZ+1));
	for(int k=0 ; k<=NZ ; k++)
	for(int j=0 ; j<=NY ; j++)
	for(int i=0 ; i<=NX ; i++)
	
	{
		
		pore=Solid_rank0[i*(NY+1)*(NZ+1)+j*(NZ+1)+k];
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
	



}

        MPI_Bcast(loc_por,NX+1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&sum,1,MPI_INT,0,MPI_COMM_WORLD);
        
   
         
	nx_pre=0;nx_aft=0;sum_nx=0;
	ave_nx=(double)sum/(mpi_size);
	porosity=(double)sum2/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1));
	disp[0]=0;
	bufloc[0]=0;   //===========
	
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
	  
	        //======================
	        bufsize[i]=nx_g[i]*(NY+1)*(NZ+1);
	        bufloc[i+1]=bufloc[i]+bufsize[i];
	        //=======================
	        
	        }
	  
	  nx_g[mpi_size-1]=(NX+1)-disp[mpi_size-1];
	  bufsize[mpi_size-1]=nx_g[mpi_size-1]*(NY+1)*(NZ+1);//======
	  
	  nx_l=nx_g[rank];
	  
	

//===========================================	  
	    Solid = new int**[nx_l];
	Psi_local = new float**[nx_l];
	for (int i=0;i<nx_l;i++)
		{
		Solid[i] = new int*[NY+1];
		Psi_local[i] = new float*[NY+1];
			for (int j=0;j<=NY;j++)
			{
			Solid[i][j]= new int[NZ+1];
			Psi_local[i][j] = new float[NZ+1];
			}
		}
	  
	  
		recv_solid = new int[nx_l*(NY+1)*(NZ+1)];
		recv_psi = new float[nx_l*(NY+1)*(NZ+1)];
//========================================
		

MPI_Barrier(MPI_COMM_WORLD);

        
        MPI_Scatterv(Solid_rank0,bufsize,bufloc,MPI_INT,recv_solid,nx_l*(NY+1)*(NZ+1),MPI_INT,0,MPI_COMM_WORLD);

	cout<<"GEOMETRY INPUT FILE PARTITIONING FOR PARALLEL READING DONE   Processor No."<<rank<<endl;
	//cout<<endl;




    
 
 
 if (rank==0)
{       
        FILE *ftest;
	ftest = fopen(filenamepsi, "r");
	ifstream fin;
	if(ftest == NULL)
	{
		cout << "\n The Concentration file (" << filenamepsi <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest);

	Psi_rank0 = new float[(NX+1)*(NY+1)*(NZ+1)];
	
	fin.open(filenamepsi);
	for(int k=0 ; k<=NZ ; k++)
	for(int j=0 ; j<=NY ; j++)
	for(int i=0 ; i<=NX ; i++)
	
	{
	
			fin >> pore2;
			Psi_rank0[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=pore2;
	
	
	}
	fin.close();
	
}


	
MPI_Barrier(MPI_COMM_WORLD);	
	
        MPI_Scatterv(Psi_rank0,bufsize,bufloc,MPI_FLOAT,recv_psi,nx_l*(NY+1)*(NZ+1),MPI_FLOAT,0,MPI_COMM_WORLD);
	        
	        
	 cout<<"CONCENTRATION FILE PARTITIONING FOR PARALLEL READING DONE   Processor No."<<rank<<endl;
//	cout<<endl;
   
MPI_Barrier(MPI_COMM_WORLD);		        
	        
	        
	if (rank==0)
{	
//        for (int i=0;i<nx_l;i++)
//	                for (int j=0;j<=NY;j++)
//	                for (int k=0;k<=NZ;k++)
//	                {
//	                        recv_solid[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=Solid_rank0[i*(NY+1)*(NZ+1)+j*(NZ+1)+k];
//	                        recv_psi[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=Psi_rank0[i*(NY+1)*(NZ+1)+j*(NZ+1)+k];
//	                }
	       delete [] Solid_rank0;	
	       delete [] Psi_rank0;  
	
	       
	
}



	for (int i=0;i<nx_l;i++)
	                for (int j=0;j<=NY;j++)
	                for (int k=0;k<=NZ;k++)
	                {
	                Solid[i][j][k]=recv_solid[i*(NY+1)*(NZ+1)+j*(NZ+1)+k];
	                Psi_local[i][j][k]=recv_psi[i*(NY+1)*(NZ+1)+j*(NZ+1)+k];
	                
	                }
	           
	 delete [] recv_psi;
	 delete [] recv_solid;
       
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

void init_Sparse_read_rock_parallel(int* Sl,int* Sr)
{	
MPI_Status status[4] ;
MPI_Request request[4];

	
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

bool mark;
int kk,ip,jp,kp,mean_l,mpi_test,s_c;

int pore;
double pore2;

	

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
			{
				Solid[i][j][k]=0;
			}

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




void init(double* rho, double** u, double** f,double* psi,double* rho_r, double* rho_b, double* rhor, double* rhob, int* SupInv)
{	
       int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	
	int nx_g[mpi_size];
	int disp[mpi_size];
	int si,sj,sm;
	
	MPI_Gather(&nx_l,1,MPI_INT,nx_g,1,MPI_INT,0,MPI_COMM_WORLD);
	
	
	if (rank==0)
		{
		disp[0]=0;
	
		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
		
		}

	MPI_Bcast(disp,mpi_size,MPI_INT,0,MPI_COMM_WORLD);

        srand((unsigned)time(0)+rank);
	
	double usqr,vsqr,rand_double;
	double c2,c4;
	
	rho0=1.0;dt=1.0/Zoom;dx=1.0/Zoom;
 
	if (lattice_v==1)
		{dx=dx_input;dt=dt_input;}

	lat_c=dx/dt;
	c_s=lat_c/sqrt(3);
	c_s2=lat_c*lat_c/3;

	c2=lat_c*lat_c;c4=c2*c2;
	
	
	niu=in_vis;
	
	//tau_f=3.0*niu/dt+0.5;
	tau_f=niu/(c_s2*dt)+0.5;
	s_v=1/tau_f;
        
	//double pr; //raduis of the obstacles
       int loc_x,loc_y,loc_z;
       
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

	if (ini_Sat<0)
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
					f[i][lm]=feq(lm,rho[i],u_tmp);	

	}
	else
	for (int i=1;i<=Count;i++)	
			
		{
			u[i][0]=inivx;
			u[i][1]=inivy;
			u[i][2]=inivz;
			u_tmp[0]=u[i][0];
			u_tmp[1]=u[i][1];
			u_tmp[2]=u[i][2];
			rand_double=(double(rand()%10000))/10000;
			if (rand_double<ini_Sat)
			        psi[i]=1;
			else
			        psi[i]=-1;
			
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
				//if (Solid[(int)(SupInv[i]/((NY+1)*(NZ+1)))][(int)((SupInv[i]%((NY+1)*(NZ+1)))/(NZ+1))][SupInv[i]%(NZ+1)]<0)
				//	f[i][lm]=0.0;
				//else
					f[i][lm]=feq(lm,rho[i],u_tmp);	
			
			if ((par_per_x>0) or (par_per_y>0) or (par_per_z>0))
			        {        
			                loc_x=(int)(SupInv[i]/((NY+1)*(NZ+1)))+disp[rank];
			                loc_y=(int)((SupInv[i]%((NY+1)*(NZ+1)))/(NZ+1));
			                loc_z=SupInv[i]%(NZ+1);
			                if (((loc_x<per_xn) or (loc_x>per_xp) or (loc_y<per_yn) or (loc_y>per_yp) or (loc_z<per_zn) or (loc_z>per_zp)) and ((ini_buf==-1) or (ini_buf==1)))
			                {        
			                        psi[i]=ini_buf;
			                        rho_r[i]=(psi[i]*rho[i]+rho[i])/2;
			                        rho_b[i]=rho[i]-rho_r[i];
			                       
			                }
			                
			                
			        }
					
	}       
	

	

	 	
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
	MPI_Status status2[12] ;
	MPI_Request request2[12];
	int mpi_test;
	int mpi_test2;
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();


double C[3];
double g_r[19],g_b[19];
double rho_0=1.0;
double lm0,lm1,cc,sum,uu;
double ux,uy,uz,nx,ny,nz;
double usqr,vsqr,eu,ef,cospsi,s_other;
double F_hat[19],GuoF[19],f_eq[19],u_tmp[3];
double m_l[19],m_inv_l[19];
int i,j,m,ind_S;
int interi,interj,interk,ip,jp,kp;
double c2,c4;
double delta_rho=0.1;
const double c_l=lat_c;
int mi;

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

	double* recvl= new double[Gcl[rank]*5];
	double* recvr = new double[Gcr[rank]*5];


	double* sendl_psi = new double[Gcl[rank]];
	double* sendr_psi = new double[Gcr[rank]];
	double* recvl_rhor= new double[Gcl[rank]];
	double* recvr_rhor = new double[Gcr[rank]];
	double* recvl_rhob= new double[Gcl[rank]];
	double* recvr_rhob = new double[Gcr[rank]];

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


		for (int ka=0;ka<Gcr[mpi_size-1];ka++)
		        {
		                sendl_rhor[ka]=0;sendl_rhob[ka]=0;
				for (int kb=0;kb<5;kb++)
				sendl[ka*5+kb]=0;
		                
		        }
		 for (int ka=0;ka<Gcl[rank+1];ka++)
		        {
		                sendr_rhor[ka]=0;sendr_rhob[ka]=0;
		                for (int kb=0;kb<5;kb++)
				sendr[ka*5+kb]=0;
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


			for (int ka=0;ka<Gcr[rank-1];ka++)
		        {
		                sendl_rhor[ka]=0;sendl_rhob[ka]=0;
				for (int kb=0;kb<5;kb++)
				sendl[ka*5+kb]=0;
		                
		        }
		        for (int ka=0;ka<Gcl[0];ka++)
		        {
		                sendr_rhor[ka]=0;sendr_rhob[ka]=0;
				for (int kb=0;kb<5;kb++)
				sendr[ka*5+kb]=0;
		                
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

			for (int ka=0;ka<Gcr[rank-1];ka++)
		        {
		                sendl_rhor[ka]=0;sendl_rhob[ka]=0;
				for (int kb=0;kb<5;kb++)
				sendl[ka*5+kb]=0;
		                
		        }
		        for (int ka=0;ka<Gcl[rank+1];ka++)
		        {
		                sendr_rhor[ka]=0;sendr_rhob[ka]=0;
				for (int kb=0;kb<5;kb++)
				sendr[ka*5+kb]=0;
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
		if (in_psi_BC>0)
		{	     
			interi=i+e[tmpi][0];
			if ((psi_xn>0) and (rank==0) and (interi<0))
			        interi=0;
			if ((psi_xp>0) and (rank==mpi_size-1) and (interi>=nx_l))
			        interi=nx_l-1;
			
			interj=j+e[tmpi][1];
			if (psi_yn-1>0)
			        {if (interj<0) {interj=0;}}
			else
			        {if (interj<0) {interj=NY;}}
			 
			if (psi_yp>0)
			        {if (interj>NY) {interj=NY;}}
			else
			        {if (interj>NY) {interj=0;}}
			
			
			interk=m+e[tmpi][2];
			if (psi_zn>0)
			        {if (interk<0) {interk=0;}}
			else
			        {if (interk<0) {interk=NZ;}}
			
			if (psi_zp>0)
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

		if ((sqrt((rho_r[ci]-rho_b[ci])*(rho_r[ci]-rho_b[ci]))>=0.9) and (ind_S==1))
		{C[0]=0;C[1]=0;C[2]=0;}

		//C[0]=0;C[1]=0;C[2]=0;
			cc=sqrt(C[0]*C[0]+C[1]*C[1]+C[2]*C[2]);
			if (cc>0)
			        {nx=C[0]/cc;ny=C[1]/cc;nz=C[2]/cc;}
			else
			        {nx=0;ny=0;nz=0;}


			//=================FORCE TERM_GUO=========================================

/*
			for (int k=0;k<19;k++)
			{	
			lm0=((elat[k][0]-u[ci][0])*gxs+(elat[k][1]-u[ci][1])*gys+(elat[k][2]-u[ci][2])*gzs)/c_s2;
			lm1=(elat[k][0]*u[ci][0]+elat[k][1]*u[ci][1]+elat[k][2]*u[ci][2])*(elat[k][0]*gxs+elat[k][1]*gys+elat[k][2]*gzs)/(c_s2*c_s2);
			GuoF[k]=w[k]*(lm0+lm1);
			//GuoF[k]=0.0;
			}

*/


//==========================
lm0=((+0.000*c_l-u[ci][0])*gxs+(+0.000*c_l-u[ci][1])*gys+(+0.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(+0.000*c_l*gxs+(+0.000*c_l)*gys+(+0.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[0]=w[0]*(lm0+lm1);

lm0=((+1.000*c_l-u[ci][0])*gxs+(+0.000*c_l-u[ci][1])*gys+(+0.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(+1.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(+1.000*c_l*gxs+(+0.000*c_l)*gys+(+0.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[1]=w[1]*(lm0+lm1);

lm0=((-1.000*c_l-u[ci][0])*gxs+(+0.000*c_l-u[ci][1])*gys+(+0.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(-1.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(-1.000*c_l*gxs+(+0.000*c_l)*gys+(+0.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[2]=w[2]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gxs+(+1.000*c_l-u[ci][1])*gys+(+0.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(+1.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(+0.000*c_l*gxs+(+1.000*c_l)*gys+(+0.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[3]=w[3]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gxs+(-1.000*c_l-u[ci][1])*gys+(+0.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(-1.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(+0.000*c_l*gxs+(-1.000*c_l)*gys+(+0.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[4]=w[4]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gxs+(+0.000*c_l-u[ci][1])*gys+(+1.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(+1.000*c_l)*u[ci][2])*(+0.000*c_l*gxs+(+0.000*c_l)*gys+(+1.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[5]=w[5]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gxs+(+0.000*c_l-u[ci][1])*gys+(-1.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(-1.000*c_l)*u[ci][2])*(+0.000*c_l*gxs+(+0.000*c_l)*gys+(-1.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[6]=w[6]*(lm0+lm1);

lm0=((+1.000*c_l-u[ci][0])*gxs+(+1.000*c_l-u[ci][1])*gys+(+0.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(+1.000*c_l*u[ci][0]+(+1.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(+1.000*c_l*gxs+(+1.000*c_l)*gys+(+0.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[7]=w[7]*(lm0+lm1);

lm0=((-1.000*c_l-u[ci][0])*gxs+(-1.000*c_l-u[ci][1])*gys+(+0.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(-1.000*c_l*u[ci][0]+(-1.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(-1.000*c_l*gxs+(-1.000*c_l)*gys+(+0.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[8]=w[8]*(lm0+lm1);

lm0=((+1.000*c_l-u[ci][0])*gxs+(-1.000*c_l-u[ci][1])*gys+(+0.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(+1.000*c_l*u[ci][0]+(-1.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(+1.000*c_l*gxs+(-1.000*c_l)*gys+(+0.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[9]=w[9]*(lm0+lm1);

lm0=((-1.000*c_l-u[ci][0])*gxs+(+1.000*c_l-u[ci][1])*gys+(+0.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(-1.000*c_l*u[ci][0]+(+1.000*c_l)*u[ci][1]+(+0.000*c_l)*u[ci][2])*(-1.000*c_l*gxs+(+1.000*c_l)*gys+(+0.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[10]=w[10]*(lm0+lm1);

lm0=((+1.000*c_l-u[ci][0])*gxs+(+0.000*c_l-u[ci][1])*gys+(+1.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(+1.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(+1.000*c_l)*u[ci][2])*(+1.000*c_l*gxs+(+0.000*c_l)*gys+(+1.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[11]=w[11]*(lm0+lm1);

lm0=((-1.000*c_l-u[ci][0])*gxs+(+0.000*c_l-u[ci][1])*gys+(-1.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(-1.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(-1.000*c_l)*u[ci][2])*(-1.000*c_l*gxs+(+0.000*c_l)*gys+(-1.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[12]=w[12]*(lm0+lm1);

lm0=((+1.000*c_l-u[ci][0])*gxs+(+0.000*c_l-u[ci][1])*gys+(-1.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(+1.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(-1.000*c_l)*u[ci][2])*(+1.000*c_l*gxs+(+0.000*c_l)*gys+(-1.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[13]=w[13]*(lm0+lm1);

lm0=((-1.000*c_l-u[ci][0])*gxs+(+0.000*c_l-u[ci][1])*gys+(+1.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(-1.000*c_l*u[ci][0]+(+0.000*c_l)*u[ci][1]+(+1.000*c_l)*u[ci][2])*(-1.000*c_l*gxs+(+0.000*c_l)*gys+(+1.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[14]=w[14]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gxs+(+1.000*c_l-u[ci][1])*gys+(+1.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(+1.000*c_l)*u[ci][1]+(+1.000*c_l)*u[ci][2])*(+0.000*c_l*gxs+(+1.000*c_l)*gys+(+1.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[15]=w[15]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gxs+(-1.000*c_l-u[ci][1])*gys+(-1.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(-1.000*c_l)*u[ci][1]+(-1.000*c_l)*u[ci][2])*(+0.000*c_l*gxs+(-1.000*c_l)*gys+(-1.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[16]=w[16]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gxs+(+1.000*c_l-u[ci][1])*gys+(-1.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(+1.000*c_l)*u[ci][1]+(-1.000*c_l)*u[ci][2])*(+0.000*c_l*gxs+(+1.000*c_l)*gys+(-1.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[17]=w[17]*(lm0+lm1);

lm0=((+0.000*c_l-u[ci][0])*gxs+(-1.000*c_l-u[ci][1])*gys+(+1.000*c_l-u[ci][2])*gzs)/c_s2;
lm1=(+0.000*c_l*u[ci][0]+(-1.000*c_l)*u[ci][1]+(+1.000*c_l)*u[ci][2])*(+0.000*c_l*gxs+(-1.000*c_l)*gys+(+1.000*c_l)*gzs)/(c_s2*c_s2);
GuoF[18]=w[18]*(lm0+lm1);

//====================

			
			//=====================equilibrium of moment=================================
			ux=u[ci][0];
			uy=u[ci][1];
			uz=u[ci][2];
			
			for(int k=0;k<19;k++)
				meq[k]=0;

			
	//========================================================================================		
			meq[0]=rho[ci];meq[3]=rho_0*ux;meq[5]=rho_0*uy;meq[7]=rho_0*uz;
			meq[1]=rho_0*(ux*ux+uy*uy+uz*uz)+CapA*cc;
			meq[9]=rho_0*(2*ux*ux-uy*uy-uz*uz)+0.5*CapA*cc*(2*nx*nx-ny*ny-nz*nz);
			meq[11]=rho_0*(uy*uy-uz*uz)+0.5*CapA*cc*(ny*ny-nz*nz);
			meq[13]=rho_0*ux*uy+0.5*CapA*cc*(nx*ny);
			meq[14]=rho_0*uy*uz+0.5*CapA*cc*(ny*nz);
			meq[15]=rho_0*ux*uz+0.5*CapA*cc*(nx*nz);
	//========================================================================================
	//		meq[0]=var_rho;meq[3]=rho_0*ux;meq[5]=rho_0*uy;meq[7]=rho_0*uz;
	//		meq[1]=-CapA*cc;
	//		meq[9]=0.5*CapA*cc*(2*nx*nx-ny*ny-nz*nz);
	//		meq[11]=0.5*CapA*cc*(ny*ny-nz*nz);
	//		meq[13]=0.5*CapA*cc*(nx*ny);
	//		meq[14]=0.5*CapA*cc*(ny*nz);
	//		meq[15]=0.5*CapA*cc*(nx*nz);

	//=======================================================================================		
			
			s_v=niu_g+(psi[ci]+1.0)/2.0*(niu_l-niu_g);
			s_v=1.0/(s_v/(c_s2*dt)+0.5);
			s_other=8*(2-s_v)/(8-s_v);
			
	//cout<<"@@@@@@@@@   "<<s_v<<"  "<<C[0]<<"   "<<C[1]<<"  "<<C[2]<<endl;
	//cout<<"@@@@@@@@@   "<<s_v<<"  "<<ux<<"   "<<uy<<"  "<<uz<<"  "<<rho_r[ci]<<" "<<rho_b[ci]<<endl;
	
	
	S[1]=s_v;S[2]=s_v;S[4]=s_other;S[6]=s_other;S[8]=s_other;S[9]=s_v;
	S[10]=s_v;S[11]=s_v;S[12]=s_v;S[13]=s_v;S[14]=s_v;S[15]=s_v;S[16]=s_other;
	S[17]=s_other;S[18]=s_other;



			//============================================================================

			/*
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
					}
			//============================================================================
			*/


//==========================
m_l[0]=+1.000*1.0*f[ci][0]+1.00000000000000000*1.0*f[ci][1]+1.00000000000000000*1.0*f[ci][2]+1.00000000000000000*1.0*f[ci][3]+1.00000000000000000*1.0*f[ci][4]+1.00000000000000000*1.0*f[ci][5]+1.00000000000000000*1.0*f[ci][6]+1.00000000000000000*1.0*f[ci][7]+1.00000000000000000*1.0*f[ci][8]+1.00000000000000000*1.0*f[ci][9]+1.00000000000000000*1.0*f[ci][10]+1.00000000000000000*1.0*f[ci][11]+1.00000000000000000*1.0*f[ci][12]+1.00000000000000000*1.0*f[ci][13]+1.00000000000000000*1.0*f[ci][14]+1.00000000000000000*1.0*f[ci][15]+1.00000000000000000*1.0*f[ci][16]+1.00000000000000000*1.0*f[ci][17]+1.00000000000000000*1.0*f[ci][18];

F_hat[0]=+1.000*1.0*GuoF[0]+1.00000000000000000*1.0*GuoF[1]+1.00000000000000000*1.0*GuoF[2]+1.00000000000000000*1.0*GuoF[3]+1.00000000000000000*1.0*GuoF[4]+1.00000000000000000*1.0*GuoF[5]+1.00000000000000000*1.0*GuoF[6]+1.00000000000000000*1.0*GuoF[7]+1.00000000000000000*1.0*GuoF[8]+1.00000000000000000*1.0*GuoF[9]+1.00000000000000000*1.0*GuoF[10]+1.00000000000000000*1.0*GuoF[11]+1.00000000000000000*1.0*GuoF[12]+1.00000000000000000*1.0*GuoF[13]+1.00000000000000000*1.0*GuoF[14]+1.00000000000000000*1.0*GuoF[15]+1.00000000000000000*1.0*GuoF[16]+1.00000000000000000*1.0*GuoF[17]+1.00000000000000000*1.0*GuoF[18];

F_hat[0]*=(1-0.5*S[0]);
m_l[0]=m_l[0]-S[0]*(m_l[0]-meq[0])+dt*F_hat[0];
//=======================================

m_l[1]=-1.00000000000000000*c_l*c_l*f[ci][0]+0.00000000000000000*c_l*c_l*f[ci][1]+0.00000000000000000*c_l*c_l*f[ci][2]+0.00000000000000000*c_l*c_l*f[ci][3]+0.00000000000000000*c_l*c_l*f[ci][4]+0.00000000000000000*c_l*c_l*f[ci][5]+0.00000000000000000*c_l*c_l*f[ci][6]+1.00000000000000000*c_l*c_l*f[ci][7]+1.00000000000000000*c_l*c_l*f[ci][8]+1.00000000000000000*c_l*c_l*f[ci][9]+1.00000000000000000*c_l*c_l*f[ci][10]+1.00000000000000000*c_l*c_l*f[ci][11]+1.00000000000000000*c_l*c_l*f[ci][12]+1.00000000000000000*c_l*c_l*f[ci][13]+1.00000000000000000*c_l*c_l*f[ci][14]+1.00000000000000000*c_l*c_l*f[ci][15]+1.00000000000000000*c_l*c_l*f[ci][16]+1.00000000000000000*c_l*c_l*f[ci][17]+1.00000000000000000*c_l*c_l*f[ci][18];

F_hat[1]=-1.00000000000000000*c_l*c_l*GuoF[0]+0.00000000000000000*c_l*c_l*GuoF[1]+0.00000000000000000*c_l*c_l*GuoF[2]+0.00000000000000000*c_l*c_l*GuoF[3]+0.00000000000000000*c_l*c_l*GuoF[4]+0.00000000000000000*c_l*c_l*GuoF[5]+0.00000000000000000*c_l*c_l*GuoF[6]+1.00000000000000000*c_l*c_l*GuoF[7]+1.00000000000000000*c_l*c_l*GuoF[8]+1.00000000000000000*c_l*c_l*GuoF[9]+1.00000000000000000*c_l*c_l*GuoF[10]+1.00000000000000000*c_l*c_l*GuoF[11]+1.00000000000000000*c_l*c_l*GuoF[12]+1.00000000000000000*c_l*c_l*GuoF[13]+1.00000000000000000*c_l*c_l*GuoF[14]+1.00000000000000000*c_l*c_l*GuoF[15]+1.00000000000000000*c_l*c_l*GuoF[16]+1.00000000000000000*c_l*c_l*GuoF[17]+1.00000000000000000*c_l*c_l*GuoF[18];

F_hat[1]*=(1-0.5*S[1]);
m_l[1]=m_l[1]-S[1]*(m_l[1]-meq[1])+dt*F_hat[1];
//=======================================

m_l[2]=+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][0]-2.00000000000000000*c_l*c_l*c_l*c_l*f[ci][1]-2.00000000000000000*c_l*c_l*c_l*c_l*f[ci][2]-2.00000000000000000*c_l*c_l*c_l*c_l*f[ci][3]-2.00000000000000000*c_l*c_l*c_l*c_l*f[ci][4]-2.00000000000000000*c_l*c_l*c_l*c_l*f[ci][5]-2.00000000000000000*c_l*c_l*c_l*c_l*f[ci][6]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][7]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][8]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][9]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][10]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][11]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][12]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][13]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][14]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][15]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][16]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][17]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][18];

F_hat[2]=+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[0]-2.00000000000000000*c_l*c_l*c_l*c_l*GuoF[1]-2.00000000000000000*c_l*c_l*c_l*c_l*GuoF[2]-2.00000000000000000*c_l*c_l*c_l*c_l*GuoF[3]-2.00000000000000000*c_l*c_l*c_l*c_l*GuoF[4]-2.00000000000000000*c_l*c_l*c_l*c_l*GuoF[5]-2.00000000000000000*c_l*c_l*c_l*c_l*GuoF[6]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[7]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[8]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[9]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[10]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[11]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[12]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[13]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[14]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[15]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[16]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[17]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[18];

F_hat[2]*=(1-0.5*S[2]);
m_l[2]=m_l[2]-S[2]*(m_l[2]-meq[2])+dt*F_hat[2];
//=======================================

m_l[3]=+0.00000000000000000*c_l*f[ci][0]+1.00000000000000000*c_l*f[ci][1]-1.00000000000000000*c_l*f[ci][2]+0.00000000000000000*c_l*f[ci][3]+0.00000000000000000*c_l*f[ci][4]+0.00000000000000000*c_l*f[ci][5]+0.00000000000000000*c_l*f[ci][6]+1.00000000000000000*c_l*f[ci][7]-1.00000000000000000*c_l*f[ci][8]+1.00000000000000000*c_l*f[ci][9]-1.00000000000000000*c_l*f[ci][10]+1.00000000000000000*c_l*f[ci][11]-1.00000000000000000*c_l*f[ci][12]+1.00000000000000000*c_l*f[ci][13]-1.00000000000000000*c_l*f[ci][14]+0.00000000000000000*c_l*f[ci][15]+0.00000000000000000*c_l*f[ci][16]+0.00000000000000000*c_l*f[ci][17]+0.00000000000000000*c_l*f[ci][18];

F_hat[3]=+0.00000000000000000*c_l*GuoF[0]+1.00000000000000000*c_l*GuoF[1]-1.00000000000000000*c_l*GuoF[2]+0.00000000000000000*c_l*GuoF[3]+0.00000000000000000*c_l*GuoF[4]+0.00000000000000000*c_l*GuoF[5]+0.00000000000000000*c_l*GuoF[6]+1.00000000000000000*c_l*GuoF[7]-1.00000000000000000*c_l*GuoF[8]+1.00000000000000000*c_l*GuoF[9]-1.00000000000000000*c_l*GuoF[10]+1.00000000000000000*c_l*GuoF[11]-1.00000000000000000*c_l*GuoF[12]+1.00000000000000000*c_l*GuoF[13]-1.00000000000000000*c_l*GuoF[14]+0.00000000000000000*c_l*GuoF[15]+0.00000000000000000*c_l*GuoF[16]+0.00000000000000000*c_l*GuoF[17]+0.00000000000000000*c_l*GuoF[18];

F_hat[3]*=(1-0.5*S[3]);
m_l[3]=m_l[3]-S[3]*(m_l[3]-meq[3])+dt*F_hat[3];
//=======================================

m_l[4]=+0.00000000000000000*c_l*c_l*c_l*f[ci][0]-2.00000000000000000*c_l*c_l*c_l*f[ci][1]+2.00000000000000000*c_l*c_l*c_l*f[ci][2]+0.00000000000000000*c_l*c_l*c_l*f[ci][3]+0.00000000000000000*c_l*c_l*c_l*f[ci][4]+0.00000000000000000*c_l*c_l*c_l*f[ci][5]+0.00000000000000000*c_l*c_l*c_l*f[ci][6]+1.00000000000000000*c_l*c_l*c_l*f[ci][7]-1.00000000000000000*c_l*c_l*c_l*f[ci][8]+1.00000000000000000*c_l*c_l*c_l*f[ci][9]-1.00000000000000000*c_l*c_l*c_l*f[ci][10]+1.00000000000000000*c_l*c_l*c_l*f[ci][11]-1.00000000000000000*c_l*c_l*c_l*f[ci][12]+1.00000000000000000*c_l*c_l*c_l*f[ci][13]-1.00000000000000000*c_l*c_l*c_l*f[ci][14]+0.00000000000000000*c_l*c_l*c_l*f[ci][15]+0.00000000000000000*c_l*c_l*c_l*f[ci][16]+0.00000000000000000*c_l*c_l*c_l*f[ci][17]+0.00000000000000000*c_l*c_l*c_l*f[ci][18];

F_hat[4]=+0.00000000000000000*c_l*c_l*c_l*GuoF[0]-2.00000000000000000*c_l*c_l*c_l*GuoF[1]+2.00000000000000000*c_l*c_l*c_l*GuoF[2]+0.00000000000000000*c_l*c_l*c_l*GuoF[3]+0.00000000000000000*c_l*c_l*c_l*GuoF[4]+0.00000000000000000*c_l*c_l*c_l*GuoF[5]+0.00000000000000000*c_l*c_l*c_l*GuoF[6]+1.00000000000000000*c_l*c_l*c_l*GuoF[7]-1.00000000000000000*c_l*c_l*c_l*GuoF[8]+1.00000000000000000*c_l*c_l*c_l*GuoF[9]-1.00000000000000000*c_l*c_l*c_l*GuoF[10]+1.00000000000000000*c_l*c_l*c_l*GuoF[11]-1.00000000000000000*c_l*c_l*c_l*GuoF[12]+1.00000000000000000*c_l*c_l*c_l*GuoF[13]-1.00000000000000000*c_l*c_l*c_l*GuoF[14]+0.00000000000000000*c_l*c_l*c_l*GuoF[15]+0.00000000000000000*c_l*c_l*c_l*GuoF[16]+0.00000000000000000*c_l*c_l*c_l*GuoF[17]+0.00000000000000000*c_l*c_l*c_l*GuoF[18];

F_hat[4]*=(1-0.5*S[4]);
m_l[4]=m_l[4]-S[4]*(m_l[4]-meq[4])+dt*F_hat[4];
//=======================================

m_l[5]=+0.00000000000000000*c_l*f[ci][0]+0.00000000000000000*c_l*f[ci][1]+0.00000000000000000*c_l*f[ci][2]+1.00000000000000000*c_l*f[ci][3]-1.00000000000000000*c_l*f[ci][4]+0.00000000000000000*c_l*f[ci][5]+0.00000000000000000*c_l*f[ci][6]+1.00000000000000000*c_l*f[ci][7]-1.00000000000000000*c_l*f[ci][8]-1.00000000000000000*c_l*f[ci][9]+1.00000000000000000*c_l*f[ci][10]+0.00000000000000000*c_l*f[ci][11]+0.00000000000000000*c_l*f[ci][12]+0.00000000000000000*c_l*f[ci][13]+0.00000000000000000*c_l*f[ci][14]+1.00000000000000000*c_l*f[ci][15]-1.00000000000000000*c_l*f[ci][16]+1.00000000000000000*c_l*f[ci][17]-1.00000000000000000*c_l*f[ci][18];

F_hat[5]=+0.00000000000000000*c_l*GuoF[0]+0.00000000000000000*c_l*GuoF[1]+0.00000000000000000*c_l*GuoF[2]+1.00000000000000000*c_l*GuoF[3]-1.00000000000000000*c_l*GuoF[4]+0.00000000000000000*c_l*GuoF[5]+0.00000000000000000*c_l*GuoF[6]+1.00000000000000000*c_l*GuoF[7]-1.00000000000000000*c_l*GuoF[8]-1.00000000000000000*c_l*GuoF[9]+1.00000000000000000*c_l*GuoF[10]+0.00000000000000000*c_l*GuoF[11]+0.00000000000000000*c_l*GuoF[12]+0.00000000000000000*c_l*GuoF[13]+0.00000000000000000*c_l*GuoF[14]+1.00000000000000000*c_l*GuoF[15]-1.00000000000000000*c_l*GuoF[16]+1.00000000000000000*c_l*GuoF[17]-1.00000000000000000*c_l*GuoF[18];

F_hat[5]*=(1-0.5*S[5]);
m_l[5]=m_l[5]-S[5]*(m_l[5]-meq[5])+dt*F_hat[5];
//=======================================

m_l[6]=+0.00000000000000000*c_l*c_l*c_l*f[ci][0]+0.00000000000000000*c_l*c_l*c_l*f[ci][1]+0.00000000000000000*c_l*c_l*c_l*f[ci][2]-2.00000000000000000*c_l*c_l*c_l*f[ci][3]+2.00000000000000000*c_l*c_l*c_l*f[ci][4]+0.00000000000000000*c_l*c_l*c_l*f[ci][5]+0.00000000000000000*c_l*c_l*c_l*f[ci][6]+1.00000000000000000*c_l*c_l*c_l*f[ci][7]-1.00000000000000000*c_l*c_l*c_l*f[ci][8]-1.00000000000000000*c_l*c_l*c_l*f[ci][9]+1.00000000000000000*c_l*c_l*c_l*f[ci][10]+0.00000000000000000*c_l*c_l*c_l*f[ci][11]+0.00000000000000000*c_l*c_l*c_l*f[ci][12]+0.00000000000000000*c_l*c_l*c_l*f[ci][13]+0.00000000000000000*c_l*c_l*c_l*f[ci][14]+1.00000000000000000*c_l*c_l*c_l*f[ci][15]-1.00000000000000000*c_l*c_l*c_l*f[ci][16]+1.00000000000000000*c_l*c_l*c_l*f[ci][17]-1.00000000000000000*c_l*c_l*c_l*f[ci][18];

F_hat[6]=+0.00000000000000000*c_l*c_l*c_l*GuoF[0]+0.00000000000000000*c_l*c_l*c_l*GuoF[1]+0.00000000000000000*c_l*c_l*c_l*GuoF[2]-2.00000000000000000*c_l*c_l*c_l*GuoF[3]+2.00000000000000000*c_l*c_l*c_l*GuoF[4]+0.00000000000000000*c_l*c_l*c_l*GuoF[5]+0.00000000000000000*c_l*c_l*c_l*GuoF[6]+1.00000000000000000*c_l*c_l*c_l*GuoF[7]-1.00000000000000000*c_l*c_l*c_l*GuoF[8]-1.00000000000000000*c_l*c_l*c_l*GuoF[9]+1.00000000000000000*c_l*c_l*c_l*GuoF[10]+0.00000000000000000*c_l*c_l*c_l*GuoF[11]+0.00000000000000000*c_l*c_l*c_l*GuoF[12]+0.00000000000000000*c_l*c_l*c_l*GuoF[13]+0.00000000000000000*c_l*c_l*c_l*GuoF[14]+1.00000000000000000*c_l*c_l*c_l*GuoF[15]-1.00000000000000000*c_l*c_l*c_l*GuoF[16]+1.00000000000000000*c_l*c_l*c_l*GuoF[17]-1.00000000000000000*c_l*c_l*c_l*GuoF[18];

F_hat[6]*=(1-0.5*S[6]);
m_l[6]=m_l[6]-S[6]*(m_l[6]-meq[6])+dt*F_hat[6];
//=======================================

m_l[7]=+0.00000000000000000*c_l*f[ci][0]+0.00000000000000000*c_l*f[ci][1]+0.00000000000000000*c_l*f[ci][2]+0.00000000000000000*c_l*f[ci][3]+0.00000000000000000*c_l*f[ci][4]+1.00000000000000000*c_l*f[ci][5]-1.00000000000000000*c_l*f[ci][6]+0.00000000000000000*c_l*f[ci][7]+0.00000000000000000*c_l*f[ci][8]+0.00000000000000000*c_l*f[ci][9]+0.00000000000000000*c_l*f[ci][10]+1.00000000000000000*c_l*f[ci][11]-1.00000000000000000*c_l*f[ci][12]-1.00000000000000000*c_l*f[ci][13]+1.00000000000000000*c_l*f[ci][14]+1.00000000000000000*c_l*f[ci][15]-1.00000000000000000*c_l*f[ci][16]-1.00000000000000000*c_l*f[ci][17]+1.00000000000000000*c_l*f[ci][18];

F_hat[7]=+0.00000000000000000*c_l*GuoF[0]+0.00000000000000000*c_l*GuoF[1]+0.00000000000000000*c_l*GuoF[2]+0.00000000000000000*c_l*GuoF[3]+0.00000000000000000*c_l*GuoF[4]+1.00000000000000000*c_l*GuoF[5]-1.00000000000000000*c_l*GuoF[6]+0.00000000000000000*c_l*GuoF[7]+0.00000000000000000*c_l*GuoF[8]+0.00000000000000000*c_l*GuoF[9]+0.00000000000000000*c_l*GuoF[10]+1.00000000000000000*c_l*GuoF[11]-1.00000000000000000*c_l*GuoF[12]-1.00000000000000000*c_l*GuoF[13]+1.00000000000000000*c_l*GuoF[14]+1.00000000000000000*c_l*GuoF[15]-1.00000000000000000*c_l*GuoF[16]-1.00000000000000000*c_l*GuoF[17]+1.00000000000000000*c_l*GuoF[18];

F_hat[7]*=(1-0.5*S[7]);
m_l[7]=m_l[7]-S[7]*(m_l[7]-meq[7])+dt*F_hat[7];
//=======================================

m_l[8]=+0.00000000000000000*c_l*c_l*c_l*f[ci][0]+0.00000000000000000*c_l*c_l*c_l*f[ci][1]+0.00000000000000000*c_l*c_l*c_l*f[ci][2]+0.00000000000000000*c_l*c_l*c_l*f[ci][3]+0.00000000000000000*c_l*c_l*c_l*f[ci][4]-2.00000000000000000*c_l*c_l*c_l*f[ci][5]+2.00000000000000000*c_l*c_l*c_l*f[ci][6]+0.00000000000000000*c_l*c_l*c_l*f[ci][7]+0.00000000000000000*c_l*c_l*c_l*f[ci][8]+0.00000000000000000*c_l*c_l*c_l*f[ci][9]+0.00000000000000000*c_l*c_l*c_l*f[ci][10]+1.00000000000000000*c_l*c_l*c_l*f[ci][11]-1.00000000000000000*c_l*c_l*c_l*f[ci][12]-1.00000000000000000*c_l*c_l*c_l*f[ci][13]+1.00000000000000000*c_l*c_l*c_l*f[ci][14]+1.00000000000000000*c_l*c_l*c_l*f[ci][15]-1.00000000000000000*c_l*c_l*c_l*f[ci][16]-1.00000000000000000*c_l*c_l*c_l*f[ci][17]+1.00000000000000000*c_l*c_l*c_l*f[ci][18];

F_hat[8]=+0.00000000000000000*c_l*c_l*c_l*GuoF[0]+0.00000000000000000*c_l*c_l*c_l*GuoF[1]+0.00000000000000000*c_l*c_l*c_l*GuoF[2]+0.00000000000000000*c_l*c_l*c_l*GuoF[3]+0.00000000000000000*c_l*c_l*c_l*GuoF[4]-2.00000000000000000*c_l*c_l*c_l*GuoF[5]+2.00000000000000000*c_l*c_l*c_l*GuoF[6]+0.00000000000000000*c_l*c_l*c_l*GuoF[7]+0.00000000000000000*c_l*c_l*c_l*GuoF[8]+0.00000000000000000*c_l*c_l*c_l*GuoF[9]+0.00000000000000000*c_l*c_l*c_l*GuoF[10]+1.00000000000000000*c_l*c_l*c_l*GuoF[11]-1.00000000000000000*c_l*c_l*c_l*GuoF[12]-1.00000000000000000*c_l*c_l*c_l*GuoF[13]+1.00000000000000000*c_l*c_l*c_l*GuoF[14]+1.00000000000000000*c_l*c_l*c_l*GuoF[15]-1.00000000000000000*c_l*c_l*c_l*GuoF[16]-1.00000000000000000*c_l*c_l*c_l*GuoF[17]+1.00000000000000000*c_l*c_l*c_l*GuoF[18];

F_hat[8]*=(1-0.5*S[8]);
m_l[8]=m_l[8]-S[8]*(m_l[8]-meq[8])+dt*F_hat[8];
//=======================================

m_l[9]=+0.00000000000000000*c_l*c_l*f[ci][0]+2.00000000000000000*c_l*c_l*f[ci][1]+2.00000000000000000*c_l*c_l*f[ci][2]-1.00000000000000000*c_l*c_l*f[ci][3]-1.00000000000000000*c_l*c_l*f[ci][4]-1.00000000000000000*c_l*c_l*f[ci][5]-1.00000000000000000*c_l*c_l*f[ci][6]+1.00000000000000000*c_l*c_l*f[ci][7]+1.00000000000000000*c_l*c_l*f[ci][8]+1.00000000000000000*c_l*c_l*f[ci][9]+1.00000000000000000*c_l*c_l*f[ci][10]+1.00000000000000000*c_l*c_l*f[ci][11]+1.00000000000000000*c_l*c_l*f[ci][12]+1.00000000000000000*c_l*c_l*f[ci][13]+1.00000000000000000*c_l*c_l*f[ci][14]-2.00000000000000000*c_l*c_l*f[ci][15]-2.00000000000000000*c_l*c_l*f[ci][16]-2.00000000000000000*c_l*c_l*f[ci][17]-2.00000000000000000*c_l*c_l*f[ci][18];

F_hat[9]=+0.00000000000000000*c_l*c_l*GuoF[0]+2.00000000000000000*c_l*c_l*GuoF[1]+2.00000000000000000*c_l*c_l*GuoF[2]-1.00000000000000000*c_l*c_l*GuoF[3]-1.00000000000000000*c_l*c_l*GuoF[4]-1.00000000000000000*c_l*c_l*GuoF[5]-1.00000000000000000*c_l*c_l*GuoF[6]+1.00000000000000000*c_l*c_l*GuoF[7]+1.00000000000000000*c_l*c_l*GuoF[8]+1.00000000000000000*c_l*c_l*GuoF[9]+1.00000000000000000*c_l*c_l*GuoF[10]+1.00000000000000000*c_l*c_l*GuoF[11]+1.00000000000000000*c_l*c_l*GuoF[12]+1.00000000000000000*c_l*c_l*GuoF[13]+1.00000000000000000*c_l*c_l*GuoF[14]-2.00000000000000000*c_l*c_l*GuoF[15]-2.00000000000000000*c_l*c_l*GuoF[16]-2.00000000000000000*c_l*c_l*GuoF[17]-2.00000000000000000*c_l*c_l*GuoF[18];

F_hat[9]*=(1-0.5*S[9]);
m_l[9]=m_l[9]-S[9]*(m_l[9]-meq[9])+dt*F_hat[9];
//=======================================

m_l[10]=+0.00000000000000000*c_l*c_l*c_l*c_l*f[ci][0]-2.00000000000000000*c_l*c_l*c_l*c_l*f[ci][1]-2.00000000000000000*c_l*c_l*c_l*c_l*f[ci][2]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][3]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][4]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][5]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][6]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][7]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][8]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][9]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][10]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][11]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][12]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][13]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][14]-2.00000000000000000*c_l*c_l*c_l*c_l*f[ci][15]-2.00000000000000000*c_l*c_l*c_l*c_l*f[ci][16]-2.00000000000000000*c_l*c_l*c_l*c_l*f[ci][17]-2.00000000000000000*c_l*c_l*c_l*c_l*f[ci][18];

F_hat[10]=+0.00000000000000000*c_l*c_l*c_l*c_l*GuoF[0]-2.00000000000000000*c_l*c_l*c_l*c_l*GuoF[1]-2.00000000000000000*c_l*c_l*c_l*c_l*GuoF[2]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[3]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[4]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[5]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[6]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[7]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[8]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[9]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[10]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[11]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[12]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[13]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[14]-2.00000000000000000*c_l*c_l*c_l*c_l*GuoF[15]-2.00000000000000000*c_l*c_l*c_l*c_l*GuoF[16]-2.00000000000000000*c_l*c_l*c_l*c_l*GuoF[17]-2.00000000000000000*c_l*c_l*c_l*c_l*GuoF[18];

F_hat[10]*=(1-0.5*S[10]);
m_l[10]=m_l[10]-S[10]*(m_l[10]-meq[10])+dt*F_hat[10];
//=======================================

m_l[11]=+0.00000000000000000*c_l*c_l*f[ci][0]+0.00000000000000000*c_l*c_l*f[ci][1]+0.00000000000000000*c_l*c_l*f[ci][2]+1.00000000000000000*c_l*c_l*f[ci][3]+1.00000000000000000*c_l*c_l*f[ci][4]-1.00000000000000000*c_l*c_l*f[ci][5]-1.00000000000000000*c_l*c_l*f[ci][6]+1.00000000000000000*c_l*c_l*f[ci][7]+1.00000000000000000*c_l*c_l*f[ci][8]+1.00000000000000000*c_l*c_l*f[ci][9]+1.00000000000000000*c_l*c_l*f[ci][10]-1.00000000000000000*c_l*c_l*f[ci][11]-1.00000000000000000*c_l*c_l*f[ci][12]-1.00000000000000000*c_l*c_l*f[ci][13]-1.00000000000000000*c_l*c_l*f[ci][14]+0.00000000000000000*c_l*c_l*f[ci][15]+0.00000000000000000*c_l*c_l*f[ci][16]+0.00000000000000000*c_l*c_l*f[ci][17]+0.00000000000000000*c_l*c_l*f[ci][18];

F_hat[11]=+0.00000000000000000*c_l*c_l*GuoF[0]+0.00000000000000000*c_l*c_l*GuoF[1]+0.00000000000000000*c_l*c_l*GuoF[2]+1.00000000000000000*c_l*c_l*GuoF[3]+1.00000000000000000*c_l*c_l*GuoF[4]-1.00000000000000000*c_l*c_l*GuoF[5]-1.00000000000000000*c_l*c_l*GuoF[6]+1.00000000000000000*c_l*c_l*GuoF[7]+1.00000000000000000*c_l*c_l*GuoF[8]+1.00000000000000000*c_l*c_l*GuoF[9]+1.00000000000000000*c_l*c_l*GuoF[10]-1.00000000000000000*c_l*c_l*GuoF[11]-1.00000000000000000*c_l*c_l*GuoF[12]-1.00000000000000000*c_l*c_l*GuoF[13]-1.00000000000000000*c_l*c_l*GuoF[14]+0.00000000000000000*c_l*c_l*GuoF[15]+0.00000000000000000*c_l*c_l*GuoF[16]+0.00000000000000000*c_l*c_l*GuoF[17]+0.00000000000000000*c_l*c_l*GuoF[18];

F_hat[11]*=(1-0.5*S[11]);
m_l[11]=m_l[11]-S[11]*(m_l[11]-meq[11])+dt*F_hat[11];
//=======================================

m_l[12]=+0.00000000000000000*c_l*c_l*c_l*c_l*f[ci][0]+0.00000000000000000*c_l*c_l*c_l*c_l*f[ci][1]+0.00000000000000000*c_l*c_l*c_l*c_l*f[ci][2]-1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][3]-1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][4]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][5]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][6]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][7]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][8]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][9]+1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][10]-1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][11]-1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][12]-1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][13]-1.00000000000000000*c_l*c_l*c_l*c_l*f[ci][14]+0.00000000000000000*c_l*c_l*c_l*c_l*f[ci][15]+0.00000000000000000*c_l*c_l*c_l*c_l*f[ci][16]+0.00000000000000000*c_l*c_l*c_l*c_l*f[ci][17]+0.00000000000000000*c_l*c_l*c_l*c_l*f[ci][18];

F_hat[12]=+0.00000000000000000*c_l*c_l*c_l*c_l*GuoF[0]+0.00000000000000000*c_l*c_l*c_l*c_l*GuoF[1]+0.00000000000000000*c_l*c_l*c_l*c_l*GuoF[2]-1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[3]-1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[4]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[5]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[6]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[7]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[8]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[9]+1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[10]-1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[11]-1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[12]-1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[13]-1.00000000000000000*c_l*c_l*c_l*c_l*GuoF[14]+0.00000000000000000*c_l*c_l*c_l*c_l*GuoF[15]+0.00000000000000000*c_l*c_l*c_l*c_l*GuoF[16]+0.00000000000000000*c_l*c_l*c_l*c_l*GuoF[17]+0.00000000000000000*c_l*c_l*c_l*c_l*GuoF[18];

F_hat[12]*=(1-0.5*S[12]);
m_l[12]=m_l[12]-S[12]*(m_l[12]-meq[12])+dt*F_hat[12];
//=======================================

m_l[13]=+0.00000000000000000*c_l*c_l*f[ci][0]+0.00000000000000000*c_l*c_l*f[ci][1]+0.00000000000000000*c_l*c_l*f[ci][2]+0.00000000000000000*c_l*c_l*f[ci][3]+0.00000000000000000*c_l*c_l*f[ci][4]+0.00000000000000000*c_l*c_l*f[ci][5]+0.00000000000000000*c_l*c_l*f[ci][6]+1.00000000000000000*c_l*c_l*f[ci][7]+1.00000000000000000*c_l*c_l*f[ci][8]-1.00000000000000000*c_l*c_l*f[ci][9]-1.00000000000000000*c_l*c_l*f[ci][10]+0.00000000000000000*c_l*c_l*f[ci][11]+0.00000000000000000*c_l*c_l*f[ci][12]+0.00000000000000000*c_l*c_l*f[ci][13]+0.00000000000000000*c_l*c_l*f[ci][14]+0.00000000000000000*c_l*c_l*f[ci][15]+0.00000000000000000*c_l*c_l*f[ci][16]+0.00000000000000000*c_l*c_l*f[ci][17]+0.00000000000000000*c_l*c_l*f[ci][18];

F_hat[13]=+0.00000000000000000*c_l*c_l*GuoF[0]+0.00000000000000000*c_l*c_l*GuoF[1]+0.00000000000000000*c_l*c_l*GuoF[2]+0.00000000000000000*c_l*c_l*GuoF[3]+0.00000000000000000*c_l*c_l*GuoF[4]+0.00000000000000000*c_l*c_l*GuoF[5]+0.00000000000000000*c_l*c_l*GuoF[6]+1.00000000000000000*c_l*c_l*GuoF[7]+1.00000000000000000*c_l*c_l*GuoF[8]-1.00000000000000000*c_l*c_l*GuoF[9]-1.00000000000000000*c_l*c_l*GuoF[10]+0.00000000000000000*c_l*c_l*GuoF[11]+0.00000000000000000*c_l*c_l*GuoF[12]+0.00000000000000000*c_l*c_l*GuoF[13]+0.00000000000000000*c_l*c_l*GuoF[14]+0.00000000000000000*c_l*c_l*GuoF[15]+0.00000000000000000*c_l*c_l*GuoF[16]+0.00000000000000000*c_l*c_l*GuoF[17]+0.00000000000000000*c_l*c_l*GuoF[18];

F_hat[13]*=(1-0.5*S[13]);
m_l[13]=m_l[13]-S[13]*(m_l[13]-meq[13])+dt*F_hat[13];
//=======================================

m_l[14]=+0.00000000000000000*c_l*c_l*f[ci][0]+0.00000000000000000*c_l*c_l*f[ci][1]+0.00000000000000000*c_l*c_l*f[ci][2]+0.00000000000000000*c_l*c_l*f[ci][3]+0.00000000000000000*c_l*c_l*f[ci][4]+0.00000000000000000*c_l*c_l*f[ci][5]+0.00000000000000000*c_l*c_l*f[ci][6]+0.00000000000000000*c_l*c_l*f[ci][7]+0.00000000000000000*c_l*c_l*f[ci][8]+0.00000000000000000*c_l*c_l*f[ci][9]+0.00000000000000000*c_l*c_l*f[ci][10]+0.00000000000000000*c_l*c_l*f[ci][11]+0.00000000000000000*c_l*c_l*f[ci][12]+0.00000000000000000*c_l*c_l*f[ci][13]+0.00000000000000000*c_l*c_l*f[ci][14]+1.00000000000000000*c_l*c_l*f[ci][15]+1.00000000000000000*c_l*c_l*f[ci][16]-1.00000000000000000*c_l*c_l*f[ci][17]-1.00000000000000000*c_l*c_l*f[ci][18];

F_hat[14]=+0.00000000000000000*c_l*c_l*GuoF[0]+0.00000000000000000*c_l*c_l*GuoF[1]+0.00000000000000000*c_l*c_l*GuoF[2]+0.00000000000000000*c_l*c_l*GuoF[3]+0.00000000000000000*c_l*c_l*GuoF[4]+0.00000000000000000*c_l*c_l*GuoF[5]+0.00000000000000000*c_l*c_l*GuoF[6]+0.00000000000000000*c_l*c_l*GuoF[7]+0.00000000000000000*c_l*c_l*GuoF[8]+0.00000000000000000*c_l*c_l*GuoF[9]+0.00000000000000000*c_l*c_l*GuoF[10]+0.00000000000000000*c_l*c_l*GuoF[11]+0.00000000000000000*c_l*c_l*GuoF[12]+0.00000000000000000*c_l*c_l*GuoF[13]+0.00000000000000000*c_l*c_l*GuoF[14]+1.00000000000000000*c_l*c_l*GuoF[15]+1.00000000000000000*c_l*c_l*GuoF[16]-1.00000000000000000*c_l*c_l*GuoF[17]-1.00000000000000000*c_l*c_l*GuoF[18];

F_hat[14]*=(1-0.5*S[14]);
m_l[14]=m_l[14]-S[14]*(m_l[14]-meq[14])+dt*F_hat[14];
//=======================================

m_l[15]=+0.00000000000000000*c_l*c_l*f[ci][0]+0.00000000000000000*c_l*c_l*f[ci][1]+0.00000000000000000*c_l*c_l*f[ci][2]+0.00000000000000000*c_l*c_l*f[ci][3]+0.00000000000000000*c_l*c_l*f[ci][4]+0.00000000000000000*c_l*c_l*f[ci][5]+0.00000000000000000*c_l*c_l*f[ci][6]+0.00000000000000000*c_l*c_l*f[ci][7]+0.00000000000000000*c_l*c_l*f[ci][8]+0.00000000000000000*c_l*c_l*f[ci][9]+0.00000000000000000*c_l*c_l*f[ci][10]+1.00000000000000000*c_l*c_l*f[ci][11]+1.00000000000000000*c_l*c_l*f[ci][12]-1.00000000000000000*c_l*c_l*f[ci][13]-1.00000000000000000*c_l*c_l*f[ci][14]+0.00000000000000000*c_l*c_l*f[ci][15]+0.00000000000000000*c_l*c_l*f[ci][16]+0.00000000000000000*c_l*c_l*f[ci][17]+0.00000000000000000*c_l*c_l*f[ci][18];

F_hat[15]=+0.00000000000000000*c_l*c_l*GuoF[0]+0.00000000000000000*c_l*c_l*GuoF[1]+0.00000000000000000*c_l*c_l*GuoF[2]+0.00000000000000000*c_l*c_l*GuoF[3]+0.00000000000000000*c_l*c_l*GuoF[4]+0.00000000000000000*c_l*c_l*GuoF[5]+0.00000000000000000*c_l*c_l*GuoF[6]+0.00000000000000000*c_l*c_l*GuoF[7]+0.00000000000000000*c_l*c_l*GuoF[8]+0.00000000000000000*c_l*c_l*GuoF[9]+0.00000000000000000*c_l*c_l*GuoF[10]+1.00000000000000000*c_l*c_l*GuoF[11]+1.00000000000000000*c_l*c_l*GuoF[12]-1.00000000000000000*c_l*c_l*GuoF[13]-1.00000000000000000*c_l*c_l*GuoF[14]+0.00000000000000000*c_l*c_l*GuoF[15]+0.00000000000000000*c_l*c_l*GuoF[16]+0.00000000000000000*c_l*c_l*GuoF[17]+0.00000000000000000*c_l*c_l*GuoF[18];

F_hat[15]*=(1-0.5*S[15]);
m_l[15]=m_l[15]-S[15]*(m_l[15]-meq[15])+dt*F_hat[15];
//=======================================

m_l[16]=+0.00000000000000000*c_l*c_l*c_l*f[ci][0]+0.00000000000000000*c_l*c_l*c_l*f[ci][1]+0.00000000000000000*c_l*c_l*c_l*f[ci][2]+0.00000000000000000*c_l*c_l*c_l*f[ci][3]+0.00000000000000000*c_l*c_l*c_l*f[ci][4]+0.00000000000000000*c_l*c_l*c_l*f[ci][5]+0.00000000000000000*c_l*c_l*c_l*f[ci][6]+1.00000000000000000*c_l*c_l*c_l*f[ci][7]-1.00000000000000000*c_l*c_l*c_l*f[ci][8]+1.00000000000000000*c_l*c_l*c_l*f[ci][9]-1.00000000000000000*c_l*c_l*c_l*f[ci][10]-1.00000000000000000*c_l*c_l*c_l*f[ci][11]+1.00000000000000000*c_l*c_l*c_l*f[ci][12]-1.00000000000000000*c_l*c_l*c_l*f[ci][13]+1.00000000000000000*c_l*c_l*c_l*f[ci][14]+0.00000000000000000*c_l*c_l*c_l*f[ci][15]+0.00000000000000000*c_l*c_l*c_l*f[ci][16]+0.00000000000000000*c_l*c_l*c_l*f[ci][17]+0.00000000000000000*c_l*c_l*c_l*f[ci][18];

F_hat[16]=+0.00000000000000000*c_l*c_l*c_l*GuoF[0]+0.00000000000000000*c_l*c_l*c_l*GuoF[1]+0.00000000000000000*c_l*c_l*c_l*GuoF[2]+0.00000000000000000*c_l*c_l*c_l*GuoF[3]+0.00000000000000000*c_l*c_l*c_l*GuoF[4]+0.00000000000000000*c_l*c_l*c_l*GuoF[5]+0.00000000000000000*c_l*c_l*c_l*GuoF[6]+1.00000000000000000*c_l*c_l*c_l*GuoF[7]-1.00000000000000000*c_l*c_l*c_l*GuoF[8]+1.00000000000000000*c_l*c_l*c_l*GuoF[9]-1.00000000000000000*c_l*c_l*c_l*GuoF[10]-1.00000000000000000*c_l*c_l*c_l*GuoF[11]+1.00000000000000000*c_l*c_l*c_l*GuoF[12]-1.00000000000000000*c_l*c_l*c_l*GuoF[13]+1.00000000000000000*c_l*c_l*c_l*GuoF[14]+0.00000000000000000*c_l*c_l*c_l*GuoF[15]+0.00000000000000000*c_l*c_l*c_l*GuoF[16]+0.00000000000000000*c_l*c_l*c_l*GuoF[17]+0.00000000000000000*c_l*c_l*c_l*GuoF[18];

F_hat[16]*=(1-0.5*S[16]);
m_l[16]=m_l[16]-S[16]*(m_l[16]-meq[16])+dt*F_hat[16];
//=======================================

m_l[17]=+0.00000000000000000*c_l*c_l*c_l*f[ci][0]+0.00000000000000000*c_l*c_l*c_l*f[ci][1]+0.00000000000000000*c_l*c_l*c_l*f[ci][2]+0.00000000000000000*c_l*c_l*c_l*f[ci][3]+0.00000000000000000*c_l*c_l*c_l*f[ci][4]+0.00000000000000000*c_l*c_l*c_l*f[ci][5]+0.00000000000000000*c_l*c_l*c_l*f[ci][6]-1.00000000000000000*c_l*c_l*c_l*f[ci][7]+1.00000000000000000*c_l*c_l*c_l*f[ci][8]+1.00000000000000000*c_l*c_l*c_l*f[ci][9]-1.00000000000000000*c_l*c_l*c_l*f[ci][10]+0.00000000000000000*c_l*c_l*c_l*f[ci][11]+0.00000000000000000*c_l*c_l*c_l*f[ci][12]+0.00000000000000000*c_l*c_l*c_l*f[ci][13]+0.00000000000000000*c_l*c_l*c_l*f[ci][14]+1.00000000000000000*c_l*c_l*c_l*f[ci][15]-1.00000000000000000*c_l*c_l*c_l*f[ci][16]+1.00000000000000000*c_l*c_l*c_l*f[ci][17]-1.00000000000000000*c_l*c_l*c_l*f[ci][18];

F_hat[17]=+0.00000000000000000*c_l*c_l*c_l*GuoF[0]+0.00000000000000000*c_l*c_l*c_l*GuoF[1]+0.00000000000000000*c_l*c_l*c_l*GuoF[2]+0.00000000000000000*c_l*c_l*c_l*GuoF[3]+0.00000000000000000*c_l*c_l*c_l*GuoF[4]+0.00000000000000000*c_l*c_l*c_l*GuoF[5]+0.00000000000000000*c_l*c_l*c_l*GuoF[6]-1.00000000000000000*c_l*c_l*c_l*GuoF[7]+1.00000000000000000*c_l*c_l*c_l*GuoF[8]+1.00000000000000000*c_l*c_l*c_l*GuoF[9]-1.00000000000000000*c_l*c_l*c_l*GuoF[10]+0.00000000000000000*c_l*c_l*c_l*GuoF[11]+0.00000000000000000*c_l*c_l*c_l*GuoF[12]+0.00000000000000000*c_l*c_l*c_l*GuoF[13]+0.00000000000000000*c_l*c_l*c_l*GuoF[14]+1.00000000000000000*c_l*c_l*c_l*GuoF[15]-1.00000000000000000*c_l*c_l*c_l*GuoF[16]+1.00000000000000000*c_l*c_l*c_l*GuoF[17]-1.00000000000000000*c_l*c_l*c_l*GuoF[18];

F_hat[17]*=(1-0.5*S[17]);
m_l[17]=m_l[17]-S[17]*(m_l[17]-meq[17])+dt*F_hat[17];
//=======================================

m_l[18]=+0.00000000000000000*c_l*c_l*c_l*f[ci][0]+0.00000000000000000*c_l*c_l*c_l*f[ci][1]+0.00000000000000000*c_l*c_l*c_l*f[ci][2]+0.00000000000000000*c_l*c_l*c_l*f[ci][3]+0.00000000000000000*c_l*c_l*c_l*f[ci][4]+0.00000000000000000*c_l*c_l*c_l*f[ci][5]+0.00000000000000000*c_l*c_l*c_l*f[ci][6]+0.00000000000000000*c_l*c_l*c_l*f[ci][7]+0.00000000000000000*c_l*c_l*c_l*f[ci][8]+0.00000000000000000*c_l*c_l*c_l*f[ci][9]+0.00000000000000000*c_l*c_l*c_l*f[ci][10]+1.00000000000000000*c_l*c_l*c_l*f[ci][11]-1.00000000000000000*c_l*c_l*c_l*f[ci][12]-1.00000000000000000*c_l*c_l*c_l*f[ci][13]+1.00000000000000000*c_l*c_l*c_l*f[ci][14]-1.00000000000000000*c_l*c_l*c_l*f[ci][15]+1.00000000000000000*c_l*c_l*c_l*f[ci][16]+1.00000000000000000*c_l*c_l*c_l*f[ci][17]-1.00000000000000000*c_l*c_l*c_l*f[ci][18];

F_hat[18]=+0.00000000000000000*c_l*c_l*c_l*GuoF[0]+0.00000000000000000*c_l*c_l*c_l*GuoF[1]+0.00000000000000000*c_l*c_l*c_l*GuoF[2]+0.00000000000000000*c_l*c_l*c_l*GuoF[3]+0.00000000000000000*c_l*c_l*c_l*GuoF[4]+0.00000000000000000*c_l*c_l*c_l*GuoF[5]+0.00000000000000000*c_l*c_l*c_l*GuoF[6]+0.00000000000000000*c_l*c_l*c_l*GuoF[7]+0.00000000000000000*c_l*c_l*c_l*GuoF[8]+0.00000000000000000*c_l*c_l*c_l*GuoF[9]+0.00000000000000000*c_l*c_l*c_l*GuoF[10]+1.00000000000000000*c_l*c_l*c_l*GuoF[11]-1.00000000000000000*c_l*c_l*c_l*GuoF[12]-1.00000000000000000*c_l*c_l*c_l*GuoF[13]+1.00000000000000000*c_l*c_l*c_l*GuoF[14]-1.00000000000000000*c_l*c_l*c_l*GuoF[15]+1.00000000000000000*c_l*c_l*c_l*GuoF[16]+1.00000000000000000*c_l*c_l*c_l*GuoF[17]-1.00000000000000000*c_l*c_l*c_l*GuoF[18];

F_hat[18]*=(1-0.5*S[18]);
m_l[18]=m_l[18]-S[18]*(m_l[18]-meq[18])+dt*F_hat[18];
//=======================================


//==========================
m_inv_l[0]=+((double)0X1.5555555555555P-2)*1.0*m_l[0]+((double)-0X1P-1)*1.0/(c_l*c_l)*m_l[1]+((double)0X1.5555555555555P-3)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)0X0P+0)*1.0/c_l*m_l[3]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)0X0P+0)*1.0/(c_l)*m_l[5]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)0X0P+0)*1.0/c_l*m_l[7]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[9]+((double)0X0P+0)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[11]+((double)0X0P+0)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[1]=+((double)0X1.C71C71C71C71CP-5)*1.0*m_l[0]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[1]+((double)-0X1.C71C71C71C71EP-5)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)0X1.5555555555555P-3)*1.0/c_l*m_l[3]+((double)-0X1.5555555555556P-3)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)0X0P+0)*1.0/(c_l)*m_l[5]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)0X0P+0)*1.0/c_l*m_l[7]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)0X1.5555555555556P-4)*1.0/(c_l*c_l)*m_l[9]+((double)-0X1.5555555555554P-4)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[11]+((double)0X0P+0)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[2]=+((double)0X1.C71C71C71C71DP-5)*1.0*m_l[0]+((double)-0X1P-56)*1.0/(c_l*c_l)*m_l[1]+((double)-0X1.C71C71C71C71DP-5)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)-0X1.5555555555555P-3)*1.0/c_l*m_l[3]+((double)0X1.5555555555556P-3)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)0X0P+0)*1.0/(c_l)*m_l[5]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)0X0P+0)*1.0/c_l*m_l[7]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)0X1.5555555555556P-4)*1.0/(c_l*c_l)*m_l[9]+((double)-0X1.5555555555554P-4)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[11]+((double)0X0P+0)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[3]=+((double)0X1.C71C71C71C71DP-5)*1.0*m_l[0]+((double)0X1P-56)*1.0/(c_l*c_l)*m_l[1]+((double)-0X1.C71C71C71C71EP-5)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)0X0P+0)*1.0/c_l*m_l[3]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)0X1.5555555555554P-3)*1.0/(c_l)*m_l[5]+((double)-0X1.5555555555556P-3)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)0X0P+0)*1.0/c_l*m_l[7]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)-0X1.5555555555555P-5)*1.0/(c_l*c_l)*m_l[9]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)0X1.FFFFFFFFFFFFEP-4)*1.0/(c_l*c_l)*m_l[11]+((double)-0X1.0000000000001P-3)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)0X1P-56)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[4]=+((double)0X1.C71C71C71C71DP-5)*1.0*m_l[0]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[1]+((double)-0X1.C71C71C71C71AP-5)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)0X0P+0)*1.0/c_l*m_l[3]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)-0X1.5555555555555P-3)*1.0/(c_l)*m_l[5]+((double)0X1.5555555555556P-3)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)0X0P+0)*1.0/c_l*m_l[7]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)-0X1.5555555555555P-5)*1.0/(c_l*c_l)*m_l[9]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)0X1.FFFFFFFFFFFFEP-4)*1.0/(c_l*c_l)*m_l[11]+((double)-0X1.0000000000001P-3)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)-0X1P-56)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[5]=+((double)0X1.C71C71C71C718P-5)*1.0*m_l[0]+((double)-0X1.8P-56)*1.0/(c_l*c_l)*m_l[1]+((double)-0X1.C71C71C71C71DP-5)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)0X0P+0)*1.0/c_l*m_l[3]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)0X0P+0)*1.0/(c_l)*m_l[5]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)0X1.5555555555555P-3)*1.0/c_l*m_l[7]+((double)-0X1.5555555555556P-3)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)-0X1.5555555555558P-5)*1.0/(c_l*c_l)*m_l[9]+((double)0X1.5555555555552P-5)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)-0X1.FFFFFFFFFFFFDP-4)*1.0/(c_l*c_l)*m_l[11]+((double)0X1.0000000000002P-3)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[6]=+((double)0X1.C71C71C71C71AP-5)*1.0*m_l[0]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[1]+((double)-0X1.C71C71C71C71DP-5)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)0X0P+0)*1.0/c_l*m_l[3]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)-0X1.5555555555554P-58)*1.0/(c_l)*m_l[5]+((double)-0X1.5555555555554P-59)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)-0X1.5555555555555P-3)*1.0/c_l*m_l[7]+((double)0X1.5555555555555P-3)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)-0X1.5555555555556P-5)*1.0/(c_l*c_l)*m_l[9]+((double)0X1.5555555555554P-5)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)-0X1.FFFFFFFFFFFFEP-4)*1.0/(c_l*c_l)*m_l[11]+((double)0X1P-3)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)-0X1P-55)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[7]=+((double)0X1.C71C71C71C71BP-6)*1.0*m_l[0]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l)*m_l[1]+((double)0X1.C71C71C71C71BP-7)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)0X1.5555555555555P-4)*1.0/c_l*m_l[3]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)0X1.5555555555555P-4)*1.0/(c_l)*m_l[5]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)0X0P+0)*1.0/c_l*m_l[7]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l)*m_l[9]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)0X1.0000000000001P-4)*1.0/(c_l*c_l)*m_l[11]+((double)0X1.0000000000001P-4)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X1P-2)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)-0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[8]=+((double)0X1.C71C71C71C71BP-6)*1.0*m_l[0]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l)*m_l[1]+((double)0X1.C71C71C71C71BP-7)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)-0X1.5555555555555P-4)*1.0/c_l*m_l[3]+((double)-0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)-0X1.5555555555555P-4)*1.0/(c_l)*m_l[5]+((double)-0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)0X0P+0)*1.0/c_l*m_l[7]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l)*m_l[9]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)0X1.0000000000001P-4)*1.0/(c_l*c_l)*m_l[11]+((double)0X1.0000000000001P-4)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X1P-2)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)-0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[9]=+((double)0X1.C71C71C71C71BP-6)*1.0*m_l[0]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l)*m_l[1]+((double)0X1.C71C71C71C71BP-7)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)0X1.5555555555555P-4)*1.0/c_l*m_l[3]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)-0X1.5555555555555P-4)*1.0/(c_l)*m_l[5]+((double)-0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)0X0P+0)*1.0/c_l*m_l[7]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l)*m_l[9]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)0X1.0000000000001P-4)*1.0/(c_l*c_l)*m_l[11]+((double)0X1.0000000000001P-4)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)-0X1P-2)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[10]=+((double)0X1.C71C71C71C71BP-6)*1.0*m_l[0]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l)*m_l[1]+((double)0X1.C71C71C71C71BP-7)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)-0X1.5555555555555P-4)*1.0/c_l*m_l[3]+((double)-0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)0X1.5555555555555P-4)*1.0/(c_l)*m_l[5]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)0X0P+0)*1.0/c_l*m_l[7]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l)*m_l[9]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)0X1.0000000000001P-4)*1.0/(c_l*c_l)*m_l[11]+((double)0X1.0000000000001P-4)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)-0X1P-2)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)-0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)-0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[11]=+((double)0X1.C71C71C71C71DP-6)*1.0*m_l[0]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l)*m_l[1]+((double)0X1.C71C71C71C71DP-7)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)0X1.5555555555555P-4)*1.0/c_l*m_l[3]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)0X0P+0)*1.0/(c_l)*m_l[5]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)0X1.5555555555555P-4)*1.0/c_l*m_l[7]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l)*m_l[9]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)-0X1.0000000000001P-4)*1.0/(c_l*c_l)*m_l[11]+((double)-0X1.0000000000001P-4)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)0X1P-2)*1.0/(c_l*c_l)*m_l[15]+((double)-0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[12]=+((double)0X1.C71C71C71C71DP-6)*1.0*m_l[0]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l)*m_l[1]+((double)0X1.C71C71C71C71DP-7)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)-0X1.5555555555555P-4)*1.0/c_l*m_l[3]+((double)-0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)0X0P+0)*1.0/(c_l)*m_l[5]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)-0X1.5555555555555P-4)*1.0/c_l*m_l[7]+((double)-0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l)*m_l[9]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)-0X1.0000000000001P-4)*1.0/(c_l*c_l)*m_l[11]+((double)-0X1.0000000000001P-4)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)0X1P-2)*1.0/(c_l*c_l)*m_l[15]+((double)0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)-0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[13]=+((double)0X1.C71C71C71C71DP-6)*1.0*m_l[0]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l)*m_l[1]+((double)0X1.C71C71C71C71DP-7)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)0X1.5555555555555P-4)*1.0/c_l*m_l[3]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)0X0P+0)*1.0/(c_l)*m_l[5]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)-0X1.5555555555555P-4)*1.0/c_l*m_l[7]+((double)-0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l)*m_l[9]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)-0X1.0000000000001P-4)*1.0/(c_l*c_l)*m_l[11]+((double)-0X1.0000000000001P-4)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)-0X1P-2)*1.0/(c_l*c_l)*m_l[15]+((double)-0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)-0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[14]=+((double)0X1.C71C71C71C71DP-6)*1.0*m_l[0]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l)*m_l[1]+((double)0X1.C71C71C71C71DP-7)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)-0X1.5555555555555P-4)*1.0/c_l*m_l[3]+((double)-0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)0X0P+0)*1.0/(c_l)*m_l[5]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)0X1.5555555555555P-4)*1.0/c_l*m_l[7]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l)*m_l[9]+((double)0X1.5555555555555P-6)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)-0X1.0000000000001P-4)*1.0/(c_l*c_l)*m_l[11]+((double)-0X1.0000000000001P-4)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[14]+((double)-0X1P-2)*1.0/(c_l*c_l)*m_l[15]+((double)0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[15]=+((double)0X1.C71C71C71C71CP-6)*1.0*m_l[0]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l)*m_l[1]+((double)0X1.C71C71C71C71CP-7)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)0X0P+0)*1.0/c_l*m_l[3]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)0X1.5555555555555P-4)*1.0/(c_l)*m_l[5]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)0X1.5555555555555P-4)*1.0/c_l*m_l[7]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)-0X1.5555555555554P-5)*1.0/(c_l*c_l)*m_l[9]+((double)-0X1.5555555555554P-5)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[11]+((double)0X0P+0)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)0X1P-2)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)-0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[16]=+((double)0X1.C71C71C71C71AP-6)*1.0*m_l[0]+((double)0X1.5555555555557P-5)*1.0/(c_l*c_l)*m_l[1]+((double)0X1.C71C71C71C71AP-7)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)0X0P+0)*1.0/c_l*m_l[3]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)-0X1.5555555555555P-4)*1.0/(c_l)*m_l[5]+((double)-0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)-0X1.5555555555555P-4)*1.0/c_l*m_l[7]+((double)-0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)-0X1.5555555555554P-5)*1.0/(c_l*c_l)*m_l[9]+((double)-0X1.5555555555554P-5)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[11]+((double)0X0P+0)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)0X1P-2)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)-0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[17]=+((double)0X1.C71C71C71C71CP-6)*1.0*m_l[0]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l)*m_l[1]+((double)0X1.C71C71C71C71CP-7)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)0X0P+0)*1.0/c_l*m_l[3]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)0X1.5555555555555P-4)*1.0/(c_l)*m_l[5]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)-0X1.5555555555555P-4)*1.0/c_l*m_l[7]+((double)-0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)-0X1.5555555555556P-5)*1.0/(c_l*c_l)*m_l[9]+((double)-0X1.5555555555556P-5)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[11]+((double)0X0P+0)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)-0X1P-2)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[18];

m_inv_l[18]=+((double)0X1.C71C71C71C71CP-6)*1.0*m_l[0]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l)*m_l[1]+((double)0X1.C71C71C71C71CP-7)*1.0/(c_l*c_l*c_l*c_l)*m_l[2]+((double)-0X0P+0)*1.0/c_l*m_l[3]+((double)-0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[4]+((double)-0X1.5555555555555P-4)*1.0/(c_l)*m_l[5]+((double)-0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[6]+((double)0X1.5555555555555P-4)*1.0/c_l*m_l[7]+((double)0X1.5555555555555P-5)*1.0/(c_l*c_l*c_l)*m_l[8]+((double)-0X1.5555555555556P-5)*1.0/(c_l*c_l)*m_l[9]+((double)-0X1.5555555555556P-5)*1.0/(c_l*c_l*c_l*c_l)*m_l[10]+((double)-0X0P+0)*1.0/(c_l*c_l)*m_l[11]+((double)-0X0P+0)*1.0/(c_l*c_l*c_l*c_l)*m_l[12]+((double)-0X0P+0)*1.0/(c_l*c_l)*m_l[13]+((double)-0X1P-2)*1.0/(c_l*c_l)*m_l[14]+((double)0X0P+0)*1.0/(c_l*c_l)*m_l[15]+((double)-0X0P+0)*1.0/(c_l*c_l*c_l)*m_l[16]+((double)-0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[17]+((double)-0X1P-3)*1.0/(c_l*c_l*c_l)*m_l[18];

//====================
		//for (int mis=0; mis<19; mis++)
                //        f[ci][mis]=m_inv_l[mis];

		for (int mi=0; mi<19; mi++)
			{
			        //mi=SWAPE_INV[mis];
			 
			sum=m_inv_l[mi];
			//==============================================================================
				
			
			
			//F[ci][mi]=0;
			ip=i+e[mi][0];
			jp=j+e[mi][1];if (jp<0) {jp=NY;}; if (jp>NY) {jp=0;};
			kp=m+e[mi][2];if (kp<0) {kp=NZ;}; if (kp>NZ) {kp=0;};


			if (ip<0) 
				if (Sl[jp*(NZ+1)+kp]>0)
				{
				//sendl[(Sl[jp*(NZ+1)+kp]-1)*5+FLN[mi]]=sum;
				//sendl_rhob[Sl[jp*(NZ+1)+kp]-1]+=g_b[lm];
				//cout<<g_r[lm]<<"    1"<<endl;
				sendl[(Sl[jp*(NZ+1)+kp]-1)*5+FLN[mi]]=sum;
				 //if (mis<=9)
				 //        f[ci][mi]=f[ci][LR[mi]];
				 //        f[ci][LR[mi]]=sum;
				}
				else
				{
				F[ci][LR[mi]]=sum;
				//if (mis<=9) 
			         //         f[ci][mi]=f[ci][LR[mi]];
				//f[ci][LR[mi]]=sum;
				}
					
						
					
					
			if (ip>=nx_l)
				if (Sr[jp*(NZ+1)+kp]>0)
				{
				//sendr[(Sr[jp*(NZ+1)+kp]-1)*5+FRP[mi]]=sum;
				//sendr_rhob[Sr[jp*(NZ+1)+kp]-1]+=g_b[lm];
				//cout<<g_r[lm]<<"    2"<<endl;
				sendr[(Sr[jp*(NZ+1)+kp]-1)*5+FRP[mi]]=sum;
				//if (mis<=9)
				 //        f[ci][mi]=f[ci][LR[mi]];
				 //        f[ci][LR[mi]]=sum;
				}
				else
				{
				F[ci][LR[mi]]=sum;
				//     if (mis<=9) 
				//                f[ci][mi]=f[ci][LR[mi]];
				//f[ci][LR[mi]]=sum;      
				}

			if ((ip>=0) and (ip<nx_l)) 
				if (Solid[ip][jp][kp]>0)
				{
				F[Solid[ip][jp][kp]][mi]=sum;
				//if (mis<=9)
				//        
				//        {
				          
				 //                       
				//                       f[ci][mi]=f[ci][LR[mi]];
				//                        f[ci][LR[mi]]=f[Solid[ip][jp][kp]][mi];
				//                        f[Solid[ip][jp][kp]][mi]=sum;
				        
				 //                
				//
				 //               }
				}
				else
				{
				F[ci][LR[mi]]=sum;
				//if (mis<=9) 
				//                f[ci][mi]=f[ci][LR[mi]];
				        
				//        f[ci][LR[mi]]=sum;
				}
		//=======================G streaming=================================================
		//for(int lm=0;lm<19;lm++)
                //{
                 eu=elat[mi][0]*u[ci][0]+elat[mi][1]*u[ci][1]+elat[mi][2]*u[ci][2];
                 //g_r[mi]=w[mi]*rho_r[ci]*(1+3*eu/c2);
                 //g_b[mi]=w[mi]*rho_b[ci]*(1+3*eu/c2);

		 g_r[mi]=w[mi]*rho_r[ci]*(1+3*eu/c2+4.5*eu*eu/c4-1.5*uu/c2); 
		//if (g_r[mi]<0) 
		//cout<<mi<<"  "<<g_r[mi]<<"  "<<n<<"  "<<rho_r[ci]<<endl;
                 g_b[mi]=w[mi]*rho_b[ci]*(1+3*eu/c2+4.5*eu*eu/c4-1.5*uu/c2);
		
			
		//	cout<<" "<<g_r[lm]<<" "<<g_b[lm]<<"  the number "<<n<<"  vector "<<lm<<endl;
                 }

		
                 
           if (cc>0)
           for(int kk=1;kk<19;kk+=2)
                {
                //ef=elat[kk][0]*C[0]+elat[kk][1]*C[1]+elat[kk][2]*C[2];
                ef=e[kk][0]*C[0]+e[kk][1]*C[1]+e[kk][2]*C[2];
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

	
	MPI_Waitall(12,request2, status2);

	MPI_Testall(12,request2,&mpi_test2,status2);	
		
			for(i=1;i<=Gcl[rank];i++)
			        {
			
			        rhor[i]+=recvl_rhor[i-1];
			        rhob[i]+=recvl_rhob[i-1];
				
				for (int lm=0;lm<5;lm++)
					if (recvl[(i-1)*5+lm]>0)
				        F[i][RP[lm]]=recvl[(i-1)*5+lm];
			        }
			for(j=Count-Gcr[rank]+1;j<=Count;j++)
			        {
				rhor[j]+=recvr_rhor[j-(Count-Gcr[rank]+1)];
				rhob[j]+=recvr_rhob[j-(Count-Gcr[rank]+1)];
				for (int lm=0;lm<5;lm++)
					if (recvr[(j-(Count-Gcr[rank]+1))*5+lm]>0)
			        	F[j][LN[lm]]=recvr[(j-(Count-Gcr[rank]+1))*5+lm];
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



void comput_macro_variables( double* rho,double** u,double** u0,double** f,double** F,double* rho_r, double* rho_b, double* rhor, double* rhob, double* psi,int* SupInv)
{
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	srand((unsigned)time(0)+rank);
	double rand_double,psi_rand;

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
				
				rho_r[i]=rhor[i];
				rho_b[i]=rhob[i];
				rhor[i]=0;
				rhob[i]=0;
				u[i][0]=(u[i][0]+dt*gxs)/rho[i];
				u[i][1]=(u[i][1]+dt*gys)/rho[i];
				u[i][2]=(u[i][2]+dt*gzs)/rho[i];
				
				
				psi[i]=(rho_r[i]-rho_b[i])/(rho_r[i]+rho_b[i]);
				
			}
			
		
	if (in_psi_BC==1)
                {
                   if ((psi_xn>0) and (rank==0))
			if (psi_xn==1)    
                           {
                           for(int j=0;j<=NY;j++)
                                   for (int k=0;k<=NZ;k++)
				   if (Solid[0][j][k]>0)
                                   {
                                   rho_r[Solid[0][j][k]]=(Psi_local[0][j][k]*1.0+1.0)/2;
                                   rho_b[Solid[0][j][k]]=1.0-rho_r[Solid[0][j][k]];
                                   psi[Solid[0][j][k]]=Psi_local[0][j][k];
				//cout<<Psi_local[0][j][k]<<"  "<<j<<"  "<<k<<endl;
                                   
                                   
                                   }
                           }
			else
			if (psi_xn==2)
			{
                           for(int j=0;j<=NY;j++)
                                   for (int k=0;k<=NZ;k++)
				   if ((Solid[0][j][k]>0) and (Solid[1][j][k]))
                                   {
                                   psi[Solid[0][j][k]]=psi[Solid[1][j][k]];
                                   rho_r[Solid[0][j][k]]=(psi[Solid[0][j][k]]*1.0+1.0)/2;
                                   rho_b[Solid[0][j][k]]=1.0-rho_r[Solid[0][j][k]];
                                   
                                   }
                        }
			else
			if (psi_xn==3)
			{
			for(int j=0;j<=NY;j++)
                                   for (int k=0;k<=NZ;k++)
				   if (Solid[0][j][k]>0)
                                   {
				   rand_double=(double(rand()%10000))/10000;
					if (rand_double<ini_Sat)
			        		psi_rand=1;
					else
			        		psi_rand=-1;
			
                                   rho_r[Solid[0][j][k]]=(psi_rand*1.0+1.0)/2;
                                   rho_b[Solid[0][j][k]]=1.0-rho_r[Solid[0][j][k]];
                                   psi[Solid[0][j][k]]=psi_rand;
                                   
                                   }
			}




                        
                      if ((psi_xp>0) and (rank==mpi_size-1))    
                        if (psi_xp==1)   
			{
                           for(int j=0;j<=NY;j++)
                                   for (int k=0;k<=NZ;k++)
				   if (Solid[nx_l-1][j][k]>0)
                                   {
                                   rho_r[Solid[nx_l-1][j][k]]=(Psi_local[nx_l-1][j][k]*1.0+1.0)/2;
                                   rho_b[Solid[nx_l-1][j][k]]=1.0-rho_r[Solid[nx_l-1][j][k]];
                                   psi[Solid[nx_l-1][j][k]]=Psi_local[nx_l-1][j][k];
                                   
                                   }
                           }  
			else
			if (psi_xp==2)
			{
                           for(int j=0;j<=NY;j++)
                                   for (int k=0;k<=NZ;k++)
				   if ((Solid[nx_l-1][j][k]>0) and (Solid[nx_l-2][j][k]>0))
                                   {
                              
                                   psi[Solid[nx_l-1][j][k]]=psi[Solid[nx_l-2][j][k]];
                                   rho_r[Solid[nx_l-1][j][k]]=(psi[Solid[nx_l-1][j][k]]*1.0+1.0)/2;
                                   rho_b[Solid[nx_l-1][j][k]]=1.0-rho_r[Solid[nx_l-1][j][k]];
                                   }
                           } 
			else 
			if (psi_xp==3)	
                        {
						
			for(int j=0;j<=NY;j++)
                             for (int k=0;k<=NZ;k++)
				   if (Solid[nx_l-1][j][k]>0)
                                   {
				rand_double=(double(rand()%10000))/10000;
					if (rand_double<ini_Sat)
			        		psi_rand=1;
					else
			        		psi_rand=-1;	
                                   rho_r[Solid[nx_l-1][j][k]]=(psi_rand*1.0+1.0)/2;
                                   rho_b[Solid[nx_l-1][j][k]]=1.0-rho_r[Solid[nx_l-1][j][k]];
                                   psi[Solid[nx_l-1][j][k]]=psi_rand;
                                   
                                   }


			}
                        
                         if (psi_yn==1)    
                           {
                           for(int i=0;i<nx_l;i++)
                                   for (int k=0;k<=NZ;k++)
				   if (Solid[i][0][k]>0)
                                   {
                                   rho_r[Solid[i][0][k]]=(Psi_local[i][0][k]*1.0+1.0)/2;
                                   rho_b[Solid[i][0][k]]=1.0-rho_r[Solid[i][0][k]];
                                   psi[Solid[i][0][k]]=Psi_local[i][0][k];
                                   
                                   }
                           }
			else 
			 if (psi_yn==2)    
                           {
                           for(int i=0;i<nx_l;i++)
                                   for (int k=0;k<=NZ;k++)
				   if ((Solid[i][0][k]>0) and (Solid[i][1][k]>0))
                                   {
                                   psi[Solid[i][0][k]]=psi[Solid[i][1][k]];
                                   rho_r[Solid[i][0][k]]=(psi[Solid[i][0][k]]*1.0+1.0)/2;
                                   rho_b[Solid[i][0][k]]=1.0-rho_r[Solid[i][0][k]];
                                   }
                           }
			else
			if (psi_yn==3)
			{
                           for(int i=0;i<nx_l;i++)
                                   for (int k=0;k<=NZ;k++)
				   if (Solid[i][0][k]>0)
                                   {

					rand_double=(double(rand()%10000))/10000;
					if (rand_double<ini_Sat)
			        		psi_rand=1;
					else
			        		psi_rand=-1;
                                   rho_r[Solid[i][0][k]]=(psi_rand*1.0+1.0)/2;
                                   rho_b[Solid[i][0][k]]=1.0-rho_r[Solid[i][0][k]];
                                   psi[Solid[i][0][k]]=psi_rand;
                                   
                                   }
                           }
 
                        
                        if (psi_yp==1)  
                           {
                           for(int i=0;i<nx_l;i++)
                                   for (int k=0;k<=NZ;k++)
					if (Solid[i][NY][k]>0)
                                   {
                                   rho_r[Solid[i][NY][k]]=(Psi_local[i][NY][k]*1.0+1.0)/2;
                                   rho_b[Solid[i][NY][k]]=1.0-rho_r[Solid[i][NY][k]];
                                   psi[Solid[i][NY][k]]=Psi_local[i][NY][k];
                                   
                                   }
                           }  
			else
			if (psi_yp==2)  
                           {
                           for(int i=0;i<nx_l;i++)
                                   for (int k=0;k<=NZ;k++)
					if ((Solid[i][NY][k]>0) and (Solid[i][NY-1][k]>0))
                                   {
                                   
                                   psi[Solid[i][NY][k]]=psi[Solid[i][NY-1][k]];
                                   rho_r[Solid[i][NY][k]]=(psi[Solid[i][NY][k]]*1.0+1.0)/2;
                                   rho_b[Solid[i][NY][k]]=1.0-rho_r[Solid[i][NY][k]];
                                   }
                           } 
				else
				if (psi_yp==3)  
				{
				for(int i=0;i<nx_l;i++)
                                   for (int k=0;k<=NZ;k++)
					if (Solid[i][NY][k]>0)
                                   {

					rand_double=(double(rand()%10000))/10000;
					if (rand_double<ini_Sat)
			        		psi_rand=1;
					else
			        		psi_rand=-1;
                                   rho_r[Solid[i][NY][k]]=(psi_rand*1.0+1.0)/2;
                                   rho_b[Solid[i][NY][k]]=1.0-rho_r[Solid[i][NY][k]];
                                   psi[Solid[i][NY][k]]=psi_rand;
                                   
                                   }
				} 


                        
                        if (psi_zn==1)    
                           {
                           for(int i=0;i<nx_l;i++)
                                   for (int j=0;j<=NY;j++)
					if (Solid[i][j][0]>0)
                                   {
                                   rho_r[Solid[i][j][0]]=(Psi_local[i][j][0]*1.0+1.0)/2;
                                   rho_b[Solid[i][j][0]]=1.0-rho_r[Solid[i][j][0]];
                                   psi[Solid[i][j][0]]=Psi_local[i][j][0];
                                   
                                   }
                           } 
			else
				 if (psi_zn==2)    
                           {
                           for(int i=0;i<nx_l;i++)
                                   for (int j=0;j<=NY;j++)
					if ((Solid[i][j][0]>0) and (Solid[i][j][1]>0))
                                   {
                                   
                                   psi[Solid[i][j][0]]=psi[Solid[i][j][1]];
                                   rho_r[Solid[i][j][0]]=(psi[Solid[i][j][0]]*1.0+1.0)/2;
                                   rho_b[Solid[i][j][0]]=1.0-rho_r[Solid[i][j][0]];
                                   }
                           }
			else 
			if (psi_zn==3)
				{
				for(int i=0;i<nx_l;i++)
                                   for (int j=0;j<=NY;j++)
					if (Solid[i][j][0]>0)
                                   {
					rand_double=(double(rand()%10000))/10000;
					if (rand_double<ini_Sat)
			        		psi_rand=1;
					else
			        		psi_rand=-1;
                                   rho_r[Solid[i][j][0]]=(psi_rand*1.0+1.0)/2;
                                   rho_b[Solid[i][j][0]]=1.0-rho_r[Solid[i][j][0]];
                                   psi[Solid[i][j][0]]=psi_rand;
                                   
                                   }
				}
                        


                    if (psi_zp==1)    
                           {
                           for(int i=0;i<nx_l;i++)
                                   for (int j=0;j<=NY;j++)
					if (Solid[i][j][NZ]>0)
                                   {
                                   rho_r[Solid[i][j][NZ]]=(Psi_local[i][j][NZ]*1.0+1.0)/2;
                                   rho_b[Solid[i][j][NZ]]=1.0-rho_r[Solid[i][j][NZ]];
                                   psi[Solid[i][j][NZ]]=Psi_local[i][j][NZ];
                                   
                                   }
                           } 
		else
			if (psi_zp==2)    
                           {
                           for(int i=0;i<nx_l;i++)
                                   for (int j=0;j<=NY;j++)
					if ((Solid[i][j][NZ]>0) and (Solid[i][j][NZ-1]>0))
                                   {
                                   
                                   psi[Solid[i][j][NZ]]=psi[Solid[i][j][NZ-1]];
                                   rho_r[Solid[i][j][NZ]]=(psi[Solid[i][j][NZ]]*1.0+1.0)/2;
                                   rho_b[Solid[i][j][NZ]]=1.0-rho_r[Solid[i][j][NZ]];
                                   }
                           } 
			else
			if (psi_zp==3)  
			   {
                           for(int i=0;i<nx_l;i++)
                                   for (int j=0;j<=NY;j++)
					if (Solid[i][j][NZ]>0)
                                   {
					rand_double=(double(rand()%10000))/10000;
					if (rand_double<ini_Sat)
			        		psi_rand=1;
					else
			        		psi_rand=-1;
                                   rho_r[Solid[i][j][NZ]]=(psi_rand*1.0+1.0)/2;
                                   rho_b[Solid[i][j][NZ]]=1.0-rho_r[Solid[i][j][NZ]];
                                   psi[Solid[i][j][NZ]]=psi_rand;
                                   
                                   }
                           }   
                        
                        
                }
	
                     
			
	//MPI_Barrier(MPI_COMM_WORLD); 

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
double m[19];


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


//Non-equilibrium boundary condition   M Matrix
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
	out<<"SCALARS sample_scalars double"<<endl;
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
	out<<"SCALARS sample_scalars double"<<endl;
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
	out<<"SCALARS sample_scalars double"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;

        for(int k=0;k<NZ0;k++)
      		for(int j=0; j<NY0; j++)
			for(int i=0;i<NX0;i++)
				out<<setprecision(preci)<<rbuf_psi[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]<<endl;

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
	out<<"SCALARS sample_scalars double"<<endl;
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
				out<<setprecision(preci)<<psi[Solid[i][j][k]]<<endl;
			else
				out<<setprecision(preci)<<psi_0<<endl;
			
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
				out<<setprecision(preci)<<rece[i*NY0*NZ0+j*NZ0+k]<<endl;	

		
		out.close();
		}
		
		
		if (rank==0)
			delete [] rece;
		
	}

	delete [] send;

}




void Comput_Perm_LOCAL(double* psi,double** u,double* Per_l,double* Per_g,int PerDIr)
{


	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();


	
	double Perm_l[3];
	double Perm_g[3];
	
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
		


	if (rank==mpi_size-1)
	{
		
	for (int j=0;j<=NY;j++)
		for (int k=0;k<=NZ;k++)
		if (Solid[nx_l-2][j][k]>0)
		{
			if (psi[Solid[nx_l-2][j][k]]>=rel_perm_psi)
			{
			Q_l[0]+=u[Solid[nx_l-2][j][k]][0];
			Q_l[1]+=u[Solid[nx_l-2][j][k]][1];
			Q_l[2]+=u[Solid[nx_l-2][j][k]][2];
			}
			if (psi[Solid[nx_l-2][j][k]]<=-rel_perm_psi)
			{
			Q_g[0]+=u[Solid[nx_l-2][j][k]][0];
			Q_g[1]+=u[Solid[nx_l-2][j][k]][1];
			Q_g[2]+=u[Solid[nx_l-2][j][k]][2];
			}
		}

	Perm_l[0]=Q_l[0]/((per_yp-per_yn+1)*(per_zp-per_zn+1))*(niu_l)/(gx+dp);
	Perm_l[1]=Q_l[1]/((per_xp-per_xn+1)*(per_zp-per_zn+1))*(niu_l)/(gy+dp);
	Perm_l[2]=Q_l[2]/((per_xp-per_xn+1)*(per_yp-per_yn+1))*(niu_l)/(gz+dp);


	Perm_g[0]=Q_g[0]/((per_yp-per_yn+1)*(per_zp-per_zn+1))*(niu_g)/(gx+dp);
	Perm_g[1]=Q_g[1]/((per_xp-per_xn+1)*(per_zp-per_zn+1))*(niu_g)/(gy+dp);
	Perm_g[2]=Q_g[2]/((per_xp-per_xn+1)*(per_yp-per_yn+1))*(niu_g)/(gz+dp);
	
	}
	
	//cout<<Q_l[0]<<"   "<<Q_g[0]<<endl;

	MPI_Barrier(MPI_COMM_WORLD);

	MPI_Bcast(Perm_l,3,MPI_DOUBLE,mpi_size-1,MPI_COMM_WORLD);
	MPI_Bcast(Perm_g,3,MPI_DOUBLE,mpi_size-1,MPI_COMM_WORLD);

	
		Per_l[0]=Perm_l[0];Per_g[0]=Perm_g[0];
		Per_l[1]=Perm_l[1];Per_g[1]=Perm_g[1];
		Per_l[2]=Perm_l[2];Per_g[2]=Perm_g[2];
	
	

}




double Comput_Perm(double* psi,double** u,double* Per_l,double* Per_g,int PerDIr, int* SupInv)
{


	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	int nx_g[mpi_size];
	int disp[mpi_size];
	int si,sj,sm;
	
	MPI_Gather(&nx_l,1,MPI_INT,nx_g,1,MPI_INT,0,MPI_COMM_WORLD);
	
	
	if (rank==0)
		{
		disp[0]=0;
	
		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
		
		}

	MPI_Bcast(disp,mpi_size,MPI_INT,0,MPI_COMM_WORLD);

//	if (rank==2)
//		for (int i=0;i<mpi_size;i++)
//			cout<<disp[i]<<"   "<<i<<endl;


//	if (rank==1)
//		cout<<per_zp<<" "<<per_zn<<endl;

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
	        if (psi[i]>=rel_perm_psi)
		        {
		                Q_l[0]+=u[i][0];
		                Q_l[1]+=u[i][1];
		                Q_l[2]+=u[i][2];
		        }
		if (psi[i]<=-rel_perm_psi)
		        {
                                Q_g[0]+=u[i][0];
                                Q_g[1]+=u[i][1];
                                Q_g[2]+=u[i][2];
		        
		        
		        }
		}


	}
	else
	for (int i=1;i<=Count;i++)
		{
	        if (psi[i]>=rel_perm_psi)
		        {
		                Q_l[0]+=u[i][0];
		                Q_l[1]+=u[i][1];
		                Q_l[2]+=u[i][2];
		        }
		if (psi[i]<=-rel_perm_psi)
		        {
                                Q_g[0]+=u[i][0];
                                Q_g[1]+=u[i][1];
                                Q_g[2]+=u[i][2];
		        
		        
		        }
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

		Perm_l[0]=Q_l[0]/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1))*(niu_l)/(gx+dp);
		Perm_l[1]=Q_l[1]/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1))*(niu_l)/(gy+dp);
		Perm_l[2]=Q_l[2]/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1))*(niu_l)/(gz+dp);

		Perm_g[0]=Q_g[0]/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1))*(niu_g)/(gx+dp);
		Perm_g[1]=Q_g[1]/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1))*(niu_g)/(gy+dp);
		Perm_g[2]=Q_g[2]/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1))*(niu_g)/(gz+dp);


		vxl=(Q_l[0]+Q_g[0])/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1));
		vyl=(Q_l[1]+Q_g[1])/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1));
		vzl=(Q_l[2]+Q_g[2])/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1));
		Capillary=sqrt(vxl*vxl+vyl*vyl+vzl*vzl)/porosity*niu_l/CapA;   

		u_ave=sqrt((vxl/porosity)*(vxl/porosity)+(vyl/porosity)*(vyl/porosity)+(vzl/porosity)*(vzl/porosity));

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



double Comput_Saturation(double* psi,int*** Solid,int* SupInv)
{
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	double S_l,S_g;
	
	int nx_g[mpi_size];
	int disp[mpi_size];
	int si,sj,sm;
	
	MPI_Gather(&nx_l,1,MPI_INT,nx_g,1,MPI_INT,0,MPI_COMM_WORLD);
	
	
	if (rank==0)
		{
		disp[0]=0;
	
		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
		
		}

	MPI_Bcast(disp,mpi_size,MPI_INT,0,MPI_COMM_WORLD);
	
double *rbuf_l,*rbuf_g;

	rbuf_l=new double[mpi_size];
	rbuf_g=new double[mpi_size];

	S_l=0;S_g=0;

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
	
			if (psi[i]>=0) 
			S_l+=1;
			else
			S_g+=1;
		}
		
		
	}
	else
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
	
	
	//cout<<porosity<<endl;
	S_l=S_l/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1)*porosity);
	S_g=S_g/((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1)*porosity);
	}

	
	delete [] rbuf_l;
	delete [] rbuf_g;
	
	
	return (S_l);
			

}




void Backup_init(double* rho, double** u, double** f,double* psi,double* rho_r, double* rho_b, double* rhor, double* rhob)
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
	
	ostringstream name5;
	name5<<pfix<<"LBM_checkpoint_rhor_"<<mode_backup_ini<<"."<<rank<<".bin_input";
 	ostringstream name4;
	name4<<pfix<<"LBM_checkpoint_f_"<<mode_backup_ini<<"."<<rank<<".bin_input";
	ostringstream name3;
	name3<<pfix<<"LBM_checkpoint_psi_"<<mode_backup_ini<<"."<<rank<<".bin_input";
 	ostringstream name2;
	name2<<pfix<<"LBM_checkpoint_velocity_"<<mode_backup_ini<<"."<<rank<<".bin_input";
	ostringstream name;
	name<<pfix<<"LBM_checkpoint_rho_"<<mode_backup_ini<<"."<<rank<<".bin_input";
	ostringstream name6;
	name6<<pfix<<"LBM_checkpoint_rhob_"<<mode_backup_ini<<"."<<rank<<".bin_input";
	 
	fstream fin;
	fin.open(name.str().c_str(),ios::in);
	if (fin.fail())
	        {
	        cout<<"\n file open error on" << name.str().c_str()<<endl;
	        exit(-1);
	        }
	fin.read((char *)(&rho[0]), sizeof(double)*(Count+1));
 
       fin.close();
       
   
	fin.open(name2.str().c_str(),ios::in);
	if (fin.fail())
	        {
	        cout<<"\n file open error on" << name2.str().c_str()<<endl;
	        exit(-1);
	        }
	fin.read((char *)(&u[0][0]), sizeof(double)*(Count+1)*3);
      
       fin.close();
       
       fin.open(name3.str().c_str(),ios::in);
       if (fin.fail())
	        {
	        cout<<"\n file open error on" << name3.str().c_str()<<endl;
	        exit(-1);
	        }
       fin.read((char *)(&psi[0]), sizeof(double)*(Count+1));
     
       fin.close();
       
       fin.open(name4.str().c_str(),ios::in);
       if (fin.fail())
	        {
	        cout<<"\n file open error on" << name4.str().c_str()<<endl;
	        exit(-1);
	        }
	fin.read((char *)(&f[0][0]), sizeof(double)*(Count+1)*19);
      
       fin.close();
       
       fin.open(name5.str().c_str(),ios::in);
       if (fin.fail())
	        {
	        cout<<"\n file open error on" << name5.str().c_str()<<endl;
	        exit(-1);
	        }
	fin.read((char *)(&rho_r[0]), sizeof(double)*(Count+1));
        
       fin.close();

       fin.open(name6.str().c_str(),ios::in);
       if (fin.fail())
	        {
	        cout<<"\n file open error on" << name6.str().c_str()<<endl;
	        exit(-1);
	        }
	fin.read((char *)(&rho_b[0]), sizeof(double)*(Count+1));
        
       fin.close();
       
       //cout<<psi[3]<<endl;
       
       
        	for(int i=1;i<=Count;i++)
        	{
        	       
			rhor[i]=0;
			rhob[i]=0;
		}
  
       
      
	
	
	
		       if (stab==1)
				{gxs=0;gys=0;gzs=0;}
			else
				{gxs=gx;gys=gy;gzs=gz;}
			
	
	 	
}




void Backup(int m,double* rho,double* psi, double** u, double** f, double* rho_r, double* rho_b)
{

        int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	
	ostringstream name;
	name<<outputfile<<"LBM_checkpoint_velocity_"<<m<<"."<<rank<<".bin_input";
	ofstream out;
	out.open(name.str().c_str());
	
        out.write((char *)(&u[0][0]),sizeof(double)*(Count+1)*3);	
			
	out.close();
	

	
	
	ostringstream name2;
	name2<<outputfile<<"LBM_checkpoint_rho_"<<m<<"."<<rank<<".bin_input";
	ofstream out2;
	out2.open(name2.str().c_str());
	

	out2.write((char *)(&rho[0]), sizeof(double)*(Count+1));	
			
	out2.close();
	

	
	ostringstream name3;
	name3<<outputfile<<"LBM_checkpoint_psi_"<<m<<"."<<rank<<".bin_input";
	//ofstream out;
	out.open(name3.str().c_str());
	
	
	out.write((char *)(&psi[0]), sizeof(double)*(Count+1));	
	
	out.close();
	
	
	ostringstream name4;
	name4<<outputfile<<"LBM_checkpoint_f_"<<m<<"."<<rank<<".bin_input";
	//ofstream out;
	out.open(name4.str().c_str());

        out.write((char *)(&f[0][0]), sizeof(double)*(Count+1)*19);      
        out.close();
        
	
        ostringstream name5;
	name5<<outputfile<<"LBM_checkpoint_rhor_"<<m<<"."<<rank<<".bin_input";
	//ofstream out2;
	out2.open(name5.str().c_str());
	

	out2.write((char *)(&rho_r[0]), sizeof(double)*(Count+1));	
			
	out2.close();

	ostringstream name6;
	name6<<outputfile<<"LBM_checkpoint_rhob_"<<m<<"."<<rank<<".bin_input";
	//ofstream out2;
	out2.open(name6.str().c_str());
	

	out2.write((char *)(&rho_b[0]), sizeof(double)*(Count+1));	
			
	out2.close();

}



