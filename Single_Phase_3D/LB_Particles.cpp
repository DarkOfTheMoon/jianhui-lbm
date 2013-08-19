//Use INPUT_LB_RW.dat
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string.h>
#include<math.h> 
# include "mpi.h"

//======PARMETIS===============
//#include<parmetis.h>  
//#include "/home/jy810/LBM/source/CODE/jianhui-lbm/parmetis-4.0.2/include/parmetis.h"
//mpic++ SINGLE_PHASE_MPI_Software_Spars_3DPartition.cpp -lparmetis -lmetis -o paratest

//=============================


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
double pre_u_ave=1.0;


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

void Partition_Solid(int***);

void periodic_streaming(double** ,double** ,int* ,int***,int*, int*,double*, double**);

void periodic_streaming_MR(double** ,double** ,int* ,int*** ,int* ,int* ,double* ,double** );

void standard_bounceback_boundary(int,double**);

void collision(double*,double** ,double** ,double** , int* ,int***,int*, int*);

void collision_nnf(double*,double** ,double** ,double** , int* ,int***,int*, int*);

void comput_macro_variables( double* ,double**,double** ,double** ,double**  ,int* ,int***);

double Error(double** ,double** ,double*, double*);

void boundary_velocity(int,double,int , double ,int ,double ,int ,double ,int ,double ,int ,double ,double** ,double**,double* ,double** ,int*** );

void boundary_pressure(int ,double ,int , double ,int ,double ,int ,double ,int ,double ,int ,double ,double** ,double**,double** ,double* ,int*** );

void output_velocity(int ,double* ,double** ,int ,int ,int ,int ,int*** );

void output_velocity_compact(int m,double* rho,double** u,int MirX,int MirY,int MirZ,int mir,int*** Solid);	

void output_density(int ,double* ,int ,int ,int ,int ,int*** );	

void Geometry(int*** );	

void output_velocity_b(int ,double* ,double** ,int ,int ,int ,int ,int*** );

void output_density_b(int ,double* ,int ,int ,int ,int ,int*** );	

void Geometry_Par(int*** );

void Backup(int ,double* ,double**, double**);

double Comput_Perm(double** u,double*,int,int*);

double S[19];

void Comput_MI(double[19][19], double[19][19]);

int inverse(mat &a);

inline double feq(int,double, double[3]);

void Suppliment(int*,int***);

void Backup_init(double* rho, double** u, double** f, char[128], char[128],char[128]);

void Parallelize_Geometry();

void Partition_Solid_SELF(int***);

void Comput_Grop_Perm(double** ,double* ,int ,int* );

void output_velocity_for_solute(int ,double* ,double** ,int ,int,int,int,int*** );

void Comput_Perm_LOCAL(double** ,double* ,int);


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
int Zoom,vel_sol;


int wr_per,pre_xp,pre_xn,pre_yp,pre_yn,pre_zp,pre_zn,fre_backup,lattice_v;
int vel_xp,vel_xn,vel_yp,vel_yn,vel_zp,vel_zn,Sub_BC,Out_Mode,mode_backup_ini;
double in_vis,p_xp,p_xn,p_yp,p_yn,p_zp,p_zn,dx_input,dt_input;
double inivx,inivy,inivz,v_xp,v_xn,v_yp,v_yn,v_zp,v_zn;
double error_perm;
int par_per_x,par_per_y,par_per_z,per_xp,per_xn,per_yp,per_yn,per_zp,per_zn;


int n_gperm1,n_gperm2,n_gperm3,n_gperm4,gperm,loc_perm;
int size_gperm1,size_gperm2,size_gperm3,size_gperm4;
int c0_gperm1,c1_gperm1,c2_gperm1,c0_gperm2,c1_gperm2,c2_gperm2,c0_gperm3,c1_gperm3,c2_gperm3,c0_gperm4,c1_gperm4,c2_gperm4;

char outputfile[128]="./";
int NCHAR=128;
	char     filename[128], dummy[128+1], backup_rho[128], backup_velocity[128],backup_f[128];
	int      dummyInt;
	
int*** Solid;	
int*** Solid2;
double** u;



char pfix[128];
int decbin;
double Permia_LOCAL[3]={0.0,0.0,0.0};


//==============hybrid lb dispersion===============
int updatesss=0;

//=================================================

//============MPI==TRANSFER==INITIALIZATION===============

double* sendl;
double* sendr;

	double* recvl;
	double* recvr;
	
	
int* sumss;             //fluid nodes number of every partition            [procn+1]	
int* bufinfo;                                   //number of nodes of the partirion that need communiate with current processors
                                                        // 0= no contact with current processor, >0 number of nodes that need to communicate with current processor.
                                                         //start from 1 size procn+1       designed for 19 components for f function transfer
int* com_ind;                          //commu nodes indexs (partition no.)   com_ind[0,start from 0] size new int[com_n]
int* com_loc;                           //mpi commu different nodes starting locations in buffet  arrays                                                         
int** nei_loc;                           //index of 18 neibourghs, nei_loc[3][0] is the first neighbour (e[1][]) of node 3
//-------------------
int* coor;      //start from 1, int [sumss[procind]+1];
//------------------

double** bufsend;                       //send buffet for f functions bufsend[comm_index][number]
double** bufrecv;                       //recv buffet for f functions
//-----------------
int** buflocsend;                       //exchange commu info, used to locate the data after when it is received from MPI communications
int** buflocrecv;                       //two digital combination index*19+ls ls is the direction of 18 vectors

int* bclx;
int* bcly;
int* bclz;
int* bcrx;
int* bcry;
int* bcrz;
int bclxn,bclyn,bclzn,bcrxn,bcryn,bcrzn;
 int com_n;

//==================================================

//============NNF_module==============================
double* niu_vn;
int nnf_module;
double nnf_m,nnf_n;
//==================================================



int main(int argc , char *argv [])
{	

MPI :: Init (argc , argv );
MPI_Status status ;

double start , finish,remain,elaps;

int rank = MPI :: COMM_WORLD . Get_rank ();
int para_size=MPI :: COMM_WORLD . Get_size ();
//**********************
int mpi_size=para_size;
//***********************


int dif,ts,th,tm;
int tse,the,tme;
double st1,st2; 

	string pfix2;

        strcpy(pfix,"./");
        if (argc>2)
		{
		//strcpy(pfix2,argv[2]);
		//updatesss=argv[2];
		updatesss=atoi(argv[2]);
                //strcpy(pfix,argv[2])
		}	
	
	
	//cout<<updatesss<<"	&&&&&&	"<<endl;
        

	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();



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
	fin >> loc_perm;					fin.getline(dummy, NCHAR);
	fin >> par_per_x >> par_per_y >>par_per_z;	fin.getline(dummy, NCHAR);
	fin >> per_xp >> per_xn;			fin.getline(dummy, NCHAR);
	fin >> per_yp >> per_yn;			fin.getline(dummy, NCHAR);
	fin >> per_zp >> per_zn;			fin.getline(dummy, NCHAR);
							fin.getline(dummy, NCHAR);
	fin >> fre_backup;                        	fin.getline(dummy, NCHAR);
	fin >>mode_backup_ini;                		fin.getline(dummy, NCHAR);

	fin >> vel_sol;                                	fin.getline(dummy, NCHAR);
	

	fin.getline(dummy, NCHAR);
	fin >> st1;					fin.getline(dummy, NCHAR);
	fin >> st2;	fin.getline(dummy, NCHAR);
	fin >> size_gperm1>>size_gperm2>>size_gperm3>>size_gperm4; fin.getline(dummy, NCHAR);
	fin >> c0_gperm1>>c1_gperm1>>c2_gperm1; 	fin.getline(dummy, NCHAR);	
	fin >> c0_gperm2>>c1_gperm2>>c2_gperm2; 	fin.getline(dummy, NCHAR);
	fin >> c0_gperm3>>c1_gperm3>>c2_gperm3; 	fin.getline(dummy, NCHAR);
	fin >> c0_gperm4>>c1_gperm4>>c2_gperm4; 	fin.getline(dummy, NCHAR);
							        fin.getline(dummy, NCHAR);
	fin >> decbin;                                fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin >>nnf_module;                        fin.getline(dummy, NCHAR);
	fin >> nnf_m >> nnf_n;                fin.getline(dummy, NCHAR);

fin.close();
	
	//cout<<nnf_n<<"    asdfa "<<endl;
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
	MPI_Bcast(&dt_input,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&vel_sol,1,MPI_INT,0,MPI_COMM_WORLD);
	
	MPI_Bcast(&par_per_x,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&par_per_y,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&par_per_z,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&per_yp,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&per_xp,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&per_yn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&per_zp,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&per_zn,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&per_xn,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Bcast(&gperm,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&n_gperm1,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&n_gperm2,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&n_gperm3,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&n_gperm4,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&size_gperm1,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&size_gperm2,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&size_gperm3,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&size_gperm4,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&c0_gperm1,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&c1_gperm1,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&c2_gperm1,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&c0_gperm2,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&c1_gperm2,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&c2_gperm2,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&c0_gperm3,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&c1_gperm3,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&c2_gperm3,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&c0_gperm4,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&c1_gperm4,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&c2_gperm4,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&loc_perm,1,MPI_INT,0,MPI_COMM_WORLD);

	MPI_Bcast(&nnf_module,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&nnf_n,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&nnf_m,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
	MPI_Bcast(&st1,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&st2,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	Par_Geo=0;

if (mirX==1)
	NX=NX*2+1;
if (mirY==1)
	NY=NY*2+1;

	
      
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


	double* Permia;
	double* rho;
	
	double**f;
	double**F;
	double**u0;
	int* SupInv;

	
	int*  Sl;
	int*  Sr;

	
	



        Parallelize_Geometry();
        
	
        
        MPI_Barrier(MPI_COMM_WORLD);
        
	

	//***************************************************
	//WARRING: SPARSE MATRIX STARTS FROM INDEX 1 NOT 0!!!
	//***************************************************

	Permia = new double[3];
	rho = new double[Count+1];
	
	if (nnf_module==1)
	        niu_vn = new double[Count+1];
	
	
	/*
	u = new double*[Count+1];
	u[0] = new double[(Count+1)*3];
	for (int i=1;i<=Count;i++) 
	        u[i] = u[i-1]+3;
	*/

	
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
	
	
	//SupInv = new int[Count+1];


	Comput_MI(M,MI);
	
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	if (Out_Mode==1)
		Geometry_Par(Solid);

	if ((freVe>=0) or (freDe>=0))
	{
	
	
		Geometry(Solid);
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
strcat(FileName2,"Permeability_error_local_Perm.txt");
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


 
	for(n=0;n<=n_max;n++)
	{
	
	if (nnf_module==0)
	        collision(rho,u,f,F,SupInv,Solid,Sl,Sr);
	else
	        collision_nnf(rho,u,f,F,SupInv,Solid,Sl,Sr);//cout<<"markcollision"<<endl;
//cout<<"@@@@@@@@@@@   "<<n<<endl;
	
	if ((1-pre_xp)*(1-pre_xn)*(1-pre_yp)*(1-pre_yn)*(1-pre_zp)*(1-pre_zn)==0)
		boundary_pressure(pre_xp,p_xp,pre_xn,p_xn,pre_yp,p_yp,pre_yn,p_yn,pre_zp,p_zp,pre_zn,p_zn,f,F,u,rho,Solid);

	
	if ((1-vel_xp)*(1-vel_xn)*(1-vel_yp)*(1-vel_yn)*(1-vel_zp)*(1-vel_zn)==0)
		boundary_velocity(vel_xp,v_xp,vel_xn,v_xn,vel_yp,v_yp,vel_yn,v_yn,vel_zp,v_zp,vel_zn,v_zn,f,F,rho,u,Solid);

//cout<<"@@@@@@@@@@@   "<<n<<endl;
  		comput_macro_variables(rho,u,u0,f,F,SupInv,Solid); 


 
  		
	
	if(n%freRe==0)
		{       
			
			

			pre_u_ave=u_ave;
			 error=Error(u,u0,&u_max,&u_ave);if (u_max>=10.0)	U_max_ref+=1;
			error_perm=Comput_Perm(u,Permia,PerDir,SupInv); 
			
			
			MPI_Bcast(&u_ave,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			MPI_Bcast(&error,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
			
			//cout<<error<<"                ccccccccccccc      "<<u_ave<<"        "<<rank<<endl;
			
			
			
			 if (rank==0)
			{ 
			    
			//ofstream fin(FileName,ios::app);          
			finish = MPI_Wtime();
			
			remain=(n_max-n)*((finish-start)/n);
			
			th=int(remain/3600);
			tm=int((remain-th*3600)/60);
			ts=int(remain-(th*3600+tm*60));

			elaps=finish-start;
			the=int(elaps/3600);
			tme=int((elaps-the*3600)/60);
			tse=int(elaps-(the*3600+tme*60));



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
			cout<<"The average velocity is: "<<setprecision(6)<<u_ave<<endl;
		//===============================================================================================
			cout<<"The max relative error of uv is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			cout<<"Elapsed time is "<<the<<"h"<<tme<<"m"<<tse<<"s"<<endl;
			cout<<"The expected completion time is "<<th<<"h"<<tm<<"m"<<ts<<"s"<<endl;
			cout<<endl;
			}
			
			//if ((Out_Mode==1) and (abs((u_ave2-u_ave)/(u_ave2))<1e-7))
			//cout<<st1<<"	"<<st2<<"	"<<loc_perm<<" "<<error<<" "<<abs((pre_u_ave-u_ave)/(pre_u_ave))<<endl; 
			if ((loc_perm==1) and (error<st2) and (abs((pre_u_ave-u_ave)/(pre_u_ave))<st1))
				{
				//cout<<"@@@@@@@@@@@@@@@@@@@@@"<<endl;
				output_velocity_compact(n,rho,u,mirX,mirY,mirZ,mir,Solid);
				n=n_max+1;				
				}

			
			if ((freDe>0) and (n%freDe==0))
					output_density(n,rho,mirX,mirY,mirZ,mir,Solid);
				
			
			if ((freVe>0) and (n%freVe==0))
					output_velocity(n,rho,u,mirX,mirY,mirZ,mir,Solid);
				
				
			
			 
			if(error!=error) {cout<<"PROGRAM STOP"<<endl;break;};
			if(U_max_ref>=5) {cout<<"PROGRAM STOP DUE TO HIGH VELOCITY"<<endl;break;}
		}	
	}

	if (fre_backup>=0)
			        Backup(n_max,rho,u,f);
	
/*
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
    int nx=NX+1;
    int ny=NY+1;
    int nz=NZ+1;
    
    int ii,jj,kk,pore2;
    int ip,jp,kp;
    
   if (par_per_x==0)
		{per_xn=0;per_xp=NX;}
	if (par_per_y==0)
		{per_yn=0;per_yp=NY;}
	if (par_per_z==0)
		{per_zn=0;per_zp=NZ;}
                                
    int procind=MPI :: COMM_WORLD . Get_rank ()+1;                               //currentprocessor index, start from 1
    int procn=MPI :: COMM_WORLD . Get_size ();                                  //total processor number
    int neib=0;
   
            bufinfo=new int[procn+1];
            for (int i=0;i<=procn;i++)
                    bufinfo[i]=0;   
            
            
    com_n=0;                                                                    //mpi commu numbers     number of neighbour partitions which need communication
    int tmpint;
    int proc_com[procn+1];                                              //index convert proc index---->commu index in current processor
    for (int i=0;i<=procn;i++)
        proc_com[i]=0;
        
      //-------------------
      int* sumtmp;
	int upx,upy,upz;
	double updoux,updouy,updouz;
	int*** Solid3;
      //-------------------  
        
      Solid = new int**[nx];
      Solid2 = new int**[nx];
	Solid3 = new int**[nx];
	
	
	for (int i=0;i<nx;i++)				///*********
		Solid[i]=new int*[ny],Solid2[i]=new int*[ny],Solid3[i]=new int*[ny];

	Solid[0][0]=new int[nx*ny*nz],Solid2[0][0]=new int[nx*ny*nz],Solid3[0][0]=new int[nx*ny*nz];

	
 	for (int i=1;i<ny;i++)
               Solid[0][i]=Solid[0][i-1]+nz,Solid2[0][i]=Solid2[0][i-1]+nz,Solid3[0][i]=Solid3[0][i-1]+nz;
       
       for (int i=1;i<nx;i++)
       {
               Solid[i][0]=Solid[i-1][0]+ny*nz,Solid2[i][0]=Solid2[i-1][0]+ny*nz,Solid3[i][0]=Solid3[i-1][0]+ny*nz;
               for (int j=1;j<ny;j++)
                       Solid[i][j]=Solid[i][j-1]+nz,Solid2[i][j]=Solid2[i][j-1]+nz,Solid3[i][j]=Solid3[i][j-1]+nz;
       }	
	
      
      for(int k=0 ; k<=NZ ; k++)
	for(int j=0 ; j<=NY ; j++)
	for(int i=0 ; i<=NX ; i++)
		Solid[i][j][k]=0,Solid2[i][j][k]=0,Solid3[i][j][k]=0;
      
   porosity=0.0;
   
 
    int* recv_solid;
    int pore;
   int pre_sum=0;
   
   
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
	
	fin.open(filename);
	for(int k=0 ; k<=NZ ; k++)
	for(int j=0 ; j<=NY ; j++)
	for(int i=0 ; i<=NX ; i++)
	
	{
		
			fin >> pore;
			
	
		
		Solid3[i][j][k]=pore;
		
	}
	fin.close();

	}
	else
	{
	fstream fin;
	fin.open(filename,ios::in);
	if (fin.fail())
	        {
	        cout<<"\n file open error on " << filename<<endl;
	        exit(-1);
	        }
	
	fin.read((char *)(Solid3[0][0]), sizeof(int)*(NX+1)*(NY+1)*(NZ+1));
	//fin.read((char *)(Solid3[0][0]), sizeof(bool)*(NX+1)*(NY+1)*(NZ+1));
	
	fin.close();
	}

	//***********************update part 1**********************************
	
	for(int k=0 ; k<=NZ ; k++)
	for(int j=0 ; j<=NY ; j++)
	for(int i=0 ; i<=NX ; i++)
	{
		Solid[i][j][k]=Solid3[i][j][k];
		if (Solid3[i][j][k]==0)
		        pre_sum++;
	}
	
	cout<<pre_sum<<"           %%%%%     previous velocity number"<<endl;

	if (updatesss==1)
	{
	FILE *ftest;
	ftest = fopen("update.txt", "r");

	if(ftest == NULL)
		{
		cout << "\n The pore geometry file ( update.txt ) does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(-1);
		}
		fclose(ftest);
	
	fin.open("update.txt");
	while (!fin.eof())
	{
	fin >> tmpint >> upx >> upy >> upz;fin.getline(dummy, NCHAR,'\n');
	cout<<upx<<"	"<<upy<<"	&&&&&&&&&&&&& "<<rank<<endl;
		Solid[upx][upy][upz] = tmpint;
	}
	
	fin.close();
	}
	//********************************************************************

}
  
        MPI_Bcast(&pre_sum,1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(Solid[0][0],(NX+1)*(NY+1)*(NZ+1),MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(Solid3[0][0],(NX+1)*(NY+1)*(NZ+1),MPI_INT,0,MPI_COMM_WORLD);

        for (int k=0;k<=NZ;k++)
		for (int j=0;j<=NY;j++)
        		for (int i=0;i<=NX;i++)
                
			{
	                
			if ((Solid[i][j][k]==0) and (i>=per_xn) and (i<=per_xp) and (j>=per_yn) and (j<=per_yp) and (k>=per_zn) and (k<=per_zp))
			        porosity+=1.0;
			}
	        sumss=new int [procn+1];
	for (int i=0;i<=procn;i++)
	        sumss[i]=0;
	
	porosity=porosity/(double)((per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1));

	
	
	Partition_Solid_SELF(Solid);




	//cout<<porosity<<endl;
	
	for(int k=0 ; k<nz ; k++)			
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)			
	        if (Solid[i][j][k]>0)
	        {
	                sumss[Solid[i][j][k]]++;
	                Solid2[i][j][k]=sumss[Solid[i][j][k]];
	                
	
	
	                //=======calculate neibough numbers==========
	                if (Solid[i][j][k]==procind)
	                for (int ls=1;ls<19;ls++)
	                {
	                        ii=i+e[ls][0];
	                                if (ii<0)
	                                        ii=nx-1;
	                                if (ii>=nx)
	                                        ii=0;
	                                
	                                        
	                        jj=j+e[ls][1];
	                                if (jj<0)
	                                        jj=ny-1;
	                                if (jj>=ny)
	                                        jj=0;
	                                
	                        kk=k+e[ls][2];
	                                       if (kk<0)
	                                        kk=nz-1;
	                                        if (kk>=nz)
	                                        kk=0; 
	                        
	                                if ((Solid[ii][jj][kk]>0) and (Solid[ii][jj][kk]!=procind))
	                                        {
	                                                bufinfo[Solid[ii][jj][kk]]++;
	                                        //cout<<procind<<"        "<<Solid[ii][jj][kk]<<endl;
	                                        }
	                                
	                                        
	                
	                }
	        //=======================================================
	        }
	        else
	                Solid2[i][j][k]=0;
	       
	        //======coordinate of nodes===========
	        coor = new int[sumss[procind]+1];
	                for(int k=0 ; k<nz ; k++)			
	                for(int j=0 ; j<ny ; j++)
	                for(int i=0 ; i<nx ; i++)
	                        if (Solid[i][j][k]==procind)
	                                coor[Solid2[i][j][k]]=i*ny*nz+j*nz+k;
	        //=============================
	        
	        
	        for (int i=1;i<=procn;i++)
	        {
	                //cout<<bufinfo[i]<<endl;
	                if (bufinfo[i]>0)
	                        com_n++;
	        }
	        com_ind=new int[com_n];
	                
	        tmpint=0;
	        for (int i=1;i<=procn;i++)
	                if (bufinfo[i]>0)
	                {com_ind[tmpint]=i;proc_com[i]=tmpint;tmpint++;}
	                
	
	
	bufsend = new double* [com_n];
	        for (int i=0;i<com_n;i++)
	                bufsend[i] = new double[bufinfo[com_ind[i]]];
	        
	bufrecv = new double* [com_n];
	        for (int i=0;i<com_n;i++)
	                bufrecv[i] = new double[bufinfo[com_ind[i]]];
	        
	        
	nei_loc= new int*[sumss[procind]+1];
	nei_loc[0] =new int[(sumss[procind]+1)*19];
	        for (int i=1;i<=sumss[procind];i++)
	                nei_loc[i] = nei_loc[i-1]+19;
	  
	        buflocsend = new int*[com_n];
	        for (int i=0;i<com_n;i++)
	                buflocsend[i] = new int[bufinfo[com_ind[i]]];
	        
	        buflocrecv = new int*[com_n];
	        for (int i=0;i<com_n;i++)
	                buflocrecv[i] = new int[bufinfo[com_ind[i]]];
	        
	        sumtmp = new int[com_n];
	                for(int i=0;i<com_n;i++)
	                        sumtmp[i]=0;
	                
	                
	
	        for (int ci=1;ci<=sumss[procind];ci++)
	                {
	                      ii=(int)(coor[ci]/(ny*nz));
	                      jj=(int)((coor[ci]%(ny*nz))/nz);
	                      kk=(int)(coor[ci]%nz);
	                      
	                      for (int mi=0; mi<19; mi++)
			{
			
			
			        ip=ii+e[mi][0];if (ip<0) {ip=nx-1;};if (ip>=nx) {ip=0;};
			        jp=jj+e[mi][1];if (jp<0) {jp=ny-1;}; if (jp>=ny) {jp=0;};
			        kp=kk+e[mi][2];if (kp<0) {kp=nz-1;}; if (kp>=nz) {kp=0;};
			        
			        if (Solid[ip][jp][kp]==procind)
			                nei_loc[ci][mi]=Solid2[ip][jp][kp];
			        else
			                if (Solid[ip][jp][kp]==0)
			                        nei_loc[ci][mi]=0;
			                else
			                {
			                       // nei_loc[ci][mi]=-Solid[ip][jp][kp];
			                       nei_loc[ci][mi]=-proc_com[Solid[ip][jp][kp]]-1;
			                        buflocsend[proc_com[Solid[ip][jp][kp]]][sumtmp[proc_com[Solid[ip][jp][kp]]]]=Solid2[ip][jp][kp]*19+mi;
			                      
			                        sumtmp[proc_com[Solid[ip][jp][kp]]]++;
			                       
			                }
			                
			}
			                
			                
	                    
	                }
	                
	                
	     
	                                               
	       MPI_Status status[com_n*2] ;
	       MPI_Request request[com_n*2];         
	       int mpi_test=procn;
	       
	       for (int i=0;i<com_n;i++)
	       
	               {
	                       MPI_Isend(buflocsend[i],bufinfo[com_ind[i]], MPI_INT, com_ind[i]-1, (procind-1)*procn+com_ind[i]-1, MPI_COMM_WORLD,&request[2*i]);
	                       
	                       MPI_Irecv(buflocrecv[i],bufinfo[com_ind[i]], MPI_INT, com_ind[i]-1, (com_ind[i]-1)*procn+procind-1, MPI_COMM_WORLD,&request[2*i+1]);		
	               }
	               
	               
	               
	      //-------------BCs--------------------------------------------------------         
	               bclxn=0;bclyn=0;bclzn=0;bcrxn=0;bcryn=0;bcrzn=0;
	                
	                 for(int k=0; k<nz ; k++)			
	                 for(int j=0 ; j<ny ; j++)
	                 {
	                         if (Solid[0][j][k]==procind)
	                                 bclxn++;
	                         if (Solid[nx-1][j][k]==procind)
	                                 bcrxn++;
	                 }
	                 if (bclxn>0)
	                         bclx=new int[bclxn];
	                 if (bcrxn>0)
	                         bcrx= new int[bcrxn];
	                 bclxn=0;bcrxn=0;
	                 for(int k=0; k<nz ; k++)			
	                 for(int j=0 ; j<ny ; j++)
	                 {
	                         if (Solid[0][j][k]==procind)
	                                 bclx[bclxn]=Solid2[0][j][k],bclxn++;
	                         
	                         if (Solid[nx-1][j][k]==procind)
	                                 bcrx[bcrxn]=Solid2[nx-1][j][k],bcrxn++;
	                 }
	                 
	                 
	                 
	               for(int k=0 ; k<nz ; k++)			
	               for(int i=0 ; i<nx ; i++)  
	                 {
	                         if (Solid[i][0][k]==procind)
	                                 bclyn++;
	                         if (Solid[i][ny-1][k]==procind)
	                                 bcryn++;
	                 }
	                 if (bclyn>0)
	                         bcly=new int[bclyn];
	                 if (bcryn>0)
	                         bcry= new int[bcryn];
	                 bclyn=0;bcryn=0;
	                for(int k=0 ; k<nz ; k++)			
	                        for(int i=0 ; i<nx ; i++)  
	                 {
	                         if (Solid[i][0][k]==procind)
	                                 bcly[bclyn]=Solid2[i][0][k],bclyn++;
	                         
	                         if (Solid[i][ny-1][k]==procind)
	                                 bcry[bcryn]=Solid2[i][ny-1][k],bcryn++;
	                 }   
	                 
	                 
	                 for(int j=0 ; j<ny ; j++)     
	                 for(int i=0 ; i<nx ; i++)  
	                 {
	                         if (Solid[i][j][0]==procind)
	                                 bclzn++;
	                         if (Solid[i][j][nz-1]==procind)
	                                 bcrzn++;
	                 }
	                 if (bclzn>0)
	                         bclz=new int[bclzn];
	                 if (bcrzn>0)
	                         bcrz= new int[bcrzn];
	                 bclzn=0;bcrzn=0;
	               for(int j=0 ; j<ny ; j++)     
	                 for(int i=0 ; i<nx ; i++) 
	                 {
	                         if (Solid[i][j][0]==procind)
	                                 bclz[bclzn]=Solid2[i][j][0],bclzn++;
	                         
	                         if (Solid[i][j][nz-1]==procind)
	                                 bcrz[bcrzn]=Solid2[i][j][nz-1],bcrzn++;
	                 }   
	                 
	      //---------------------------------------------------------------------------           
	                      
	                 
      		MPI_Waitall(2*com_n,request, status);
      		MPI_Testall(2*com_n,request,&mpi_test,status);
	      Count=sumss[procind];  


	//--------------hybrid lb dispersion------------------------
	u = new double*[Count+1];
	u[0] = new double[(Count+1)*3];
	for (int i=1;i<=Count;i++) 
	        u[i] = u[i-1]+3;

	//delete [] Solid_rank0;   
	   for (int i=0;i<com_n;i++)
	           delete [] buflocsend[i];
	   delete [] buflocsend;
	   
	   delete [] sumtmp;
	   
	   
	  // delete [] Solid2[0][0];
		//for (int i=0;i<nx;i++)
		//	delete [] Solid2[i];
		//delete [] Solid2;
    	
	


	if (updatesss==1)
	{
		if (rank==0)
		{
		cout<<"READING VELOCITY DATA"<<endl;
		cout<<endl;
		}


	int tmpsum[mpi_size];
	for (int i=0;i<mpi_size;i++)
		tmpsum[i]=3;	

	const int root_rank=0;
	
	double rho_0=1.0;
	
	
	MPI_Status status;
	MPI_Request request;

	double* rbuf_v;



	int nx_g[mpi_size];
	int disp[mpi_size];
	
	for (int i=0;i<mpi_size;i++)
		nx_g[i]=(sumss[i+1]+1)*3;


	

		disp[0]=0;
	
		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
		
	
	//if (rank==root_rank)
	//	rbuf_v = new double[disp[mpi_size-1]+nx_g[mpi_size-1]];
	rbuf_v = new double[pre_sum*3];
	
	int NX0=NX+1;
	int NY0=NY+1;
	int NZ0=NZ+1;
	//MPI_Gatherv(u[0],nx_g[rank],MPI_DOUBLE,rbuf_v,nx_g,disp,MPI_DOUBLE,root_rank,MPI_COMM_WORLD);

	ostringstream namevel;
	namevel<<"vel.bin";

	if (rank==0)
	{
	fstream fin;
	fin.open(namevel.str().c_str(),ios::in);
       if (fin.fail())
	        {
	        cout<<"\n file open error on " << namevel.str().c_str()<<endl;
	        exit(-1);
	        }
	        
	        
       fin.read((char *)(&rbuf_v[0]), sizeof(double)*(pre_sum*3));
        	
       fin.close();
	}
	//cout<<"@@@@@@@@@@@@"<<endl;

	MPI_Bcast(rbuf_v,pre_sum*3,MPI_DOUBLE,0,MPI_COMM_WORLD);
	



	tmpint=0;pore=0;
	for(int k=0 ; k<NZ+1 ; k++)			
	         for(int j=0 ; j<NY+1 ; j++)
	                for(int i=0 ; i<NX+1 ; i++)
			if (Solid3[i][j][k]==0)
			{
			
				if (Solid[i][j][k]==rank+1)
				{
				u[Solid2[i][j][k]][0]=rbuf_v[tmpint];
				//u[pore][1]=rbuf_v[tmpint+1];
				u[Solid2[i][j][k]][1]=rbuf_v[tmpint+1];
				u[Solid2[i][j][k]][2]=rbuf_v[tmpint+2];
				//pore++;
				}
			tmpint+=3;		
			}


	delete [] rbuf_v;

	cout<<pre_sum*3<<"        "<<tmpint<<"        previous sum velocity and reading velocity numbers"<<endl;
	}
	

	//***********************update part 2**********************************
	
	
	FILE *ftest;
	ftest = fopen("update.txt", "r");

	if(ftest == NULL)
		{
		cout << "\n The pore geometry file (" << filename <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(-1);
		}
		fclose(ftest);
	fstream fin;
	fin.open("update.txt");
	while (!fin.eof())
	{
	fin >> tmpint >> upx >> upy >> upz >> updoux >> updouy >> updouz;fin.getline(dummy, NCHAR);
	if ((Solid[upx][upy][upz]==rank+1) and (tmpint==0))
		{
			u[Solid2[upx][upy][upz]][0]=updoux;
			u[Solid2[upx][upy][upz]][1]=updouy;
			u[Solid2[upx][upy][upz]][2]=updouz;
		}
	}
	
	fin.close();
	//********************************************************************
	

	
	//----------------------------------------------------------

	delete [] Solid2[0][0];
		for (int i=0;i<nx;i++)
			delete [] Solid2[i];
		delete [] Solid2;

	delete [] Solid3[0][0];
		for (int i=0;i<nx;i++)
			delete [] Solid3[i];
		delete [] Solid3;
	  
	   if (rank>0)
	   {
	      
	           delete [] Solid[0][0];
			for (int i=0;i<nx;i++)
			delete [] Solid[i];
		delete [] Solid;
	   }
	   
}



void init(double* rho, double** u, double** f,int*** Solid)
{	
      

//	int rank = MPI :: COMM_WORLD . Get_rank ();
//	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	

	Zoom=1;
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
			if (updatesss!=1)
			{
			u[i][0]=inivx;
			u[i][1]=inivy;
			u[i][2]=inivz;
			}
			u_tmp[0]=u[i][0];
			u_tmp[1]=u[i][1];
			u_tmp[2]=u[i][2];

			rho[i]=1.0;
			
			if (nnf_module==1)
			niu_vn[i]=niu;
			

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


inline double feq(int k,double rho, double u[3])
{
	double eu,uv,feq;
        double c2,c4;

	c2=lat_c*lat_c;c4=c2*c2;
	eu=(elat[k][0]*u[0]+elat[k][1]*u[1]+elat[k][2]*u[2]);
	uv=(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);// WITH FORCE TERM:GRAVITY IN X DIRECTION
	feq=w[k]*rho*(1.0+3.0*eu/c2+4.5*eu*eu/c4-1.5*uv/c2);
	return feq;
}





void collision(double* rho,double** u,double** f,double** F,int* SupInv,int*** Solid, int* Sl, int* Sr)
{

	MPI_Status status[com_n*2] ;
	MPI_Request request[com_n*2];
	int mpi_test;
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	int procn=mpi_size; 
	int procind=rank+1;
	
	int testl1,testl2;
	int sumtmp[com_n];
	for (int i=0;i<com_n;i++)
	        sumtmp[i]=0;
	
	
double lm0,lm1,sum;
double usqr,vsqr;
double F_hat[19],GuoF[19],f_eq[19],u_tmp[3];
double m_l[19],m_inv_l[19];
int i,j,m,ip,jp,kp;


const double c_l=lat_c;	


	
	for(int ci=1;ci<=Count;ci++)	
	

		{	


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
			
			if (nei_loc[ci][mi]>0)
			        F[nei_loc[ci][mi]][mi]=m_inv_l[mi];
			else
	                                   if (nei_loc[ci][mi]==0)
	                                            F[ci][LR[mi]]=m_inv_l[mi];
	                                    else
	                                    {
	                                            bufsend[-nei_loc[ci][mi]-1][sumtmp[-nei_loc[ci][mi]-1]]=m_inv_l[mi];
	                                            sumtmp[-nei_loc[ci][mi]-1]++;
	                                    }
			        
			        
			        
			        
			        
			
			
			}
			//============================================================================
			
			
			}

  for (int i=0;i<com_n;i++)
	{
	MPI_Isend(bufsend[i],bufinfo[com_ind[i]], MPI_DOUBLE, com_ind[i]-1, (procind-1)*procn+com_ind[i]-1, MPI_COMM_WORLD,&request[2*i]);
	                       
	MPI_Irecv(bufrecv[i],bufinfo[com_ind[i]], MPI_DOUBLE, com_ind[i]-1, (com_ind[i]-1)*procn+procind-1, MPI_COMM_WORLD,&request[2*i+1]);		
	}

	               MPI_Waitall(2*com_n,request, status);
	               MPI_Testall(2*com_n,request,&mpi_test,status);

//cout<<"@@@@@@@@@@@@@@@1	"<<rank<<endl;	               
	               for (int i=0;i<com_n;i++)
	                       {
	                               for (int j=0;j<bufinfo[com_ind[i]];j++)
	                                       {
	                                               testl1=(int)(buflocrecv[i][j]/19);
	                                               testl2=(int)(buflocrecv[i][j]%19);
	                                               //cout<<testarr[testl1*19+testl2]<<"  "<<procind<<"   "<<testl1<<"        "<<testl2<<"        "<<buflocrecv[i][j]<<endl;
	                                               //cout<<buflocsend[i][j]<<"  ";
	                                               //testarr[testl1*19+testl2]=1;
	                                               F[testl1][testl2]=bufrecv[i][j];
						
	                                       }
	                       }
	




}



void collision_nnf(double* rho,double** u,double** f,double** F,int* SupInv,int*** Solid, int* Sl, int* Sr)
{

	MPI_Status status[4] ;
	MPI_Request request[4];
	int mpi_test;
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	int procn=mpi_size; 
	int procind=rank+1;


	int testl1,testl2;
	int sumtmp[com_n];
	for (int i=0;i<com_n;i++)
	        sumtmp[i]=0;


double lm0,lm1,sum;
double usqr,vsqr;
double F_hat[19],GuoF[19],f_eq[19],u_tmp[3];
double m_l[19],m_inv_l[19];
int i,j,m,ip,jp,kp;

//===============NNF module=======================
double pxx,pww,pxy,pyz,pxz;
double sxx,syy,szz,sxy,sxz,syz;
double s_other,s_vn;

double feq1;
//=============================================



const double c_l=lat_c;	

	
	
	for(int ci=1;ci<=Count;ci++)	
	

		{	
				
			
	


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
			
		

//==========================
m_l[0]=+1.000*f[ci][0]+1.000*f[ci][1]+1.000*f[ci][2]+1.000*f[ci][3]+1.000*f[ci][4]+1.000*f[ci][5]+1.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]+1.000*f[ci][9]+1.000*f[ci][10]+1.000*f[ci][11]+1.000*f[ci][12]+1.000*f[ci][13]+1.000*f[ci][14]+1.000*f[ci][15]+1.000*f[ci][16]+1.000*f[ci][17]+1.000*f[ci][18];

F_hat[0]=+1.000*GuoF[0]+1.000*GuoF[1]+1.000*GuoF[2]+1.000*GuoF[3]+1.000*GuoF[4]+1.000*GuoF[5]+1.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]+1.000*GuoF[9]+1.000*GuoF[10]+1.000*GuoF[11]+1.000*GuoF[12]+1.000*GuoF[13]+1.000*GuoF[14]+1.000*GuoF[15]+1.000*GuoF[16]+1.000*GuoF[17]+1.000*GuoF[18];

meq[0]=+1.000*f_eq[0]+1.000*f_eq[1]+1.000*f_eq[2]+1.000*f_eq[3]+1.000*f_eq[4]+1.000*f_eq[5]+1.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]+1.000*f_eq[9]+1.000*f_eq[10]+1.000*f_eq[11]+1.000*f_eq[12]+1.000*f_eq[13]+1.000*f_eq[14]+1.000*f_eq[15]+1.000*f_eq[16]+1.000*f_eq[17]+1.000*f_eq[18];

//F_hat[0]*=(1-0.5*S[0]);
//m_l[0]=m_l[0]-S[0]*(m_l[0]-meq[0])+dt*F_hat[0];
//=======================================

m_l[1]=-30.000*f[ci][0]-11.000*f[ci][1]-11.000*f[ci][2]-11.000*f[ci][3]-11.000*f[ci][4]-11.000*f[ci][5]-11.000*f[ci][6]+8.000*f[ci][7]+8.000*f[ci][8]+8.000*f[ci][9]+8.000*f[ci][10]+8.000*f[ci][11]+8.000*f[ci][12]+8.000*f[ci][13]+8.000*f[ci][14]+8.000*f[ci][15]+8.000*f[ci][16]+8.000*f[ci][17]+8.000*f[ci][18];

F_hat[1]=-30.000*GuoF[0]-11.000*GuoF[1]-11.000*GuoF[2]-11.000*GuoF[3]-11.000*GuoF[4]-11.000*GuoF[5]-11.000*GuoF[6]+8.000*GuoF[7]+8.000*GuoF[8]+8.000*GuoF[9]+8.000*GuoF[10]+8.000*GuoF[11]+8.000*GuoF[12]+8.000*GuoF[13]+8.000*GuoF[14]+8.000*GuoF[15]+8.000*GuoF[16]+8.000*GuoF[17]+8.000*GuoF[18];

meq[1]=-30.000*f_eq[0]-11.000*f_eq[1]-11.000*f_eq[2]-11.000*f_eq[3]-11.000*f_eq[4]-11.000*f_eq[5]-11.000*f_eq[6]+8.000*f_eq[7]+8.000*f_eq[8]+8.000*f_eq[9]+8.000*f_eq[10]+8.000*f_eq[11]+8.000*f_eq[12]+8.000*f_eq[13]+8.000*f_eq[14]+8.000*f_eq[15]+8.000*f_eq[16]+8.000*f_eq[17]+8.000*f_eq[18];

//F_hat[1]*=(1-0.5*S[1]);
//m_l[1]=m_l[1]-S[1]*(m_l[1]-meq[1])+dt*F_hat[1];
//=======================================

m_l[2]=+12.000*f[ci][0]-4.000*f[ci][1]-4.000*f[ci][2]-4.000*f[ci][3]-4.000*f[ci][4]-4.000*f[ci][5]-4.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]+1.000*f[ci][9]+1.000*f[ci][10]+1.000*f[ci][11]+1.000*f[ci][12]+1.000*f[ci][13]+1.000*f[ci][14]+1.000*f[ci][15]+1.000*f[ci][16]+1.000*f[ci][17]+1.000*f[ci][18];

F_hat[2]=+12.000*GuoF[0]-4.000*GuoF[1]-4.000*GuoF[2]-4.000*GuoF[3]-4.000*GuoF[4]-4.000*GuoF[5]-4.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]+1.000*GuoF[9]+1.000*GuoF[10]+1.000*GuoF[11]+1.000*GuoF[12]+1.000*GuoF[13]+1.000*GuoF[14]+1.000*GuoF[15]+1.000*GuoF[16]+1.000*GuoF[17]+1.000*GuoF[18];

meq[2]=+12.000*f_eq[0]-4.000*f_eq[1]-4.000*f_eq[2]-4.000*f_eq[3]-4.000*f_eq[4]-4.000*f_eq[5]-4.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]+1.000*f_eq[9]+1.000*f_eq[10]+1.000*f_eq[11]+1.000*f_eq[12]+1.000*f_eq[13]+1.000*f_eq[14]+1.000*f_eq[15]+1.000*f_eq[16]+1.000*f_eq[17]+1.000*f_eq[18];

//F_hat[2]*=(1-0.5*S[2]);
//m_l[2]=m_l[2]-S[2]*(m_l[2]-meq[2])+dt*F_hat[2];
//=======================================

m_l[3]=+0.000*f[ci][0]+1.000*f[ci][1]-1.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+1.000*f[ci][7]-1.000*f[ci][8]+1.000*f[ci][9]-1.000*f[ci][10]+0.000*f[ci][11]+0.000*f[ci][12]+0.000*f[ci][13]+0.000*f[ci][14]+1.000*f[ci][15]-1.000*f[ci][16]+1.000*f[ci][17]-1.000*f[ci][18];

F_hat[3]=+0.000*GuoF[0]+1.000*GuoF[1]-1.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+1.000*GuoF[7]-1.000*GuoF[8]+1.000*GuoF[9]-1.000*GuoF[10]+0.000*GuoF[11]+0.000*GuoF[12]+0.000*GuoF[13]+0.000*GuoF[14]+1.000*GuoF[15]-1.000*GuoF[16]+1.000*GuoF[17]-1.000*GuoF[18];

meq[3]=+0.000*f_eq[0]+1.000*f_eq[1]-1.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+1.000*f_eq[7]-1.000*f_eq[8]+1.000*f_eq[9]-1.000*f_eq[10]+0.000*f_eq[11]+0.000*f_eq[12]+0.000*f_eq[13]+0.000*f_eq[14]+1.000*f_eq[15]-1.000*f_eq[16]+1.000*f_eq[17]-1.000*f_eq[18];

//F_hat[3]*=(1-0.5*S[3]);
//m_l[3]=m_l[3]-S[3]*(m_l[3]-meq[3])+dt*F_hat[3];
//=======================================

m_l[4]=+0.000*f[ci][0]-4.000*f[ci][1]+4.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+1.000*f[ci][7]-1.000*f[ci][8]+1.000*f[ci][9]-1.000*f[ci][10]+0.000*f[ci][11]+0.000*f[ci][12]+0.000*f[ci][13]+0.000*f[ci][14]+1.000*f[ci][15]-1.000*f[ci][16]+1.000*f[ci][17]-1.000*f[ci][18];

F_hat[4]=+0.000*GuoF[0]-4.000*GuoF[1]+4.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+1.000*GuoF[7]-1.000*GuoF[8]+1.000*GuoF[9]-1.000*GuoF[10]+0.000*GuoF[11]+0.000*GuoF[12]+0.000*GuoF[13]+0.000*GuoF[14]+1.000*GuoF[15]-1.000*GuoF[16]+1.000*GuoF[17]-1.000*GuoF[18];

meq[4]=+0.000*f_eq[0]-4.000*f_eq[1]+4.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+1.000*f_eq[7]-1.000*f_eq[8]+1.000*f_eq[9]-1.000*f_eq[10]+0.000*f_eq[11]+0.000*f_eq[12]+0.000*f_eq[13]+0.000*f_eq[14]+1.000*f_eq[15]-1.000*f_eq[16]+1.000*f_eq[17]-1.000*f_eq[18];

//F_hat[4]*=(1-0.5*S[4]);
//m_l[4]=m_l[4]-S[4]*(m_l[4]-meq[4])+dt*F_hat[4];
//=======================================

m_l[5]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+1.000*f[ci][3]-1.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]-1.000*f[ci][9]-1.000*f[ci][10]+1.000*f[ci][11]-1.000*f[ci][12]+1.000*f[ci][13]-1.000*f[ci][14]+0.000*f[ci][15]+0.000*f[ci][16]+0.000*f[ci][17]+0.000*f[ci][18];

F_hat[5]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+1.000*GuoF[3]-1.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]-1.000*GuoF[9]-1.000*GuoF[10]+1.000*GuoF[11]-1.000*GuoF[12]+1.000*GuoF[13]-1.000*GuoF[14]+0.000*GuoF[15]+0.000*GuoF[16]+0.000*GuoF[17]+0.000*GuoF[18];

meq[5]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+1.000*f_eq[3]-1.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]-1.000*f_eq[9]-1.000*f_eq[10]+1.000*f_eq[11]-1.000*f_eq[12]+1.000*f_eq[13]-1.000*f_eq[14]+0.000*f_eq[15]+0.000*f_eq[16]+0.000*f_eq[17]+0.000*f_eq[18];

//F_hat[5]*=(1-0.5*S[5]);
//m_l[5]=m_l[5]-S[5]*(m_l[5]-meq[5])+dt*F_hat[5];
//=======================================

m_l[6]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]-4.000*f[ci][3]+4.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]-1.000*f[ci][9]-1.000*f[ci][10]+1.000*f[ci][11]-1.000*f[ci][12]+1.000*f[ci][13]-1.000*f[ci][14]+0.000*f[ci][15]+0.000*f[ci][16]+0.000*f[ci][17]+0.000*f[ci][18];

F_hat[6]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]-4.000*GuoF[3]+4.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]-1.000*GuoF[9]-1.000*GuoF[10]+1.000*GuoF[11]-1.000*GuoF[12]+1.000*GuoF[13]-1.000*GuoF[14]+0.000*GuoF[15]+0.000*GuoF[16]+0.000*GuoF[17]+0.000*GuoF[18];

meq[6]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]-4.000*f_eq[3]+4.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]-1.000*f_eq[9]-1.000*f_eq[10]+1.000*f_eq[11]-1.000*f_eq[12]+1.000*f_eq[13]-1.000*f_eq[14]+0.000*f_eq[15]+0.000*f_eq[16]+0.000*f_eq[17]+0.000*f_eq[18];

//F_hat[6]*=(1-0.5*S[6]);
//m_l[6]=m_l[6]-S[6]*(m_l[6]-meq[6])+dt*F_hat[6];
//=======================================

m_l[7]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+1.000*f[ci][5]-1.000*f[ci][6]+0.000*f[ci][7]+0.000*f[ci][8]+0.000*f[ci][9]+0.000*f[ci][10]+1.000*f[ci][11]+1.000*f[ci][12]-1.000*f[ci][13]-1.000*f[ci][14]+1.000*f[ci][15]+1.000*f[ci][16]-1.000*f[ci][17]-1.000*f[ci][18];

F_hat[7]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+1.000*GuoF[5]-1.000*GuoF[6]+0.000*GuoF[7]+0.000*GuoF[8]+0.000*GuoF[9]+0.000*GuoF[10]+1.000*GuoF[11]+1.000*GuoF[12]-1.000*GuoF[13]-1.000*GuoF[14]+1.000*GuoF[15]+1.000*GuoF[16]-1.000*GuoF[17]-1.000*GuoF[18];

meq[7]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+1.000*f_eq[5]-1.000*f_eq[6]+0.000*f_eq[7]+0.000*f_eq[8]+0.000*f_eq[9]+0.000*f_eq[10]+1.000*f_eq[11]+1.000*f_eq[12]-1.000*f_eq[13]-1.000*f_eq[14]+1.000*f_eq[15]+1.000*f_eq[16]-1.000*f_eq[17]-1.000*f_eq[18];

//F_hat[7]*=(1-0.5*S[7]);
//m_l[7]=m_l[7]-S[7]*(m_l[7]-meq[7])+dt*F_hat[7];
//=======================================

m_l[8]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]-4.000*f[ci][5]+4.000*f[ci][6]+0.000*f[ci][7]+0.000*f[ci][8]+0.000*f[ci][9]+0.000*f[ci][10]+1.000*f[ci][11]+1.000*f[ci][12]-1.000*f[ci][13]-1.000*f[ci][14]+1.000*f[ci][15]+1.000*f[ci][16]-1.000*f[ci][17]-1.000*f[ci][18];

F_hat[8]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]-4.000*GuoF[5]+4.000*GuoF[6]+0.000*GuoF[7]+0.000*GuoF[8]+0.000*GuoF[9]+0.000*GuoF[10]+1.000*GuoF[11]+1.000*GuoF[12]-1.000*GuoF[13]-1.000*GuoF[14]+1.000*GuoF[15]+1.000*GuoF[16]-1.000*GuoF[17]-1.000*GuoF[18];

meq[8]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]-4.000*f_eq[5]+4.000*f_eq[6]+0.000*f_eq[7]+0.000*f_eq[8]+0.000*f_eq[9]+0.000*f_eq[10]+1.000*f_eq[11]+1.000*f_eq[12]-1.000*f_eq[13]-1.000*f_eq[14]+1.000*f_eq[15]+1.000*f_eq[16]-1.000*f_eq[17]-1.000*f_eq[18];

//F_hat[8]*=(1-0.5*S[8]);
//m_l[8]=m_l[8]-S[8]*(m_l[8]-meq[8])+dt*F_hat[8];
//=======================================

m_l[9]=+0.000*f[ci][0]+2.000*f[ci][1]+2.000*f[ci][2]-1.000*f[ci][3]-1.000*f[ci][4]-1.000*f[ci][5]-1.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]+1.000*f[ci][9]+1.000*f[ci][10]-2.000*f[ci][11]-2.000*f[ci][12]-2.000*f[ci][13]-2.000*f[ci][14]+1.000*f[ci][15]+1.000*f[ci][16]+1.000*f[ci][17]+1.000*f[ci][18];

F_hat[9]=+0.000*GuoF[0]+2.000*GuoF[1]+2.000*GuoF[2]-1.000*GuoF[3]-1.000*GuoF[4]-1.000*GuoF[5]-1.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]+1.000*GuoF[9]+1.000*GuoF[10]-2.000*GuoF[11]-2.000*GuoF[12]-2.000*GuoF[13]-2.000*GuoF[14]+1.000*GuoF[15]+1.000*GuoF[16]+1.000*GuoF[17]+1.000*GuoF[18];

meq[9]=+0.000*f_eq[0]+2.000*f_eq[1]+2.000*f_eq[2]-1.000*f_eq[3]-1.000*f_eq[4]-1.000*f_eq[5]-1.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]+1.000*f_eq[9]+1.000*f_eq[10]-2.000*f_eq[11]-2.000*f_eq[12]-2.000*f_eq[13]-2.000*f_eq[14]+1.000*f_eq[15]+1.000*f_eq[16]+1.000*f_eq[17]+1.000*f_eq[18];

//F_hat[9]*=(1-0.5*S[9]);
//m_l[9]=m_l[9]-S[9]*(m_l[9]-meq[9])+dt*F_hat[9];
//=======================================

m_l[10]=+0.000*f[ci][0]-4.000*f[ci][1]-4.000*f[ci][2]+2.000*f[ci][3]+2.000*f[ci][4]+2.000*f[ci][5]+2.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]+1.000*f[ci][9]+1.000*f[ci][10]-2.000*f[ci][11]-2.000*f[ci][12]-2.000*f[ci][13]-2.000*f[ci][14]+1.000*f[ci][15]+1.000*f[ci][16]+1.000*f[ci][17]+1.000*f[ci][18];

F_hat[10]=+0.000*GuoF[0]-4.000*GuoF[1]-4.000*GuoF[2]+2.000*GuoF[3]+2.000*GuoF[4]+2.000*GuoF[5]+2.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]+1.000*GuoF[9]+1.000*GuoF[10]-2.000*GuoF[11]-2.000*GuoF[12]-2.000*GuoF[13]-2.000*GuoF[14]+1.000*GuoF[15]+1.000*GuoF[16]+1.000*GuoF[17]+1.000*GuoF[18];

meq[10]=+0.000*f_eq[0]-4.000*f_eq[1]-4.000*f_eq[2]+2.000*f_eq[3]+2.000*f_eq[4]+2.000*f_eq[5]+2.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]+1.000*f_eq[9]+1.000*f_eq[10]-2.000*f_eq[11]-2.000*f_eq[12]-2.000*f_eq[13]-2.000*f_eq[14]+1.000*f_eq[15]+1.000*f_eq[16]+1.000*f_eq[17]+1.000*f_eq[18];

//F_hat[10]*=(1-0.5*S[10]);
//m_l[10]=m_l[10]-S[10]*(m_l[10]-meq[10])+dt*F_hat[10];
//=======================================

m_l[11]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+1.000*f[ci][3]+1.000*f[ci][4]-1.000*f[ci][5]-1.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]+1.000*f[ci][9]+1.000*f[ci][10]+0.000*f[ci][11]+0.000*f[ci][12]+0.000*f[ci][13]+0.000*f[ci][14]-1.000*f[ci][15]-1.000*f[ci][16]-1.000*f[ci][17]-1.000*f[ci][18];

F_hat[11]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+1.000*GuoF[3]+1.000*GuoF[4]-1.000*GuoF[5]-1.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]+1.000*GuoF[9]+1.000*GuoF[10]+0.000*GuoF[11]+0.000*GuoF[12]+0.000*GuoF[13]+0.000*GuoF[14]-1.000*GuoF[15]-1.000*GuoF[16]-1.000*GuoF[17]-1.000*GuoF[18];

meq[11]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+1.000*f_eq[3]+1.000*f_eq[4]-1.000*f_eq[5]-1.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]+1.000*f_eq[9]+1.000*f_eq[10]+0.000*f_eq[11]+0.000*f_eq[12]+0.000*f_eq[13]+0.000*f_eq[14]-1.000*f_eq[15]-1.000*f_eq[16]-1.000*f_eq[17]-1.000*f_eq[18];

//F_hat[11]*=(1-0.5*S[11]);
//m_l[11]=m_l[11]-S[11]*(m_l[11]-meq[11])+dt*F_hat[11];
//=======================================

m_l[12]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]-2.000*f[ci][3]-2.000*f[ci][4]+2.000*f[ci][5]+2.000*f[ci][6]+1.000*f[ci][7]+1.000*f[ci][8]+1.000*f[ci][9]+1.000*f[ci][10]+0.000*f[ci][11]+0.000*f[ci][12]+0.000*f[ci][13]+0.000*f[ci][14]-1.000*f[ci][15]-1.000*f[ci][16]-1.000*f[ci][17]-1.000*f[ci][18];

F_hat[12]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]-2.000*GuoF[3]-2.000*GuoF[4]+2.000*GuoF[5]+2.000*GuoF[6]+1.000*GuoF[7]+1.000*GuoF[8]+1.000*GuoF[9]+1.000*GuoF[10]+0.000*GuoF[11]+0.000*GuoF[12]+0.000*GuoF[13]+0.000*GuoF[14]-1.000*GuoF[15]-1.000*GuoF[16]-1.000*GuoF[17]-1.000*GuoF[18];

meq[12]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]-2.000*f_eq[3]-2.000*f_eq[4]+2.000*f_eq[5]+2.000*f_eq[6]+1.000*f_eq[7]+1.000*f_eq[8]+1.000*f_eq[9]+1.000*f_eq[10]+0.000*f_eq[11]+0.000*f_eq[12]+0.000*f_eq[13]+0.000*f_eq[14]-1.000*f_eq[15]-1.000*f_eq[16]-1.000*f_eq[17]-1.000*f_eq[18];

//F_hat[12]*=(1-0.5*S[12]);
//m_l[12]=m_l[12]-S[12]*(m_l[12]-meq[12])+dt*F_hat[12];
//=======================================

m_l[13]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+1.000*f[ci][7]-1.000*f[ci][8]-1.000*f[ci][9]+1.000*f[ci][10]+0.000*f[ci][11]+0.000*f[ci][12]+0.000*f[ci][13]+0.000*f[ci][14]+0.000*f[ci][15]+0.000*f[ci][16]+0.000*f[ci][17]+0.000*f[ci][18];

F_hat[13]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+1.000*GuoF[7]-1.000*GuoF[8]-1.000*GuoF[9]+1.000*GuoF[10]+0.000*GuoF[11]+0.000*GuoF[12]+0.000*GuoF[13]+0.000*GuoF[14]+0.000*GuoF[15]+0.000*GuoF[16]+0.000*GuoF[17]+0.000*GuoF[18];

meq[13]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+1.000*f_eq[7]-1.000*f_eq[8]-1.000*f_eq[9]+1.000*f_eq[10]+0.000*f_eq[11]+0.000*f_eq[12]+0.000*f_eq[13]+0.000*f_eq[14]+0.000*f_eq[15]+0.000*f_eq[16]+0.000*f_eq[17]+0.000*f_eq[18];

//F_hat[13]*=(1-0.5*S[13]);
//m_l[13]=m_l[13]-S[13]*(m_l[13]-meq[13])+dt*F_hat[13];
//=======================================

m_l[14]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+0.000*f[ci][7]+0.000*f[ci][8]+0.000*f[ci][9]+0.000*f[ci][10]+1.000*f[ci][11]-1.000*f[ci][12]-1.000*f[ci][13]+1.000*f[ci][14]+0.000*f[ci][15]+0.000*f[ci][16]+0.000*f[ci][17]+0.000*f[ci][18];

F_hat[14]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+0.000*GuoF[7]+0.000*GuoF[8]+0.000*GuoF[9]+0.000*GuoF[10]+1.000*GuoF[11]-1.000*GuoF[12]-1.000*GuoF[13]+1.000*GuoF[14]+0.000*GuoF[15]+0.000*GuoF[16]+0.000*GuoF[17]+0.000*GuoF[18];

meq[14]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+0.000*f_eq[7]+0.000*f_eq[8]+0.000*f_eq[9]+0.000*f_eq[10]+1.000*f_eq[11]-1.000*f_eq[12]-1.000*f_eq[13]+1.000*f_eq[14]+0.000*f_eq[15]+0.000*f_eq[16]+0.000*f_eq[17]+0.000*f_eq[18];

//F_hat[14]*=(1-0.5*S[14]);
//m_l[14]=m_l[14]-S[14]*(m_l[14]-meq[14])+dt*F_hat[14];
//=======================================

m_l[15]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+0.000*f[ci][7]+0.000*f[ci][8]+0.000*f[ci][9]+0.000*f[ci][10]+0.000*f[ci][11]+0.000*f[ci][12]+0.000*f[ci][13]+0.000*f[ci][14]+1.000*f[ci][15]-1.000*f[ci][16]-1.000*f[ci][17]+1.000*f[ci][18];

F_hat[15]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+0.000*GuoF[7]+0.000*GuoF[8]+0.000*GuoF[9]+0.000*GuoF[10]+0.000*GuoF[11]+0.000*GuoF[12]+0.000*GuoF[13]+0.000*GuoF[14]+1.000*GuoF[15]-1.000*GuoF[16]-1.000*GuoF[17]+1.000*GuoF[18];

meq[15]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+0.000*f_eq[7]+0.000*f_eq[8]+0.000*f_eq[9]+0.000*f_eq[10]+0.000*f_eq[11]+0.000*f_eq[12]+0.000*f_eq[13]+0.000*f_eq[14]+1.000*f_eq[15]-1.000*f_eq[16]-1.000*f_eq[17]+1.000*f_eq[18];

//F_hat[15]*=(1-0.5*S[15]);
//m_l[15]=m_l[15]-S[15]*(m_l[15]-meq[15])+dt*F_hat[15];
//=======================================

m_l[16]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+1.000*f[ci][7]-1.000*f[ci][8]+1.000*f[ci][9]-1.000*f[ci][10]+0.000*f[ci][11]+0.000*f[ci][12]+0.000*f[ci][13]+0.000*f[ci][14]-1.000*f[ci][15]+1.000*f[ci][16]-1.000*f[ci][17]+1.000*f[ci][18];

F_hat[16]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+1.000*GuoF[7]-1.000*GuoF[8]+1.000*GuoF[9]-1.000*GuoF[10]+0.000*GuoF[11]+0.000*GuoF[12]+0.000*GuoF[13]+0.000*GuoF[14]-1.000*GuoF[15]+1.000*GuoF[16]-1.000*GuoF[17]+1.000*GuoF[18];

meq[16]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+1.000*f_eq[7]-1.000*f_eq[8]+1.000*f_eq[9]-1.000*f_eq[10]+0.000*f_eq[11]+0.000*f_eq[12]+0.000*f_eq[13]+0.000*f_eq[14]-1.000*f_eq[15]+1.000*f_eq[16]-1.000*f_eq[17]+1.000*f_eq[18];

//F_hat[16]*=(1-0.5*S[16]);
//m_l[16]=m_l[16]-S[16]*(m_l[16]-meq[16])+dt*F_hat[16];
//=======================================

m_l[17]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]-1.000*f[ci][7]-1.000*f[ci][8]+1.000*f[ci][9]+1.000*f[ci][10]+1.000*f[ci][11]-1.000*f[ci][12]+1.000*f[ci][13]-1.000*f[ci][14]+0.000*f[ci][15]+0.000*f[ci][16]+0.000*f[ci][17]+0.000*f[ci][18];

F_hat[17]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]-1.000*GuoF[7]-1.000*GuoF[8]+1.000*GuoF[9]+1.000*GuoF[10]+1.000*GuoF[11]-1.000*GuoF[12]+1.000*GuoF[13]-1.000*GuoF[14]+0.000*GuoF[15]+0.000*GuoF[16]+0.000*GuoF[17]+0.000*GuoF[18];

meq[17]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]-1.000*f_eq[7]-1.000*f_eq[8]+1.000*f_eq[9]+1.000*f_eq[10]+1.000*f_eq[11]-1.000*f_eq[12]+1.000*f_eq[13]-1.000*f_eq[14]+0.000*f_eq[15]+0.000*f_eq[16]+0.000*f_eq[17]+0.000*f_eq[18];

//F_hat[17]*=(1-0.5*S[17]);
//m_l[17]=m_l[17]-S[17]*(m_l[17]-meq[17])+dt*F_hat[17];
//=======================================

m_l[18]=+0.000*f[ci][0]+0.000*f[ci][1]+0.000*f[ci][2]+0.000*f[ci][3]+0.000*f[ci][4]+0.000*f[ci][5]+0.000*f[ci][6]+0.000*f[ci][7]+0.000*f[ci][8]+0.000*f[ci][9]+0.000*f[ci][10]-1.000*f[ci][11]-1.000*f[ci][12]+1.000*f[ci][13]+1.000*f[ci][14]+1.000*f[ci][15]+1.000*f[ci][16]-1.000*f[ci][17]-1.000*f[ci][18];

F_hat[18]=+0.000*GuoF[0]+0.000*GuoF[1]+0.000*GuoF[2]+0.000*GuoF[3]+0.000*GuoF[4]+0.000*GuoF[5]+0.000*GuoF[6]+0.000*GuoF[7]+0.000*GuoF[8]+0.000*GuoF[9]+0.000*GuoF[10]-1.000*GuoF[11]-1.000*GuoF[12]+1.000*GuoF[13]+1.000*GuoF[14]+1.000*GuoF[15]+1.000*GuoF[16]-1.000*GuoF[17]-1.000*GuoF[18];

meq[18]=+0.000*f_eq[0]+0.000*f_eq[1]+0.000*f_eq[2]+0.000*f_eq[3]+0.000*f_eq[4]+0.000*f_eq[5]+0.000*f_eq[6]+0.000*f_eq[7]+0.000*f_eq[8]+0.000*f_eq[9]+0.000*f_eq[10]-1.000*f_eq[11]-1.000*f_eq[12]+1.000*f_eq[13]+1.000*f_eq[14]+1.000*f_eq[15]+1.000*f_eq[16]-1.000*f_eq[17]-1.000*f_eq[18];

//F_hat[18]*=(1-0.5*S[18]);
//m_l[18]=m_l[18]-S[18]*(m_l[18]-meq[18])+dt*F_hat[18];



//===============NNF module=======================
if (n>10)
{
s_vn=1.0/(niu_vn[ci]/(c_s2*dt)+0.5);
//if (n==0)
//cout<<niu_vn[ci]<<"           &&& "	<<s_vn<<endl;

//==============MRT===========================
/*
pxx=m_l[9]/3;
pww=m_l[11];
pxy=m_l[13];
pyz=m_l[14];
pxz=m_l[15];

sxx=-(1-s_vn/2.0)*(1.0/3.0*m_l[1]+pxx-u[ci][0]*u[ci][0]);
syy=-(1-s_vn/2.0)*(1.0/3.0*m_l[1]-0.5*pxx+0.5*pww-u[ci][1]*u[ci][1]);
szz=-(1-s_vn/2.0)*(1.0/3.0*m_l[1]-0.5*pxx-0.5*pww-u[ci][2]*u[ci][2]);
sxy=-(1-s_vn/2.0)*(pxy-u[ci][0]*u[ci][1]);
syz=-(1-s_vn/2.0)*(pyz-u[ci][1]*u[ci][2]);
sxz=-(1-s_vn/2.0)*(pxz-u[ci][0]*u[ci][2]);

sxx=(sxx+rho[ci]/3.0)/(niu_vn[ci]*rho[ci])*0.5;
syy=(syy+rho[ci]/3.0)/(niu_vn[ci]*rho[ci])*0.5;
szz=(szz+rho[ci]/3.0)/(niu_vn[ci]*rho[ci])*0.5;
sxy=sxy/(niu_vn[ci]*rho[ci])*0.5;
syz=syz/(niu_vn[ci]*rho[ci])*0.5;
sxz=sxz/(niu_vn[ci]*rho[ci])*0.5;
*/
//============================================



//============SRT================
sxx=0;syy=0;szz=0;sxy=0;syz=0;sxz=0;
for (int sint=0;sint<19;sint++)
        {       
                feq1=f[ci][sint]-f_eq[sint];
                sxx+=e[sint][0]*e[sint][0]*feq1;
                syy+=e[sint][1]*e[sint][1]*feq1;
                szz+=e[sint][2]*e[sint][2]*feq1;
                sxy+=e[sint][0]*e[sint][1]*feq1;
                syz+=e[sint][1]*e[sint][2]*feq1;
                sxz+=e[sint][0]*e[sint][2]*feq1;    
        }
       
        
        sxx=-s_vn/(2*rho[ci]*c_s2*dt)*sxx;
        syy=-s_vn/(2*rho[ci]*c_s2*dt)*syy;
        szz=-s_vn/(2*rho[ci]*c_s2*dt)*szz;
        sxy=-s_vn/(2*rho[ci]*c_s2*dt)*sxy;
        syz=-s_vn/(2*rho[ci]*c_s2*dt)*syz;
        sxz=-s_vn/(2*rho[ci]*c_s2*dt)*sxz;
//==================================




        niu_vn[ci]=nnf_m*exp((nnf_n-1)*log(sqrt(sxx*sxx+syy*syy+szz*szz+2*sxy*sxy+2*syz*syz+2*sxz*sxz)));
        
        tau_f=niu_vn[ci]/(c_s2*dt)+0.5;
	s_vn=1/tau_f;
	s_other=8*(2-s_vn)/(8-s_vn);







	S[1]=s_vn;
	S[2]=s_vn;
	S[4]=s_other;
	S[6]=s_other;
	S[8]=s_other;
	S[9]=s_vn;
	S[10]=s_vn;
	S[11]=s_vn;
	S[12]=s_vn;
	S[13]=s_vn;
	S[14]=s_vn;
	S[15]=s_vn;
	S[16]=s_other;
	S[17]=s_other;
	S[18]=s_other;

	//=============================================
}



//========================================RELAXATION====================
F_hat[0]*=(1-0.5*S[0]);
F_hat[1]*=(1-0.5*S[1]);
F_hat[2]*=(1-0.5*S[2]);
F_hat[3]*=(1-0.5*S[3]);
F_hat[4]*=(1-0.5*S[4]);
F_hat[5]*=(1-0.5*S[5]);
F_hat[6]*=(1-0.5*S[6]);
F_hat[7]*=(1-0.5*S[7]);
F_hat[8]*=(1-0.5*S[8]);
F_hat[9]*=(1-0.5*S[9]);
F_hat[10]*=(1-0.5*S[10]);
F_hat[11]*=(1-0.5*S[11]);
F_hat[12]*=(1-0.5*S[12]);
F_hat[13]*=(1-0.5*S[13]);
F_hat[14]*=(1-0.5*S[14]);
F_hat[15]*=(1-0.5*S[15]);
F_hat[16]*=(1-0.5*S[16]);
F_hat[17]*=(1-0.5*S[17]);
F_hat[18]*=(1-0.5*S[18]);

m_l[0]=m_l[0]-S[0]*(m_l[0]-meq[0])+dt*F_hat[0];
m_l[1]=m_l[1]-S[1]*(m_l[1]-meq[1])+dt*F_hat[1];
m_l[2]=m_l[2]-S[2]*(m_l[2]-meq[2])+dt*F_hat[2];
m_l[3]=m_l[3]-S[3]*(m_l[3]-meq[3])+dt*F_hat[3];
m_l[4]=m_l[4]-S[4]*(m_l[4]-meq[4])+dt*F_hat[4];
m_l[5]=m_l[5]-S[5]*(m_l[5]-meq[5])+dt*F_hat[5];
m_l[6]=m_l[6]-S[6]*(m_l[6]-meq[6])+dt*F_hat[6];
m_l[7]=m_l[7]-S[7]*(m_l[7]-meq[7])+dt*F_hat[7];
m_l[8]=m_l[8]-S[8]*(m_l[8]-meq[8])+dt*F_hat[8];
m_l[9]=m_l[9]-S[9]*(m_l[9]-meq[9])+dt*F_hat[9];
m_l[10]=m_l[10]-S[10]*(m_l[10]-meq[10])+dt*F_hat[10];
m_l[11]=m_l[11]-S[11]*(m_l[11]-meq[11])+dt*F_hat[11];
m_l[12]=m_l[12]-S[12]*(m_l[12]-meq[12])+dt*F_hat[12];
m_l[13]=m_l[13]-S[13]*(m_l[13]-meq[13])+dt*F_hat[13];
m_l[14]=m_l[14]-S[14]*(m_l[14]-meq[14])+dt*F_hat[14];
m_l[15]=m_l[15]-S[15]*(m_l[15]-meq[15])+dt*F_hat[15];
m_l[16]=m_l[16]-S[16]*(m_l[16]-meq[16])+dt*F_hat[16];
m_l[17]=m_l[17]-S[17]*(m_l[17]-meq[17])+dt*F_hat[17];
m_l[18]=m_l[18]-S[18]*(m_l[18]-meq[18])+dt*F_hat[18];
//====================================================================






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
			
			if (nei_loc[ci][mi]>0)
			        F[nei_loc[ci][mi]][mi]=m_inv_l[mi];
			else
	                                   if (nei_loc[ci][mi]==0)
	                                            F[ci][LR[mi]]=m_inv_l[mi];
	                                    else
	                                    {
	                                            bufsend[-nei_loc[ci][mi]-1][sumtmp[-nei_loc[ci][mi]-1]]=m_inv_l[mi];
	                                            sumtmp[-nei_loc[ci][mi]-1]++;
	                                    }
			        
			        
			        
			        
			        
			
			
			}
			//============================================================================
			
			
			}


	for (int i=0;i<com_n;i++)
	{
	MPI_Isend(bufsend[i],bufinfo[com_ind[i]], MPI_DOUBLE, com_ind[i]-1, (procind-1)*procn+com_ind[i]-1, MPI_COMM_WORLD,&request[2*i]);
	                       
	MPI_Irecv(bufrecv[i],bufinfo[com_ind[i]], MPI_DOUBLE, com_ind[i]-1, (com_ind[i]-1)*procn+procind-1, MPI_COMM_WORLD,&request[2*i+1]);		
	}


	
	MPI_Waitall(4,request, status);

	MPI_Testall(4,request,&mpi_test,status);	



		
			for (int i=0;i<com_n;i++)
	                       {
	                               for (int j=0;j<bufinfo[com_ind[i]];j++)
	                                       {
	                                               testl1=(int)(buflocrecv[i][j]/19);
	                                               testl2=(int)(buflocrecv[i][j]%19);
	                                             
	                                               F[testl1][testl2]=bufrecv[i][j];
						
	                                       }
	                       }

                        


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
					rho[i]+=f[i][k];//cout<<f[i][k]<<"  ";
					u[i][0]+=elat[k][0]*f[i][k];
					u[i][1]+=elat[k][1]*f[i][k];
					u[i][2]+=elat[k][2]*f[i][k];
					}
				
				//cout<<endl;
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
if ((Sub_BC==0) or (Sub_BC==1))
{

if ((yp-1)*(yn-1)==0)
	{
	if (yn==1)
	for (int i=0;i<bclyn;i++)
		for (int ks=0;ks<Q;ks++)
		F[bcly[i]][ks]=feq(ks,1.0,u_yn); 
	if (yp==1)
	for (int i=0;i<bcryn;i++)
		for (int ks=0;ks<Q;ks++)
		F[bcry[i]][ks]=feq(ks,1.0,u_yp); 
	}
 
/*
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
*/

if ((zp-1)*(zn-1)==0)	
	{
	if (zn==1)
	for (int i=0;i<bclzn;i++)
		for (int ks=0;ks<Q;ks++)
		F[bclz[i]][ks]=feq(ks,1.0,u_zn); 

	if (zp==1)
	for (int i=0;i<bcrzn;i++)
		for (int ks=0;ks<Q;ks++)
		F[bcrz[i]][ks]=feq(ks,1.0,u_zp); 
	}
 


/*	
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
		for (int ks=0;ks<Q;ks++)
		{
		if ((zp==1) && (Solid[i][j][NZ]>0))    
		        F[Solid[i][j][NZ]][ks]=feq(ks,1.0,u_zp); 
		if ((zn==1) && (Solid[i][j][0]>0))
		        F[Solid[i][j][0]][ks]=feq(ks,1.0,u_zn);
		}

*/

if ((xp-1)*(xn-1)==0)
	{
	if (xn==1)
	for (int i=0;i<bclxn;i++)
		for (int ks=0;ks<Q;ks++)
		F[bclx[i]][ks]=feq(ks,1.0,u_xn); 

	if (xp==1)
	for (int i=0;i<bcrxn;i++)
		for (int ks=0;ks<Q;ks++)
		F[bcrx[i]][ks]=feq(ks,1.0,u_xp); 
	}


/*
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
			
	
*/

}




//Non-Equilibrium boundary condition the distribution of BC points is updated using static pressure: rho=1.0
if (Sub_BC==2)
{
if (yp==1)
for (int i=0;i<bcryn;i++)
	for (int ks=0;ks<Q;ks++)
	if (e[ks][1]<0)
	F[bcry[i]][ks]=feq(LR[ks],1.0,u_yp)-F[bcry[i]][LR[ks]]+feq(ks,1.0,u_yp);

if (yn==1)
for (int i=0;i<bclyn;i++)
	for (int ks=0;ks<Q;ks++)
	if (e[ks][1]>0)
	F[bcly[i]][ks]=feq(LR[ks],1.0,u_yn)-F[bcly[i]][LR[ks]]+feq(ks,1.0,u_yn);


if (xp==1)
for (int i=0;i<bcrzn;i++)
	for (int ks=0;ks<Q;ks++)
	if (e[ks][2]<0)
	F[bcrz[i]][ks]=feq(LR[ks],1.0,u_zp)-F[bcrz[i]][LR[ks]]+feq(ks,1.0,u_zp);

if (xn==1)
for (int i=0;i<bclzn;i++)
	for (int ks=0;ks<Q;ks++)
	if (e[ks][2]>0)
	F[bclz[i]][ks]=feq(LR[ks],1.0,u_zn)-F[bclz[i]][LR[ks]]+feq(ks,1.0,u_zn);

if (zp==1)
for (int i=0;i<bcrxn;i++)
	for (int ks=0;ks<Q;ks++)
	if (e[ks][0]<0)
	F[bcrx[i]][ks]=feq(LR[ks],1.0,u_xp)-F[bcrx[i]][LR[ks]]+feq(ks,1.0,u_xp);

if (zn==1)
for (int i=0;i<bclxn;i++)
	for (int ks=0;ks<Q;ks++)
	if (e[ks][0]>0)
	F[bclx[i]][ks]=feq(LR[ks],1.0,u_xn)-F[bclx[i]][LR[ks]]+feq(ks,1.0,u_xn);


	
}

//Non-Equilibrium boundary condition use neighbouring points' value to update the BC points
if (Sub_BC==3)
{
/*
const int e[19][3]=
{{0,0,0},{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};
*/

if (yp==1)
for (int i=0;i<bcryn;i++)
	for (int ks=0;ks<Q;ks++)
	if (nei_loc[bcry[i]][4]>0)
	F[bcry[i]][ks]=feq(ks,rho[nei_loc[bcry[i]][4]],u_yp);
if (yn==1)
for (int i=0;i<bclyn;i++)
	for (int ks=0;ks<Q;ks++)
	if (nei_loc[bcly[i]][3]>0)
	F[bcly[i]][ks]=feq(ks,rho[nei_loc[bcly[i]][3]],u_yn);
if (xp==1)
for (int i=0;i<bcrxn;i++)
	for (int ks=0;ks<Q;ks++)
	if (nei_loc[bcrx[i]][2]>0)
	F[bcrx[i]][ks]=feq(ks,rho[nei_loc[bcrx[i]][2]],u_xp);
if (xn==1)
for (int i=0;i<bclxn;i++)
	for (int ks=0;ks<Q;ks++)
	if (nei_loc[bclx[i]][1]>0)
	F[bclx[i]][ks]=feq(ks,rho[nei_loc[bclx[i]][1]],u_xn);
if (zp==1)
for (int i=0;i<bcrzn;i++)
	for (int ks=0;ks<Q;ks++)
	if (nei_loc[bcrz[i]][6]>0)
	F[bcrz[i]][ks]=feq(ks,rho[nei_loc[bcrz[i]][6]],u_zp);
if (zn==1)
for (int i=0;i<bclzn;i++)
	for (int ks=0;ks<Q;ks++)
	if (nei_loc[bclz[i]][5]>0)
	F[bclz[i]][ks]=feq(ks,rho[nei_loc[bclz[i]][5]],u_zn);


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
if (yp==1)
for (int i=0;i<bcryn;i++)
	for (int ks=0;ks<Q;ks++)
		F[bcry[i]][ks]=feq(ks,rho_yp,u_ls);
if (yn==1)	
for (int i=0;i<bclyn;i++)
	for (int ks=0;ks<Q;ks++)
		F[bcly[i]][ks]=feq(ks,rho_yn,u_ls);
if (xp==1)
for (int i=0;i<bcrxn;i++)
	for (int ks=0;ks<Q;ks++)
		F[bcrx[i]][ks]=feq(ks,rho_xp,u_ls);
if (xn==1)	
for (int i=0;i<bclxn;i++)
	for (int ks=0;ks<Q;ks++)
		F[bclx[i]][ks]=feq(ks,rho_xn,u_ls);
if (zp==1)
for (int i=0;i<bcrzn;i++)
	for (int ks=0;ks<Q;ks++)
		F[bcrz[i]][ks]=feq(ks,rho_zp,u_ls);
if (zn==1)	
for (int i=0;i<bclzn;i++)
	for (int ks=0;ks<Q;ks++)
		F[bclz[i]][ks]=feq(ks,rho_zn,u_ls);
		
}
       


//Equilibruim boundary condition: Use the value of neighbouring points values to update the distributions of BC points
 if (Sub_BC==1)
{


if (yp==1)
for (int i=0;i<bcryn;i++)
	if (nei_loc[bcry[i]][4]>0)
	for (int ks=0;ks<Q;ks++)
	{
	u_ls[0]=u[nei_loc[bcry[i]][4]][0];
	u_ls[1]=u[nei_loc[bcry[i]][4]][1];
	u_ls[2]=u[nei_loc[bcry[i]][4]][2];
	F[bcry[i]][ks]=feq(ks,rho_yp,u_ls);
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bcry[i]][ks]=feq(ks,rho_yp,u_ls);
		}			

if (yn==1)	
for (int i=0;i<bclyn;i++)
	if (nei_loc[bcly[i]][3]>0)
	for (int ks=0;ks<Q;ks++)
	{
	u_ls[0]=u[nei_loc[bcly[i]][3]][0];
	u_ls[1]=u[nei_loc[bcly[i]][3]][1];
	u_ls[2]=u[nei_loc[bcly[i]][3]][2];
	F[bcly[i]][ks]=feq(ks,rho_yn,u_ls);
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bcry[i]][ks]=feq(ks,rho_yn,u_ls);
		}

if (xp==1)
for (int i=0;i<bcrxn;i++)
	if (nei_loc[bcrx[i]][2]>0)
	for (int ks=0;ks<Q;ks++)
	{
	u_ls[0]=u[nei_loc[bcrx[i]][2]][0];
	u_ls[1]=u[nei_loc[bcrx[i]][2]][1];
	u_ls[2]=u[nei_loc[bcrx[i]][2]][2];
	F[bcrx[i]][ks]=feq(ks,rho_xp,u_ls);
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bcry[i]][ks]=feq(ks,rho_xp,u_ls);
		}

if (xn==1)
for (int i=0;i<bclxn;i++)
	if (nei_loc[bclx[i]][1]>0)
	for (int ks=0;ks<Q;ks++)
	{
	u_ls[0]=u[nei_loc[bclx[i]][1]][0];
	u_ls[1]=u[nei_loc[bclx[i]][1]][1];
	u_ls[2]=u[nei_loc[bclx[i]][1]][2];
	F[bclx[i]][ks]=feq(ks,rho_xn,u_ls);
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bcry[i]][ks]=feq(ks,rho_xn,u_ls);
		}


if (zp==1)
for (int i=0;i<bcrzn;i++)
	if (nei_loc[bcrz[i]][6]>0)
	for (int ks=0;ks<Q;ks++)
	{
	u_ls[0]=u[nei_loc[bcrz[i]][6]][0];
	u_ls[1]=u[nei_loc[bcrz[i]][6]][1];
	u_ls[2]=u[nei_loc[bcrz[i]][6]][2];
	F[bcrz[i]][ks]=feq(ks,rho_zp,u_ls);
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bcry[i]][ks]=feq(ks,rho_zp,u_ls);
		}


if (zn==1)
for (int i=0;i<bclzn;i++)
	if (nei_loc[bclz[i]][5]>0)
	for (int ks=0;ks<Q;ks++)
	{
	u_ls[0]=u[nei_loc[bclz[i]][5]][0];
	u_ls[1]=u[nei_loc[bclz[i]][5]][1];
	u_ls[2]=u[nei_loc[bclz[i]][5]][2];
	F[bclz[i]][ks]=feq(ks,rho_zn,u_ls);
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bcry[i]][ks]=feq(ks,rho_zn,u_ls);
		}



}	  


//Non-equilibrium boundary condition  
if (Sub_BC==2)
{


if (yp==1)
for (int i=0;i<bcryn;i++)
	if (nei_loc[bcry[i]][4]>0)
	{
	for (int mi=0; mi<19; mi++)
		{
		m_l[mi]=0;
		for (int mj=0; mj<19; mj++)
			m_l[mi]+=M[mi][mj]*f[nei_loc[bcry[i]][4]][mj];
		}
		m_l[0]=rho_yp;
		for (int mi=0; mi<19; mi++)
			{
			F[bcry[i]][mi]=0;
			for (int mj=0; mj<19; mj++)
				F[bcry[i]][mi]+=MI[mi][mj]*m_l[mj];
			}
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bcry[i]][ks]=feq(ks,rho_yp,u_ls);
		}	


if (yn==1)
for (int i=0;i<bclyn;i++)
	if (nei_loc[bcly[i]][3]>0)
	{
	for (int mi=0; mi<19; mi++)
		{
		m_l[mi]=0;
		for (int mj=0; mj<19; mj++)
			m_l[mi]+=M[mi][mj]*f[nei_loc[bcly[i]][3]][mj];
		}
		m_l[0]=rho_yn;
		for (int mi=0; mi<19; mi++)
			{
			F[bcly[i]][mi]=0;
			for (int mj=0; mj<19; mj++)
				F[bcly[i]][mi]+=MI[mi][mj]*m_l[mj];
			}
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bcly[i]][ks]=feq(ks,rho_yn,u_ls);
		}



if (xp==1)
for (int i=0;i<bcrxn;i++)
	if (nei_loc[bcrx[i]][2]>0)
	{
	for (int mi=0; mi<19; mi++)
		{
		m_l[mi]=0;
		for (int mj=0; mj<19; mj++)
			m_l[mi]+=M[mi][mj]*f[nei_loc[bcrx[i]][2]][mj];
		}
		m_l[0]=rho_xp;
		for (int mi=0; mi<19; mi++)
			{
			F[bcrx[i]][mi]=0;
			for (int mj=0; mj<19; mj++)
				F[bcrx[i]][mi]+=MI[mi][mj]*m_l[mj];
			}
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bcrx[i]][ks]=feq(ks,rho_xp,u_ls);
		}	


if (xn==1)
for (int i=0;i<bclxn;i++)
	if (nei_loc[bclx[i]][1]>0)
	{
	for (int mi=0; mi<19; mi++)
		{
		m_l[mi]=0;
		for (int mj=0; mj<19; mj++)
			m_l[mi]+=M[mi][mj]*f[nei_loc[bclx[i]][1]][mj];
		}
		m_l[0]=rho_xn;
		for (int mi=0; mi<19; mi++)
			{
			F[bclx[i]][mi]=0;
			for (int mj=0; mj<19; mj++)
				F[bclx[i]][mi]+=MI[mi][mj]*m_l[mj];
			}
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bclx[i]][ks]=feq(ks,rho_xn,u_ls);
		}



if (zp==1)
for (int i=0;i<bcrzn;i++)
	if (nei_loc[bcrz[i]][6]>0)
	{
	for (int mi=0; mi<19; mi++)
		{
		m_l[mi]=0;
		for (int mj=0; mj<19; mj++)
			m_l[mi]+=M[mi][mj]*f[nei_loc[bcrz[i]][6]][mj];
		}
		m_l[0]=rho_zp;
		for (int mi=0; mi<19; mi++)
			{
			F[bcrz[i]][mi]=0;
			for (int mj=0; mj<19; mj++)
				F[bcrz[i]][mi]+=MI[mi][mj]*m_l[mj];
			}
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bcrz[i]][ks]=feq(ks,rho_zp,u_ls);
		}	

if (zn==1)
for (int i=0;i<bclzn;i++)
	if (nei_loc[bclz[i]][5]>0)
	{
	for (int mi=0; mi<19; mi++)
		{
		m_l[mi]=0;
		for (int mj=0; mj<19; mj++)
			m_l[mi]+=M[mi][mj]*f[nei_loc[bclz[i]][5]][mj];
		}
		m_l[0]=rho_zn;
		for (int mi=0; mi<19; mi++)
			{
			F[bclz[i]][mi]=0;
			for (int mj=0; mj<19; mj++)
				F[bclz[i]][mi]+=MI[mi][mj]*m_l[mj];
			}
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bclz[i]][ks]=feq(ks,rho_zn,u_ls);
		}




}		

//Non-equilibrium boundary condition
if (Sub_BC==3)
{



if (yp==1)
for (int i=0;i<bcryn;i++)
	if (nei_loc[bcry[i]][4]>0)
	for (int ks=0;ks<Q;ks++)
	{
	u_ls[0]=u[nei_loc[bcry[i]][4]][0];
	u_ls[1]=u[nei_loc[bcry[i]][4]][1];
	u_ls[2]=u[nei_loc[bcry[i]][4]][2];
	F[bcry[i]][ks]=feq(ks,rho_yp,u_ls)+f[nei_loc[bcry[i]][4]][ks]-feq(ks,rho[nei_loc[bcry[i]][4]],u[nei_loc[bcry[i]][4]]);
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bcry[i]][ks]=feq(ks,rho_yp,u_ls);
		}


if (yn==1)
for (int i=0;i<bclyn;i++)
	if (nei_loc[bcly[i]][3]>0)
	for (int ks=0;ks<Q;ks++)
	{
	u_ls[0]=u[nei_loc[bcly[i]][3]][0];
	u_ls[1]=u[nei_loc[bcly[i]][3]][1];
	u_ls[2]=u[nei_loc[bcly[i]][3]][2];
	F[bcly[i]][ks]=feq(ks,rho_yn,u_ls)+f[nei_loc[bcly[i]][4]][ks]-feq(ks,rho[nei_loc[bcly[i]][3]],u[nei_loc[bcly[i]][3]]);
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bcly[i]][ks]=feq(ks,rho_yn,u_ls);
		}

if (xp==1)
for (int i=0;i<bcrxn;i++)
	if (nei_loc[bcrx[i]][2]>0)
	for (int ks=0;ks<Q;ks++)
	{
	u_ls[0]=u[nei_loc[bcrx[i]][2]][0];
	u_ls[1]=u[nei_loc[bcrx[i]][2]][1];
	u_ls[2]=u[nei_loc[bcrx[i]][2]][2];
	F[bcrx[i]][ks]=feq(ks,rho_xp,u_ls)+f[nei_loc[bcrx[i]][4]][ks]-feq(ks,rho[nei_loc[bcrx[i]][2]],u[nei_loc[bcrx[i]][2]]);
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bcrx[i]][ks]=feq(ks,rho_xp,u_ls);
		}


if (xn==1)
for (int i=0;i<bclxn;i++)
	if (nei_loc[bclx[i]][1]>0)
	for (int ks=0;ks<Q;ks++)
	{
	u_ls[0]=u[nei_loc[bclx[i]][1]][0];
	u_ls[1]=u[nei_loc[bclx[i]][1]][1];
	u_ls[2]=u[nei_loc[bclx[i]][1]][2];
	F[bclx[i]][ks]=feq(ks,rho_xn,u_ls)+f[nei_loc[bclx[i]][4]][ks]-feq(ks,rho[nei_loc[bclx[i]][1]],u[nei_loc[bclx[i]][1]]);
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bclx[i]][ks]=feq(ks,rho_xn,u_ls);
		}



if (zp==1)
for (int i=0;i<bcrzn;i++)
	if (nei_loc[bcrz[i]][6]>0)
	for (int ks=0;ks<Q;ks++)
	{
	u_ls[0]=u[nei_loc[bcrz[i]][6]][0];
	u_ls[1]=u[nei_loc[bcrz[i]][6]][1];
	u_ls[2]=u[nei_loc[bcrz[i]][6]][2];
	F[bcrz[i]][ks]=feq(ks,rho_zp,u_ls)+f[nei_loc[bcrz[i]][4]][ks]-feq(ks,rho[nei_loc[bcrz[i]][6]],u[nei_loc[bcrz[i]][6]]);
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bcrz[i]][ks]=feq(ks,rho_zp,u_ls);
		}


if (zn==1)
for (int i=0;i<bclzn;i++)
	if (nei_loc[bclz[i]][5]>0)
	for (int ks=0;ks<Q;ks++)
	{
	u_ls[0]=u[nei_loc[bclz[i]][5]][0];
	u_ls[1]=u[nei_loc[bclz[i]][5]][1];
	u_ls[2]=u[nei_loc[bclz[i]][5]][2];
	F[bclz[i]][ks]=feq(ks,rho_zn,u_ls)+f[nei_loc[bclz[i]][4]][ks]-feq(ks,rho[nei_loc[bclz[i]][5]],u[nei_loc[bclz[i]][5]]);
	}
	else
		for (int ks=0;ks<Q;ks++)
		{
		u_ls[0]=0.0;
		u_ls[1]=0.0;u_ls[2]=0.0;
		F[bclz[i]][ks]=feq(ks,rho_zn,u_ls);
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

		//cout<<temp1<<"	@@@@@@@@@@"<<endl;

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
	
	//MPI_Bcast(&error_in,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//cout<<error_in<<"          wwwwwwwwwwwwwww      "<<rank<<endl;
	
	return(error_in);

}




//OUTPUT SUBROUTAINS:
//ALL THE OUTPUTS ARE TRANSFERED TO PROCESSOR 0, AND EXPORT TO DAT FILE BY PROCESSOR 0


void Geometry(int*** Solid)	
{	
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	ostringstream name;
	name<<outputfile<<"LBM_Geometry"<<".vtk";
	if (rank==0)
	{

	ofstream out;
	out.open(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Geometry"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<NX+1<<"         "<<NY+1<<"         "<<NZ+1<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<(NX+1)*(NY+1)*(NZ+1)<<endl;
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;

	

	for(int k=0;k<NZ+1;k++)
		for(int j=0; j<NY+1; j++)
			for(int i=0;i<NX+1;i++)
			if (Solid[i][j][k]>0)
				out<<0<<" ";
			else
				out<<1<<" ";


	out.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	
	
		
}


void Geometry_Par(int*** Solid)	
{	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	ostringstream name;
	name<<outputfile<<"LBM_Geometry_Mesh"<<".vtk";
	if (rank==0)
	{

	ofstream out;
	out.open(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Geometry"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<NX+1<<"         "<<NY+1<<"         "<<NZ+1<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<(NX+1)*(NY+1)*(NZ+1)<<endl;
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;

	

	for(int k=0;k<NZ+1;k++)
		for(int j=0; j<NY+1; j++)
			for(int i=0;i<NX+1;i++)
			out<<Solid[i][j][k]<<" ";


	out.close();
	}
	MPI_Barrier(MPI_COMM_WORLD);
	
	
}


void output_velocity(int m,double* rho,double** u,int MirX,int MirY,int MirZ,int mir,int*** Solid)	
{
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	const int mpi_size=MPI :: COMM_WORLD . Get_size ();
	int procind=rank+1;
	int procn=mpi_size;

	int tmpsum[mpi_size];
	for (int i=0;i<mpi_size;i++)
		tmpsum[i]=3;	

	const int root_rank=0;
	
	double rho_0=1.0;
	
	
	MPI_Status status;
	MPI_Request request;

	double* rbuf_v;



	int nx_g[mpi_size];
	int disp[mpi_size];
	
	for (int i=0;i<mpi_size;i++)
		nx_g[i]=(sumss[i+1]+1)*3;


	

		disp[0]=0;
	
		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
		
	
	if (rank==root_rank)
		rbuf_v = new double[disp[mpi_size-1]+nx_g[mpi_size-1]];

	
	int NX0=NX+1;
	int NY0=NY+1;
	int NZ0=NZ+1;
	MPI_Gatherv(u[0],nx_g[rank],MPI_DOUBLE,rbuf_v,nx_g,disp,MPI_DOUBLE,root_rank,MPI_COMM_WORLD);



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
	out<<"DIMENSIONS         "<<NX0<<"         "<<NY0<<"         "<<NZ0<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;

	out<<"POINT_DATA     "<<NX0*NY0*NZ0<<endl;
	out<<"VECTORS sample_vectors double"<<endl;
	out<<endl;

	for(int k=0 ; k<NZ0 ; k++)			
	         for(int j=0 ; j<NY0 ; j++)
	                for(int i=0 ; i<NX0 ; i++)
			if (Solid[i][j][k]==0)
				out<<0<<" "<<0<<" "<<0<<endl;
			else
				{
				out<<rbuf_v[disp[Solid[i][j][k]-1]+tmpsum[Solid[i][j][k]-1]]<<" ";
				out<<rbuf_v[disp[Solid[i][j][k]-1]+tmpsum[Solid[i][j][k]-1]+1]<<" ";
				out<<rbuf_v[disp[Solid[i][j][k]-1]+tmpsum[Solid[i][j][k]-1]+2]<<endl;
				tmpsum[Solid[i][j][k]-1]+=3;
				}
	
	out.close();

	}
        MPI_Barrier(MPI_COMM_WORLD);
	
	if (rank==root_rank)
	delete [] rbuf_v;
		
}

void output_velocity_compact(int m,double* rho,double** u,int MirX,int MirY,int MirZ,int mir,int*** Solid)	
{
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	const int mpi_size=MPI :: COMM_WORLD . Get_size ();
	int procind=rank+1;
	int procn=mpi_size;

	int tmpsum[mpi_size];
	for (int i=0;i<mpi_size;i++)
		tmpsum[i]=3;	

	const int root_rank=0;
	
	double rho_0=1.0;
	
	
	MPI_Status status;
	MPI_Request request;

	double* rbuf_v;
	double* rbuf_v2;



	int nx_g[mpi_size];
	int disp[mpi_size];
	
	for (int i=0;i<mpi_size;i++)
		nx_g[i]=(sumss[i+1]+1)*3;


	

		disp[0]=0;
	
		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
		
	
	if (rank==root_rank)
		{
		rbuf_v = new double[disp[mpi_size-1]+nx_g[mpi_size-1]];
		rbuf_v2 = new double[disp[mpi_size-1]+nx_g[mpi_size-1]];
		}

	
	int NX0=NX+1;
	int NY0=NY+1;
	int NZ0=NZ+1;
	MPI_Gatherv(u[0],nx_g[rank],MPI_DOUBLE,rbuf_v,nx_g,disp,MPI_DOUBLE,root_rank,MPI_COMM_WORLD);

	int sumtmp=0;

	ostringstream name;
	name<<"vel.bin";
	if (rank==root_rank)
	{

	ofstream out;
	out.open(name.str().c_str());
	//out<<"# vtk DataFile Version 2.0"<<endl;
	//out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Velocity"<<endl;
	//out<<"ASCII"<<endl;
	//out<<"DATASET STRUCTURED_POINTS"<<endl;
	//out<<"DIMENSIONS         "<<NX0<<"         "<<NY0<<"         "<<NZ0<<endl;
	//out<<"ORIGIN 0 0 0"<<endl;
	//out<<"SPACING 1 1 1"<<endl;

	//out<<"POINT_DATA     "<<NX0*NY0*NZ0<<endl;
	//out<<"VECTORS sample_vectors double"<<endl;
	//out<<endl;

	

	for(int k=0 ; k<NZ0 ; k++)			
	         for(int j=0 ; j<NY0 ; j++)
	                for(int i=0 ; i<NX0 ; i++)
			if (Solid[i][j][k]>0)
				{
				rbuf_v2[sumtmp]=rbuf_v[disp[Solid[i][j][k]-1]+tmpsum[Solid[i][j][k]-1]];
				rbuf_v2[sumtmp+1]=rbuf_v[disp[Solid[i][j][k]-1]+tmpsum[Solid[i][j][k]-1]+1];
				rbuf_v2[sumtmp+2]=rbuf_v[disp[Solid[i][j][k]-1]+tmpsum[Solid[i][j][k]-1]+2];
				tmpsum[Solid[i][j][k]-1]+=3;
				sumtmp+=3;
				}
	out.write((char *)(&rbuf_v2[0]), sizeof(double)*(disp[mpi_size-1]+nx_g[mpi_size-1]));
	out.close();

	}
        MPI_Barrier(MPI_COMM_WORLD);
	
        cout<<sumtmp<<"         update export sum of velocity"<<endl;
        
	if (rank==root_rank)
		{
		delete [] rbuf_v;
		delete [] rbuf_v2;
		}
		
}

void output_density(int m,double* rho,int MirX,int MirY,int MirZ,int mir,int*** Solid)	
{
    
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	int procind=rank+1;
	int procn=mpi_size;

	int tmpsum[mpi_size];
	for (int i=0;i<mpi_size;i++)
		tmpsum[i]=1;	

	const int root_rank=0;
	
	double rho_0=1.0;
	
	
	MPI_Status status;
	MPI_Request request;

	double* rbuf_v;



	int nx_g[mpi_size];
	int disp[mpi_size];
	
	for (int i=0;i<mpi_size;i++)
		nx_g[i]=sumss[i+1]+1;


	

		disp[0]=0;
	
		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
		
	
	if (rank==root_rank)
		rbuf_v = new double[disp[mpi_size-1]+nx_g[mpi_size-1]];

	
	
	
	int NX0=NX+1;
	int NY0=NY+1;
	int NZ0=NZ+1;
	MPI_Gatherv(rho,nx_g[rank],MPI_DOUBLE,rbuf_v,nx_g,disp,MPI_DOUBLE,root_rank,MPI_COMM_WORLD);

        


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
	out<<"DIMENSIONS         "<<NX0<<"         "<<NY0<<"         "<<NZ0<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<NX0*NY0*NZ0<<endl;
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;

	for(int k=0 ; k<NZ0 ; k++)			
	         for(int j=0 ; j<NY0 ; j++)
	                for(int i=0 ; i<NX0 ; i++)
			if (Solid[i][j][k]==0)
				out<<rho0<<" ";
			else
				{
				out<<rbuf_v[disp[Solid[i][j][k]-1]+tmpsum[Solid[i][j][k]-1]]<<" ";
				tmpsum[Solid[i][j][k]-1]++;
				}
	
	out.close();

	}
        MPI_Barrier(MPI_COMM_WORLD);
	
	if (rank==root_rank)
	delete [] rbuf_v;

}


void Backup(int m,double* rho,double** u, double** f)
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
	
	
	
	ostringstream name4;
	name4<<outputfile<<"LBM_checkpoint_f_"<<m<<"."<<rank<<".bin_input";
	//ofstream out;
	out.open(name4.str().c_str());
	
	
        
        out.write((char *)(&f[0][0]), sizeof(double)*(Count+1)*19);
                
        
        //cout<<"fffssssssssss     "<<f[6][0]<<endl;
        out.close();
        
	
	
}

double Comput_Perm(double** u,double* Permia,int PerDIr,int* SupInv)
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	//int ls_int=0;
	
	
	int si,sj,sm;
	
	double avex,avey,avez;
	
	




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
			dp=abs(p_xp-p_xn)*c_s2/(per_xp-per_xn+1)/dx;break;
		case 2:
			dp=abs(p_yp-p_yn)*c_s2/(per_yp-per_yn+1)/dx;break;
		case 3:
			dp=abs(p_zp-p_zn)*c_s2/(per_zp-per_zn+1)/dx;break;
		default:
			dp=abs(p_xp-p_xn)*c_s2/(per_xp-per_xn+1)/dx;
		}

		
	if ((par_per_x-1)*(par_per_y-1)*(par_per_z-1)==0)	
	for (int i=1;i<=Count;i++)
	{
		si=(int)(coor[i]/((NY+1)*(NZ+1)));
		sj=(int)((coor[i]%((NY+1)*(NZ+1)))/(NZ+1));
		sm=(int)(coor[i]%(NZ+1)); 
		//si+=disp[rank];
		//if (rank==1)
		//cout<<rank<<"  "<<si<<" "<<sj<<" "<<sm<<endl;
		//cout<<si<<"  "<<sj<<"  "<<sm<<endl;
		//cout<<rank<<"        "<<per_yp<<"        "<<per_yn<<endl;
		//ls_int++;
		if ((si>=per_xn) and (si<=per_xp) and (sj>=per_yn) and (sj<=per_yp) and (sm>=per_zn) and (sm<=per_zp))
		{
		  //ls_int++;     
	        Q[0]+=u[i][0];
		Q[1]+=u[i][1];
		Q[2]+=u[i][2];
		}
		  //   else
		     //           cout<<si<<"  "<<sj<<"  "<<sm<<endl;


	}
	else
	for (int i=1;i<=Count;i++)
	       {
		Q[0]+=u[i][0];
		Q[1]+=u[i][1];
		Q[2]+=u[i][2];
		} 

	

	MPI_Barrier(MPI_COMM_WORLD);

	     //   cout<<ls_int<<"        "<<Q[0]<<"        "<<Q[1]<<"        "<<Q[2]<<"        "<<rank<<endl;

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

		//cout<<(per_xp-per_xn+1)*(per_yp-per_yn+1)*(per_zp-per_zn+1)*porosity<<endl;
		//cout<<"SUM= "<<Q[0]<<"        "<<Q[1]<<"        "<<Q[2]<<endl;
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
	name4<<pfix<<"LBM_checkpoint_f_"<<mode_backup_ini<<"."<<rank<<".bin_input";
	
 	ostringstream name2;
	name2<<pfix<<"LBM_checkpoint_velocity_"<<mode_backup_ini<<"."<<rank<<".bin_input";
	
	ostringstream name;
	name<<pfix<<"LBM_checkpoint_rho_"<<mode_backup_ini<<"."<<rank<<".bin_input";
	
	// cout<<name2.str().c_str()<<"         mmmmmmmmm"<<endl;
	 
	fstream fin;
	fin.open(name.str().c_str(),ios::in);
	if (fin.fail())
	        {
	        cout<<"\n file open error on " << name.str().c_str()<<endl;
	        exit(-1);
	        }
	
	        fin.read((char *)(&rho[0]), sizeof(double)*(Count+1));
	        
        	//for(int i=1;i<=Count;i++)
        	   //     fin >> rho[i];
  
       fin.close();
       
   
	fin.open(name2.str().c_str(),ios::in);
	if (fin.fail())
	        {
	        cout<<"\n file open error on " << name2.str().c_str()<<endl;
	        exit(-1);
	        }
	        
	fin.read((char *)(&u[0][0]), sizeof(double)*(Count+1)*3);
        //	for(int i=1;i<=Count;i++)
        //	        fin >> u[i][0] >> u[i][1] >> u[i][2];
  
       fin.close();
       
      
       fin.open(name4.str().c_str(),ios::in);
       if (fin.fail())
	        {
	        cout<<"\n file open error on " << name4.str().c_str()<<endl;
	        exit(-1);
	        }
	        
	        
       fin.read((char *)(&f[0][0]), sizeof(double)*(Count+1)*19);
        	//for(int i=1;i<=Count;i++)
        	//        fin >> f[i][0] >> f[i][1] >> f[i][2] >> f[i][3] >> f[i][4] >> f[i][5] >>f[i][6] >> f[i][7] >> f[i][8] >> f[i][9] >> f[i][10] >> f[i][11] >> f[i][12] >> f[i][13] >> f[i][14] >> f[i][15] >> f[i][16] >>f[i][17] >> f[i][18];
  
       fin.close();
       
      
	//cout<<f[6][0]<<"       aaaaaa"<<endl;

	 	
}



void output_velocity_for_solute(int m,double* rho,double** u,int MirX,int MirY,int MirZ,int mir,int*** Solid)	
{
  	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	

	double* lsu= new double[Count*3];
		for (int i=1;i<=Count;i++)
			for (int j=0;j<3;j++)
			lsu[(i-1)*3+j]=u[i][j];


	ostringstream name;
	name<<outputfile<<"Velocity_for_solute_"<<m<<"."<<rank<<".bin";
	ofstream out;
	out.open(name.str().c_str());


        
			        out.write((char *)lsu,sizeof(double)*Count*3);
	
	out.close();
		
	delete [] lsu;
	

		
}

/*
void Partition_Solid(int*** Solid)
{


int rank = MPI :: COMM_WORLD . Get_rank ();
int para_size=MPI :: COMM_WORLD . Get_size ();
	MPI_Comm comm;
	MPI_Comm_dup(MPI_COMM_WORLD, &comm);

	int sum=0;
	int sum3,sum_rec;
	int ii,jj,kk;
	int mesh_par=para_size;

	int nx=NX+1;
	int ny=NY+1;
	int nz=NZ+1;

	int* rbuf_result;

	for(int k=0 ; k<nz ; k++)				
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)				
	

	
		{	
			
			
			
			
			if (Solid[i][j][k] == 0)	{sum++;Solid[i][j][k] = sum-1;}
			else
				{Solid[i][j][k] = -1;}
			
		
			
			
		}
	
	sum_rec=sum;
	



	idx_t* vtxdist=NULL;
	idx_t *xadj, *adjncy;
        vtxdist = new idx_t[para_size+1];
		


        int ls1,ls2;
        int *proc_size;
        proc_size = new int[para_size];	
        
        ls1=sum_rec%para_size;
        ls2=(int)sum_rec/para_size;
        
        for (int i=0;i<para_size;i++)
                proc_size[i]=ls2;
        for(int i=0;i<ls1;i++)
                proc_size[i]++;
        
    


        
        vtxdist[0]=0;
        for (int i=1;i<=para_size;i++)
        {        vtxdist[i]=vtxdist[i-1]+proc_size[i-1];
                
        }
        
        
        
        xadj = new idx_t[proc_size[rank]+1];
        sum=0;
        
        
        for(int k=0 ; k<nz ; k++)			
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)	
	{
	        if ((Solid[i][j][k]>=vtxdist[rank]) and (Solid[i][j][k]<vtxdist[rank+1]))
	                for (int ls=1;ls<19;ls++)
		{
		ii=i+e[ls][0];
		jj=j+e[ls][1];
		kk=k+e[ls][2];	
		
		
	
		
		if ((ii>=0) and (ii<nx) and (jj>=0) and (jj<ny) and (kk>=0) and (kk<nz) and (Solid[ii][jj][kk]>=0))
			sum++;
		}
		
	}
        
      adjncy = new idx_t[sum];
      xadj[0]=0;
      
      sum=0;sum3=0;
      for(int k=0 ; k<nz ; k++)			
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)	
	{
	        if ((Solid[i][j][k]>=vtxdist[rank]) and (Solid[i][j][k]<vtxdist[rank+1]))
	        {
	                for (int ls=1;ls<19;ls++)
	                {
		ii=i+e[ls][0];
		jj=j+e[ls][1];
		kk=k+e[ls][2];	
		
	
		
		if ((ii>=0) and (ii<nx) and (jj>=0) and (jj<ny) and (kk>=0) and (kk<nz) and (Solid[ii][jj][kk]>=0))
		        {adjncy[sum]=Solid[ii][jj][kk];sum++;}
		        }
		        
		        
		        
		sum3++;
		xadj[sum3]=sum;
		
		
		}
		
	}


	

	

	idx_t *vwgt=NULL;
	
	idx_t *adjwgt=NULL;
	idx_t wgtflag=0;
	idx_t numflag=0;
	idx_t ncon=1;
	idx_t nparts=mesh_par;
	real_t *tpwgts;
	real_t *ubvec;
	idx_t options[10];
	idx_t edgecut;
	idx_t *part;
	part = new idx_t[proc_size[rank]];
	//vwgt = new idx_t[proc_size[rank]*ncon];

	
	tpwgts = new real_t[ncon*nparts];
	for (int i=0;i<ncon*nparts;i++)
	        tpwgts[i]=1.0/(real_t)nparts;
	
	ubvec = new real_t[ncon];
	for (int i=0;i<ncon;i++)
	        ubvec[i]=1.05;
	
	//vwgt = new idx_t[proc_size[rank]];
	//	for (int i=0;i<proc_size[rank];i++)
	//	vwgt[i]=10;

	//adjwgt = new idx_t[sum];
	//	for (int i=0;i<sum;i++)
	//	adjwgt[i]=1;

	
	options[0] = 0;
	
	if (rank==0)
	cout<<"Start Partition "<<endl;
	


	MPI_Barrier(MPI_COMM_WORLD);
	

	
      ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt,adjwgt, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);

	if (rank==0)
	cout<<"Partition Complete"<<endl;

	int* part_int;
	int* disp;
	int* sta_parts;

	part_int = new int[proc_size[rank]];

	for (int i=0;i<proc_size[rank];i++)
		part_int[i]=part[i];

	if (rank==0)
		rbuf_result = new int[sum_rec];
	
	disp = new int[para_size];
	for (int i=0;i<para_size;i++)
		disp[i]=vtxdist[i];



	MPI_Gatherv(part_int,proc_size[rank],MPI_INT,rbuf_result,proc_size,disp,MPI_INT,0,MPI_COMM_WORLD);

	if (rank==0)
	{
	sta_parts = new int[mesh_par];
	for (int i=0;i<mesh_par;i++)
		sta_parts[i]=0;

	sum=0;
	for (int k=0;k<nz;k++)
	for (int j=0;j<ny;j++)
	for (int i=0;i<nx;i++)
		if (Solid[i][j][k]>=0)
			{Solid[i][j][k]=rbuf_result[sum]+1;sta_parts[rbuf_result[sum]]++;sum++;}
		else	
			{Solid[i][j][k]=0;}

	
	for (int i=0;i<mesh_par;i++)
		cout<<sta_parts[i]<<"	"<<i<<"		"<<sta_parts[i]-(int)sum/mesh_par<<"	"<<(double)(sta_parts[i]-(int)sum/mesh_par)/sta_parts[i]<<endl;
	}

	MPI_Bcast(Solid[0][0],nx*ny*nz,MPI_INT,0,MPI_COMM_WORLD);

	
	delete [] proc_size;
	delete [] xadj;
	delete [] adjncy;
	delete [] part;
	delete [] tpwgts;
	delete [] ubvec;
	delete [] part_int;
	delete [] disp;

	if (rank==0)
		delete [] rbuf_result;

	
	
}
*/

void Partition_Solid_SELF(int*** Solid)
{
	int nx=NX+1;
	int ny=NY+1;
	int nz=NZ+1;
	cout<<endl;
	cout<<"MESH PARTITION INITIALIZATION START"<<endl;

int sum=0;
int mesh_par=MPI :: COMM_WORLD . Get_size ();
for(int k=0 ; k<nz ; k++)			
	for(int j=0 ; j<ny ; j++)	
	for(int i=0 ; i<nx ; i++)

		if (Solid[i][j][k]==0)
			{sum++;}
		else
			{Solid[i][j][k]=-1;}


int nxref,nyref,nzref;
int divnumori;
divnumori=mesh_par;
int divnum;
int oddval=0;
int evennum=0;
double oddpor;
int dir;
int dint;
int sumin,numgeonum;
int* sum_loc;


int *nnx,*nny,*nnz,*npx,*npy,*npz;
nnx=new int[divnumori];
nny=new int[divnumori];
nnz=new int[divnumori];
npx=new int[divnumori];
npy=new int[divnumori];
npz=new int[divnumori];




cout<<divnumori<<endl;
divnum=divnumori;
while (divnum%2==0)
        {evennum++;divnum=divnum/2;}
oddval=divnum;
divnum=divnumori;

cout<<evennum<<"         "<<oddval<<endl;
 

	sum_loc = new int[divnum+1];
	for (int i=0;i<=divnum;i++)
	        sum_loc[i]=0;
		
	
	nxref=nx;nyref=ny;nzref=nz;
	oddpor=sum/oddval;
	if (oddval>1)
	        {
	                //oddpor=sum/oddval;
	                //cout<<oddpor<<endl;
	               if (nx>ny)
	                       if (nx>nz)
	                       dir=1;
	                       else
	                               dir=3;
                       else
	                               if (ny>nz)
	                               dir=2;
	                               else
	                               dir=3;
	                   //cout<<dir<<endl;
	                if (dir==1)
	                        for (int sn=1;sn<=oddval;sn++)
	                        {
	                                sumin=0;
	                                dint=0;
	                                while ((sumin<oddpor*sn) and (dint<nx))
	                                {
	                                for (int j=0;j<ny;j++)
	                                        for (int k=0;k<nz;k++)
	                                        {
	                                                if ((Solid[dint][j][k]<sn) and (Solid[dint][j][k]>=0))
	                                                        {
	                                                                sumin++;
	                                                                if (Solid[dint][j][k]==0)
	                                                                        Solid[dint][j][k]=sn;
	                                                        }
	                                                
	                                        }
	                                        dint++;
	                                }
	                                
	                             nxref=nx/oddval;   
	                                
	                                
	                        }
	                        
	                        
	                   if (dir==2)
	                        for (int sn=1;sn<=oddval;sn++)
	                        {
	                                sumin=0;
	                                dint=0;
	                                while ((sumin<oddpor*sn) and (dint<ny))
	                                {
	                                for (int i=0;i<nx;i++)
	                                        for (int k=0;k<nz;k++)
	                                        {
	                                                if ((Solid[i][dint][k]<sn) and (Solid[i][dint][k]>=0))
	                                                        {
	                                                                sumin++;
	                                                                if (Solid[i][dint][k]==0)
	                                                                        Solid[i][dint][k]=sn;
	                                                        }
	                                                
	                                        }
	                                        dint++;
	                                }
	                                
	                                
	                              nyref=ny/oddval;     
	                                
	                        }
	                             
	                 if (dir==3)
	                        for (int sn=1;sn<=oddval;sn++)
	                        {
	                                sumin=0;
	                                dint=0;
	                                while ((sumin<oddpor*sn) and (dint<nz))
	                                {
	                                for (int i=0;i<nx;i++)
	                                        for (int j=0;j<ny;j++)
	                                        {
	                                                if ((Solid[i][j][dint]<sn) and (Solid[i][j][dint]>=0))
	                                                        {
	                                                                sumin++;
	                                                                if (Solid[i][j][dint]==0)
	                                                                        Solid[i][j][dint]=sn;
	                                                        }
	                                                
	                                        }
	                                        dint++;
	                                }
	                                
	                           nzref=nz/oddval;     
	                                
	                                
	                        }       
	                
	        }
	        else
	                {
	                  for (int k=0;k<nz;k++)
	              for (int j=0;j<ny;j++)
	              for (int i=0;i<nx;i++)
	                      if (Solid[i][j][k]==0)
	                              Solid[i][j][k]=1;
	                      else
	                              Solid[i][j][k]=-1;
	                }
	
	
	      for (int k=0;k<nz;k++)
	              for (int j=0;j<ny;j++)
	              for (int i=0;i<nx;i++)
	             if (Solid[i][j][k]<0)
	             Solid[i][j][k]=0;
	
	
	
	numgeonum=oddval;
	cout<<numgeonum<<endl;
	//cout<<"@@@@@@@@@@@"<<endl;
	for (int sn=1;sn<=evennum;sn++)
	        {
	                
	                oddpor=oddpor/2;
	                for (int i=0;i<=divnum;i++)
	                        sum_loc[i]=0;
	                
	          if (nxref>nyref)
	                       if (nxref>nzref)
	                       dir=1;
	                       else
	                               dir=3;
                       else
	                               if (nyref>nzref)
	                               dir=2;
	                               else
	                               dir=3;      
	                
	             //-----------------------------
	            // cout<<dir<<"            eeeeeeeeeeeee"<<endl;
	             //cout<<Solid[88][0][0]<<"        afadsfasdf         "<<endl;
	             if (dir==1)
	                     {
	                         for (int i=0;i<nx;i++)
	                                 {
	                                         for (int j=0;j<ny;j++)
	                                         for (int k=0;k<nz;k++)
	                                         if (Solid[i][j][k]>0)
	                                                 if (sum_loc[Solid[i][j][k]]>=0)
	                                                 {sum_loc[Solid[i][j][k]]++;Solid[i][j][k]*=-1;}
	                                         
	             
	                                   
	                             
	                               for (int si=0;si<=divnum;si++)
	                                       if (sum_loc[si]>oddpor)
	                                               sum_loc[si]=-1;
	                                       
	                                 }
	                       nxref=nxref/2;          
	                     }
	                     
	            if (dir==2)
	                     {
	                         for (int j=0;j<ny;j++)
	                                 {
	                                         for (int i=0;i<nx;i++)
	                                         for (int k=0;k<nz;k++)
	                                         if (Solid[i][j][k]>0)
	                                                 if (sum_loc[Solid[i][j][k]]>=0)
	                                                 {sum_loc[Solid[i][j][k]]++;Solid[i][j][k]*=-1;}
	                                                 //cout<<Solid[i][j][k]<<"                bbbbbbbbb        "<<i<<"        "<<j<<"        "<<k<<endl;
	                                         
	                                         
	                             
	                               for (int siz=1;siz<=divnum;siz++)               
	                                       if (sum_loc[siz]>oddpor)
	                                               sum_loc[siz]=-1;
	                                 
	                                       
	                                 }
	                       nyref=nyref/2;          
	                     }         
	                     
	                if (dir==3)
	                     {
	                         for (int k=0;k<nz;k++)
	                                 {
	                                         for (int i=0;i<nx;i++)
	                                         for (int j=0;j<ny;j++)
	                                         if (Solid[i][j][k]>0)
	                                                 if (sum_loc[Solid[i][j][k]]>=0)
	                                                 {sum_loc[Solid[i][j][k]]++;Solid[i][j][k]*=-1;}
	                                
	                             
	                               for (int si=0;si<=divnum;si++)
	                                       if (sum_loc[si]>oddpor)
	                                               sum_loc[si]=-1;
	                                       
	                                 }
	                      nzref=nzref/2;           
	                     }              
	         numgeonum*=2;
	         cout<<numgeonum<<endl;

	         for (int k=0;k<nz;k++)
	              for (int j=0;j<ny;j++)
	              for (int i=0;i<nx;i++)
	              {
	                      if (Solid[i][j][k]>0)
	                                   Solid[i][j][k]=Solid[i][j][k]*2;
	                           if (Solid[i][j][k]<0)
	                                          Solid[i][j][k]=-Solid[i][j][k]*2-1;
	                
	                
	                
	        }
	
	  }
	    
	
	//============decomposition complete=======================
	
	//==================================================
	
	

	






	

	

}
