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


void Read_Rock(int***,double*, double*,double*,char[128], char[128]);

void tests();

void init_Sparse_read_rock_parallel(int*,int*);

void init(double*, double**, double**, double*, double*, double*, double*, double*, int*);

void collision(double*,double** ,double** ,double** ,double*, double*, double*, double*, double*,int* ,int***,int*, int*);

void comput_macro_variables( double* ,double**,double** ,double** ,double**, double*, double*, int* ,int***);

double Error(double** ,double** ,double*, double*);

void boundary_velocity(int,double,int,double,int,double,int,double,int,double,int,double,double**,double**,double*,double**,int***,double*);

void boundary_pressure(int,double,int,double,int,double,int,double,int,double,int,double,double**,double**,double**,double*,int***,double*);

void output_velocity(int ,double* ,double** ,int ,int ,int ,int ,int*** );

void output_density(int ,double* ,int ,int ,int ,int ,int*** );	

void Geometry(int*** );	

void output_velocity_b(int ,double* ,double** ,int ,int ,int ,int ,int*** );

void output_density_b(int ,double* ,int ,int ,int ,int ,int*** );	

void Geometry_b(int*** );

void Backup(int ,double* ,double**, double**);

double S[19];

void Comput_MI(double[19][19], double[19][19]);

int inverse(mat &a);

double feq(int,double, double[3],double);

void Suppliment(int*,int***);

void Parallelize_Geometry();

void Backup_init(double* rho, double** u, double** f, double*, double*,double*, double*, double*,char[128], char[128],char[128],int*);

int e[19][3]=
{{0,0,0},{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};

double elat[19][3]=
{{0,0,0},{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};

double w[19]={1.0/3.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};

int LR[19]={0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15};
int FRP[19]={0,0,0,0,0,0,0,1,0,2,0,0,0,0,0,3,0,4,0};
int FLN[19]={0,0,0,0,0,0,0,0,1,0,2,0,0,0,0,0,3,0,4};
int RP[5]={1,7,9,15,17};
int LN[5]={2,8,10,16,18};



int n,nx_l,n_max,in_BC,freRe,freDe,freVe,Par_Geo,Par_nx,Par_ny,Par_nz,lattice_v;
int Zoom;


int wr_per,pre_xp,pre_xn,pre_yp,pre_yn,pre_zp,pre_zn,fre_backup;
int vel_xp,vel_xn,vel_yp,vel_yn,vel_zp,vel_zn,Sub_BC,Out_Mode,mode_backup_ini;
double in_vis,p_xp,p_xn,p_yp,p_yn,p_zp,p_zn;
double inivx,inivy,inivz,v_xp,v_xn,v_yp,v_yn,v_zp,v_zn;
double error_perm,F_epsilon,c_s,c_s2,dx_input,dt_input,lat_c;

char outputfile[128]="./";
int NCHAR=128;
	char     filename[128], perfile[128], dummy[128+1], backup_rho[128], backup_velocity[128],backup_f[128];
	int      dummyInt;

	
int*** Solid;
	float* Per_Int;
	float* Por_Int;
	
	
int main(int argc , char *argv [])
{	

MPI :: Init (argc , argv );
MPI_Status status ;

double start , finish,remain;

int rank = MPI :: COMM_WORLD . Get_rank ();
int para_size=MPI :: COMM_WORLD . Get_size ();

int dif,ts,th,tm;

int tse,the,tme;
double elaps;

 

	MPI_Barrier(MPI_COMM_WORLD);
	start = MPI_Wtime();

	if (rank==0)
	{
	ifstream fin(argv[1]);
	                                                        fin.getline(dummy, NCHAR);
	fin >> filename;				fin.getline(dummy, NCHAR);
	fin >> perfile;					fin.getline(dummy, NCHAR);
	fin >> NX >> NY >> NZ;				fin.getline(dummy, NCHAR);
	fin >> n_max;					fin.getline(dummy, NCHAR);
	fin >> reso;					fin.getline(dummy, NCHAR);
	fin >> in_BC;					fin.getline(dummy, NCHAR);
	fin >> F_epsilon;				fin.getline(dummy, NCHAR);
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
	fin >> freRe;					fin.getline(dummy, NCHAR);
	fin >> Out_Mode;				fin.getline(dummy, NCHAR);
	fin >> freVe;					fin.getline(dummy, NCHAR);
	fin >> freDe;					fin.getline(dummy, NCHAR);
							fin.getline(dummy, NCHAR);
	fin >> lattice_v >> dx_input >> dt_input;	fin.getline(dummy, NCHAR);
	fin >> outputfile;				fin.getline(dummy, NCHAR);
	fin >> Sub_BC;					fin.getline(dummy, NCHAR);
	fin >> fre_backup;                        	fin.getline(dummy, NCHAR);
	fin >>mode_backup_ini;                		fin.getline(dummy, NCHAR);
	fin >> backup_rho;                        	fin.getline(dummy, NCHAR);
	fin >> backup_velocity;                		fin.getline(dummy, NCHAR);
	fin >> backup_f;                        	fin.getline(dummy, NCHAR);
	
	
	fin.close();
	
	cout<<backup_f<<endl;
	
	
	NX=NX-1;NY=NY-1;NZ=NZ-1;
	}

	
	
	
	MPI_Bcast(&filename,128,MPI_CHAR,0,MPI_COMM_WORLD);MPI_Bcast(&perfile,128,MPI_CHAR,0,MPI_COMM_WORLD);
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
	MPI_Bcast(&inivz,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&backup_f,128,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&freRe,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&freVe,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&freDe,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&backup_velocity,128,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&mir,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&in_vis,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Zoom,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&outputfile,128,MPI_CHAR,0,MPI_COMM_WORLD);
	MPI_Bcast(&fre_backup,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&mode_backup_ini,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&Sub_BC,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&Out_Mode,1,MPI_INT,0,MPI_COMM_WORLD);
	MPI_Bcast(&backup_rho,128,MPI_CHAR,0,MPI_COMM_WORLD);MPI_Bcast(&F_epsilon,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	
      	MPI_Bcast(&lattice_v,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&dx_input,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&dt_input,1,MPI_DOUBLE,0,MPI_COMM_WORLD);


int U_max_ref=0;

mirX=0;
mirY=0;
mirZ=0;
mir=1;
Zoom=1;

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


	nx_l=(int)((NX+1)/para_size);
	dif=(NX+1)-nx_l*para_size;
	
	if (rank>para_size-1-dif)
		nx_l+=1;



	double* rho;
	double** u;
	double**f;
	double**F;
	double**u0;
	int* SupInv;
	double* forcex;
	double* forcey;
	double* forcez;
	
	


	double* PerC;
	double* PorC;


	
	int*  Sl;
	int*  Sr;

	 Parallelize_Geometry();
	
	


	Sl = new int[(NY+1)*(NZ+1)];
	Sr = new int[(NY+1)*(NZ+1)];


//	Per_Int = new double[((NX+1)/Zoom)*((NY+1)/Zoom)*((NZ+1)/Zoom)];
//	Por_Int = new double[((NX+1)/Zoom)*((NY+1)/Zoom)*((NZ+1)/Zoom)];


	

//	Solid = new int**[nx_l];
//	for (int i=0;i<nx_l;i++)
//		{
//		Solid[i] = new int*[NY+1];
//			for (int j=0;j<=NY;j++)
//			Solid[i][j]= new int[NZ+1];
//		}


	
	init_Sparse_read_rock_parallel(Sl,Sr);
	
	
	

	

	//***************************************************
	//WARRING: SPARSE MATRIX STARTS FROM INDEX 1 NOT 0!!!
	//***************************************************

	
	rho = new double[Count+1];
	forcex = new double[Count+1];
	forcey = new double[Count+1];
	forcez = new double[Count+1];
	
	PerC = new double[Count+1];
	PorC = new double[Count+1];

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
	init(rho,u,f,PerC, PorC, forcex, forcey,forcez, SupInv);
else
        Backup_init(rho,u,f,PerC,PorC, forcex,forcey, forcez,backup_rho,backup_velocity,backup_f,SupInv);

if (rank==0)
		cout<<"Porosity= "<<porosity<<endl;


	delete [] Per_Int;
	delete [] Por_Int;



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
	
	
	collision(rho,u,f,F,PerC, PorC,forcex,forcey,forcez,SupInv,Solid,Sl,Sr);

	
	if ((1-pre_xp)*(1-pre_xn)*(1-pre_yp)*(1-pre_yn)*(1-pre_zp)*(1-pre_zn)==0)
		boundary_pressure(pre_xp,p_xp,pre_xn,p_xn,pre_yp,p_yp,pre_yn,p_yn,pre_zp,p_zp,pre_zn,p_zn,f,F,u,rho,Solid,PorC);

	
	if ((1-vel_xp)*(1-vel_xn)*(1-vel_yp)*(1-vel_yn)*(1-vel_zp)*(1-vel_zn)==0)
		boundary_velocity(vel_xp,v_xp,vel_xn,v_xn,vel_yp,v_yp,vel_yn,v_yn,vel_zp,v_zp,vel_zn,v_zn,f,F,rho,u,Solid,PorC);

  		comput_macro_variables(rho,u,u0,f,F,PerC, PorC,SupInv,Solid); 

	
	
	if(n%freRe==0)
		{       
			
			

			if (rank==0)
			{
			 ofstream fin(FileName,ios::out);       
			 fin<<"The"<<n-freRe<<"th computation result:"<<endl;
			Re=u_ave*(NY+1)/(1.0/3.0*(1/s_v-0.5));
			fin<<"The Maximum velocity is: "<<setprecision(6)<<u_max<<"   Re="<<Re<<"     Courant Number="<<u_max*dt/dx<<endl;
			fin<<"The max relative error of velocity is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			fin<<"Elapsed time is "<< the<<"h"<<tme<<"m"<<tse<<"s"<<endl;
			fin<<"The expected completion time is "<<th<<"h"<<tm<<"m"<<ts<<"s"<<endl;
			fin<<endl;  
			fin.close();
			}

			error=Error(u,u0,&u_max,&u_ave);if (u_max>=10.0)	U_max_ref+=1;

			     
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
		//============================================================================================
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


		
			ofstream finf3(FileName3,ios::app);
			finf3<<gx<<endl;
			finf3.close();
		

			cout<<"The"<<n<<"th computation result:"<<endl;
			cout<<"The Density of point(NX/2,NY/2,NZ/2) is: "<<setprecision(6)
				<<rho[Solid[(nx_l/2)][NY/2][NZ/2]]<<endl;
		//=============================================================================================

		//==============================================================================================

		//==============================================================================================
			
			Re=u_ave*(NY+1)/(1.0/3.0*(1/s_v-0.5));
			cout<<"The Maximum velocity is: "<<setprecision(6)<<u_max<<"   Re="<<Re<<"     Courant Number="<<u_max*dt/dx<<endl;
			
		//===============================================================================================
			cout<<"The max relative error of uv is: "
				<<setiosflags(ios::scientific)<<error<<endl;
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
				
			if ((fre_backup>=0) and (n%fre_backup==0) and (n>0))
			        Backup(n,rho,u,f);
			        
			if((error!=error) and (n>100)) {cout<<"PROGRAM STOP"<<endl;break;};
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
	delete [] forcex;
	delete [] forcey;
	delete [] forcez;

	delete [] PerC;
	delete [] PorC;

	delete [] SupInv;
	delete [] Sl;
	delete [] Sr;


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

void Parallelize_Geometry()
{
        int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	MPI_Status status;
	MPI_Request request;

int nx_g[mpi_size];
int disp[mpi_size];
int pore;
float pore2,pore3;
int loc_por[NX+1];
int sum=0;
double ave_nx;
int nx_pre,nx_aft,n_i,sum_nx;

//====================
int* Solid_rank0;
float* Per_rank0;
float* Por_rank0;
int* recv_solid;
float* recv_per;
float* recv_por;

int bufsize[mpi_size];
int bufloc[mpi_size];
//=====================


	for (int i=0;i<=NX;i++)
	        loc_por[i]=0;
	
	
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
		        }
	}
	fin.close();
}

        MPI_Bcast(loc_por,NX+1,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&sum,1,MPI_INT,0,MPI_COMM_WORLD);
        
   
         
	nx_pre=0;nx_aft=0;sum_nx=0;
	ave_nx=(double)sum/(mpi_size);
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
	for (int i=0;i<nx_l;i++)
	{
		Solid[i] = new int*[NY+1];
			for (int j=0;j<=NY;j++)
			Solid[i][j]= new int[NZ+1];
	}
		
	  Per_rank0 = new float[(NX+1)*(NY+1)*(NZ+1)];
	  Por_rank0 = new float[(NX+1)*(NY+1)*(NZ+1)];
	  
		recv_solid = new int[nx_l*(NY+1)*(NZ+1)];
		recv_per = new float[nx_l*(NY+1)*(NZ+1)];
		recv_por = new float[nx_l*(NY+1)*(NZ+1)];
//========================================
		

MPI_Barrier(MPI_COMM_WORLD);

        
        MPI_Scatterv(Solid_rank0,bufsize,bufloc,MPI_INT,recv_solid,nx_l*(NY+1)*(NZ+1),MPI_INT,0,MPI_COMM_WORLD);

	cout<<"GEOMETRY INPUT FILE PARTITIONING FOR PARALLEL READING DONE   Processor No."<<rank<<endl;
	//cout<<endl;




    
 
 
 if (rank==0)
{       
        FILE *ftest;
	ftest = fopen(perfile, "r");
	ifstream fin;
	if(ftest == NULL)
	{
		cout << "\n The Concentration file (" << perfile <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest);

	Per_rank0 = new float[(NX+1)*(NY+1)*(NZ+1)];
	
	fin.open(perfile);
	for(int k=0 ; k<=NZ ; k++)
	for(int j=0 ; j<=NY ; j++)
	for(int i=0 ; i<=NX ; i++)
	
	{
	
			fin >> pore2 >> pore3;
			Per_rank0[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=pore2;
			Por_rank0[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=pore3;
			
	
	
	}
	fin.close();
	
}


	
MPI_Barrier(MPI_COMM_WORLD);	
	
        MPI_Scatterv(Per_rank0,bufsize,bufloc,MPI_FLOAT,recv_per,nx_l*(NY+1)*(NZ+1),MPI_FLOAT,0,MPI_COMM_WORLD);
	MPI_Scatterv(Por_rank0,bufsize,bufloc,MPI_FLOAT,recv_por,nx_l*(NY+1)*(NZ+1),MPI_FLOAT,0,MPI_COMM_WORLD);
	        
	 cout<<"PERMEABILITY AND POROSITY FILE PARTITIONING FOR PARALLEL READING DONE   Processor No."<<rank<<endl;
//	cout<<endl;
   
MPI_Barrier(MPI_COMM_WORLD);		        
	        
	        
	if (rank==0)
{	

	       delete [] Solid_rank0;	
	       delete [] Per_rank0;  
	       delete [] Por_rank0;
	
}


        Per_Int = new float[nx_l*(NY+1)*(NZ+1)];
        Por_Int = new float[nx_l*(NY+1)*(NZ+1)];

	for (int i=0;i<nx_l;i++)
	                for (int j=0;j<=NY;j++)
	                for (int k=0;k<=NZ;k++)
	                {
	                Solid[i][j][k]=recv_solid[i*(NY+1)*(NZ+1)+j*(NZ+1)+k];
	                Per_Int[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=recv_per[i*(NY+1)*(NZ+1)+j*(NZ+1)+k];
	                Por_Int[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=recv_por[i*(NY+1)*(NZ+1)+j*(NZ+1)+k];
	                
	                }
	           
	 delete [] recv_per;
	 delete [] recv_solid;
	 delete [] recv_por;
       
}



void init_Sparse_read_rock_parallel(int* Sl,int* Sr)
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



double feq(int k,double rho,double u[3],double epsi)
{

	double eu,uv,feq;
        double c2,c4;
	c2=lat_c*lat_c;c4=c2*c2;
	eu=(elat[k][0]*u[0]+elat[k][1]*u[1]+elat[k][2]*u[2]);
	uv=(u[0]*u[0]+u[1]*u[1]+u[2]*u[2]);
	feq=w[k]*rho*(1.0+3.0*eu/c2+4.5*eu*eu/(c4*epsi)-1.5*uv/(c2*epsi));
	if (n==0)
	//cout<<epsi<<endl;
	return feq;

}



void init(double* rho, double** u, double** f,double* PerC, double* PorC,  double* forcex, double* forcey, double* forcez, int* SupInv)
{	

      
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	int* nx_g = new int[mpi_size];
	int* disp = new int[mpi_size];
	int ip,jp,kp,ijk;
	
	MPI_Gather(&nx_l,1,MPI_INT,nx_g,1,MPI_INT,0,MPI_COMM_WORLD);
	
	if (rank==0)
	{
		disp[0]=0;
		

		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
	}

	MPI_Bcast(disp,mpi_size,MPI_INT,0,MPI_COMM_WORLD);

	
			
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
	double uu;

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


			
	for (int i=1;i<=Count;i++)	
				
		{ 	
				
			//cout<<i<<"    @@@@@@@@@@@@@@    "<<endl;

			ip=(int)(SupInv[i]/((NY+1)*(NZ+1)));
			jp=(int)((SupInv[i]%((NY+1)*(NZ+1)))/(NZ+1));
			kp=(int)(SupInv[i]%(NZ+1));

			//ip=(ip-ip%Zoom)/Zoom;
			//jp=(jp-jp%Zoom)/Zoom;
			//kp=(kp-kp%Zoom)/Zoom;
 			//cout<<ip<<"   "<<jp<<"  "<<kp<<"    "<<disp[rank]<<endl;


			u[i][0]=inivx;
			u[i][1]=inivy;
			u[i][2]=inivz;
			u_tmp[0]=u[i][0];
			u_tmp[1]=u[i][1];
			u_tmp[2]=u[i][2];

			rho[i]=1.0;
			
			//ijk=ip*((NY+1)/Zoom)*((NZ+1)/Zoom)+jp*((NZ+1)/Zoom)+kp;
			//cout<<ijk<<"      @@@@@@@@      "<<Per_Int[10700]<<endl;


			PerC[i]=Per_Int[ip*(NY+1)*(NZ+1)+jp*(NZ+1)+kp];
			PorC[i]=Por_Int[ip*(NY+1)*(NZ+1)+jp*(NZ+1)+kp];
			
			//if (PorC[i]==0)
			 //       cout<<ip<<"  "<<jp<<"  "<<kp<<"  "<<rank<<endl;
			
			//***********************************************************************

			
			uu=sqrt(u[i][0]*u[i][0]+u[i][1]*u[i][1]+u[i][2]*u[i][2]);
			forcex[i]=-PorC[i]*niu/PerC[i]*u[i][0]-PorC[i]*F_epsilon*u[i][0]*uu/(sqrt(PerC[i]))+PorC[i]*gx;
			forcey[i]=-PorC[i]*niu/PerC[i]*u[i][1]-PorC[i]*F_epsilon*u[i][1]*uu/(sqrt(PerC[i]))+PorC[i]*gy;
			forcez[i]=-PorC[i]*niu/PerC[i]*u[i][2]-PorC[i]*F_epsilon*u[i][2]*uu/(sqrt(PerC[i]))+PorC[i]*gz;







			//***********************************************************************


			//INITIALIZATION OF m and f

			for (int lm=0;lm<19;lm++)
			
					f[i][lm]=feq(lm,rho[i],u_tmp,PorC[i]);
				

		
		

	}

	
	delete [] nx_g;
	delete [] disp;

	 	
}




void collision(double* rho,double** u,double** f,double** F, double* PerC,double* PorC,double* forcex, double* forcey, double* forcez, int* SupInv,int*** Solid, int* Sl, int* Sr)
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
double usqr,vsqr,uu;
double F_hat[19],GuoF[19],f_eq[19],u_tmp[3];
double m_l[19];
int i,j,m,ip,jp,kp;

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
			

			uu=sqrt(u[ci][0]*u[ci][0]+u[ci][1]*u[ci][1]+u[ci][2]*u[ci][2]);
			forcex[ci]=-PorC[ci]*in_vis/PerC[ci]*u[ci][0]-PorC[ci]*F_epsilon*u[ci][0]*uu/(sqrt(PerC[ci]))+PorC[ci]*gx;
			forcey[ci]=-PorC[ci]*in_vis/PerC[ci]*u[ci][1]-PorC[ci]*F_epsilon*u[ci][1]*uu/(sqrt(PerC[ci]))+PorC[ci]*gy;
			forcez[ci]=-PorC[ci]*in_vis/PerC[ci]*u[ci][2]-PorC[ci]*F_epsilon*u[ci][2]*uu/(sqrt(PerC[ci]))+PorC[ci]*gz;

			

			//=================FORCE TERM_GUO=========================================
			for (int k=0;k<19;k++)
			{	
			lm0=((elat[k][0]-u[ci][0]/PorC[ci])*forcex[ci]+(elat[k][1]-u[ci][1]/PorC[ci])*forcey[ci]+(elat[k][2]-u[ci][2]/PorC[ci])*forcez[ci])/c_s2;
			lm1=(elat[k][0]*u[ci][0]+elat[k][1]*u[ci][1]+elat[k][2]*u[ci][2])*(elat[k][0]*forcex[ci]+elat[k][1]*forcey[ci]+elat[k][2]*forcez[ci])/(c_s2*c_s2);
			lm1/=PorC[ci];
			GuoF[k]=w[k]*(lm0+lm1);
			//GuoF[k]=0.0;
			}


			
			//=====================equilibrium of moment=================================
			u_tmp[0]=u[ci][0];
			u_tmp[1]=u[ci][1];
			u_tmp[2]=u[ci][2];

			for(int k=0;k<19;k++)
				{
				f_eq[k]=feq(k,rho[ci],u_tmp,PorC[ci]);
				}
			
			

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
m_l[0]=+1.000*f[ci][0]+1.00000000000000000*f[ci][1]+1.00000000000000000*f[ci][2]+1.00000000000000000*f[ci][3]+1.00000000000000000*f[ci][4]+1.00000000000000000*f[ci][5]+1.00000000000000000*f[ci][6]+1.00000000000000000*f[ci][7]+1.00000000000000000*f[ci][8]+1.00000000000000000*f[ci][9]+1.00000000000000000*f[ci][10]+1.00000000000000000*f[ci][11]+1.00000000000000000*f[ci][12]+1.00000000000000000*f[ci][13]+1.00000000000000000*f[ci][14]+1.00000000000000000*f[ci][15]+1.00000000000000000*f[ci][16]+1.00000000000000000*f[ci][17]+1.00000000000000000*f[ci][18];

F_hat[0]=+1.000*GuoF[0]+1.00000000000000000*GuoF[1]+1.00000000000000000*GuoF[2]+1.00000000000000000*GuoF[3]+1.00000000000000000*GuoF[4]+1.00000000000000000*GuoF[5]+1.00000000000000000*GuoF[6]+1.00000000000000000*GuoF[7]+1.00000000000000000*GuoF[8]+1.00000000000000000*GuoF[9]+1.00000000000000000*GuoF[10]+1.00000000000000000*GuoF[11]+1.00000000000000000*GuoF[12]+1.00000000000000000*GuoF[13]+1.00000000000000000*GuoF[14]+1.00000000000000000*GuoF[15]+1.00000000000000000*GuoF[16]+1.00000000000000000*GuoF[17]+1.00000000000000000*GuoF[18];

meq[0]=+1.000*f_eq[0]+1.00000000000000000*f_eq[1]+1.00000000000000000*f_eq[2]+1.00000000000000000*f_eq[3]+1.00000000000000000*f_eq[4]+1.00000000000000000*f_eq[5]+1.00000000000000000*f_eq[6]+1.00000000000000000*f_eq[7]+1.00000000000000000*f_eq[8]+1.00000000000000000*f_eq[9]+1.00000000000000000*f_eq[10]+1.00000000000000000*f_eq[11]+1.00000000000000000*f_eq[12]+1.00000000000000000*f_eq[13]+1.00000000000000000*f_eq[14]+1.00000000000000000*f_eq[15]+1.00000000000000000*f_eq[16]+1.00000000000000000*f_eq[17]+1.00000000000000000*f_eq[18];

F_hat[0]*=(1-0.5*S[0]);
m_l[0]=m_l[0]-S[0]*(m_l[0]-meq[0])+dt*F_hat[0];
//=======================================

m_l[1]=-30.00000000000000000*f[ci][0]-11.00000000000000000*f[ci][1]-11.00000000000000000*f[ci][2]-11.00000000000000000*f[ci][3]-11.00000000000000000*f[ci][4]-11.00000000000000000*f[ci][5]-11.00000000000000000*f[ci][6]+8.00000000000000000*f[ci][7]+8.00000000000000000*f[ci][8]+8.00000000000000000*f[ci][9]+8.00000000000000000*f[ci][10]+8.00000000000000000*f[ci][11]+8.00000000000000000*f[ci][12]+8.00000000000000000*f[ci][13]+8.00000000000000000*f[ci][14]+8.00000000000000000*f[ci][15]+8.00000000000000000*f[ci][16]+8.00000000000000000*f[ci][17]+8.00000000000000000*f[ci][18];

F_hat[1]=-30.00000000000000000*GuoF[0]-11.00000000000000000*GuoF[1]-11.00000000000000000*GuoF[2]-11.00000000000000000*GuoF[3]-11.00000000000000000*GuoF[4]-11.00000000000000000*GuoF[5]-11.00000000000000000*GuoF[6]+8.00000000000000000*GuoF[7]+8.00000000000000000*GuoF[8]+8.00000000000000000*GuoF[9]+8.00000000000000000*GuoF[10]+8.00000000000000000*GuoF[11]+8.00000000000000000*GuoF[12]+8.00000000000000000*GuoF[13]+8.00000000000000000*GuoF[14]+8.00000000000000000*GuoF[15]+8.00000000000000000*GuoF[16]+8.00000000000000000*GuoF[17]+8.00000000000000000*GuoF[18];

meq[1]=-30.00000000000000000*f_eq[0]-11.00000000000000000*f_eq[1]-11.00000000000000000*f_eq[2]-11.00000000000000000*f_eq[3]-11.00000000000000000*f_eq[4]-11.00000000000000000*f_eq[5]-11.00000000000000000*f_eq[6]+8.00000000000000000*f_eq[7]+8.00000000000000000*f_eq[8]+8.00000000000000000*f_eq[9]+8.00000000000000000*f_eq[10]+8.00000000000000000*f_eq[11]+8.00000000000000000*f_eq[12]+8.00000000000000000*f_eq[13]+8.00000000000000000*f_eq[14]+8.00000000000000000*f_eq[15]+8.00000000000000000*f_eq[16]+8.00000000000000000*f_eq[17]+8.00000000000000000*f_eq[18];

F_hat[1]*=(1-0.5*S[1]);
m_l[1]=m_l[1]-S[1]*(m_l[1]-meq[1])+dt*F_hat[1];
//=======================================

m_l[2]=+12.00000000000000000*f[ci][0]-4.00000000000000000*f[ci][1]-4.00000000000000000*f[ci][2]-4.00000000000000000*f[ci][3]-4.00000000000000000*f[ci][4]-4.00000000000000000*f[ci][5]-4.00000000000000000*f[ci][6]+1.00000000000000000*f[ci][7]+1.00000000000000000*f[ci][8]+1.00000000000000000*f[ci][9]+1.00000000000000000*f[ci][10]+1.00000000000000000*f[ci][11]+1.00000000000000000*f[ci][12]+1.00000000000000000*f[ci][13]+1.00000000000000000*f[ci][14]+1.00000000000000000*f[ci][15]+1.00000000000000000*f[ci][16]+1.00000000000000000*f[ci][17]+1.00000000000000000*f[ci][18];

F_hat[2]=+12.00000000000000000*GuoF[0]-4.00000000000000000*GuoF[1]-4.00000000000000000*GuoF[2]-4.00000000000000000*GuoF[3]-4.00000000000000000*GuoF[4]-4.00000000000000000*GuoF[5]-4.00000000000000000*GuoF[6]+1.00000000000000000*GuoF[7]+1.00000000000000000*GuoF[8]+1.00000000000000000*GuoF[9]+1.00000000000000000*GuoF[10]+1.00000000000000000*GuoF[11]+1.00000000000000000*GuoF[12]+1.00000000000000000*GuoF[13]+1.00000000000000000*GuoF[14]+1.00000000000000000*GuoF[15]+1.00000000000000000*GuoF[16]+1.00000000000000000*GuoF[17]+1.00000000000000000*GuoF[18];

meq[2]=+12.00000000000000000*f_eq[0]-4.00000000000000000*f_eq[1]-4.00000000000000000*f_eq[2]-4.00000000000000000*f_eq[3]-4.00000000000000000*f_eq[4]-4.00000000000000000*f_eq[5]-4.00000000000000000*f_eq[6]+1.00000000000000000*f_eq[7]+1.00000000000000000*f_eq[8]+1.00000000000000000*f_eq[9]+1.00000000000000000*f_eq[10]+1.00000000000000000*f_eq[11]+1.00000000000000000*f_eq[12]+1.00000000000000000*f_eq[13]+1.00000000000000000*f_eq[14]+1.00000000000000000*f_eq[15]+1.00000000000000000*f_eq[16]+1.00000000000000000*f_eq[17]+1.00000000000000000*f_eq[18];

F_hat[2]*=(1-0.5*S[2]);
m_l[2]=m_l[2]-S[2]*(m_l[2]-meq[2])+dt*F_hat[2];
//=======================================

m_l[3]=+0.00000000000000000*f[ci][0]+1.00000000000000000*f[ci][1]-1.00000000000000000*f[ci][2]+0.00000000000000000*f[ci][3]+0.00000000000000000*f[ci][4]+0.00000000000000000*f[ci][5]+0.00000000000000000*f[ci][6]+1.00000000000000000*f[ci][7]-1.00000000000000000*f[ci][8]+1.00000000000000000*f[ci][9]-1.00000000000000000*f[ci][10]+0.00000000000000000*f[ci][11]+0.00000000000000000*f[ci][12]+0.00000000000000000*f[ci][13]+0.00000000000000000*f[ci][14]+1.00000000000000000*f[ci][15]-1.00000000000000000*f[ci][16]+1.00000000000000000*f[ci][17]-1.00000000000000000*f[ci][18];

F_hat[3]=+0.00000000000000000*GuoF[0]+1.00000000000000000*GuoF[1]-1.00000000000000000*GuoF[2]+0.00000000000000000*GuoF[3]+0.00000000000000000*GuoF[4]+0.00000000000000000*GuoF[5]+0.00000000000000000*GuoF[6]+1.00000000000000000*GuoF[7]-1.00000000000000000*GuoF[8]+1.00000000000000000*GuoF[9]-1.00000000000000000*GuoF[10]+0.00000000000000000*GuoF[11]+0.00000000000000000*GuoF[12]+0.00000000000000000*GuoF[13]+0.00000000000000000*GuoF[14]+1.00000000000000000*GuoF[15]-1.00000000000000000*GuoF[16]+1.00000000000000000*GuoF[17]-1.00000000000000000*GuoF[18];

meq[3]=+0.00000000000000000*f_eq[0]+1.00000000000000000*f_eq[1]-1.00000000000000000*f_eq[2]+0.00000000000000000*f_eq[3]+0.00000000000000000*f_eq[4]+0.00000000000000000*f_eq[5]+0.00000000000000000*f_eq[6]+1.00000000000000000*f_eq[7]-1.00000000000000000*f_eq[8]+1.00000000000000000*f_eq[9]-1.00000000000000000*f_eq[10]+0.00000000000000000*f_eq[11]+0.00000000000000000*f_eq[12]+0.00000000000000000*f_eq[13]+0.00000000000000000*f_eq[14]+1.00000000000000000*f_eq[15]-1.00000000000000000*f_eq[16]+1.00000000000000000*f_eq[17]-1.00000000000000000*f_eq[18];

F_hat[3]*=(1-0.5*S[3]);
m_l[3]=m_l[3]-S[3]*(m_l[3]-meq[3])+dt*F_hat[3];
//=======================================

m_l[4]=+0.00000000000000000*f[ci][0]-4.00000000000000000*f[ci][1]+4.00000000000000000*f[ci][2]+0.00000000000000000*f[ci][3]+0.00000000000000000*f[ci][4]+0.00000000000000000*f[ci][5]+0.00000000000000000*f[ci][6]+1.00000000000000000*f[ci][7]-1.00000000000000000*f[ci][8]+1.00000000000000000*f[ci][9]-1.00000000000000000*f[ci][10]+0.00000000000000000*f[ci][11]+0.00000000000000000*f[ci][12]+0.00000000000000000*f[ci][13]+0.00000000000000000*f[ci][14]+1.00000000000000000*f[ci][15]-1.00000000000000000*f[ci][16]+1.00000000000000000*f[ci][17]-1.00000000000000000*f[ci][18];

F_hat[4]=+0.00000000000000000*GuoF[0]-4.00000000000000000*GuoF[1]+4.00000000000000000*GuoF[2]+0.00000000000000000*GuoF[3]+0.00000000000000000*GuoF[4]+0.00000000000000000*GuoF[5]+0.00000000000000000*GuoF[6]+1.00000000000000000*GuoF[7]-1.00000000000000000*GuoF[8]+1.00000000000000000*GuoF[9]-1.00000000000000000*GuoF[10]+0.00000000000000000*GuoF[11]+0.00000000000000000*GuoF[12]+0.00000000000000000*GuoF[13]+0.00000000000000000*GuoF[14]+1.00000000000000000*GuoF[15]-1.00000000000000000*GuoF[16]+1.00000000000000000*GuoF[17]-1.00000000000000000*GuoF[18];

meq[4]=+0.00000000000000000*f_eq[0]-4.00000000000000000*f_eq[1]+4.00000000000000000*f_eq[2]+0.00000000000000000*f_eq[3]+0.00000000000000000*f_eq[4]+0.00000000000000000*f_eq[5]+0.00000000000000000*f_eq[6]+1.00000000000000000*f_eq[7]-1.00000000000000000*f_eq[8]+1.00000000000000000*f_eq[9]-1.00000000000000000*f_eq[10]+0.00000000000000000*f_eq[11]+0.00000000000000000*f_eq[12]+0.00000000000000000*f_eq[13]+0.00000000000000000*f_eq[14]+1.00000000000000000*f_eq[15]-1.00000000000000000*f_eq[16]+1.00000000000000000*f_eq[17]-1.00000000000000000*f_eq[18];

F_hat[4]*=(1-0.5*S[4]);
m_l[4]=m_l[4]-S[4]*(m_l[4]-meq[4])+dt*F_hat[4];
//=======================================

m_l[5]=+0.00000000000000000*f[ci][0]+0.00000000000000000*f[ci][1]+0.00000000000000000*f[ci][2]+1.00000000000000000*f[ci][3]-1.00000000000000000*f[ci][4]+0.00000000000000000*f[ci][5]+0.00000000000000000*f[ci][6]+1.00000000000000000*f[ci][7]+1.00000000000000000*f[ci][8]-1.00000000000000000*f[ci][9]-1.00000000000000000*f[ci][10]+1.00000000000000000*f[ci][11]-1.00000000000000000*f[ci][12]+1.00000000000000000*f[ci][13]-1.00000000000000000*f[ci][14]+0.00000000000000000*f[ci][15]+0.00000000000000000*f[ci][16]+0.00000000000000000*f[ci][17]+0.00000000000000000*f[ci][18];

F_hat[5]=+0.00000000000000000*GuoF[0]+0.00000000000000000*GuoF[1]+0.00000000000000000*GuoF[2]+1.00000000000000000*GuoF[3]-1.00000000000000000*GuoF[4]+0.00000000000000000*GuoF[5]+0.00000000000000000*GuoF[6]+1.00000000000000000*GuoF[7]+1.00000000000000000*GuoF[8]-1.00000000000000000*GuoF[9]-1.00000000000000000*GuoF[10]+1.00000000000000000*GuoF[11]-1.00000000000000000*GuoF[12]+1.00000000000000000*GuoF[13]-1.00000000000000000*GuoF[14]+0.00000000000000000*GuoF[15]+0.00000000000000000*GuoF[16]+0.00000000000000000*GuoF[17]+0.00000000000000000*GuoF[18];

meq[5]=+0.00000000000000000*f_eq[0]+0.00000000000000000*f_eq[1]+0.00000000000000000*f_eq[2]+1.00000000000000000*f_eq[3]-1.00000000000000000*f_eq[4]+0.00000000000000000*f_eq[5]+0.00000000000000000*f_eq[6]+1.00000000000000000*f_eq[7]+1.00000000000000000*f_eq[8]-1.00000000000000000*f_eq[9]-1.00000000000000000*f_eq[10]+1.00000000000000000*f_eq[11]-1.00000000000000000*f_eq[12]+1.00000000000000000*f_eq[13]-1.00000000000000000*f_eq[14]+0.00000000000000000*f_eq[15]+0.00000000000000000*f_eq[16]+0.00000000000000000*f_eq[17]+0.00000000000000000*f_eq[18];

F_hat[5]*=(1-0.5*S[5]);
m_l[5]=m_l[5]-S[5]*(m_l[5]-meq[5])+dt*F_hat[5];
//=======================================

m_l[6]=+0.00000000000000000*f[ci][0]+0.00000000000000000*f[ci][1]+0.00000000000000000*f[ci][2]-4.00000000000000000*f[ci][3]+4.00000000000000000*f[ci][4]+0.00000000000000000*f[ci][5]+0.00000000000000000*f[ci][6]+1.00000000000000000*f[ci][7]+1.00000000000000000*f[ci][8]-1.00000000000000000*f[ci][9]-1.00000000000000000*f[ci][10]+1.00000000000000000*f[ci][11]-1.00000000000000000*f[ci][12]+1.00000000000000000*f[ci][13]-1.00000000000000000*f[ci][14]+0.00000000000000000*f[ci][15]+0.00000000000000000*f[ci][16]+0.00000000000000000*f[ci][17]+0.00000000000000000*f[ci][18];

F_hat[6]=+0.00000000000000000*GuoF[0]+0.00000000000000000*GuoF[1]+0.00000000000000000*GuoF[2]-4.00000000000000000*GuoF[3]+4.00000000000000000*GuoF[4]+0.00000000000000000*GuoF[5]+0.00000000000000000*GuoF[6]+1.00000000000000000*GuoF[7]+1.00000000000000000*GuoF[8]-1.00000000000000000*GuoF[9]-1.00000000000000000*GuoF[10]+1.00000000000000000*GuoF[11]-1.00000000000000000*GuoF[12]+1.00000000000000000*GuoF[13]-1.00000000000000000*GuoF[14]+0.00000000000000000*GuoF[15]+0.00000000000000000*GuoF[16]+0.00000000000000000*GuoF[17]+0.00000000000000000*GuoF[18];

meq[6]=+0.00000000000000000*f_eq[0]+0.00000000000000000*f_eq[1]+0.00000000000000000*f_eq[2]-4.00000000000000000*f_eq[3]+4.00000000000000000*f_eq[4]+0.00000000000000000*f_eq[5]+0.00000000000000000*f_eq[6]+1.00000000000000000*f_eq[7]+1.00000000000000000*f_eq[8]-1.00000000000000000*f_eq[9]-1.00000000000000000*f_eq[10]+1.00000000000000000*f_eq[11]-1.00000000000000000*f_eq[12]+1.00000000000000000*f_eq[13]-1.00000000000000000*f_eq[14]+0.00000000000000000*f_eq[15]+0.00000000000000000*f_eq[16]+0.00000000000000000*f_eq[17]+0.00000000000000000*f_eq[18];

F_hat[6]*=(1-0.5*S[6]);
m_l[6]=m_l[6]-S[6]*(m_l[6]-meq[6])+dt*F_hat[6];
//=======================================

m_l[7]=+0.00000000000000000*f[ci][0]+0.00000000000000000*f[ci][1]+0.00000000000000000*f[ci][2]+0.00000000000000000*f[ci][3]+0.00000000000000000*f[ci][4]+1.00000000000000000*f[ci][5]-1.00000000000000000*f[ci][6]+0.00000000000000000*f[ci][7]+0.00000000000000000*f[ci][8]+0.00000000000000000*f[ci][9]+0.00000000000000000*f[ci][10]+1.00000000000000000*f[ci][11]+1.00000000000000000*f[ci][12]-1.00000000000000000*f[ci][13]-1.00000000000000000*f[ci][14]+1.00000000000000000*f[ci][15]+1.00000000000000000*f[ci][16]-1.00000000000000000*f[ci][17]-1.00000000000000000*f[ci][18];

F_hat[7]=+0.00000000000000000*GuoF[0]+0.00000000000000000*GuoF[1]+0.00000000000000000*GuoF[2]+0.00000000000000000*GuoF[3]+0.00000000000000000*GuoF[4]+1.00000000000000000*GuoF[5]-1.00000000000000000*GuoF[6]+0.00000000000000000*GuoF[7]+0.00000000000000000*GuoF[8]+0.00000000000000000*GuoF[9]+0.00000000000000000*GuoF[10]+1.00000000000000000*GuoF[11]+1.00000000000000000*GuoF[12]-1.00000000000000000*GuoF[13]-1.00000000000000000*GuoF[14]+1.00000000000000000*GuoF[15]+1.00000000000000000*GuoF[16]-1.00000000000000000*GuoF[17]-1.00000000000000000*GuoF[18];

meq[7]=+0.00000000000000000*f_eq[0]+0.00000000000000000*f_eq[1]+0.00000000000000000*f_eq[2]+0.00000000000000000*f_eq[3]+0.00000000000000000*f_eq[4]+1.00000000000000000*f_eq[5]-1.00000000000000000*f_eq[6]+0.00000000000000000*f_eq[7]+0.00000000000000000*f_eq[8]+0.00000000000000000*f_eq[9]+0.00000000000000000*f_eq[10]+1.00000000000000000*f_eq[11]+1.00000000000000000*f_eq[12]-1.00000000000000000*f_eq[13]-1.00000000000000000*f_eq[14]+1.00000000000000000*f_eq[15]+1.00000000000000000*f_eq[16]-1.00000000000000000*f_eq[17]-1.00000000000000000*f_eq[18];

F_hat[7]*=(1-0.5*S[7]);
m_l[7]=m_l[7]-S[7]*(m_l[7]-meq[7])+dt*F_hat[7];
//=======================================

m_l[8]=+0.00000000000000000*f[ci][0]+0.00000000000000000*f[ci][1]+0.00000000000000000*f[ci][2]+0.00000000000000000*f[ci][3]+0.00000000000000000*f[ci][4]-4.00000000000000000*f[ci][5]+4.00000000000000000*f[ci][6]+0.00000000000000000*f[ci][7]+0.00000000000000000*f[ci][8]+0.00000000000000000*f[ci][9]+0.00000000000000000*f[ci][10]+1.00000000000000000*f[ci][11]+1.00000000000000000*f[ci][12]-1.00000000000000000*f[ci][13]-1.00000000000000000*f[ci][14]+1.00000000000000000*f[ci][15]+1.00000000000000000*f[ci][16]-1.00000000000000000*f[ci][17]-1.00000000000000000*f[ci][18];

F_hat[8]=+0.00000000000000000*GuoF[0]+0.00000000000000000*GuoF[1]+0.00000000000000000*GuoF[2]+0.00000000000000000*GuoF[3]+0.00000000000000000*GuoF[4]-4.00000000000000000*GuoF[5]+4.00000000000000000*GuoF[6]+0.00000000000000000*GuoF[7]+0.00000000000000000*GuoF[8]+0.00000000000000000*GuoF[9]+0.00000000000000000*GuoF[10]+1.00000000000000000*GuoF[11]+1.00000000000000000*GuoF[12]-1.00000000000000000*GuoF[13]-1.00000000000000000*GuoF[14]+1.00000000000000000*GuoF[15]+1.00000000000000000*GuoF[16]-1.00000000000000000*GuoF[17]-1.00000000000000000*GuoF[18];

meq[8]=+0.00000000000000000*f_eq[0]+0.00000000000000000*f_eq[1]+0.00000000000000000*f_eq[2]+0.00000000000000000*f_eq[3]+0.00000000000000000*f_eq[4]-4.00000000000000000*f_eq[5]+4.00000000000000000*f_eq[6]+0.00000000000000000*f_eq[7]+0.00000000000000000*f_eq[8]+0.00000000000000000*f_eq[9]+0.00000000000000000*f_eq[10]+1.00000000000000000*f_eq[11]+1.00000000000000000*f_eq[12]-1.00000000000000000*f_eq[13]-1.00000000000000000*f_eq[14]+1.00000000000000000*f_eq[15]+1.00000000000000000*f_eq[16]-1.00000000000000000*f_eq[17]-1.00000000000000000*f_eq[18];

F_hat[8]*=(1-0.5*S[8]);
m_l[8]=m_l[8]-S[8]*(m_l[8]-meq[8])+dt*F_hat[8];
//=======================================

m_l[9]=+0.00000000000000000*f[ci][0]+2.00000000000000000*f[ci][1]+2.00000000000000000*f[ci][2]-1.00000000000000000*f[ci][3]-1.00000000000000000*f[ci][4]-1.00000000000000000*f[ci][5]-1.00000000000000000*f[ci][6]+1.00000000000000000*f[ci][7]+1.00000000000000000*f[ci][8]+1.00000000000000000*f[ci][9]+1.00000000000000000*f[ci][10]-2.00000000000000000*f[ci][11]-2.00000000000000000*f[ci][12]-2.00000000000000000*f[ci][13]-2.00000000000000000*f[ci][14]+1.00000000000000000*f[ci][15]+1.00000000000000000*f[ci][16]+1.00000000000000000*f[ci][17]+1.00000000000000000*f[ci][18];

F_hat[9]=+0.00000000000000000*GuoF[0]+2.00000000000000000*GuoF[1]+2.00000000000000000*GuoF[2]-1.00000000000000000*GuoF[3]-1.00000000000000000*GuoF[4]-1.00000000000000000*GuoF[5]-1.00000000000000000*GuoF[6]+1.00000000000000000*GuoF[7]+1.00000000000000000*GuoF[8]+1.00000000000000000*GuoF[9]+1.00000000000000000*GuoF[10]-2.00000000000000000*GuoF[11]-2.00000000000000000*GuoF[12]-2.00000000000000000*GuoF[13]-2.00000000000000000*GuoF[14]+1.00000000000000000*GuoF[15]+1.00000000000000000*GuoF[16]+1.00000000000000000*GuoF[17]+1.00000000000000000*GuoF[18];

meq[9]=+0.00000000000000000*f_eq[0]+2.00000000000000000*f_eq[1]+2.00000000000000000*f_eq[2]-1.00000000000000000*f_eq[3]-1.00000000000000000*f_eq[4]-1.00000000000000000*f_eq[5]-1.00000000000000000*f_eq[6]+1.00000000000000000*f_eq[7]+1.00000000000000000*f_eq[8]+1.00000000000000000*f_eq[9]+1.00000000000000000*f_eq[10]-2.00000000000000000*f_eq[11]-2.00000000000000000*f_eq[12]-2.00000000000000000*f_eq[13]-2.00000000000000000*f_eq[14]+1.00000000000000000*f_eq[15]+1.00000000000000000*f_eq[16]+1.00000000000000000*f_eq[17]+1.00000000000000000*f_eq[18];

F_hat[9]*=(1-0.5*S[9]);
m_l[9]=m_l[9]-S[9]*(m_l[9]-meq[9])+dt*F_hat[9];
//=======================================

m_l[10]=+0.00000000000000000*f[ci][0]-4.00000000000000000*f[ci][1]-4.00000000000000000*f[ci][2]+2.00000000000000000*f[ci][3]+2.00000000000000000*f[ci][4]+2.00000000000000000*f[ci][5]+2.00000000000000000*f[ci][6]+1.00000000000000000*f[ci][7]+1.00000000000000000*f[ci][8]+1.00000000000000000*f[ci][9]+1.00000000000000000*f[ci][10]-2.00000000000000000*f[ci][11]-2.00000000000000000*f[ci][12]-2.00000000000000000*f[ci][13]-2.00000000000000000*f[ci][14]+1.00000000000000000*f[ci][15]+1.00000000000000000*f[ci][16]+1.00000000000000000*f[ci][17]+1.00000000000000000*f[ci][18];

F_hat[10]=+0.00000000000000000*GuoF[0]-4.00000000000000000*GuoF[1]-4.00000000000000000*GuoF[2]+2.00000000000000000*GuoF[3]+2.00000000000000000*GuoF[4]+2.00000000000000000*GuoF[5]+2.00000000000000000*GuoF[6]+1.00000000000000000*GuoF[7]+1.00000000000000000*GuoF[8]+1.00000000000000000*GuoF[9]+1.00000000000000000*GuoF[10]-2.00000000000000000*GuoF[11]-2.00000000000000000*GuoF[12]-2.00000000000000000*GuoF[13]-2.00000000000000000*GuoF[14]+1.00000000000000000*GuoF[15]+1.00000000000000000*GuoF[16]+1.00000000000000000*GuoF[17]+1.00000000000000000*GuoF[18];

meq[10]=+0.00000000000000000*f_eq[0]-4.00000000000000000*f_eq[1]-4.00000000000000000*f_eq[2]+2.00000000000000000*f_eq[3]+2.00000000000000000*f_eq[4]+2.00000000000000000*f_eq[5]+2.00000000000000000*f_eq[6]+1.00000000000000000*f_eq[7]+1.00000000000000000*f_eq[8]+1.00000000000000000*f_eq[9]+1.00000000000000000*f_eq[10]-2.00000000000000000*f_eq[11]-2.00000000000000000*f_eq[12]-2.00000000000000000*f_eq[13]-2.00000000000000000*f_eq[14]+1.00000000000000000*f_eq[15]+1.00000000000000000*f_eq[16]+1.00000000000000000*f_eq[17]+1.00000000000000000*f_eq[18];

F_hat[10]*=(1-0.5*S[10]);
m_l[10]=m_l[10]-S[10]*(m_l[10]-meq[10])+dt*F_hat[10];
//=======================================

m_l[11]=+0.00000000000000000*f[ci][0]+0.00000000000000000*f[ci][1]+0.00000000000000000*f[ci][2]+1.00000000000000000*f[ci][3]+1.00000000000000000*f[ci][4]-1.00000000000000000*f[ci][5]-1.00000000000000000*f[ci][6]+1.00000000000000000*f[ci][7]+1.00000000000000000*f[ci][8]+1.00000000000000000*f[ci][9]+1.00000000000000000*f[ci][10]+0.00000000000000000*f[ci][11]+0.00000000000000000*f[ci][12]+0.00000000000000000*f[ci][13]+0.00000000000000000*f[ci][14]-1.00000000000000000*f[ci][15]-1.00000000000000000*f[ci][16]-1.00000000000000000*f[ci][17]-1.00000000000000000*f[ci][18];

F_hat[11]=+0.00000000000000000*GuoF[0]+0.00000000000000000*GuoF[1]+0.00000000000000000*GuoF[2]+1.00000000000000000*GuoF[3]+1.00000000000000000*GuoF[4]-1.00000000000000000*GuoF[5]-1.00000000000000000*GuoF[6]+1.00000000000000000*GuoF[7]+1.00000000000000000*GuoF[8]+1.00000000000000000*GuoF[9]+1.00000000000000000*GuoF[10]+0.00000000000000000*GuoF[11]+0.00000000000000000*GuoF[12]+0.00000000000000000*GuoF[13]+0.00000000000000000*GuoF[14]-1.00000000000000000*GuoF[15]-1.00000000000000000*GuoF[16]-1.00000000000000000*GuoF[17]-1.00000000000000000*GuoF[18];

meq[11]=+0.00000000000000000*f_eq[0]+0.00000000000000000*f_eq[1]+0.00000000000000000*f_eq[2]+1.00000000000000000*f_eq[3]+1.00000000000000000*f_eq[4]-1.00000000000000000*f_eq[5]-1.00000000000000000*f_eq[6]+1.00000000000000000*f_eq[7]+1.00000000000000000*f_eq[8]+1.00000000000000000*f_eq[9]+1.00000000000000000*f_eq[10]+0.00000000000000000*f_eq[11]+0.00000000000000000*f_eq[12]+0.00000000000000000*f_eq[13]+0.00000000000000000*f_eq[14]-1.00000000000000000*f_eq[15]-1.00000000000000000*f_eq[16]-1.00000000000000000*f_eq[17]-1.00000000000000000*f_eq[18];

F_hat[11]*=(1-0.5*S[11]);
m_l[11]=m_l[11]-S[11]*(m_l[11]-meq[11])+dt*F_hat[11];
//=======================================

m_l[12]=+0.00000000000000000*f[ci][0]+0.00000000000000000*f[ci][1]+0.00000000000000000*f[ci][2]-2.00000000000000000*f[ci][3]-2.00000000000000000*f[ci][4]+2.00000000000000000*f[ci][5]+2.00000000000000000*f[ci][6]+1.00000000000000000*f[ci][7]+1.00000000000000000*f[ci][8]+1.00000000000000000*f[ci][9]+1.00000000000000000*f[ci][10]+0.00000000000000000*f[ci][11]+0.00000000000000000*f[ci][12]+0.00000000000000000*f[ci][13]+0.00000000000000000*f[ci][14]-1.00000000000000000*f[ci][15]-1.00000000000000000*f[ci][16]-1.00000000000000000*f[ci][17]-1.00000000000000000*f[ci][18];

F_hat[12]=+0.00000000000000000*GuoF[0]+0.00000000000000000*GuoF[1]+0.00000000000000000*GuoF[2]-2.00000000000000000*GuoF[3]-2.00000000000000000*GuoF[4]+2.00000000000000000*GuoF[5]+2.00000000000000000*GuoF[6]+1.00000000000000000*GuoF[7]+1.00000000000000000*GuoF[8]+1.00000000000000000*GuoF[9]+1.00000000000000000*GuoF[10]+0.00000000000000000*GuoF[11]+0.00000000000000000*GuoF[12]+0.00000000000000000*GuoF[13]+0.00000000000000000*GuoF[14]-1.00000000000000000*GuoF[15]-1.00000000000000000*GuoF[16]-1.00000000000000000*GuoF[17]-1.00000000000000000*GuoF[18];

meq[12]=+0.00000000000000000*f_eq[0]+0.00000000000000000*f_eq[1]+0.00000000000000000*f_eq[2]-2.00000000000000000*f_eq[3]-2.00000000000000000*f_eq[4]+2.00000000000000000*f_eq[5]+2.00000000000000000*f_eq[6]+1.00000000000000000*f_eq[7]+1.00000000000000000*f_eq[8]+1.00000000000000000*f_eq[9]+1.00000000000000000*f_eq[10]+0.00000000000000000*f_eq[11]+0.00000000000000000*f_eq[12]+0.00000000000000000*f_eq[13]+0.00000000000000000*f_eq[14]-1.00000000000000000*f_eq[15]-1.00000000000000000*f_eq[16]-1.00000000000000000*f_eq[17]-1.00000000000000000*f_eq[18];

F_hat[12]*=(1-0.5*S[12]);
m_l[12]=m_l[12]-S[12]*(m_l[12]-meq[12])+dt*F_hat[12];
//=======================================

m_l[13]=+0.00000000000000000*f[ci][0]+0.00000000000000000*f[ci][1]+0.00000000000000000*f[ci][2]+0.00000000000000000*f[ci][3]+0.00000000000000000*f[ci][4]+0.00000000000000000*f[ci][5]+0.00000000000000000*f[ci][6]+1.00000000000000000*f[ci][7]-1.00000000000000000*f[ci][8]-1.00000000000000000*f[ci][9]+1.00000000000000000*f[ci][10]+0.00000000000000000*f[ci][11]+0.00000000000000000*f[ci][12]+0.00000000000000000*f[ci][13]+0.00000000000000000*f[ci][14]+0.00000000000000000*f[ci][15]+0.00000000000000000*f[ci][16]+0.00000000000000000*f[ci][17]+0.00000000000000000*f[ci][18];

F_hat[13]=+0.00000000000000000*GuoF[0]+0.00000000000000000*GuoF[1]+0.00000000000000000*GuoF[2]+0.00000000000000000*GuoF[3]+0.00000000000000000*GuoF[4]+0.00000000000000000*GuoF[5]+0.00000000000000000*GuoF[6]+1.00000000000000000*GuoF[7]-1.00000000000000000*GuoF[8]-1.00000000000000000*GuoF[9]+1.00000000000000000*GuoF[10]+0.00000000000000000*GuoF[11]+0.00000000000000000*GuoF[12]+0.00000000000000000*GuoF[13]+0.00000000000000000*GuoF[14]+0.00000000000000000*GuoF[15]+0.00000000000000000*GuoF[16]+0.00000000000000000*GuoF[17]+0.00000000000000000*GuoF[18];

meq[13]=+0.00000000000000000*f_eq[0]+0.00000000000000000*f_eq[1]+0.00000000000000000*f_eq[2]+0.00000000000000000*f_eq[3]+0.00000000000000000*f_eq[4]+0.00000000000000000*f_eq[5]+0.00000000000000000*f_eq[6]+1.00000000000000000*f_eq[7]-1.00000000000000000*f_eq[8]-1.00000000000000000*f_eq[9]+1.00000000000000000*f_eq[10]+0.00000000000000000*f_eq[11]+0.00000000000000000*f_eq[12]+0.00000000000000000*f_eq[13]+0.00000000000000000*f_eq[14]+0.00000000000000000*f_eq[15]+0.00000000000000000*f_eq[16]+0.00000000000000000*f_eq[17]+0.00000000000000000*f_eq[18];

F_hat[13]*=(1-0.5*S[13]);
m_l[13]=m_l[13]-S[13]*(m_l[13]-meq[13])+dt*F_hat[13];
//=======================================

m_l[14]=+0.00000000000000000*f[ci][0]+0.00000000000000000*f[ci][1]+0.00000000000000000*f[ci][2]+0.00000000000000000*f[ci][3]+0.00000000000000000*f[ci][4]+0.00000000000000000*f[ci][5]+0.00000000000000000*f[ci][6]+0.00000000000000000*f[ci][7]+0.00000000000000000*f[ci][8]+0.00000000000000000*f[ci][9]+0.00000000000000000*f[ci][10]+1.00000000000000000*f[ci][11]-1.00000000000000000*f[ci][12]-1.00000000000000000*f[ci][13]+1.00000000000000000*f[ci][14]+0.00000000000000000*f[ci][15]+0.00000000000000000*f[ci][16]+0.00000000000000000*f[ci][17]+0.00000000000000000*f[ci][18];

F_hat[14]=+0.00000000000000000*GuoF[0]+0.00000000000000000*GuoF[1]+0.00000000000000000*GuoF[2]+0.00000000000000000*GuoF[3]+0.00000000000000000*GuoF[4]+0.00000000000000000*GuoF[5]+0.00000000000000000*GuoF[6]+0.00000000000000000*GuoF[7]+0.00000000000000000*GuoF[8]+0.00000000000000000*GuoF[9]+0.00000000000000000*GuoF[10]+1.00000000000000000*GuoF[11]-1.00000000000000000*GuoF[12]-1.00000000000000000*GuoF[13]+1.00000000000000000*GuoF[14]+0.00000000000000000*GuoF[15]+0.00000000000000000*GuoF[16]+0.00000000000000000*GuoF[17]+0.00000000000000000*GuoF[18];

meq[14]=+0.00000000000000000*f_eq[0]+0.00000000000000000*f_eq[1]+0.00000000000000000*f_eq[2]+0.00000000000000000*f_eq[3]+0.00000000000000000*f_eq[4]+0.00000000000000000*f_eq[5]+0.00000000000000000*f_eq[6]+0.00000000000000000*f_eq[7]+0.00000000000000000*f_eq[8]+0.00000000000000000*f_eq[9]+0.00000000000000000*f_eq[10]+1.00000000000000000*f_eq[11]-1.00000000000000000*f_eq[12]-1.00000000000000000*f_eq[13]+1.00000000000000000*f_eq[14]+0.00000000000000000*f_eq[15]+0.00000000000000000*f_eq[16]+0.00000000000000000*f_eq[17]+0.00000000000000000*f_eq[18];

F_hat[14]*=(1-0.5*S[14]);
m_l[14]=m_l[14]-S[14]*(m_l[14]-meq[14])+dt*F_hat[14];
//=======================================

m_l[15]=+0.00000000000000000*f[ci][0]+0.00000000000000000*f[ci][1]+0.00000000000000000*f[ci][2]+0.00000000000000000*f[ci][3]+0.00000000000000000*f[ci][4]+0.00000000000000000*f[ci][5]+0.00000000000000000*f[ci][6]+0.00000000000000000*f[ci][7]+0.00000000000000000*f[ci][8]+0.00000000000000000*f[ci][9]+0.00000000000000000*f[ci][10]+0.00000000000000000*f[ci][11]+0.00000000000000000*f[ci][12]+0.00000000000000000*f[ci][13]+0.00000000000000000*f[ci][14]+1.00000000000000000*f[ci][15]-1.00000000000000000*f[ci][16]-1.00000000000000000*f[ci][17]+1.00000000000000000*f[ci][18];

F_hat[15]=+0.00000000000000000*GuoF[0]+0.00000000000000000*GuoF[1]+0.00000000000000000*GuoF[2]+0.00000000000000000*GuoF[3]+0.00000000000000000*GuoF[4]+0.00000000000000000*GuoF[5]+0.00000000000000000*GuoF[6]+0.00000000000000000*GuoF[7]+0.00000000000000000*GuoF[8]+0.00000000000000000*GuoF[9]+0.00000000000000000*GuoF[10]+0.00000000000000000*GuoF[11]+0.00000000000000000*GuoF[12]+0.00000000000000000*GuoF[13]+0.00000000000000000*GuoF[14]+1.00000000000000000*GuoF[15]-1.00000000000000000*GuoF[16]-1.00000000000000000*GuoF[17]+1.00000000000000000*GuoF[18];

meq[15]=+0.00000000000000000*f_eq[0]+0.00000000000000000*f_eq[1]+0.00000000000000000*f_eq[2]+0.00000000000000000*f_eq[3]+0.00000000000000000*f_eq[4]+0.00000000000000000*f_eq[5]+0.00000000000000000*f_eq[6]+0.00000000000000000*f_eq[7]+0.00000000000000000*f_eq[8]+0.00000000000000000*f_eq[9]+0.00000000000000000*f_eq[10]+0.00000000000000000*f_eq[11]+0.00000000000000000*f_eq[12]+0.00000000000000000*f_eq[13]+0.00000000000000000*f_eq[14]+1.00000000000000000*f_eq[15]-1.00000000000000000*f_eq[16]-1.00000000000000000*f_eq[17]+1.00000000000000000*f_eq[18];

F_hat[15]*=(1-0.5*S[15]);
m_l[15]=m_l[15]-S[15]*(m_l[15]-meq[15])+dt*F_hat[15];
//=======================================

m_l[16]=+0.00000000000000000*f[ci][0]+0.00000000000000000*f[ci][1]+0.00000000000000000*f[ci][2]+0.00000000000000000*f[ci][3]+0.00000000000000000*f[ci][4]+0.00000000000000000*f[ci][5]+0.00000000000000000*f[ci][6]+1.00000000000000000*f[ci][7]-1.00000000000000000*f[ci][8]+1.00000000000000000*f[ci][9]-1.00000000000000000*f[ci][10]+0.00000000000000000*f[ci][11]+0.00000000000000000*f[ci][12]+0.00000000000000000*f[ci][13]+0.00000000000000000*f[ci][14]-1.00000000000000000*f[ci][15]+1.00000000000000000*f[ci][16]-1.00000000000000000*f[ci][17]+1.00000000000000000*f[ci][18];

F_hat[16]=+0.00000000000000000*GuoF[0]+0.00000000000000000*GuoF[1]+0.00000000000000000*GuoF[2]+0.00000000000000000*GuoF[3]+0.00000000000000000*GuoF[4]+0.00000000000000000*GuoF[5]+0.00000000000000000*GuoF[6]+1.00000000000000000*GuoF[7]-1.00000000000000000*GuoF[8]+1.00000000000000000*GuoF[9]-1.00000000000000000*GuoF[10]+0.00000000000000000*GuoF[11]+0.00000000000000000*GuoF[12]+0.00000000000000000*GuoF[13]+0.00000000000000000*GuoF[14]-1.00000000000000000*GuoF[15]+1.00000000000000000*GuoF[16]-1.00000000000000000*GuoF[17]+1.00000000000000000*GuoF[18];

meq[16]=+0.00000000000000000*f_eq[0]+0.00000000000000000*f_eq[1]+0.00000000000000000*f_eq[2]+0.00000000000000000*f_eq[3]+0.00000000000000000*f_eq[4]+0.00000000000000000*f_eq[5]+0.00000000000000000*f_eq[6]+1.00000000000000000*f_eq[7]-1.00000000000000000*f_eq[8]+1.00000000000000000*f_eq[9]-1.00000000000000000*f_eq[10]+0.00000000000000000*f_eq[11]+0.00000000000000000*f_eq[12]+0.00000000000000000*f_eq[13]+0.00000000000000000*f_eq[14]-1.00000000000000000*f_eq[15]+1.00000000000000000*f_eq[16]-1.00000000000000000*f_eq[17]+1.00000000000000000*f_eq[18];

F_hat[16]*=(1-0.5*S[16]);
m_l[16]=m_l[16]-S[16]*(m_l[16]-meq[16])+dt*F_hat[16];
//=======================================

m_l[17]=+0.00000000000000000*f[ci][0]+0.00000000000000000*f[ci][1]+0.00000000000000000*f[ci][2]+0.00000000000000000*f[ci][3]+0.00000000000000000*f[ci][4]+0.00000000000000000*f[ci][5]+0.00000000000000000*f[ci][6]-1.00000000000000000*f[ci][7]-1.00000000000000000*f[ci][8]+1.00000000000000000*f[ci][9]+1.00000000000000000*f[ci][10]+1.00000000000000000*f[ci][11]-1.00000000000000000*f[ci][12]+1.00000000000000000*f[ci][13]-1.00000000000000000*f[ci][14]+0.00000000000000000*f[ci][15]+0.00000000000000000*f[ci][16]+0.00000000000000000*f[ci][17]+0.00000000000000000*f[ci][18];

F_hat[17]=+0.00000000000000000*GuoF[0]+0.00000000000000000*GuoF[1]+0.00000000000000000*GuoF[2]+0.00000000000000000*GuoF[3]+0.00000000000000000*GuoF[4]+0.00000000000000000*GuoF[5]+0.00000000000000000*GuoF[6]-1.00000000000000000*GuoF[7]-1.00000000000000000*GuoF[8]+1.00000000000000000*GuoF[9]+1.00000000000000000*GuoF[10]+1.00000000000000000*GuoF[11]-1.00000000000000000*GuoF[12]+1.00000000000000000*GuoF[13]-1.00000000000000000*GuoF[14]+0.00000000000000000*GuoF[15]+0.00000000000000000*GuoF[16]+0.00000000000000000*GuoF[17]+0.00000000000000000*GuoF[18];

meq[17]=+0.00000000000000000*f_eq[0]+0.00000000000000000*f_eq[1]+0.00000000000000000*f_eq[2]+0.00000000000000000*f_eq[3]+0.00000000000000000*f_eq[4]+0.00000000000000000*f_eq[5]+0.00000000000000000*f_eq[6]-1.00000000000000000*f_eq[7]-1.00000000000000000*f_eq[8]+1.00000000000000000*f_eq[9]+1.00000000000000000*f_eq[10]+1.00000000000000000*f_eq[11]-1.00000000000000000*f_eq[12]+1.00000000000000000*f_eq[13]-1.00000000000000000*f_eq[14]+0.00000000000000000*f_eq[15]+0.00000000000000000*f_eq[16]+0.00000000000000000*f_eq[17]+0.00000000000000000*f_eq[18];

F_hat[17]*=(1-0.5*S[17]);
m_l[17]=m_l[17]-S[17]*(m_l[17]-meq[17])+dt*F_hat[17];
//=======================================

m_l[18]=+0.00000000000000000*f[ci][0]+0.00000000000000000*f[ci][1]+0.00000000000000000*f[ci][2]+0.00000000000000000*f[ci][3]+0.00000000000000000*f[ci][4]+0.00000000000000000*f[ci][5]+0.00000000000000000*f[ci][6]+0.00000000000000000*f[ci][7]+0.00000000000000000*f[ci][8]+0.00000000000000000*f[ci][9]+0.00000000000000000*f[ci][10]-1.00000000000000000*f[ci][11]-1.00000000000000000*f[ci][12]+1.00000000000000000*f[ci][13]+1.00000000000000000*f[ci][14]+1.00000000000000000*f[ci][15]+1.00000000000000000*f[ci][16]-1.00000000000000000*f[ci][17]-1.00000000000000000*f[ci][18];

F_hat[18]=+0.00000000000000000*GuoF[0]+0.00000000000000000*GuoF[1]+0.00000000000000000*GuoF[2]+0.00000000000000000*GuoF[3]+0.00000000000000000*GuoF[4]+0.00000000000000000*GuoF[5]+0.00000000000000000*GuoF[6]+0.00000000000000000*GuoF[7]+0.00000000000000000*GuoF[8]+0.00000000000000000*GuoF[9]+0.00000000000000000*GuoF[10]-1.00000000000000000*GuoF[11]-1.00000000000000000*GuoF[12]+1.00000000000000000*GuoF[13]+1.00000000000000000*GuoF[14]+1.00000000000000000*GuoF[15]+1.00000000000000000*GuoF[16]-1.00000000000000000*GuoF[17]-1.00000000000000000*GuoF[18];

meq[18]=+0.00000000000000000*f_eq[0]+0.00000000000000000*f_eq[1]+0.00000000000000000*f_eq[2]+0.00000000000000000*f_eq[3]+0.00000000000000000*f_eq[4]+0.00000000000000000*f_eq[5]+0.00000000000000000*f_eq[6]+0.00000000000000000*f_eq[7]+0.00000000000000000*f_eq[8]+0.00000000000000000*f_eq[9]+0.00000000000000000*f_eq[10]-1.00000000000000000*f_eq[11]-1.00000000000000000*f_eq[12]+1.00000000000000000*f_eq[13]+1.00000000000000000*f_eq[14]+1.00000000000000000*f_eq[15]+1.00000000000000000*f_eq[16]-1.00000000000000000*f_eq[17]-1.00000000000000000*f_eq[18];

F_hat[18]*=(1-0.5*S[18]);
m_l[18]=m_l[18]-S[18]*(m_l[18]-meq[18])+dt*F_hat[18];
//=======================================






			//==============================================================================




			// ==================   f=M_-1m matrix calculation and streaming =============================
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




void comput_macro_variables( double* rho,double** u,double** u0,double** f,double** F,double* PerC, double* PorC,int* SupInv,int*** Solid)
{
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();


	double c0,c1,uu;



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
				

				u[i][0]=(u[i][0]+dt*PorC[i]*rho[i]*gx/2)/rho[i];
				u[i][1]=(u[i][1]+dt*PorC[i]*rho[i]*gy/2)/rho[i];
				u[i][2]=(u[i][2]+dt*PorC[i]*rho[i]*gz/2)/rho[i];

				c0=0.5*(1+PorC[i]*dt/2*in_vis/PerC[i]);
				c1=PorC[i]*dt/2*F_epsilon/(sqrt(PerC[i]));

				uu=sqrt(u[i][0]*u[i][0]+u[i][1]*u[i][1]+u[i][2]*u[i][2]);
				u[i][0]=u[i][0]/(c0+sqrt(c0*c0+c1*uu));
				u[i][1]=u[i][1]/(c0+sqrt(c0*c0+c1*uu));
				u[i][2]=u[i][2]/(c0+sqrt(c0*c0+c1*uu));
				
		
			}


	
	//MPI_Barrier(MPI_COMM_WORLD); 

}



void boundary_velocity(int xp,double v_xp,int xn, double v_xn,int yp,double v_yp,int yn,double v_yn,int zp,double v_zp,int zn,double v_zn,double** f,double** F,double* rho,double** u,int*** Solid,double* PorC)

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
		        F[Solid[i][NY][k]][ks]=feq(ks,1.0,u_yp,PorC[Solid[i][NY][k]]);
		if ((yn==1) && (Solid[i][0][k]>0))
		       F[Solid[i][0][k]][ks]=feq(ks,1.0,u_yn,PorC[Solid[i][0][k]]);      
		}

if ((zp-1)*(zn-1)==0)		
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
		for (int ks=0;ks<Q;ks++)
		{
		if ((zp==1) && (Solid[i][j][NZ]>0))    
		        F[Solid[i][j][NZ]][ks]=feq(ks,1.0,u_zp,PorC[Solid[i][j][NZ]]); 
		if ((zn==1) && (Solid[i][j][0]>0))
		        F[Solid[i][j][0]][ks]=feq(ks,1.0,u_zn,PorC[Solid[i][j][0]]);
		}


if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
			F[Solid[nx_l-1][j][k]][ks]=feq(ks,1.0,u_xp,PorC[Solid[nx_l-1][j][k]]);
			



if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
			F[Solid[0][j][k]][ks]=feq(ks,1.0,u_xn,PorC[Solid[0][j][k]]);
			
	

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
		                 F[Solid[i][NY][k]][ks]=feq(ks,1.0,u_yp,PorC[Solid[i][NY][k]]);
		if ((yn==1) && (Solid[i][0][k]>0))
		        //if (Solid[i][1][k]>0)
		        //         F[Solid[i][0][k]][ks]=feq(ks,rho[Solid[i][1][k]],u_yn);
		        // else
		                 F[Solid[i][0][k]][ks]=feq(ks,1.0,u_yn,PorC[Solid[i][0][k]]);      
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
		                F[Solid[i][j][NZ]][ks]=feq(ks,1.0,u_zp,PorC[Solid[i][j][NZ]]); 
		if ((zn==1) && (Solid[i][j][0]>0))
		        //if (Solid[i][j][1]>0)
		       //         F[Solid[i][j][0]][ks]=feq(ks,rho[Solid[i][j][1]],u_zn);
		       // else
		                F[Solid[i][j][0]][ks]=feq(ks,1.0,u_zn,PorC[Solid[i][j][0]]);
		}


if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
		        //if (Solid[nx_l-2][j][k]>0)
			//        F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho[Solid[nx_l-2][j][k]],u_xp);
			//else
			         F[Solid[nx_l-1][j][k]][ks]=feq(ks,1.0,u_xp,PorC[Solid[nx_l-1][j][k]]);
			



if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
		        //if (Solid[1][j][k]>0)
		        //        F[Solid[0][j][k]][ks]=feq(ks,rho[Solid[1][j][k]],u_xn);
		        //else
		                F[Solid[0][j][k]][ks]=feq(ks,1.0,u_xn,PorC[Solid[0][j][k]]);
			
	

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
				F[Solid[i][NY][k]][ks]=feq(LR[ks],1.0,u_yp,PorC[Solid[i][NY][k]])-F[Solid[i][NY][k]][LR[ks]]-feq(ks,1.0,u_yp,PorC[Solid[i][NY][k]]);

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
					F[Solid[i][0][k]][ks]=feq(LR[ks],1.0,u_yn,PorC[Solid[i][0][k]])-F[Solid[i][0][k]][LR[ks]]+feq(ks,1.0,u_yn,PorC[Solid[i][0][k]]);

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
					F[Solid[i][j][NZ]][ks]=feq(LR[ks],1.0,u_zp,PorC[Solid[i][j][NZ]])-F[Solid[i][j][NZ]][LR[ks]]+feq(ks,1.0,u_zp,PorC[Solid[i][j][NZ]]);
		
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
					F[Solid[i][j][0]][ks]=feq(LR[ks],1.0,u_zn,PorC[Solid[i][j][0]])-F[Solid[i][j][0]][LR[ks]]+feq(ks,1.0,u_zn,PorC[Solid[i][j][0]]);
	
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
					F[Solid[nx_l-1][j][k]][ks]=feq(LR[ks],1.0,u_xp,PorC[Solid[nx_l-1][j][k]])-F[Solid[nx_l-1][j][k]][LR[ks]]+feq(ks,1.0,u_xp,PorC[Solid[nx_l-1][j][k]]);
		
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
					F[Solid[0][j][k]][ks]=feq(LR[ks],1.0,u_xn,PorC[Solid[0][j][k]])-F[Solid[0][j][k]][LR[ks]]+feq(ks,1.0,u_xn,PorC[Solid[0][j][k]]);
		
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
				//	F[Solid[i][NY][k]][ks]=feq(ks,rho[Solid[i][NY-1][k]],u_yp,PorC[Solid[i][NY][k]])+f[Solid[i][NY-1][k]][ks]-feq(ks,rho[Solid[i][NY-1][k]],u[Solid[i][NY-1][k]],PorC[Solid[i][NY-1][k]]);
				//else
					F[Solid[i][NY][k]][ks]=feq(ks,1.0,u_yp,PorC[Solid[i][NY][k]]);

			}


if (yn==1)
for (int i=0;i<nx_l;i++)
	for(int k=0;k<=NZ;k++)
	if (Solid[i][0][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
				//if (Solid[i][1][k]>0) 
				//	F[Solid[i][0][k]][ks]=feq(ks,rho[Solid[i][1][k]],u_yn,PorC[Solid[i][0][k]])+f[Solid[i][1][k]][ks]-feq(ks,rho[Solid[i][1][k]],u[Solid[i][1][k]],PorC[Solid[i][1][k]]);
				//else
					F[Solid[i][0][k]][ks]=feq(ks,1.0,u_yn,PorC[Solid[i][0][k]]);

			}



if (zp==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][NZ]>0)
		for (int ks=0;ks<Q;ks++)
			{
				//if (Solid[i][j][NZ-1]>0) 
				//	F[Solid[i][j][NZ]][ks]=feq(ks,rho[Solid[i][j][NZ-1]],u_zp,PorC[Solid[i][j][NZ]])+f[Solid[i][j][NZ-1]][ks]-feq(ks,rho[Solid[i][j][NZ-1]],u[Solid[i][j][NZ-1]],PorC[Solid[i][j][NZ-1]]);
				//else
					F[Solid[i][j][NZ]][ks]=feq(ks,1.0,u_zp,PorC[Solid[i][j][NZ]]);
		
			}



if (zn==1)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
	if (Solid[i][j][0]>0)
		for (int ks=0;ks<Q;ks++)
			{
				//if (Solid[i][j][1]>0) 
				//	F[Solid[i][j][0]][ks]=feq(ks,rho[Solid[i][j][1]],u_zn,PorC[Solid[i][j][0]])+f[Solid[i][j][1]][ks]-feq(ks,rho[Solid[i][j][1]],u[Solid[i][j][1]],PorC[Solid[i][j][1]]);
				//else
					F[Solid[i][j][0]][ks]=feq(ks,1.0,u_zn,PorC[Solid[i][j][1]]);
	
			}




if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
	if (Solid[nx_l-1][k][j]>0)
		for (int ks=0;ks<Q;ks++)
			{
				//if (Solid[nx_l-2][j][k]>0) 
				//	F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho[Solid[nx_l-2][j][k]],u_xp,PorC[Solid[nx_l-1][j][k]])+f[Solid[nx_l-2][j][k]][ks]-feq(ks,rho[Solid[nx_l-2][j][k]],u[Solid[nx_l-2][j][k]],PorC[Solid[nx_l-2][j][k]]);
				//else
					F[Solid[nx_l-1][j][k]][ks]=feq(ks,1.0,u_xp,PorC[Solid[nx_l-1][j][k]]);
		
			}



if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
	if (Solid[0][j][k]>0)
		for (int ks=0;ks<Q;ks++)
			{
				//if (Solid[1][j][k]>0) 
				//	F[Solid[0][j][k]][ks]=feq(ks,rho[Solid[1][j][k]],u_xn,PorC[Solid[0][j][k]])+f[Solid[1][j][k]][ks]-feq(ks,rho[Solid[1][j][k]],u[Solid[1][j][k]],u_xn,PorC[Solid[1][j][k]]);
				//else
					F[Solid[0][j][k]][ks]=feq(ks,1.0,u_xn,PorC[Solid[0][j][k]]);
		
			}
	
}


}



void boundary_pressure(int xp,double rho_xp,int xn, double rho_xn,int yp,double rho_yp,int yn,double rho_yn,int zp,double rho_zp,int zn,double rho_zn,double** f,double** F,double** u,double* rho,int*** Solid, double* PorC)


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
		        F[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls,PorC[Solid[i][NY][k]]);
		if ((yn==1) && (Solid[i][0][k]>0))
		        F[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls,PorC[Solid[i][0][k]]);
		
		}
	        
if ((zp-1)*(zn-1)==0)
for (int i=0;i<nx_l;i++)
	for(int j=0;j<=NY;j++)
		for (int ks=0;ks<Q;ks++)
		{
		if ((zp==1) && (Solid[i][j][NZ]>0))
		        F[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls,PorC[Solid[i][j][NZ]]);
		if ((zn==1) && (Solid[i][j][0]>0))
		        F[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls,PorC[Solid[i][j][0]]);
		}

		
if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
        for (int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
		if (Solid[nx_l-1][j][k]>0)
			F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls,PorC[Solid[nx_l-1][j][k]]);
			

if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)		
		if (Solid[0][j][k]>0)
			F[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls,PorC[Solid[0][j][k]]);
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
					F[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls,PorC[Solid[i][NY][k]]);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls,PorC[Solid[i][NY][k]]);
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
					F[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls,PorC[Solid[i][0][k]]);
					}
			else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls,PorC[Solid[i][0][k]]);
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
					F[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls,PorC[Solid[i][j][NZ]]);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls,PorC[Solid[i][j][NZ]]);
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
					F[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls,PorC[Solid[i][j][0]]);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls,PorC[Solid[i][j][0]]);
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
					F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls,PorC[Solid[nx_l-1][j][k]]);
					}
				else	
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;
					u_ls[2]=0.0;
					F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls,PorC[Solid[nx_l-1][j][k]]);
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
					F[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls,PorC[Solid[0][j][k]]);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;
					u_ls[2]=0.0;
					F[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls,PorC[Solid[0][j][k]]);
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
						m_l[mi]+=M[mi][mj]*f[Solid[i][NY-1][k]][mj];					}

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
				F[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls,PorC[Solid[i][NY][k]]);
			
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
						m_l[mi]+=M[mi][mj]*f[Solid[i][1][k]][mj];					}

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
				F[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls,PorC[Solid[i][0][k]]);
			
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
						m_l[mi]+=M[mi][mj]*f[Solid[i][j][NZ-1]][mj];					}

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
				F[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls,PorC[Solid[i][j][NZ]]);
			
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
						m_l[mi]+=M[mi][mj]*f[Solid[i][j][1]][mj];					}

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
				F[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls,PorC[Solid[i][j][0]]);
			
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
						m_l[mi]+=M[mi][mj]*f[Solid[nx_l-2][j][k]][mj];					}

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
				F[Solid[nx_l-1][j][NZ]][ks]=feq(ks,rho_xp,u_ls,PorC[Solid[nx_l-1][j][NZ]]);
			
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
						m_l[mi]+=M[mi][mj]*f[Solid[1][j][k]][mj];					}

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
				F[Solid[0][j][NZ]][ks]=feq(ks,rho_xn,u_ls,PorC[Solid[0][j][NZ]]);
			
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
					F[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls,PorC[Solid[i][NY][k]])+f[Solid[i][NY-1][k]][ks]-feq(ks,rho[Solid[i][NY-1][k]],u[Solid[i][NY-1][k]],PorC[Solid[i][NY-1][k]]);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[i][NY][k]][ks]=feq(ks,rho_yp,u_ls,PorC[Solid[i][NY][k]]);
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
					F[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls,PorC[Solid[i][0][k]])+f[Solid[i][1][k]][ks]-feq(ks,rho[Solid[i][1][k]],u[Solid[i][1][k]],PorC[Solid[i][1][k]]);
					}
			else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[i][0][k]][ks]=feq(ks,rho_yn,u_ls,PorC[Solid[i][0][k]]);
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
				F[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls,PorC[Solid[i][j][NZ]])+f[Solid[i][j][NZ-1]][ks]-feq(ks,rho[Solid[i][j][NZ-1]],u[Solid[i][j][NZ-1]],PorC[Solid[i][j][NZ-1]]);
				}
				else
				{
				u_ls[0]=0.0;
				u_ls[1]=0.0;u_ls[2]=0.0;
				F[Solid[i][j][NZ]][ks]=feq(ks,rho_zp,u_ls,PorC[Solid[i][j][NZ]]);
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
			F[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls,PorC[Solid[i][j][0]])+f[Solid[i][j][1]][ks]-feq(ks,rho[Solid[i][j][1]],u[Solid[i][j][1]],PorC[Solid[i][j][1]]);
			}
			else
			{
			u_ls[0]=0.0;
			u_ls[1]=0.0;u_ls[2]=0.0;
			F[Solid[i][j][0]][ks]=feq(ks,rho_zn,u_ls,PorC[Solid[i][j][0]]);
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
					F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls,PorC[Solid[nx_l-1][j][k]])+f[Solid[nx_l-2][j][k]][ks]-feq(ks,rho[Solid[nx_l-2][j][k]],u[Solid[nx_l-2][j][k]],PorC[Solid[nx_l-1][j][k]]);
					//F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls);
					}
				else	
					{
					u_ls[0]=0.0;u_ls[1]=0.0;u_ls[2]=0.0;
					F[Solid[nx_l-1][j][k]][ks]=feq(ks,rho_xp,u_ls,PorC[Solid[nx_l-1][j][k]]);
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
					F[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls,PorC[Solid[0][j][k]])+f[Solid[1][j][k]][ks]-feq(ks,rho[Solid[1][j][k]],u[Solid[1][j][k]],PorC[Solid[1][j][k]]);
					//f[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls);
					}
				else
					{
					u_ls[0]=0.0;
					u_ls[1]=0.0;
					u_ls[2]=0.0;
					F[Solid[0][j][k]][ks]=feq(ks,rho_xn,u_ls,PorC[Solid[0][j][k]]);
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




void Backup_init(double* rho, double** u, double** f, double* PerC, double* PorC, double* forcex, double* forcey, double* forcez, char backup_rho[128], char backup_velocity[128], char backup_f[128],int* SupInv)
{	
      
  
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	int* nx_g = new int[mpi_size];
	int* disp = new int[mpi_size];
	int ip,jp,kp,ijk;
	
	MPI_Gather(&nx_l,1,MPI_INT,nx_g,1,MPI_INT,0,MPI_COMM_WORLD);
	
	if (rank==0)
	{
		disp[0]=0;
		

		for (int i=1;i<mpi_size;i++)
			disp[i]=disp[i-1]+nx_g[i-1];
	}

	MPI_Bcast(disp,mpi_size,MPI_INT,0,MPI_COMM_WORLD);

	
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
	double uu;


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
        	        fin >> f[i][0] >> f[i][1] >> f[i][2] >> f[i][3] >> f[i][4] >> f[i][5] >>f[i][6] >> f[i][7] >> f[i][8] >> f[i][9] >> f[i][10] >> f[i][11] >> f[i][12] >> f[i][13] >> f[i][14] >> f[i][15] >> f[i][16] >>f[i][17] >> f[i][18] ;
       fin.close();
       
     MPI_Barrier(MPI_COMM_WORLD);
			
	for (int i=1;i<=Count;i++)	
				
		{ 	
				
			

			ip=(int)(SupInv[i]/((NY+1)*(NZ+1)));
			jp=(int)((SupInv[i]%((NY+1)*(NZ+1)))/(NZ+1));
			kp=(int)(SupInv[i]%(NZ+1));

			ip=(ip-ip%Zoom)/Zoom;
			jp=(jp-jp%Zoom)/Zoom;
			kp=(kp-kp%Zoom)/Zoom;
 			

			PerC[i]=Per_Int[ip*(NY+1)*(NZ+1)+jp*(NZ+1)+kp];
			PorC[i]=Por_Int[ip*(NY+1)*(NZ+1)+jp*(NZ+1)+kp];
			
			

			
			//***********************************************************************

			
			uu=sqrt(u[i][0]*u[i][0]+u[i][1]*u[i][1]+u[i][2]*u[i][2]);
			forcex[i]=-PorC[i]*niu/PerC[i]*u[i][0]-PorC[i]*F_epsilon*u[i][0]*uu/(sqrt(PerC[i]))+PorC[i]*gx;
			forcey[i]=-PorC[i]*niu/PerC[i]*u[i][1]-PorC[i]*F_epsilon*u[i][1]*uu/(sqrt(PerC[i]))+PorC[i]*gy;
			forcez[i]=-PorC[i]*niu/PerC[i]*u[i][2]-PorC[i]*F_epsilon*u[i][2]*uu/(sqrt(PerC[i]))+PorC[i]*gz;







			//***********************************************************************


			
				

		
		

	}

	
	delete [] nx_g;
	delete [] disp;

	 	
}
