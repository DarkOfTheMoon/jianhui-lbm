#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
# include "mpi.h"


#define MAXN 100
#define eps 1e-12
#define zero(x) (fabs(x)<eps)

struct mat{
    int n,m;
    double data[MAXN][MAXN];
};



//D2Q9 STANDARD LATTICE MRT LATTICE BOLTZMANN METHOD
//UNIFORM MESH, LATTICE VELOCITY 1


using namespace std;        


void init(double***, double****, double****,double[19],int[19][3],bool***, double***, double***, double***, double[19], int, int, int,double,double,double,double,int,int,int,double*,char[128],double,double,double);

void cylinder_creation(int, int, double);

void collision(double***, double****, double****,double[19],int[19][3],double[19][19], double[19][19], double[19],double***, double***, double***, int, int, int,bool***);
   
double Comput_Perm(double****,double* ,int, int ,int ,double ,double ,double ,double ,int ,int , int,double,int );

void periodic_streaming(double****,double****,bool***,int[19][3],int,int,int);

void comput_macro_variables_IMR(double***, double****,double****,double****,int[19][3],bool***, double***, double***, double***, int,int, int,int,double,double*,double*,double*);

void comput_macro_variables(double***, double****,double****,double****,int[19][3],bool***, double***, double***, double***, int,int, int);

void standard_bounceback_boundary(int,int,int, double****);

void evolution(double***, double****,double****,double****);

double Error(double****,double****,int,int,int,double*,double*);

void output_velocity(int, double****,int,int,int);

void output_velocity_Vector(int, double****,int,int,int,int,int,int,int);

void output_density(int, double***,int,int,int,int,int,int,int);

void boundary_pressure(int,double,int, double,int,double,int,double ,int,double,int,double,double****,double****,double***,int, int,int, double[19],int[19][3],bool***);

void boundary_velocity(int,double,int, double,int,double,int,double ,int,double,int,double,double****,double***,int, int,int, double[19],int[19][3],bool***);

void Read_Rock(bool*** ,int,int,int,int,int,int,double*,char[128]);



double S[19];
void Comput_MI(double[19][19], double[19][19]);
int inverse(mat &a);
double feq(int,double, double[3],double[19],int[19][3]);
void Geometry(bool***, int,int,int ,int,int,int,int);	
int Par_Geo,Par_nx,Par_ny,Par_nz,Zoom;

int PerDir;

char outputfile[128]="./";

int main (int argc , char *argv [])
{
MPI :: Init (argc , argv );
MPI_Status status ;

int rank = MPI :: COMM_WORLD . Get_rank ();
int para_size=MPI :: COMM_WORLD . Get_size ();

int NX,NY,NZ,freRe,freVe,freDe;
double gx,gy,gz,reso;

int mirX=0;
int mirY=0;
int mirZ=0;
int mir=0;
const int NCHAR=128;
int in_IMR,wr_per,pre_xp,pre_xn,pre_yp,pre_yn,pre_zp,pre_zn;
int vel_xp,vel_xn,vel_yp,vel_yn,vel_zp,vel_zn;
double in_vis,p_xp,p_xn,p_yp,p_yn,p_zp,p_zn;
double inivx,inivy,inivz,v_xp,v_xn,v_yp,v_yn,v_zp,v_zn;
int vtk_time,n_max;


	char     filename[NCHAR], dummy[NCHAR+1];
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
	fin.close();
	

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
	








const int Q=19;          
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



if ( (NX+1)%(MPI :: COMM_WORLD . Get_size ())!=0) 
	{
	cout<<NX+1<<"  "<<MPI :: COMM_WORLD . Get_size ()<<endl;
	cout<<"DOMAIN SCALE CAN NOT BE DIVIDED BY THE NUMBER OF PROCESSORS"<<endl;
	exit (0);
	}

	
double porosity;	




double s_v=2.0/(6*in_vis+1);
//double s_v=2.0/(6*0.1+1.0);
double v_max;

double u_ave;

/*
double M[19][19]=
{{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
{-30,-11,-11,-11,-11,-11,-11,8,8,8,8,8,8,8,8,8,8,8,8},
{12,-4,-4,-4,-4,-4,-4,1,1,1,1,1,1,1,1,1,1,1,1},
{0,1,-1,0,0,0,0,1,-1,1,-1,0,0,0,0,1,-1,1,-1},
{0,-4,4,0,0,0,0,1,-1,1,-1,0,0,0,0,1,-1,1,-1},
{0,0,0,1,-1,0,0,1,1,-1,-1,1,-1,1,-1,0,0,0,0},
{0,0,0,-4,4,0,0,4,4,-4,-4,4,-4,4,-4,0,0,0,0},
{0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1},
{0,0,0,0,0,-4,4,0,0,0,0,1,1,-1,-1,1,1,-1,-1},
{0,2,2,-1,-1,-1,-1,1,1,1,1,-2,-2,-2,-2,1,1,1,1},
{0,-4,-4,5,5,5,5,4,4,4,4,-2,-2,-2,-2,1,1,1,1},
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
double* Permia;
double m[19];
double meq[19];

double uMax=0.08;
double Re=300.0;
//forcex=gx;forcey=gy;
//niu=Umax*(NX+1)/Re;
double niu=uMax*(NY+1)/Re;
double tau_f=3.0*niu+0.5;
//tau_f=1;
//s_v=1.0/tau_f;
//s_v=0.6;

if (Zoom>1)	
	reso=reso/Zoom;


int n;
int e[19][3]=
{{0,0,0},{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};
double w[19]={1.0/3.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};

	double error;
	double error_perm;

	double ****u;
	double ****u0;
	double ****f;
	double ***rho;
	double ****F;
	bool*** Solid;
	double*** forcex;
	double*** forcey;
	double*** forcez;

  // Allocate memory
	Permia = new double[3];
  	u = new double***[(NX+1)/para_size];
	u0 = new double***[(NX+1)/para_size];
	f = new double***[(NX+1)/para_size];
	rho = new double**[(NX+1)/para_size];
	F = new double***[(NX+1)/para_size];
	Solid = new bool**[(NX+1)/para_size];
	forcex = new double**[(NX+1)/para_size];
	forcey = new double**[(NX+1)/para_size];
	forcez = new double**[(NX+1)/para_size];

  for (int i = 0; i < (NX+1)/para_size; i++) 
{

    	u[i] = new double**[NY+1];
	u0[i] = new double**[NY+1];
	f[i] = new double**[NY+1];
	F[i] = new double**[NY+1];
	rho[i]=new double*[NY+1];
	Solid[i]=new bool*[NY+1];
	forcex[i] = new double*[NY+1];
	forcey[i] = new double*[NY+1];
	forcez[i] = new double*[NY+1];



	
    for (int j = 0; j < NY+1; j++)
	{
      	u[i][j] = new double*[NZ+1];
	u0[i][j] = new double*[NZ+1];
	f[i][j] = new double*[NZ+1];
	rho[i][j]= new double[NZ+1];
	F[i][j] = new double*[NZ+1];
	Solid[i][j] = new bool[NZ+1];
	forcex[i][j]= new double[NZ+1];
	forcey[i][j]= new double[NZ+1];
	forcez[i][j]= new double[NZ+1];


	for (int k=0;k<NZ+1;k++)
		{
		u[i][j][k] = new double[3];
		u0[i][j][k] = new double[3];
		f[i][j][k] = new double[19];
		F[i][j][k] = new double[19];
		}
	}
  }


Comput_MI(M,MI);



	using namespace std;
	init(rho,u,f,w,e,Solid,forcex,forcey,forcez,S,NX,NY,NZ,s_v,gx,gy,gz,mirX,mirY,mirZ,&porosity,filename,inivx,inivy,inivz);
	Geometry(Solid,NX,NY,NZ,mirX,mirY,mirZ,mir);
	
	if (rank==0)
		cout<<"Porosity= "<<porosity<<endl;

char FileName[128]="Results.txt";
char FileName2[128]="Permeability.txt";
char FileName3[128]="bodyforce.txt";

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

	collision(rho,u,f,w,e,M,MI,S,forcex,forcey,forcez,NX,NY,NZ,Solid);

	periodic_streaming(f,F,Solid,e,NX,NY,NZ);

	if ((1-pre_xp)*(1-pre_xn)*(1-pre_yp)*(1-pre_yn)*(1-pre_zp)*(1-pre_zn)==0)
		boundary_pressure(pre_xp,p_xp,pre_xn,p_xn,pre_yp,p_yp,pre_yn,p_yn,pre_zp,p_zp,pre_zn,p_zn,f,u,rho,NX,NY,NZ,w,e,Solid);
	
	if ((1-vel_xp)*(1-vel_xn)*(1-vel_yp)*(1-vel_yn)*(1-vel_zp)*(1-vel_zn)==0)
		boundary_velocity(vel_xp,v_xp,vel_xn,v_xn,vel_yp,v_yp,vel_yn,v_yn,vel_zp,v_zp,vel_zn,v_zn,f,rho,NX,NY,NZ,w,e,Solid);

	if (in_IMR==1)
  		comput_macro_variables_IMR(rho,u,u0,f,e,Solid,forcex,forcey,forcez,NX,NY,NZ,n,porosity,&gx,&gy,&gz);
	else
		comput_macro_variables(rho,u,u0,f,e,Solid,forcex,forcey,forcez,NX,NY,NZ);

		if(n%freRe==0)
		{       
			error=Error(u,u0,NX,NY,NZ,&v_max,&u_ave);if (v_max>=1.0)	U_max_ref+=1;
			error_perm=Comput_Perm(u,Permia,NX,NY,NZ,s_v,gx,gy,gz,mirX,mirY,mirZ,porosity,PerDir);
			if (rank==0)
			{
			ofstream fin(FileName,ios::app);
			fin<<"The"<<n<<"th computation result:"<<endl;
		//=============================================================================================
			fin<<"The permiability is: "<<Permia[0]*reso*reso*100<<", "<<Permia[1]*reso*reso*100<<", "<<Permia[2]*reso*reso*100<<endl;
			fin<<"The relative error of permiability computing is: "<<error_perm<<endl;
		//==============================================================================================

		//==============================================================================================
			//cout<<"The Density of point(NX/2,NY/2,NZ/2) is: "<<setprecision(6)
			//	<<rho[int((NX+1)/para_size/2)][NY/2][NZ/2]<<endl;
			Re=u_ave*(NY+1)/(1.0/3.0*(1/s_v-0.5));
			fin<<"The Maximum velocity is: "<<setprecision(6)<<v_max<<"   Re="<<Re<<endl;
			
		//===============================================================================================
			fin<<"The max relative error of velocity is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			fin<<endl;
			fin.close();


		if (wr_per==1)
			{
			ofstream finfs(FileName2,ios::app);
			switch(PerDir)
				{
				case 1:
				finfs<<Permia[0]*reso*reso*100<<endl;break;
				case 2:
				finfs<<Permia[1]*reso*reso*100<<endl;break;
				case 3:
				finfs<<Permia[2]*reso*reso*100<<endl;break;
				default:
				finfs<<Permia[0]*reso*reso*100<<endl;break;
				}
			finfs.close();
			}

		
			ofstream finf3(FileName3,ios::app);
			finf3<<gx<<endl;
			finf3.close();
		

			cout<<"The"<<n<<"th computation result:"<<endl;
		//=============================================================================================
			cout<<"The permiability is: "<<Permia[0]*reso*reso*100<<", "<<Permia[1]*reso*reso*100<<", "<<Permia[2]*reso*reso*100<<endl;
			cout<<"The relative error of permiability computing is: "<<error_perm<<endl;
		//==============================================================================================

		//==============================================================================================
			//cout<<"The Density of point(NX/2,NY/2,NZ/2) is: "<<setprecision(6)
			//	<<rho[int((NX+1)/para_size/2)][NY/2][NZ/2]<<endl;
			Re=u_ave*(NY+1)/(1.0/3.0*(1/s_v-0.5));
			cout<<"The Maximum velocity is: "<<setprecision(6)<<v_max<<"   Re="<<Re<<endl;
			
		//===============================================================================================
			cout<<"The max relative error of uv is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			cout<<endl;
                        }

			if ((freDe>=0) and (n%freDe==0))
			output_density(n,rho,NX,NY,NZ,mirX,mirY,mirZ,mir);
			if ((freVe>=0) and (n%freVe==0))
			output_velocity_Vector(n,u,NX,NY,NZ,mirX,mirY,mirZ,mir);
						

				if(error!=error) {cout<<"PROGRAM STOP"<<endl;break;};
				if(U_max_ref>=5) {cout<<"PROGRAM STOP DUE TO HIGH VELOCITY"<<endl;break;}
			
		}
	}




 for (int i = 0; i < (NX+1)/para_size; i++) 
	{
    	for (int j = 0; j < NY+1; j++)
		{
		for (int k=0;k<NZ+1;k++)
      			{
			delete [] u[i][j][k];
			delete [] u0[i][j][k];
			delete [] f[i][j][k];
			delete [] F[i][j][k];
			}

    		delete [] u[i][j];
		delete [] u0[i][j];
		delete [] f[i][j];
		delete [] F[i][j];
		delete [] rho[i][j];
		delete [] Solid[i][j];
		delete [] forcex[i][j];
		delete [] forcey[i][j];
		delete [] forcez[i][j];


  		}
  	delete [] u[i];
	delete [] u0[i];
	delete [] f[i];
	delete [] F[i];
	delete [] rho[i];
	delete [] Solid[i];
	delete [] forcex[i];
	delete [] forcey[i];
	delete [] forcez[i];	





	}

	delete [] u;
	delete [] u0;
	delete [] f;
	delete [] F;
	delete [] rho;
	delete [] Solid;
	delete [] forcex;
	delete [] forcey;
	delete [] forcez;	

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




double feq(int k,double rho, double u[3],double w[19],int e[19][3])
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


void init(double*** rho, double**** u, double**** f,double w[19], int e[19][3],bool*** Solid, double*** forcex, double*** forcey, double*** forcez, double S[19], int NX,int NY, int NZ, double s_v, double gx,double gy,double gz,int mirX,int mirY,int mirZ,double* porosity,char poreFileName[128],double inivx,double inivy,double inivz)
{	
        //=================
        //double res[9][9];
        //=================
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	//double Cylinder_r=8;
	double usqr,vsqr;
	double dx=1.0;
	double dy=1.0;
	double Lx=dx*double(NX);
	double Ly=dy*double(NY);
	double dt=1;//sqrt(3);
	double c=dx/dt;
	double rho0=1.0;
 	double uMax=0.00;
	
       
	
	double pr; //raduis of the obstacles
        

	//srand(time(0));                        //RANDOM NUMBER GENERATION SEED
	//cylinder_creation(50,103,Cylinder_r);

	double u_tmp[3];

	double s_other=8*(2-s_v)/(8-s_v);
	
	//s_other=s_v;	

	
/*
	S[0]=0;
	S[1]=s_other;
	S[2]=s_other;
	S[3]=0;
	S[4]=s_other;
	S[5]=0;
	S[6]=s_other;
	S[7]=0;
	S[8]=s_other;
	S[9]=s_v;
	S[10]=s_other;
	S[11]=s_v;
	S[12]=s_other;
	S[13]=s_v;
	S[14]=s_v;
	S[15]=s_v;
	S[16]=s_other;
	S[17]=s_other;
	S[18]=s_other;
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

	if (!(poreFileName=="NOSOLID"))
		Read_Rock(Solid,NX,NY,NZ,mirX,mirY,mirZ,porosity,poreFileName);
	else
		{
		for(int i=0;i<(NX+1)/mpi_size;i++)	
			for(int j=0;j<=NY;j++)
				for(int k=0;k<=NZ;k++)
				Solid[i][j][k]=0;
		}

	
for(int i=0;i<(NX+1)/mpi_size;i++)	
	for(int j=0;j<=NY;j++)
		for(int k=0;k<=NZ;k++)		
		{    	
			

                        
			u[i][j][k][0]=inivx;
			u[i][j][k][1]=inivy;
			u[i][j][k][2]=inivz;
			u_tmp[0]=u[i][j][k][0];
			u_tmp[1]=u[i][j][k][1];
			u_tmp[2]=u[i][j][k][2];

			rho[i][j][k]=1.0;
			
			
			

			//***********************************************************************

			forcex[i][j][k]=gx;
			forcey[i][j][k]=gy;
			forcez[i][j][k]=gz;


			




			//***********************************************************************


			//INITIALIZATION OF m and f

			for (int lm=0;lm<19;lm++)
				{
				f[i][j][k][lm]=feq(lm,rho[i][j][k],u_tmp,w,e);
				//if (Solid[i][j][k]) {f[i][j][k][lm]=0.0;} 
				}

		
		

	}

		

	 	
}





void periodic_streaming(double**** f,double**** F,bool*** Solid,int e[19][3],int NX,int NY,int NZ)
{

	MPI_Status status[4] ;
	MPI_Request request[4];

	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	int Buf_Size=(NY+1)*(NZ+1)*19;
	
	double* swapl= new double[Buf_Size];
	double* swapr= new double[Buf_Size];

	double* swapl_send= new double[Buf_Size];
	double* swapr_send= new double[Buf_Size];


	


for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
		for(int l=0;l<19;l++)
		{
		swapl_send[j*(NZ+1)*19+k*19+l]=f[0][j][k][l];
		swapr_send[j*(NZ+1)*19+k*19+l]=f[(NX+1)/mpi_size-1][j][k][l];
		}
MPI_Barrier(MPI_COMM_WORLD);

if (rank==0)
		{
		MPI_Isend(swapr_send, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(swapl_send, (NY+1)*(NZ+1)*19, MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(swapr , (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(swapl, (NY+1)*(NZ+1)*19, MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );
		
		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(swapl_send, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(swapr_send, (NY+1)*(NZ+1)*19, MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(swapl, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(swapr, (NY+1)*(NZ+1)*19, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request[3]);
			
			}
			else
			{
			MPI_Isend(swapl_send, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(swapr_send, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(swapl, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(swapr, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);
			
			};

	
	MPI_Waitall(4,request, status);


int ip,jp,kp;
for(int i=0;i<(NX+1)/mpi_size;i++)	
	for(int j=0;j<=NY;j++)
		for(int k=0;k<=NZ;k++)
			//if (!Solid[i][j][k])
                        for(int lm=0;lm<19;lm++)


			{       

				ip=i-e[lm][0];
				jp=j-e[lm][1];if (jp<0) {jp=NY;}; if (jp>NY) {jp=0;};
				kp=k-e[lm][2];if (kp<0) {kp=NZ;}; if (kp>NZ) {kp=0;};
				if (ip<0)
					{
					F[i][j][k][lm]=swapl[jp*(NZ+1)*19+kp*19+lm];
					}
				if (ip>(NX+1)/mpi_size-1)
					{
					F[i][j][k][lm]=swapr[jp*(NZ+1)*19+kp*19+lm];
					}
				if ((ip>=0) and (ip<=(NX+1)/mpi_size-1))
					{
					F[i][j][k][lm]=f[ip][jp][kp][lm];
					}
					
			}
			



	for(int i=0;i<(NX+1)/mpi_size;i++)	
		for(int j=0;j<=NY;j++)
			for(int k=0;k<=NZ;k++)
				for(int lm=0;lm<19;lm++)
                		f[i][j][k][lm]=F[i][j][k][lm];

	delete [] swapl;
	delete [] swapr;
	delete [] swapl_send;
	delete [] swapr_send;


}




void collision(double*** rho,double**** u,double**** f,double w[19], int e[19][3],double M[19][19], double MI[19][19], double S[19],double*** forcex, double*** forcey, double*** forcez,int NX,int NY,int NZ,bool*** Solid)
{

double lm0,lm1;
double usqr,vsqr;
double F_hat[19],GuoF[19],f_eq[19],u_tmp[3];
double m_l[19];
double meq[19];


	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	for(int i=0;i<(NX+1)/mpi_size;i++)	
		for(int j=0;j<=NY;j++)
			for(int m=0;m<=NZ;m++)
			
			if (!Solid[i][j][m])     // LIQUID PHASE SEPERATION
                        //WARNNING:THE SOUND SPEED HERE IS sqrt(RT)=1/sqrt(3)

			{
			

			//=================FORCE TERM_GUO=========================================
			for (int k=0;k<19;k++)
			{	
			lm0=((e[k][0]-u[i][j][m][0])*forcex[i][j][m]+(e[k][1]-u[i][j][m][1])*forcey[i][j][m]+(e[k][2]-u[i][j][m][2])*forcez[i][j][m])*3;
			lm1=(e[k][0]*u[i][j][m][0]+e[k][1]*u[i][j][m][1]+e[k][2]*u[i][j][m][2])*(e[k][0]*forcex[i][j][m]+e[k][1]*forcey[i][j][m]+e[k][2]*forcez[i][j][m])*9;
			GuoF[k]=w[k]*(lm0+lm1);
			//GuoF[k]=0.0;
			}


			
			//=====================equilibrium of moment=================================
			u_tmp[0]=u[i][j][m][0];
			u_tmp[1]=u[i][j][m][1];
			u_tmp[2]=u[i][j][m][2];

			for(int k=0;k<19;k++)
				{
				f_eq[k]=feq(k,rho[i][j][m],u_tmp,w,e);
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
						m_l[mi]+=M[mi][mj]*f[i][j][m][mj];
						F_hat[mi]+=M[mi][mj]*GuoF[mj];
						}
					F_hat[mi]*=(1-0.5*S[mi]);
					}
			//============================================================================



				
			
				
				for (int sk=0;sk<19;sk++)
				//m[sk]=m[sk]-S[sk]*(m[sk]-meq[sk]);
				//m[sk]=m[sk]-S[sk]*(m[sk]-meq[sk])+(1-S[sk]*F_hat[sk]/2);
				m_l[sk]=m_l[sk]-S[sk]*(m_l[sk]-meq[sk])+F_hat[sk];

			//}

			// ==================   f=M_-1m matrix calculation  =============================
				for (int mi=0; mi<19; mi++)
					{f[i][j][m][mi]=0;
					for (int mj=0; mj<19; mj++)
						f[i][j][m][mi]+=MI[mi][mj]*m_l[mj];
					}
			//============================================================================
			
			
			}
                        else    
				standard_bounceback_boundary(i,j,m,f);
	
				   //BOUNCEBACK BOUNDARY CONDITION
			
		

}

void Read_Rock(bool*** Solid,int NX,int NY,int NZ,int mirX,int mirY, int mirZ,double* porosity,char poreFileName[128])
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
int i, j, k;

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

	
	

	if (Zoom<=1)
	{
	for (i=0;i<(NX+1)/mpi_size;i++)
		for (j=0;j<(NY+1);j++)
			for (k=0;k<(NZ+1);k++)
			if (Solid_Int[(rank*(NX+1)/mpi_size+i)*(NY+1)*(NZ+1)+j*(NZ+1)+k]==1)
				Solid[i][j][k]=1;
			else
				Solid[i][j][k]=0;
	}
	else
	{
	for (i=0;i<(nx)/mpi_size;i++)
		for (j=0;j<(ny);j++)
			for (k=0;k<(nz);k++)
			{			
			if (Solid_Int[(rank*(nx)/mpi_size+i)*(ny)*(nz)+j*(nz)+k]==1)
				{				
				for (int iz=0;iz<=Zoom-1;iz++)
					for (int jz=0;jz<=Zoom-1;jz++)
						for (int kz=0;kz<=Zoom-1;kz++)
						Solid[Zoom*i+iz][Zoom*j+jz][Zoom*k+kz]=1;	
				}			
				
			else
				{
				for (int iz=0;iz<=Zoom-1;iz++)
					for (int jz=0;jz<=Zoom-1;jz++)
						for (int kz=0;kz<=Zoom-1;kz++)
						Solid[Zoom*i+iz][Zoom*j+jz][Zoom*k+kz]=0;
				}
			}


	}
	

	delete [] Solid_Int;


}







void comput_macro_variables_IMR( double*** rho,double**** u,double**** u0,double**** f,int e[19][3],bool*** Solid, double*** forcex, double*** forcey, double*** forcez, int NX,int NY,int NZ,int n,double porosity,double* gx,double* gy, double* gz)
{
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	double rho0=1.0;
	double dp[3],rhok;
	double *rbuf;
	rbuf=new double[mpi_size*3];

	for(int i=0;i<(NX+1)/mpi_size;i++)	
		for(int j=0;j<=NY;j++)
			for(int m=0;m<=NZ;m++)
                   
			{ 

		

			if (!Solid[i][j][m]) // NOT INTERIOR SOLID NODE
			{
				
				rhok=rho[i][j][m];
				u0[i][j][m][0]=u[i][j][m][0];
				u0[i][j][m][1]=u[i][j][m][1];
				u0[i][j][m][2]=u[i][j][m][2];
				rho[i][j][m]=0;
				u[i][j][m][0]=0;
				u[i][j][m][1]=0;
				u[i][j][m][2]=0;
	
				for(int k=0;k<19;k++)
					{
					//f[i][j][k]=F[i][j][k];
					rho[i][j][m]+=f[i][j][m][k];
					u[i][j][m][0]+=e[k][0]*f[i][j][m][k];
					u[i][j][m][1]+=e[k][1]*f[i][j][m][k];
					u[i][j][m][2]+=e[k][2]*f[i][j][m][k];
					}
			

				
					

				u[i][j][m][0]=(u[i][j][m][0]+forcex[i][j][m]/2)/rho[i][j][m];
				u[i][j][m][1]=(u[i][j][m][1]+forcey[i][j][m]/2)/rho[i][j][m];
				u[i][j][m][2]=(u[i][j][m][2]+forcez[i][j][m]/2)/rho[i][j][m];

				dp[0]+=forcex[i][j][m]-(u[i][j][m][0]-u0[i][j][m][0])*rhok;
				dp[1]+=forcey[i][j][m]-(u[i][j][m][1]-u0[i][j][m][1])*rhok;
				dp[2]+=forcez[i][j][m]-(u[i][j][m][2]-u0[i][j][m][2])*rhok;
				
			}
			else
			{		u[i][j][m][0]=0.0; 
					u[i][j][m][1]=0.0;
					u[i][j][m][2]=0.0;
			}

                        
			}

	if ((n%20==0) and (n>100))
	{
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
	for(int i=0;i<(NX+1)/mpi_size;i++)	
		for(int j=0;j<=NY;j++)
			for(int m=0;m<=NZ;m++)
			{
		switch(PerDir)
		{
		case 1:
			forcex[i][j][m]=dp[0];break;
		case 2:
			forcey[i][j][m]=dp[1];break;
		case 3:
			forcez[i][j][m]=dp[2];break;
		default:
			forcex[i][j][m]=dp[0];
		}

			}
	*gx=dp[0];
	*gy=dp[1];
	*gz=dp[2];
	
	}
	delete [] rbuf;

}



void comput_macro_variables( double*** rho,double**** u,double**** u0,double**** f,int e[19][3],bool*** Solid, double*** forcex, double*** forcey, double*** forcez, int NX,int NY,int NZ)
{
	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	double rho0=1.0;

	for(int i=0;i<(NX+1)/mpi_size;i++)	
		for(int j=0;j<=NY;j++)
			for(int m=0;m<=NZ;m++)
                   
			{ 



			if (!Solid[i][j][m]) // NOT INTERIOR SOLID NODE
			{
				u0[i][j][m][0]=u[i][j][m][0];
				u0[i][j][m][1]=u[i][j][m][1];
				u0[i][j][m][2]=u[i][j][m][2];
				rho[i][j][m]=0;
				u[i][j][m][0]=0;
				u[i][j][m][1]=0;
				u[i][j][m][2]=0;
	
				for(int k=0;k<19;k++)
					{
					//f[i][j][k]=F[i][j][k];
					rho[i][j][m]+=f[i][j][m][k];
					u[i][j][m][0]+=e[k][0]*f[i][j][m][k];
					u[i][j][m][1]+=e[k][1]*f[i][j][m][k];
					u[i][j][m][2]+=e[k][2]*f[i][j][m][k];
					}
			

				
					

				u[i][j][m][0]=(u[i][j][m][0]+forcex[i][j][m]/2)/rho[i][j][m];
				u[i][j][m][1]=(u[i][j][m][1]+forcey[i][j][m]/2)/rho[i][j][m];
				u[i][j][m][2]=(u[i][j][m][2]+forcez[i][j][m]/2)/rho[i][j][m];
				
			}
			else
			{		u[i][j][m][0]=0.0; 
					u[i][j][m][1]=0.0;
					u[i][j][m][2]=0.0;
			}

                        
			}

}



void standard_bounceback_boundary(int it, int jt,int kt, double**** f)
{
	
	double tmp;
			tmp = f[it][jt][kt][1];f[it][jt][kt][1] = f[it][jt][kt][2];f[it][jt][kt][2] = tmp;
			tmp = f[it][jt][kt][3];f[it][jt][kt][3] = f[it][jt][kt][4];f[it][jt][kt][4] = tmp;
                        tmp = f[it][jt][kt][5];f[it][jt][kt][5] = f[it][jt][kt][6];f[it][jt][kt][6] = tmp;
			tmp = f[it][jt][kt][7];f[it][jt][kt][7] = f[it][jt][kt][10];f[it][jt][kt][10] = tmp;
			tmp = f[it][jt][kt][8];f[it][jt][kt][8] = f[it][jt][kt][9];f[it][jt][kt][9] = tmp;
			tmp = f[it][jt][kt][11];f[it][jt][kt][11] = f[it][jt][kt][14];f[it][jt][kt][14] = tmp;
                        tmp = f[it][jt][kt][12];f[it][jt][kt][12] = f[it][jt][kt][13];f[it][jt][kt][13] = tmp;
			tmp = f[it][jt][kt][15];f[it][jt][kt][15] = f[it][jt][kt][18];f[it][jt][kt][18] = tmp;
			tmp = f[it][jt][kt][16];f[it][jt][kt][16] = f[it][jt][kt][17];f[it][jt][kt][17] = tmp;


			

}




void boundary_velocity(int xp,double v_xp,int xn, double v_xn,int yp,double v_yp,int yn,double v_yn,int zp,double v_zp,int zn,double v_zn,double**** f,double*** rho,int NX, int NY,int NZ, double w[19], int e[19][3],bool*** Solid)

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
for (int i=0;i<=(NX+1)/mpi_size-1;i++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
			{
			if (!Solid[i][NY][k])
			f[i][NY][k][ks]=feq(ks,rho[i][NY-1][k],u_yp,w,e);
			}


if (yn==1)
for (int i=0;i<=(NX+1)/mpi_size-1;i++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
			{
			if (!Solid[i][0][k])
			f[i][0][k][ks]=feq(ks,rho[i][1][k],u_yn,w,e);
			}



if (zp==1)
for (int i=0;i<=(NX+1)/mpi_size-1;i++)
	for(int j=0;j<=NY;j++)
		for (int ks=0;ks<Q;ks++)
			{
			if (!Solid[i][j][NZ])
			f[i][j][NZ][ks]=feq(ks,rho[i][j][NZ-1],u_zp,w,e);
			}



if (zn==1)
for (int i=0;i<=(NX+1)/mpi_size-1;i++)
	for(int j=0;j<=NY;j++)
		for (int ks=0;ks<Q;ks++)
			{
			if (!Solid[i][j][0])
			f[i][j][0][ks]=feq(ks,rho[i][j][1],u_zn,w,e);
			}




if ((xp==1) && (rank==mpi_size-1))
for (int j=0;j<=NY;j++)
	for (int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
			{
			if (!Solid[(NX+1)/mpi_size-1][k][j])
			f[(NX+1)/mpi_size-1][j][k][ks]=feq(ks,rho[(NX+1)/mpi_size-2][j][k],u_xp,w,e);
			}



if ((xn==1) && (rank==0))
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
		for (int ks=0;ks<Q;ks++)
			{
			if (!Solid[0][j][k])
			f[0][j][k][ks]=feq(ks,rho[1][j][k],u_xn,w,e);
			}
	


}



void boundary_pressure(int xp,double rho_xp,int xn, double rho_xn,int yp,double rho_yp,int yn,double rho_yn,int zp,double rho_zp,int zn,double rho_zn,double**** f,double**** u,double*** rho,int NX, int NY,int NZ, double w[19], int e[19][3],bool*** Solid)

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
for (int i=0;i<=(NX+1)/mpi_size-1;i++)
	for(int k=0;k<=NZ;k++)
	{ 	
		u_ls[0]=u[i][NY-1][k][0];
		u_ls[1]=u[i][NY-1][k][1];
		u_ls[2]=u[i][NY-1][k][2];
		for (int ks=0;ks<Q;ks++)
			{
			if (!Solid[i][NY][k])
			f[i][NY][k][ks]=feq(ks,rho_yp,u_ls,w,e);
			}
	}


if (yn==1)
for (int i=0;i<=(NX+1)/mpi_size-1;i++)
	for(int k=0;k<=NZ;k++)
	{ 	
		u_ls[0]=u[i][1][k][0];
		u_ls[1]=u[i][1][k][1];
		u_ls[2]=u[i][1][k][2];
		for (int ks=0;ks<Q;ks++)
			{
			if (!Solid[i][0][k])
			f[i][0][k][ks]=feq(ks,rho_yn,u_ls,w,e);
			}
	}


if (zp==1)
for (int i=0;i<=(NX+1)/mpi_size-1;i++)
	for(int j=0;j<=NY;j++)
	{ 	
		u_ls[0]=u[i][j][NZ-1][0];
		u_ls[1]=u[i][j][NZ-1][1];
		u_ls[2]=u[i][j][NZ-1][2];
		for (int ks=0;ks<Q;ks++)
			{
			if (!Solid[i][j][NZ])
			f[i][j][NZ][ks]=feq(ks,rho_zp,u_ls,w,e);
			}
	}



if (zn==1)
for (int i=0;i<=(NX+1)/mpi_size-1;i++)
	for(int j=0;j<=NY;j++)
	{ 	
		u_ls[0]=u[i][j][1][0];
		u_ls[1]=u[i][j][1][1];
		u_ls[2]=u[i][j][1][2];
		for (int ks=0;ks<Q;ks++)
			{
			if (!Solid[i][j][0])
			f[i][j][0][ks]=feq(ks,rho_zn,u_ls,w,e);
			}
	}



if ((xp==1) && (rank==mpi_size-1))
      for (int j=0;j<=NY;j++)
		for (int k=0;k<=NZ;k++)
	{	
		u_ls[0]=u[(NX+1)/mpi_size-2][j][k][0];
		u_ls[1]=u[(NX+1)/mpi_size-2][j][k][1];
		u_ls[2]=u[(NX+1)/mpi_size-2][j][k][2];
		for (int ks=0;ks<Q;ks++)
			{
			if (!Solid[(NX+1)/mpi_size-1][j][k])
			f[(NX+1)/mpi_size-1][j][k][ks]=feq(ks,rho_xp,u_ls,w,e);
			//f[(NX+1)/mpi_size-1][j][k][ks]=feq(ks,rho[(NX+1)/mpi_size-2][j][k],u_ls,w,e);
			}

	}


if ((xn==1) && (rank==0))
      for (int j=0;j<=NY;j++)
		for(int k=0;k<=NZ;k++)
	{	
		u_ls[0]=u[1][j][k][0];
		u_ls[1]=u[1][j][k][1];
		u_ls[2]=u[1][j][k][1];
		for (int ks=0;ks<Q;ks++)
			{
			if (!Solid[0][j][k])
			f[0][j][k][ks]=feq(ks,rho_xn,u_ls,w,e);
			//f[0][j][k][ks]=feq(ks,rho[1][j][k],u_ls,w,e);
			}
	

	}

}




void output_velocity(int m,double**** u,int NX,int NY,int NZ)	
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	const int root_rank=0;
	int Gather_Count=((NX+1)/mpi_size)*(NY+1)*(NZ+1)*3;
	double* u_storage= new double[Gather_Count];
	double* rbuf_u;

	
	if (rank==root_rank)
		rbuf_u= new double[(NX+1)*(NY+1)*(NZ+1)*3];
	
	for (int i=0;i<(NX+1)/mpi_size;i++)
		for (int j=0;j<=NY;j++)
			for(int k=0;k<=NZ;k++)
			{
			u_storage[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3]=u[i][j][k][0];
			u_storage[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+1]=u[i][j][k][1];
			u_storage[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+2]=u[i][j][k][2];
			}

	


	MPI_Gather(u_storage,Gather_Count,MPI_DOUBLE,rbuf_u,Gather_Count,MPI_DOUBLE,root_rank,MPI_COMM_WORLD);



	if (rank==root_rank)
	{
	ostringstream name;
	name<<"LBM_velocity_x"<<m<<".vtk";
	ofstream out(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Velocity-X"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<NX+1<<"         "<<NY+1<<"         "<<NZ+1<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<(NX+1)*(NY+1)*(NZ+1)<<endl;
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	for(int k=0;k<=NZ;k++)
      		for(int j=0; j<=NY; j++)
			for(int i=0;i<=NX;i++)
        
			out<<"		"<<rbuf_u[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3]<<endl;
	
	ostringstream name2;
	name2<<"LBM_velocity_y"<<m<<".vtk";
	ofstream out2(name2.str().c_str());
	out2<<"# vtk DataFile Version 2.0"<<endl;
	out2<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Velocity-X"<<endl;
	out2<<"ASCII"<<endl;
	out2<<"DATASET STRUCTURED_POINTS"<<endl;
	out2<<"DIMENSIONS         "<<NX+1<<"         "<<NY+1<<"         "<<NZ+1<<endl;
	out2<<"ORIGIN 0 0 0"<<endl;
	out2<<"SPACING 1 1 1"<<endl;
	out2<<"POINT_DATA     "<<(NX+1)*(NY+1)*(NZ+1)<<endl;
	out2<<"SCALARS sample_scalars float"<<endl;
	out2<<"LOOKUP_TABLE default"<<endl;
	for(int k=0;k<=NZ;k++)
      		for(int j=0; j<=NY; j++)
			for(int i=0;i<=NX;i++)
        
			out2<<"		"<<rbuf_u[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+1]<<endl;

	
	ostringstream name3;
	name3<<"LBM_velocity_z"<<m<<".vtk";
	ofstream out3(name3.str().c_str());
	out3<<"# vtk DataFile Version 2.0"<<endl;
	out3<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Velocity-X"<<endl;
	out3<<"ASCII"<<endl;
	out3<<"DATASET STRUCTURED_POINTS"<<endl;
	out3<<"DIMENSIONS         "<<NX+1<<"         "<<NY+1<<"         "<<NZ+1<<endl;
	out3<<"ORIGIN 0 0 0"<<endl;
	out3<<"SPACING 1 1 1"<<endl;
	out3<<"POINT_DATA     "<<(NX+1)*(NY+1)*(NZ+1)<<endl;
	out3<<"SCALARS sample_scalars float"<<endl;
	out3<<"LOOKUP_TABLE default"<<endl;

	for(int k=0;k<=NZ;k++)
      		for(int j=0; j<=NY; j++)
			for(int i=0;i<=NX;i++)
			
			out3<<"		"<<rbuf_u[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+2]<<endl;
	


	}


	delete [] u_storage;
	if (rank==root_rank)
		delete [] rbuf_u;

		
}


void output_velocity_Vector(int m,double**** u,int NX,int NY,int NZ,int MirX,int MirY,int MirZ,int mir)	
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	const int root_rank=mpi_size-1;
	int Gather_Count=((NX+1)/mpi_size)*(NY+1)*(NZ+1)*3;
	double* u_storage= new double[Gather_Count];
	double* rbuf_u;

	
	if (rank==root_rank)
		rbuf_u= new double[(NX+1)*(NY+1)*(NZ+1)*3];
	
	for (int i=0;i<(NX+1)/mpi_size;i++)
		for (int j=0;j<=NY;j++)
			for(int k=0;k<=NZ;k++)
			{
			u_storage[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3]=u[i][j][k][0];
			u_storage[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+1]=u[i][j][k][1];
			u_storage[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+2]=u[i][j][k][2];
			}

	


	MPI_Gather(u_storage,Gather_Count,MPI_DOUBLE,rbuf_u,Gather_Count,MPI_DOUBLE,root_rank,MPI_COMM_WORLD);
	
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
        		out<<rbuf_u[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3]<<" "<<rbuf_u[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+1]<<" "<<rbuf_u[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+2]<<" "<<endl;
			//out<<endl;
			}

	out.close();
			


	}

	delete [] u_storage;
	if (rank==root_rank)
		delete [] rbuf_u;

}





void Geometry(bool*** Solid,int NX,int NY,int NZ,int MirX,int MirY,int MirZ,int mir)	
{	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	const int root_rank=0;
	int Gather_Count=((NX+1)/mpi_size)*(NY+1)*(NZ+1);

	int* Solid_storage= new int[Gather_Count];
	int* rbuf;

	for(int i=0;i<(NX+1)/mpi_size;i++)
		for(int j=0;j<=NY;j++)
			for(int k=0;k<=NZ;k++)
			if (Solid[i][j][k])
				Solid_storage[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=1;
			else
				Solid_storage[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=0;

	if (rank==root_rank)
		rbuf= new int[(NX+1)*(NY+1)*(NZ+1)];
	
	MPI_Gather(Solid_storage,Gather_Count,MPI_INT,rbuf,Gather_Count,MPI_INT,root_rank,MPI_COMM_WORLD);
	
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

		
}



void output_density(int m,double*** rho, int NX,int NY, int NZ,int MirX,int MirY,int MirZ,int mir)	
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	const int mpi_size=MPI :: COMM_WORLD . Get_size ();
	const int root_rank=mpi_size-1;
	int Gather_Count=((NX+1)/mpi_size)*(NY+1)*(NZ+1);

	double* rho_storage= new double[Gather_Count];
	double* rbuf_rho;

	double rho0=0.0;

	if (rank==root_rank)
		rbuf_rho= new double[(NX+1)*(NY+1)*(NZ+1)];
	
	for (int i=0;i<(NX+1)/mpi_size;i++)
		for (int j=0;j<=NY;j++)
			for(int k=0;k<=NZ;k++)
			{
			rho_storage[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=rho[i][j][k];
			}

	


	MPI_Gather(rho_storage,Gather_Count,MPI_DOUBLE,rbuf_rho,Gather_Count,MPI_DOUBLE,root_rank,MPI_COMM_WORLD);
	
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

	

	delete [] rho_storage;
	
	if (rank==root_rank)
		delete [] rbuf_rho;

		
}



double Comput_Perm(double**** u,double* Permia,int NX,int NY,int NZ,double s_v,double gx,double gy,double gz,int mirX,int mirY, int mirZ,double porosity,int PerDIr)
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	double *rbuf;
	rbuf=new double[mpi_size*3];
	double Perm[3];
	double error;
	double Q[3]={0.0,0.0,0.0};

for(int i=0; i<(NX+1)/mpi_size; i++)

	for(int j=0;j<=NY;j++)
		for(int k=0;k<=NZ;k++)
		{
		Q[0]+=u[i][j][k][0];
		Q[1]+=u[i][j][k][1];
		Q[2]+=u[i][j][k][2];
		}

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

		Perm[0]=Q[0]/((NX+1)*(NY+1)*(NZ+1))*(1.0/3.0*(1/s_v-0.5))/gx;
		Perm[1]=Q[1]/((NX+1)*(NY+1)*(NZ+1))*(1.0/3.0*(1/s_v-0.5))/gy;
		Perm[2]=Q[2]/((NX+1)*(NY+1)*(NZ+1))*(1.0/3.0*(1/s_v-0.5))/gz;
		
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

double Error(double**** u,double**** u0,int NX,int NY,int NZ,double *v_max,double* u_average)
{	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	double *rbuf,*um,*uave,u_compt;
	double temp1,temp2,temp3;
	temp1=0;
	temp2=0;
	temp3=0;
	double error;
	double u_max;
	rbuf=new double[mpi_size];
	um = new double[mpi_size];
	uave = new double[mpi_size];
	u_max=0;
	*v_max=0;

for(int i=1; i<(NX+1)/mpi_size; i++)
	for(int j=1;j<NY;j++)
		for(int k=1;k<NZ;k++)
		{	
			temp1 += (u[i][j][k][0]-u0[i][j][k][0])*(u[i][j][k][0]-u0[i][j][k][0])+(u[i][j][k][1]-u0[i][j][k][1])*(u[i][j][k][1]-u0[i][j][k][1])+(u[i][j][k][2]-u0[i][j][k][2])*(u[i][j][k][2]-u0[i][j][k][2]);
			temp2 += u[i][j][k][0]*u[i][j][k][0]+u[i][j][k][1]*u[i][j][k][1]+u[i][j][k][2]*u[i][j][k][2];	
			temp3+=sqrt(u[i][j][k][0]*u[i][j][k][0]+u[i][j][k][1]*u[i][j][k][1]+u[i][j][k][2]*u[i][j][k][2]);
			if (u[i][j][k][0]*u[i][j][k][0]+u[i][j][k][1]*u[i][j][k][1]+u[i][j][k][2]*u[i][j][k][2]>u_max)
				u_max=u[i][j][k][0]*u[i][j][k][0]+u[i][j][k][1]*u[i][j][k][1]+u[i][j][k][2]*u[i][j][k][2];
		}
		temp1=sqrt(temp1);
		temp2=sqrt(temp2);
		error=temp1/(temp2+1e-30);
		u_max=sqrt(u_max);
		
	
	MPI_Gather(&error,1,MPI_DOUBLE,rbuf,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
		
	MPI_Gather(&u_max,1,MPI_DOUBLE,um,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	MPI_Gather(&temp3,1,MPI_DOUBLE,uave,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

	u_compt=0;
	if (rank==0)
	    for (int i=0;i<mpi_size;i++)
		{
		if (rbuf[i]>error)
			error=rbuf[i];
		if (um[i]>*v_max)
			*v_max=um[i];
		u_compt+=uave[i];

		}

	u_compt/=(NX+1)*(NY+1)*(NZ+1);
	*u_average=u_compt;
	
	delete [] rbuf;
	delete [] um;
	delete [] uave;

	return(error);

}



