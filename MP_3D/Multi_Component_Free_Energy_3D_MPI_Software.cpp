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


void init(double***, double****, double****,double***, double****, double[19],int[19][3],bool***, double***, double***, double***, double[19], int, int, int,double,double,double,double,int,int,int,double,double);

void cylinder_creation(int, int, double);


void collision(double***, double****, double****,double****, double***, double****, double****, double[19],int[19][3],double[19][19], double[19][19], double[19],double***, double***, double***, int, int, int,bool***,double,double,int);   


void periodic_streaming(double***, double****, double****,double****, double***, double****, double****, double[19],int[19][3],double[19][19], double[19][19], double[19],double***, double***, double***, int, int, int,bool***,double,double,int);

void comput_macro_variables(double***, double****,double****,double****,double****, double***, int[19][3],bool***, double***, double***, double***,int,int, int,double);

void standard_bounceback_boundary(int,int,int, double****,double****);

void evolution(double***, double****,double****,double****);

double Error(double****,double****,int,int,int,double*,double*);

void output_velocity(int, double****,int,int,int);

void output_velocity_Vector(int, double****,int,int,int);

void output_density(int, double***,int,int,int);

void output_psi(int, double***,int,int,int,bool ***);

void Read_Rock(bool*** ,int,int,int,int,int,int,double*, char[128]);

void Read_psi(double***,int,int,int,int,int, int,char[128]);

double S[19];

void Comput_MI(double[19][19], double[19][19]);
int inverse(mat &a);

void Geometry(bool***, int,int,int );	

double Comput_Perm(double***,double****,double*,double*,int,int,int,double,double ,double ,int );

double Comput_Saturation(double*** ,bool*** ,int,int,int);

int mirX=0;
int mirY=0;
int mirZ=0;


double ca=0.04;
double kappa,CM;

int in_IMR,wr_per,pre_xp,pre_xn,pre_yp,pre_yn,pre_zp,pre_zn,Zoom;
int vel_xp,vel_xn,vel_yp,vel_yn,vel_zp,vel_zn,freRe,freVe,freDe,frePsi,PerDir,mir;
double in_vis,p_xp,p_xn,p_yp,p_yn,p_zp,p_zn,reso,porosity,Re_l,Re_g;
double inivx,inivy,inivz,v_xp,v_xn,v_yp,v_yn,v_zp,v_zn,S_l,S_g;
int vtk_time,n_max;
int Par_Geo,Par_nx,Par_ny,Par_nz;
double gx,gy,gz,niu_l,niu_g,ContactAngle_parameter,Permeability;

char filename[128],filenamepsi[128],outputfile[128];


int main (int argc , char * argv [])
{
MPI :: Init (argc , argv );
MPI_Status status ;

int rank = MPI :: COMM_WORLD . Get_rank ();
int para_size=MPI :: COMM_WORLD . Get_size ();



const int Q=19;          
int NX,NY,NZ;		



double s_v=1.0;
double v_max,u_ave,error_Per;
double Per_l[3],Per_g[3];

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

double m[19];
double meq[19];

int NCHAR=128;

char     dummy[NCHAR+1];
	int      dummyInt;
if (rank==0)
	{
	ifstream fin(argv[1]);
	
							fin.getline(dummy, NCHAR);
	fin >> filename;				fin.getline(dummy, NCHAR);
	fin >> filenamepsi;				fin.getline(dummy, NCHAR);
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
	fin >> niu_l;					fin.getline(dummy, NCHAR);
	fin >> niu_g;					fin.getline(dummy, NCHAR);
	fin >> ContactAngle_parameter;			fin.getline(dummy, NCHAR);
	fin >> kappa;					fin.getline(dummy, NCHAR);
	fin >> CM;					fin.getline(dummy, NCHAR);
	fin >> inivx >> inivy >> inivz;			fin.getline(dummy, NCHAR);
	fin >> Permeability;				fin.getline(dummy, NCHAR);
							fin.getline(dummy, NCHAR);
	fin >> wr_per;					fin.getline(dummy, NCHAR);
	fin >> PerDir;					fin.getline(dummy, NCHAR);
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
	

	NX=NX-1;NY=NY-1;NZ=NZ-1;
	}

	MPI_Bcast(&filename,128,MPI_CHAR,0,MPI_COMM_WORLD);MPI_Bcast(&filenamepsi,128,MPI_CHAR,0,MPI_COMM_WORLD);
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
	MPI_Bcast(&mir,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&niu_l,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&niu_g,1,MPI_DOUBLE,0,MPI_COMM_WORLD); MPI_Bcast(&ContactAngle_parameter,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&frePsi,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&kappa,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&CM,1,MPI_DOUBLE,0,MPI_COMM_WORLD);MPI_Bcast(&Permeability,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	MPI_Bcast(&Zoom,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&outputfile,128,MPI_CHAR,0,MPI_COMM_WORLD);




if (mirX)
	NX=NX*2-1;
if (mirY)
	NY=NY*2-1;
if (mirZ)
	NZ=NZ*2-1;


if (Zoom>1)
	{	
	NX=(NX+1)*Zoom-1;
	NY=(NY+1)*Zoom-1;
	NZ=(NZ+1)*Zoom-1;
	}


int n;
int e[19][3]=
{{0,0,0},{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};
double w[19]={1.0/3.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/18.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0,1.0/36.0};

	double error;

	double ****u;
	double ****u0;
	double ****f;
	double ***rho;
	double ****F;
	bool*** Solid;
	double*** forcex;
	double*** forcey;
	double*** forcez;
	double ****g;
	double ****Fg;
	double ***psi;


if ( (NX+1)%(MPI :: COMM_WORLD . Get_size ())!=0) exit (0);

  // Allocate memory
  	u = new double***[(NX+1)/para_size];
	u0 = new double***[(NX+1)/para_size];
	f = new double***[(NX+1)/para_size];
	rho = new double**[(NX+1)/para_size];
	F = new double***[(NX+1)/para_size];
	Solid = new bool**[(NX+1)/para_size];
	forcex = new double**[(NX+1)/para_size];
	forcey = new double**[(NX+1)/para_size];
	forcez = new double**[(NX+1)/para_size];
	psi = new double**[(NX+1)/para_size];
	g = new double***[(NX+1)/para_size];
	Fg =  new double***[(NX+1)/para_size];


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
	psi[i]=new double*[NY+1];
	g[i]=new double**[NY+1];
	Fg[i]=new double**[NY+1];


	
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
	psi[i][j]=new double[NZ+1];
	g[i][j]=new double*[NZ+1];
	Fg[i][j]=new double*[NZ+1];

	for (int k=0;k<NZ+1;k++)
		{
		u[i][j][k] = new double[3];
		u0[i][j][k] = new double[3];
		f[i][j][k] = new double[19];
		F[i][j][k] = new double[19];
		g[i][j][k] = new double[19];
		Fg[i][j][k]= new double[19];
		}
	}
  }


Comput_MI(M,MI);
n=0;


	using namespace std;
	init(rho,u,f,psi,g,w,e,Solid,forcex,forcey,forcez,S,NX,NY,NZ,s_v,gx,gy,gz,mirX,mirY,mirZ,niu_l,niu_g);Geometry(Solid,NX,NY,NZ);


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

	collision(rho,u,f,g,psi,F,Fg,w,e,M,MI,S,forcex,forcey,forcez,NX,NY,NZ,Solid,niu_l,niu_g,n);
	
	periodic_streaming(rho,u,f,g,psi,F,Fg,w,e,M,MI,S,forcex,forcey,forcez,NX,NY,NZ,Solid,niu_l,niu_g,n);//cout<<"mark"<<endl;


  	comput_macro_variables(rho,u,u0,f,g,psi,e,Solid,forcex,forcey,forcez,NX,NY,NZ,ContactAngle_parameter); 

		if(n%freRe==0)
		{       
			error=Error(u,u0,NX,NY,NZ,&v_max,&u_ave);
			error_Per=Comput_Perm(psi,u,Per_l,Per_g,NX,NY,NZ,gx,gy,gz,PerDir);
			S_l=Comput_Saturation(psi,Solid,NX,NY,NZ);
			
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
			fin<<"The relative permeability of component 1 is "<<Per_l[0]*reso*reso*100/Permeability<<", "<<Per_l[1]*reso*reso*100/Permeability<<", "<<Per_l[2]*reso*reso*100/Permeability<<endl;
			fin<<"The relative permeability of component 2 is "<<Per_g[0]*reso*reso*100/Permeability<<", "<<Per_g[0]*reso*reso*100/Permeability<<", "<<Per_g[2]*reso*reso*100/Permeability<<endl;
			fin<<"Satuation of Component 1: "<<S_l<<", "<<"The satuation of Component 2: "<<1-S_l<<endl;
			fin<<"The relative error of permiability computing is: "<<error_Per<<endl;
			fin<<endl;
			fin.close();

/*
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
*/

//============================================================================================================

			cout<<"The"<<n<<"th computation result:"<<endl;
			cout<<"The Density of point(NX/2,NY/2,NZ/2) is: "<<setprecision(6)
				<<rho[int((NX+1)/para_size/2)][NY/2][NZ/2]<<endl;
			
			cout<<"The Maximum velocity is: "<<setprecision(6)<<v_max<<"   Re_l="<<Re_l<<"   Re_g="<<Re_g<<endl;
			
			cout<<"The max relative error of uv is: "
				<<setiosflags(ios::scientific)<<error<<endl;
			cout<<"The relative permeability of component 1 is "<<Per_l[0]*reso*reso*100/Permeability<<", "<<Per_l[1]*reso*reso*100/Permeability<<", "<<Per_l[2]*reso*reso*100/Permeability<<endl;
			cout<<"The relative permeability of component 2 is "<<Per_g[0]*reso*reso*100/Permeability<<", "<<Per_g[0]*reso*reso*100/Permeability<<", "<<Per_g[2]*reso*reso*100/Permeability<<endl;
			cout<<"Satuation of Component 1: "<<S_l<<", "<<"The satuation of Component 2: "<<1-S_l<<endl;
			cout<<"The relative error of permiability computing is: "<<error_Per<<endl;
			cout<<endl;
                        }

			if(n>=0)
			{
				if ((freDe>=0) and (n%freDe==0))
						output_density(n,rho,NX,NY,NZ);
				if ((frePsi>=0) and (n%frePsi==0))
						output_psi(n,psi,NX,NY,NZ,Solid);
				if ((freVe>=0) and (n%freVe==0))
						output_velocity_Vector(n,u,NX,NY,NZ);
					
				if(error!=error) {cout<<"PROGRAM STOP"<<endl;break;};
			}

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
			delete [] Fg[i][j][k];
			delete [] g[i][j][k];
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
		delete [] Fg[i][j];
		delete [] g[i][j];
		delete [] psi[i][j];


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
	delete [] Fg[i];
	delete [] g[i];
	delete [] psi[i];




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
	delete [] Fg;
	delete [] g;
	delete [] psi;

	
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




void init(double*** rho, double**** u, double**** f,double*** psi, double**** g,double w[19], int e[19][3],bool*** Solid, double*** forcex, double*** forcey, double*** forcez, double S[19], int NX,int NY, int NZ, double s_v, double gx,double gy,double gz,int mirX,int mirY,int mirZ,double niu_l,double niu_g)
{	
        //=================
        //double res[9][9];
        //=================
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	double Cylinder_r=7;
	double usqr,vsqr;
	double dx=1.0;
	double dy=1.0;
	double Lx=dx*double(NX);
	double Ly=dy*double(NY);
	double dt=1;//sqrt(3);
	double c=dx/dt;
	double rho0=1.0;
 	double uMax=0.00;
	double Re=100.0;
	//forcex=gx;forcey=gy;
	//niu=U*Lx/Re;
	double niu=uMax*Cylinder_r*2/Re;
	double tau_f=3.0*niu+0.5;
	
	double pr; //raduis of the obstacles
        

	

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


	if (!(filename=="NOSOLID"))
		
		Read_Rock(Solid,NX,NY,NZ,mirX,mirY,mirZ,&porosity,filename);
		
	else
		{
		for(int i=0;i<(NX+1)/mpi_size;i++)	
			for(int j=0;j<=NY;j++)
				for(int k=0;k<=NZ;k++)
				Solid[i][j][k]=0;
		}

	Read_psi(psi,NX,NY,NZ,mirX,mirY,mirZ,filenamepsi);

	
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
			
			


			s_v=niu_g+(psi[i][j][k]+1.0)/2.0*(niu_l-niu_g);
			s_v=1.0/(3*s_v+0.5);
			
			S[9]=s_v;S[11]=s_v;S[13]=s_v;S[14]=s_v;S[15]=s_v;




			//***********************************************************************


		

		
		

	}

		
		


	 	
}





void periodic_streaming(double*** rho,double**** u,double**** f,double**** g, double*** psi,double**** F, double**** Fg,double w[19], int e[19][3],double M[19][19], double MI[19][19], double S[19],double*** forcex, double*** forcey, double*** forcez,int NX,int NY,int NZ,bool*** Solid,double niu_l,double niu_g,int n)
{

	MPI_Status status[4] ;
	MPI_Request request[4];

	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	int Buf_Size=(NY+1)*(NZ+1)*19;
	
	double* swapl= new double[Buf_Size];
	double* swapr= new double[Buf_Size];
	double* swapl_g= new double[Buf_Size];
	double* swapr_g= new double[Buf_Size];
	


	
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
		for(int l=0;l<19;l++)
		{
		swapl[j*(NZ+1)*19+k*19+l]=f[0][j][k][l];
		swapr[j*(NZ+1)*19+k*19+l]=f[(NX+1)/mpi_size-1][j][k][l];
		swapl_g[j*(NZ+1)*19+k*19+l]=g[0][j][k][l];
		swapr_g[j*(NZ+1)*19+k*19+l]=g[(NX+1)/mpi_size-1][j][k][l];
		}
	


if (rank==0)
		{
		MPI_Isend(swapr, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(swapl, (NY+1)*(NZ+1)*19, MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(swapr , (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(swapl, (NY+1)*(NZ+1)*19, MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );

		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(swapl, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(swapr, (NY+1)*(NZ+1)*19, MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(swapl, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(swapr, (NY+1)*(NZ+1)*19, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request[3]);
			
			}
			else
			{
			MPI_Isend(swapl, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(swapr, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(swapl, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(swapr, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);
			
			};

	
	MPI_Waitall(4,request, status);



if (rank==0)
		{
		MPI_Isend(swapr_g, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(swapl_g, (NY+1)*(NZ+1)*19, MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(swapr_g, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(swapl_g, (NY+1)*(NZ+1)*19, MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );

		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(swapl_g, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(swapr_g, (NY+1)*(NZ+1)*19, MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(swapl_g, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(swapr_g, (NY+1)*(NZ+1)*19, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request[3]);
			
			}
			else
			{
			MPI_Isend(swapl_g, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(swapr_g, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(swapl_g, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(swapr_g, (NY+1)*(NZ+1)*19, MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);
			
			};

	
	MPI_Waitall(4,request, status);



int ip,jp,kp;
for(int i=0;i<(NX+1)/mpi_size;i++)	
	for(int j=0;j<=NY;j++)
		for(int k=0;k<=NZ;k++)
			//if (!interior_node(i,j))
                        for(int lm=0;lm<19;lm++)


			{       

				ip=i-e[lm][0];
				jp=j-e[lm][1];if (jp<0) {jp=NY;}; if (jp>NY) {jp=0;};
				kp=k-e[lm][2];if (kp<0) {kp=NZ;}; if (kp>NZ) {kp=0;};
				if (ip<0)
					{
					F[i][j][k][lm]=swapl[jp*(NZ+1)*19+kp*19+lm];
					Fg[i][j][k][lm]=swapl_g[jp*(NZ+1)*19+kp*19+lm];
					}
				if (ip>(NX+1)/mpi_size-1)
					{
					F[i][j][k][lm]=swapr[jp*(NZ+1)*19+kp*19+lm];
					Fg[i][j][k][lm]=swapr_g[jp*(NZ+1)*19+kp*19+lm];
					}
				if ((ip>=0) and (ip<=(NX+1)/mpi_size-1))
					{
					F[i][j][k][lm]=f[ip][jp][kp][lm];
					Fg[i][j][k][lm]=g[ip][jp][kp][lm];
					}
					
			}
			



	for(int i=0;i<(NX+1)/mpi_size;i++)	
		for(int j=0;j<=NY;j++)
			for(int k=0;k<=NZ;k++)
				for(int lm=0;lm<19;lm++)
				{                		
				f[i][j][k][lm]=F[i][j][k][lm];
				g[i][j][k][lm]=Fg[i][j][k][lm];
				}

	delete [] swapl;
	delete [] swapr;
	delete [] swapl_g;
	delete [] swapr_g;



}


void collision(double*** rho,double**** u,double**** f,double**** g, double*** psi,double**** F, double**** Fg,double w[19], int e[19][3],double M[19][19], double MI[19][19], double S[19],double*** forcex, double*** forcey, double*** forcez,int NX,int NY,int NZ,bool*** Solid,double niu_l,double niu_g,int n)
{

	MPI_Status status[4] ;
	MPI_Request request[4];

	
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	
	
	
	double* swapl_r = new double[(NY+1)*(NZ+1)];
	double* swapr_r = new double[(NY+1)*(NZ+1)];
	double* swapl_psi = new double[(NY+1)*(NZ+1)];
	double* swapr_psi = new double[(NY+1)*(NZ+1)];


double lm0,lm1,s_v;
double A,B,C,D,E;
double F_hat[19],GuoF[19],f_eq[19];
double m_l[19],meq[19];;
double feq[19],feq_g[19],Fi[19];

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

	
for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
		{
		swapl_r[j*(NZ+1)+k]=rho[0][j][k];
		swapr_r[j*(NZ+1)+k]=rho[(NX+1)/mpi_size-1][j][k];
		swapl_psi[j*(NZ+1)+k]=psi[0][j][k];
		swapr_psi[j*(NZ+1)+k]=psi[(NX+1)/mpi_size-1][j][k];
		}


if (rank==0)
		{
		MPI_Isend(swapr_r, (NY+1)*(NZ+1), MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(swapl_r, (NY+1)*(NZ+1), MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(swapr_r , (NY+1)*(NZ+1), MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(swapl_r, (NY+1)*(NZ+1), MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );

		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(swapl_r, (NY+1)*(NZ+1), MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(swapr_r, (NY+1)*(NZ+1), MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(swapl_r, (NY+1)*(NZ+1), MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(swapr_r, (NY+1)*(NZ+1), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request[3]);
			
			}
			else
			{
			MPI_Isend(swapl_r, (NY+1)*(NZ+1), MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(swapr_r, (NY+1)*(NZ+1), MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(swapl_r, (NY+1)*(NZ+1), MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(swapr_r, (NY+1)*(NZ+1), MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);
			
			};

	
	MPI_Waitall(4,request, status);

if (rank==0)
		{
		MPI_Isend(swapr_psi, (NY+1)*(NZ+1), MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(swapl_psi, (NY+1)*(NZ+1), MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(swapr_psi , (NY+1)*(NZ+1), MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(swapl_psi, (NY+1)*(NZ+1), MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );

		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(swapl_psi, (NY+1)*(NZ+1), MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(swapr_psi, (NY+1)*(NZ+1), MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(swapl_psi, (NY+1)*(NZ+1), MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(swapr_psi, (NY+1)*(NZ+1), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request[3]);
			
			}
			else
			{
			MPI_Isend(swapl_psi, (NY+1)*(NZ+1), MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(swapr_psi, (NY+1)*(NZ+1), MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(swapl_psi, (NY+1)*(NZ+1), MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(swapr_psi, (NY+1)*(NZ+1), MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);
			
			};

	
	MPI_Waitall(4,request, status);



double term1,term2,term3,tempvar;
double lambda,miu,niu_temp;  	//NEED TO BE DETERMINED  *******
int interi,interj,interk;

	for(int i=0;i<(NX+1)/mpi_size;i++)	
		for(int j=0;j<=NY;j++)
			for(int m=0;m<=NZ;m++)
			
			if ((!Solid[i][j][m]) or (n==0))     // LIQUID PHASE SEPERATION
                        //WARNNING:THE SOUND SPEED HERE IS sqrt(RT)=1/sqrt(3)

			{
			


//=====================================BULK PRESSURE DEFINITION=====================================================

//==============================NUMERICAL DERIVATIVES COMPUTING=(BINARY SYSTEM)============================
		par[1]=0;par[0]=0;par[2]=0;delxy=0;par_r[1]=0;par_r[0]=0;par_r[2]=0;delxy_r=0;
		for (int tmpi=0;tmpi<19;tmpi++)
			{

			//-------------------PERIODIC BOUNDARY CONDITION---------------------------
			interi=i+e[tmpi][0];
			interj=j+e[tmpi][1];if (interj<0) {interj=NY;}; if (interj>NY) {interj=0;};
			interk=m+e[tmpi][2];if (interk<0) {interk=NZ;}; if (interk>NZ) {interk=0;};
			//-------------------------------------------------------------------------
			if (interi<0)
			{
			par_r[0]+=px[tmpi]*swapl_r[interj*(NZ+1)+interk];
			par_r[1]+=py[tmpi]*swapl_r[interj*(NZ+1)+interk];
			par_r[2]+=pz[tmpi]*swapl_r[interj*(NZ+1)+interk];

			par[0]+=px[tmpi]*swapl_psi[interj*(NZ+1)+interk];
			par[1]+=py[tmpi]*swapl_psi[interj*(NZ+1)+interk];
			par[2]+=pz[tmpi]*swapl_psi[interj*(NZ+1)+interk];

			delxy+=lap[tmpi]*swapl_psi[interj*(NZ+1)+interk];
			delxy_r+=lap[tmpi]*swapl_r[interj*(NZ+1)+interk];
			}

			if (interi>(NX+1)/mpi_size-1)
			{
			par_r[0]+=px[tmpi]*swapr_r[interj*(NZ+1)+interk];
			par_r[1]+=py[tmpi]*swapr_r[interj*(NZ+1)+interk];
			par_r[2]+=pz[tmpi]*swapr_r[interj*(NZ+1)+interk];

			par[0]+=px[tmpi]*swapr_psi[interj*(NZ+1)+interk];
			par[1]+=py[tmpi]*swapr_psi[interj*(NZ+1)+interk];
			par[2]+=pz[tmpi]*swapr_psi[interj*(NZ+1)+interk];

			delxy+=lap[tmpi]*swapr_psi[interj*(NZ+1)+interk];
			delxy_r+=lap[tmpi]*swapr_r[interj*(NZ+1)+interk];
			}

			if ((interi>=0) and (interi<=(NX+1)/mpi_size-1))
			{
			par_r[0]+=px[tmpi]*rho[interi][interj][interk];
			par_r[1]+=py[tmpi]*rho[interi][interj][interk];
			par_r[2]+=pz[tmpi]*rho[interi][interj][interk];

			par[0]+=px[tmpi]*psi[interi][interj][interk];
			par[1]+=py[tmpi]*psi[interi][interj][interk];
			par[2]+=pz[tmpi]*psi[interi][interj][interk];

			delxy+=lap[tmpi]*psi[interi][interj][interk];
			delxy_r+=lap[tmpi]*rho[interi][interj][interk];
			}
			
			
			}

	tempvar=psi[i][j][m]*psi[i][j][m];
	p0=rho[i][j][m]/3.0+ca*(-0.5*tempvar+0.75*tempvar*tempvar);
	lambda=0.0;
	miu=ca*(-psi[i][j][m]+tempvar*psi[i][j][m])-kappa*delxy;

//===========================================================================================================

//==============================EQUILIBRIUM FUNCTION COMPUTING================================================




	for (int k=1;k<19;k++)
		{
		feq[k]=0;feq_g[k]=0;
		
		feq[k]+=p0-kappa*psi[i][j][m]*delxy+(e[k][0]*u[i][j][m][0]+e[k][1]*u[i][j][m][1]+e[k][2]*u[i][j][m][2])*rho[i][j][m];
		
		term1=0;term2=0;term3=0;
		
		for(int ii=0;ii<=2;ii++)
			for (int jj=0;jj<=2;jj++)
			{
		term1+=(e[k][ii]*e[k][jj]-(1.0/3.0)*(delta[ii][jj]))*(rho[i][j][m]*u[i][j][m][ii]*u[i][j][m][jj]+lambda*(u[i][j][m][ii]*par_r[jj]+u[i][j][m][jj]*par_r[ii]+delta[ii][jj]*(u[i][j][m][0]*par_r[0]+u[i][j][m][1]*par_r[1]+u[i][j][m][2]*par_r[2])));

			term2+=1.5*(e[k][ii]*e[k][jj]-1.0/3.0*delta[ii][jj])*psi[i][j][m]*u[i][j][m][ii]*u[i][j][m][jj];
		
			}
		feq_g[k]=w[k]*3.0*(2.0*CM*miu+e[k][0]*psi[i][j][m]*u[i][j][m][0]+e[k][1]*psi[i][j][m]*u[i][j][m][1]+e[k][2]*psi[i][j][m]*u[i][j][m][2]+term2);

		term1*=3.0/2.0;
		term2=kappa*(ome[k][0]*par[0]*par[0]+ome[k][1]*par[1]*par[1]+ome[k][2]*par[2]*par[2]+ome[k][3]*par[0]*par[1]+ome[k][4]*par[1]*par[2]+ome[k][5]*par[2]*par[0]);

		
		feq[k]=(term1+feq[k])*w[k]*3.0+term2;	
		
		}

	feq[0]=rho[i][j][m]-(feq[1]+feq[2]+feq[3]+feq[4]+feq[5]+feq[6]+feq[7]+feq[8]+feq[9]+feq[10]+feq[11]+feq[12]+feq[13]+feq[14]+feq[15]+feq[16]+feq[17]+feq[18]);

	feq_g[0]=psi[i][j][m]-(feq_g[1]+feq_g[2]+feq_g[3]+feq_g[4]+feq_g[5]+feq_g[6]+feq_g[7]+feq_g[8]+feq_g[9]+feq_g[10]+feq_g[11]+feq_g[12]+feq_g[13]+feq_g[14]+feq_g[15]+feq_g[16]+feq_g[17]+feq_g[18]);

	

	//=========================================================================================================

//==========================================EQILIBRIUM FUNCTION COMPUTING FOR G FUNCTION===========================
//***THIS PART IS FOR G FUNCTION RELAXATION TIME =1 , FOR OTHER RELAXATION TIME,MODIFICATION SHOULD BE INVOLVED****





			//=================FORCE TERM_GUO=========================================
			for (int k=0;k<19;k++)
			{	
			lm0=((e[k][0]-u[i][j][m][0])*forcex[i][j][m]+(e[k][1]-u[i][j][m][1])*forcey[i][j][m]+(e[k][2]-u[i][j][m][2])*forcez[i][j][m])*3;
			lm1=(e[k][0]*u[i][j][m][0]+e[k][1]*u[i][j][m][1]+e[k][2]*u[i][j][m][2])*(e[k][0]*forcex[i][j][m]+e[k][1]*forcey[i][j][m]+e[k][2]*forcez[i][j][m])*9;
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

			s_v=niu_g+(psi[i][j][m]+1.0)/2.0*(niu_l-niu_g);
			s_v=1.0/(3*s_v+0.5);

			S[9]=s_v;S[11]=s_v;S[13]=s_v;S[14]=s_v;S[15]=s_v;

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


			g[i][j][m][0]=feq_g[0];g[i][j][m][1]=feq_g[1];g[i][j][m][2]=feq_g[2];g[i][j][m][3]=feq_g[3];
			g[i][j][m][4]=feq_g[4];g[i][j][m][5]=feq_g[5];g[i][j][m][6]=feq_g[6];g[i][j][m][7]=feq_g[7];
			g[i][j][m][8]=feq_g[8];g[i][j][m][9]=feq_g[9];g[i][j][m][10]=feq_g[10];g[i][j][m][11]=feq_g[11];
			g[i][j][m][12]=feq_g[12];g[i][j][m][13]=feq_g[13];g[i][j][m][14]=feq_g[14];g[i][j][m][15]=feq_g[15];
			g[i][j][m][16]=feq_g[16];g[i][j][m][17]=feq_g[17];g[i][j][m][18]=feq_g[18];

			if (n==0)
			{
			f[i][j][m][0]=feq[0];f[i][j][m][1]=feq[1];f[i][j][m][2]=feq[2];f[i][j][m][3]=feq[3];
			f[i][j][m][4]=feq[4];f[i][j][m][5]=feq[5];f[i][j][m][6]=feq[6];f[i][j][m][7]=feq[7];
			f[i][j][m][8]=feq[8];f[i][j][m][9]=feq[9];f[i][j][m][10]=feq[10];f[i][j][m][11]=feq[11];
			f[i][j][m][12]=feq[12];f[i][j][m][13]=feq[13];f[i][j][m][14]=feq[14];f[i][j][m][15]=feq[15];
			f[i][j][m][16]=feq[16];f[i][j][m][17]=feq[17];f[i][j][m][18]=feq[18];
			}



			}
                        else    
				
                        		{
					standard_bounceback_boundary(i,j,m,f,g);
					};    //BOUNCEBACK BOUNDARY CONDITION









	delete [] swapl_r;
	delete [] swapr_r;
	delete [] swapl_psi;
	delete [] swapr_psi;







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
	//double porosity;

	// Calculate Porosity
	int nNodes = 0;
	for(i=0 ; i<nx*ny*nz ; i++)
		if(Solid_Int[i] == 0) nNodes++;

	*porosity = (double)nNodes / (nx*ny*nz);

}


	MPI_Bcast(porosity,1,MPI_DOUBLE,0,MPI_COMM_WORLD);

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



void Read_psi(double*** psi,int NX,int NY,int NZ,int mirX,int mirY, int mirZ,char poreFileName[128])
{


int rank = MPI :: COMM_WORLD . Get_rank ();
int mpi_size=MPI :: COMM_WORLD . Get_size ();



int nx0=NX+1;
int ny0=NY+1;
int nz0=NZ+1;


int nx=NX+1;
int ny=NY+1;
int nz=NZ+1;

double* psi_Int;
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

psi_Int = new double[nx*ny*nz];

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
			fin >> pore ;
			if( pore == -1.0 || pore == 1.0) break;
			 
		}
		psi_Int[i*ny*nz+j*nz+k] = pore; 
		
	}
	fin.close();

	// Mirroring the rock
	if(mirX==1){
		for(i=nx0 ; i<nx ; i++)
		for(j=0   ; j<ny ; j++)
		for(k=0   ; k<nz ; k++)
				psi_Int[i*ny*nz+j*nz+k] = psi_Int[(nx-i-1)*ny*nz+j*nz+k];
	}

	if(mirY==1){
		for(j=ny0 ; j<ny ; j++)
		for(i=0   ; i<nx ; i++)
		for(k=0   ; k<nz ; k++)
				psi_Int[i*ny*nz+j*nz+k] = psi_Int[i*ny*nz+(ny-j-1)*nz+k];
	}

	if(mirZ==1){
		for(k=nz0 ; k<nz ; k++)
		for(j=0   ; j<ny ; j++)
		for(i=0   ; i<nx ; i++)
				psi_Int[i*ny*nz+j*nz+k] = psi_Int[i*ny*nz+j*nz+nz-k-1];
	}
	//double porosity;

	
}


	
	MPI_Bcast(psi_Int,nx*ny*nz,MPI_DOUBLE,0,MPI_COMM_WORLD);

	if (Zoom<=1)
	{
	for (i=0;i<(NX+1)/mpi_size;i++)
		for (j=0;j<(NY+1);j++)
			for (k=0;k<(NZ+1);k++)
			psi[i][j][k]=psi_Int[(rank*(NX+1)/mpi_size+i)*(NY+1)*(NZ+1)+j*(NZ+1)+k];

	}
	else
	{
	for (i=0;i<(nx)/mpi_size;i++)
		for (j=0;j<(ny);j++)
			for (k=0;k<(nz);k++)			
				for (int iz=0;iz<=Zoom-1;iz++)
					for (int jz=0;jz<=Zoom-1;jz++)
						for (int kz=0;kz<=Zoom-1;kz++)
						psi[Zoom*i+iz][Zoom*j+jz][Zoom*k+kz]=psi_Int[(rank*(nx)/mpi_size+i)*(ny)*(nz)+j*(nz)+k];	
							
				
	
			


	}

	

	delete [] psi_Int;


}



void comput_macro_variables( double*** rho,double**** u,double**** u0,double**** f,double**** g,double*** psi,int e[19][3],bool*** Solid, double*** forcex, double*** forcey, double*** forcez, int NX,int NY,int NZ,double ContactAngle_parameter)
{
	MPI_Status status[4] ;
	MPI_Request request[4];
	int interi,interj,interk;
	bool markmac;
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();
	double rho0=1.0;


	int* swapl= new int[(NY+1)*(NZ+1)];
	int* swapr= new int[(NY+1)*(NZ+1)];
	double* swapl_psi = new double[(NY+1)*(NZ+1)];
	double* swapr_psi = new double[(NY+1)*(NZ+1)];

for (int j=0;j<=NY;j++)
	for(int k=0;k<=NZ;k++)
		{
		if (Solid[0][j][k]==1)
		swapl[j*(NZ+1)+k]=1;
		else
		swapl[j*(NZ+1)+k]=0;
		swapl_psi[j*(NZ+1)+k]=psi[0][j][k];

		if (Solid[(NX+1)/mpi_size-1][j][k]==1)
		swapr[j*(NZ+1)+k]=1;
		else
		swapl[j*(NZ+1)+k]=0;
		swapr_psi[j*(NZ+1)+k]=psi[(NX+1)/mpi_size-1][j][k];
		}


if (rank==0)
		{
		MPI_Isend(swapr, (NY+1)*(NZ+1), MPI_INT, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(swapl, (NY+1)*(NZ+1), MPI_INT, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(swapr , (NY+1)*(NZ+1), MPI_INT, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(swapl, (NY+1)*(NZ+1), MPI_INT, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );
		
		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(swapl, (NY+1)*(NZ+1), MPI_INT, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(swapr, (NY+1)*(NZ+1), MPI_INT, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(swapl, (NY+1)*(NZ+1), MPI_INT, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(swapr, (NY+1)*(NZ+1), MPI_INT, 0, 0, MPI_COMM_WORLD,&request[3]);
			
			}
			else
			{
			MPI_Isend(swapl, (NY+1)*(NZ+1), MPI_INT, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(swapr, (NY+1)*(NZ+1), MPI_INT, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(swapl, (NY+1)*(NZ+1), MPI_INT, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);
      			MPI_Irecv(swapr, (NY+1)*(NZ+1), MPI_INT, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);
			
			};

	
	MPI_Waitall(4,request, status);


if (rank==0)
		{
		MPI_Isend(swapr_psi, (NY+1)*(NZ+1), MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[0]);
      		MPI_Isend(swapl_psi, (NY+1)*(NZ+1), MPI_DOUBLE, mpi_size-1, rank*2, MPI_COMM_WORLD,&request[1]);
		MPI_Irecv(swapr_psi , (NY+1)*(NZ+1), MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[2]);		
      		MPI_Irecv(swapl_psi, (NY+1)*(NZ+1), MPI_DOUBLE, mpi_size-1, (mpi_size-1)*2+1, MPI_COMM_WORLD,&request[3] );
		
		}
		else
		if (rank==mpi_size-1)
			{
			MPI_Isend(swapl_psi, (NY+1)*(NZ+1), MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(swapr_psi, (NY+1)*(NZ+1), MPI_DOUBLE, 0, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(swapl_psi, (NY+1)*(NZ+1), MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2] );
      			MPI_Irecv(swapr_psi, (NY+1)*(NZ+1), MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,&request[3]);
			
			}
			else
			{
			MPI_Isend(swapl_psi, (NY+1)*(NZ+1), MPI_DOUBLE, rank-1, rank*2, MPI_COMM_WORLD,&request[0]);
      			MPI_Isend(swapr_psi, (NY+1)*(NZ+1), MPI_DOUBLE, rank+1, rank*2+1, MPI_COMM_WORLD,&request[1]);
			MPI_Irecv(swapl_psi, (NY+1)*(NZ+1), MPI_DOUBLE, rank-1, (rank-1)*2+1, MPI_COMM_WORLD,&request[2]);

      			MPI_Irecv(swapr_psi, (NY+1)*(NZ+1), MPI_DOUBLE, rank+1, (rank+1)*2, MPI_COMM_WORLD,&request[3]);
			
			};

	
	MPI_Waitall(4,request, status);

	for(int i=0;i<(NX+1)/mpi_size;i++)	
		for(int j=0;j<=NY;j++)
			for(int m=0;m<=NZ;m++)
                   
			{ 



			if (!Solid[i][j][m]) // NOT INTERIOR SOLID NODE
			{
				u0[i][j][m][0]=u[i][j][m][0];
				u0[i][j][m][1]=u[i][j][m][1];
				u0[i][j][m][2]=u[i][j][m][2];
				rho[i][j][m]=0;psi[i][j][m]=0;
				u[i][j][m][0]=0;
				u[i][j][m][1]=0;
				u[i][j][m][2]=0;
	
				for(int k=0;k<19;k++)
					{
					//f[i][j][k]=F[i][j][k];
					rho[i][j][m]+=f[i][j][m][k];
					psi[i][j][m]+=g[i][j][m][k];
					u[i][j][m][0]+=e[k][0]*f[i][j][m][k];
					u[i][j][m][1]+=e[k][1]*f[i][j][m][k];
					u[i][j][m][2]+=e[k][2]*f[i][j][m][k];
					}
			

				
					

				u[i][j][m][0]=(u[i][j][m][0]+forcex[i][j][m]/2)/rho[i][j][m];
				u[i][j][m][1]=(u[i][j][m][1]+forcey[i][j][m]/2)/rho[i][j][m];
				u[i][j][m][2]=(u[i][j][m][2]+forcez[i][j][m]/2)/rho[i][j][m];
				
			}
			else
			{		u[i][j][m][0]=0; 
					u[i][j][m][1]=0;
					u[i][j][m][2]=0;
					psi[i][j][m]=0;

				markmac=0;
						for (int k=1;k<18;k++)
						{
						interi=i+e[k][0];
						interj=j+e[k][1];if (interj<0) {interj=NY;}; if (interj>NY) {interj=0;};
						interk=m+e[k][2];if (interk<0) {interk=NZ;}; if (interk>NZ) {interk=0;};
						
							
						if ((interi<0) && (swapl[interj*(NZ+1)+interk]==0) && (markmac==0))
							{
							
							psi[i][j][m]=swapl_psi[interj*(NZ+1)+interk]-2*ContactAngle_parameter;	
						
							markmac=1;
							}
						if ((interi>(NX+1)/mpi_size-1) && (swapr[interj*(NZ+1)+interk]==0) && (markmac==0))
							{
							
							psi[i][j][m]=swapr_psi[interj*(NZ+1)+interk]-2*ContactAngle_parameter;	
						
							markmac=1;
							}
						
						if ((interi>=0) && (interi<=(NX+1)/mpi_size-1) && (Solid[interi][interj][interk]==0) && (markmac==0))
							{
							
							psi[i][j][m]=psi[interi][interj][interk]-2*ContactAngle_parameter;	
						
							markmac=1;
							}



						}
			}

                        
			}
	delete [] swapl;
	delete [] swapr;
	delete [] swapl_psi;
	delete [] swapr_psi;

}




void standard_bounceback_boundary(int it, int jt,int kt, double**** f, double**** g)
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

			tmp = g[it][jt][kt][1];g[it][jt][kt][1] = g[it][jt][kt][2];g[it][jt][kt][2] = tmp;
			tmp = g[it][jt][kt][3];g[it][jt][kt][3] = g[it][jt][kt][4];g[it][jt][kt][4] = tmp;
                        tmp = g[it][jt][kt][5];g[it][jt][kt][5] = g[it][jt][kt][6];g[it][jt][kt][6] = tmp;
			tmp = g[it][jt][kt][7];g[it][jt][kt][7] = g[it][jt][kt][10];g[it][jt][kt][10] = tmp;
			tmp = g[it][jt][kt][8];g[it][jt][kt][8] = g[it][jt][kt][9];g[it][jt][kt][9] = tmp;
			tmp = g[it][jt][kt][11];g[it][jt][kt][11] = g[it][jt][kt][14];g[it][jt][kt][14] = tmp;
                        tmp = g[it][jt][kt][12];g[it][jt][kt][12] = g[it][jt][kt][13];g[it][jt][kt][13] = tmp;
			tmp = g[it][jt][kt][15];g[it][jt][kt][15] = g[it][jt][kt][18];g[it][jt][kt][18] = tmp;
			tmp = g[it][jt][kt][16];g[it][jt][kt][16] = g[it][jt][kt][17];g[it][jt][kt][17] = tmp;
			

}




void output_velocity_Vector(int m,double**** u,int NX,int NY,int NZ)	
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
	name<<outputfile<<"LBM_velocity_Vector_"<<m<<".vtk";
	ofstream out(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Velocity"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<NX+1<<"         "<<NY+1<<"         "<<NZ+1<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<(NX+1)*(NY+1)*(NZ+1)<<endl;
	out<<"VECTORS sample_vectors double"<<endl;
	out<<endl;
	//out<<"LOOKUP_TABLE default"<<endl;
	for(int k=0;k<=NZ;k++)
      		for(int j=0; j<=NY; j++)
			{
			for(int i=0;i<=NX;i++)
        		out<<rbuf_u[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3]<<" "<<rbuf_u[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+1]<<" "<<rbuf_u[i*(NY+1)*(NZ+1)*3+j*(NZ+1)*3+k*3+2]<<" "<<endl;
			//out<<endl;
			}
			


	}

	delete [] u_storage;
	if (rank==root_rank)
		delete [] rbuf_u;

}





void Geometry(bool*** Solid,int NX,int NY,int NZ)	
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
	

	if (rank==root_rank)
	{
	ostringstream name;
	name<<outputfile<<"LBM_Geometry"<<".vtk";
	ofstream out(name.str().c_str());
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
for(int k=0;k<=NZ;k++)
        for(int j=0; j<=NY; j++)
		for(int i=0;i<=NX;i++)
			//for(int k=0;k<=NZ;k++)
			out<<"		"<<rbuf[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]<<endl;
	}
		
	delete [] Solid_storage;
	if (rank==root_rank)
		delete [] rbuf;

		
}



void output_density(int m,double*** rho, int NX,int NY, int NZ)	
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
	
	if (rank==root_rank)
	{
	
	ostringstream name;
	name<<outputfile<<"LBM_Density_"<<m<<".vtk";
	ofstream out(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
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
			out<<"		"<<rbuf_rho[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]<<endl;
	}

	

	delete [] rho_storage;
	
	if (rank==root_rank)
		delete [] rbuf_rho;

		
}




void output_psi(int m,double*** psi, int NX,int NY, int NZ,bool*** Solid)	
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	const int mpi_size=MPI :: COMM_WORLD . Get_size ();
	const int root_rank=mpi_size-1;
	int Gather_Count=((NX+1)/mpi_size)*(NY+1)*(NZ+1);

	double* psi_storage= new double[Gather_Count];
	double* rbuf_psi;

	double rho0=0.0;

	if (rank==root_rank)
		rbuf_psi= new double[(NX+1)*(NY+1)*(NZ+1)];
	
	for (int i=0;i<(NX+1)/mpi_size;i++)
		for (int j=0;j<=NY;j++)
			for(int k=0;k<=NZ;k++)
			{
			if (Solid[i][j][k]==1)
			psi_storage[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=0.0;
			else
			psi_storage[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]=psi[i][j][k];
			}

	


	MPI_Gather(psi_storage,Gather_Count,MPI_DOUBLE,rbuf_psi,Gather_Count,MPI_DOUBLE,root_rank,MPI_COMM_WORLD);
	
	if (rank==root_rank)
	{
	
	ostringstream name;
	name<<outputfile<<"LBM_psi_"<<m<<".vtk";
	ofstream out(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
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
			out<<"		"<<rbuf_psi[i*(NY+1)*(NZ+1)+j*(NZ+1)+k]<<endl;
	}

	

	delete [] psi_storage;
	
	if (rank==root_rank)
		delete [] rbuf_psi;

		
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
			temp3 += sqrt(u[i][j][k][0]*u[i][j][k][0]+u[i][j][k][1]*u[i][j][k][1]+u[i][j][k][2]*u[i][j][k][2]);
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

double Comput_Perm(double*** psi,double**** u,double* Per_l,double* Per_g,int NX,int NY,int NZ,double gx,double gy,double gz,int PerDIr)
{

	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	double *rbuf_l,*rbuf_g;
	rbuf_l=new double[mpi_size*3];
	rbuf_g=new double[mpi_size*3];

	double Perm_l[3];
	double Perm_g[3];
	double error;
	double Q_l[3]={0.0,0.0,0.0};
	double Q_g[3]={0.0,0.0,0.0};

for(int i=0; i<(NX+1)/mpi_size; i++)
	for(int j=0;j<=NY;j++)
		for(int k=0;k<=NZ;k++)
		if (psi[i][j][k]>=0)
		{
		Q_l[0]+=u[i][j][k][0];
		Q_l[1]+=u[i][j][k][1];
		Q_l[2]+=u[i][j][k][2];
		}
		else	
		{
		Q_g[0]+=u[i][j][k][0];
		Q_g[1]+=u[i][j][k][1];
		Q_g[2]+=u[i][j][k][2];
		}


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


double Comput_Saturation(double*** psi,bool*** Solid, int NX,int NY,int NZ)
{
	int rank = MPI :: COMM_WORLD . Get_rank ();
	int mpi_size=MPI :: COMM_WORLD . Get_size ();

	double S_l,S_g;
	

	
double *rbuf_l,*rbuf_g;

	rbuf_l=new double[mpi_size];
	rbuf_g=new double[mpi_size];

	S_l=0;S_g=0;
	for(int i=0; i<(NX+1)/mpi_size; i++)
		for(int j=0;j<=NY;j++)
			for(int k=0;k<=NZ;k++)
			{				
			if ((psi[i][j][k]>=0) and (Solid[i][j][k]==0))
			S_l+=1;
			if ((psi[i][j][k]<0) and (Solid[i][j][k]==0))
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

	//MPI_Bcast(&S_l,1,MPI_DOUBLE,0,MPI_COMM_WORLD);
	//MPI_Bcast(&S_g,1,MPI_DOUBLE,0,MPI_COMM_WORLD);


	
	delete [] rbuf_l;
	delete [] rbuf_g;
	
	
	return (S_l);
			

}

