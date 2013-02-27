#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include <math.h>


using namespace std; 

int main (int argc , char * argv [])
{
int NCHAR=128;
	char     filename[128], dummy[128+1];
	int      dummyInt;


int sym_x;
int sym_y;
int sym_z;
int add_buf_x_n;
int add_buf_y_n;
int add_buf_z_n;

int add_buf_x_p;
int add_buf_y_p;
int add_buf_z_p;
int add_porous_plate; //0=OFF, 1=fine plate,pore size1, pore size2, 3=posr size3
int porous_position; //-1=defualt position,end of the geometry, or give a positive value
int Zoom; //1,2,3,4...


int nx,ny,nz;
int expvtk,expdat,bindat,fil,geo_mod;
int dir;

char poreFileName[128];
char poreFileNameMET[128];
char poreFileNameVTK[128];
char poreFileNameOut[128];

int mesh_par;
int par_n;
char poreFileNamePar[128];
int expparvtk;
int par_bin;
	
ifstream fins(argv[1]);
			
							fins.getline(dummy, NCHAR);
	fins >> poreFileName;				fins.getline(dummy, NCHAR);
	fins >> nx>> ny>>nz;				fins.getline(dummy, NCHAR);
	fins >> poreFileNameVTK;			fins.getline(dummy, NCHAR);
	fins >> poreFileNameOut;			fins.getline(dummy, NCHAR);
	fins >> expvtk;					fins.getline(dummy, NCHAR);
	fins >> expdat;					fins.getline(dummy, NCHAR);
	fins >> bindat;					fins.getline(dummy, NCHAR);
							fins.getline(dummy, NCHAR);
	fins >> fil;					fins.getline(dummy, NCHAR);
							fins.getline(dummy, NCHAR);
	fins >> geo_mod;					fins.getline(dummy, NCHAR);
	fins >> dir;					fins.getline(dummy, NCHAR);
	fins >> sym_x >> sym_y >> sym_z;			fins.getline(dummy, NCHAR);
	fins >> add_buf_x_n>>add_buf_y_n>>add_buf_z_n;	fins.getline(dummy, NCHAR);
	fins >> add_buf_x_p>>add_buf_y_p>>add_buf_z_p;	fins.getline(dummy, NCHAR);
	fins >> add_porous_plate;			fins.getline(dummy, NCHAR);
	fins >> porous_position;				fins.getline(dummy, NCHAR);
	fins >> Zoom;					fins.getline(dummy, NCHAR);
							fins.getline(dummy, NCHAR);
	fins >> mesh_par;				fins.getline(dummy, NCHAR);
	fins >> poreFileNameMET;				fins.getline(dummy, NCHAR);
								fins.getline(dummy, NCHAR);
	fins >> par_n;						fins.getline(dummy, NCHAR);
	fins >> par_bin;				fins.getline(dummy, NCHAR);
	fins >> expparvtk;					fins.getline(dummy, NCHAR);




int ii,jj,kk,pore2;

//=============
int sumss[11];
//=============

int pls[7][3]={{0,2,2},{2,0,2},{2,2,0},{0,0,0},{0,0,0},{0,0,0}};

int nx1,ny1,nz1;
nx1=(nx+pls[dir][0])*(sym_x+1)+add_buf_x_n+add_buf_x_p;
ny1=(ny+pls[dir][1])*(sym_y+1)+add_buf_y_n+add_buf_y_p;
nz1=(nz+pls[dir][2])*(sym_z+1)+add_buf_z_n+add_buf_z_p;

if (geo_mod==1)
{
nx=nx1*Zoom;
ny=ny1*Zoom;
nz=nz1*Zoom;

cout<<nx<<"	"<<ny<<"	"<<nz<<endl;
}


char poreFileName2;


const int e[18][3]=
{{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};


int sum=0;
int sum2=0;
int* Solid_rank0;
 
int*** Solid;
double pore;
	


	Solid = new int**[nx];	///*********
	
	for (int i=0;i<nx;i++)				///*********
		Solid[i]=new int*[ny];

	Solid[0][0]=new int[nx*ny*nz];

	
 	for (int i=1;i<ny;i++)
               Solid[0][i]=Solid[0][i-1]+nz;
       
       for (int i=1;i<nx;i++)
       {
               Solid[i][0]=Solid[i-1][0]+ny*nz;
               for (int j=1;j<ny;j++)
                       Solid[i][j]=Solid[i][j-1]+nz;
       }

		
	

if (bindat==0)
{
sum=1;
	FILE *ftest;
	ifstream fin;
	
	ftest = fopen(poreFileNameOut, "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file  (" << poreFileNameOut <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest);

	fin.open(poreFileNameOut);

	for(int k=0 ; k<nz ; k++)				///*********
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)				///*********


	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore;
			
			//if (pore == 0.0)	{Solid[ci-1][cj-1][ck-1] = 0;}
			if (pore == 0.0)	{Solid[i][j][k] = sum;sum++;}
			
			//if (pore == 1.0) 	{Solid[ci-1][cj-1][ck-1] = 1;sum++;}
			if (pore == 1.0) 	{Solid[i][j][k] = 0;}
			
		
			
			
		}
		
	fin.close();
}
else
{
	sum=1;

	fstream fin;
	fin.open(poreFileNameOut,ios::in);
	if (fin.fail())
	        {
	        cout<<"\n file open error on " << poreFileNameOut<<endl;
	        exit(-1);
	        }
	Solid_rank0 = new int[nx*ny*nz];
	fin.read((char *)(&Solid_rank0[0]), sizeof(int)*nx*ny*nz);

	for(int k=0 ; k<nz ; k++)				///*********
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)				///*********


	                         
		{	
			
		
			
			//if (pore == 0.0)	{Solid[ci-1][cj-1][ck-1] = 0;}
			if (Solid_rank0[i*ny*nz+j*nz+k] == 0)	{Solid[i][j][k] = sum;sum++;}
			else
			
			if (Solid_rank0[i*ny*nz+j*nz+k]==1) 	{Solid[i][j][k] = 0;}
			
		
			
			
		}
	fin.close();
	
}

	ostringstream namemet;
	namemet<<poreFileNameMET<<".part."<<par_n;
FILE *ftest;
	ifstream fin;		
ftest = fopen(namemet.str().c_str(), "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file (" << namemet.str().c_str() <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest);

	fin.open(namemet.str().c_str());


	sum2=0;
	for(int k=0 ; k<nz ; k++)				///*********
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)
		if (Solid[i][j][k]>0)
			{
			fin >> pore2;		
			Solid[i][j][k]=pore2+1;
			}



	fin.close();

	
	
	
	
	if (expparvtk==1)
	{
	ostringstream name;
	name<<"Partitions.vtk";
	ofstream out;
	out.open(name.str().c_str());

	
	
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"binary"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<nz<<"         "<<ny<<"         "<<nx<<endl;       ///*********
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<nx*ny*nz<<endl;				///*********
	out<<"SCALARS sample_scalars int"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	
	out.write((char *)(&Solid[0][0][0]), sizeof(int)*nx*ny*nz); 
	/*
	for (int k=0;k<nz;k++)
	{
	for (int j=0;j<ny;j++)
	for (int i=0;i<nx;i++)
		out<<Solid[i][j][k]<<" ";
	}
	}
	*/
	
	

	
	out.close();
	}

	cout<<"Start writing Partitioned Geomtry DAT file"<<endl;
	cout<<nx<<"	"<<ny<<"	"<<nz<<endl;
	ostringstream name2;
	name2<<poreFileName<<"."<<par_n<<"."<<nx<<"x"<<ny<<"x"<<nz;
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
	ofstream out2;
	out2.open(name2.str().c_str());
	
	if (par_bin==1)
	out2.write((char *)(&Solid[0][0][0]), sizeof(int)*nx*ny*nz); 
	else
	for (int k=0;k<nz;k++)
	{
		//cout<<k<<endl;
		for (int j=0;j<ny;j++)
		for (int i=0;i<nx;i++)
			out2<<Solid[i][j][k]<<" ";
	}
	
	out2.close();

	cout<<"DAT file ouput complete"<<endl;
	cout<<endl;


}


