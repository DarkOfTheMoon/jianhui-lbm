#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>

using namespace std; 
      
int main (int argc , char * argv [])
{

int nx=128;
int ny=128;
int nz=128;
int dir=3;

int pls[4][3]={{0,2,2},{2,0,2},{2,2,0},{0,0,0}};

int sum=0;

bool*** Solid_Int;
char poreFileName[128]="128.all";

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


double pore;
	int i, j, k;
	
	Solid_Int = new bool**[nx];			
	for (i=0;i<nx;i++)				
		{
		Solid_Int[i]=new bool*[ny+pls[dir][1]];		///*********
		for (j=0;j<ny;j++)			///*********
			{
			Solid_Int[i][j]= new bool[nz+pls[dir][2]];///*********
			for (k=0;k<nz;k++)			///*********
				Solid_Int[i][j][k]= 0;
			}
		}
	
	
		while(!fin.eof())
		{	
			fin >> i >> j>> k>> pore;
			Solid_Int[i-1][j-1][k-1]=pore;
		}
		
	
	fin.close();


//***********************************************************************************************
	//cout<<sum<<endl;

	ostringstream name;
	name<<"Bentheimer"<<".vtk";
	ofstream out;
	out.open(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<nx+pls[dir][0]<<"         "<<ny+pls[dir][1]<<"         "<<nz+pls[dir][2]<<endl;       ///*********
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<(nx+pls[dir][0])*(ny+pls[dir][1])*(nz+pls[dir][2])<<endl;				///*********
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;

	for(k=0 ; k<nz+pls[dir][2] ; k++)						///*********
		for(j=0 ; j<ny+pls[dir][1] ; j++)					///*********
			for(i=0 ; i<nx+pls[dir][0] ; i++)				///*********
			out<<Solid_Int[i][j][k]<<endl;


	out.close();


//***********************************************************************************************



//=======================================================================================
dir=0;                 //%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int div=2; 			//%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
int nx_d=nx/div;
int ny_d=ny/div;
int nz_d=nz/div;
int sym_x=0;
int sym_y=0;
int sym_z=0;
int add_buf_x_n=5;
int add_buf_y_n=0;
int add_buf_z_n=0;

int add_buf_x_p=0;
int add_buf_y_p=0;
int add_buf_z_p=0;
bool*** Solid;

	
	Solid = new bool**[(nx_d+pls[dir][0])*(sym_x+1)];			
	for (i=0;i<(nx_d+pls[dir][0])*(sym_x+1);i++)				
		{
		Solid[i]=new bool*[(ny_d+pls[dir][1])*(sym_y+1)];		
		for (j=0;j<(ny_d+pls[dir][1])*(sym_y+1);j++)			
			{
			Solid[i][j]= new bool[(nz_d+pls[dir][2])*(sym_z+1)];
			for (k=0;k<(nz_d+pls[dir][2])*(sym_z+1);k++)			
				Solid[i][j][k]= 1;
			}
		}

//cout<<nx_d<<endl;

for (int i_x=1;i_x<=div;i_x++)
for (int j_x=1;j_x<=div;j_x++)
for (int k_x=1;k_x<=div;k_x++)
	{
	for (i=1;i<=nx_d;i++)
		for (j=1;j<=ny_d;j++)
			for (k=1;k<=nz_d;k++)
			Solid[i-1+pls[dir][0]/2][j-1+pls[dir][1]/2][k-1+pls[dir][2]/2]=Solid_Int[i+(i_x-1)*nx_d-1][j+(j_x-1)*ny_d-1][k+(k_x-1)*nz_d-1];
	
		
		
		
		
if (sym_x==1)
	for(k=0 ; k<nz_d+pls[dir][2] ; k++)				///*********
	for(j=0 ; j<ny_d+pls[dir][1] ; j++)
	for(i=0 ; i<nx_d+pls[dir][0] ; i++)				///*********
		Solid[nx_d+pls[dir][0]+i][j][k]=Solid[nx_d+pls[dir][0]-1-i][j][k];


if (sym_y==1)
	for(k=0 ; k<nz_d+pls[dir][2] ; k++)				///*********
	for(j=0 ; j<ny_d+pls[dir][1] ; j++)
	for(i=0 ; i<nx_d+pls[dir][0] ; i++)				///*********
		Solid[i][ny_d+pls[dir][1]+j][k]=Solid[i][ny_d+pls[dir][1]-1-j][k];


if (sym_z==1)
	for(k=0 ; k<nz_d+pls[dir][2] ; k++)				///*********
	for(j=0 ; j<ny_d+pls[dir][1] ; j++)
	for(i=0 ; i<nx_d+pls[dir][0] ; i++)				///*********
		Solid[i][j][nz_d+pls[dir][2]+k]=Solid[i][j][nz_d+pls[dir][2]-1-k];


	//ostringstream name;
	name.str("");
	name<<"Bentheimer_"<<(i_x-1)*div*div+(j_x-1)*div+(k_x-1)<<".vtk";
	ofstream out;
	out.open(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	//out<<"DIMENSIONS         "<<nx_d+pls[dir][0]<<"         "<<ny_d+pls[dir][1]<<"         "<<nz_d+pls[dir][2]<<endl; 
	out<<"DIMENSIONS         "<<(nx_d+pls[dir][0])*(sym_x+1)+add_buf_x_n+add_buf_x_p<<"         "<<(ny_d+pls[dir][1])*(sym_y+1)+add_buf_y_n+add_buf_y_p<<"         "<<(nz_d+pls[dir][2])*(sym_z+1)+add_buf_z_n+add_buf_z_p<<endl;       ///*********
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	//out<<"POINT_DATA     "<<(nx_d+pls[dir][0])*(ny_d+pls[dir][1])*(nz_d+pls[dir][2])<<endl;			
	out<<"POINT_DATA     "<<((nx_d+pls[dir][0])*(sym_x+1)+add_buf_x_n+add_buf_x_p)*((ny_d+pls[dir][1])*(sym_y+1)+add_buf_y_n+add_buf_y_p)*((nz_d+pls[dir][2])*(sym_z+1)+add_buf_z_n+add_buf_z_p)<<endl;		
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;

//	for(k=0 ; k<nz_d+pls[dir][2] ; k++)						///*********
//		for(j=0 ; j<ny_d+pls[dir][1] ; j++)					///*********
//			for(i=0 ; i<nx_d+pls[dir][0] ; i++)				///*********
//			out<<Solid[i][j][k]<<endl;

		for(k=0 ; k<(nz_d+pls[dir][2])*(sym_z+1)+add_buf_z_n+add_buf_z_p ; k++)						///*********
		for(j=0 ; j<(ny_d+pls[dir][1])*(sym_y+1)+add_buf_y_n+add_buf_y_p ; j++)					///*********
			for(i=0 ; i<(nx_d+pls[dir][0])*(sym_x+1)+add_buf_x_n+add_buf_x_p ; i++)				///*********		
			if ((i>=add_buf_x_n) and (i<(nx_d+pls[dir][0])*(sym_x+1)+add_buf_x_n) and (j>=add_buf_y_n) and (j<(ny_d+pls[dir][1])*(sym_y+1)+add_buf_y_n) and (k>=add_buf_z_n) and (k<(nz_d+pls[dir][2])*(sym_z+1)+add_buf_z_n))
			out<<Solid[i-add_buf_x_n][j-add_buf_y_n][k-add_buf_z_n]<<endl;
			else
			out<<0.0<<endl;
	out.close();

			
	}

//=========================================================================================














}


