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

int nx=90;
int ny=90;
int nz=90;
int dir=2;

int pls[4][3]={{0,2,2},{2,0,2},{2,2,0},{0,0,0}};

int sum=0;

bool*** Solid_Int;
char poreFileName[128]="maxd20-3-3.dat";

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
	
	Solid_Int = new bool**[nx+pls[dir][0]];			///*********
	for (i=0;i<nx+pls[dir][0];i++)				///*********
		{
		Solid_Int[i]=new bool*[ny+pls[dir][1]];		///*********
		for (j=0;j<ny+pls[dir][1];j++)			///*********
			{
			Solid_Int[i][j]= new bool[nz+pls[dir][2]];///*********
			for (k=0;k<nz+pls[dir][2];k++)			///*********
				Solid_Int[i][j][k]= 1;
			}
		}
	
	// Reading pore geometry
	for(k=pls[dir][2]/2 ; k<nz+pls[dir][2]/2 ; k++)				///*********
	for(j=pls[dir][1]/2 ; j<ny+pls[dir][1]/2 ; j++)
	for(i=pls[dir][0]/2 ; i<nx+pls[dir][0]/2 ; i++)				///*********
	{
		while(true)
		{	
			fin >> pore;
			if( pore == 0.0 || pore == 1.0) break;
		}
		if (pore == 0.0)	{Solid_Int[i][j][k] = 0;sum++;}
		//else			Solid_Int[i][j][k] = 1;
		if (pore == 1.0) 	{Solid_Int[i][j][k] = 1;sum++;}
	}
	fin.close();


	//cout<<sum<<endl;

	ostringstream name;
	name<<"maxd20-3-3_x_NewBC"<<".vtk";
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
			out<<"		"<<Solid_Int[i][j][k]<<endl;


	out.close();


}


