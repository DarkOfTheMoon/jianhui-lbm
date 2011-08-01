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

int nx=200;
int ny=200;
int nz=200;
int dir=1;

int pls[4][3]={{0,2,2},{2,0,2},{2,2,0},{0,0,0}};

int sum=0;

bool*** Solid_Int;
char poreFileName[128]="maxd9-3-3.dat";

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
	
	Solid_Int = new bool**[nx+2];			///*********
	for (i=0;i<nx+2;i++)				///*********
		{
		Solid_Int[i]=new bool*[ny+2];		///*********
		for (j=0;j<ny+2;j++)			///*********
			{
			Solid_Int[i][j]= new bool[nz+2];///*********
			for (k=0;k<nz+2;k++)			///*********
				Solid_Int[i][j][k]= 0;
			}
		}
	
	// Reading pore geometry
	for(k=1 ; k<nz+1 ; k++)				///*********
	for(j=1 ; j<ny+1 ; j++)
	for(i=1 ; i<nx+1 ; i++)				///*********
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

	if (pls[dir][0]>0)
		for (j=0;j<ny+2;j++)
			for (k=0;k<nz+2;k++)
			{
			Solid_Int[0][j][k]=1;
			Solid_Int[nx+1][j][k]=1;
			}


	if (pls[dir][1]>0)
		for (i=0;i<nx+2;i++)
			for (k=0;k<nz+2;k++)
			{
			Solid_Int[i][0][k]=1;
			Solid_Int[i][ny+1][k]=1;
			}



	if (pls[dir][2]>0)
		for (j=0;j<ny+2;j++)
			for (i=0;i<nx+2;i++)
			{
			Solid_Int[i][j][0]=1;
			Solid_Int[i][j][nz+1]=1;
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
	out<<"DIMENSIONS         "<<nx+2<<"         "<<ny+2<<"         "<<nz+2<<endl;       ///*********
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<(nx+2)*(ny+2)*(nz+2)<<endl;				///*********
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;

	for(k=0 ; k<nz+2 ; k++)						///*********
		for(j=0 ; j<ny+2 ; j++)					///*********
			for(i=0 ; i<nx+2 ; i++)				///*********
			out<<"		"<<Solid_Int[i][j][k]<<endl;


	out.close();


}


