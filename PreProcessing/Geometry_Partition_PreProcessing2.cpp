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

int nx=90;
int ny=90;
int nz=90;
int ii,jj,kk,pore2;

//=============
int sumss[11];
//=============



char poreFileName[128]="ftb.dat";
char poreFileName2[128]="20-3-3.graph.part.10";
char poreFileNameVTK[128]="20-3-3.vtk";

const int e[18][3]=
{{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};


int sum=0;
int sum2=0;

int*** Solid;
double pore;


	FILE *ftest;
	ifstream fin;
	
	ftest = fopen(poreFileName, "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file  (" << poreFileName <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest);

	fin.open(poreFileName);



	Solid = new int**[nx];
	
	
		
	for (int i=0; i<nx;i++)
	{
	       Solid[i] = new int*[ny];
	       for (int j=0;j<ny;j++)
	       {
	               Solid[i][j] = new int[nz];
	               for (int k=0;k<nz;k++)
	                       Solid[i][j][k] = 0;
	       }
	}
		
	


sum=1;

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
		
ftest = fopen(poreFileName2, "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file (" << poreFileName2 <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest);

	fin.open(poreFileName2);


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

	/*
	sum2=0;
	for (int i=0;i<=10;i++)
	        sumss[i]=0;
	
	fin.open(poreFileName2);
	while (!fin.eof())
	{        
	        fin>>pore2;
	        sum2++;
	        for (int i=0;i<=10;i++)
	                if (pore2==i)
	                        sumss[i]++;
	}

	cout<<sum<<"                "<<sum2<<endl;
	for (int i=0;i<=10;i++)
	        cout<<i<<"                "<<sumss[i]<<endl;
	cout<<endl;
	*/
	
	
	
	
	ostringstream name;
	name<<poreFileNameVTK;
	ofstream out;
	out.open(name.str().c_str());

	
	
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<nx<<"         "<<ny<<"         "<<nz<<endl;       ///*********
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<nx*ny*nz<<endl;				///*********
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	
	
	for (int k=0;k<nz;k++)
	{
	for (int j=0;j<ny;j++)
	for (int i=0;i<nx;i++)
		out<<Solid[i][j][k]<<" ";
	}
	
	
	

	
	out.close();


}


