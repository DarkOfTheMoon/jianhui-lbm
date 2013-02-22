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
int ii,jj,kk;

char poreFileName[128]="maxd20-3-3.dat";
char poreFileNameVTK[128]="20-3-3.graph";

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
		cout << "\n The pore geometry file (" << poreFileName <<
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
		
	
	for(int k=0 ; k<nz ; k++)				
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)
		if (Solid[i][j][k]>0)
		for (int ls=0;ls<6;ls++)
		{
		ii=i+e[ls][0];
		jj=j+e[ls][1];
		kk=k+e[ls][2];
		if ((ii>=0) and (ii<nx) and (jj>=0) and (jj<ny) and (kk>=0) and (kk<nz))
			if (Solid[ii][jj][kk]>0)
				sum2++;
		}	
	
	sum2=sum2/2;


	cout<<sum-1<<"			"<<sum2<<endl;
	ostringstream name;
	name<<poreFileNameVTK;
	ofstream out;
	out.open(name.str().c_str());
	out<<sum-1<<"			"<<sum2<<endl;
	
	
	
	
	for (int k=0;k<nz;k++)
	for (int j=0;j<ny;j++)
	for (int i=0;i<nx;i++)
	if (Solid[i][j][k]>0)
	{
	for (int ls=0;ls<6;ls++)
		{
		ii=i+e[ls][0];
		jj=j+e[ls][1];
		kk=k+e[ls][2];	
		if ((ii>=0) and (ii<nx) and (jj>=0) and (jj<ny) and (kk>=0) and (kk<nz))
			if (Solid[ii][jj][kk]>0)
				out<<Solid[ii][jj][kk]<<" ";
		}

	out<<endl;
	}
	

	
	out.close();


}


