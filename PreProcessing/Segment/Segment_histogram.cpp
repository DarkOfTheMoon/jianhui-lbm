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

int nx=300;
int ny=300;
int nz=1;

char poreFileName[128]="300x300_3phaseket-filtxt.txt";
char poreFileNameVTK[128]="segment.vtk";
char poreFileNameOut[128]="segment.dat";
//output VTK file,0=no, 1=yes
int VTK_OUT=1;
int local_min=10; //half of local min
//===========VTK AND OUT ARE ALL WRITTEN IN BINARY FORMAT===============
//===========================================================





int histo[256];
double*** Solid;
double*** cas;
double*** cas2;
double*** cas3;
int*** seg;
double his1[256];
double his2[256];
double his3[256];
double his4[256];
double his5[256];



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


int pore;
	int i, j, k,ci,cj,ck;
	
	
	Solid = new double**[nx];
	cas = new double**[nx];cas2 = new double**[nx];cas3 = new double**[nx];
	seg= new int**[nx];
		
	for (i=0; i<nx;i++)
	{
	       Solid[i] = new double*[ny];
		cas[i] = new double*[ny];cas2[i] = new double*[ny];cas3[i] = new double*[ny];
		seg[i] = new int*[ny];

	       for (j=0;j<ny;j++)
	       {
	               Solid[i][j] = new double[nz];
			cas[i][j] = new double[nz];cas2[i][j] = new double[nz];cas3[i][j] = new double[nz];
			seg[i][j] = new int[nz];


	               for (k=0;k<nz;k++)
				{
	                       	Solid[i][j][k] = 0;
				cas[i][j][k] = 0;cas2[i][j][k] = 0;cas3[i][j][k] = 0;
				seg[i][j][k] = 0;
				}
	       }
	}
		
	

for (i=0;i<256;i++)
	{histo[i]=0;his1[i]=0.0;his4[i]=0.0;}



	for(k=0 ; k<nz ; k++)				///*********
	for(j=0 ; j<ny ; j++)
	for(i=0 ; i<nx ; i++)				///*********


	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore;
			
			//if (pore == 0.0)	{Solid[ci-1][cj-1][ck-1] = 0;}
			//if (pore == 0.0)	{Solid[i][j][k] = 0;}
			
			//if (pore == 1.0) 	{Solid[ci-1][cj-1][ck-1] = 1;sum++;}
			//if (pore == 1.0) 	{Solid[i][j][k] = 1;sum++;}
			
			Solid[i][j][k]=(double) pore;
			
			seg[i][j][k]=pore;
			
			histo[pore]+=1;
		}
		
	fin.close();
		
	for (i=1;i<255;i++)
		{
		his1[i]=histo[i+1]-histo[i];
		his2[i]=(histo[i+1]-histo[i-1])/2;
		his3[i]=(histo[i+1]+histo[i-1]-2*histo[i]);
		}
	
	//for (i=1;i<255;i++)
	for (i=2;i<254;i++)
		his1[i]=(his2[i-1]+his2[i]+his2[i+1])/3;
		his1[i]=(his2[i-2]+his2[i-1]+his2[i]+his2[i+1]+his2[i+2])/5;

double max_val=0.0;	
double max_val2=0.0;
double max1,max2,max3;
int sum1;


	
	for (i=0;i<256;i++)
		{
		if (histo[i]>max_val)
			max_val=histo[i];
		//if (histo[255-i]>max_val2)
		//	max_val2=histo[i];
		his3[i]=max_val;//his4[i]=max_val2;
		}


	for (i=local_min;i<256-local_min;i++)
		{
		sum1=7000000;
		for (j=i-local_min;j<i+local_min;j++)
			if (histo[j]<sum1)
				sum1=histo[j];	

		his4[i]=sum1;
		}

	if (VTK_OUT==1)
	{
	cout<<"Start writing VTK file"<<endl;
	cout<<nx<<"         "<<ny<<"         "<<nz<<endl;

	ostringstream name;
	name<<poreFileNameVTK;
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
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
	out<<"SCALARS sample_scalars int"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	
	
	for (k=0;k<nz;k++)
	{
	cout<<k<<endl;
	for (j=0;j<ny;j++)
	for (i=0;i<nx;i++)
		out<<seg[i][j][k]<<" ";
	}
	

	//out.write((char *)(&seg[0][0][0]), sizeof(int)*nx*ny*nz); 

	out.close();

	cout<<"VTK file ouput COMPLETE"<<endl;
	}


	//=================================================================
	cout<<"Start writing DAT file"<<endl;
	cout<<nx<<"         "<<ny<<"         "<<nz<<endl;
	ostringstream name2;
	name2<<poreFileNameOut;
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
	ofstream out2;
	out2.open(name2.str().c_str());
	
	
	out2.write((char *)(&seg[0][0][0]), sizeof(int)*nx*ny*nz); 

	out2.close();

	cout<<"DAT file ouput complete"<<endl;
	



	ostringstream name3;
	name3<<"histo_test.dat";
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
	ofstream out3;
	out3.open(name3.str().c_str());
	
	for (i=0;i<256;i++)
		out3<<histo[i]<<endl;
	
	out3.close();

	name3.str("");
	name3<<"histo_test2.dat";
	out3.open(name3.str().c_str());
	for (i=0;i<256;i++)
		out3<<his2[i]<<endl;
	
	out3.close();


	name3.str("");
	name3<<"histo_test3.dat";
	out3.open(name3.str().c_str());
	for (i=0;i<256;i++)
		out3<<his1[i]<<endl;
	
	out3.close();

	name3.str("");
	name3<<"histo_test4.dat";
	out3.open(name3.str().c_str());
	for (i=0;i<256;i++)
		out3<<his4[i]<<endl;
	
	out3.close();
			
}


