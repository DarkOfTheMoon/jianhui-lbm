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

int nx=220;
int ny=200;
int nz=200;

int Phase_Filter_Visual=0;	//1= defualt setting for filter phases for visualization


int sym_x=1;
int sym_y=0;
int sym_z=0;

int xn=20; 		//20,200,20,200,20,200
int xp=200;
int yn=20;
int yp=200;
int zn=20;
int zp=200;


int input_vtk=1;	//0=NO,1=YES
int output=0;		//export processed file? 
int output_vtk=1;
	
char poreFileName[128]="psi2_LBM_psi_450000.vtk";

char outputFileName[128]="sys_cut_phase_0.50.vtk";



//======================================================

if (Phase_Filter_Visual==1)
	{
	sym_x=0;sym_y=0;sym_z=0;
	xn=0;xp=nx;yn=0;yp=ny;zn=0;zp=nz;
	input_vtk=1;output=1;output_vtk=1;
	}


int NCHAR=128;

char     dummy[128+1];

int nx2,ny2,nz2;
nx2=xp-xn;
ny2=yp-yn;
nz2=zp-zn;







double*** Solid_Int;
double*** Solid;


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
	int i, j, k,ci,cj,ck;
	
	Solid_Int = new double**[nx2*(sym_x+1)];	///*********
	Solid = new double**[nx2];
	for (i=0;i<nx2*(sym_x+1);i++)				///*********
		{
		Solid_Int[i]=new double*[ny2*(sym_y+1)];		///*********
		for (j=0;j<ny2*(sym_y+1);j++)			///*********
			{
			Solid_Int[i][j]= new double[nz2*(sym_z+1)];///*********
			for (k=0;k<nz2*(sym_z+1);k++)			///*********
				Solid_Int[i][j][k]= 1;
			      
			}
		}
		
	for (i=0; i<nx2;i++)
	{
	       Solid[i] = new double*[ny2];
	       for (j=0;j<ny2;j++)
	       {
	               Solid[i][j] = new double[nz2];
	               for (k=0;k<nz2;k++)
	                       Solid[i][j][k] = 0;
	       }
	}
		
		

	if (input_vtk==1)
	{
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	}

	for(k=0 ; k<nz ; k++)				///*********
	for(j=0 ; j<ny ; j++)
	for(i=0 ; i<nx ; i++)				///*********

	 
	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore;
			if ((i>=xn) and (i<xp) and (j>=yn) and (j<yp) and (k>=zn) and (k<zp))
				Solid[i-xn][j-yn][k-zn]=pore;
			
			//if (pore == 0.0)	{Solid[ci-1][cj-1][ck-1] = 0;}
			//if (pore == 0.0)	{Solid[i][j][k] = 0;}
			
			//if (pore == 1.0) 	{Solid[ci-1][cj-1][ck-1] = 1;sum++;}
			//if (pore == 1.0) 	{Solid[i][j][k] = 1;sum++;}
			
			
			
			
		}
		
	fin.close();


	cout<<endl;	
	cout<<"File Reading Complete"<<endl;
	//cout<<"Porosity = "<<1-(double(sum)/(nx2*ny2*nz2))<<endl;	
	//cout<<sum<<endl;	
		
	int sum=0;
	int sum2=0;	
	// Reading pore geometry
	for(k=0 ; k<nz2 ; k++)				///*********
	for(j=0 ; j<ny2 ; j++)
	for(i=0 ; i<nx2 ; i++)				///*********

	{
	        Solid_Int[i][j][k]=Solid[i][j][k];

		if ((Solid[i][j][k]>0.00001) or (Solid[i][j][k]<-0.0001))
			sum++;
		
		if (Solid[i][j][k]>0.2)
			sum2++;
	        
	        
	}

	cout<<"Saturation of component 1 is "<<(double)sum2/sum<<endl;

	if (sym_x==1)
	for(k=0 ; k<nz2 ; k++)				///*********
	for(j=0 ; j<ny2 ; j++)
	for(i=0 ; i<nx2 ; i++)				///*********
		Solid_Int[nx2+i][j][k]=Solid_Int[nx2-1-i][j][k];


if (sym_y==1)
	for(k=0 ; k<nz2; k++)				///*********
	for(j=0 ; j<ny2; j++)
	for(i=0 ; i<nx2; i++)				///*********
		Solid_Int[i][ny2+j][k]=Solid_Int[i][ny2-1-j][k];


if (sym_z==1)
	for(k=0 ; k<nz2 ; k++)				///*********
	for(j=0 ; j<ny2; j++)
	for(i=0 ; i<nx2; i++)				///*********
		Solid_Int[i][j][nz2+k]=Solid_Int[i][j][nz2-1-k];



	if (output==1)
	{
	//cout<<sum<<endl;
	cout<<"Start writting output file"<<endl;
	cout<<endl;

	ostringstream name;
	name<<outputFileName;
	ofstream out;
	out.open(name.str().c_str());
	
	if (output_vtk==1)
	{
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<nx2*(sym_x+1)<<"         "<<ny2*(sym_y+1)<<"         "<<nz2*(sym_z+1)<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<nx2*(sym_x+1)*ny2*(sym_y+1)*nz2*(sym_z+1)<<endl;				///*********
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	}
	
	
	for(k=0 ; k<nz2*(sym_z+1) ; k++)						///*********
		for(j=0 ; j<ny2*(sym_y+1) ; j++)					///*********
			for(i=0 ; i<nx2*(sym_x+1); i++)				///*********		
			if (Solid_Int[i][j][k]>=0.0)
				out<<1.0<<endl;
			else
				//if (Solid_Int[i][j][k]<0)
				out<<-1.0<<endl;
				//else
				//out<<0.0<<endl;
				
				//out<<Solid_Int[i][j][k]<<endl;
			



	out.close();
	}

cout<<nx2*(sym_x+1)<<"         "<<ny2*(sym_y+1)<<"         "<<nz2*(sym_z+1)<<endl;

}


