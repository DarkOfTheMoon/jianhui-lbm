#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<math.h> 

using namespace std; 
      
int main (int argc , char * argv [])
{
int nx=90;
int ny=90;
int nz=90;
int Zoom=4; //1,2,3,4...
int sum=0;

int*** Solid_Int;
bool*** Solid;
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

int pore;
	int i, j, k,ci,cj,ck;
	
	Solid_Int = new int**[nx*Zoom];	///*********
	
	for (i=0;i<nx*Zoom;i++)				///*********
		Solid_Int[i]=new int*[ny*Zoom];

	Solid_Int[0][0]=new int[nx*ny*nz*Zoom*Zoom*Zoom];

	
 	for (int i=1;i<ny*Zoom;i++)
               Solid_Int[0][i]=Solid_Int[0][i-1]+ny*Zoom;
       
       for (int i=1;i<nx*Zoom;i++)
       {
               Solid_Int[i][0]=Solid_Int[i-1][0]+ny*nz*Zoom*Zoom;
               for (int j=1;j<ny*Zoom;j++)
                       Solid_Int[i][j]=Solid_Int[i][j-1]+nz*Zoom;
       }




		
	
		

		
	for(k=0 ; k<nz ; k++)				///*********
	for(j=0 ; j<ny ; j++)
	for(i=0 ; i<nx ; i++)				///*********


	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore;
			
			//if (pore == 0.0)	{Solid[ci-1][cj-1][ck-1] = 0;}
			if (pore == 0)	
				{
				for (ci=0;ci<Zoom;ci++)
				for (cj=0;cj<Zoom;cj++)
				for (ck=0;ck<Zoom;ck++) 
					Solid_Int[i*Zoom+ci][j*Zoom+cj][k*Zoom+ck] = 0;
				}
			
			//if (pore == 1.0) 	{Solid[ci-1][cj-1][ck-1] = 1;sum++;}
			//if (pore == 1) 	{Solid[i][j][k] = 1;sum++;}
			
			if (pore == 1)	
				{sum++;
				for (ci=0;ci<Zoom;ci++)
				for (cj=0;cj<Zoom;cj++)
				for (ck=0;ck<Zoom;ck++) 
					Solid_Int[i*Zoom+ci][j*Zoom+cj][k*Zoom+ck] = 1;
				}
			
			
		}
		
	fin.close();
		
	
	cout<<"Porosity = "<<1-(double(sum)/(nx*ny*nz))<<endl;	
	//cout<<sum<<endl;	

	ostringstream name;
	//name<<poreFileName<<"_sym.dat";
	name<<"output2.vtk";
	ofstream out;
	out.open(name.str().c_str());
	
	
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	//out<<"binary"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<nx*Zoom<<"         "<<ny*Zoom<<"         "<<nz*Zoom<<endl;       ///*********
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<nx*ny*nz*Zoom*Zoom*Zoom<<endl;				///*********
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	

	// out.write((char *)(&Solid_Int[0][0][0]), sizeof(int)*nx*ny*nz*Zoom*Zoom*Zoom); 
		for(k=0 ; k<nz*Zoom ; k++)
		{cout<<k<<endl;						///*********
		for(j=0 ; j<ny*Zoom ; j++)					///*********
			for(i=0 ; i<nx*Zoom; i++)				///*********		
			out<<Solid_Int[i][j][k]<<" ";
		}
	



	out.close();

cout<<nx*Zoom<<"         "<<ny*Zoom<<"         "<<nz*Zoom<<endl;




}
