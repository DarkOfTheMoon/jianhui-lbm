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

int nx=310;
int ny=300;
int nz=300;
char poreFileName[128]="LV60_310_300_300.dat";
char poreFileName2[128]="LV60_cut";
int partition=3;
int read_mod=1;  //0=decimal, 1=binary
int add_solid=1;	//0=no,1=yes
int add_buf_x_n=0;

int add_buf_x_p=0;





int num_par;
num_par=partition*partition;






int sum=0;
int nxs=nx+add_buf_x_n+add_buf_x_p;

int* Solid;
int*** Solid2;
double pore;
	int i, j, k,ci,cj,ck;
Solid = new int[nx*ny*nz];

Solid2 = new int**[nxs];	///*********
	
	for (i=0;i<nxs;i++)				///*********
		Solid2[i]=new int*[ny/partition+add_solid*2];

	Solid2[0][0]=new int[nxs*(ny/partition+add_solid*2)*(nz/partition+add_solid*2)];

	
 	for (int i=1;i<ny/partition+add_solid*2;i++)
               Solid2[0][i]=Solid2[0][i-1]+nz/partition;
       
       for (int i=1;i<nxs;i++)
       {
               Solid2[i][0]=Solid2[i-1][0]+(ny/partition+add_solid*2)*(nz/partition+add_solid*2);
               for (int j=1;j<ny/partition+add_solid*2;j++)
                       Solid2[i][j]=Solid2[i][j-1]+nz/partition+add_solid*2;
       }



if (read_mod==0)
{
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
	
/*
fstream fin;
	fin.open(filename,ios::in);
	if (fin.fail())
	        {
	        cout<<"\n file open error on " << filename<<endl;
	        exit(-1);
	        }
	Solid_rank0 = new int[(NX+1)*(NY+1)*(NZ+1)];
	fin.read((char *)(&Solid_rank0[0]), sizeof(int)*(NX+1)*(NY+1)*(NZ+1));
*/

	
	
		
	
		
		
	for(k=0 ; k<nz ; k++)				///*********
	for(j=0 ; j<ny ; j++)
	for(i=0 ; i<nx ; i++)				///*********


	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore;
			
			//if (pore == 0.0)	{Solid[ci-1][cj-1][ck-1] = 0;}
			if (pore == 0.0)	{Solid[i*ny*nz+j*nz+k] = 0;}
			
			//if (pore == 1.0) 	{Solid[ci-1][cj-1][ck-1] = 1;sum++;}
			if (pore == 1.0) 	{Solid[i*ny*nz+j*nz+k] = 1;sum++;}
			
			
			
			
		}
		
	fin.close();
		
	cout<<"Geometry Reading Complete"<<endl;
	cout<<"Porosity = "<<1-(double(sum)/(nx*ny*nz))<<endl;	
	//cout<<sum<<endl;	
		
}		
else
{
	fstream fin;
	fin.open(poreFileName,ios::in);
	if (fin.fail())
	        {
	        cout<<"\n file open error on " <<poreFileName <<endl;
	        exit(-1);
	        }
	
	fin.read((char *)(&Solid[0]), sizeof(int)*nx*ny*nz);	
	cout<<"Geometry Reading Complete"<<endl;
}


	ostringstream name;
	for (int xi=0;xi<partition;xi++)
		for (int xj=0;xj<partition;xj++)
		{
		
		for (k=0;k<nz/partition+add_solid*2;k++)
			for (i=add_buf_x_n;i<nx+add_buf_x_n;i++)
			{
			Solid2[i][0][k]=1;
			Solid2[i][ny/partition+add_solid*2-1][k]=1;
			}
		
		for (j=0;j<ny/partition+add_solid*2;j++)
			for (i=add_buf_x_n;i<nx+add_buf_x_n;i++)
			{
			Solid2[i][j][0]=1;
			Solid2[i][j][nz/partition+add_solid*2-1]=1;
			}
		
		for (k=0;k<nz/partition+add_solid*2;k++)
			for (j=0;j<ny/partition+add_solid*2;j++)
			{	
			for (i=0;i<add_buf_x_n;i++)
			Solid2[i][j][k]=0;
			for (i=nx;i<nx+add_buf_x_p;i++)
			Solid2[i][j][k]=0;
			}

		for (k=0;k<nz/partition;k++)
			for (j=0;j<ny/partition;j++)
				for (i=0;i<nx;i++)
				Solid2[i+add_buf_x_n][j+add_solid][k+add_solid]=Solid[i*ny*nz+(j+partition/partition*xi)*nz+(k+partition/partition*xj)];


		
		
		name.str("");
		name<<poreFileName2<<xi*partition+xj<<".vtk";
		ofstream out;
		out.open(name.str().c_str());
		out<<"# vtk DataFile Version 2.0"<<endl;
		out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
		out<<"binary"<<endl;
		out<<"DATASET STRUCTURED_POINTS"<<endl;
		out<<"DIMENSIONS         "<<nz/partition+add_solid*2<<"         "<<ny/partition+add_solid*2<<"         "<<nxs<<endl;       ///*********
		out<<"ORIGIN 0 0 0"<<endl;
		out<<"SPACING 1 1 1"<<endl;
		out<<"POINT_DATA     "<<nxs*(ny/partition+add_solid*2)*(nz/partition+add_solid*2)<<endl;				///*********
		out<<"SCALARS sample_scalars int"<<endl;
		out<<"LOOKUP_TABLE default"<<endl;
	
	
	//cout<<Solid[2][0][0]<<"   fffff  "<<endl;
	//for(k=0+xj*nz/partition ; k<(xj+1)*nz/partition ; k++)						///*********
	//	for(j=0+xi*ny/partition ; j<(xi+1)*ny/partition ; j++)					///*********
	//		for(i=0 ; i<nx/partition ; i++)	
			out.write((char *)(&Solid2[0][0][0]), sizeof(int)*nxs*(ny/partition+add_solid*2)*(nz/partition+add_solid*2)); 

	//out.close();

	//cout<<"VTK file ouput COMPLETE"<<endl;
			//cout<<i<<" "<<j<<" "<<k<<endl;		
			//out<<Solid[i*ny*nz+j*nz+k]<<endl;
			

	//out.write((char *)(&Solid2[0]), sizeof(int)*nx1*ny1*nz1*Zoom*Zoom*Zoom); 

	out.close();
	cout<<"VTK File "<<xi*partition+xj<<" COMPLETE"<<endl;
	
		name.str("");
		name<<poreFileName2<<xi*partition+xj<<"_"<<nxs<<"_"<<ny/partition+add_solid*2<<"_"<<nz/partition+add_solid*2<<".cutdat";
		out.open(name.str().c_str());
		out.write((char *)(&Solid2[0][0][0]), sizeof(int)*nxs*(ny/partition+add_solid*2)*(nz/partition+add_solid*2));
		out.close();
		cout<<"DAT File "<<xi*partition+xj<<" COMPLETE"<<endl;
	


	
		}
	//cout<<sum<<endl;
/*
	ostringstream name;name.str("");
	name<<poreFileName<<i<<".inputdat";
	name<<poreFileName<<"_2cut.vtk";
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
	
	
	//cout<<Solid[2][0][0]<<"   fffff  "<<endl;
	for(k=0 ; k<nz ; k++)						///*********
		for(j=0 ; j<ny ; j++)					///*********
			for(i=0 ; i<nx ; i++)	
			{
			//cout<<i<<" "<<j<<" "<<k<<endl;		
			out<<Solid[i*ny*nz+j*nz+k]<<endl;
			}




	out.close();

*/

}


