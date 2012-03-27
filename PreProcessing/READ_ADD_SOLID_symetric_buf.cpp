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

int nx=87;
int ny=87;
int nz=87;
int dir=3;
int sym_x=0;
int sym_y=0;
int sym_z=0;
int add_buf_x_n=0;
int add_buf_y_n=2;
int add_buf_z_n=4;

int add_buf_x_p=0;
int add_buf_y_p=0;
int add_buf_z_p=0;

int Zoom=1; //1,2,3,4...
char poreFileName[128]="20-3-3.dat";
char poreFileNameVTK[128]="20-3-3.vtk";
char poreFileNameOut[128]="20-3-3_x.dat";
//output VTK file,0=no, 1=yes
int VTK_OUT=1;
//===========VTK AND OUT ARE ALL WRITTEN IN BINARY FORMAT===============
//===========================================================





int pls[4][3]={{0,2,2},{2,0,2},{2,2,0},{0,0,0}};

int sum=0;



int nx1,ny1,nz1;
nx1=(nx+pls[dir][0])*(sym_x+1)+add_buf_x_n+add_buf_x_p;
ny1=(ny+pls[dir][1])*(sym_y+1)+add_buf_y_n+add_buf_y_p;
nz1=(nz+pls[dir][2])*(sym_z+1)+add_buf_z_n+add_buf_z_p;


bool*** Solid_Int;
bool*** Solid;
int*** Solid_Int2;



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
	
	Solid_Int = new bool**[(nx+pls[dir][0])*(sym_x+1)];	///*********
	Solid = new bool**[nx];
	for (i=0;i<(nx+pls[dir][0])*(sym_x+1);i++)				///*********
		{
		Solid_Int[i]=new bool*[(ny+pls[dir][1])*(sym_y+1)];		///*********
		for (j=0;j<(ny+pls[dir][1])*(sym_y+1);j++)			///*********
			{
			Solid_Int[i][j]= new bool[(nz+pls[dir][2])*(sym_z+1)];///*********
			for (k=0;k<(nz+pls[dir][2])*(sym_z+1);k++)			///*********
				Solid_Int[i][j][k]= 1;
			      
			}
		}
		
	for (i=0; i<nx;i++)
	{
	       Solid[i] = new bool*[ny];
	       for (j=0;j<ny;j++)
	       {
	               Solid[i][j] = new bool[nz];
	               for (k=0;k<nz;k++)
	                       Solid[i][j][k] = 0;
	       }
	}
		
	Solid_Int2 = new int**[nx1*Zoom];	///*********
	
	for (i=0;i<nx1*Zoom;i++)				///*********
		Solid_Int2[i]=new int*[ny1*Zoom];

	Solid_Int2[0][0]=new int[nx1*ny1*nz1*Zoom*Zoom*Zoom];

	
 	for (int i=1;i<ny1*Zoom;i++)
               Solid_Int2[0][i]=Solid_Int2[0][i-1]+nz1*Zoom;
       
       for (int i=1;i<nx1*Zoom;i++)
       {
               Solid_Int2[i][0]=Solid_Int2[i-1][0]+ny1*nz1*Zoom*Zoom;
               for (int j=1;j<ny1*Zoom;j++)
                       Solid_Int2[i][j]=Solid_Int2[i][j-1]+nz1*Zoom;
       }
	





	for(k=0 ; k<nz ; k++)				///*********
	for(j=0 ; j<ny ; j++)
	for(i=0 ; i<nx ; i++)				///*********


	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore;
			
			//if (pore == 0.0)	{Solid[ci-1][cj-1][ck-1] = 0;}
			if (pore == 0.0)	{Solid[i][j][k] = 0;}
			
			//if (pore == 1.0) 	{Solid[ci-1][cj-1][ck-1] = 1;sum++;}
			if (pore == 1.0) 	{Solid[i][j][k] = 1;sum++;}
			
			
			
			
		}
		
	fin.close();
		
	
	cout<<"Porosity = "<<1-(double(sum)/(nx*ny*nz))<<endl;	
	//cout<<sum<<endl;	
		
		
	// Reading pore geometry
	for(k=pls[dir][2]/2 ; k<nz+pls[dir][2]/2 ; k++)				///*********
	for(j=pls[dir][1]/2 ; j<ny+pls[dir][1]/2 ; j++)
	for(i=pls[dir][0]/2 ; i<nx+pls[dir][0]/2 ; i++)				///*********

	{
	        Solid_Int[i][j][k]=Solid[i-pls[dir][0]/2][j-pls[dir][1]/2][k-pls[dir][2]/2];
	        
	        
	}

	if (sym_x==1)
	for(k=0 ; k<nz+pls[dir][2] ; k++)				///*********
	for(j=0 ; j<ny+pls[dir][1] ; j++)
	for(i=0 ; i<nx+pls[dir][0] ; i++)				///*********
		Solid_Int[nx+pls[dir][0]+i][j][k]=Solid_Int[nx+pls[dir][0]-1-i][j][k];


if (sym_y==1)
	for(k=0 ; k<nz+pls[dir][2] ; k++)				///*********
	for(j=0 ; j<ny+pls[dir][1] ; j++)
	for(i=0 ; i<nx+pls[dir][0] ; i++)				///*********
		Solid_Int[i][ny+pls[dir][1]+j][k]=Solid_Int[i][ny+pls[dir][1]-1-j][k];


if (sym_z==1)
	for(k=0 ; k<nz+pls[dir][2] ; k++)				///*********
	for(j=0 ; j<ny+pls[dir][1] ; j++)
	for(i=0 ; i<nx+pls[dir][0] ; i++)				///*********
		Solid_Int[i][j][nz+pls[dir][2]+k]=Solid_Int[i][j][nz+pls[dir][2]-1-k];

	
	for(k=0 ; k<(nz+pls[dir][2])*(sym_z+1)+add_buf_z_n+add_buf_z_p; k++)						
		for(j=0 ; j<(ny+pls[dir][1])*(sym_y+1)+add_buf_y_n+add_buf_y_p; j++)					
			for(i=0 ; i<(nx+pls[dir][0])*(sym_x+1)+add_buf_x_n+add_buf_x_p; i++)				
			if ((i>=add_buf_x_n) and (i<(nx+pls[dir][0])*(sym_x+1)+add_buf_x_n) and (j>=add_buf_y_n) and (j<(ny+pls[dir][1])*(sym_y+1)+add_buf_y_n) and (k>=add_buf_z_n) and (k<(nz+pls[dir][2])*(sym_z+1)+add_buf_z_n))
			for (ci=0;ci<Zoom;ci++)
				for (cj=0;cj<Zoom;cj++)
				for (ck=0;ck<Zoom;ck++) 
					Solid_Int2[i*Zoom+ci][j*Zoom+cj][k*Zoom+ck] = Solid_Int[i-add_buf_x_n][j-add_buf_y_n][k-add_buf_z_n];
			else
			for (ci=0;ci<Zoom;ci++)
				for (cj=0;cj<Zoom;cj++)
				for (ck=0;ck<Zoom;ck++) 
					Solid_Int2[i*Zoom+ci][j*Zoom+cj][k*Zoom+ck] = 0;



	if (VTK_OUT==1)
	{
	cout<<"Start writing VTK file"<<endl;
	cout<<nx1*Zoom<<"         "<<ny1*Zoom<<"         "<<nz1*Zoom<<endl;

	ostringstream name;
	name<<poreFileNameVTK;
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
	ofstream out;
	out.open(name.str().c_str());
	
	
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"binary"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<nz1*Zoom<<"         "<<ny1*Zoom<<"         "<<nx1*Zoom<<endl;       ///*********
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<nx1*ny1*nz1*Zoom*Zoom*Zoom<<endl;				///*********
	out<<"SCALARS sample_scalars int"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	
	/*
	for (k=0;k<nz1*Zoom;k++)
	{
	cout<<k<<endl;
	for (j=0;j<ny1*Zoom;j++)
	for (i=0;i<nx1*Zoom;i++)
		out<<Solid_Int2[i][j][k]<<" ";
	}
	*/

	out.write((char *)(&Solid_Int2[0][0][0]), sizeof(int)*nx1*ny1*nz1*Zoom*Zoom*Zoom); 

	out.close();

	cout<<"VTK file ouput COMPLETE"<<endl;
	}


	//=================================================================
	cout<<"Start writing DAT file"<<endl;
	cout<<nx1*Zoom<<"         "<<ny1*Zoom<<"         "<<nz1*Zoom<<endl;
	ostringstream name2;
	name2<<poreFileNameOut;
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
	ofstream out2;
	out2.open(name2.str().c_str());
	
	
	out2.write((char *)(&Solid_Int2[0][0][0]), sizeof(int)*nx1*ny1*nz1*Zoom*Zoom*Zoom); 

	out2.close();

	cout<<"DAT file ouput complete"<<endl;
	

}


