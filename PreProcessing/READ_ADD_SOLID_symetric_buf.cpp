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

int nx=350;
int ny=350;
int nz=350;
int dir=3;
int sym_x=0;
int sym_y=0;
int sym_z=0;
int add_buf_x_n=0;
int add_buf_y_n=0;
int add_buf_z_n=0;

int add_buf_x_p=0;
int add_buf_y_p=0;
int add_buf_z_p=0;

int pls[4][3]={{0,2,2},{2,0,2},{2,2,0},{0,0,0}};

int sum=0;

/*Fg = new double**[num_psi];
	for (int i=0;i<num_psi;i++)
	        Fg[i]= new double*[Count+1];
       Fg[0][0]= new double[num_psi*(Count+1)*7];
       
       for (int i=1;i<=Count;i++)
               Fg[0][i]=Fg[0][i-1]+7;
       
       for (int i=1;i<num_psi;i++)
       {
               Fg[i][0]=Fg[i-1][0]+(Count+1)*7;
               for (int j=1;j<=Count;j++)
                       Fg[i][j]=Fg[i][j-1]+7;
       }
*/

int*** Solid_Int;
bool*** Solid;

char poreFileName[128]="Guitingnb9.dat";

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
	
	Solid_Int = new int**[(nx+pls[dir][0])*(sym_x+1)];	///*********
	Solid = new bool**[nx];
	for (i=0;i<(nx+pls[dir][0])*(sym_x+1);i++)				///*********
		Solid_Int[i]=new int*[(ny+pls[dir][1])*(sym_y+1)];

	Solid_Int[0][0]=new int[(nx+pls[dir][0])*(sym_x+1)*(ny+pls[dir][1])*(sym_y+1)*(nz+pls[dir][2])*(sym_z+1)];

	
 	for (int i=1;i<(ny+pls[dir][1])*(sym_y+1);i++)
               Solid_Int[0][i]=Solid_Int[0][i-1]+(nz+pls[dir][2])*(sym_z+1);
       
       for (int i=1;i<(nx+pls[dir][0])*(sym_x+1);i++)
       {
               Solid_Int[i][0]=Solid_Int[i-1][0]+(ny+pls[dir][1])*(sym_y+1)*(nz+pls[dir][2])*(sym_z+1);
               for (int j=1;j<(ny+pls[dir][1])*(sym_y+1);j++)
                       Solid_Int[i][j]=Solid_Int[i][j-1]+(nz+pls[dir][2])*(sym_z+1);
       }


/*
	Solid_Int = new int**[(nx+pls[dir][0])*(sym_x+1)];	///*********
	Solid = new int**[nx];
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
*/

		
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




	//cout<<sum<<endl;

	ostringstream name;
	//name<<poreFileName<<"_sym.dat";
	name<<"Guitingnb9_x_sym_700x350x350_7.7.vtk";
	ofstream out;
	out.open(name.str().c_str());
	
	
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<(nx+pls[dir][0])*(sym_x+1)+add_buf_x_n+add_buf_x_p<<"         "<<(ny+pls[dir][1])*(sym_y+1)+add_buf_y_n+add_buf_y_p<<"         "<<(nz+pls[dir][2])*(sym_z+1)+add_buf_z_n+add_buf_z_p<<endl;       ///*********
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<((nx+pls[dir][0])*(sym_x+1)+add_buf_x_n+add_buf_x_p)*((ny+pls[dir][1])*(sym_y+1)+add_buf_y_n+add_buf_y_p)*((nz+pls[dir][2])*(sym_z+1)+add_buf_z_n+add_buf_z_p)<<endl;				///*********
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	


	
	for(k=0 ; k<(nz+pls[dir][2])*(sym_z+1)+add_buf_z_n+add_buf_z_p ; k++)						///*********
		for(j=0 ; j<(ny+pls[dir][1])*(sym_y+1)+add_buf_y_n+add_buf_y_p ; j++)					///*********
			for(i=0 ; i<(nx+pls[dir][0])*(sym_x+1)+add_buf_x_n+add_buf_x_p ; i++)				///*********		
			if ((i>=add_buf_x_n) and (i<(nx+pls[dir][0])*(sym_x+1)+add_buf_x_n) and (j>=add_buf_y_n) and (j<(ny+pls[dir][1])*(sym_y+1)+add_buf_y_n) and (k>=add_buf_z_n) and (k<(nz+pls[dir][2])*(sym_z+1)+add_buf_z_n))
			out<<Solid_Int[i-add_buf_x_n][j-add_buf_y_n][k-add_buf_z_n]<<endl;
			else
			out<<0.0<<endl;



	out.close();

cout<<(nx+pls[dir][0])*(sym_x+1)+add_buf_x_n+add_buf_x_p<<"         "<<(ny+pls[dir][1])*(sym_y+1)+add_buf_y_n+add_buf_y_p<<"         "<<(nz+pls[dir][2])*(sym_z+1)+add_buf_z_n+add_buf_z_p<<endl;

}


