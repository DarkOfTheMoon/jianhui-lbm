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

int nx_read=640;
int ny_read=320;
int nz_read=320;

int nx=320;
int ny=320;
int nz=320;

//=====source data option======
int source_data_opt=1;	//1=OWN,2=Keehm
double bodyf=2e-6;
double vis=0.05;
//=============================

double dx=0.009;	//resolution mm
const int sub_n=6;

int sub_num[sub_n]={4,3,3,2,2,2};
int sub_size[sub_n]={80,120,160,180,200,230};
int sub_size2[sub_n]={80,120,160,180,200,230};
int sub_size3[sub_n]={80,120,160,180,200,230};


	//-------------AUTO MODEL FOR SUBSIZE Y AND Z-------------------
	for (int i=0;i<sub_n;i++)
		{
		sub_size2[i]=(int)sub_size[i]*ny/nx;cout<<sub_size[i]<<"	"<<sub_size2[i]<<endl;
		sub_size3[i]=(int)sub_size[i]*nz/nx;cout<<sub_size[i]<<"	"<<sub_size2[i]<<endl;
		}
	
	//--------------------------------------------------------------


int input_vtk=1;	//0=NO,1=YES
int pgDir=1;
	
char poreFileName[128]="port_perm_LBM_velocity_Vector_120000.vtk";	//velocity
char poreFileName_geo[128]="port_perm_LBM_Geometry.vtk";		//geometry


char ouput_prefix[128]="port_perm_";


//======================================================

double pgBB[3];
double nu=1.0/6.0;
double	pgMag		 = 1.0/6.0*1e-3;	

	pgBB[0] = pgBB[1] = pgBB[2] = 0.0;

	if     (pgDir == 3) pgBB[2] = pgMag/2.0;
	else if(pgDir == 2) pgBB[1] = pgMag/2.0;
	else                pgBB[0] = pgMag/2.0;


int NCHAR=128;
char     dummy[128+1];

int inter_sub[sub_n];



for (int i=0;i<sub_n;i++)
	{
	if (sub_num[i]>1)
	inter_sub[i]=(nx-sub_size[i])/(sub_num[i]-1);
	else
		inter_sub[i]=0;
	//cout<<i<<"	"<<inter_sub[i]<<endl;
	}



double**** vel;
bool*** Solid;


	FILE *ftest;
	ifstream fin;
	
	/*
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
*/

double pore1,pore2,pore3;
	int i, j, k;
	
double pore;
	
	
	
	Solid = new bool**[nx];
	
		
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


	vel = new double***[nx];
		for (int i=0;i<nx;i++)
		{
		vel[i] = new double**[ny];
			for (int j=0;j<ny;j++)
			{
			vel[i][j] = new double*[nz];
				for (int k=0;k<nz;k++)
				{
				vel[i][j][k] = new double[3];
				vel[i][j][k][0]=0.0;
				vel[i][j][k][1]=0.0;
				vel[i][j][k][2]=0.0;
				}
			}
		}



	ftest = fopen(poreFileName_geo, "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file (" << poreFileName <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest);

	fin.open(poreFileName_geo);
	
	cout<<"Start Reading Geometry File"<<endl;

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

	for(k=0 ; k<nz_read ; k++)				///*********
	for(j=0 ; j<ny_read ; j++)
	for(i=0 ; i<nx_read ; i++)				///*********

	 
	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore;
			if ((i<nx) and (j<ny) and (k<nz))
			if (pore==1.0)
				{Solid[i][j][k]=1;}
			else
				{Solid[i][j][k]=0;}
			
			
			
		}
		
	fin.close();

	cout<<endl;	
	cout<<"Geometry File Reading Complete"<<endl;
	cout<<endl;




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

	cout<<"Start Reading Velocity File	"<<endl;
		

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

	for(k=0 ; k<nz_read ; k++)				///*********
	for(j=0 ; j<ny_read ; j++)
	for(i=0 ; i<nx_read ; i++)				///*********

	 
	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore1>>pore2>>pore3;
			if ((i<nx) and (j<ny) and (k<nz))
			if (Solid[i][j][k]==0) 
			{
			vel[i][j][k][0]=pore1;
			vel[i][j][k][1]=pore2;
			vel[i][j][k][2]=pore3;
			}
			
			
			
			
		}
		
	fin.close();


	cout<<endl;	
	cout<<"Velocity File Reading Complete"<<endl;
	cout<<endl;



int st_i,st_j,st_k;
int sum=0;
double vx,vy,vz;
double permX,permY,permZ;
double factor;
/*
	const int sub_n=4;

int sub_num[sub_n]={4,3,2,2};
int sub_size[sub_n]={80,100,160,200};

inter_sub

*/
	ostringstream name;
	ofstream out;


	for (int pri_ind=0;pri_ind<sub_n;pri_ind++)
	{
				name.str("");
				name<<ouput_prefix<<"Sub_Sample_Perm_"<<sub_size[pri_ind]<<"x"<<sub_size2[pri_ind]<<"x"<<sub_size3[pri_ind]<<".sta_dat";
				out.open(name.str().c_str());
		for (int li=0;li<sub_num[pri_ind];li++)
			for (int lj=0;lj<sub_num[pri_ind];lj++)
				for (int lk=0;lk<sub_num[pri_ind];lk++)
				{
				sum=0;vx=0.0;vy=0.0;vz=0.0;
				st_i=li*inter_sub[pri_ind];
				st_j=lj*inter_sub[pri_ind];
				st_k=lk*inter_sub[pri_ind];
			
				
				for (int i=st_i;i<st_i+sub_size[pri_ind];i++)
					for (int j=st_j;j<st_j+sub_size2[pri_ind];j++)
						for (int k=st_k;k<st_k+sub_size3[pri_ind];k++)
						{
						
						if (Solid[i][j][k]==0)
							sum++;
						//permX = 1e9*(Q[0]/(nx*ny*nz)+factor*pgBB[0])*nu*dx*dx/pgMag
						vx+=vel[i][j][k][0];
						vy+=vel[i][j][k][1];
						vz+=vel[i][j][k][2];


						}
				factor=(double)sum/(sub_size[pri_ind]*sub_size2[pri_ind]*sub_size3[pri_ind]);
				//cout<<sum<<"	"<<sub_size[pri_ind]*sub_size[pri_ind]*sub_size[pri_ind]<<endl;
		if (source_data_opt==2)
		{
		permX = 1e9*(vx/(sub_size[pri_ind]*sub_size2[pri_ind]*sub_size3[pri_ind])+factor*pgBB[0])*nu*dx*dx/pgMag;
		permY = 1e9*(vy/(sub_size[pri_ind]*sub_size2[pri_ind]*sub_size3[pri_ind])+factor*pgBB[1])*nu*dx*dx/pgMag;  
		permZ = 1e9*(vz/(sub_size[pri_ind]*sub_size2[pri_ind]*sub_size3[pri_ind])+factor*pgBB[2])*nu*dx*dx/pgMag;
		}
		else
		{
		permX=vx/(sub_size[pri_ind]*sub_size2[pri_ind]*sub_size3[pri_ind])*vis/bodyf*dx*dx*1e9;
		permY=vy/(sub_size[pri_ind]*sub_size2[pri_ind]*sub_size3[pri_ind])*vis/bodyf*dx*dx*1e9;
		permZ=vz/(sub_size[pri_ind]*sub_size2[pri_ind]*sub_size3[pri_ind])*vis/bodyf*dx*dx*1e9;
		}
				
				
				if     (pgDir == 3) out<<factor<<" "<<permZ<<endl;
				else if(pgDir == 2) out<<factor<<" "<<permY<<endl;
				else                out<<factor<<" "<<permX<<endl;
							//out<<factor<<" "<<permX<<" "<<permY<<" "<<permZ<<endl;
				}
			out.close();

	}

}


