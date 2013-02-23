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

//============
int procn=10;           //total processor number
int* sumss;             //fluid nodes number of every partition            [procn+1]
                                
int procind=6+1;                        //current processor index, start from 1
int neib=0;
int* bufinfo;                                   //number of nodes of the partirion that need communiate with current processors
                                                        // 0= no contact with current processor, >0 number of nodes that need to communicate with current processor.
                                                         //start from 1 size procn+1       designed for 19 components for f function transfer
                                                         
                                                         
                                                         
                                                        
        bufinfo=new int[procn+1];
for (int i=0;i<=procn;i++)
        bufinfo[i]=0;                           

int com_n=0;                            //mpi commu numbers     number of neighbour partitions which need communication
int com_memo_sum=0;             //mpi commu buffets length, length of buffet array 
int* com_ind;                          //commu nodes indexs (partition no.)   com_ind[0,start from 0] size new int[com_n]
int* com_loc;                           //mpi commu different nodes starting locations in buffet  arrays
int tmpint;

//-------------------
int* coor;      //start from 1, int [sumss[procind]+1];
//------------------

double** bufsend;
double** bufrecv;

//============

//Geo_par.dat file
//0=solid, 1,2,3.... indicates different partitions




char poreFileName[128]="geo_par.dat";
char poreFileName2[128]="20-3-3.graph.part.10";
char poreFileNameVTK[128]="20-3-3.vtk";

const int e[19][3]=
{{0,0,0},{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};


int sum=0;
int sum2=0;

int*** Solid;
int*** Solid2;

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
	Solid2 = new int**[nx];
	
		
	for (int i=0; i<nx;i++)
	{
	       Solid[i] = new int*[ny];
	       Solid2[i] = new int*[ny];
	       for (int j=0;j<ny;j++)
	       {
	               Solid[i][j] = new int[nz];
	                Solid2[i][j] = new int[nz];
	               for (int k=0;k<nz;k++)
	                       Solid[i][j][k] = 0,Solid2[i][j][k] = 0;
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
			
			Solid[i][j][k]=pore;
			
		
			
			
		}
		
	fin.close();

	sumss=new int [procn+1];
	for (int i=0;i<=procn;i++)
	        sumss[i]=0;
	
	for(int k=0 ; k<nz ; k++)			
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)			
	        if (Solid[i][j][k]>0)
	        {
	                sumss[Solid[i][j][k]]++;
	                Solid2[i][j][k]=sumss[Solid[i][j][k]];
	                
	                
	                //=======calculate neibough numbers==========
	                if (Solid[i][j][k]==procind)
	                for (int ls=1;ls<19;ls++)
	                {
	                        ii=i+e[ls][0];
	                                if (ii<0)
	                                        ii=nx;
	                                if (ii>=nx)
	                                        ii=0;
	                                
	                                        
	                        jj=j+e[ls][1];
	                                if (jj<0)
	                                        jj=ny;
	                                if (jj>=ny)
	                                        jj=0;
	                                
	                        kk=k+e[ls][2];
	                                       if (kk<0)
	                                        kk=nz;
	                                        if (kk>=nz)
	                                        kk=0; 
	                        
	                                if ((Solid[ii][jj][kk]>0) and (Solid[ii][jj][kk]!=procind))
	                                        {
	                                                bufinfo[Solid[ii][jj][kk]]++;
	                                        //cout<<procind<<"        "<<Solid[ii][jj][kk]<<endl;
	                                        }
	                                
	                                        
	                
	                }
	        //=======================================================
	        }
	        else
	                Solid2[i][j][k]=0;
	       
	        //======coordinate of nodes===========
	        coor = new int[sumss[procind]+1];
	                for(int k=0 ; k<nz ; k++)			
	                for(int j=0 ; j<ny ; j++)
	                for(int i=0 ; i<nx ; i++)
	                        if (Solid[i][j][k]==procind)
	                                coor[Solid2[i][j][k]]=i*ny*nz+j*ny+k;
	        //=============================
	        
	        
	        for (int i=1;i<=procn;i++)
	        {
	                //cout<<bufinfo[i]<<endl;
	                if (bufinfo[i]>0)
	                        com_n++;
	        }
	        com_ind=new int[com_n];
	                
	        tmpint=0;
	        for (int i=1;i<=procn;i++)
	                if (bufinfo[i]>0)
	                {com_ind[tmpint]=i;tmpint++;com_memo_sum+=bufinfo[i];}
	                
	cout<<com_n<<endl;
	for (int i=0;i<com_n;i++)
	        cout<<com_ind[i]<<endl;
	
	bufsend = new double* [com_n];
	        for (int i=0;i<com_n;i++)
	                bufsend[i] = new double[bufinfo[com_ind[i]]];
	        
	bufrecv = new double* [com_n];
	        for (int i=0;i<com_n;i++)
	                bufrecv[i] = new double[bufinfo[com_ind[i]]];
	        
	        
	        
	        
	
	
}


