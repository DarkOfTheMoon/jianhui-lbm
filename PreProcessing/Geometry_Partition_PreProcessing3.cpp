#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include <math.h>
# include "mpi.h"

using namespace std; 

int main (int argc , char * argv [])
{
MPI :: Init (argc , argv );



int nx=90;
int ny=90;
int nz=90;
int ii,jj,kk,pore2;
int ip,jp,kp;
//============
int procn=MPI :: COMM_WORLD . Get_size ();           //total processor number
int* sumss;             //fluid nodes number of every partition            [procn+1]
                                
int procind=MPI :: COMM_WORLD . Get_rank ()+1;                        //current processor index, start from 1
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
int** nei_loc;                           //index of 18 neibourghs, nei_loc[3][0] is the first neighbour (e[1][]) of node 3

//-------------------
int* coor;      //start from 1, int [sumss[procind]+1];
//------------------

double** bufsend;                       //send buffet for f functions bufsend[comm_index][number]
double** bufrecv;                       //recv buffet for f functions
//-----------------
int** buflocsend;                       //exchange commu info, used to locate the data after when it is received from MPI communications
int** buflocrecv;                       //two digital combination index*19+ls ls is the direction of 18 vectors
//----------------

int proc_com[procn+1];                  //index convert proc index---->commu index in current processor
for (int i=0;i<=procn;i++)
        proc_com[i]=0;

//-------------------
int* sumtmp;
//-------------------
//============

//Geo_par.dat file
//0=solid, 1,2,3.... indicates different partitions

char poreFileName[128]="geo_par.dat";
char poreFileName2[128]="20-3-3.graph.part.10";
char poreFileNameVTK[128]="20-3-3.vtk";

const int e[19][3]=
{{0,0,0},{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};

const int LR[19]={0,2,1,4,3,6,5,10,9,8,7,14,13,12,11,18,17,16,15};

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
	                                        ii=nx-1;
	                                if (ii>=nx)
	                                        ii=0;
	                                
	                                        
	                        jj=j+e[ls][1];
	                                if (jj<0)
	                                        jj=ny-1;
	                                if (jj>=ny)
	                                        jj=0;
	                                
	                        kk=k+e[ls][2];
	                                       if (kk<0)
	                                        kk=nz-1;
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
	                                coor[Solid2[i][j][k]]=i*ny*nz+j*nz+k;
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
	                {com_ind[tmpint]=i;proc_com[i]=tmpint;tmpint++;com_memo_sum+=bufinfo[i];}
	                
	//cout<<com_n<<endl;
	//for (int i=0;i<com_n;i++)
	//        cout<<com_ind[i]<<endl;
	
	
	
	
	
	
	bufsend = new double* [com_n];
	        for (int i=0;i<com_n;i++)
	                bufsend[i] = new double[bufinfo[com_ind[i]]];
	        
	bufrecv = new double* [com_n];
	        for (int i=0;i<com_n;i++)
	                bufrecv[i] = new double[bufinfo[com_ind[i]]];
	        
	        
	nei_loc= new int*[sumss[procind]+1];
	nei_loc[0] =new int[(sumss[procind]+1)*19];
	        for (int i=1;i<=sumss[procind];i++)
	                nei_loc[i] = nei_loc[i-1]+19;
	  
	        buflocsend = new int*[com_n];
	        for (int i=0;i<com_n;i++)
	                buflocsend[i] = new int[bufinfo[com_ind[i]]];
	        
	        buflocrecv = new int*[com_n];
	        for (int i=0;i<com_n;i++)
	                buflocrecv[i] = new int[bufinfo[com_ind[i]]];
	        
	        sumtmp = new int[com_n];
	                for(int i=0;i<com_n;i++)
	                        sumtmp[i]=0;
	                
	                
	        //tmpint=0;
	        for (int ci=1;ci<=sumss[procind];ci++)
	                {
	                      ii=(int)(coor[ci]/(ny*nz));
	                      jj=(int)((coor[ci]%(ny*nz))/nz);
	                      kk=(int)(coor[ci]%nz);
	                      
	                      for (int mi=0; mi<19; mi++)
			{
			
			
			        ip=ii+e[mi][0];if (ip<0) {ip=nx-1;};if (ip>=nx) {ip=0;};
			        jp=jj+e[mi][1];if (jp<0) {jp=ny-1;}; if (jp>=ny) {jp=0;};
			        kp=kk+e[mi][2];if (kp<0) {kp=nz-1;}; if (kp>=nz) {kp=0;};
			        
			        if (Solid[ip][jp][kp]==procind)
			                nei_loc[ci][mi]=Solid2[ip][jp][kp];
			        else
			                if (Solid[ip][jp][kp]==0)
			                        nei_loc[ci][mi]=0;
			                else
			                {
			                        nei_loc[ci][mi]=-Solid[ip][jp][kp];
			                        //cout<<nei_loc[ci][mi]<<endl;
			                        buflocsend[proc_com[Solid[ip][jp][kp]]][sumtmp[proc_com[Solid[ip][jp][kp]]]]=Solid2[ip][jp][kp]*19+mi;
			                        //cout
			                        //<<buflocsend[proc_com[Solid[ip][jp][kp]]][sumtmp[proc_com[Solid[ip][jp][kp]]]]
			                        //<<endl;
			                        
			                        sumtmp[proc_com[Solid[ip][jp][kp]]]++;
			                        //if (Solid[ip][jp][kp]==6)        tmpint++;
			                }
			                
			}
			                
			                
	                    
	                }
	                
	                
	                
	                //for (int i=0;i<com_n;i++)
	                //        cout<<sumtmp[i]<<"        "<<bufinfo[com_ind[i]]<<"        "<<procind<<"@"<<i<<endl;
	                
	        //for (int i=0;i<com_n;i++)
	           //                   for (int j=0;j<bufinfo[com_ind[i]];j++)
	              //                       cout<<buflocsend[i][j]<<"        "<<i<<"        "<<j<<"        "<<procind<<endl;
	                                               
	                                               
	                                               
	       MPI_Status status[com_n*2] ;
	       MPI_Request request[com_n*2];         
	       int mpi_test=procn;
	       
	       for (int i=0;i<com_n;i++)
	       
	               {
	                       MPI_Isend(buflocsend[i],bufinfo[com_ind[i]], MPI_INT, com_ind[i]-1, (procind-1)*procn+com_ind[i]-1, MPI_COMM_WORLD,&request[2*i]);
	                       
	                       MPI_Irecv(buflocrecv[i],bufinfo[com_ind[i]], MPI_INT, com_ind[i]-1, (com_ind[i]-1)*procn+procind-1, MPI_COMM_WORLD,&request[2*i+1]);		
	               }
      		
      		MPI_Waitall(2*com_n,request, status);
      		MPI_Testall(2*com_n,request,&mpi_test,status);

      		
      		int* testarr;
      		int testl1,testl2;
      		testarr=new int[(sumss[procind]+1)*19];
      		       
      		for (int ci=0;ci<(sumss[procind]+1)*19;ci++)
	               testarr[ci]=0;
	       
	       for (int ci=1;ci<=sumss[procind];ci++)
	               {
	                    for (int mi=0; mi<19; mi++)
	                            if (nei_loc[ci][mi]>0)
	                                    testarr[nei_loc[ci][mi]*19+mi]=1;
	                            else
	                                   if (nei_loc[ci][mi]==0)
	                                            testarr[ci*19+LR[mi]]=1;
	                            
	                                    
	                                    
	                                    
	               }
	               
	               
	               for (int i=0;i<com_n;i++)
	                       {
	                               for (int j=0;j<bufinfo[com_ind[i]];j++)
	                                       {
	                                               testl1=(int)(buflocrecv[i][j]/19);
	                                               testl2=(int)(buflocrecv[i][j]%19);
	                                               //cout<<testarr[testl1*19+testl2]<<"  "<<procind<<"   "<<testl1<<"        "<<testl2<<"        "<<buflocrecv[i][j]<<endl;
	                                               //cout<<buflocsend[i][j]<<"  ";
	                                               testarr[testl1*19+testl2]=1;
	                                       }
	                       }
	                      
	                       
	                      
	                       for (int ci=19;ci<(sumss[procind]+1)*19;ci++)
	                       {        //cout<<testarr[ci]<<endl;
	                               if (testarr[ci]<1)
	                                       cout<<ci<<"        "<<procind<<endl;
	                       }

	                       
	        
	      //cout<<endl;
	      //for (int i=0;i<com_n;i++)
	      //         cout<<bufinfo[com_ind[i]]<<"                "<<sumtmp[i]<<"                "<<com_ind[i]<<endl;
	        
	MPI :: Finalize ();
}


