/*
 * Copyright 1997, Regents of the University of Minnesota
 *
 * main.c
 *
 * This is the entry point of the ILUT
 *
 * Started 10/19/95
 * George
 *
 * $Id: parmetis.c,v 1.5 2003/07/30 21:18:54 karypis Exp $
 *
 */

//#include <parmetisbin.h>sstream



#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<math.h>

#include<parmetis.h>

using namespace std; 


/*************************************************************************
* Let the game begin
**************************************************************************/
int main(int argc, char *argv[])
{


        
MPI_Init(&argc, &argv);
//=================
 MPI_Comm comm;
 //comm=MPI_COMM_WORLD;

	MPI_Comm_dup(MPI_COMM_WORLD, &comm);
 
 //=================
 
int rank = MPI :: COMM_WORLD . Get_rank ();
int para_size=MPI :: COMM_WORLD . Get_size ();


int NCHAR=128;
	char     filename[128], dummy[128+1];
	int      dummyInt;



int sym_x;
int sym_y;
int sym_z;
int add_buf_x_n;
int add_buf_y_n;
int add_buf_z_n;

int add_buf_x_p;
int add_buf_y_p;
int add_buf_z_p;
int add_porous_plate; //0=OFF, 1=fine plate,pore size1, pore size2, 3=posr size3
int porous_position; //-1=defualt position,end of the geometry, or give a positive value
int Zoom; //1,2,3,4...


int nx,ny,nz;
int expvtk,expdat,bindat,fil,geo_mod;
int dir;

int sum_rec;


const int e[18][3]=
{{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};

char poreFileName[128];
char poreFileNameMET[128];
char poreFileNameVTK[128];
char poreFileNameOut[128];
int mark,ii,jj,kk,loop,sum2;

int sum=0;
int sum3;
int mesh_par;
int partition_vtk;
int reci,recj,reck;

if (rank==0)
{
ifstream fins(argv[1]);
			
							fins.getline(dummy, NCHAR);
	fins >> poreFileName;				fins.getline(dummy, NCHAR);
	fins >> nx>> ny>>nz;				fins.getline(dummy, NCHAR);
	fins >> poreFileNameVTK;			fins.getline(dummy, NCHAR);
	fins >> poreFileNameOut;			fins.getline(dummy, NCHAR);
	fins >> expvtk;					fins.getline(dummy, NCHAR);
	fins >> expdat;					fins.getline(dummy, NCHAR);
	fins >> bindat;					fins.getline(dummy, NCHAR);
							fins.getline(dummy, NCHAR);
	fins >> fil;					fins.getline(dummy, NCHAR);
							fins.getline(dummy, NCHAR);
	fins >> geo_mod;					fins.getline(dummy, NCHAR);
	fins >> dir;					fins.getline(dummy, NCHAR);
	fins >> sym_x >> sym_y >> sym_z;			fins.getline(dummy, NCHAR);
	fins >> add_buf_x_n>>add_buf_y_n>>add_buf_z_n;	fins.getline(dummy, NCHAR);
	fins >> add_buf_x_p>>add_buf_y_p>>add_buf_z_p;	fins.getline(dummy, NCHAR);
	fins >> add_porous_plate;			fins.getline(dummy, NCHAR);
	fins >> porous_position;				fins.getline(dummy, NCHAR);
	fins >> Zoom;					fins.getline(dummy, NCHAR);
							fins.getline(dummy, NCHAR);
	fins >> mesh_par;				fins.getline(dummy, NCHAR);
	fins >> poreFileNameMET;				fins.getline(dummy, NCHAR);
	fins >> partition_vtk;				fins.getline(dummy, NCHAR);



fins.close();	

}


       MPI_Bcast(&poreFileName,128,MPI_CHAR,0,MPI_COMM_WORLD);MPI_Bcast(&poreFileNameVTK,128,MPI_CHAR,0,MPI_COMM_WORLD); 
       MPI_Bcast(&poreFileNameOut,128,MPI_CHAR,0,MPI_COMM_WORLD);
       
       MPI_Bcast(&nx,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&ny,1,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(&nz,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&expvtk,1,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(&expdat,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&bindat,1,MPI_INT,0,MPI_COMM_WORLD);
       MPI_Bcast(&fil,1,MPI_INT,0,MPI_COMM_WORLD);MPI_Bcast(&mesh_par,1,MPI_INT,0,MPI_COMM_WORLD);

       
       
int*** Solid;
double pore;
FILE *ftest;
	ifstream fin;
	
if (rank==0)
{
	
	
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
}

	
	Solid = new int**[nx];	///*********
	
	for (int i=0;i<nx;i++)				///*********
		Solid[i]=new int*[ny];

	Solid[0][0]=new int[nx*ny*nz];

	
 	for (int i=1;i<ny;i++)
               Solid[0][i]=Solid[0][i-1]+nz;
       
       for (int i=1;i<nx;i++)
       {
               Solid[i][0]=Solid[i-1][0]+ny*nz;
               for (int j=1;j<ny;j++)
                       Solid[i][j]=Solid[i][j-1]+nz;
       }


	/*
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

	*/


	if (rank==0)
	{
	cout<<"Start reading source geometry file"<<endl;	
	for(int k=0 ; k<nz ; k++)				///*********
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)				///*********
	

	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore;
			
			//if (pore == 0.0)	{Solid[ci-1][cj-1][ck-1] = 0;}
			if (pore == 0)	{sum++;Solid[i][j][k] = sum-1;}
			else
			//if (pore == 1.0) 	{Solid[ci-1][cj-1][ck-1] = 1;sum++;}
			//if (pore == 1.0) 	
				{Solid[i][j][k] = -1;}
			
		
			
			
		}
	//cout<<"Porosity = "<<(double(sum)/(nx*ny*nz))<<endl;	
	sum_rec=sum;
	fin.close();

        }
        
        MPI_Bcast(Solid[0][0],nx*ny*nz,MPI_INT,0,MPI_COMM_WORLD);
        MPI_Bcast(&sum_rec,1,MPI_INT,0,MPI_COMM_WORLD);
        
        
        idx_t* vtxdist=NULL;
	idx_t *xadj, *adjncy;
        vtxdist = new idx_t[para_size+1];
	//vtxdist = (idx_t*)malloc(sizeof(idx_t)*(para_size+1));
	//vtxdist = imalloc(para_size,"main: vtxdist");

        //cout<<sizeof(idx_t)<<"		"<<sizeof(int)<<"	******"<<endl;


        int ls1,ls2;
        int *proc_size;
        proc_size = new int[para_size];
        
        ls1=sum_rec%para_size;
        ls2=(int)sum_rec/para_size;
        
        for (int i=0;i<para_size;i++)
                proc_size[i]=ls2;
        for(int i=0;i<ls1;i++)
                proc_size[i]++;
        
        //if (rank==0)
        //        for (int i=0;i<para_size;i++)
        //        cout<<proc_size[i]<<"		"<<rank<<endl;
        


        
        vtxdist[0]=0;
        for (int i=1;i<=para_size;i++)
        {        vtxdist[i]=vtxdist[i-1]+proc_size[i-1];
                //if (rank==4)
                //        cout<<vtxdist[i]<<"	"<<vtxdist[i-1]<<"	"<<i<<"		"<<proc_size[i-1]<<endl;
        }
        
        
        
        xadj = new idx_t[proc_size[rank]+1];
        sum=0;
        
        
        for(int k=0 ; k<nz ; k++)			
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)	
	{
	        if ((Solid[i][j][k]>=vtxdist[rank]) and (Solid[i][j][k]<vtxdist[rank+1]))
	                for (int ls=0;ls<18;ls++)
		{
		ii=i+e[ls][0];
		jj=j+e[ls][1];
		kk=k+e[ls][2];	
		
		
		//================
		/*
		if (ii>=nx) ii=0;
		if (ii<0) ii=nx-1;
		if (jj>=ny) jj=0;
		if (jj<0) jj=ny-1;
		if (kk>=nz) kk=0;
		if (kk<0) kk=nz-1;
		*/
		//================
		
		if ((ii>=0) and (ii<nx) and (jj>=0) and (jj<ny) and (kk>=0) and (kk<nz) and (Solid[ii][jj][kk]>=0))
			sum++;
		}
		
	}
        
      adjncy = new idx_t[sum];
      xadj[0]=0;
      
      sum=0;sum3=0;
      for(int k=0 ; k<nz ; k++)			
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)	
	{
	        if ((Solid[i][j][k]>=vtxdist[rank]) and (Solid[i][j][k]<vtxdist[rank+1]))
	        {
	                for (int ls=0;ls<18;ls++)
	                {
		ii=i+e[ls][0];
		jj=j+e[ls][1];
		kk=k+e[ls][2];	
		
		
		//================
		/*
		if (ii>=nx) ii=0;
		if (ii<0) ii=nx-1;
		if (jj>=ny) jj=0;
		if (jj<0) jj=ny-1;
		if (kk>=nz) kk=0;
		if (kk<0) kk=nz-1;
		*/
		//================
		
		if ((ii>=0) and (ii<nx) and (jj>=0) and (jj<ny) and (kk>=0) and (kk<nz) and (Solid[ii][jj][kk]>=0))
		        {adjncy[sum]=Solid[ii][jj][kk];sum++;}
		        }
		        
		        
		        
		sum3++;
		xadj[sum3]=sum;
		
		
		}
		
	}


	//if (rank==2)
	//	for (int i=0;i<proc_size[rank]+1;i++)
	//		cout<<xadj[i]<<"	"<<adjncy[i]<<"		"<<vtxdist[rank+1]<<endl;

	//cout<<"@@@@@@@@@@@@	"<<endl;

	

	idx_t *vwgt=NULL;
	
	idx_t *adjwgt=NULL;
	idx_t wgtflag=0;
	idx_t numflag=0;
	idx_t ncon=20;
	idx_t nparts=mesh_par;
	real_t *tpwgts;
	real_t *ubvec;
	idx_t options[10];
	idx_t edgecut;
	idx_t *part;
	part = new idx_t[proc_size[rank]];
	//vwgt = new idx_t[proc_size[rank]*ncon];

	
	tpwgts = new real_t[ncon*nparts];
	for (int i=0;i<ncon*nparts;i++)
	        tpwgts[i]=1.0/(real_t)nparts;
	
	ubvec = new real_t[ncon];
	for (int i=0;i<ncon;i++)
	        ubvec[i]=1.05;
	
	//vwgt = new idx_t[proc_size[rank]];
	//	for (int i=0;i<proc_size[rank];i++)
	//	vwgt[i]=1;

	//adjwgt = new idx_t[sum];
	//	for (int i=0;i<sum;i++)
	//	adjwgt[i]=1;

	//cout<<vtxdist[0]<<" "<<vtxdist[1]<<" "<<vtxdist[2]<<"	"<<rank<<"	"<<&vtxdist[0]<<"   "<<&vtxdist[1]<<"	"<<&vtxdist[2]<<" "<<&vtxdist[3]<<endl;
	//cout<<xadj[0]<<" "<<xadj[1]<<" "<<xadj[2]<<"	"<<rank<<endl;
	options[0] = 0;
	
	MPI_Barrier(MPI_COMM_WORLD);
	
	
	
      ParMETIS_V3_PartKway(vtxdist, xadj, adjncy, vwgt,adjwgt, &wgtflag, &numflag, &ncon, &nparts, tpwgts, ubvec, options, &edgecut, part, &comm);
	
	
        MPI :: Finalize ();
        
        
    
}
