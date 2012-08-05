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
int nxref,nyref,nzref;
int divnumori=3;
int divnum;
int oddval=0;
int evennum=0;
double oddpor;
int dir;
int dint;
int sumin,numgeonum;
int* sum_loc;


char poreFileName[128]="maxd20-3-3.dat";
char poreFileNameVTK[128]="20-3-3.vtk";
divnum=divnumori;
while (divnum%2==0)
        {evennum++;divnum=divnum/2;}
oddval=divnum;
divnum=divnumori;

cout<<evennum<<"         "<<oddval<<endl;
 


int sum=0;


int*** Solid;
double pore;


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



	Solid = new int**[nx];
	sum_loc = new int[divnum+1];
	for (int i=0;i<=divnum;i++)
	        sum_loc[i]=0;
		
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
		
	




	for(int k=0 ; k<nz ; k++)				///*********
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)				///*********


	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore;
			
			//if (pore == 0.0)	{Solid[ci-1][cj-1][ck-1] = 0;}
			if (pore == 0.0)	{Solid[i][j][k] = 0;sum++;}
			
			//if (pore == 1.0) 	{Solid[ci-1][cj-1][ck-1] = 1;sum++;}
			if (pore == 1.0) 	{Solid[i][j][k] = -1;}
			
		
			
			
		}
		
	fin.close();
		
	
	cout<<"Porosity = "<<(double(sum)/(nx*ny*nz))<<endl;	
	
	//cout<<evennum<<"         "<<oddval<<endl;
	nxref=nx;nyref=ny;nzref=nz;
	oddpor=sum/oddval;
	if (oddval>1)
	        {
	                //oddpor=sum/oddval;
	                //cout<<oddpor<<endl;
	               if (nx>ny)
	                       if (nx>nz)
	                       dir=1;
	                       else
	                               dir=3;
                       else
	                               if (ny>nz)
	                               dir=2;
	                               else
	                               dir=3;
	                   //cout<<dir<<endl;
	                if (dir==1)
	                        for (int sn=1;sn<=oddval;sn++)
	                        {
	                                sumin=0;
	                                dint=0;
	                                while ((sumin<oddpor*sn) and (dint<nx))
	                                {
	                                for (int j=0;j<ny;j++)
	                                        for (int k=0;k<nz;k++)
	                                        {
	                                                if ((Solid[dint][j][k]<sn) and (Solid[dint][j][k]>=0))
	                                                        {
	                                                                sumin++;
	                                                                if (Solid[dint][j][k]==0)
	                                                                        Solid[dint][j][k]=sn;
	                                                        }
	                                                
	                                        }
	                                        dint++;
	                                }
	                                
	                             nxref=nx/oddval;   
	                                
	                                
	                        }
	                        
	                        
	                   if (dir==2)
	                        for (int sn=1;sn<=oddval;sn++)
	                        {
	                                sumin=0;
	                                dint=0;
	                                while ((sumin<oddpor*sn) and (dint<ny))
	                                {
	                                for (int i=0;i<nx;i++)
	                                        for (int k=0;k<nz;k++)
	                                        {
	                                                if ((Solid[i][dint][k]<sn) and (Solid[i][dint][k]>=0))
	                                                        {
	                                                                sumin++;
	                                                                if (Solid[i][dint][k]==0)
	                                                                        Solid[i][dint][k]=sn;
	                                                        }
	                                                
	                                        }
	                                        dint++;
	                                }
	                                
	                                
	                              nyref=ny/oddval;     
	                                
	                        }
	                             
	                 if (dir==3)
	                        for (int sn=1;sn<=oddval;sn++)
	                        {
	                                sumin=0;
	                                dint=0;
	                                while ((sumin<oddpor*sn) and (dint<nz))
	                                {
	                                for (int i=0;i<nx;i++)
	                                        for (int j=0;j<ny;j++)
	                                        {
	                                                if ((Solid[i][j][dint]<sn) and (Solid[i][j][dint]>=0))
	                                                        {
	                                                                sumin++;
	                                                                if (Solid[i][j][dint]==0)
	                                                                        Solid[i][j][dint]=sn;
	                                                        }
	                                                
	                                        }
	                                        dint++;
	                                }
	                                
	                           nzref=nz/oddval;     
	                                
	                                
	                        }       
	                
	        }
	        else
	                {
	                  for (int k=0;k<nz;k++)
	              for (int j=0;j<ny;j++)
	              for (int i=0;i<nx;i++)
	                      if (Solid[i][j][k]==0)
	                              Solid[i][j][k]=1;
	                      else
	                              Solid[i][j][k]=-1;
	                }
	
	
	      for (int k=0;k<nz;k++)
	              for (int j=0;j<ny;j++)
	              for (int i=0;i<nx;i++)
	             if (Solid[i][j][k]<0)
	             Solid[i][j][k]=0;
	
	
	
	numgeonum=oddval;
	cout<<numgeonum<<endl;
	for (int sn=1;sn<=evennum;sn++)
	        {
	                
	                oddpor=oddpor/2;
	                for (int i=0;i<=divnum;i++)
	                        sum_loc[i]=0;
	                
	          if (nxref>nyref)
	                       if (nxref>nzref)
	                       dir=1;
	                       else
	                               dir=3;
                       else
	                               if (nyref>nzref)
	                               dir=2;
	                               else
	                               dir=3;      
	                
	             //-----------------------------
	            // cout<<dir<<"            eeeeeeeeeeeee"<<endl;
	             //cout<<Solid[88][0][0]<<"        afadsfasdf         "<<endl;
	             if (dir==1)
	                     {
	                         for (int i=0;i<nx;i++)
	                                 {
	                                         for (int j=0;j<ny;j++)
	                                         for (int k=0;k<nz;k++)
	                                         if (Solid[i][j][k]>0)
	                                                 if (sum_loc[Solid[i][j][k]]>=0)
	                                                 {sum_loc[Solid[i][j][k]]++;Solid[i][j][k]*=-1;}
	                                         
	             
	                                   
	                             
	                               for (int si=0;si<=divnum;si++)
	                                       if (sum_loc[si]>oddpor)
	                                               sum_loc[si]=-1;
	                                       
	                                 }
	                       nxref=nxref/2;          
	                     }
	                     
	            if (dir==2)
	                     {
	                         for (int j=0;j<ny;j++)
	                                 {
	                                         for (int i=0;i<nx;i++)
	                                         for (int k=0;k<nz;k++)
	                                         if (Solid[i][j][k]>0)
	                                                 if (sum_loc[Solid[i][j][k]]>=0)
	                                                 {sum_loc[Solid[i][j][k]]++;Solid[i][j][k]*=-1;}
	                                                 //cout<<Solid[i][j][k]<<"                bbbbbbbbb        "<<i<<"        "<<j<<"        "<<k<<endl;
	                                         
	                                         
	                             
	                               for (int siz=1;siz<=divnum;siz++)               
	                                       if (sum_loc[siz]>oddpor)
	                                               sum_loc[siz]=-1;
	                                 
	                                       
	                                 }
	                       nyref=nyref/2;          
	                     }         
	                     
	                if (dir==3)
	                     {
	                         for (int k=0;k<nz;k++)
	                                 {
	                                         for (int i=0;i<nx;i++)
	                                         for (int j=0;j<ny;j++)
	                                         if (Solid[i][j][k]>0)
	                                                 if (sum_loc[Solid[i][j][k]]>=0)
	                                                 {sum_loc[Solid[i][j][k]]++;Solid[i][j][k]*=-1;}
	                                
	                             
	                               for (int si=0;si<=divnum;si++)
	                                       if (sum_loc[si]>oddpor)
	                                               sum_loc[si]=-1;
	                                       
	                                 }
	                      nzref=nzref/2;           
	                     }              
	         numgeonum*=2;
	         cout<<numgeonum<<endl;
	         for (int k=0;k<nz;k++)
	              for (int j=0;j<ny;j++)
	              for (int i=0;i<nx;i++)
	              {
	                      if (Solid[i][j][k]>0)
	                                   Solid[i][j][k]=Solid[i][j][k]*2;
	                           if (Solid[i][j][k]<0)
	                                          Solid[i][j][k]=-Solid[i][j][k]*2-1;
	                
	                
	                
	        }
	
	  }
	    
	
	//============decomposition complete=======================
	
	//==================================================
	
	
	
	
	
	
	
	
	
	cout<<"Start writing VTK file"<<endl;
	

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
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	
	
	for (int k=0;k<nz;k++)
	{
	//cout<<k<<endl;
	for (int j=0;j<ny;j++)
	for (int i=0;i<nx;i++)
		out<<Solid[i][j][k]<<" ";
	}
	

	//out.write((char *)(&Solid_Int2[0][0][0]), sizeof(int)*nx1*ny1*nz1*Zoom*Zoom*Zoom); 

	out.close();

	cout<<"VTK file ouput COMPLETE"<<endl;
	
}


