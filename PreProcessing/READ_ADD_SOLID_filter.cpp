#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>

using namespace std; 
const int e[7][3]=
{{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
 const int e1[6][3]=
 {{0,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
 const int e2[6][3]=
 {{0,0,0},{1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
 const int e3[6][3]=
 {{0,0,0},{1,0,0},{-1,0,0},{0,-1,0},{0,0,1},{0,0,-1}};
 const int e4[6][3]=
 {{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,0,1},{0,0,-1}};
const int e5[6][3]=
 {{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,-1}};
const int e6[6][3]=
 {{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1}};
 
int main (int argc , char * argv [])
{

int nx=300;
int ny=300;
int nz=300;
int dir=3; 	//0,1,2 x,y,z add extra solid boundaries, 
		//3,no BCs, 4,5,6: solid the most outside to form a BC,flow dirc:x,y,z
int sym_x=1;
int sym_y=0;
int sym_z=0;

int add_buf_x_n=0;
int add_buf_y_n=0;
int add_buf_z_n=0;

int add_buf_x_p=0;
int add_buf_y_p=0;
int add_buf_z_p=0;
int add_porous_plate=0;
int porous_position=206; //-1=defualt position,end of the geometry, or give a positive value
int Zoom=1; //1,2,3,4...
int Filter=1;
char poreFileName[128]="MtGambier_nb5.dat_cut.dat";
char poreFileNameVTK[128]="MtGambier_nb5.vtk";
char poreFileNameOut[128]="MtGambier_f600x300x300_9.dat";
char Non_conective_poreFileNameOut[128]="check.vtk";
//output VTK file,0=no, 1=yes
int VTK_OUT=0;
int vtk_check=0;
//===========VTK AND OUT ARE ALL WRITTEN IN BINARY FORMAT===============
//===========================================================





int pls[7][3]={{0,2,2},{2,0,2},{2,2,0},{0,0,0},{0,0,0},{0,0,0}};

int sum=0;
int sum2;
int si,sj,sk;


int nx1,ny1,nz1;
nx1=(nx+pls[dir][0])*(sym_x+1)+add_buf_x_n+add_buf_x_p;
ny1=(ny+pls[dir][1])*(sym_y+1)+add_buf_y_n+add_buf_y_p;
nz1=(nz+pls[dir][2])*(sym_z+1)+add_buf_z_n+add_buf_z_p;


int*** Solid_Int;
int*** Solid;
int*** Solid2;
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
	
	Solid_Int = new int**[(nx+pls[dir][0])*(sym_x+1)];	///*********
	Solid = new int**[nx];
	for (i=0;i<(nx+pls[dir][0])*(sym_x+1);i++)				///*********
		{
		Solid_Int[i]=new int*[(ny+pls[dir][1])*(sym_y+1)];		///*********
		for (j=0;j<(ny+pls[dir][1])*(sym_y+1);j++)			///*********
			{
			Solid_Int[i][j]= new int[(nz+pls[dir][2])*(sym_z+1)];///*********
			for (k=0;k<(nz+pls[dir][2])*(sym_z+1);k++)			///*********
				Solid_Int[i][j][k]= 1;
			      
			}
		}
		
	for (i=0; i<nx;i++)
	{
	       Solid[i] = new int*[ny];
	       for (j=0;j<ny;j++)
	       {
	               Solid[i][j] = new int[nz];
	               for (k=0;k<nz;k++)
	                       Solid[i][j][k] = 0;
	       }
	}
	
	Solid2 = new int**[nx];
	for (i=0; i<nx;i++)
	{
	       Solid2[i] = new int*[ny];
	       for (j=0;j<ny;j++)
	       {
	               Solid2[i][j] = new int[nz];
	               for (k=0;k<nz;k++)
	                       Solid2[i][j][k] = 0;
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
			if (pore == 0.0)	{Solid[i][j][k] = 0;Solid2[i][j][k] = 0;}
			
			//if (pore == 1.0) 	{Solid[ci-1][cj-1][ck-1] = 1;sum++;}
			if (pore == 1.0) 	{Solid[i][j][k] = 1;sum++;Solid2[i][j][k] = 1;}
			
		
			
			
		}
		
	fin.close();
	cout<<"Porosity = "<<1-(double(sum)/(nx*ny*nz))<<endl;	
	
	
	
	if (Filter==1)
	{
	//==============FILTER    ============================================
	for(i=0 ; i<nx ; i++)
	for(k=0 ; k<nz ; k++)			
	for(j=0 ; j<ny ; j++)
	if (((i==0) or (i==nx-1) or (j==0) or (j==ny-1) or (k==0) or (k==nz-1)) and (Solid2[i][j][k]==0))
	        Solid2[i][j][k]=-1;
	
	sum=0;sum2=-1;
	
	while (sum!=sum2)
	{
	        sum2=sum;
	for(i=0 ; i<nx ; i++)
	for(k=0 ; k<nz ; k++)			
	for(j=0 ; j<ny ; j++)
	if (Solid2[i][j][k]==0)
	        for (int mi=0;mi<7;mi++)
	        {
	            si=i+e[mi][0];
	                     sj=j+e[mi][1];
	                     sk=k+e[mi][2];
	                     if ((si>=0) and (si<nx) and (sj>=0) and (sj<ny) and (sk>=0) and (sk<nz))
	                              if (Solid2[si][sj][sk]==-1)
	                             {Solid2[i][j][k]=-1;sum++;mi=7;}    
	        }
	       cout<<sum<<endl; 
	}      
	       
	
	/*
	for(i=1 ; i<nx ; i++)
	for(k=0 ; k<nz ; k++)			
	for(j=0 ; j<ny ; j++)

	{
	 if ((Solid2[i-1][j][k]==-1) and (Solid2[i][j][k]<1))
	             for (int mi=0;mi<6;mi++)
	             {
	                     si=i+e1[mi][0];
	                     sj=j+e1[mi][1];
	                     sk=k+e1[mi][2];
	                     if ((si>=0) and (si<nx) and (sj>=0) and (sj<ny) and (sk>=0) and (sk<nz))
	                             if (Solid2[si][sj][sk]<1)
	             Solid2[si][sj][sk]=-1;    
	         }
	         
	}
	
	
	for(k=0 ; k<nz ; k++)			
	for(j=0 ; j<ny ; j++)
	if (Solid2[nx-1][j][k]<=0)
	        Solid2[nx-1][j][k]=-1;
	
	for(i=nx-2 ; i>=0 ; i--)	
	for(k=0 ; k<nz ; k++)			
	for(j=0 ; j<ny ; j++)
	
	{
	 if ((Solid2[i+1][j][k]==-1) and (Solid2[i][j][k]<1))
	             for (int mi=0;mi<6;mi++)
	             {
	                     si=i+e2[mi][0];
	                     sj=j+e2[mi][1];
	                     sk=k+e2[mi][2];
	                     if ((si>=0) and (si<nx) and (sj>=0) and (sj<ny) and (sk>=0) and (sk<nz))
	                             if (Solid2[si][sj][sk]<1)
	             Solid2[si][sj][sk]=-1;    
	         }
	         
	}
	
	
	
	for(k=0 ; k<nz ; k++)			
	for(i=0 ; i<nx ; i++)
	if (Solid2[i][0][k]<=0)
	        Solid2[i][0][k]=-1;
	
	for(j=1 ; j<ny ; j++)
	for(k=0 ; k<nz ; k++)			
	for(i=0 ; i<nx ; i++)	
	{
	 if ((Solid2[i][j-1][k]==-1) and (Solid2[i][j][k]<1))
	             for (int mi=0;mi<6;mi++)
	             {
	                     si=i+e3[mi][0];
	                     sj=j+e3[mi][1];
	                     sk=k+e3[mi][2];
	                     if ((si>=0) and (si<nx) and (sj>=0) and (sj<ny) and (sk>=0) and (sk<nz))
	                             if (Solid2[si][sj][sk]<1)
	             Solid2[si][sj][sk]=-1;    
	         }
	         
	}
	
	
	for(k=0 ; k<nz ; k++)			
	for(i=0 ; i<nx ; i++)
	if (Solid2[i][ny-1][k]<=0)
	        Solid2[i][ny-1][k]=-1;
	for(j=ny-2 ; j>=0 ; j--)
	for(k=0 ; k<nz ; k++)			
	for(i=0 ; i<nx ; i++)	
	
	{
	 if ((Solid2[i][j+1][k]==-1) and (Solid2[i][j][k]<1))
	             for (int mi=0;mi<6;mi++)
	             {
	                     si=i+e4[mi][0];
	                     sj=j+e4[mi][1];
	                     sk=k+e4[mi][2];
	                     if ((si>=0) and (si<nx) and (sj>=0) and (sj<ny) and (sk>=0) and (sk<nz))
	                             if (Solid2[si][sj][sk]<1)
	             Solid2[si][sj][sk]=-1;    
	         }
	         
	}
	
	
	
	for(j=0 ; j<ny ; j++)
	for(i=0 ; i<nx ; i++)	
	if (Solid2[i][j][0]<=0)
	        Solid2[i][j][0]=-1;
	
	for(k=1 ; k<nz ; k++)	
	for(j=0 ; j<ny ; j++)
	for(i=0 ; i<nx ; i++)	
	{
	 if ((Solid2[i][j][k-1]==-1) and (Solid2[i][j][k]<1))
	             for (int mi=0;mi<6;mi++)
	             {
	                     si=i+e5[mi][0];
	                     sj=j+e5[mi][1];
	                     sk=k+e5[mi][2];
	                     if ((si>=0) and (si<nx) and (sj>=0) and (sj<ny) and (sk>=0) and (sk<nz))
	                             if (Solid2[si][sj][sk]<1)
	             Solid2[si][sj][sk]=-1;    
	         }
	         
	}
	
	
	
	for(j=0 ; j<ny ; j++)
	for(i=0 ; i<nx ; i++)	
	if (Solid2[i][j][nz-1]<=0)
	        Solid2[i][j][nz-1]=-1;
	for(k=nz-2 ; k>=0 ; k--)			
	for(j=0 ; j<ny ; j++)
	for(i=0 ; i<nx ; i++)	
	
	{
	 if ((Solid2[i][j][k+1]==-1) and (Solid2[i][j][k]<1))
	             for (int mi=0;mi<6;mi++)
	             {
	                     si=i+e6[mi][0];
	                     sj=j+e6[mi][1];
	                     sk=k+e6[mi][2];
	                     if ((si>=0) and (si<nx) and (sj>=0) and (sj<ny) and (sk>=0) and (sk<nz))
	                             if (Solid2[si][sj][sk]<1)
	             Solid2[si][sj][sk]=-1;    
	         }
	         
	}
	
	
	*/
	
	
	//============================================================
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	sum=0;
	for(k=0 ; k<nz ; k++)			
	for(j=0 ; j<ny ; j++)
	for(i=0 ; i<nx ; i++)
        {	
	if (Solid2[i][j][k]==0)
	        Solid[i][j][k]=1;
	if (Solid[i][j][k]==1)
	        sum++;
	}
	
	cout<<"Porosity2 = "<<1-(double(sum)/(nx*ny*nz))<<endl;	
	//cout<<sum<<endl;
	}
	//=======================FILTER============================
	
	
	
	
	
	
	if (dir==4)
	{
		for (int i=0;i<nx;i++)
			{
			for (int j=0;j<ny;j++)
				{
				Solid[i][j][0]=1;
				Solid[i][j][nz-1]=1;
				}
			for (int k=0;k<nz;k++)
				{
				Solid[i][0][k]=1;
				Solid[i][ny-1][k]=1;
				}
			}
	}	

if (dir==5)
	{
		for (int j=0;j<ny;j++)
			{
			for (int i=0;i<nx;i++)
				{
				Solid[i][j][0]=1;
				Solid[i][j][nz-1]=1;
				}
			for (int k=0;k<nz;k++)
				{
				Solid[0][j][k]=1;
				Solid[nx-1][j][k]=1;
				}
			}
	}	


if (dir==6)
	{
		for (int k=0;k<nz;k++)
			{
			for (int j=0;j<ny;j++)
				{
				Solid[0][j][k]=1;
				Solid[nx-1][j][k]=1;
				}
			for (int i=0;i<nx;i++)
				{
				Solid[i][0][k]=1;
				Solid[i][ny-1][k]=1;
				}
			}
	}	







		
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


	if (add_porous_plate==1)
		for(int k=0 ; k<nz1*Zoom ; k=k++)				///*********
		for(int j=0 ; j<ny1*Zoom; j=j++)
			if (((k%2==1) and (j%2==0)) or ((k%2==0) and (j%2==1)))
			if (porous_position<0)
				Solid_Int2[nx1*Zoom-1][j][k]=1;
			else
				Solid_Int2[porous_position][j][k]=1;



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
	
	
	if (vtk_check==1)
	        {
	//================================================================
	cout<<"Start non-conective pore vtk file"<<endl;
	cout<<nx<<"         "<<ny<<"         "<<nz<<endl;
	ostringstream name3;
	name3<<Non_conective_poreFileNameOut;
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
	ofstream out3;
	out3.open(name3.str().c_str());
	out3<<"# vtk DataFile Version 2.0"<<endl;
	out3<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out3<<"ASCII"<<endl;
	out3<<"DATASET STRUCTURED_POINTS"<<endl;
	out3<<"DIMENSIONS         "<<nx<<"         "<<ny<<"         "<<nz<<endl;       ///*********
	out3<<"ORIGIN 0 0 0"<<endl;
	out3<<"SPACING 1 1 1"<<endl;
	out3<<"POINT_DATA     "<<nx*ny*nz<<endl;				///*********
	out3<<"SCALARS sample_scalars float"<<endl;
	out3<<"LOOKUP_TABLE default"<<endl;
	
	
	for (k=0;k<nz;k++)
	for (j=0;j<ny;j++)
	for (i=0;i<nx;i++)
		out3<<Solid2[i][j][k]<<" ";
	
	

	out3.close();

	cout<<"non-conective pore vtk  file ouput complete"<<endl;
	        }
	        
	        
	        
}


