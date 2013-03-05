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
//ifstream fin(argv[1]);
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


	
	cout<<"Start reading source geometry file"<<endl;	
	for(int k=0 ; k<nz ; k++)				///*********
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)				///*********
	

	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore;
			
			//if (pore == 0.0)	{Solid[ci-1][cj-1][ck-1] = 0;}
			if (pore == 0)	{Solid[i][j][k] = 0;sum++;reci=i;recj=j;reck=k;}
			else
			//if (pore == 1.0) 	{Solid[ci-1][cj-1][ck-1] = 1;sum++;}
			//if (pore == 1.0) 	
				{Solid[i][j][k] = 1;}
			
		
			
			
		}
	cout<<"Porosity = "<<(double(sum)/(nx*ny*nz))<<endl;	
	sum_rec=sum;
	fin.close();

	if (fil==1)
	{
	cout<<endl;
	cout<<"=================================================================\n";
	cout<<"		    Preprocessing Module-Remove isolated pores\n";
	cout<<"	Jianhui Yang - All right reserved	OCT2012	\n";
	cout<<"=================================================================\n";
	cout<<endl;
	
	
	sum2=sum;

	sum3=sum;

	for (int i=0;i<nx;i++)
		for (int j=0;j<ny;j++)
			for (int k=0;k<nz;k++)
			if ((i==0) and (Solid[i][j][k]==0))
				{Solid[i][j][k]=2;sum--;}
			else
				if ((i==nx-1) and (Solid[i][j][k]==0))
					{Solid[i][j][k]=3;sum--;}
				else
					if ((j==0) and (Solid[i][j][k]==0))
						{Solid[i][j][k]=4;sum--;}
					else
						if ((j==ny-1) and (Solid[i][j][k]==0))
							{Solid[i][j][k]=5;sum--;}
						else
							if ((k==0) and (Solid[i][j][k]==0))
								{Solid[i][j][k]=6;sum--;}
							else
								if ((k==nz-1) and (Solid[i][j][k]==0))
									{Solid[i][j][k]=7;sum--;}
	

	
	loop=0;mark=1;
	while (mark>0)
	{
	loop++;cout<<loop<<"	The residule sum="<<sum<<"	ratio="<<double(sum)/double(sum2)<<endl;
	mark=0;
	for (int i=0;i<nx;i++)
		for (int j=0;j<ny;j++)
			for (int k=0;k<nz;k++)
			if (Solid[i][j][k]>1)
				for (int ls=0;ls<18;ls++)
				{
				ii=i+e[ls][0];
				jj=j+e[ls][1];
				kk=k+e[ls][2];
				
				if ((ii>=0) and (ii<nx) and (jj>=0) and (jj<ny) and (kk>=0) and (kk<nz))
				{
					if (Solid[ii][jj][kk]==0)
						{Solid[ii][jj][kk]=Solid[i][j][k];sum--;mark++;if (Solid[i][j][k]==10) sum3--;}
					if ((Solid[ii][jj][kk]>1) and (Solid[ii][jj][kk]!=Solid[i][j][k]))
						{
						if (Solid[i][j][k]<10)
							sum3--;
						Solid[i][j][k]=10;
						
						}	
				}			


				}

				
	}
	

	
	
	cout<<"The residule sum="<<sum<<"	ratio="<<double(sum)/double(sum2)<<endl;
	cout<<"Initial seed for boundary seraching, sum="<<sum3<<"   ratio="<<double(sum3)/double(sum2)<<endl;


	cout<<endl;
	cout<<"Start searching isolated pores on boundaries"<<endl;
	//sum2=sum;
	
	loop=0;mark=1;
	while (mark>0)
	{
	loop++;cout<<loop<<"	The residule sum="<<sum3<<"	ratio="<<double(sum3)/double(sum2)<<endl;
	mark=0;
	for (int i=0;i<nx;i++)
		for (int j=0;j<ny;j++)
			for (int k=0;k<nz;k++)
			if (Solid[i][j][k]==10)
				for (int ls=0;ls<18;ls++)
				{
				ii=i+e[ls][0];
				jj=j+e[ls][1];
				kk=k+e[ls][2];
				
				if ((ii>=0) and (ii<nx) and (jj>=0) and (jj<ny) and (kk>=0) and (kk<nz) and (Solid[ii][jj][kk]>1) and (Solid[ii][jj][kk]<10))
				{Solid[ii][jj][kk]=10;sum3--;mark++;}

				}

				
	}

	cout<<"The residule sum="<<sum3<<"	ratio="<<double(sum3)/double(sum2)<<endl;


	for(int k=0 ; k<nz ; k++)				
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)	
	if (Solid[i][j][k]==10)
		Solid[i][j][k]=0;
	else	
		Solid[i][j][k]=1;			
	//============decomposition complete=======================
	
	//==================================================
	
	}


bool*** Solid_Int;

int*** Solid_Int2;





if (geo_mod==1)
{

int pls[7][3]={{0,2,2},{2,0,2},{2,2,0},{0,0,0},{0,0,0},{0,0,0}};

int sum=0;



int nx1,ny1,nz1;
nx1=(nx+pls[dir][0])*(sym_x+1)+add_buf_x_n+add_buf_x_p;
ny1=(ny+pls[dir][1])*(sym_y+1)+add_buf_y_n+add_buf_y_p;
nz1=(nz+pls[dir][2])*(sym_z+1)+add_buf_z_n+add_buf_z_p;







double pore;
	int i, j, k,ci,cj,ck;
	
	Solid_Int = new bool**[(nx+pls[dir][0])*(sym_x+1)];	///*********
	
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
			
	if (add_porous_plate==2)
		for(int k=0 ; k<nz1*Zoom ; k=k++)				///*********
		for(int j=0 ; j<ny1*Zoom; j=j++)
			if ((k%3==0) or (j%3==2))
			if (porous_position<0)
				Solid_Int2[nx1*Zoom-1][j][k]=1;
			else
				Solid_Int2[porous_position][j][k]=1;
			
	if (add_porous_plate==3)
		for(int k=0 ; k<nz1*Zoom ; k=k++)				///*********
		for(int j=0 ; j<ny1*Zoom; j=j++)
			if ((k%4==0) or (j%4==2))
			if (porous_position<0)
				Solid_Int2[nx1*Zoom-1][j][k]=1;
			else
				Solid_Int2[porous_position][j][k]=1;

	/*
	for (int i=0;i<nx;i++)
		{
		for (int j=0;j<ny;j++)
			delete [] Solid[i][j];
		delete [] Solid[i];
		}
	delete [] Solid;
	*/


	delete [] Solid[0];
	
	nx=nx1*Zoom;
	ny=ny1*Zoom;
	nz=nz1*Zoom;

	Solid = new int**[nx1*Zoom];	///*********
	
	for (i=0;i<nx1*Zoom;i++)				///*********
		Solid[i]=new int*[ny1*Zoom];

	Solid[0][0]=new int[nx1*ny1*nz1*Zoom*Zoom*Zoom];

	
 	for (int i=1;i<ny1*Zoom;i++)
               Solid[0][i]=Solid[0][i-1]+nz1*Zoom;
       
       for (int i=1;i<nx1*Zoom;i++)
       {
               Solid[i][0]=Solid[i-1][0]+ny1*nz1*Zoom*Zoom;
               for (int j=1;j<ny1*Zoom;j++)
                       Solid[i][j]=Solid[i][j-1]+nz1*Zoom;
       }
	
	for(int k=0 ; k<nz ; k++)				
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)	
		Solid[i][j][k]=Solid_Int2[i][j][k];

	
	
	


}


if (expvtk==1)
	{
	cout<<"Start writing VTK file"<<endl;
	cout<<nx<<"         "<<ny<<"         "<<nz<<endl;

	ostringstream name;
	name<<poreFileNameVTK;
	ofstream out;
	out.open(name.str().c_str());
	
	//if (geo_mod==1)
	//{
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"binary"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<nz<<"         "<<ny<<"         "<<nx<<endl;       ///*********
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<nx*ny*nz<<endl;				///*********
	out<<"SCALARS sample_scalars int"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	out.write((char *)(&Solid[0][0][0]), sizeof(int)*nx*ny*nz); 
	/*
	}
	else
	{
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
	for (int j=0;j<ny;j++)
	for (int i=0;i<nx;i++)
		out<<Solid[i][j][k]<<" ";
	}
	}
	*/

	out.close();

	cout<<"VTK file ouput COMPLETE"<<endl;
	}




if (expdat==1)
{



	cout<<"Start writing DAT file"<<endl;
	cout<<nx<<"	"<<ny<<"	"<<nz<<endl;
	ostringstream name2;
	name2<<poreFileNameOut;
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
	ofstream out2;
	out2.open(name2.str().c_str());
	
	if (bindat==1)
	out2.write((char *)(&Solid[0][0][0]), sizeof(int)*nx*ny*nz); 
	else
	for (int k=0;k<nz;k++)
	{
		//cout<<k<<endl;
		for (int j=0;j<ny;j++)
		for (int i=0;i<nx;i++)
			out2<<Solid[i][j][k]<<" ";
	}
	
	out2.close();

	cout<<"DAT file ouput complete"<<endl;
	cout<<endl;

}
	
	



if (mesh_par==1)
{

	cout<<endl;
	cout<<"MESH PARTITION INITIALIZATION START"<<endl;

int ii,jj,kk;
sum2=0;

sum=1;


	cout<<sum<<endl;
	for(int k=0 ; k<nz ; k++)				///*********
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)				///*********


	//while (!fin.eof())                                        //**********
		{	
			
			
			
			if (Solid[i][j][k] == 0.0)	{Solid[i][j][k] = sum;sum++;}
			else
			
			if (Solid[i][j][k]== 1.0) 	{Solid[i][j][k] = 0;}
			
		
			
			
		}

	
	for(int k=0 ; k<nz ; k++)				
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)
		if (Solid[i][j][k]>0)
		{//cout<<i<<"	"<<j<<"		"<<k<<endl;
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
		
		if ((ii>=0) and (ii<nx) and (jj>=0) and (jj<ny) and (kk>=0) and (kk<nz))
			if (Solid[ii][jj][kk]>0)
				{sum2++;}
		}	
		}
	
	sum2=sum2/2;

	cout<<endl;
	cout<<"EXPORT METIS INPUT DATA FILE"<<endl;
	cout<<sum-1<<"			"<<sum2<<endl;
	ostringstream name3;
	name3<<poreFileNameMET;
	ofstream out;
	out.open(name3.str().c_str());
	out<<sum-1<<"			"<<sum2<<endl;
	
	
	
	
	for (int k=0;k<nz;k++)
	for (int j=0;j<ny;j++)
	for (int i=0;i<nx;i++)
	if (Solid[i][j][k]>0)
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
		
		if ((ii>=0) and (ii<nx) and (jj>=0) and (jj<ny) and (kk>=0) and (kk<nz))
			if (Solid[ii][jj][kk]>0)
				out<<Solid[ii][jj][kk]<<" ";
		}

	out<<endl;
	}
	cout<<"EXPORT COMPLETE"<<endl;

	
	out.close();




	

}	
	
	

if (mesh_par>1)
{

	cout<<endl;
	cout<<"MESH PARTITION INITIALIZATION START"<<endl;

sum=0;
for(int k=0 ; k<nz ; k++)			
	for(int j=0 ; j<ny ; j++)	
	for(int i=0 ; i<nx ; i++)

		if (Solid[i][j][k]==0)
			{sum++;}
		else
			{Solid[i][j][k]=-1;}


int nxref,nyref,nzref;
int divnumori;
divnumori=mesh_par;
int divnum;
int oddval=0;
int evennum=0;
double oddpor;
int dir;
int dint;
int sumin,numgeonum;
int* sum_loc;


int *nnx,*nny,*nnz,*npx,*npy,*npz;
nnx=new int[divnumori];
nny=new int[divnumori];
nnz=new int[divnumori];
npx=new int[divnumori];
npy=new int[divnumori];
npz=new int[divnumori];




cout<<divnumori<<endl;
divnum=divnumori;
while (divnum%2==0)
        {evennum++;divnum=divnum/2;}
oddval=divnum;
divnum=divnumori;

cout<<evennum<<"         "<<oddval<<endl;
 

	sum_loc = new int[divnum+1];
	for (int i=0;i<=divnum;i++)
	        sum_loc[i]=0;
		
	
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
	//cout<<"@@@@@@@@@@@"<<endl;
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
	
	




	if (partition_vtk==1)
	{
	cout<<endl;
	cout<<"Start decomposed VTK file"<<endl;
	cout<<nx<<"         "<<ny<<"         "<<nz<<endl;

	ostringstream name;
	name<<"Partition.vtk";
	ofstream out;
	out.open(name.str().c_str());
	
	//if (geo_mod==1)
	//{
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"binary"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<nz<<"         "<<ny<<"         "<<nx<<endl;       ///*********
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<nx*ny*nz<<endl;				///*********
	out<<"SCALARS sample_scalars int"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	out.write((char *)(&Solid[0][0][0]), sizeof(int)*nx*ny*nz); 


	out.close();

	cout<<"Decomposed VTK file ouput COMPLETE"<<endl;
	}

	
	
	cout<<endl;
	cout<<"Start writing MESH DAT file"<<endl;
	cout<<nx<<"	"<<ny<<"	"<<nz<<endl;
	ostringstream name3;
	name3<<poreFileName<<"_"<<mesh_par<<".dat";
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
	ofstream out3;
	out3.open(name3.str().c_str());
	
	if (bindat==1)
	out3.write((char *)(&Solid[0][0][0]), sizeof(int)*nx*ny*nz); 
	else
	for (int k=0;k<nz;k++)
	{
		//cout<<k<<endl;
		for (int j=0;j<ny;j++)
		for (int i=0;i<nx;i++)
			out3<<Solid[i][j][k]<<" ";
	}
	
	out3.close();

	cout<<"DAT file ouput complete"<<endl;
	cout<<endl;

	






	

}	

	
	
	
	
	
	
}


