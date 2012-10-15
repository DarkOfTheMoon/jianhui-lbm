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





int nx,ny,nz;



const int e[18][3]=
{{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};

char poreFileName[128];

char poreFileNameVTK[128];
char poreFileNameOut[128];
int mark,ii,jj,kk,loop,sum2;

int sum=0;

ifstream fins(argv[1]);
							fins.getline(dummy, NCHAR);
	fins >> poreFileName;				fins.getline(dummy, NCHAR);
	fins >> nx>> ny>>nz;				fins.getline(dummy, NCHAR);
	fins >> poreFileNameVTK;				fins.getline(dummy, NCHAR);
	fins >> poreFileNameOut;				fins.getline(dummy, NCHAR);


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
		
	cout<<endl;
	cout<<"=================================================================\n";
	cout<<"		    Preprocessing Module-Remove isolated pores\n";
	cout<<"	Jianhui Yang - All right reserved	OCT2012	\n";
	cout<<"=================================================================\n";
	cout<<endl;
	cout<<"Start reading source geometry file"<<endl;





	for(int k=0 ; k<nz ; k++)				///*********
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)				///*********


	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore;
			
			//if (pore == 0.0)	{Solid[ci-1][cj-1][ck-1] = 0;}
			if (pore == 0.0)	{Solid[i][j][k] = 0;sum++;}
			else
			//if (pore == 1.0) 	{Solid[ci-1][cj-1][ck-1] = 1;sum++;}
			//if (pore == 1.0) 	
				{Solid[i][j][k] = 1;}
			
		
			
			
		}
		
	fin.close();
		
	
	cout<<"Porosity = "<<(double(sum)/(nx*ny*nz))<<endl;	
	
	sum2=sum;



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
						{Solid[ii][jj][kk]=Solid[i][j][k];sum--;mark++;}
					if ((Solid[ii][jj][kk]>1) and (Solid[ii][jj][kk]!=Solid[i][j][k]))
						{Solid[i][j][k]=10;}	
				}			


				}

				
	}
			
	
	
	cout<<"The residule sum="<<sum<<"	ratio="<<double(sum)/double(sum2)<<endl;


	cout<<endl;
	cout<<"Start searching isolated pores on boundaries"<<endl;
	//sum2=sum;
	
	loop=0;mark=1;
	while (mark>0)
	{
	loop++;cout<<loop<<"	The residule sum="<<sum<<"	ratio="<<double(sum)/double(sum2)<<endl;
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
				{Solid[ii][jj][kk]=10;sum--;mark++;}

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
		if (Solid[i][j][k]==1)
		out<<"2 ";
			else
			if (Solid[i][j][k]==10)
			out<<"0 ";
			else
			out<<"1 ";
			
			
		//out<<Solid[i][j][k]<<" ";
	}
	

	//out.write((char *)(&Solid_Int2[0][0][0]), sizeof(int)*nx1*ny1*nz1*Zoom*Zoom*Zoom); 

	out.close();

	cout<<"VTK file ouput COMPLETE"<<endl;


	cout<<"Start writing DAT file"<<endl;
	cout<<nx<<"         "<<ny<<"         "<<nz<<endl;
	ostringstream name2;
	name2<<poreFileNameOut;
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
	ofstream out2;
	out2.open(name2.str().c_str());


	//for (int i=0;i<nx;i++)
	//for (int j=0;j<ny;j++)

	for (int k=0;k<nz;k++)
	{
	//cout<<k<<endl;
	for (int j=0;j<ny;j++)
	for (int i=0;i<nx;i++)
		if (Solid[i][j][k]==10)
			out2<<"0 ";
		else
			out2<<"1 ";
		//out<<Solid[i][j][k]<<" ";
	}
	
	
	//out2.write((char *)(&Solid_Int2[0][0][0]), sizeof(int)*nx1*ny1*nz1*Zoom*Zoom*Zoom); 

	out2.close();

	cout<<"DAT file ouput complete"<<endl;
	
	
}


