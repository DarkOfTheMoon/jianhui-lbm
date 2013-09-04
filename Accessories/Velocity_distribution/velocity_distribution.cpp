#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include <math.h>

#include <cstdlib>
#include <algorithm>


using namespace std; 

int main (int argc , char * argv [])
{
//ifstream fin(argv[1]);
int NCHAR=128;
	char     filename[128], dummy[128+1];
	int      dummyInt;

const int max_cluster=1000000;
double crival=0.0000001;


int nx,ny,nz;
int vtk_format;


const int e[18][3]=
{{1,0,0,},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,1,0},{1,-1,0},{-1,-1,0},{0,1,1},
{0,-1,1},{0,1,-1},{0,-1,-1},{1,0,1},{-1,0,1},{1,0,-1},{-1,0,-1}};

char poreFileName[128];

char poreFileNameVTK[128];
char poreFileNameOut[128];
int mark,ii,jj,kk,loop,sum2;

int sum=0;
int sum1=0;
int pore_sum;
int sum3;
int sn[max_cluster];


int* sn2;
sn2=new int[max_cluster];

int phase_ind;
int exp_vtk;
double* xres;
double* xind;
int local_res=50;
double loc_res=double(9.0/local_res);


int lnx,rnx,lny,rny,lnz,rnz;

for (int i=0;i<max_cluster;i++)
	sn[i]=0;


ifstream fins(argv[1]);
							fins.getline(dummy, NCHAR);
	fins >> poreFileName;				fins.getline(dummy, NCHAR);
	fins >> nx>> ny>>nz;				fins.getline(dummy, NCHAR);
	fins >> poreFileNameVTK;				fins.getline(dummy, NCHAR);
	fins >> vtk_format;				fins.getline(dummy, NCHAR);
	fins >> phase_ind;				fins.getline(dummy, NCHAR);
	fins >> lnx>>rnx;                                fins.getline(dummy, NCHAR);
	fins >> lny>>rny;                                fins.getline(dummy, NCHAR);
	fins >> lnz>>rnz;                                fins.getline(dummy, NCHAR);
	fins >> crival;				 fins.getline(dummy, NCHAR);
fins.close();	


//cout<<phase_ind<<endl;


xres = new double[phase_ind];
xind = new double[phase_ind];
for (int i=0;i<phase_ind;i++)
        xres[i]=0.0,xind[i]=0.0;


double*** Solid;
double pore,pore1,pore2,pore3;


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



	Solid = new double**[nx];
	
		
	for (int i=0; i<nx;i++)
	{
	       Solid[i] = new double*[ny];
	       for (int j=0;j<ny;j++)
	       {
	               Solid[i][j] = new double[nz];
	               for (int k=0;k<nz;k++)
	                       Solid[i][j][k] = 0;
	       }
	}
		
	cout<<endl;
	cout<<"=================================================================\n";
	cout<<"		    		Velocity Distribution Programme\n";
	cout<<"	Jianhui Yang - All right reserved	JAN2013	\n";
	cout<<"=================================================================\n";
	cout<<endl;
	cout<<"Start reading source geometry file"<<endl;

	
	if (vtk_format==1)
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

double vave=0.0;
double vmax=0.0;

	for(int k=0 ; k<nz ; k++)				///*********
	for(int j=0 ; j<ny ; j++)
	for(int i=0 ; i<nx ; i++)				///*********


	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore1 >> pore2 >> pore3;
			pore=sqrt(pore1*pore1+pore2*pore2+pore3*pore3);
			//if (pore == 0.0)	{Solid[ci-1][cj-1][ck-1] = 0;}
			if (pore > crival)	{Solid[i][j][k] = pore;        
			                        if ((i>=lnx) and (i<rnx) and (j>=lny) and (j<rny) and (k>=lnz) and (k<rnz)) 
			                                        {
			                                        sum++;vave+=pore;
			                                        if (pore>vmax) 
			                                                vmax=pore;
			                                        }
			                                   }
			else
				{Solid[i][j][k] = 0.0;}
			
			
		
			
			
		}
		
	fin.close();
		
	pore_sum=sum;
	cout<<"porosity="<<double(pore_sum)/double((rnx-lnx)*(rny-lny)*(rnz-lnz))<<endl;
	//cout<<"non wetting saturation="<<double(sum1)/double(sum)<<endl;
	cout<<"READING COMPLETE"<<endl;	
	vave/=(double)sum;
	
	cout<<"Average Velocity="<<vave<<endl;
	cout<<"Maximum Velocity="<<vmax<<endl;

	sum3=0;

	
	
	for (int lsi=0;lsi<phase_ind;lsi++)
		xind[lsi]=10.0/double(phase_ind)*(lsi);
		
	sum3=0;
	
	double interval;
	double cal_max=vave*5;
	//interval=vmax/(phase_ind-1);
	//interval=cal_max/(phase_ind-1);
	

	for (int i=lnx;i<rnx;i++)
		for (int j=lny;j<rny;j++)
			for (int k=lnz;k<rnz;k++)
			if (Solid[i][j][k]>crival)
						for (int lsi=0;lsi<phase_ind-1;lsi++)
				          	if ((Solid[i][j][k]>xind[lsi]*vave) and (Solid[i][j][k]<=vave*xind[lsi+1]))
				          		{xres[lsi]+=1.0;sum3++;
							//	cout<<i<<" "<<j<<" "<<k<<endl;
							}
				
				
				
			
	//	cout<<"@@@@@@@@@@@@@@@"<<endl;		
	for (int ls=0;ls<phase_ind;ls++)
	        xres[ls]=xres[ls]/(double)sum;
	
	
	cout<<"Statistical portion="<<double(sum3)/sum<<endl;
	        //cout<<sum<<endl;
	
				//cout<<xres[0]<<"                @@@@@@@@@@@@@@@@"<<endl;
		ostringstream name3;
		name3<<"velocity_distribution.txt";

		ofstream out3;
		out3.open(name3.str().c_str());
		//cout<<"@@@@@@@@@    "<<phase_ind<<endl;
		//cout<<xres[0]<<"                @@@@@@@@@@@@@@@@"<<endl;
		sum2=0;
			for (int lsi=0;lsi<phase_ind;lsi++)   
		        out3<<xind[lsi]<<"        "<<xres[lsi]<<endl;
	
		out3.close();





	
}


