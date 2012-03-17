#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>






const int depth=6;
const int nx=300;
const int ny=300;
const int nz=300;
using namespace std; 
      
int main (int argc , char * argv [])
{


int size[depth]={300,150,120,100,60,50};
char poreFileName[128]="LV60.dat";
char outputname[128]="prosity_output.dat";

int num[depth];
double ave[depth];
double min[depth];
double max[depth];
double sigma[depth];

int ls1,ls2,ls3,sta2,sta3;
double** porosity;


//================
if (ny>nx)
        sta2=nx;
else
        sta2=ny;

if (nz<sta2)
        sta2=nz;

sta3=sta2/10;
//=================
double* porosity2;
double* porosity3;
porosity2 = new double[sta2];
porosity3 = new double[sta2];
for (int i=0;i<sta2;i++)
{
        porosity2[i]=0.0;
        porosity3[i]=0.0;
       
}

//cout<<sta2<<"         "<<sta3<<endl;





for (int i=0;i<depth;i++)
        {
                ls1=nx/size[i];
                ls2=ny/size[i]; 
                ls3=nz/size[i];
                if (ls1<ls2)
                        num[i]=ls1;
                else    
                        num[i]=ls2;
                
                if (ls3<num[i])
                        num[i]=ls3;
                //cout<<num[i]<<endl;
        }

int sum=0;


bool*** Solid;


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
	
	
	Solid = new bool**[nx];
	porosity = new double*[depth];
		
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
		
		for(int i=0;i<depth;i++)
		        porosity[i]= new double[num[depth-1]*num[depth-1]*num[depth-1]];
		
		for (int i=0;i<depth;i++)
		        for (int j=0;j<num[depth-1]*num[depth-1]*num[depth-1];j++)
		                porosity[i][j]=0.0;
		        
		        
		        
		        //for (int i=0;i<depth;i++)
		           //    cout<<num[i]<<endl;
		        
		
		
	for(k=0 ; k<nz ; k++)				///*********
	for(j=0 ; j<ny ; j++)
	for(i=0 ; i<nx ; i++)				///*********


	
		{	
			
			fin >> pore;
			Solid[i][j][k]=pore;
			if (pore==1.0)
			{
			        sum+=1;
			        
			        for (int ii=sta3;ii<sta2;ii++)
			                {
			                        if ((i<ii) and (j<ii) and (k<ii))
			                                porosity2[ii]+=1.0;
			                }
			               
			        for (int ii=0;ii<depth;ii++)
			                {
			                        ls1=i/size[ii];
			                        ls2=j/size[ii];
			                        ls3=k/size[ii];
			                        if ((ls1<num[ii]) and (ls2<num[ii]) and (ls3<num[ii]))
			                                porosity[ii][ls1*num[ii]*num[ii]+ls2*num[ii]+ls3]+=1.0;
			                        
			                        
			                }
			}
			
			
			
		}
		
	fin.close();
		
	for (int i=0;i<depth;i++)
	{        max[i]=0.0;
	        ave[i]=0.0;
	        min[i]=1.0;
	        sigma[i]=0.0;
	}
	
	
	for (int i=0;i<depth;i++)
		        for (int j=0;j<num[depth-1]*num[depth-1]*num[depth-1];j++)
		        porosity[i][j]=1-porosity[i][j]/(double)(size[i]*size[i]*size[i]);
	
	for (int i=0;i<depth;i++)
	                {
		        for (int j=0;j<num[i]*num[i]*num[i];j++)
		        {
		                if (porosity[i][j]<min[i])
		                        min[i]=porosity[i][j];
		                if (porosity[i][j]>max[i])
		                        max[i]=porosity[i][j];
		                
		                ave[i]+=porosity[i][j];
		                sigma[i]+=porosity[i][j]*porosity[i][j];
		                
		  
		        } 
		        
		        ave[i]/=num[i]*num[i]*num[i];
		        sigma[i]/=num[i]*num[i]*num[i];
		        sigma[i]-=ave[i]*ave[i];
		        sigma[i]=sqrt(sigma[i]);
		        }
		        
	
		        for (int i=0;i<sta2;i++)
		                porosity2[i]=1-porosity2[i]/(i*i*i);
		        
		        
		        
		        
	cout<<"Porosity = "<<1-(double(sum)/(nx*ny*nz))<<endl;	
	cout<<porosity[1][0]<<endl;
	//cout<<sum<<endl;	
	
	ostringstream name;
	name<<outputname;
	ofstream out;
	out.open(name.str().c_str());	
	
	for (int i=0;i<depth;i++)
		        for (int j=0;j<num[i]*num[i]*num[i];j++)
		        out<<size[i]<<" "<<porosity[i][j]<<endl;
		
		out.close();
	
	
	ostringstream name2;
	name2<<outputname<<"_sta";

	out.open(name2.str().c_str());	
	
	for (int i=0;i<depth;i++)
		        //for (int j=0;j<num[i]*num[i]*num[i];j++)
		        out<<size[i]<<" "<<ave[i]<<"  "<<min[i]<<"  "<<max[i]<<" "<<sigma[i]<<endl;
		
		out.close();
		
		
	ostringstream name3;
	name3<<outputname<<"_sta2";

	out.open(name3.str().c_str());	
	
	for (int i=0;i<depth;i++)
		        if (num[i]>1)
		        out<<size[i]<<"  "<<sigma[i]<<endl;
		
		out.close();
		
	ostringstream name4;
	name4<<outputname<<"_sub";

	out.open(name4.str().c_str());	
	
	for (int i=sta3;i<sta2;i++)
		        //for (int j=0;j<num[i]*num[i]*num[i];j++)
		        out<<i<<"  "<<porosity2[i]<<endl;
		
		out.close();	
		

		
}


