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
int nx=256;
int ny=130;
int nz=130;
srand((unsigned)time(0));
double resa;
double sum;

ostringstream name;
	name<<"phase"<<"_x.dat";
	ofstream out;
	out.open(name.str().c_str());
	sum=0;
	
	for(int k=0 ; k<nz; k++)						///*********
		for(int j=0 ; j<ny; j++)					///*********
			for(int i=0 ; i<nx; i++)
			{
			resa=(double(rand()%10000))/10000;
			if (resa>0.6)
			        {out<<1.0<<endl;sum+=1;}
			else
			out<<-1.0<<endl;
			}
	out.close();
	sum/=(nx*ny*nz);
	cout<<"Initial satuation is "<<sum<<endl;
}
