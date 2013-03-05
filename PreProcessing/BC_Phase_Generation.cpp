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
int nx=30;
int ny=100;
int nz=3;


ostringstream name;
	name<<"BC"<<".dat";
	ofstream out;
	out.open(name.str().c_str());

	for(int k=0 ; k<nz; k++)						///*********
		for(int j=0 ; j<ny; j++)					///*********
			for(int i=0 ; i<nx; i++)
			//if (i<=3)
		        if ((i==0) or (i==nx-1))
			out<<1.0<<endl;
			else
			out<<0.0<<endl;
	out.close();
	
	
ostringstream name2;
	name2<<"phase"<<"_x.dat";
	ofstream out2;
	out2.open(name2.str().c_str());

	for(int k=0 ; k<nz; k++)						///*********
		for(int j=0 ; j<ny; j++)					///*********
			for(int i=0 ; i<nx; i++)
			//if (i<=3)
		        if ((i>=3) and (i<=5))
			out2<<1.0<<endl;
			else
			out2<<-1.0<<endl;
	out2.close();
}
