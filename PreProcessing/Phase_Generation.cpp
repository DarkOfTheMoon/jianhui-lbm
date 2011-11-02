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

ostringstream name;
	name<<"phase"<<"_x.dat";
	ofstream out;
	out.open(name.str().c_str());

	for(int k=0 ; k<nz; k++)						///*********
		for(int j=0 ; j<ny; j++)					///*********
			for(int i=0 ; i<nx; i++)
			//if (i<=3)
		        if ((i>=3) and (i<=5))
			out<<1.0<<endl;
			else
			out<<0.0<<endl;
	out.close();
}
