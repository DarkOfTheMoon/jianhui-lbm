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



int total_number=5;
char prefix[128]="rel_";
char poreFileName[128]="Summary.outdat";



float var1,var2,var3;	
var1=0;var2=0;var3=0;
char chd[128];
int testd;
ostringstream name;
ostringstream name2;
name2<<poreFileName;
        FILE *ftest;
        FILE *fp;
	ifstream fin;
	ofstream out;
int NCHAR=128;	
char dummy[128+1];
ofstream fins;	
	fins.open(poreFileName,ios::trunc);
	fins.close();
	
	
	
for (int i=1;i<=total_number;i++)
{
        
        name.str("");
	name<<prefix<<i<<"_Results.txt";
	
	ftest = fopen(name.str().c_str(), "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file (" << poreFileName <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	else
	        {
	                
	                fin.open(name.str().c_str());
	                fin.getline(dummy, NCHAR);        //       The30800th computation result:     
	                fin.getline(dummy, NCHAR);        //The Maximum velocity is: 0.0120153   Re_l=1.52508   Re_g=1.52508
	                fin.getline(dummy, NCHAR);           // Courant Number=0.0120153	 Capillary Num=0.0076254
	                fin.getline(dummy, NCHAR);        //The max relative error of velocity is: 7.971677e-01
	                fin>> chd >> chd>>chd>>chd >> chd>>chd>>chd>>var2;fin.getline(dummy, NCHAR);         //   The relative permeability of component 1 is -3.992098e+00, -3.079636e-03, 1.535446e-07
	                fin >>chd >> chd>>chd>>chd >> chd>>chd>>chd>> var3;fin.getline(dummy, NCHAR);        //The relative permeability of component 2 is 2.143594e+02, 1.349168e-02, -5.743910e-05
	                fin.getline(dummy, NCHAR);            //The LOCAL relative permeability of component 1 is -1.108770e+01, -8.600400e-03, -3.606825e-08
	                fin.getline(dummy, NCHAR);        //The LOCAL relative permeability of component 2 is 2.432098e+02, 1.871785e-02, -4.650501e-06
	                 fin>> chd>>chd>>chd >> chd>>var1;fin.getline(dummy, NCHAR);    //  Satuation of Component 1: 8.510638e-02, The satuation of Component 2: 9.148936e-0
	                fin.getline(dummy, NCHAR);        //The relative error of permiability computing is: 3.423860e-03
	                fin.getline(dummy, NCHAR);            //Elapsed time is 0h3m22s
	                fin.getline(dummy, NCHAR);        //The expected completion time is 1h46m20s

	                ofstream out(poreFileName,ios::app);
	                out<<var1<<" "<<var2<<" "<<var3<<endl;
	                out.close();
	                
	                fin.close();
	                
	                
	            
	                
	        }
	
}

}


