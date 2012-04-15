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
char prefix[128]="rel_";    //prefix<<i<<"_Results.txt"
char poreFileName[128]="Summary.outdat";
char poreFileName2[128]="Summary_Least_Square.outdat";


float var1,var2,var3,var4,var5;	
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
	
	fins.open(poreFileName2,ios::trunc);
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
	                fin>> chd >> chd>>chd>>chd >> chd>>chd>>chd>>var2;fin.getline(dummy, NCHAR);        //   The relative permeability of component 1 is -3.992098e+00, -3.079636e-03, 1.535446e-07
	                fin >>chd >> chd>>chd>>chd >> chd>>chd>>chd>> var3;fin.getline(dummy, NCHAR);         //   The relative permeability of component 2 is -3.992098e+00, -3.079636e-03, 1.535446e-07
	                fin>> chd >> chd>>chd>>chd >> chd>>chd>>chd>>chd>>chd>>var4;fin.getline(dummy, NCHAR);            //The LEAST SQUARED relative permeability of component 1 is 1.594603e-01
	                fin >>chd >> chd>>chd>>chd >> chd>>chd>>chd>> chd>>chd>>var5;fin.getline(dummy, NCHAR);        //The LEAST SQUARED relative permeability of component 2 is 4.766752e-01
	                fin>> chd>>chd>>chd >> chd>>var1;fin.getline(dummy, NCHAR);    //  Satuation of Component 1: 8.510638e-02, The satuation of Component 2: 9.148936e-0
	                fin.getline(dummy, NCHAR);        //The relative error of permiability computing is: 3.423860e-03
	                fin.getline(dummy, NCHAR);            //Elapsed time is 0h3m22s
	                fin.getline(dummy, NCHAR);        //The expected completion time is 1h46m20s
			fin.getline(dummy, NCHAR);
			
	                fin.getline(dummy, NCHAR);        //       The30800th computation result:     
	                fin.getline(dummy, NCHAR);        //The Maximum velocity is: 0.0120153   Re_l=1.52508   Re_g=1.52508
	                fin.getline(dummy, NCHAR);           // Courant Number=0.0120153	 Capillary Num=0.0076254
	                fin.getline(dummy, NCHAR);        //The max relative error of velocity is: 7.971677e-01
	                fin.getline(dummy, NCHAR);         //   The relative permeability of component 1 is -3.992098e+00, -3.079636e-03, 1.535446e-07
	                fin.getline(dummy, NCHAR);        //The relative permeability of component 2 is 2.143594e+02, 1.349168e-02, -5.743910e-05
	                fin.getline(dummy, NCHAR);            //The LEAST SQUARED relative permeability of component 1 is 1.594603e-01
	                fin.getline(dummy, NCHAR);        //The LEAST SQUARED relative permeability of component 2 is 4.766752e-01
	                fin.getline(dummy, NCHAR);    //  Satuation of Component 1: 8.510638e-02, The satuation of Component 2: 9.148936e-0
	                fin.getline(dummy, NCHAR);        //The relative error of permiability computing is: 3.423860e-03
	                fin.getline(dummy, NCHAR);            //Elapsed time is 0h3m22s
	                fin.getline(dummy, NCHAR);        //The expected completion time is 1h46m20s

			
	                ofstream out(poreFileName,ios::app);
	                out<<var1<<" "<<var2<<" "<<var3<<endl;
	                out.close();
	                
	            
	                
	                ofstream out2(poreFileName2,ios::app);
	                out2<<var1<<" "<<var4<<" "<<var5<<endl;
	                out2.close();
	                
	                fin.close();
	            
	                
	        }
	
}

}


