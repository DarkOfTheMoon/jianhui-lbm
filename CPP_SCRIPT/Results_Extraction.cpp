/*

TO USE THE CODE
g++ Results_Extraction.cpp -o test
./test rel_ 2000

then the code will extract datas from rel_1_Reslts.txt, rel_2_results.txt ..... rel_10_Results.txt
from the most recent 2000 time step


*/



#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string.h>



using namespace std; 
      
int main (int argc , char * argv [])
{
const int total_number=3;
int indexs[total_number]={2,3,5};


int i;

const int imax=300000;
int num;
       int le_num;
       if (argc>2)
               le_num=atoi(argv[2]);
       else
               le_num=1500;

char prefix[128];    //prefix<<i<<"_Results.txt"
strcpy(prefix,argv[1]);

char poreFileName[128]="Summary.outdat";
char poreFileNameles[128]="Summary_Least_Square.outdat";
char poreFileNameave[128]="Summary_Ave.outdat";

char rel1file[128];
char rel2file[128];

float rel1[imax];
float rel2[imax];
char ls[128];
char ls2[128];
char index[128];
       

float var1,var2,var3,var4,var5;	
float ave1,ave2;

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
	
	fins.open(poreFileNameles,ios::trunc);
	fins.close();
	
	fins.open(poreFileNameave,ios::trunc);
	fins.close();
	
for (int ige=0;ige<total_number;ige++)
{
        i=indexs[ige];
        //cout<<i<<endl;
        
        
        name.str("");
	name<<prefix<<i<<"_Results.txt";
	
	ftest = fopen(name.str().c_str(), "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file (" << name.str().c_str() <<
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

			/*
	                ofstream out(poreFileName,ios::app);
	                out<<var1<<" "<<var2<<" "<<var3<<endl;
	                out.close();
	                
	            
	                
	                ofstream out2(poreFileNameles,ios::app);
	                out2<<var1<<" "<<var4<<" "<<var5<<endl;
	                out2.close();
	                */
	                fin.close();
	            
	                
	        }
	        
	        
	        
	        
	 name.str("");
	name<<prefix<<i<<"_Relative_Permeability_Component1.txt";
	
	ftest = fopen(name.str().c_str(), "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file (" <<name.str().c_str() <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	
	                fclose(ftest);
	             strcpy(rel1file,prefix);
	             sprintf(index, "%d", i);
	             strcat(rel1file,index);
	             strcat(rel1file,"_Relative_Permeability_Component1.txt");
	             fin.open(rel1file);
	
	
	num=0;
	while ((!fin.eof()) and (num<imax))
		{	
			fin >>ls; fin.getline(dummy, NCHAR);
			rel1[num] =   atof(ls);
			//cout<<num<<"   @     "<<rel1[num]<<" "<<ls<<" "<<rel1file<<endl;
			
			num++;
			
		}
	fin.close();   
	
	        
	        
	        
	        //=================================
	name.str("");
	name<<prefix<<i<<"_Relative_Permeability_Component2.txt";
	
	ftest = fopen(name.str().c_str(), "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file (" << name.str().c_str() <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	
	                fclose(ftest);
	             strcpy(rel2file,prefix);
	             sprintf(index, "%d", i);
	             strcat(rel2file,index);
	             strcat(rel2file,"_Relative_Permeability_Component2.txt");
	          
	             fin.open(rel2file);
	
	
	num=0;
	while ((!fin.eof()) and (num<imax))
		{	
			fin >>ls; fin.getline(dummy, NCHAR);
			rel2[num] =   atof(ls);
			//cout<<num<<"        "<<rel1[num]<<"   "<<name.str().c_str()<<endl;
			num++;
			
		}
	fin.close();   
	       
	        
	        
	        
	        
	        double tb1,yb1,x11,x01,lst1,lsb1;
	  double tb2,yb2,x12,x02,lst2,lsb2;
	        ave1=0;ave2=0; tb1=0;yb1=0;tb2=0;yb2=0;
	for (int j=num-le_num;j<num;j++)
	        {
	                ave1+=rel1[j];ave2+=rel2[j];
	                tb1+=j;tb2+=j;
	                //cout<<rel1[j]<<endl;
	        }
	        tb1/=le_num;tb2/=le_num;
	        ave1/=le_num;ave2/=le_num;
	        yb1=ave1;yb2=ave2;
	        lst1=0;lsb1=0;lst2=0;lsb2=0;
	        for (int j=num-le_num;j<num;j++)
	        {
	                lst1+=(j-tb1)*(rel1[j]-yb1);
	                lsb1+=(j-tb1)*(j-tb1);
	                lst2+=(j-tb2)*(rel2[j]-yb2);
	                lsb2+=(j-tb2)*(j-tb2);
	                
	        }
	        
	                                x11=lst1/lsb1;
	                          x01=yb1-x11*tb1;
	                          x12=lst2/lsb2;
	                          x02=yb2-x12*tb2;
	                          
	                          var4=x01+x11*(num-1);
	                          var5=x02+x12*(num-1);
	                          
	               ofstream out(poreFileName,ios::app);
	                out<<var1<<" "<<var2<<" "<<var3<<endl;
	                out.close();
	                
	            
	                
	                ofstream out2(poreFileNameles,ios::app);
	                out2<<var1<<" "<<var4<<" "<<var5<<endl;
	                out2.close();
	                
	                
	                ofstream out3(poreFileNameave,ios::app);
	                out3<<var1<<" "<<ave1<<" "<<ave2<<endl;
	                out3.close();
	                
	                          
	                          
	        
	        
	
}

}


