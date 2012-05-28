/*

TO USE THE CODE:
compile the code with g++ Least...cpp -o test
run the code 
./test rel_34_ 5000
rel_34_ is the prefix, 5000 is the length of least square fitting

then use 
sh new_plot.sh rel_34_
to visulize the resutls

*/





#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string.h>

using namespace std; 
      
int main(int argc , char *argv [])
{
        
       const int imax=200000;
       int le_num;
       if (argc>2)
               le_num=atoi(argv[2]);
       else
               le_num=1500;
       
       
         char poreFileName1[128];
       // cout<<argv[1]<<endl;
        strcpy(poreFileName1,argv[1]);
        
        strcat(poreFileName1,"Relative_Permeability_Component1.txt");
        
         char poreFileName2[128];
         strcpy(poreFileName2,argv[1]);
        strcat(poreFileName2,"Relative_Permeability_Component2.txt");
        int NCHAR=128;char dummy[128+1];
        float rel1[imax];
        float rel2[imax];
        char ls[128];
        int num;
        

	FILE *ftest;
	ifstream fin;
	
	ftest = fopen(poreFileName1, "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file (" << poreFileName1 <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest);

	fin.open(poreFileName1);
	
	
	num=0;
	while ((!fin.eof()) and (num<imax))
		{	
			fin >>ls; fin.getline(dummy, NCHAR);
			rel1[num] =   atof(ls);
			//cout<<num<<"        "<<rel1[num]<<endl;
			num++;
			
		}
	fin.close();
	
	ftest = fopen(poreFileName2, "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file (" << poreFileName2 <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest);

	fin.open(poreFileName2);
	for (int i=0;i<num;i++)
	        {
	            fin >>ls; fin.getline(dummy, NCHAR);
			rel2[i] =   atof(ls);
			//cout<<num<<"        "<<rel1[num]<<endl;
			
	                
	        }
	        
	  double tb1,yb1,x11,x01,lst1,lsb1;
	  double tb2,yb2,x12,x02,lst2,lsb2;
	  double* rel_new1;
	  double* rel_new2;
	  
	  rel_new1= new double[num-le_num];
	  rel_new2=new double[num-le_num];
	  for (int i=0;i<num-le_num;i++)
	          {
	                  tb1=0;yb1=0;tb2=0;yb2=0;
	                  for (int j=i;j<i+le_num;j++)
	                          {tb1+=j;yb1+=rel1[j];tb2+=j;yb2+=rel2[j];}
	                  tb1/=le_num;yb1/=le_num;
	                  lst1=0;lsb1=0;lst2=0;lsb2=0;
	                   tb2/=le_num;yb2/=le_num;
	               
	                  for (int j=i;j<i+le_num;j++)
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
	                   rel_new1[i]=x01+x11*i;  
	                   rel_new2[i]=x02+x12*i;
	                          
	          }
	        
	
	ostringstream name;
	name<<argv[1]<<"new.dat";
	ofstream out;
	out.open(name.str().c_str());
	
	for (int i=0;i<num-le_num;i++)
	        out<<rel_new1[i]<<" "<<rel_new2[i]<<endl;
	
	
	out.close();
	
	 
	 
	 
	 
	 
	cout<<num<<"        "<<le_num<<endl;
	
	
	
		
		








}
