#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>

using namespace std;       
int main(int argc , char *argv [])
{	
//FileStream   file1=new   FileStream("test_file"); 
        //ofstream fs("binary_file.bin");
        //fstream fs("binary_file.bin",ios::out|ios::binary|ios::app);
        fstream fs;
        fs.open("binary_file.bin",ios::out|ios::binary);
        double i[3]={3.2567,3.4578422,7.434665};
        double ia[3][2]={{2.3521,12.356674},{256.35321,1.22222},{3.1111,5.1111}};
        
        
      fs.write((char *)(ia),sizeof(double)*6);
      fs.close();
      
      
      //fs.open("binary_file.bin",ios::app|ios::binary);
      //fs.write((char *)(ia),sizeof(ia)*3);
      
      
     // cout<<sizeof(i)<<endl;
      //fs.close();
      
      double j[6];
      double ja[3][2];
      fstream file;
      
      file.open("binary_file.bin",ios::in);
      
      file.read((char *)(&ja[0][0]),sizeof(double)*6);
      //file.read((char *)(ja),sizeof(double)*3);
      
      cout<<setiosflags(ios::scientific)<<j[0]<<"  "<<j[1]<<"  "<<j[2]<<endl;
      cout<<endl;
      cout<<ja[0][0]<<"  "<<ja[1][1]<<"  "<<ja[2][0]<<endl;
      cout<<setprecision(10)<<j[3]<<"  "<<j[4]<<"  "<<j[5]<<endl;
      file.close();
      
      
        double num1=267.435346;
        
        char* a;
        //char* b="^-??p??9p?";
        char b[128];
        //b="^-??p??9p?";
        
        
        a=(char *)(&num1);
        memcpy (b,&num1,8);
        //cout<<a<<endl;
        cout<<"16 digital= "<<hex<<num1<<endl;
        
        double* num2=(double *)(a);
        
        cout<<*num2*3<<"   "<<*num2<<endl;
        
        sprintf(b,"%A\n",num1);
        
        double n=0X1.0B6F72D5E071CP+8;
        cout<<"n="<<n<<"         "<<b<<endl;
}



