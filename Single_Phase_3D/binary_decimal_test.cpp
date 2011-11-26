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
        ofstream fs("binary_file.bin");
        
        double i[3]={3.2567,3.4578422,7.434665};
      
      fs.write((char *)(&i),sizeof(i)*3);
     // cout<<sizeof(i)<<endl;
      fs.close();
      
      double j[3];
      fstream file;
      
      file.open("binary_file.bin",ios_base::in);
      
      file.read((char *)(&j),sizeof(i)*3);
      
      cout<<setiosflags(ios::scientific)<<j[0]<<"  "<<j[1]<<"  "<<j[2]<<endl;
      cout<<setprecision(10)<<j[0]<<"  "<<j[1]<<"  "<<j[2]<<endl;
      file.close();
      
      
        double num1=267.435346;
        
        char* a;
        
        a=(char *)(&num1);
        
        //cout<<a<<endl;
        
        double* num2=(double *)(a);
        
        cout<<*num2*3<<"   "<<*num2<<endl;
}



