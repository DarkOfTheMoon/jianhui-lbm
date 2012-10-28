#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<cstring>

//=================
#include<stdio.h>
#include<netdb.h>
#include <time.h> 
//=================


using namespace std; 
      
int main (int argc , char * argv [])
{
char   name[128]; 
hostent*   pHost; 
gethostname(name,   128);
int date_time[6];//year,month,day,hour,minutes,second


cout<<"host name is: "<<name<<endl;


time_t rawtime;  
struct tm * timeinfo;  

time ( &rawtime );  
timeinfo = localtime ( &rawtime );  
printf ( "The current time is: %d:%d:%d", timeinfo->tm_hour,timeinfo->tm_min,timeinfo->tm_sec );  
cout<<endl;
printf("today is %4d.%02d.%02d  %02d:%02d:%02d\n", timeinfo->tm_year+1900, 
timeinfo->tm_mon+1, timeinfo->tm_mday, timeinfo->tm_hour, timeinfo->tm_min, timeinfo->tm_sec); 
cout<<endl;

date_time[0]=timeinfo->tm_year+1900;
date_time[1]=timeinfo->tm_mon+1;
date_time[2]=timeinfo->tm_mday;
date_time[3]=timeinfo->tm_hour;
date_time[4]=timeinfo->tm_min;
date_time[5]=timeinfo->tm_sec;





//====================================
/*
char a[128]="";
char s1[128] = "Hello";
strcpy(a,"Hello");

cout<<strcmp(a,s1)<<endl;

//cout<<sizeof(name)<<endl;
*/
//====================================


char output1[128]="licence1.dat";
char output2[128]="licence2.dat";

        ostringstream name1;
	name1<<output1;
	
	ofstream out;
	out.open(name1.str().c_str());
	out.write((char *)(&name), sizeof(char)*128); 
	out.close();

	
	
	name1.str("");
	name1<<output2;
	out.open(name1.str().c_str());
	out.write((char *)(&date_time[0]), sizeof(int)*6); 
	out.close();

	
	
	//------------receive and verify------------------------
	int ver_date_time[6];
	char ver_host_name[128];
	
	fstream fin;
	fin.open("licence1.dat",ios::in);
	if (fin.fail())
	        {
	        cout<<"\n file open error on licence1.dat"<<endl;
	        exit(-1);
	        }
	
	fin.read((char *)(ver_host_name), sizeof(char)*128);
	fin.close();
	
	fin.open("licence2.dat",ios::in);
	if (fin.fail())
	        {
	        cout<<"\n file open error on licence2.dat"<<endl;
	        exit(-1);
	        }
	
	fin.read((char *)(ver_date_time), sizeof(int)*6);
	fin.close();
	
	if (strcmp(name,ver_host_name)==0)
	        cout<<"Licence check OK"<<endl;
	else
	        cout<<"unvalid licence file, please check"<<endl;
	
	
	cout<<ver_date_time[0]<<endl;
	
	
	

}


