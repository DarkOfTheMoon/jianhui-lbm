#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>

#include<stdio.h>
#include<netdb.h>
#include <time.h> 



using namespace std; 
      
int main (int argc , char * argv [])
{
char   name[128]; 
hostent*   pHost; 
gethostname(name,   128);

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

}


