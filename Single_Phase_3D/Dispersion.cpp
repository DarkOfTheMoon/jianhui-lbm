#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<math.h> 
#define _USE_MATH_DEFINES

using namespace std;   

int nx=100;
int ny=16;
int nz=3;
int par_num=16*3*3;
int max_n=1000;
double delta_t=0.5;
const double tor=1e-7;

bool*** Solid;
double**** V;
double** x;
double t;
double d_m=0.0013;


int main(int argc , char *argv [])
{     
        
double phi,theta;        
srand(time(NULL));
const double PI = std::atan(1.0)*4;



Solid= new bool**[nx];
V= new double***[nx];
x=new double*[par_num];

for (int i=0;i<par_num;i++)
        {
        x[i]= new double[3];
        
        }
for (int i=0;i<nx;i++)
        {
        Solid[i]= new bool*[ny];
        V[i]= new double**[ny];
                for (int j=0;j<ny;j++)
                {
                        Solid[i][j]=new bool[nz];
                        V[i][j]=new  double*[nz];
                                for (int k=0;k<nz;k++)
                                        V[i][j][k]=new double[3];
                }
        }
        
        
        

        
        for (int i=0;i<nx;i++)
                for (int j=0;j<ny;j++)
                for (int k=0;k<nz;k++)
                {
                        V[i][j][k][0]=0.0;
                         V[i][j][k][1]=0.0;
                        V[i][j][k][2]=0.0;
                        
                        
                if ((j==0) or (j==ny))
                        Solid[i][j][k]=0;
                else
                        Solid[i][j][k]=0;
                
                }
                
       
        for (int j=0;j<ny;j++)
                for (int k=0;k<nz;k++)
                for (int m=0;m<3;m++)
                        {
                                 x[j*3*3+k*3+m][0]=50+ (double(rand()%10000))/10000*0.5;
                                 x[j*3*3+k*3+m][1]=j+ (double(rand()%10000))/10000*0.5;
                                 x[j*3*3+k*3+m][2]=k+ (double(rand()%10000))/10000*0.5;
                        
                        }
                        
             double lsx,lsy,lsz;
             int lx,ly,lz;
         for (int n=0;n<max_n;n++)
                 {
                 t=n*delta_t;
                 for (int pn=0;pn<par_num;pn++)
                         {
                         lx=round(x[pn][0]);
                          ly=round(x[pn][1]); 
                          lz=round(x[pn][2]);
                          lsx=x[pn][0];lsy=x[pn][1];lsz=x[pn][2];
                          
                         if ((abs(x[pn][0]-lx)<tor) or (abs(x[pn][1]-ly)<tor) or (abs(x[pn][2]-lz)<tor))
                                 {
                                         do
                                         {
                                         phi=(double(rand()%10000))/10000*PI;
                                         theta=(double(rand()%10000))/10000*2*PI;
                                         lsx=x[pn][0]+sqrt(6*d_m*delta_t)*sin(phi)*cos(theta);
                                         lsy=x[pn][1]+sqrt(6*d_m*delta_t)*sin(phi)*sin(theta);
                                         lsz=x[pn][2]+sqrt(6*d_m*delta_t)*cos(theta);
                                         
                                         }while (!Solid[(int)round(lsx)][(int)round(lsy)][(int)round(lsz)]);
                                        
                                        x[pn][0]=lsx;x[pn][1]=lsy;x[pn][2]=lsz;
                                 
                                 }
                                 else
                                         {
                                         phi=(double(rand()%10000))/10000*PI;
                                         theta=(double(rand()%10000))/10000*2*PI;
                                         lsx=x[pn][0]+delta_t*V[lx][ly][lz][0]+sqrt(6*d_m*delta_t)*sin(phi)*cos(theta);
                                         lsy=x[pn][1]+delta_t*V[lx][ly][lz][1]+sqrt(6*d_m*delta_t)*sin(phi)*sin(theta);
                                         lsz=x[pn][2]+delta_t*V[lx][ly][lz][2]+sqrt(6*d_m*delta_t)*cos(theta);
                                         
                                         
                                         
                                         
                                         
                                         
                                         
                                         
                                         }
                         
                         
                         
                         }
                 
                 
                 
                 }
               
                
 
         cout<<abs(-3.3)<<"   "<<round(3.4999999)<<"   "<<sin(0.25*PI)<<endl;
         cout<<PI<<endl;
         
         
}


