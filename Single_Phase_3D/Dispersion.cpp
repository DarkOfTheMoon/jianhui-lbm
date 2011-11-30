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
int par_num=16*3*100;
int max_n=120000;
double delta_t=0.5;
//const double tor=1e-7;

bool*** Solid;
int** Dec;
double**** V;
double** x;
double t;
double d_m=0.0013;
double* psi_s;
double*** psi_g;

int main(int argc , char *argv [])
{     
        
double phi,theta;        
srand(time(NULL));
const double PI = std::atan(1.0)*4;



Solid= new bool**[nx];
psi_g= new double**[nx];
Dec = new int*[par_num];
V= new double***[nx];
x=new double*[par_num];
psi_s= new double[nx];

for (int i=0;i<par_num;i++)
        {
        x[i]= new double[3];
        Dec[i]= new int [3];
        
        }
for (int i=0;i<nx;i++)
        {
        Solid[i]= new bool*[ny];
      psi_g[i]= new double*[ny];
        V[i]= new double**[ny];
                for (int j=0;j<ny;j++)
                {
                        Solid[i][j]=new bool[nz];
                        psi_g[i][j]= new double[nz];
                        V[i][j]=new  double*[nz];
                                for (int k=0;k<nz;k++)
                                {
                                        V[i][j][k]=new double[3];
                                   
                                }
                }
        }
        
        
        
for (int is=0;is<par_num;is++)
                 psi_s[is]=0.0;
        
        for (int i=0;i<nx;i++)
                for (int j=0;j<ny;j++)
                for (int k=0;k<nz;k++)
                {
                        V[i][j][k][0]=1.0e-4;
                         V[i][j][k][1]=0.0;
                        V[i][j][k][2]=0.0;
                        
                        
                if ((j==0) or (j==ny-1))
                        Solid[i][j][k]=0;
                else
                        Solid[i][j][k]=0;
                
                }
                
  //==========INITIAL X===============     
        for (int j=0;j<ny;j++)
                for (int k=0;k<nz;k++)
                for (int m=0;m<100;m++)
                        {   
                                Dec[j*3*100+k*100+m][0]=0.0;
                                Dec[j*3*100+k*100+m][1]=0.0;
                                Dec[j*3*100+k*100+m][2]=0.0;
                                
                                 x[j*3*100+k*100+m][0]=(double)50+ (double(rand()%10000))/10000*0.5;
                                 x[j*3*100+k*100+m][1]=(double)j+ (double(rand()%10000))/10000*0.5;
                                 x[j*3*100+k*100+m][2]=(double)k+ (double(rand()%10000))/10000*0.5;
                                 //cout<<j*3*3+k*3+m<<endl;
                        }
  //==========================
  
  
             double lsx,lsy,lsz,dtx,dty,dtz,dt_l,llsx,llsy,llsz;
             int lx,ly,lz,lxp,lyp,lzp,dec0,dec1,dec2;
             for (int n=0;n<max_n;n++)
                 {
                 t=n*delta_t;
                 for (int pn=0;pn<par_num;pn++)
                         {
                         lx=(int)round(x[pn][0]);
                         ly=(int)round(x[pn][1]); 
                         lz=(int)round(x[pn][2]);
                          //lsx=x[pn][0];lsy=x[pn][1];lsz=x[pn][2];
                          
                        // if ((abs(x[pn][0]-lx)<tor) or (abs(x[pn][1]-ly)<tor) or (abs(x[pn][2]-lz)<tor))
                        
                                         phi=(double(rand()%10000))/10000*PI;
                                         theta=(double(rand()%10000))/10000*2*PI;
                                         lsx=x[pn][0]+delta_t*V[lx][ly][lz][0]+sqrt(6*d_m*delta_t)*sin(phi)*cos(theta);
                                         lsy=x[pn][1]+delta_t*V[lx][ly][lz][1]+sqrt(6*d_m*delta_t)*sin(phi)*sin(theta);
                                         lsz=x[pn][2]+delta_t*V[lx][ly][lz][2]+sqrt(6*d_m*delta_t)*cos(theta);
                                         
                                         lxp=(int)round(lsx);
                                         lyp=(int)round(lsy);
                                         lzp=(int)round(lsz);
                                         
                                   
                           if ((lxp!=lx) or (lyp!=ly) or (lzp!=lz))
                                   {
                                           dec0=0;dec1=0;dec2=0;
                                           if (lxp>=nx) {lxp=0;dec0=1;lsx-=(double)nx;}
                                           if (lxp<0) {lxp=nx-1;dec0=-1;lsx+=(double)nx;}
                                           if (lyp>=ny) {lyp=0;dec1=1;lsy-=(double)ny;}
                                           if (lyp<0) {lyp=ny-1;dec1=-1;lsy+=(double)ny;}
                                           if (lzp>=nz) {lzp=0;dec2=1;lsz-=(double)nz;}
                                           if (lzp<0) {lzp=nz-1;dec2=-1;lsz+=(double)nz;}
                                           
                                           if (Solid[lxp][lyp][lzp]==1)
                                                   {
                                                           dtx=delta_t;dty=delta_t;dtz=delta_t;
                                                   if (lxp+dec0*nx>lx)
                                                   {
                                                           if (abs(V[lx][ly][lz][0])<1e-30) 
                                                           dtx=0.0;
                                                                else
                                                           dtx=abs(((double)lx+0.5-lsx)/V[lx][ly][lz][0]);
                                                           
                                                   }
                                                   
                                                   if (lyp+dec1*ny>ly)
                                                   {
                                                           if (abs(V[lx][ly][lz][1])<1e-30) 
                                                           dty=0.0;
                                                                else
                                                           dty=abs(((double)ly+0.5-lsy)/V[lx][ly][lz][1]);
                                                   }
                                                   
                                                   
                                                   if (lzp+dec2*nz>lz)
                                                   {
                                                           if (abs(V[lx][ly][lz][2])<1e-30) 
                                                           dtz=0.0;
                                                                else
                                                           dtz=abs(((double)lz+0.5-lsz)/V[lx][ly][lz][2]);
                                                           
                                                           
                                                           
                                                   }
                                                   if (lxp+dec0*nx<lx) 
                                                   {
                                                           if (abs(V[lx][ly][lz][0])<1e-30) 
                                                           dtx=0.0;
                                                                else
                                                         dtx=abs(((double)lx-0.5-lsx)/V[lx][ly][lz][0]);
                                                   }
                                                   
                                                   
                                                   if (lyp+dec1*ny<ly)
                                                   {
                                                           if (abs(V[lx][ly][lz][1])<1e-30) 
                                                           dty=0.0;
                                                                else
                                                           dty=abs(((double)ly-0.5-lsy)/V[lx][ly][lz][1]);
                                                   }
                                                   
                                                   if (lzp+dec2*nz<lz)
                                                   {
                                                           if (abs(V[lx][ly][lz][2])<1e-30) 
                                                           dtz=0.0;
                                                                else
                                                           dtz=abs(((double)lz-0.5-lsz)/V[lx][ly][lz][2]);
                                                   }
                                                   
                                                   
                                                 
                                                   
                                                    dt_l=dtx<dty?dtx:dty;
                                                    dt_l=dt_l<dtz?dt_l:dtz;
                                                    
                                                    lsx=x[pn][0]+dt_l*V[lx][ly][lz][0];
                                                    lsy=x[pn][1]+dt_l*V[lx][ly][lz][1];
                                                    lsz=x[pn][2]+dt_l*V[lx][ly][lz][2];
                                                    
                                                    do 
                                                            {
                                                                    phi=(double(rand()%10000))/10000*PI;
                                                                    theta=(double(rand()%10000))/10000*2*PI;
                                                            llsx=lsx+sqrt(6*d_m*(delta_t-dt_l))*sin(phi)*cos(theta);
                                                            llsy=lsy+sqrt(6*d_m*(delta_t-dt_l))*sin(phi)*sin(theta);
                                                            llsz=lsz+sqrt(6*d_m*(delta_t-dt_l))*cos(theta);
                                                            }while(Solid[(int)round(lsx)][(int)round(lsy)][(int)round(lsz)]==1);
                                                    
                                                    
                                                   }
                                                   else
                                                           {
                                                           x[pn][0]=lsx;
                                                           x[pn][1]=lsy;
                                                           x[pn][2]=lsz;
                                                           
                                                           Dec[pn][0]+=dec0;
                                                           Dec[pn][1]+=dec1;
                                                           Dec[pn][2]+=dec2;
                                                           }
                                           
                                          
                                   
                                   
                                   
                                   
                                   
                                   
                                   }
                                   else
                                           {  
                                          //  dec0=0;dec1=0;dec2=0;       
                                           //if (lxp>=nx) {lxp=0;dec0=1;lsx-=(double)nx;}
                                           //if (lxp<0) {lxp=nx-1;dec0=-1;lsx+=(double)nx;}
                                           //if (lyp>=ny) {lyp=0;dec1=1;lsy-=(double)ny;}
                                           //if (lyp<0) {lyp=ny-1;dec1=-1;lsy+=(double)ny;}
                                           //if (lzp>=nz) {lzp=0;dec2=1;lsz-=(double)nz;}
                                           //if (lzp<0) {lzp=nz-1;dec2=-1;lsz+=(double)nz;}
                                           
                                          //                      Dec[pn][0]+=dec0;
                                          //                 Dec[pn][1]+=dec1;
                                          //                 Dec[pn][2]+=dec2;
                                           x[pn][0]=lsx;
                                           x[pn][1]=lsy;
                                           x[pn][2]=lsz;
                                           
                                           
                                           
                                           
                                           }
                                         
                                         
                                         
                                         
                                         
                                      
                         
                         
                         
                         }
                 
                 
                 //cout<<x[1][0]<<"  "<<x[1][1]<<"  "<<x[1][2]<<"  "<<n<<endl;
                 }
               
               //cout<<x[1][0]<<"  "<<x[1][1]<<"  "<<x[1][2]<<"   ffff"<<endl;  
 
        // cout<<abs(-3.3)<<"   "<<(int)round(-0.3)<<"   "<<sin(0.25*PI)<<endl;
         cout<<PI<<endl;
         
 cout<<x[1][0]<<"  "<<x[1][1]<<"  "<<x[1][2]<<"   ffff"<<endl;
        
         //cout<<x[1][0]<<"  "<<x[1][1]<<"  "<<x[1][2]<<"   gggggggggf"<<endl; 
         for (int i=0;i<nx;i++)
                   for (int j=0;j<ny;j++)
                for (int k=0;k<nz;k++)
                  psi_g[i][j][k]=0.0;
          
        // cout<<x[1][0]<<"  "<<x[1][1]<<"  "<<x[1][2]<<"   ggggf"<<endl; 
         char FileName[128]="";//strcpy(FileName,outputfile);
         char FileName2[128]="";//strcpy(FileName2,outputfile);
         char FileName3[128]="";//strcpy(FileName3,outputfile);
         strcat(FileName,"Psi.txt");
         strcat(FileName2,"psi_general.vtk");
         strcat(FileName3,"location_particles.txt");
         
         //cout<<x[1][0]<<"  "<<x[1][1]<<"  "<<x[1][2]<<"aaaa"<<endl;
          ofstream fin(FileName,ios::out);       
        for (int i=0;i<par_num;i++)
                {
                        lx=(int)round(x[i][0]);
                        ly=(int)round(x[i][1]);
                        lz=(int)round(x[i][2]);
                psi_s[lx]+=1.0;
                psi_g[lx][ly][lz]+=1.0;
                }
                
              for (int i=0;i<nx;i++)
                      fin<<psi_s[i]<<endl;
             fin.close();
             
             
             //cout<<x[1][0]<<"  "<<x[1][1]<<"  "<<x[1][2]<<endl;
             
             ofstream out(FileName2,ios::out); 
             
             out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<nx<<"         "<<ny<<"         "<<nz<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<nx*ny*nz<<endl;
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;

            for (int i=0;i<nx;i++)
                   for (int j=0;j<ny;j++)
                for (int k=0;k<nz;k++)
                out<<psi_g[i][j][k]<<endl;
        out.close();
         
        ofstream out2(FileName3,ios::out); 
         for (int i=0;i<par_num;i++)
                 out2<<x[i][0]<<"   "<<x[i][1]<<"   "<<x[i][2]<<endl;
         
         out2.close();
                 
         
}


