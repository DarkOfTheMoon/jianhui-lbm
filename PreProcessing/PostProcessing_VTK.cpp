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
char VTKName[128]="LBM_Geometry.vtk";
char VTKName_head[128]="LBM_psi_";
char poreFileNameOut[128]="128_out.vtk";
double tor=1.0e-1;

double*** Solid;
double*** Solid2;


int nx,ny,nz;
ostringstream name;
ostringstream name2;
ofstream fins;	
char chd[128];
int NCHAR=128;	
char dummy[128+1];
int xp,xn,yp,yn,zp,zn;



	FILE *ftest;
	ifstream fin;
	name.str("");
	name<<VTKName;
	ftest = fopen(name.str().c_str(), "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file (" <<name.str().c_str()  <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest);
	
	name2.str("");
	name2<<poreFileNameOut;	
	fins.open(name2.str().c_str(),ios::out);
	
	fin.open(name.str().c_str());
	fin.getline(dummy, NCHAR);fins<<dummy<<endl;
	fin.getline(dummy, NCHAR);fins<<dummy<<endl;
	fin.getline(dummy, NCHAR);fins<<dummy<<endl;
	fin.getline(dummy, NCHAR);fins<<dummy<<endl;
	fin.getline(dummy, NCHAR);fins<<dummy<<endl;
	fin.getline(dummy, NCHAR);fins<<dummy<<endl;
	fin.getline(dummy, NCHAR);fins<<dummy<<endl;
	fin.getline(dummy, NCHAR);fins<<dummy<<endl;
	fin.getline(dummy, NCHAR);fins<<dummy<<endl;
	fin.getline(dummy, NCHAR);fins<<dummy<<endl;
	fin.close();
	
	fin.open(name.str().c_str());
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin>> chd >>nx>>ny>>nz;fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	fin.getline(dummy, NCHAR);
	
	
	
        Solid = new double**[nx];
        Solid2 = new double**[nx];
        for (int i=0;i<nx;i++)
        {
                Solid[i]=new double*[ny];
                Solid2[i]=new double*[ny];
                for (int j=0;j<ny;j++)
                        {
                        Solid[i][j]= new double[nz];
                        Solid2[i][j]= new double[nz];
                        }
        }
        
        
	
	for (int k=0;k<nz;k++)
	        for (int j=0;j<ny;j++)
	        for (int i=0;i<nx;i++)
	        fin>>Solid[i][j][k];
	
	for (int k=0;k<nz;k++)
	        for (int j=0;j<ny;j++)
	        for (int i=0;i<nx;i++)
	        {
	                xp=i+1;if (xp>=nx) xp=0;
	                xn=i-1;if (xn<0) xn=nx-1;
	                yp=j+1;if (yp>=ny) yp=0;
	                yn=j-1;if (yn<0) yn=ny-1;
	                zp=k+1;if (zp>=nz) zp=0;
	                zn=k-1;if (zn<0) zn=nz-1;
	                
	                Solid2[i][j][k]=Solid[xp][j][k]+Solid[xn][j][k]+Solid[i][yp][k]+Solid[i][yn][k]+Solid[i][j][zn]+Solid[i][j][zp]-6*Solid[i][j][k];
	                
	                if (Solid2[i][j][k]>=tor)
	                        Solid2[i][j][k]=1;
	                else
	                        //if (Solid2[i][j][k]<=-tor)
	                       // Solid2[i][j][k]=-1;
	                        //else
	                                Solid2[i][j][k]=0;
	                                
	                                
	                                
	                                
	        }
	
	
cout<<nx<<" "<<ny<<" "<<nz<<endl;


	
	
	
	
	
	
	
	
	for (int k=0;k<nz;k++)
	        for (int j=0;j<ny;j++)
	        for (int i=0;i<nx;i++)
	        fins<<Solid2[i][j][k]<<" ";
	fins.close();
	fin.close();
	
}


