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

int total_number=4;        
 char poreFileName[128]="INPUT_128";
 ostringstream name;
for (int i=0;i<total_number;i++)
{
        
        name.str("");
	name<<poreFileName<<"."<<i<<".dat";
	ofstream out;
	out.open(name.str().c_str());
	
	out<<"128.all.dat                 :Geometry"<<endl;
        out<<"128 130 130                :nx ny nz"<<endl;
        out<<"60	                 :Maxmum time step"<<endl;
        out<<"10.0	                :x=1 (um) Resolution"<<endl;
        out<<"0	        	       :Pressure/Velocity Boundary "<<endl;
        out<<"0 0 0                           :mirroring option for x,y,z"<<endl;
        out<<"1.0e-5 0.0 0.0          :body force for x,y,z"<<endl;
        out<<"0 1.0 0 1.0144 	:Pressure Boundary in X direction"<<endl;
        out<<"0 1.0 0 1.0	:Pressure Boundary in Y direction"<<endl;
        out<<"0 1.0 0 1.0	:Pressure Boundary in Z direction"<<endl;
        out<<"0 0.0 0 0.07	:Velocity Boundary in X direction "<<endl;
        out<<"0 0.0 0 0.07	:Velocity Boundary in Y direction"<<endl;
        out<<"0 0.0 0 0.07	:Velocity Boundary in Z direction"<<endl;
        out<<"0.1            	:Viscosity"<<endl;
        out<<"0.0 0.0 0.0	:initial velocity for x,y,z"<<endl;
        out<<"=========OUTPUT==CONTROL==================="<<endl;
        out<<"1		:Permeability writing (1=yes, 0=no)"<<endl;
        out<<"1               :Direction of Permeability Calculation (1=X, 2=Y, 3=Z)"<<endl;
        out<<"10		:Frequency of results writing (interval in time steps)"<<endl;
        out<<"0               :Memory Saving mode for output subroutains "<<endl;
        out<<"-1		:Frequency of velocity field writing "<<endl;
        out<<"-1		:Frequency of density field writing"<<endl;
        out<<"0		:Mirroring part output (0=No,1=YES)"<<endl;
        out<<"==============ADVANCE==PARAMETER=========="<<endl;
        out<<"0 30 32 3 	:ONLY USE PART OF THE ORIGINAL GEOMETRY"<<endl;
        out<<"1		:FACTOR OF MESHREFINEMENT"<<endl;
        out<<"1 0.25 1.0	:Self define lattice velocity"<<endl;
        out<<"./128.		:OUTPUT PATH,DEFAULT"<<endl;
        out<<"1              :PRESSURE AND VELOCITY BOUNDARY CONDITION"<<endl;
        out<<"-1		:BACKUP FREQUENCY (-1=NO BACKUP)"<<endl;
        out<<"0               :USE BACKUP DATA (0=OFF, 1=ON)"<<endl;
        out<<"LBM_Backup_Density_4000            :BACKUP DATA FOR INTIALIZATION--DENSITY"<<endl;
        out<<"LBM_Backup_Velocity_4000           :BACKUP DATA FOR INTIALIZATION--Velocity"<<endl;
        out<<"LBM_Backup_f_4000                        :BACKUP DATA FOR DISTRIBUTION FUNCTION F"<<endl;
        out<<i<<endl;
             
        
}
}

