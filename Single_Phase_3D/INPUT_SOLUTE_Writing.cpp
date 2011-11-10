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



int total_number=10;
double g[15]={1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,2e-6,3e-6,4e-6,5e-6,6e-6,7e-6,8e-6,9e-6,1e-5};
double diffu[15]={0.05,0.02,0.01,0.005,0.002,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};
        
 char poreFileName[128]="INPUT_Solute_";
 ostringstream name;
for (int i=0;i<total_number;i++)
{
        
        name.str("");
	name<<poreFileName<<"."<<i<<".inputdat";
	ofstream out;
	out.open(name.str().c_str());
	
	out<<"=============INPUT FILE FOR 3D LATTICE BOLTZMANN MPI CODE--Solute and heat transfer======="<<endl;
	out<<"128.all.dat                 :Geometry"<<endl;
	out<<"phase.dat       	    	:Initial components distribution"<<endl;
        out<<"128 130 130                :nx ny nz"<<endl;
        out<<"60	                 :Maxmum time step"<<endl;
        out<<"10.0	                :x=1 (um) Resolution"<<endl;
        out<<"0	        	       :Pressure/Velocity Boundary "<<endl;
       	out<<"0.0 20	:Buoyancy Parameter (0= no buoyancy, reference_psi F=/rho*g+Bouyancy*(/psi-reference_psi)*g)"<<endl;
	out<<"0.0e-6			:Gravity (Only for Buoyancy, gravity*thermal expansion)"<<endl;
        out<<g[i]<<" 0.0 0.0          :body force for x,y,z"<<endl;
        out<<"0 1.0 0 1.0144 	:Pressure Boundary in X direction"<<endl;
        out<<"0 1.0 0 1.0	:Pressure Boundary in Y direction"<<endl;
        out<<"0 1.0 0 1.0	:Pressure Boundary in Z direction"<<endl;
        out<<"0 0.0 0 0.07	:Velocity Boundary in X direction "<<endl;
        out<<"0 0.0 0 0.07	:Velocity Boundary in Y direction"<<endl;
        out<<"0 0.0 0 0.07	:Velocity Boundary in Z direction"<<endl;
	out<<"0 10.0 0 30.0	:Constant concentration BC in X direction ()"<<endl;
	out<<"0 10.0 0 30.0	    		:Constant concentration BC in Y direction"<<endl;
	out<<"0 1.0 0 1.0	    		:Constant concentration BC in Z direction"<<endl;
	out<<"0 0.0 0 0.0			:Constant Diffusive Flux BC in X direction (p,n)"<<endl;
	out<<"0 0.0 0 0.0			:Constant Diffusive Flux BC in Y direction"<<endl;
	out<<"0 0.0 0 0.0			:Constant Diffusive Flux BC in Z direction"<<endl;
        out<<"0.02            	:Viscosity"<<endl;
	out<<diffu[i]<<"  		        :Diffusion coefficient"<<endl;
        out<<"0.0 0.0 0.0	:initial velocity for x,y,z"<<endl;
        out<<"=========OUTPUT==CONTROL==================="<<endl;
        out<<"1		:Permeability writing (1=yes, 0=no)"<<endl;
        out<<"1               :Direction of Permeability Calculation (1=X, 2=Y, 3=Z)"<<endl;
        out<<"10		:Frequency of results writing (interval in time steps)"<<endl;
        out<<"0               :Memory Saving mode for output subroutains "<<endl;
        out<<"-1		:Frequency of velocity field writing "<<endl;
        out<<"-1		:Frequency of density field writing"<<endl;
	out<<"-1				:Frequency of concentration writing "<<endl;
	out<<"0 100           :Dispersion statistics (0=No,1=X direction, 2=Y, 3=Z;Frequncy of writing)"<<endl;
        out<<"==============ADVANCE==PARAMETER=========="<<endl;
        out<<"0 0.25 1.0	:Self define lattice velocity"<<endl;
        out<<"./scr_"<<i<<".		:OUTPUT PATH,DEFAULT"<<endl;
        out<<"1              :PRESSURE AND VELOCITY BOUNDARY CONDITION"<<endl;
	out<<"1			:CONSTANT PSI AND FLUX BOUNDARY CONDITION OPTIONS: 0,1:Inamuro,GUO"<<endl;
	out<<"2000                 :Concentration initialization from the xth timesteps (-1=OFF)"<<endl;
	out<<"0 0 0       :Pemeability calculation Partially  "<<endl;
	out<<"90 20				:Permeability calculation partially"<<endl;
	out<<"0 100				:Permeability calculation partially "<<endl;
	out<<"0 100				:Permeability calculation partially"<<endl;
	out<<"====================BACKUP=PARAMETERS=========================="<<endl;
        out<<"-1		:BACKUP FREQUENCY (-1=NO BACKUP)"<<endl;
        out<<"0               :USE BACKUP DATA (0=OFF, 1=ON)"<<endl;
        out<<"LBM_Backup_Density_4000            :BACKUP DATA FOR INTIALIZATION--DENSITY"<<endl;
        out<<"LBM_Backup_Velocity_4000           :BACKUP DATA FOR INTIALIZATION--Velocity"<<endl;
	out<<"LBM_Backup_Concentration_500           :BACKUP DATA FOR INTIALIZATION--Velocity"<<endl;
        out<<"LBM_Backup_f_4000                        :BACKUP DATA FOR DISTRIBUTION FUNCTION F"<<endl;
	out<<"LBM_Backup_g_500                        :BACKUP DATA FOR DISTRIBUTION FUNCTION fg"<<endl;
        out<<i<<endl;
             
        
}
}

