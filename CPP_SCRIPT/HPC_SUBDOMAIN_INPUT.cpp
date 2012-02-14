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



int total_number=11;
double g[15]={1e-6,1e-6,1e-6,1e-6,1e-6,1e-6,2e-6,3e-6,4e-6,5e-6,6e-6,7e-6,8e-6,9e-6,1e-5};
double diffu[15]={0.05,0.02,0.01,0.005,0.002,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001,0.001};

double var1[11];
        
 	char poreFileName[128]="INPUT_CG_";
	char poreFileName2[128]="LV60_";
 ostringstream name;
for (int i=0;i<total_number;i++)
{
        
        name.str("");
	name<<poreFileName<<i<<".inputdat";
	ofstream out;
	out.open(name.str().c_str());
	


	//===========================================SP==SOLUTE_HEAT==============================================
	/*
	out<<"=============INPUT FILE FOR 3D LATTICE BOLTZMANN MPI CODE--Solute and heat transfer======="<<endl;
	out<<"128.all.dat                 :Geometry"<<endl;
	out<<"phase.dat       	    	:Initial components distribution"<<endl;
        out<<"128 130 130                :nx ny nz"<<endl;
        out<<"60	                 :Maxmum time step"<<endl;
        out<<"10.0	                :x=1 (um) Resolution"<<endl;
        out<<"0	        	       :Pressure/Velocity Boundary "<<endl;
       	out<<"0.0 20	:Buoyancy Parameter (0=no buoyancy, reference_psi F=/rho*g+Bouyancy*(/psi-reference_psi)*g)"<<endl;
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
        out<<"500				:intial data for velocity (for ST only)"<<endl;
	out<<"3 3 16				:Number of peaks, width,starting coordinate (only for ST_6V_peaks)"<<endl;
	out<<"======================GEOMETRY READING======================"<<endl;
	out<<"1		:Geometry Reading format, 0=decimal,1=binary"<<endl;
        out<<i<<endl;
        	
	*/


//==============================================================================================	
//========================================SINGLE_PHASE==========================================
//==============================================================================================




	

/*   
        out<<"=================SINGLE PHASE LATTICE BOLTZMANN CODE INPUT FILE==============================="<<endl;
out<<"LV60_106_128_128.dat        :Geometry filename"<<endl;
out<<"106 128 128                 :nx ny nz"<<endl;
out<<"1000	                 :Maxmum time step"<<endl;
out<<"7.249	                :x=1 (um) Resolution (for Permeability calculation, resolution/dx)"<<endl;
out<<"0	        	       :Pressure/Velocity Boundary (1=Yes, 0=No)   Convergence accelerator (1=yes, 0=no)"<<endl;
out<<"2.0e-6 0.0 0.0          :body force for x,y,z"<<endl;
out<<"0 1.0 0 1.003 	:Pressure Boundary in X direction (Format detials can be found within this file)p=c_s^2rho"<<endl;
out<<"0 1.0 0 1.0	:Pressure Boundary in Y direction (Format detials can be found within this file)"<<endl;
out<<"0 1.0 0 1.0	:Pressure Boundary in Z direction (Format detials can be found within this file)"<<endl;
out<<"0 0.0 0 0.07	:Velocity Boundary in X direction (Format detials can be found within this file)"<<endl;
out<<"0 0.0 0 0.07	:Velocity Boundary in Y direction (Format detials can be found within this file)"<<endl;
out<<"0 0.0 0 0.07	:Velocity Boundary in Z direction (Format detials can be found within this file)"<<endl;
out<<"0.1            	:Viscosity"<<endl;
out<<"0.0e-3 0.0 0.0	:initial velocity for x,y,z"<<endl;
out<<"=========OUTPUT==CONTROL==================="<<endl;
out<<"1		:Permeability writing (1=yes, 0=no)"<<endl;
out<<"1               :Direction of Permeability Calculation (1=X, 2=Y, 3=Z)"<<endl;
out<<"50		:Frequency of results writing (interval in time steps)"<<endl;
out<<"0               :Memory Saving mode for output subroutains (0=No, 1=Yes, recommended for large outputs)"<<endl;
out<<"-1		:Frequency of velocity field writing (in vtk format, -1=no velocity writing)"<<endl;
out<<"-1		:Frequency of density field writing (in vtk format, -1=no density writing)"<<endl;
out<<"==============ADVANCE==PARAMETER=========="<<endl;
out<<"0 1.0 2.0	:Self define lattice velocity: 0=DEFAULT, dx, dt ((u_x+u_y+u_z)*dt/dx<=1 Courant Number)"<<endl;
out<<"./		:OUTPUT PATH,DEFAULT: ./ (REMEMBER TO INCLUDE / AT THE END OF PATH)"<<endl;
out<<"2              :PRESSURE AND VELOCITY BOUNDARY CONDITION OPTIONS: 0,1,2,3: EBC_S,EBC_D,TOLKE_BC,NEBC_D"<<endl;
out<<"1				:LOCAL PERM CALCULATION (0=NO, 1=YES)"<<endl;
out<<"0 0 0				:Pemeability calculation Partially "<<endl; 
out<<"130 20				:Permeability calculation partially Starting point and Ending point in X"<<endl;
out<<"0 100				:Permeability calculation partially Starting point and Ending point in Y"<<endl;
out<<"0 100				:Permeability calculation partially Starting point and Ending point in Z"<<endl;
out<<"======================BACKUP OPTIONS==========================="<<endl;
out<<"-1		:BACKUP FREQUENCY (-1=NO BACKUP)"<<endl;
out<<"0               :USE BACKUP DATA (0=OFF, 1=ON)"<<endl;
out<<"-1      :Velocity output for Solute transport every x timesteps (-1=NO)"<<endl;
out<<"======================GROP PERM PARAMETER ====================="<<endl;
out<<"0		:GROP PERM OPTION (0=OFF, or Frequency of exacution)"<<endl;
out<<"1 1 2 4		:NUMBERS OF SUBDOMAINS FOR SUB-PERM 1,2,3,4"<<endl;
out<<"512 256 128 64	:SIZE OF SUBDOMAINS FOR SUB-PERM 1,2,3,4"<<endl;
out<<"0 1 1		:STARTING COORDINATES OF SUB-PERM 1"<<endl;
out<<"128 129 129	:STARTING COORDINATES OF SUB-PERM 2"<<endl;
out<<"128 129 129	:STARTING COORDINATES OF SUB-PERM 3"<<endl;
out<<"128 129 129	:STARTING COORDINATES OF SUB-PERM 4"<<endl;
out<<"======================GEOMETRY READING======================="<<endl;
out<<"1		:Geometry Reading format, 0=decimal,1=binary"<<endl;
out<<i<<endl;
*/




//=========================================================================
//============================MC_CG========================================
//=========================================================================




out<<"========3D LATTICE BOLTZMANN MPI CODE--MULTI COMPONENT CG======="<<endl;
out<<"LV60_bufn_309_302_302.dat         :Geometry filenameEmily_Berea_y_614_764_7.dat"<<endl;
out<<"phase_309_302_302.dat           	:Initial components distribution"<<endl;
out<<"309 302 302               	:nx ny nz"<<endl;
out<<"10000000		     		:Maximum time step"<<endl;
out<<"7.249	             		:dx (um) Resolution (for Permeability calculation)"<<endl;
out<<"0	        	      	:Pressure Or Velocity Boundary (1=YES, 0=No)"<<endl;
out<<"1				:Psi constant BC"<<endl;
out<<"1.0e-4 0.0 0.0	          	:body force for x,y,z"<<endl;
out<<"0 1.0 0 1.0027	    		:Pressure Boundary in X direction ()p=c_s^2*/rho"<<endl;
out<<"0 1.0 0 1.0	    		:Pressure Boundary in Y direction ()"<<endl;
out<<"0 1.0 0 1.0	    	:Pressure Boundary in Z direction (Format detials can be found within this file)"<<endl;
out<<"0 0.0 0 0.0	   	:Velocity Boundary in X direction (Format detials can be found within this file)"<<endl;
out<<"0 0.0 0 0.07	   	:Velocity Boundary in Y direction (Format detials can be found within this file)"<<endl;
out<<"0 0.0 0 0.07	   	:Velocity Boundary in Z direction (Format detials can be found within this file)"<<endl;
out<<"2 1			:Psi constant BC in X 0=OFF, 1=Ini psi, 2=Neibourghing,3=Rand"<<endl;
out<<"0 0				:Psi constant BC in Y"<<endl;
out<<"0 0				:Psi constant BC in Z"<<endl;
out<<"0.05    		       		:Viscosity (Component A, psi=1)"<<endl;
out<<"0.05    		       		:Viscosity (Component B, psi=-1)"<<endl;
out<<"0.85     		      	:Contact Angle Cos(Theta) (Positive=1 wetting, Negative -1 wetting)"<<endl;
out<<"0.7e-2    		      	:Surface tension (Kappa)"<<endl;
out<<"0.0 0.0 0.0			:initial velocity for x,y,z"<<endl;
out<<"18500	        	       	:Permeability (Single Phase mD)"<<endl;
out<<"=========OUTPUT==CONTROL==================="<<endl;
out<<"1		         	:Permeability writing (1=yes, 0=no)"<<endl;
out<<"1                       	:Direction of Permeability Calculation (1=X, 2=Y, 3=Z)"<<endl;
out<<"50		         	:Frequency of results writing (interval in time steps)"<<endl;
out<<"0               		:Memory Saving mode for output(0=No, 1=Yes,for large outputs, x,z opposite)"<<endl;
out<<"-1				:Freqency of velocity field writing (in vtk format, -1=no velocity writing)"<<endl;
out<<"-1				:Freqency of density field writing (in vtk format, -1=no density writing)"<<endl;
out<<"-1				:Freqency of concentration writing (in vtk format, -1=no density writing)"<<endl;
out<<"==============ADVANCE==PARAMETER=========="<<endl;
out<<"0 1.0 4.0	              :Self define lattice velocity: 0=DEFAULT,dx, dt ((u_x+u_y+u_z)*dt/dx<=1 Courant Number)"<<endl;
out<<"/work/jy810/Sandpack/LV60/"<<i<<"_	:OUTPUT PATH,DEFAULT: ./ (INCLUDE / AT THE END)"<<endl;
out<<"2            :PRESSURE AND VELOCITY BOUNDARY CONDITION OPTIONS: 0,1,2,3: EBC_S,EBC_D,TOLKE_BC,NEBC_D"<<endl;
out<<"1  1000	:MULTI-COMPONENT STABALIZER: (a,b) a=0=OFF, a=1=ON, BODY FORCE APPLIED AFTER b steps"<<endl;
out<<"-1               :PRESET SATUATION, 0--1 initial component distri not needed, -1=OFF(for Component A,1)"<<endl;
out<<"1				:PRESET VALUE FOR BUFFET AREA, 0=NO,1=COMP A,-1=COMP B (valid when preset satuation)"<<endl;
out<<"0.8				:Relative permeability calcualtion 0..1 (psi>=value, cal the flux for Comp1)"<<endl;
out<<"1 0 0				:Pemeability calculation Partially  (1=ON, 0=OFF)"<<endl;
out<<"308 9				:Permeability calculation partially Starting point and Ending point in X"<<endl;
out<<"0 100				:Permeability calculation partially Starting point and Ending point in Y"<<endl;
out<<"0 100				:Permeability calculation partially Starting point and Ending point in Z"<<endl;
out<<"======================BACKUP==CONTROL===================="<<endl;
out<<"-1		                   :BACKUP FREQUENCY (-1=NO BACKUP,0=Backup at the end of computation)"<<endl;
out<<20000*(i+1)<<"                              :INITIALIZATION WITH BACKUP DATA (0=OFF, 1=ON)"<<endl;
out<<"======================GEOMETRY READING======================="<<endl;
out<<"1		:Geometry Reading format, 0=decimal,1=binary"<<endl;
out<<"======================PRESSURE OR BODAY FORCE SETTING (CHANGE WITH TIME)============="<<endl;
out<<"0 1			:0=OFF,1,2,3=ON,1=X,2=Y,3=Z;|| 1=PRESSURE,2=BODY FORCE"<<endl;
out<<"-1 0.0 0.0 0.0		:ACT AT TIMESTEP n, PRESSURE N, PRESSURE P, BODYFORCE(1)"<<endl;
out<<"-1 0.0 0.0 0.0		:ACT AT TIMESTEP n, PRESSURE N, PRESSURE P, BODYFORCE(2)"<<endl;
out<<"-1 0.0 0.0 0.0		:ACT AT TIMESTEP n, PRESSURE N, PRESSURE P, BODYFORCE(3)"<<endl;
out<<"-1 0.0 0.0 0.0		:ACT AT TIMESTEP n, PRESSURE N, PRESSURE P, BODYFORCE(4)"<<endl;
out<<"-1 0.0 0.0 0.0		:ACT AT TIMESTEP n, PRESSURE N, PRESSURE P, BODYFORCE(5)"<<endl;
out<<i<<endl;

}




//======================SCRIPT WRITTING========================================
for (int i=0;i<total_number;i++)
{
        
        name.str("");
	name<<poreFileName2<<i<<".sc";
	ofstream out;
	out.open(name.str().c_str());
	
out<<"#!/bin/sh "<<endl;
out<<"#PBS -N Sandp"<<i<<endl;
out<<"#PBS -l walltime=190:00:00 "<<endl;
out<<"#PBS -l select=4:ncpus=12:icib=true "<<endl;
//out<<"#PBS -l select=6:ncpus=12 "<<endl;




out<<"module load intel-suite mpi"<<endl;

out<<"cp /work/jy810/Sandpack/LV60/*.dat ."<<endl;
out<<"cp /work/jy810/Sandpack/LV60/*_"<<i*20000<<"*.bin_input ."<<endl;

out<<endl;

out<<"mpiexec /home/jy810/Multi_Component/MC_CG "<<poreFileName<<i<<".inputdat";
out<<" /work/jy810/Sandpack/LV60/";
out<<">/work/jy810/Sandpack/LV60/"<<i<<"_report.txt"<<endl;



}

}


