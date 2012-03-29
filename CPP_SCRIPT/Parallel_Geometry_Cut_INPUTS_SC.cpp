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

int nx=300;
int ny=300;
int nz=300;
         
int sym_x=0;

int add_buf_x_n=0;


int add_buf_x_p=0;

char poreFileName[128]="LV60.dat";
char poreFileNameVTK[128]="LV60_";                             //vtk headings
char poreFileNameOut[128]="LV60_";                       //output headings
//output VTK file,0=no, 1=yes
int division_num=2;
int VTK_OUT=1;
int input_sc=1;
int total_number=11;
char poreFileNameInput[128]="INPUT_CG_";
char poreFileNameSc[128]="LV60_";
//===========VTK AND OUT ARE ALL WRITTEN IN BINARY FORMAT===============
//===========================================================
ostringstream name;
ostringstream name2;


int dir=3;     
int add_buf_y_n=0;
int add_buf_z_n=0;
int Zoom=1; //1,2,3,4...
int nums;
int add_buf_y_p=0;
int add_buf_z_p=0;
int sym_y=0;
int sym_z=0;
int per_div;
per_div=ny/division_num;
if (per_div>nz/division_num)
        per_div=nz/division_num;


int pls[4][3]={{0,2,2},{2,0,2},{2,2,0},{0,0,0}};

int sum=0;



int nx1,ny1,nz1;
nx1=(nx+pls[dir][0])*(sym_x+1)+add_buf_x_n+add_buf_x_p;
ny1=(ny+pls[dir][1])*(sym_y+1)+add_buf_y_n+add_buf_y_p;
nz1=(nz+pls[dir][2])*(sym_z+1)+add_buf_z_n+add_buf_z_p;


bool*** Solid_Int;
bool*** Solid;
int*** Solid_Int2;



	FILE *ftest;
	ifstream fin;
	
	ftest = fopen(poreFileName, "r");

	if(ftest == NULL)
	{
		cout << "\n The pore geometry file (" << poreFileName <<
			") does not exist!!!!\n";
		cout << " Please check the file\n\n";

		exit(0);
	}
	fclose(ftest);

	fin.open(poreFileName);


double pore;
	int i, j, k,ci,cj,ck;
	
	Solid_Int = new bool**[(nx+pls[dir][0])*(sym_x+1)];	///*********
	Solid = new bool**[nx];
	for (i=0;i<(nx+pls[dir][0])*(sym_x+1);i++)				///*********
		{
		Solid_Int[i]=new bool*[(ny+pls[dir][1])*(sym_y+1)];		///*********
		for (j=0;j<(ny+pls[dir][1])*(sym_y+1);j++)			///*********
			{
			Solid_Int[i][j]= new bool[(nz+pls[dir][2])*(sym_z+1)];///*********
			for (k=0;k<(nz+pls[dir][2])*(sym_z+1);k++)			///*********
				Solid_Int[i][j][k]= 1;
			      
			}
		}
		
	for (i=0; i<nx;i++)
	{
	       Solid[i] = new bool*[ny];
	       for (j=0;j<ny;j++)
	       {
	               Solid[i][j] = new bool[nz];
	               for (k=0;k<nz;k++)
	                       Solid[i][j][k] = 0;
	       }
	}
		
	Solid_Int2 = new int**[nx1*Zoom];	///*********
	
	for (i=0;i<nx1*Zoom;i++)				///*********
		Solid_Int2[i]=new int*[ny1*Zoom];

	Solid_Int2[0][0]=new int[nx1*ny1*nz1*Zoom*Zoom*Zoom];

	
 	for (int i=1;i<ny1*Zoom;i++)
               Solid_Int2[0][i]=Solid_Int2[0][i-1]+nz1*Zoom;
       
       for (int i=1;i<nx1*Zoom;i++)
       {
               Solid_Int2[i][0]=Solid_Int2[i-1][0]+ny1*nz1*Zoom*Zoom;
               for (int j=1;j<ny1*Zoom;j++)
                       Solid_Int2[i][j]=Solid_Int2[i][j-1]+nz1*Zoom;
       }
	





	for(k=0 ; k<nz ; k++)				///*********
	for(j=0 ; j<ny ; j++)
	for(i=0 ; i<nx ; i++)				///*********


	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore;
			
			//if (pore == 0.0)	{Solid[ci-1][cj-1][ck-1] = 0;}
			if (pore == 0.0)	{Solid[i][j][k] = 0;}
			
			//if (pore == 1.0) 	{Solid[ci-1][cj-1][ck-1] = 1;sum++;}
			if (pore == 1.0) 	{Solid[i][j][k] = 1;sum++;}
			
			
			
			
		}
		
	fin.close();
		
	
	cout<<"Porosity = "<<1-(double(sum)/(nx*ny*nz))<<endl;	
	//cout<<sum<<endl;	
		
		
	// Reading pore geometry
	for(k=pls[dir][2]/2 ; k<nz+pls[dir][2]/2 ; k++)				///*********
	for(j=pls[dir][1]/2 ; j<ny+pls[dir][1]/2 ; j++)
	for(i=pls[dir][0]/2 ; i<nx+pls[dir][0]/2 ; i++)				///*********

	{
	        Solid_Int[i][j][k]=Solid[i-pls[dir][0]/2][j-pls[dir][1]/2][k-pls[dir][2]/2];
	        
	        
	}

	if (sym_x==1)
	for(k=0 ; k<nz+pls[dir][2] ; k++)				///*********
	for(j=0 ; j<ny+pls[dir][1] ; j++)
	for(i=0 ; i<nx+pls[dir][0] ; i++)				///*********
		Solid_Int[nx+pls[dir][0]+i][j][k]=Solid_Int[nx+pls[dir][0]-1-i][j][k];


if (sym_y==1)
	for(k=0 ; k<nz+pls[dir][2] ; k++)				///*********
	for(j=0 ; j<ny+pls[dir][1] ; j++)
	for(i=0 ; i<nx+pls[dir][0] ; i++)				///*********
		Solid_Int[i][ny+pls[dir][1]+j][k]=Solid_Int[i][ny+pls[dir][1]-1-j][k];


if (sym_z==1)
	for(k=0 ; k<nz+pls[dir][2] ; k++)				///*********
	for(j=0 ; j<ny+pls[dir][1] ; j++)
	for(i=0 ; i<nx+pls[dir][0] ; i++)				///*********
		Solid_Int[i][j][nz+pls[dir][2]+k]=Solid_Int[i][j][nz+pls[dir][2]-1-k];

	
	for(k=0 ; k<(nz+pls[dir][2])*(sym_z+1)+add_buf_z_n+add_buf_z_p; k++)						
		for(j=0 ; j<(ny+pls[dir][1])*(sym_y+1)+add_buf_y_n+add_buf_y_p; j++)					
			for(i=0 ; i<(nx+pls[dir][0])*(sym_x+1)+add_buf_x_n+add_buf_x_p; i++)				
			if ((i>=add_buf_x_n) and (i<(nx+pls[dir][0])*(sym_x+1)+add_buf_x_n) and (j>=add_buf_y_n) and (j<(ny+pls[dir][1])*(sym_y+1)+add_buf_y_n) and (k>=add_buf_z_n) and (k<(nz+pls[dir][2])*(sym_z+1)+add_buf_z_n))
			for (ci=0;ci<Zoom;ci++)
				for (cj=0;cj<Zoom;cj++)
				for (ck=0;ck<Zoom;ck++) 
					Solid_Int2[i*Zoom+ci][j*Zoom+cj][k*Zoom+ck] = Solid_Int[i-add_buf_x_n][j-add_buf_y_n][k-add_buf_z_n];
			else
			for (ci=0;ci<Zoom;ci++)
				for (cj=0;cj<Zoom;cj++)
				for (ck=0;ck<Zoom;ck++) 
					Solid_Int2[i*Zoom+ci][j*Zoom+cj][k*Zoom+ck] = 0;

				nums=0;

for (int outk=0;outk<division_num;outk++)
        for (int outj=0;outj<division_num;outj++)
{
        
	if (VTK_OUT==1)
	{
	cout<<"Start writing VTK file  No."<<nums<<endl;
	cout<<nx1*Zoom<<"         "<<ny1*Zoom<<"         "<<nz1*Zoom<<endl;

	name.str("");
	name<<poreFileNameVTK<<nums<<".vtk";
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
	ofstream out;
	out.open(name.str().c_str());
	
	
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<nx1*Zoom<<"         "<<per_div<<"         "<<per_div<<endl;       ///*********
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<nx1*per_div*per_div*Zoom*Zoom*Zoom<<endl;				///*********
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	
	
	for (k=per_div*outk;k<per_div*(outk+1);k++)
	{
	for (j=per_div*outj;j<per_div*(outj+1);j++)
	for (i=0;i<nx1*Zoom;i++)
		out<<Solid_Int2[i][j][k]<<" ";
	}
	

	//out.write((char *)(&Solid_Int2[0][0][0]), sizeof(int)*nx1*ny1*nz1*Zoom*Zoom*Zoom); 

	//out.close();
	out.close();
	cout<<"VTK file ouput COMPLETE No."<<nums<<endl;
	}


	//=================================================================
	cout<<"Start writing DAT file No."<<nums<<endl;
	cout<<nx1*Zoom<<"         "<<per_div<<"         "<<per_div<<endl;
	name2.str("");
	name2<<poreFileNameOut<<nums<<"_"<<nx1*Zoom<<"x"<<per_div<<"x"<<per_div<<".dat";
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
	ofstream out2;
	out2.open(name2.str().c_str());
	
	for (k=per_div*outk;k<per_div*(outk+1);k++)
	{
	for (j=per_div*outj;j<per_div*(outj+1);j++)
	for (i=0;i<nx1*Zoom;i++)
		out2<<Solid_Int2[i][j][k]<<" ";
	}
	
	out2.close();
	//out2.write((char *)(&Solid_Int2[0][0][0]), sizeof(int)*nx1*ny1*nz1*Zoom*Zoom*Zoom); 

	//out2.close();

	cout<<"DAT file ouput complete  No."<<nums<<endl;nums++;
	
}


if (input_sc==1)
{
for (int geoi=0;geoi<nums;geoi++)
for (int i=1;i<=total_number;i++)
{
        
        name.str("");
	//name<<poreFileNameOut<<geoi<<"_"<<nx1*Zoom<<"x"<<per_div<<"x"<<per_div<<".dat";
	name<<poreFileNameInput<<geoi<<"_"<<i<<".inputdat";
	ofstream out;
	out.open(name.str().c_str());
	


	//===========================================SP==SOLUTE_HEAT==============================================
	/*
	out<<"=============INPUT FILE FOR 3D LATTICE BOLTZMANN MPI CODE--Solute and heat transfer======="<<endl;
	out<<poreFileNameOut<<geoi<<"_"<<nx1*Zoom<<"x"<<per_div<<"x"<<per_div<<".dat"<<endl;
	out<<"phase.dat       	    	:Initial components distribution"<<endl;
        out<<nx1*Zoom<<" "<<per_div<<" "<<per_div<<"                :nx ny nz"<<endl;
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
        out<<"./scr"<<geoi<<"_"<<i<<"_		:OUTPUT PATH,DEFAULT"<<endl;
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
out<<poreFileNameOut<<geoi<<"_"<<nx1*Zoom<<"x"<<per_div<<"x"<<per_div<<".dat"<<endl;
out<<nx1*Zoom<<" "<<per_div<<" "<<per_div<<"                 :nx ny nz"<<endl;
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
out<<"./sp"<<geoi<<"_"<<i<<"_		:OUTPUT PATH,DEFAULT: ./ (REMEMBER TO INCLUDE / AT THE END OF PATH)"<<endl;
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
out<<poreFileNameOut<<geoi<<"_"<<nx1*Zoom<<"x"<<per_div<<"x"<<per_div<<".dat"<<endl;
out<<"phase_309_302_302.dat           	:Initial components distribution"<<endl;
out<<nx1*Zoom<<" "<<per_div<<" "<<per_div<<"               	:nx ny nz"<<endl;
out<<"10000000		     		:Maximum time step"<<endl;
out<<"7.249	             		:dx (um) Resolution (for Permeability calculation)"<<endl;
out<<"0	        	      	:Pressure Or Velocity Boundary (1=YES, 0=No)"<<endl;
out<<"0				:Psi constant BC"<<endl;
out<<"7.0e-6 0.0 0.0	          	:body force for x,y,z"<<endl;
out<<"0 1.0 0 1.0027	    		:Pressure Boundary in X direction ()p=c_s^2*/rho"<<endl;
out<<"0 1.0 0 1.0	    		:Pressure Boundary in Y direction ()"<<endl;
out<<"0 1.0 0 1.0	    	:Pressure Boundary in Z direction (Format detials can be found within this file)"<<endl;
out<<"0 0.0 0 0.0	   	:Velocity Boundary in X direction (Format detials can be found within this file)"<<endl;
out<<"0 0.0 0 0.07	   	:Velocity Boundary in Y direction (Format detials can be found within this file)"<<endl;
out<<"0 0.0 0 0.07	   	:Velocity Boundary in Z direction (Format detials can be found within this file)"<<endl;
out<<"0 0			:Psi constant BC in X 0=OFF, 1=Ini psi, 2=Neibourghing,3=Rand"<<endl;
out<<"0 0				:Psi constant BC in Y"<<endl;
out<<"0 0				:Psi constant BC in Z"<<endl;
out<<"0.03    		       		:Viscosity (Component A, psi=1)"<<endl;
out<<"0.03    		       		:Viscosity (Component B, psi=-1)"<<endl;
out<<"0.87     		      	:Contact Angle Cos(Theta) (Positive=1 wetting, Negative -1 wetting)"<<endl;
out<<"0.9e-2    		      	:Surface tension (Kappa)"<<endl;
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
out<<"/work/jy810/Sandpack/LV60/rel"<<geoi<<"_"<<i<<"_	:OUTPUT PATH,DEFAULT: ./ (INCLUDE / AT THE END)"<<endl;
out<<"2            :PRESSURE AND VELOCITY BOUNDARY CONDITION OPTIONS: 0,1,2,3: EBC_S,EBC_D,TOLKE_BC,NEBC_D"<<endl;
out<<"1  1000	:MULTI-COMPONENT STABALIZER: (a,b) a=0=OFF, a=1=ON, BODY FORCE APPLIED AFTER b steps"<<endl;
out<<(double)1.0/total_number*i<<"               :PRESET SATUATION, 0--1, distri not needed, -1=OFF(for Comp A,1)"<<endl;
out<<"0				:PRESET VALUE FOR BUFFET AREA, 0=NO,1=COMP A,-1=COMP B (valid when preset satuation)"<<endl;
out<<"0.8				:Relative permeability calcualtion 0..1 (psi>=value, cal the flux for Comp1)"<<endl;
out<<"1 0 0				:Pemeability calculation Partially  (1=ON, 0=OFF)"<<endl;
out<<"308 9				:Permeability calculation partially Starting point and Ending point in X"<<endl;
out<<"0 100				:Permeability calculation partially Starting point and Ending point in Y"<<endl;
out<<"0 100				:Permeability calculation partially Starting point and Ending point in Z"<<endl;
out<<"======================BACKUP==CONTROL===================="<<endl;
out<<"-1		                   :BACKUP FREQUENCY (-1=NO BACKUP,0=Backup at the end of computation)"<<endl;
out<<"0                              :INITIALIZATION WITH BACKUP DATA (0=OFF, 1=ON)"<<endl;
out<<"======================GEOMETRY READING======================="<<endl;
out<<"1		:Geometry Reading format, 0=decimal,1=binary"<<endl;
out<<"======================PRESSURE OR BODAY FORCE SETTING (CHANGE WITH TIME)============="<<endl;
out<<"0 1			:0=OFF,1,2,3=ON,1=X,2=Y,3=Z;|| 1=PRESSURE,2=BODY FORCE"<<endl;
out<<"-1 0.0 0.0 0.0		:ACT AT TIMESTEP n, PRESSURE N, PRESSURE P, BODYFORCE(1)"<<endl;
out<<"-1 0.0 0.0 0.0		:ACT AT TIMESTEP n, PRESSURE N, PRESSURE P, BODYFORCE(2)"<<endl;
out<<"-1 0.0 0.0 0.0		:ACT AT TIMESTEP n, PRESSURE N, PRESSURE P, BODYFORCE(3)"<<endl;
out<<"-1 0.0 0.0 0.0		:ACT AT TIMESTEP n, PRESSURE N, PRESSURE P, BODYFORCE(4)"<<endl;
out<<"-1 0.0 0.0 0.0		:ACT AT TIMESTEP n, PRESSURE N, PRESSURE P, BODYFORCE(5)"<<endl;
out<<"0			:BODY FORCE APPLIED ON 1:Phase1,0:BOTH,-1:Phase2"<<endl;
out<<"=====================CAPILLARY PRESSURE MEASURMENTS================"<<endl;
out<<"0 1 2000                     :0=OFF,1,2,3=ON,1=X,2=Y,3=Z;|| 1=PRESSURE,2=BODY FORCE; time steps of 1 changement"<<endl;
out<<"1.0 0.993 6 10           :PRESSURE N,P; Chan times, term condition: no. output intervals"<<endl;                   
out<<"1.0e-6                          :Error of Saturation stable condition"<<endl;
out<<i<<endl;

}




//======================SCRIPT WRITTING========================================
for (int geoi=0;geoi<nums;geoi++)
for (int i=0;i<total_number;i++)
{
        
        name.str("");
	name<<poreFileNameSc<<geoi<<"_"<<i<<".sc";
	ofstream out;
	out.open(name.str().c_str());
	
out<<"#!/bin/sh "<<endl;
out<<"#PBS -N rel_perm"<<geoi<<"_"<<i<<endl;
out<<"#PBS -l walltime=69:00:00 "<<endl;
out<<"#PBS -l select=2:ncpus=12:icib=true "<<endl;
//out<<"#PBS -l select=2:ncpus=12 "<<endl;




out<<"module load intel-suite mpi"<<endl;

out<<"cp /work/jy810/Sandpack/LV60/*.dat ."<<endl;
out<<"cp /work/jy810/Sandpack/LV60/*.inputdat ."<<endl;

//out<<"cp /work/jy810/Sandpack/LV60/*_"<<i*20000<<"*.bin_input ."<<endl;

out<<endl;

//out<<"mpiexec /home/jy810/Single_Phase/SP "<<poreFileNameInput<<geoi<<"_"<<i<<".inputdat";
out<<"mpiexec /home/jy810/Multi_Component/MC_CG "<<poreFileNameInput<<geoi<<"_"<<i<<".inputdat";
//out<<" /work/jy810/Sandpack/LV60/";
//out<<">/work/jy810/Sandpack/LV60/"<<i<<"_report.txt"<<endl;



}
}








}


