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


double reso=5.5;		//Resolution microns/pixel
double time=2;		//Time micro-second/lu
double den=1000.0;		//Density kg/m^3


double la_dx=1.0;		//lattice dx
double la_dt=1.0;		//lattice dt
double la_den=1.0;		//lattice density



double la_ift=0.5e-2;		//lattice surface tension
double la_vis=0.05;		//lattice viscosity
double la_pressure=1.0e-5;	//lattice pressure
double la_dx2=1.0;		//lattice length
double la_dt2=1.0;		//lattice time
double la_v=1e-5;		//lattice velocity
double la_bodyf=1e-5;		//lattice force



double x0,t0,m0;

double p_dx,p_dt,p_ift,p_vis,p_pressure,p_v,p_bodyf;


x0=reso*1e-6;
t0=time*1e-6;
m0=den/la_den*x0*x0*x0;

p_v=la_v*x0/t0;
p_vis=la_vis*x0*x0/t0;


cout<<endl;
cout<<"============================="<<endl;
cout<<"lattice velocity "<<la_v<<"	Pysical unite "<<p_v<<" m/s"<<endl;
cout<<endl;

cout<<"============================="<<endl;
cout<<"lattice viscosity "<<la_vis<<"	Pysical unite "<<p_vis*1e6<<" mm^2/s, or "<<p_vis<<"m^2/s"<<endl;
cout<<endl;

cout<<"============================="<<endl;
cout<<"lattice surface tension "<<la_ift<<"	Pysical unite "<<la_ift*m0/t0/t0*1e3<<" mN/m;dyn/cm"<<endl;
cout<<endl;


cout<<"============================="<<endl;
cout<<"lattice pressure "<<la_pressure<<"	Pysical unite "<<la_pressure*m0/t0/t0/x0<<" Pa or "<<la_pressure*m0/t0/t0/x0/(6.8948e3)<<" psi"<<endl;
cout<<"or "<<la_pressure*m0/t0/t0/x0*1e-5<<" bar"<<endl;
cout<<endl;

cout<<"============================="<<endl;
cout<<"lattice body force accerlaration  "<<la_bodyf<<"	Pysical unite "<<la_bodyf*x0/t0/t0<<" N/Kg"<<endl;
cout<<endl;

}
