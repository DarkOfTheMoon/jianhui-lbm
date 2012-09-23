#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>

//D2Q9 STANDARD LATTICE MRT LATTICE BOLTZMANN METHOD
//UNIFORM MESH, LATTICE VELOCITY 1


using namespace std;        
const int Q=9;          
const int NX=64;		
const int NY=64;		
const double U=0.0;   //lib velocity (top boundary)
const double gx=0.000003; // gravity coeficient in x direction
const double gy=-0.00; // gravity coeficient in y direction	
double alfa,beta;// parameter for moment equiliburium 
double tau=1.0;
double tau_R=1.0;


const double s_e=1.0;
const double s_eps=1.0;
const double s_q=1.0;
const double s_v=1.0;
const double s_m=1.0;

const double s2_e=1.0;
const double s2_eps=1.0;
const double s2_q=1.0;
const double s2_v=1.0;//0.55
const double s2_m=1.0;


bool Solid[NX+1][NY+1];
double M[9][9]=	{{1,1,1,1,1,1,1,1,1},
		{-4,-1,-1,-1,-1,2,2,2,2},
		{4,-2,-2,-2,-2,1,1,1,1},
		{0,1,0,-1,0,1,-1,-1,1},
		{0,-2,0,2,0,1,-1,-1,1},
		{0,0,1,0,-1,1,1,-1,-1},
		{0,0,-2,0,2,1,1,-1,-1},
		{0,1,-1,1,-1,0,0,0,0},
		{0,0,0,0,0,1,-1,1,-1}};

double MI[9][9]= {{0.1111111111,-0.1111111111,0.11111111111,0,0,0,0,0,0},
		{0.1111111111,-0.0277777778,-0.0555555556,0.1666666667,-0.1666666667,0,0,0.25,0},
		{0.1111111111,-0.0277777778,-0.0555555556,0,-0.0000,0.1666666667,-0.1666666667,-0.2500,0},
		{0.1111111111,-0.0277777778,-0.0555555556,-0.1666666667,0.1666666667,0 ,0,0.2500,0},
		{0.1111111111,-0.0277777778,-0.0555555556,0,0.0000,-0.1666666667,0.1666666667,-0.2500,0},
		{0.1111111111,0.0555555556,0.0277777778,0.1666666667,0.0833333333,0.1666666667,0.0833333333,0,0.25},
		{0.1111111111,0.0555555556,0.0277777778,-0.1666666667,-0.0833333333,0.1666666667,0.0833333333,0,-0.25},
		{0.1111111111,0.0555555556,0.0277777778,-0.1666666667,-0.0833333333,-0.1666666667,-0.0833333333,0,0.25},
		{0.1111111111,0.0555555556,0.0277777778,0.1666666667,0.0833333333,-0.1666666667,-0.0833333333,0,-0.25}};


int eve[25][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1},
		{2,0},{0,2},{-2,0},{0,-2},{2,2},{-2,2},{-2,-2},{2,-2},{2,1},
		{1,2},{-1,2},{-2,1},{-2,-1},{-1,-2},{1,-2},{2,-1}};





double m[9];
double meq[9];
double m_R[9];
double meq_R[9];



int e[Q][2]={{0,0},{1,0},{0,1},{-1,0},{0,-1},{1,1},{-1,1},{-1,-1},{1,-1}};
double w[Q]={4.0/9,1.0/9,1.0/9,1.0/9,1.0/9,1.0/36,1.0/36,1.0/36,1.0/36};
double rho[NX+1][NY+1],u[NX+1][NY+1][2],u0[NX+1][NY+1][2],f[NX+1][NY+1][Q],F[NX+1][NY+1][Q];
double rho_R[NX+1][NY+1],u_R[NX+1][NY+1][2],f_R[NX+1][NY+1][Q],F_R[NX+1][NY+1][Q];
double psi[NX+1][NY+1],psi_R[NX+1][NY+1];

double rho_mac[NX+1][NY+1],u_mac[NX+1][NY+1][2];




int n;
double uMax,c,Re,dx,dy,Lx,Ly,dt,rho0,P0,tau_f,niu,error,SFx,SFy;

double Fx=0;
double Fy=0;// Interparticle force
double rho_v=1.0;
double rho_l=1.0;
double forcex[NX+1][NY+1];
double forcey[NX+1][NY+1];

double forcex_R[NX+1][NY+1];
double forcey_R[NX+1][NY+1];

void init();
void cylinder_creation(int, int, double);
void feq(int,int); //EQUILIBRIUM EQUATION COMPUTING
void collision();   
bool interior_node(int,int) ;
void periodic_streaming();
void streaming();
void comput_macro_variables();
void standard_bounceback_boundary(int,int);
void boundary_Zou_He_pressure(bool, double, bool,double, bool, double, bool, double);
void boundary_Zou_He_velocity(bool, double, bool,double, bool, double, bool, double);
void evolution();
void output(int m);
void outputpoi(int m);
void Error();
void output_matrix(int);
void output_velocity(int);
void output_velocity_x(int);
void output_velocity_y(int);
void outputpoi2(int m);
void output_density(int);
void output_density2(int);
//double S[9]={0,s_e,s_eps,s_m,s_q,s_m,s_q,s_v,s_v};
//double S_R[9]={0,s2_e,s2_eps,s2_m,s2_q,s2_m,s2_q,s2_v,s2_v};
double S[9]={s_v,s_v,s_v,s_v,s_v,s_v,s_v,s_v,s_v};
double S_R[9]={s2_v,s2_v,s2_v,s2_v,s2_v,s2_v,s2_v,s2_v,s2_v};
double G=-90;  //G value, could be changed when needed
double Gads=-47;  //Gads, parameter for Solid-Liquid force
double G_MC=2.0;
double G_MC_R=2.0;



int main()
{	


	
	

	using namespace std;
	init();output_density(0);//cout<<"initial finish"<<endl;
	for(n=0;n<=100000;n++)
	{	
		evolution();  
		if(n%50==0)
		{       
			Error();
			cout<<"The"<<n<<"th computation result:"<<endl;
			cout<<"Forcex "<<forcex[NX/2][NY/2]<<" Forcey "<<forcey[NX/2][NY/2]<<endl;
			cout<<"Forcex_R "<<forcex_R[NX/2][NY/2]<<" Forcey_R "<<forcey_R[NX/2][NY/2]<<endl;
			cout<<"The u,v of point(NX/2,NY/2) is: "<<setprecision(6)
				<<u[NX/2][NY/2][0]<<", "<<u[NX/2][NY/2][1]<<endl;
			cout<<"The max relative error of uv is: "
				<<setiosflags(ios::scientific)<<error<<endl;
                        
			if(n>=1)
			{
				if(n%500==0) 	
						{output_density(n);output_density2(n);
						//outputpoi(n);outputpoi2(n);
						 output_velocity(n);//output_matrix(n);
						output_velocity_y(n);output_velocity_x(n);
						//outputpoi2(n);
						};
				//if(error<1.0e-6) {output(n);break;};
			}
		}
	}
	return 0;
}

void init()
{	
        //=================
        //double res[9][9];
        //=================




	double usqr,vsqr;
	dx=1.0;
	dy=1.0;
	Lx=dx*double(NX);
	Ly=dy*double(NY);
	dt=1;//sqrt(3);
	c=dx/dt;
	rho0=0.0;
 	uMax=0.1;
	Re=0.01;
	//forcex=gx;forcey=gy;
	//niu=U*Lx/Re;
	niu=uMax*200/Re;
	tau_f=3.0*niu+0.5;
        tau_f=1.0;
	std::cout<<"tau_f= "<<tau_f<<endl;
	double pr; //raduis of the obstacles
        alfa=1;
	beta=-3;

	srand(time(0));                        //RANDOM NUMBER GENERATION SEED
	//cylinder_creation(30,52,9);




	 //for(int i=0;i<=NX;i++)	
	//	for(int j=0;j<=NY;j++)
	//	Solid[i][j]=0;
	
	for(int i=0;i<=NX;i++)	
		for(int j=0;j<=NY;j++)
		{    
			
			Solid[i][j]=0;


		     //pr=2.5+double(rand()%100)/100*3;
                     //=======================================================================   
                      // if (((i%30==0) or(i==NX)) && (j%25==0) && (j>=5) && (j<=NY-30))
		     //	  cylinder_creation(i,j,pr);
                     if ((j==0) or (j==NY))
		     		Solid[i][j]=1;
            

		    //==========================================================================
			

                        
			u[i][j][0]=0.000;
			u[i][j][1]=0;
			u_R[i][j][0]=0.000;
			u_R[i][j][1]=0;

			rho[i][j]=0.0;rho_R[i][j]=0.0;
			//rho[i][j]=rho_v+double(rand()%100)/100;
			//if ((j==0) or (j==NY))
			//	u[i][j][1]=uMax;
			//rho[i][j]=rho0;  //CONSTANT DENSITY INITIALIZATION
			//rho[i][j]=(double(rand()%10)>5) ? (0):(1);
			//if (sqrt((i-32)*(i-32)+(j-32)*(j-32))<=15)
				//rho[i][j]=rho0+double(rand()%100)/100;
			
			if ((j>=17) and (j<=47)) 
				rho[i][j]=rho_l;//+double(rand()%100)/100;
			else
				//rho_R[i][j]=rho_v;//+double(rand()%100)/100;
			        rho[i][j]=rho_l;
			

			//***********************************************************************

			//forcex[i][j]=0.0001;
			//forcey[i][j]=0.0;


			
			rho_mac[i][j]=rho[i][j]+rho_R[i][j]; 
			if (rho[i][j]+rho_R[i][j]!=0) //if (n==0) cout<<rho_mac[i][j]<<"  "<<n<<endl;
			{
			//u_mac[i][j][0]=(rho[i][j]*u[i][j][0]*s_m+rho_R[i][j]*u_R[i][j][0]*s2_m)/(rho[i][j]*s_m+rho_R[i][j]*s2_m);
			//u_mac[i][j][1]=(rho[i][j]*u[i][j][1]*s_m+rho_R[i][j]*u_R[i][j][1]*s2_m)/(rho[i][j]*s_m+rho_R[i][j]*s2_m);
			u_mac[i][j][0]=u[i][j][0]+u_R[i][j][0];
			u_mac[i][j][1]=u[i][j][1]+u_R[i][j][1];
			
			}
			else

			{
			u_mac[i][j][0]=0.0;
			u_mac[i][j][1]=0.0;
			}


			//***********************************************************************


			//INITIALIZATION OF m and f

			
				for (int mi=0; mi<9; mi++)
					{
					f[i][j][mi]=w[mi]*rho[i][j];
					f_R[i][j][mi]=w[mi]*rho_R[i][j];
					};
				






		//cout<<' '<<i<<' '<<j<<endl;
                //cout<<meq[0]<<' '<<meq[1]<<' '<<meq[2]<<' '<<meq[3]<< ' '<<endl;
 		//cout<<meq[4]<<' '<<meq[5]<<' '<<meq[6]<<' '<<meq[7]<< ' '<<meq[8]<<endl;
		//cout<<endl;
                                
	 	}


//for (int ik=0;ik<9;ik++) 
//	{
//	for (int jk=0;jk<9;jk++)
//	{res[ik][jk]=0;
//		for (int kk=0;kk<9;kk++)
//		   res[ik][jk]=res[ik][jk]+M[ik][kk]*MI[kk][jk];
//         cout<<res[ik][jk]<<' ';
//     
//    }
//	cout<<endl;
//	}

		
}

void feq(int ik, int jk)//,int k,double rho, double u[2])
{
	
        
	int ip,jp;
	double tmp;
	//double G=-0.0120;  //G value, could be changed when needed
	//double Gads=-47;  //Gads, parameter for Solid-Liquid force
	


	double TE[4]={9,13,25,9};


//=========================LONG==RANGE==MC==FROCE==========================================================
	
	double weight[4][5]={{1.0/3.0,1.0/12.0,0,0,0},{4.0/15.0,1.0/10.0,1.0/120.0,0,0},
				 {4.0/21.0,4.0/45.0,1.0/60.0,2.0/315.0,1.0/5040.0},
				{3.0,0.75,0.0,0.0,0.0}};
	double Fx_R,Fy_R,SFx_R,SFy_R;
	double tmpgx,tmpgy,tmp_Rgx,tmp_Rgy;

	int mod=2;
	
		Fx=0;SFx=0;
		Fy=0;SFy=0;
		Fx_R=0;SFx_R=0;
		Fy_R=0;SFy_R=0;

		for(int kl=1;kl<TE[mod];kl++)
			{       

				//=====================================================================

				ip=ik+eve[kl][0];if (ip<0) {ip=NX+1+ip;}; if (ip>NX) {ip=ip-NX-1;};
				jp=jk+eve[kl][1];if (jp<=0) {jp=1-jp;}; if (jp>=NY) {jp=NY-1-(jp-NY);};
				//jp=jk+eve[kl][1];if (jp<0) {jp=NY+1+jp;}; if (jp>NY) {jp=jp-NY-1;};
				// WARNING: THIS PART IS FOR PERIODIC BOUNDARY CONDITIONS
				// FOR OTHER BOUNDARIERS, IT SHOULD BE MODIFIED	
				//=====================================================================
				if (Solid[ip][jp]==0)
				{
				switch ((eve[kl][0]*eve[kl][0]+eve[kl][1]*eve[kl][1]))
					{
					case 1:
					
						Fx+=(weight[mod][0])*eve[kl][0]*rho_R[ip][jp];
						Fy+=(weight[mod][0])*eve[kl][1]*rho_R[ip][jp];
						Fx_R+=(weight[mod][0])*eve[kl][0]*rho[ip][jp];
						Fy_R+=(weight[mod][0])*eve[kl][1]*rho[ip][jp];

						break;

						
					case 2:
					
						Fx+=(weight[mod][1])*eve[kl][0]*rho_R[ip][jp];
						Fy+=(weight[mod][1])*eve[kl][1]*rho_R[ip][jp];
						Fx_R+=(weight[mod][1])*eve[kl][0]*rho[ip][jp];
						Fy_R+=(weight[mod][1])*eve[kl][1]*rho[ip][jp];
						break;
					case 4:
						Fx+=(weight[mod][2])*eve[kl][0]*rho_R[ip][jp];
						Fy+=(weight[mod][2])*eve[kl][1]*rho_R[ip][jp];
						Fx_R+=(weight[mod][2])*eve[kl][0]*rho[ip][jp];
						Fy_R+=(weight[mod][2])*eve[kl][1]*rho[ip][jp];

						break;

					case 5:
						Fx+=(weight[mod][3])*eve[kl][0]*rho_R[ip][jp];
						Fy+=(weight[mod][3])*eve[kl][1]*rho_R[ip][jp];
						Fx_R+=(weight[mod][3])*eve[kl][0]*rho[ip][jp];
						Fy_R+=(weight[mod][3])*eve[kl][1]*rho[ip][jp];

						break;

					case 8:
						Fx+=(weight[mod][4])*eve[kl][0]*rho_R[ip][jp];
						Fy+=(weight[mod][4])*eve[kl][1]*rho_R[ip][jp];
						Fx_R+=(weight[mod][4])*eve[kl][0]*rho[ip][jp];
						Fy_R+=(weight[mod][4])*eve[kl][1]*rho[ip][jp];

						break;

				

					}
			}
			}


		
		tmp=-G_MC*rho[ik][jk]*Fx;Fx=tmp/3;
		tmp=-G_MC*rho[ik][jk]*Fy;Fy=tmp/3;

		tmp=-G_MC_R*rho_R[ik][jk]*Fx_R;Fx_R=tmp/3;
		tmp=-G_MC_R*rho_R[ik][jk]*Fy_R;Fy_R=tmp/3;


		SFx=-Gads*rho[ik][jk]*SFx;
		SFy=-Gads*rho[ik][jk]*SFy;
		
		//if (Fx!=0)
		//cout<<Fx<<" "<<ik<<","<<jk<<endl;



	//interparicle force included.
 	forcex[ik][jk]=gx*rho[ik][jk]+Fx+SFx;forcey[ik][jk]=gy*rho[ik][jk]+Fy+SFy;
	forcex_R[ik][jk]=gx*rho_R[ik][jk]+Fx_R+SFx_R;forcey_R[ik][jk]=gy*rho_R[ik][jk]+Fy_R+SFy_R;
	/*
	if (rho[ik][jk]<rho_l*0.5) {tmpgx=gx*rho[ik][jk];tmpgy=gy*rho[ik][jk];} else {tmpgx=gx*rho_l;tmpgy=gy*rho_l;};
	if (rho_R[ik][jk]<rho_v*0.5) {tmp_Rgx=gx*rho_R[ik][jk];tmp_Rgy=gy*rho_R[ik][jk];} else {tmp_Rgx=gx*rho_v;tmp_Rgy=gy*rho_v;};

	forcex[ik][jk]=tmpgx+Fx+SFx;forcey[ik][jk]=tmpgy+Fy+SFy;
	forcex_R[ik][jk]=tmp_Rgx+Fx_R+SFx_R;forcey_R[ik][jk]=tmp_Rgy+Fy_R+SFy_R;
	*/
}





void periodic_streaming()
{


	


	int ip,jp;
	for(int i=0;i<=NX;i++)	
		for(int j=0;j<=NY;j++)
			if (!interior_node(i,j))
                        for(int k=0;k<Q;k++)
			{       
				ip=i+e[k][0];if (ip<0) {ip=NX;}; if (ip>NX) {ip=0;};
				jp=j+e[k][1];if (jp<0) {jp=NY;}; if (jp>NY) {jp=0;};
				//if ((jp>=0) && (jp<=NY))// && (ip>=0) && (ip<=NX))
					F[ip][jp][k]=f[i][j][k];
					F_R[ip][jp][k]=f_R[i][j][k];
					
			}
			//else
			//	for(int k=0;k<Q;k++)
			//		{
			//		ip=i+e[k][0];jp=j+e[k][1];
			//		if ((jp>=0) && (jp<=NY) &&(ip>=0) &&(ip<=NX))
			//		F[ip][jp][k]=f[i][j][k];
			//		};



	for(int i=0;i<=NX;i++)	
		for(int j=0;j<=NY;j++)
                		{                 
	                        f[i][j][0]=F[i][j][0];f[i][j][1]=F[i][j][1];f[i][j][2]=F[i][j][2];
				f[i][j][3]=F[i][j][3];f[i][j][4]=F[i][j][4];f[i][j][5]=F[i][j][5];
				f[i][j][6]=F[i][j][6];f[i][j][7]=F[i][j][7];f[i][j][8]=F[i][j][8];

				f_R[i][j][0]=F_R[i][j][0];f_R[i][j][1]=F_R[i][j][1];f_R[i][j][2]=F_R[i][j][2];
				f_R[i][j][3]=F_R[i][j][3];f_R[i][j][4]=F_R[i][j][4];f_R[i][j][5]=F_R[i][j][5];
				f_R[i][j][6]=F_R[i][j][6];f_R[i][j][7]=F_R[i][j][7];f_R[i][j][8]=F_R[i][j][8];


				};


}

void streaming()
{


	int ip,jp;
	for(int i=0;i<=NX;i++)	
		for(int j=0;j<=NY;j++)
 			if (!Solid[i][j])
                        for(int k=0;k<Q;k++)
			{       
				ip=i-e[k][0];
				jp=j-e[k][1];
				if ((ip>=0) && (ip<=NX) && (jp>=0) && (jp<=NY))
				F[i][j][k]=f[ip][jp][k];
				F_R[i][j][k]=f_R[ip][jp][k];
			};

	for(int i=0;i<=NX;i++)	
		for(int j=0;j<=NY;j++)
                		{                 
	                        f[i][j][0]=F[i][j][0];f[i][j][1]=F[i][j][1];f[i][j][2]=F[i][j][2];
				f[i][j][3]=F[i][j][3];f[i][j][4]=F[i][j][4];f[i][j][5]=F[i][j][5];
				f[i][j][6]=F[i][j][6];f[i][j][7]=F[i][j][7];f[i][j][8]=F[i][j][8];

				f_R[i][j][0]=F_R[i][j][0];f_R[i][j][1]=F_R[i][j][1];f_R[i][j][2]=F_R[i][j][2];
				f_R[i][j][3]=F_R[i][j][3];f_R[i][j][4]=F_R[i][j][4];f_R[i][j][5]=F_R[i][j][5];
				f_R[i][j][6]=F_R[i][j][6];f_R[i][j][7]=F_R[i][j][7];f_R[i][j][8]=F_R[i][j][8];
				};


}


bool interior_node(int ti,int tj)
{

	
	bool status_solid_interior;

	status_solid_interior=1;
	if (Solid[ti][tj]==0) 
		status_solid_interior=0; 
	else
    		for (int k=0;k<Q;k++)
			if ((ti+e[k][0]>=0) && (ti+e[k][0]<=NX) && (tj+e[k][1]>=0) && (tj+e[k][1]<=NY) && (Solid[ti+e[k][0]][tj+e[k][1]]==0))
				status_solid_interior=0;

return status_solid_interior; 



}




void collision()
{


double usqr,vsqr;
double F_hat[9],GuoF[9];
double F_hat_R[9],GuoF_R[9];
double lm0,lm1,eu,uv,uv2,eu2,feq_r[9],feq_R[9],delta_feq_r[9],delta_feq_R[9];
double uprim,vprim;



	for(int i=0;i<=NX;i++)	
		for(int j=0;j<=NY;j++)
			//if (!Solid[i][j]) //&& (i!=0) && (i!=NX))  // NOT INTERIOR SOLID NODE
								   // AND DO NOT INCLUDE THE FLUX BOUNDARY
			if (!Solid[i][j])     // LIQUID PHASE SEPERATION
                        //WARNNING:THE SOUND SPEED HERE IS sqrt(RT)=1/sqrt(3)

			{
			feq(i,j);
			//=================FORCE TERM_GUO=========================================
			for (int k=0;k<9;k++)
			{	
			lm0=((e[k][0]-u[i][j][0])*forcex[i][j]+(e[k][1]-u[i][j][1])*forcey[i][j])*3.0;
			lm1=(e[k][0]*u[i][j][0]+e[k][1]*u[i][j][1])*(e[k][0]*forcex[i][j]+e[k][1]*forcey[i][j])*9.0;
			GuoF[k]=w[k]*(lm0+lm1);

			lm0=((e[k][0]-u_R[i][j][0])*forcex_R[i][j]+(e[k][1]-u_R[i][j][1])*forcey_R[i][j])*3.0;
			lm1=(e[k][0]*u_R[i][j][0]+e[k][1]*u_R[i][j][1])*(e[k][0]*forcex_R[i][j]+e[k][1]*forcey_R[i][j])*9.0;
			GuoF_R[k]=w[k]*(lm0+lm1);

				
			//======================HSD FORCE==========================================
			
			//lm0=((e[k][0]-u[i][j][0])*forcex[i][j]+(e[k][1]-u[i][j][1])*forcey[i][j])*3.0/rho[i][j];
			//eu=e[k][0]*u[i][j][0]+e[k][1]*u[i][j][1];
			//uv=u[i][j][0]*u[i][j][0]+u[i][j][1]*u[i][j][1];
        		
			//feq_l=w[k]*rho[i][j]*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
			//GuoF[k]=lm0*feq_l;
			//=========================================================================


			//*****************************TEST**PART***********************************
			//SDF=lm0*feq_l/rho[i][j];
			//feq_l=feq_l*(1-lm0/2);
			//SDF=(1-1.0/(2.0*tau_f))*lm0*feq_l;	

			//f[i][j][k]=f[i][j][k]-GuoF[k]/2;
			//***************************************************************************

			//GuoF[k]=0.0;
			//GuoF_R[k]=0.0;


				
			}

			//*****************************TEST**PART***********************************
			
			
		
			//***************************************************************************

			//=====================equilibrium of moment=================================
			//uprim=(rho[i][j]*u[i][j][0]/s_m+rho_R[i][j]*u_R[i][j][0]/s2_m)/(rho[i][j]/s_m+rho_R[i][j]/s2_m);
			//vprim=(rho[i][j]*u[i][j][1]/s_m+rho_R[i][j]*u_R[i][j][1]/s2_m)/(rho[i][j]/s_m+rho_R[i][j]/s2_m);
			
			
			//uprim=u_mac[i][j][0];
			//vprim=u_mac[i][j][1];
			//usqr=uprim*uprim;
			//vsqr=vprim*vprim;


//---------------------------------------------------------------------------------------------
			if (rho[i][j]==0)
			{
			uprim=0.0;vprim=0.0;
			}
			else
			{
			uprim=u[i][j][0];
			vprim=u[i][j][1];
			}
			


			usqr=uprim*uprim;
			vsqr=vprim*vprim;
			uv=usqr+vsqr;
			uv2=(uprim+forcex[i][j]*rho[i][j])*(uprim+forcex[i][j]*rho[i][j])+(vprim+forcey[i][j]*rho[i][j])*(vprim+forcey[i][j]*rho[i][j]);
			
//-----------------------------------------------------------------------------------------------
			for (int ks=0; ks<9;ks++)
				{
				eu=e[ks][0]*uprim+e[ks][1]*vprim;
				eu2=e[ks][0]*(uprim+forcex[i][j]*rho[i][j])+e[ks][1]*(vprim+forcey[i][j]*rho[i][j]);
				feq_r[ks]=w[ks]*rho[i][j]*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
				delta_feq_r[ks]=w[ks]*rho[i][j]*(1.0+3.0*eu2+4.5*eu2*eu2-1.5*uv2)-feq_r[ks];
				}


                        //meq[0]=1;meq[1]=-2+3*(usqr+vsqr);meq[2]=alfa+beta*(usqr+vsqr);
			//meq[3]=uprim;meq[4]=-uprim;meq[5]=vprim;
			//meq[6]=-vprim;meq[7]=usqr-vsqr;meq[8]=uprim*vprim;

			
			//uprim=u_mac[i][j][0];
			//vprim=u_mac[i][j][1];
			//usqr=uprim*uprim;
			//vsqr=vprim*vprim;
//-------------------------------------------------------------------------------------------------------
			if (rho_R[i][j]==0)
			{
			uprim=0.0;vprim=0.0;
			}
			else
			{
			uprim=u_R[i][j][0];
			vprim=u_R[i][j][1];
			}

			usqr=uprim*uprim;
			vsqr=vprim*vprim;
			uv=usqr+vsqr;
			uv2=(uprim+forcex_R[i][j]*rho_R[i][j])*(uprim+forcex_R[i][j]*rho_R[i][j])+(vprim+forcey_R[i][j]*rho_R[i][j])*(vprim+forcey_R[i][j]*rho_R[i][j]);
//-------------------------------------------------------------------------------------------------------
		for (int ks=0; ks<9;ks++)
				{
				eu=e[ks][0]*uprim+e[ks][1]*vprim;
				eu2=e[ks][0]*(uprim+forcex_R[i][j]*rho_R[i][j])+e[ks][1]*(vprim+forcey_R[i][j]*rho_R[i][j]);
				feq_R[ks]=w[ks]*rho_R[i][j]*(1.0+3.0*eu+4.5*eu*eu-1.5*uv);
				delta_feq_R[ks]=w[ks]*rho_R[i][j]*(1.0+3.0*eu2+4.5*eu2*eu2-1.5*uv2)-feq_R[ks];
				}




			for (int sk=0;sk<9;sk++)
			{
			        f[i][j][sk]=f[i][j][sk]+(feq_r[sk]-f[i][j][sk])/tau+delta_feq_r[sk];
			        f_R[i][j][sk]=f_R[i][j][sk]+(feq_R[sk]-f_R[i][j][sk])/tau_R+delta_feq_R[sk];
			        //f[i][j][sk]=f[i][j][sk]+(feq_r[sk]-f[i][j][sk])/tau;
			        //f_R[i][j][sk]=f_R[i][j][sk]+(feq_R[sk]-f_R[i][j][sk])/tau_R;
			}
			
			}
                        else    
				
				if (!interior_node(i,j)) 
                        		{standard_bounceback_boundary(i,j);};    //BOUNCEBACK BOUNDARY CONDITION
			

		
}








void comput_macro_variables()
{
	
	double psi0=4.0;
	double eos_a=1;
	double eos_b=4;
	double eos_R=1;
	
	double Pre,bp4,bp4s,bp4c,alf;
	double T=0.07556;
	double omega=0.344;



	for(int i=0;i<=NX;i++)	
		for(int j=0;j<=NY;j++)
                   
			{ 



			if (!Solid[i][j]) // NOT INTERIOR SOLID NODE
			{
				u0[i][j][0]=u_mac[i][j][0];
				u0[i][j][1]=u_mac[i][j][1];
				rho[i][j]=0;
				u[i][j][0]=0;
				u[i][j][1]=0;
				for(int k=0;k<Q;k++)
					{
					//f[i][j][k]=F[i][j][k];
					rho[i][j]+=f[i][j][k];
					u[i][j][0]+=e[k][0]*f[i][j][k];
					u[i][j][1]+=e[k][1]*f[i][j][k];
					}
				
				if (rho[i][j]==0)
				{ 
				u[i][j][0]=0;
				u[i][j][1]=0;
				}
				else
				{
				u[i][j][0]=(u[i][j][0]+forcex[i][j]/2)/rho[i][j];
				u[i][j][1]=(u[i][j][1]+forcey[i][j]/2)/rho[i][j];
				
				//u[i][j][0]/=rho[i][j];
				//u[i][j][1]/=rho[i][j];

				
				}

				
				rho_R[i][j]=0;
				u_R[i][j][0]=0;
				u_R[i][j][1]=0;
				for(int k=0;k<Q;k++)
					{
					//f[i][j][k]=F[i][j][k];
					rho_R[i][j]+=f_R[i][j][k];
					u_R[i][j][0]+=e[k][0]*f_R[i][j][k];
					u_R[i][j][1]+=e[k][1]*f_R[i][j][k];
					}
				

				if (rho_R[i][j]==0)
				{ 
				u_R[i][j][0]=0;
				u_R[i][j][1]=0;
				}
				else
				{
				u_R[i][j][0]=(u_R[i][j][0]+forcex_R[i][j]/2)/rho_R[i][j];
				u_R[i][j][1]=(u_R[i][j][1]+forcey_R[i][j]/2)/rho_R[i][j];

				//u_R[i][j][0]/=rho_R[i][j];
				//u_R[i][j][1]/=rho_R[i][j];

				}
				

		rho_mac[i][j]=rho[i][j]+rho_R[i][j]; 
		if (rho[i][j]+rho_R[i][j]!=0) //if (n==0) cout<<rho_mac[i][j]<<"  "<<n<<endl;
		{
		//if (n==0) cout<<rho_mac[i][j]<<"  "<<n<<endl;
		//u_mac[i][j][0]=(rho[i][j]*u[i][j][0]*s_m+rho_R[i][j]*u_R[i][j][0]*s2_m)/(rho[i][j]*s_m+rho_R[i][j]*s2_m);
		//u_mac[i][j][1]=(rho[i][j]*u[i][j][1]*s_m+rho_R[i][j]*u_R[i][j][1]*s2_m)/(rho[i][j]*s_m+rho_R[i][j]*s2_m);

		u_mac[i][j][0]=u[i][j][0]+u_R[i][j][0];
		u_mac[i][j][1]=u[i][j][1]+u_R[i][j][1];

		//u_mac[i][j][0]=(rho[i][j]*u[i][j][0]*s_m+rho_R[i][j]*u_R[i][j][0]*s2_m)/(rho[i][j]+rho_R[i][j]);
		//u_mac[i][j][1]=(rho[i][j]*u[i][j][1]*s_m+rho_R[i][j]*u_R[i][j][1]*s2_m)/(rho[i][j]+rho_R[i][j]);



		}
		else
		
			{
			u_mac[i][j][0]=0.0;
			u_mac[i][j][1]=0.0;
			}









//=======================================C-S EOS========================================
//		bp4=eos_b*rho[i][j]/4.0;
//		Pre=rho[i][j]*eos_R*T*(1.0+bp4+bp4*bp4-bp4*bp4*bp4)/((1.0-bp4)*(1.0-bp4)*(1.0-bp4))-eos_a*rho[i][j]*rho[i][j];
//		psi[i][j]=Pre-rho[i][j]/3;
//======================================================================================



//==========================SHAN CHEN EOS===============================================
			//psi[i][j]=exp(1/rho[i][j]);
//======================================================================================



//============================================Enskog===================================
//		bp4=eos_b*rho[i][j]/4.0;
//		Pre=rho[i][j]*eos_R*T*(1.0+bp4+bp4*bp4-bp4*bp4*bp4)/((1.0-bp4)*(1.0-bp4)*(1.0-bp4))-eos_a*rho[i][j]*rho[i][j];
//		psi[i][j]=Pre-rho[i][j]/3;

//======================================================================================

			
//==========================VDW EOS from Sukop book=====================================
			psi[i][j]=psi0*exp(-rho0/rho[i][j]);
//======================================================================================

//*******************************************************
//	Pre=(psi0*exp(-rho0/rho[i][j]))*(psi0*exp(-rho0/rho[i][j]))*G/6+rho[i][j]/3;
//	psi[i][j]=Pre-rho[i][j]/3;
//	bp4=(1.0-eos_b*rho[i][j])/((1.0-2*eos_b*rho[i][j])*(1.0-2*eos_b*rho[i][j])*(1.0-2*eos_b*rho[i][j]));
//	bp4=1.0+5.0/8.0*eos_b*rho[i][j]+0.2869*eos_b*rho[i][j]*eos_b*rho[i][j]+0.1103*eos_b*rho[i][j]*eos_b*rho[i][j]*eos_b*rho[i][j];
//	psi[i][j]=eos_b*rho[i][j]*rho[i][j]*bp4*eos_R*T-eos_a*rho[i][j]*rho[i][j];
	
//*******************************************************


				//u[i][j][0]/=rho[i][j];
				//u[i][j][1]/=rho[i][j];


				
					
		//====================HSD FORCE TERM=======********************=================
			//	u[i][j][0]=(u[i][j][0]+forcex[i][j]*s_v/2)/rho[i][j];
			//	u[i][j][1]=(u[i][j][1]+forcey[i][j]*s_v/2)/rho[i][j];
		//==============================================================================


		/*==========================GUO FORCE TERM======================================
				if (rho[i][j]==0)
				{
				u[i][j][0]=0;
				u[i][j][1]=0;
				}
				else
				{
				u[i][j][0]=u[i][j][0]+(forcex[i][j]/2)/rho[i][j];
				u[i][j][1]=u[i][j][1]+(forcey[i][j]/2)/rho[i][j];
				}

				if (rho_R[i][j]==0)
				{
				u_R[i][j][0]=0;
				u_R[i][j][1]=0;
				}
				else
				{
				u_R[i][j][0]=u_R[i][j][0]+(forcex_R[i][j]/2)/rho_R[i][j];
				u_R[i][j][1]=u_R[i][j][1]+(forcey_R[i][j]/2)/rho_R[i][j];
				}

				
		//==============================================================================*/


			}
			else
			{		u[i][j][0]=0; 
					u[i][j][1]=0;

					u_R[i][j][0]=0; 
					u_R[i][j][1]=0;
			}

                        
			}


//cout<<forcex[36][68]<<" "<<n<<endl;
}


void standard_bounceback_boundary(int it, int jt)
{

	double tmp;
			tmp = f[it][jt][1];f[it][jt][1] = f[it][jt][3];f[it][jt][3] = tmp;
			tmp = f[it][jt][2];f[it][jt][2] = f[it][jt][4];f[it][jt][4] = tmp;
			tmp = f[it][jt][5];f[it][jt][5] = f[it][jt][7];f[it][jt][7] = tmp;
			tmp = f[it][jt][6];f[it][jt][6] = f[it][jt][8];f[it][jt][8] = tmp;

			tmp = f_R[it][jt][1];f_R[it][jt][1] = f_R[it][jt][3];f_R[it][jt][3] = tmp;
			tmp = f_R[it][jt][2];f_R[it][jt][2] = f_R[it][jt][4];f_R[it][jt][4] = tmp;
			tmp = f_R[it][jt][5];f_R[it][jt][5] = f_R[it][jt][7];f_R[it][jt][7] = tmp;
			tmp = f_R[it][jt][6];f_R[it][jt][6] = f_R[it][jt][8];f_R[it][jt][8] = tmp;



                        	
			

}



void cylinder_creation(int cx, int cy, double cr)
{

// CYLINDER CREATION FOR INITIALIZATION
// FORMAT: 
// COORDINATE X, COORDINATE Y, RADUIS OF CYLINDER
// ALL VALUES ARE IN LATTICE DIMENSION
int Tmp_r;
	Tmp_r=int(cr)+1;


//================================UP AND BOTTOM BOUNDARY================
for(int i=0;i<=NX;i++)
	{
	Solid[i][NY]=1;
	Solid[i][0]=1;
	}
//======================================================================



for (int i=cx-Tmp_r;i<=cx+Tmp_r;i++)

	for (int j=cy-Tmp_r;j<=cy+Tmp_r;j++)
		

		  
			if ((sqrt((i-cx)*(i-cx)+(j-cy)*(j-cy))<=cr) && (i>=0) && (i<=NX) && (j>=0) && (j<=NY))
				Solid[i][j]=1;

				



}

void boundary_Zou_He_velocity(bool north, double v_n, bool south,double v_s, bool east, double v_e, bool west, double v_w)

// PLEASE SPECIFY THE BOUNDARY LOCATION BEFORE USE
// FOR EXAMPLE: IF SOUTH BOUNDARY IS INCLUDE, THEN THE BOOL VARIABLE SOUTH=1, AND ALSO SPECIFY THE LOCAL RHO
// IF THIS BOUNDARY CONDITION IS NOT GOING TO BE USED, PLEASE SET THE CORRESPONDING BOOL VARIABLE AS 0
// FORMAT:
// (NORTH BOUNDARY MARK,LOCAL RHO VALUE,SOUTH BOUNDARY MARK, LOCAL RHO VALUE,EAST BOUNDARY MARK, LOCAL RHO VALUE
// WEST BOUNDARY MARK, LOCAL RHO VALUE)

{

double rho0_local,ru,vtmp;


	if (north)
      
	for (int i=0;i<=NX;i++)
	{
		rho0_local = ( f[i][NY][0]+f[i][NY][1]+f[i][NY][3]+2.0*(f[i][NY][2]+f[i][NY][5]+f[i][NY][6]))/(1+v_n);
		ru = rho0_local*v_n;
		f[i][NY][4]=f[i][NY][2]-(2.0/3.0)*ru;
		f[i][NY][7]=f[i][NY][5]-(1.0/6.0)*ru+(1.0/2.0)*(f[i][NY][1]-f[i][NY][3]);
		f[i][NY][8]=f[i][NY][6]-(1.0/6.0)*ru+(1.0/2.0)*(f[i][NY][3]-f[i][NY][1]);

	}


if (south)
      for (int i=0;i<=NX;i++)
	{
		rho0_local = ( f[i][0][0]+f[i][0][1]+f[i][0][3]+2.0*(f[i][0][4]+f[i][0][7]+f[i][0][8]))/(1-v_s);
		ru = rho0_local*v_s;
		f[i][0][2]=f[i][0][4]+(2.0/3.0)*ru;
		f[i][0][5]=f[i][0][7]+(1.0/6.0)*ru-(1.0/2.0)*(f[i][0][1]-f[i][0][3]);
		f[i][0][6]=f[i][0][8]+(1.0/6.0)*ru-(1.0/2.0)*(f[i][0][3]-f[i][0][1]);

	}


if (east)
      for (int j=1;j<NY;j++)
	{
		rho0_local = ( f[NX][j][0]+f[NX][j][2]+f[NX][j][4]+2.0*(f[NX][j][1]+f[NX][j][5]+f[NX][j][8]))/(1+v_e);
		ru = rho0_local*v_e;
		f[NX][j][3]=f[NX][j][1]-(2.0/3.0)*ru;
		f[NX][j][7]=f[NX][j][5]-(1.0/6.0)*ru+(1.0/2.0)*(f[NX][j][2]-f[NX][j][4]);
		f[NX][j][6]=f[NX][j][8]-(1.0/6.0)*ru+(1.0/2.0)*(f[NX][j][4]-f[NX][j][2]);

	}


if (west)
      for (int j=1;j<NY;j++)
	{	//vtmp=4*uMax/((NY-1)*(NY-1))*((j-1.5)*(NY-1)-(j-1.5)*(j-1.5));
		rho0_local = ( f[0][j][0]+f[0][j][2]+f[0][j][4]+2.0*(f[0][j][3]+f[0][j][7]+f[0][j][6]))/(1-v_w);
		ru = rho0_local*v_w;
		f[0][j][1]=f[0][j][3]+(2.0/3.0)*ru;
		f[0][j][5]=f[0][j][7]+(1.0/6.0)*ru-(1.0/2.0)*(f[0][j][2]-f[0][j][4]);
		f[0][j][8]=f[0][j][6]+(1.0/6.0)*ru-(1.0/2.0)*(f[0][j][4]-f[0][j][2]);

	}

}





void boundary_Zou_He_pressure(bool north, double rho_n, bool south,double rho_s, bool east, double rho_e, bool west, double rho_w)

// PLEASE SPECIFY THE BOUNDARY LOCATION BEFORE USE
// FOR EXAMPLE: IF SOUTH BOUNDARY IS INCLUDE, THEN THE BOOL VARIABLE SOUTH=1, AND ALSO SPECIFY THE LOCAL RHO
// IF THIS BOUNDARY CONDITION IS NOT GOING TO BE USED, PLEASE SET THE CORRESPONDING BOOL VARIABLE AS 0
// FORMAT:
// (NORTH BOUNDARY MARK,LOCAL RHO VALUE,SOUTH BOUNDARY MARK, LOCAL RHO VALUE,EAST BOUNDARY MARK, LOCAL RHO VALUE
// WEST BOUNDARY MARK, LOCAL RHO VALUE)

{

double ux0,uy0,ru;


	if (north)
      for (int i=0;i<=NX;i++)
	
	{ 	//cout<<"Zou_He_velocity_north"<<endl;
		uy0 = -1.0 + ( f[i][NY][0]+f[i][NY][1]+f[i][NY][3]+2.0*(f[i][NY][2]+f[i][NY][5]+f[i][NY][6]))/rho_n;
		ru = rho_n*uy0;
		f[i][NY][4]=f[i][NY][2]-(2.0/3.0)*ru;
		f[i][NY][7]=f[i][NY][5]-(1.0/6.0)*ru+(1.0/2.0)*(f[i][NY][1]-f[i][NY][3]);
		f[i][NY][8]=f[i][NY][6]-(1.0/6.0)*ru+(1.0/2.0)*(f[i][NY][3]-f[i][NY][1]);

	}


if (south)
      for (int i=0;i<=NX;i++)
	{	//cout<<"Zou_He_velocity_south"<<endl;
		uy0 = -1.0 + ( f[i][0][0]+f[i][0][1]+f[i][0][3]+2.0*(f[i][0][4]+f[i][0][7]+f[i][0][8]))/rho_s;
		ru = rho_s*uy0;
		f[i][0][2]=f[i][0][4]+(2.0/3.0)*ru;
		f[i][0][5]=f[i][0][7]+(1.0/6.0)*ru+(1.0/2.0)*(f[i][0][3]-f[i][0][1]);
		f[i][0][6]=f[i][0][8]+(1.0/6.0)*ru+(1.0/2.0)*(f[i][0][1]-f[i][0][3]);

	}


if (east)
      for (int j=1;j<NY;j++)
	{	//cout<<"Zou_He_velocity_east"<<endl;
		ux0 = -1.0 + ( f[NX][j][0]+f[NX][j][2]+f[NX][j][4]+2.0*(f[NX][j][1]+f[NX][j][5]+f[NX][j][8]))/rho_e;
		ru = rho_e*ux0;
		f[NX][j][3]=f[NX][j][1]-(2.0/3.0)*ru;
		f[NX][j][7]=f[NX][j][5]-(1.0/6.0)*ru+(1.0/2.0)*(f[NX][j][2]-f[NX][j][4]);
		f[NX][j][6]=f[NX][j][8]-(1.0/6.0)*ru+(1.0/2.0)*(f[NX][j][4]-f[NX][j][2]);

	}


if (west)
      for (int j=1;j<NY;j++)
	{	//cout<<"Zou_He_velocity_west"<<endl;
		ux0 = -1.0 + ( f[0][j][0]+f[0][j][2]+f[0][j][4]+2.0*(f[0][j][3]+f[0][j][7]+f[0][j][6]))/rho_w;
		ru = rho_w*ux0;
		f[0][j][1]=f[0][j][3]+(2.0/3.0)*ru;
		f[0][j][5]=f[0][j][7]+(1.0/6.0)*ru+(1.0/2.0)*(f[0][j][4]-f[0][j][2]);
		f[0][j][8]=f[0][j][6]+(1.0/6.0)*ru+(1.0/2.0)*(f[0][j][2]-f[0][j][4]);

	}

}

void evolution()	
{	
	comput_macro_variables(); 
	//collision_and_streaming();
	//boundary_Zou_He_pressure(1,rho_l,0,rho_v,0,0,0,0);
        //boundary_Zou_He_velocity(0,uMax,1,uMax,0,0,0,0);
	collision();//cout<<"mark"<<endl;
	//boundary_Zou_He_pressure(1,rho_l,0,rho_v,0,0,0,0);
	//boundary_Zou_He_velocity(1,uMax,0,uMax,0,0,0,0);
	periodic_streaming();
        //boundary_Zou_He_velocity(0,-uMax,0,-uMax,0,uMax,1,uMax);
	//boundary_Zou_He_pressure(0,rho_l,0,rho_l,1,1,0,0);	
        //boundaries_Guo_LR();
	//boundaries_Guo_TB();
        //periodic_streaming();	
       	//comput_macro_variables();

}

void output(int m)	
{
	ostringstream name;
	name<<"cavity_"<<m<<".dat";
	ofstream out(name.str().c_str());
	//out<<"Title=\"LBM Lid Driven Flow\"\n"
	//	<<"VARIABLES=\"X\",\"Y\",\"U\",\"V\"\n"
	//	<<"ZONE T= \"BOX\", I= "
	//	<<NX+1<<", J="<<NY+1<<", F=POINT"<<endl;
	for(int j=0; j<=NY; j++)
		for(int i=0; i<=NX; i++)
		{
			out<<double(i)/Lx<<" "<<double(j)/Ly<<" "
				<<u[i][j][0]<<" "<<u[i][j][1]<<endl;
		}
}

void output_matrix(int m)
{
	ostringstream name;
	name<<"matrix_"<<m<<".vtk";
	//ofstream out(name.str().c_str());
	ofstream out;
	out.open(name.str().c_str());
	//ofstream out;
	//out.open(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<NX+1<<"         "<<NY+1<<"         "<<1<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<(NX+1)*(NY+1)<<endl;
	out<<"SCALARS sample_scalars double"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	
	for(int i=0; i<=NX; i++)
		{
		for(int j=0; j<=NY; j++)
			out<<psi[i][j]<<" "<<endl;;
		//out<<endl;
		}
	

		
}

void output_velocity(int m)	
{
	ostringstream name;
	name<<"LBM_velocity_"<<m<<".vtk";
	//ofstream out(name.str().c_str());
	ofstream out;
	out.open(name.str().c_str());
	//ofstream out;
	//out.open(name.str().c_str());
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<NX+1<<"         "<<NY+1<<"         "<<1<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<(NX+1)*(NY+1)<<endl;
	out<<"VECTORS sample_vectors double"<<endl;
	//out<<"LOOKUP_TABLE default"<<endl;
	


        for(int j=0;j<=NY;j++)
		{
		for(int i=0; i<=NX; i++)
			{
			out<<u_mac[i][j][0]<<" "<<u_mac[i][j][1]<<" "<<0.0<<endl;
			};
		
		}


	out.close();

}



void output_velocity_x(int m)	
{
	ostringstream name;
	name<<"LBM_velocity_x"<<m<<".dat";
	ofstream out(name.str().c_str());
        for(int i=0;i<=NX;i++)
		{
		for(int j=0; j<=NY; j++)
			{
			out<<u_mac[i][j][0]<<" ";
			};
		out<<endl;
		}
}


void output_velocity_y(int m)	
{
	ostringstream name;
	name<<"LBM_velocity_y"<<m<<".dat";
	ofstream out(name.str().c_str());
        for(int i=0;i<=NX;i++)
		{
		for(int j=0; j<=NY; j++)
			{
			out<<u_mac[i][j][1]<<" ";
			};
		out<<endl;
		}
}


void outputpoi(int m)	
{
	ostringstream name;
	name<<"psi_"<<m<<".dat";
	ofstream out(name.str().c_str());
        int i=30; 
	//for(int i=0;i<=NX;i++)
	//	{
		for(int j=0; j<=NY; j++)
			{if (!Solid[i][j])
				out<<u_R[i][j][0]<<" ";
			 else
				out<<0<<" ";

			};
		out<<endl;
	//	}	



}


void outputpoi2(int m)	
{
	ostringstream name;
	name<<"psi2_"<<m<<".dat";
	ofstream out(name.str().c_str());
        int i=30;
	//for(int i=0;i<=NX;i++)
	//	{
		for(int j=0; j<=NY; j++)
			{if (!Solid[i][j])
				out<<u[i][j][0]<<" ";
			 else
				out<<0<<" ";

			};
		out<<endl;
	//	}	

ostringstream name2;
name2<<"psi2_mac_"<<m<<".dat";
	ofstream out2(name2.str().c_str());
        double tmp;
	
	//for(int i=0;i<=NX;i++)
	//	{
		for(int j=0; j<=NY; j++)
			{if (!Solid[i][j])
				{
				tmp=rho[i][j]*u[i][j][0]+rho_R[i][j]*u_R[i][j][0]+(forcex[i][j]+forcex_R[i][j])/2;
				tmp/=(rho_mac[i][j]);
				//out2<<u_mac[i][j][0]<<" ";
				out2<<tmp<<" ";
				}
			 else
				out2<<0<<" ";

			};
		out2<<endl;
	//	}	


}


void output_density(int m)	
{
	ostringstream name;
	name<<"Prous_"<<m<<".vtk";
	//ofstream out(name.str().c_str());
	ofstream out;
	out.open(name.str().c_str());
	
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<NX+1<<"         "<<NY+1<<"         "<<1<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<(NX+1)*(NY+1)<<endl;
	out<<"SCALARS sample_scalars double"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
        for(int j=0;j<=NY;j++)
		{
		for(int i=0; i<=NX; i++)
			{if (!Solid[i][j])
				out<<rho[i][j]<<" "<<endl;
			 else
				out<<rho0<<" "<<endl;

			};
		}
	out.close();
	

	ostringstream name2;
	name2<<"Prous_2_"<<m<<".vtk";
	
	out.open(name2.str().c_str());
	
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<NX+1<<"         "<<NY+1<<"         "<<1<<endl;
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<(NX+1)*(NY+1)<<endl;
	out<<"SCALARS sample_scalars double"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
        for(int j=0;j<=NY;j++)
		{
		for(int i=0; i<=NX; i++)
			{if (!Solid[i][j])
				out<<rho_R[i][j]<<" "<<endl;
			 else
				out<<rho0<<" "<<endl;

			};
	
		}

	out.close();



		
}

void output_density2(int m)	
{
	ostringstream name;
	name<<"Prous_"<<m<<".dat";
	//ofstream out(name.str().c_str());
	ofstream out;
	out.open(name.str().c_str());
	
	
        for(int i=0;i<=NX;i++)
		{
		for(int j=0; j<=NY; j++)
			{if (!Solid[i][j])
				out<<rho[i][j]<<" ";
			 else
				out<<rho0<<" ";

			};
		out<<endl;
		}
	out.close();
	

	ostringstream name2;
	name2<<"Prous_2_"<<m<<".dat";
	
	out.open(name2.str().c_str());
	
	
        for(int i=0;i<=NX;i++)
		{
		for(int j=0; j<=NY; j++)
			{if (!Solid[i][j])
				out<<rho_R[i][j]<<" ";
			 else
				out<<rho0<<" ";

			};
	out<<endl;
		}

	out.close();



		
}


void Error()
{
	double temp1,temp2;
	temp1=0;
	temp2=0;
	for(int i=1; i<NX; i++)
		for(int j=1;j<NY;j++)
		{
			temp1 += (u_mac[i][j][0]-u0[i][j][0])*(u_mac[i][j][0]-u0[i][j][0])+(u_mac[i][j][1]-u0[i][j][1])*(u_mac[i][j][1]-u0[i][j][1]);
			temp2 += u_mac[i][j][0]*u_mac[i][j][0]+u_mac[i][j][1]*u_mac[i][j][1];
		}
		temp1=sqrt(temp1);
		temp2=sqrt(temp2);
		error=temp1/(temp2+1e-30);  
}
