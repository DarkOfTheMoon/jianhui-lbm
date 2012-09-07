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

int nx=370;
int ny=370;
int nz=1;
int mode=2; //1=two phases, 2=three phases
char poreFileName[128]="2phasetiff.txt";
char poreFileNameVTK[128]="segment.vtk";
char poreFileNameOut[128]="segment.dat";
//output VTK file,0=no, 1=yes
int VTK_OUT=1;
int ImageJ_sequence=1;
int local_min=5; //half of local min
int peak_check_length=6;
int histo_find_accuricy=20;
int area_max_radius=20;
int area_max_fluc=30;
//double solid_part_shift_percentage=0.3;
//===========VTK AND OUT ARE ALL WRITTEN IN BINARY FORMAT===============
//===========================================================





int histo[256];
double*** Solid;
double*** cas;
double*** cas2;
double*** cas3;
int*** seg;
double his1[256];
double his2[256];
double his3[256];
double his4[256];
int his5[256];
int peak1,peak2,peak3;
int peak1s,peak2s,peak3s;
int peak1s_cor,peak2s_cor,peak3s_cor;
int val1,val2,val3;
int sum1,sum2;
double max1,max2,max3;

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


int pore;
	int i, j, k,ci,cj,ck;
	
	
	Solid = new double**[nx];
	cas = new double**[nx];cas2 = new double**[nx];cas3 = new double**[nx];
	seg= new int**[nx];
		
	for (i=0; i<nx;i++)
	{
	       Solid[i] = new double*[ny];
		cas[i] = new double*[ny];cas2[i] = new double*[ny];cas3[i] = new double*[ny];
		seg[i] = new int*[ny];

	       for (j=0;j<ny;j++)
	       {
	               Solid[i][j] = new double[nz];
			cas[i][j] = new double[nz];cas2[i][j] = new double[nz];cas3[i][j] = new double[nz];
			seg[i][j] = new int[nz];


	               for (k=0;k<nz;k++)
				{
	                       	Solid[i][j][k] = 0;
				cas[i][j][k] = 0;cas2[i][j][k] = 0;cas3[i][j][k] = 0;
				seg[i][j][k] = 0;
				}
	       }
	}
		
	

//for (i=0;i<256;i++)
//	{histo[i]=0;his1[i]=0.0;his2[i]=0.0;his4[i]=0.0;}



	for(k=0 ; k<nz ; k++)		
	{        
	 for (i=0;i<256;i++)
	{histo[i]=0;his1[i]=0.0;his2[i]=0.0;his4[i]=0.0;his5[i];}       
	        
	        
	for(j=0 ; j<ny ; j++)
	for(i=0 ; i<nx ; i++)				///*********


	//while (!fin.eof())                                        //**********
		{	
			//fin >> ci >> cj>> ck>>pore;
			fin >> pore;
			
			//if (pore == 0.0)	{Solid[ci-1][cj-1][ck-1] = 0;}
			//if (pore == 0.0)	{Solid[i][j][k] = 0;}
			
			//if (pore == 1.0) 	{Solid[ci-1][cj-1][ck-1] = 1;sum++;}
			//if (pore == 1.0) 	{Solid[i][j][k] = 1;sum++;}
			
			Solid[i][j][k]=(double) pore;
			
			seg[i][j][k]=pore;
			
			histo[pore]+=1;
		}
		
	//fin.close();
		
	for (i=1;i<255;i++)
		{
		//his1[i]=histo[i+1]-histo[i];
		his2[i]=(histo[i+1]-histo[i-1])/2;
		//his5[i]=(histo[i+1]+histo[i-1]-2*histo[i]);
		}
	
	for (i=1;i<255;i++)
	//for (i=2;i<254;i++)
		his1[i]=(his2[i-1]+his2[i]+his2[i+1])/3;
		//his1[i]=(his2[i-2]+his2[i-1]+his2[i]+his2[i+1]+his2[i+2])/5;


	for (i=0;i<256;i++)
		{sum1=0;
		for (j=-area_max_radius;j<=area_max_radius;j++)
			if ((i+j>=0) and (i+j<256) and (histo[i+j]>sum1))
				sum1=histo[i+j];
		his5[i]=sum1;
		}


	

double max_val=0.0;	
double max_val2=0.0;



	//----------peak search mode 1--------------------
	/*
	for (i=0;i<256;i++)
		{
		if (histo[i]>max_val)
			max_val=histo[i];
		//if (histo[255-i]>max_val2)
		//	max_val2=histo[i];
		his3[i]=max_val;//his4[i]=max_val2;
		}


		        peak1=0;peak2=0;peak3=0;
		for (i=0;i<256-peak_check_length;i++)
		          if ((his3[i]==histo[i]) and (his3[i]>0))
		          {        sum1=0;
		                  for (j=1;j<peak_check_length;j++)
		                          if (his3[j+i]!=his3[i])
		                          sum1=1;
		                  if (sum1==0)
		                          if (peak1==0)
		                                  peak1=i;
		                          else
		                                  if (peak2==0)
		                                          peak2=i;
		                                  else
		                                          peak3=i;
		          }
	*/
	//------------------------------------------------





	//-----------peak search mode2-----------------
	
	for (i=0;i<256;i++)
		{
		if (his5[i]>max_val)
			max_val=his5[i];
		//if (histo[255-i]>max_val2)
		//	max_val2=histo[i];
		his3[i]=max_val;//his4[i]=max_val2;
		}
	peak1s=0;peak2s=0;peak3s=0;
	
	for (i=peak_check_length;i<256-peak_check_length;i++)
		{sum1=0;
		for (j=-peak_check_length;j<=peak_check_length;j++)
			if (his5[i]!=his5[i+j])
				sum1=1;
		if (sum1==0)
			{//cout<<his5[i]<<"	@@@@@@@"<<endl;
			if (his5[i]>peak1s)
				{
				peak3s=peak2s;peak2s=peak1s;peak1s=his5[i];
				peak3s_cor=peak2s_cor;peak2s_cor=peak1s_cor;peak1s_cor=i;
				}
				else
				if ((his5[i]>peak2s) and (his5[i]<peak1s))
					{
					peak3s=peak2s;peak2s=his5[i];
					peak3s_cor=peak2s_cor;peak2s_cor=i;
					}
					else
					if ((his5[i]>peak3s) and (his5[i]<peak2s))
						{
						peak3s=his5[i];
						peak3s_cor=i;
						}
			}
		}
	//cout<<"asdfasdfas       "<<peak1s<<"	"<<peak2s<<"	"<<peak3s<<endl;


	if (peak3s_cor>peak2s_cor)
		{sum1=peak3s_cor;peak3s_cor=peak2s_cor;peak2s_cor=sum1;}
	if (peak2s_cor>peak1s_cor)
		{sum1=peak2s_cor;peak2s_cor=peak1s_cor;peak1s_cor=sum1;}
	if (peak3s_cor>peak2s_cor)
		{sum1=peak3s_cor;peak3s_cor=peak2s_cor;peak2s_cor=sum1;}


    	peak1=peak3s_cor;peak2=peak2s_cor;peak3=peak1s_cor;
	/*	for (i=0;i<256;i++)
		          if (((peak1s==histo[i]) or (peak2s==histo[i]) or (peak3s==histo[i])))
		          {        cout<<histo[i]<<"	"<<i<<endl;
		                          if (peak1==0)
		                                  peak1=i;
		                          else
		                                  if (peak2==0)
		                                          peak2=i;
		                                  else
		                                          peak3=i;
		          }
	*/
	//---------------------------------------------



	for (i=local_min;i<256-local_min;i++)
		{
		sum1=7000000;
		for (j=i-local_min;j<i+local_min;j++)
			if (histo[j]<sum1)
				sum1=histo[j];	

		his4[i]=sum1;
		}

		//========FIND PEAKS  =================
		/*
		        peak1=0;peak2=0;peak3=0;
		for (i=0;i<256-peak_check_length;i++)
		          if ((his3[i]==histo[i]) and (his3[i]>0))
		          {        sum1=0;
		                  for (j=1;j<peak_check_length;j++)
		                          if (his3[j+i]!=his3[i])
		                          sum1=1;
		                  if (sum1==0)
		                          if (peak1==0)
		                                  peak1=i;
		                          else
		                                  if (peak2==0)
		                                          peak2=i;
		                                  else
		                                          peak3=i;
		          }
*/
		//===================================================
		
		
		cout<<"========PEAKS==============="<<k<<endl;
		cout<<peak1<<"        "<<peak2<<"        "<<peak3<<endl;
		//cout<<"============================"<<endl;
		
		for (i=0;i<256;i++)
		        his2[i]=his4[i]-histo[i];
		
		val1=256;val2=256;val3=256;
		for (i=0;i<256;i++)
		        {
		                if ((his2[i]==0) and (abs(i-(peak1+(peak2-peak1)/2))<histo_find_accuricy))
		                        val1=i;
		                
		                if ((his2[i]==0) and (abs(i-(peak2+(peak3-peak2)/2))<histo_find_accuricy))
		                        val2=i;
		        }
		
		sum1=0;
		for (i=-local_min;i<local_min+1;i++)
			if (histo[val1+i]>sum1)
				sum1=histo[val1+i];
		//cout<<sum1<<endl;

		sum2=-1;
		for (i=val2;i>0;i--)
			if ((sum2<0) and (histo[i]<sum1))
				sum2=i;
		val1=sum2;
			

		//val1=val1+int((val2-val1)*solid_part_shift_percentage);
	
	cout<<"===========Segment Vals================"<<endl;
	cout<<val1<<"	"<<val2<<endl;
	cout<<"======================================="<<endl;
	cout<<endl;
	
	for(j=0 ; j<ny ; j++)
	for(i=0 ; i<nx ; i++)	
	if (seg[i][j][k]<val1)
	        seg[i][j][k]=0;
	        else
	                if (seg[i][j][k]<val2)
	                seg[i][j][k]=1;
	                else
	                       seg[i][j][k]=2;
		        
		        
	        
		        
		        
		        
		        
		        
	}	
	//=========================================
fin.close();



	
	if (VTK_OUT==1)
	{
	cout<<"Start writing VTK file"<<endl;
	cout<<nx<<"         "<<ny<<"         "<<nz<<endl;

	ostringstream name;
	name<<poreFileNameVTK;
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
	ofstream out;
	out.open(name.str().c_str());
	
	
	out<<"# vtk DataFile Version 2.0"<<endl;
	out<<"J.Yang Lattice Boltzmann Simulation 3D Single Phase-Solid-Density"<<endl;
	out<<"ASCII"<<endl;
	out<<"DATASET STRUCTURED_POINTS"<<endl;
	out<<"DIMENSIONS         "<<nx<<"         "<<ny<<"         "<<nz<<endl;       ///*********
	out<<"ORIGIN 0 0 0"<<endl;
	out<<"SPACING 1 1 1"<<endl;
	out<<"POINT_DATA     "<<nx*ny*nz<<endl;				///*********
	out<<"SCALARS sample_scalars float"<<endl;
	out<<"LOOKUP_TABLE default"<<endl;
	
	
	for (k=0;k<nz;k++)
	{
	cout<<k<<endl;
	for (j=0;j<ny;j++)
	for (i=0;i<nx;i++)
/*	
        if (seg[i][j][k]<val1)
	        out<<0<<" ";
	        else
	                if (seg[i][j][k]<val2)
	                out<<1<<" ";
	                else
	                        out<<2<<" ";
*/	  
		out<<seg[i][j][k]<<" ";
	}
	

	//out.write((char *)(&seg[0][0][0]), sizeof(int)*nx*ny*nz); 

	out.close();

	cout<<"VTK file ouput COMPLETE"<<endl;
	cout<<endl;
	}


	//=================================================================
	cout<<"Start writing cal DAT file"<<endl;
	cout<<nx<<"         "<<ny<<"         "<<nz<<endl;
	ostringstream name2;
	name2<<poreFileNameOut;
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
	ofstream out2;
	out2.open(name2.str().c_str());
	
	for (k=0;k<nz;k++)
	{
	cout<<k<<endl;
	for (j=0;j<ny;j++)
	{
		for (i=0;i<nx;i++)
		out2<<seg[i][j][k]<<" ";
		out2<<endl;
	}
	//out2.write((char *)(&seg[0][0][0]), sizeof(int)*nx*ny*nz); 
	}
	out2.close();
	cout<<"DAT file ouput complete"<<endl;
	
	cout<<endl;
	cout<<"----------------------"<<endl;
	cout<<"Start writing TXT file"<<endl;
	
	
	for (k=0;k<nz;k++)
	{
	cout<<k<<endl;
	name2.str("");
	name2<<"Segment_TXT_seqence_"<<k<<".txt";
	
	out2.open(name2.str().c_str());
	

	for (j=0;j<ny;j++)
	{
		for (i=0;i<nx;i++)
		out2<<seg[i][j][k]<<" ";
		out2<<endl;
	}
	out2.close();
	}
	

	cout<<"TXT files ouput complete"<<endl;
	cout<<"----------------------"<<endl;
	cout<<endl;






	ostringstream name3;
	name3<<"histo_test.dat";
	//name<<"Clashach_z_sym_196x196x388_8.946.dat";
	ofstream out3;
	out3.open(name3.str().c_str());
	
	for (i=0;i<256;i++)
		out3<<histo[i]<<endl;
	
	out3.close();

	name3.str("");
	name3<<"histo_test2.dat";
	out3.open(name3.str().c_str());
	for (i=0;i<256;i++)
		out3<<his1[i]<<endl;
	
	out3.close();


	name3.str("");
	name3<<"histo_test3.dat";
	out3.open(name3.str().c_str());
	for (i=0;i<256;i++)
		out3<<his5[i]<<endl;
	
	out3.close();

	name3.str("");
	name3<<"histo_test4.dat";
	out3.open(name3.str().c_str());
	for (i=0;i<256;i++)
		out3<<his4[i]<<endl;
	
	out3.close();
			
}


