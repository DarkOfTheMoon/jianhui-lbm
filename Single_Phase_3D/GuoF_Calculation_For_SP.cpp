
#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
#include<math.h> 


#define MAXN 100
#define eps 1e-12
#define zero(x) (fabs(x)<eps)


using namespace std;     



double elat[19][3]=
{{0,0,0},{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1},{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},{1,0,1},
{-1,0,-1},{1,0,-1},{-1,0,1},{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1}};

double M[19][19]=
{{1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1},
{-30,-11,-11,-11,-11,-11,-11,8,8,8,8,8,8,8,8,8,8,8,8},
{12,-4,-4,-4,-4,-4,-4,1,1,1,1,1,1,1,1,1,1,1,1},
{0,1,-1,0,0,0,0,1,-1,1,-1,0,0,0,0,1,-1,1,-1},
{0,-4,4,0,0,0,0,1,-1,1,-1,0,0,0,0,1,-1,1,-1},
{0,0,0,1,-1,0,0,1,1,-1,-1,1,-1,1,-1,0,0,0,0},
{0,0,0,-4,4,0,0,1,1,-1,-1,1,-1,1,-1,0,0,0,0},
{0,0,0,0,0,1,-1,0,0,0,0,1,1,-1,-1,1,1,-1,-1},
{0,0,0,0,0,-4,4,0,0,0,0,1,1,-1,-1,1,1,-1,-1},
{0,2,2,-1,-1,-1,-1,1,1,1,1,-2,-2,-2,-2,1,1,1,1},
{0,-4,-4,2,2,2,2,1,1,1,1,-2,-2,-2,-2,1,1,1,1},
{0,0,0,1,1,-1,-1,1,1,1,1,0,0,0,0,-1,-1,-1,-1},
{0,0,0,-2,-2,2,2,1,1,1,1,0,0,0,0,-1,-1,-1,-1},
{0,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,-1,-1,1},
{0,0,0,0,0,0,0,1,-1,1,-1,0,0,0,0,-1,1,-1,1},
{0,0,0,0,0,0,0,-1,-1,1,1,1,-1,1,-1,0,0,0,0},
{0,0,0,0,0,0,0,0,0,0,0,-1,-1,1,1,1,1,-1,-1}};

struct mat{
    int n,m;
    double data[MAXN][MAXN];
};


double MI[19][19];


void Comput_MI(double[19][19], double[19][19]);

int inverse(mat &a);

int main(int argc , char *argv [])
{
//char* str1[]={"hello","world","test char array"};
string str1[]={"hello","world","test char array"};
 int dec, sign; 
   int ndig = 3; 
Comput_MI(M,MI);
 
 double  value = -92.876535; 
   str1[2] = fcvt(value, ndig, &dec, &sign);
string str[19][3];
string strinv[19][19];
//cout<<str1[2]<<"    "<<dec<<"   "<<sign<<endl;

str1[2].insert(dec,".");


for (int i=0;i<19;i++)
        for (int j=0;j<3;j++)
        {
                str[i][j]=fcvt(elat[i][j],ndig,&dec,&sign);
                
                if (dec>0)
                        str[i][j].insert(dec,".");
                if (dec==0)
                        str[i][j].insert(1,".");
		if (dec<0)
			{
			for (int k=1;k<=-dec+1;k++)
				str[i][j].insert(0,"0");
			str[i][j].insert(1,".");
			}
              if (sign==1)
                      str[i][j].insert(0,"-");

		if (sign==0)
                      str[i][j].insert(0,"+");

               str[i][j]=str[i][j]+"*c_l";
                
        }
/*
for (int k=0;k<19;k++)
			{	
			lm0=((elat[k][0]-u[ci][0]/PorC[ci])*forcex[ci]+(elat[k][1]-u[ci][1]/PorC[ci])*forcey[ci]+(elat[k][2]-u[ci][2]/PorC[ci])*forcez[ci])/c_s2;
			lm1=(elat[k][0]*u[ci][0]+elat[k][1]*u[ci][1]+elat[k][2]*u[ci][2])*(elat[k][0]*forcex[ci]+elat[k][1]*forcey[ci]+elat[k][2]*forcez[ci])/(c_s2*c_s2);
			lm1/=PorC[ci];
			GuoF[k]=w[k]*(lm0+lm1);
			//GuoF[k]=0.0;



if (lattice_v==1)
	for (int i=0;i<19;i++)
		for (int j=0;j<3;j++)
		elat[i][j]=e[i][j]*lat_c;
					
*/



        ostringstream name;
	name<<"GuoF"<<".dat";
	ofstream out;
	out.open(name.str().c_str());
//	out<<"        for (int mi=0; mi<19; mi++)"<<endl;
	out<<"//=========================="<<endl;

	for (int mi=0; mi<19; mi++)
		{
		out<<"lm0=";
		out<<"(("<<str[mi][0]<<"-u[ci][0])*gxf+("<<str[mi][1]<<"-u[ci][1])*gyf+("<<str[mi][2]<<"-u[ci][2])*gzf)/c_s2;"<<endl;
		out<<"lm1=";
		out<<"("<<str[mi][0]<<"*u[ci][0]+("<<str[mi][1]<<")*u[ci][1]+("<<str[mi][2]<<")*u[ci][2])*("<<str[mi][0]<<"*gxf+("<<str[mi][1]<<")*gyf+("<<str[mi][2]<<")*gzf)/(c_s2*c_s2);"<<endl;
		//=============
		//out<<"lm1/=PorC[ci];"<<endl;
		//==============
		out<<"GuoF["<<mi<<"]=w["<<mi<<"]*(lm0+lm1);"<<endl;
		out<<endl;
		

		}
out<<"//===================="<<endl;


/*
	
	//sum=0;	
	//		for (int mj=0; mj<19; mj++)
	//			sum+=MI[mi][mj]*m_l[mj];
	

	out<<"//===================="<<endl;
	for (int mi=0; mi<19; mi++)
		{out<<"sum["<<mi<<"]=";
		for (int mj=0; mj<19; mj++)
			out<<strinv[mi][mj]<<"*m_l["<<mj<<"]";
		out<<";"<<endl;
		out<<endl;
		}

*/



	out.close();				
					
					
}






int inverse(mat &a){
    double t;
    int i,j,k,is[MAXN],js[MAXN];
    if(a.n!=a.m) return 0;
    for(k=0;k<a.n;k++){
        for(t=0,i=k;i<a.n;i++)
            for(j=k;j<a.n;j++)
                if(fabs(a.data[i][j])>t)
                    t=fabs(a.data[is[k]=i][js[k]=j]);
        if(zero(t)) return 0;
        if(is[k]!=k)
            for(j=0;j<a.n;j++)
                t=a.data[k][j],a.data[k][j]=a.data[is[k]][j],a.data[is[k]][j]=t;
        if(js[k]!=k)
            for(i=0;i<a.n;i++)
                t=a.data[i][k],a.data[i][k]=a.data[i][js[k]],a.data[i][js[k]]=t;
        a.data[k][k]=1/a.data[k][k];
        for(j=0;j<a.n;j++)
            if(j!=k)
                a.data[k][j]*=a.data[k][k];
        for(i=0;i<a.n;i++)
            if(i!=k)
                for(j=0;j<a.n;j++)
                    if(j!=k)
                        a.data[i][j]-=a.data[i][k]*a.data[k][j];
        for(i=0;i<a.n;i++)
            if(i!=k)
                a.data[i][k]*=-a.data[k][k];
    }
    for(k=a.n-1;k>=0;k--){
        for(j=0;j<a.n;j++)
            if(js[k]!=k)
                t=a.data[k][j],a.data[k][j]=a.data[js[k]][j],a.data[js[k]][j]=t;
        for(i=0;i<a.n;i++)
            if(is[k]!=k)
                t=a.data[i][k],a.data[i][k]=a.data[i][is[k]],a.data[i][is[k]]=t;
    }
    return 1;
}

void Comput_MI(double M[19][19], double MI[19][19])
{

double mim[19][19];

mat a;
    int i,j;
    int n_s=19;
        for(int i=0;i<n_s;i++)
            for(int j=0;j<n_s;j++)
                a.data[i][j]=M[i][j];
	a.m=a.n=n_s;
        if(inverse(a))
            for(int i=0;i<n_s;i++)
                for(int j=0;j<n_s;j++)
                    MI[i][j]=a.data[i][j];
               
            
        else
            puts("NO");


}
