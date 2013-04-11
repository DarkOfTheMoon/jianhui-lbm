#include <stdio.h>

#include <dirent.h>

#include<iostream>
#include<cmath>
#include<cstdlib>
#include<iomanip>
#include<fstream>
#include<sstream>
#include<string>
 using namespace std; 


int main( int argc, char *argv[] )

{


string nametest;
int where;


         struct dirent *pDirEntry = NULL;

         DIR          *pDir      = NULL;

         if( (pDir = opendir("/home/jy810/LBM/source/CODE/jianhui-lbm/PreProcessing/")) == NULL )

         {

                   printf("opendir failed!\n");

                   return 1;

         }

         else

         {

                   while( pDirEntry = readdir(pDir) )

                   {
			nametest=pDirEntry->d_name;
			where = nametest.find(".vtk");cout<<where<<endl;
                            printf("索引节点:%d\t 偏移量：%d\t 文件名长：%d\t文件类型：%d\t 文件名：%s\n",

                            pDirEntry->d_ino, pDirEntry->d_off,

                            pDirEntry->d_reclen,pDirEntry->d_type,pDirEntry->d_name);

                   }

                   closedir(pDir);

                   return 0;    

         }       

}
