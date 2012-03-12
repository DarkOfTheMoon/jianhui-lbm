function test
close all; clear; clc
% The following command will find all the files in the directory which their extention is *.BMP.
count = size(dir(['*.bmp']), 1)
FNAME = dir(['*.bmp'])
C=[];
for j=1:count
imgname = FNAME(j).name
% The following command will read the files obtained from the previous step.
Imod=imread(imgname);
C=[C;Imod]; 
end
% The last step will save the processed images as an ASCII file. This file will be the essential input for simulators.
 dlmwrite('Filename.txt',C,'delimiter', '\t'); clear imgname; 
end

