input=imread('SD-Final.bmp','bmp');

%===================================================
%input=input/255;%====ONLY FOR PNG FILE READING======
%===================================================

%fid = fopen('Rocks.dat','wt');
%if (fid==1)
%   error('Cannot open image file...press CTRL-C to exit ');
%end

for i=1:5

input1=[input,fliplr(input)];
input2=flipud(input1);
inputf=[input1;input2];
[x,y]=size(inputf);

input=inputf;

end





%for i=1:x
%    for j=1:y
 %       fprintf(fid,'%1d ',inputf(i,j));
 %   end
 %   fprintf(fid,'\n');
%end


%fclose(fid);
