input=imread('circle_arms_bmp1.png','png');
zlen=13;

%===================================================
input=input/255;%====ONLY FOR PNG FILE READING======
%===================================================

fid = fopen('circle_arms_geometry.dat','wt');
if (fid==1)
  error('Cannot open image file...press CTRL-C to exit ');
end

inputf=input(:,:,1);
inputf=inputf';
[x,y]=size(inputf);




% input=input./input;
% a=input(:,:,1);

% for i=1:5
% 
% input1=[input,fliplr(input)];
% input2=flipud(input1);
% inputf=[input1;input2];
% [x,y]=size(inputf);
% 
% input=inputf;
% 
% end




for k=1:zlen
for j=1:y
   for i=1:x
       if (((k==1) || (k==zlen)) && (i>141))
           fprintf(fid,'%1d\n',1);
       else
       fprintf(fid,'%1d\n',inputf(i,j));
       end
   end
  
end
end

x
y
zlen



fclose(fid);
