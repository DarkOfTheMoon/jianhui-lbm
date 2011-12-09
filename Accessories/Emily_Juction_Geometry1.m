a=10;
b=20;
c=10;
d=30;
r=3;
r2=4;


A=ones(d*2,a+b+c);

for i=1:d*2
    for j=c-r:c+r
        A(i,j)=0;
    end
end


for j=1:a+b+c
    for i=d-r:d+r
        A(i,j)=0;
    end
end

for j=a+c:a+b+c
    for i=1:d*2
        A(i,j)=0;
    end
end


for j=1:c
    for i=1:3
        A(i,j)=0;
        A(2*d+1-i,j)=0;
    end
end

for i=4:d-r-1
    for j=c-r:c-r+0  %change pertubation
        A(i,j)=1;
    end
end


% =================GEOMETRY2=====================
for i=d-r2:d+r2
    for j=c-r2:c+r2
        A(i,j)=0;
    end
end

fid = fopen('Junction.dat','wt');

zlen=round(r/8*5*2+2);
for k=1:zlen
    for j=1:a+b+c
        for i=1:d*2
            if (((k==1) || (k==zlen)) && (j<a+c))
                     fprintf(fid,'%1d\n',1);
            else
                fprintf(fid,'%1d\n',A(i,j));
            end
            
        end
    end
end



x=d*2
y=a+b+c
z=r*2+1+2
