clear
figure

area_v=[];


sizes=[70,100,130,180,240,290];
file1='HP2_Sub_Sample_Perm_70x30x29.sta_dat';
file2='HP2_Sub_Sample_Perm_100x43x42.sta_dat';
file3='HP2_Sub_Sample_Perm_130x56x55.sta_dat';
file4='HP2_Sub_Sample_Perm_180x78x76.sta_dat';
file5='HP2_Sub_Sample_Perm_240x104x101.sta_dat';
file6='HP2_Sub_Sample_Perm_290x126x122.sta_dat';

ALL_size=350;
All_Perm=5333;
All_Porosity=0.826037;
%==================================================
a1=load(file1);

xx=a1(:,1);
yy=a1(:,2);
[Perm_42,IX_42]=sort(a1',2);
a1=a1';
x=xx;y=yy;
k = convhull(x,y);
plot(x(k),y(k),'b-',x,y,'*b');
[K,V] = convhull(x,y);area_v=[area_v,V];
hold on
plot (a1(1,IX_42(1,:)),a1(2,IX_42(1,:)),'-bs');
hold on
%===================================================

%==================================================
a1=load(file2);

xx=a1(:,1);
yy=a1(:,2);
[Perm_42,IX_42]=sort(a1',2);
a1=a1';
x=xx;y=yy;
k = convhull(x,y);
plot(x(k),y(k),'b-',x,y,'*b');
[K,V] = convhull(x,y);area_v=[area_v,V];
hold on
plot (a1(1,IX_42(1,:)),a1(2,IX_42(1,:)),'-bs');
hold on
%===================================================

%==================================================
a1=load(file3);

xx=a1(:,1);
yy=a1(:,2);
[Perm_42,IX_42]=sort(a1',2);
a1=a1';
x=xx;y=yy;
k = convhull(x,y);
plot(x(k),y(k),'b-',x,y,'*b');
[K,V] = convhull(x,y);area_v=[area_v,V];
hold on
plot (a1(1,IX_42(1,:)),a1(2,IX_42(1,:)),'-bs');
hold on
%===================================================


%==================================================
a1=load(file4);

xx=a1(:,1);
yy=a1(:,2);
[Perm_42,IX_42]=sort(a1',2);
a1=a1';
x=xx;y=yy;
k = convhull(x,y);
plot(x(k),y(k),'b-',x,y,'*b');
[K,V] = convhull(x,y);area_v=[area_v,V];
hold on
plot (a1(1,IX_42(1,:)),a1(2,IX_42(1,:)),'-bs');
hold on
%===================================================

%==================================================
a1=load(file5);

xx=a1(:,1);
yy=a1(:,2);
[Perm_42,IX_42]=sort(a1',2);
a1=a1';
x=xx;y=yy;
k = convhull(x,y);
plot(x(k),y(k),'b-',x,y,'*b');
[K,V] = convhull(x,y);area_v=[area_v,V];
hold on
plot (a1(1,IX_42(1,:)),a1(2,IX_42(1,:)),'-bs');
hold on
%===================================================

%==================================================
a1=load(file6);

xx=a1(:,1);
yy=a1(:,2);
[Perm_42,IX_42]=sort(a1',2);
a1=a1';
x=xx;y=yy;
k = convhull(x,y);
plot(x(k),y(k),'b-',x,y,'*b');
[K,V] = convhull(x,y);area_v=[area_v,V];
hold on
plot (a1(1,IX_42(1,:)),a1(2,IX_42(1,:)),'-bs');
hold on
%===================================================


plot (All_Porosity,All_Perm,'-ro');

hold off


figure
plot(sizes,area_v,'r^');


%plot (perm_90(2,IX_90(2,:)),perm_90(1,IX_90(2,:)),'kd');