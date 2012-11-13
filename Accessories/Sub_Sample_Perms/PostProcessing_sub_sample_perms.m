clear
figure

area_v=[];


sizes=[80,120,160,180,200,230];
file1='port_z_Sub_Sample_Perm_80.sta_dat';
file2='port_z_Sub_Sample_Perm_120.sta_dat';
file3='port_z_Sub_Sample_Perm_160.sta_dat';
file4='port_z_Sub_Sample_Perm_180.sta_dat';
file5='port_z_Sub_Sample_Perm_200.sta_dat';
file6='port_z_Sub_Sample_Perm_230.sta_dat';


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


hold off


figure
plot(sizes,area_v,'r^');


%plot (perm_90(2,IX_90(2,:)),perm_90(1,IX_90(2,:)),'kd');