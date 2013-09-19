clear
figure

area_v=[];


sizes=[120,160,220,290,350,400];
file1='1r_Sub_Sample_Perm_120x60x60_125.sta_dat';
file2='1r_Sub_Sample_Perm_160x80x80_64.sta_dat';
file3='1r_Sub_Sample_Perm_220x110x110_27.sta_dat';
file4='1r_Sub_Sample_Perm_290x145x145_8.sta_dat';
file5='1r_Sub_Sample_Perm_350x175x175_8.sta_dat';
file6='1r_Sub_Sample_Perm_400x200x200_8.sta_dat';

ALL_size=230;
All_Perm=0.166;
All_Porosity=0.23427;

component=2;  %component 1 or 2
%==================================================
a1=load(file1);

xx=a1(:,1);
yy=a1(:,component+1);
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
yy=a1(:,component+1);
[Perm_42,IX_42]=sort(a1',2);
a1=a1';
x=xx;y=yy;
k = convhull(x,y);
plot(x(k),y(k),'k-',x,y,'^k');
[K,V] = convhull(x,y);area_v=[area_v,V];
hold on
plot (a1(1,IX_42(1,:)),a1(2,IX_42(1,:)),'-ks');
hold on
%===================================================

%==================================================
a1=load(file3);

xx=a1(:,1);
yy=a1(:,component+1);
[Perm_42,IX_42]=sort(a1',2);
a1=a1';
x=xx;y=yy;
k = convhull(x,y);
plot(x(k),y(k),'g-',x,y,'sg');
[K,V] = convhull(x,y);area_v=[area_v,V];
hold on
plot (a1(1,IX_42(1,:)),a1(2,IX_42(1,:)),'-gs');
hold on
%===================================================


%==================================================
a1=load(file4);

xx=a1(:,1);
yy=a1(:,component+1);
[Perm_42,IX_42]=sort(a1',2);
a1=a1';
x=xx;y=yy;
k = convhull(x,y);
plot(x(k),y(k),'y-',x,y,'dy');
[K,V] = convhull(x,y);area_v=[area_v,V];
hold on
plot (a1(1,IX_42(1,:)),a1(2,IX_42(1,:)),'-ys');
hold on
%===================================================

%==================================================
a1=load(file5);

xx=a1(:,1);
yy=a1(:,component+1);
[Perm_42,IX_42]=sort(a1',2);
a1=a1';
x=xx;y=yy;
k = convhull(x,y);
plot(x(k),y(k),'c-',x,y,'xc');
[K,V] = convhull(x,y);area_v=[area_v,V];
hold on
plot (a1(1,IX_42(1,:)),a1(2,IX_42(1,:)),'-cs');
hold on
%===================================================

%==================================================
a1=load(file6);

xx=a1(:,1);
yy=a1(:,component+1);
[Perm_42,IX_42]=sort(a1',2);
a1=a1';
x=xx;y=yy;
k = convhull(x,y);
plot(x(k),y(k),'r-',x,y,'>r');
[K,V] = convhull(x,y);area_v=[area_v,V];
hold on
plot (a1(1,IX_42(1,:)),a1(2,IX_42(1,:)),'-rs');
hold on
%===================================================


plot (All_Porosity,All_Perm,'-ro');

hold off


figure
plot(sizes,area_v,'r^');


%plot (perm_90(2,IX_90(2,:)),perm_90(1,IX_90(2,:)),'kd');