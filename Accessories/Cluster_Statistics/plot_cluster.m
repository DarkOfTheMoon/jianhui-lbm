a=load('20-3-3.st');
loglog(a(:,1),a(:,4),'r^');
amin=min(a(:,1));
amax=max(a(:,1));

fx=log(a(:,1));
fy=log(a(:,4));
af=polyfit(fx,fy,1);
tau=-af(1)+2

% fz=polyval(af,fx);
% fx=exp(fx);
% fz=exp(fz);
% hold on 
% loglog (fx,fz,'k*');

% hold on
% loglog(a(:,1),a(:,5),'b^');

% 
% 
% 
% hold on
% x=amin:amax;
% y=x.^(-2.12);
% loglog(x,y,'-b');
 hold off
