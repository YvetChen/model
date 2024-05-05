
clc
clear all %#ok<*CLALL>
close all

%The unit is mm.

a=0.5;% The semimajor axis length of the ellipse
mu=0.01;% Aspect ratio
b=a*mu;% The semiminor axis length of the ellipse
numc=135;% Number of circles
nume=120;% Number of ellipses
r=0.6;% Radius of a circle
P=25e6;

kg=75.47e9;
ug=32.22e9;

kd=33.4e9;
ud=15.8e9;

w=2*pi*1e6;
n=1e-3; % The viscosity of water

kf=2.18e9;

stp=1+3*kg/(4*ug);
ctp=kg*(3*kg+4*ug)/(pi*mu*ug*(3*kg+ug));
tp=(pi*a*b*nume+pi*r*r*numc)/(25*50);% Porosity

c0=pi*a*b*nume/(25*50);% Initial porosity of the conpliant holes

syms  kh1


m=1-stp/kg;
l=ctp*c0;
g=-ctp*P;
j=stp-kd;
f=kh1*(1-(stp*P/kg)-l*exp(g/kh1))-(kd+stp*P)==0;

kh=vpasolve(f,kh1);
kh=double(kh);


c=c0*exp(-ctp*P*(1/kh));

kmf1=1/kh+1/(1/(1/kd-1/kh)+3*i*w*n/(8*c*mu*mu));

kmf=1/kmf1;

umf1=1/ud-(4/15)*(1/kd-1/kmf);

umf=1/umf1;

ks1=1/kg+tp*(1/kf-1/kg)/(1+tp*(1/kf-1/kg)/(1/kmf-1/kg));

ks=1/ks1;

us=umf;

den=2750;
vs0=sqrt(us/den);
vp0=sqrt(ks/den+(4/3)*vs0*vs0);


sat = fopen('bouk modulus.txt','wt+');

fprintf(sat,'ks=%g[Pa]\n',ks); 
fprintf(sat,'us=%g[Pa]\n',us);  
fprintf(sat,'c0=%g\n',c0); 
fprintf(sat,'kh=%g[Pa]\n',kh);
fprintf(sat,'kmf=%g[Pa]\n',kmf);
fprintf(sat,'c=%g\n',c);
fprintf(sat,'c0=%g\n',c0);
fprintf(sat,'tp=%g\n',tp);
fprintf(sat,'vp0=%g\n',vp0);
       
fclose(sat);


