function h = core2(x,c,rpa,time)
% Function that computes phi1, phi2 and jLi equations
% for both initialization and integration.

global cden0 F I kpl kslow kfast t0 nxcat nx rtf rpaden0 rmax vl lsep  L sa;

nxsep = nx -nxcat;

if exist('c','var')
    ces = c(1:nxsep);
    cec = c(nxsep+1:nx);
    [count,usig] = fngrowth(vl,lsep,time,nxsep);
    rsa = rpa';
else
    ces = cden0(1:nxsep);
    cec = cden0(nxsep+1:nx);
    [count,usig]= fngrowth(vl,lsep,0,nxsep);
rsa =  rpaden0';
time = 0;
end



aa = (0.5*pi*rsa.^2/rmax^3); 
% theta = exp(-0.02009*time);
theta = 0.5-0.5*tanh((time-1050)/100);
kan = theta*kslow + (1-theta)*kfast;


% ioa = F*ka*(cec.^0.5).*usig;
ioa = F*kan*(cec.^0.5).*usig;
ioc = F*kpl*sqrt(ces(1));
epa = 1-(pi/6)*(rsa/rmax).^3;         % porosity of anode (function of radius of particle)

ep  = [epa' ones(1,nxcat)];
aac = [aa' zeros(1,nxcat)];
ce = [ces cec];

kaps = 1e2*K2(1e-3*ce);
%%

%% Variable partition
v_cell = x(1,1);
p2 = x(2:nx+1,1);
p11 = x(nx+2:nx+nxsep+1,1);
jn_a = x(nx+nxsep+2:nx+2*nxsep+1);

jn = [jn_a' zeros(1,nxcat)];

%%

%% Equations

lncex = dss010(0,L,nx,log(ce));
lncexx = dss010(0,L,nx,lncex);

p2x = [0 zeros(1,nx-1)];
p2xx = dss044(0,L,nx,p2,p2x,2,1);

h1 = p2xx + F*aac./(kaps.*ep.^1.5).*jn - 2*rtf*(1-t0)*lncexx;
h1(nx) = 2*ioc*sinh(0.5/(rtf)*(v_cell-p2(nx)));

p11x = [-I/sa zeros(1,nxsep-1)];
p11xx = dss044(0,lsep,nxsep,p11,p11x,1,2);

h2 = p11xx - F*aa'/(sa*epa'.^1.5).*jn_a';
h2(1) = p11(1);

cjn1 =2*ioa'.*sinh(0.5/(rtf)*(p11-p2(1:nxsep)));
h3 = jn_a' - cjn1';

ijn1 = F*trapz(linspace(0,lsep,nxsep),cjn1'.*aa');
% ijn1 = F*aa(end)*trapz(linspace(0,lsep,nxsep),cjn1);
trg = I;
h4 =	ijn1 - trg;


h =[h1 h2 h3 h4]';
end
%------------------------------end core()---------------------------------- 	