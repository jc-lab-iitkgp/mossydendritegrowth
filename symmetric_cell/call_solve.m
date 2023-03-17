clearvars
close all
clc
global ce0 cden0 F D I kpl kslow kfast omega t0 nx nxcat dx L lsep lcat rtf rpaden0 rmax sa vl;

alpha    = 0.5;               % Activity coefficient for both reaction.
ce0      = 2000;              % Initial electrolyte concentration (mol/m3)
D        = 7.5e-11;	          % mˆ2/s Diffusion Coeff. Electrolyte
F        = 96487;             % Farady Constants. (A*s/mol)
I        = -10;               % Discharging Current density(+ve). (A/m^2)
kpl      = 3.5e-8;            % reaction rate constant for epitaxial plating/dissoution. (m/s)
kslow       = 8e-12;          % reaction rate constant for anode reaction. (m/s) slow (root-growth)
kfast    = 2.756e-10;         % reaction rate constant for anode reaction. (m/s) fast (tip-growth)
omega    = 1.31e-5;           % Molar volume of Li (m3/mol) 
R        = 8.3143;            % Universal gas constant. (J/mol.K)
sa       = 1000;              % Conductivity of lithium metal
t0       = 0.363;             % Transference Number
T        = 300;               % Tempearature of the entire cell (K) 
rpa0     = 0.05e-6;           %	Initial radius of spherical dendrites in Anode(m)
rmax     = 0.5e-6;
rtf      = R*T/F;

%============================ Cell Geometry ===============================
L     = 750e-6 - 4.85e-6;            % Thickness of separator region. (m)
lcat = L/2;
lsep = L-lcat;
%========================== Grid Parameters ===============================
nx = 200;
nxcat = nx/2;
nxsep = nx - nxcat;
dx = L/(nx-1);
dxcat = lcat/(nxcat-1);

x = linspace(0,L,nx);
x_a = linspace(0,lsep,nxsep);
vl       = -8*(I/F*omega);   % Velocity of dendrite growth (m/s)
% ========================== Time Domain ==================================
tpla     = 3600;                     % Simulation time plating (sec)
tden     = 1800;                     % Simulation time plating (sec)


% --------------------------- (Initialization)-----------------------------
cden0 = -I*(1-t0)/(F*D)*(x-0.5*L) + ce0;
rpaden0 = rpa0*ones(1,nxcat);

for ik=1:1
[count, usig] = fngrowth(vl,lsep,0,nxsep);

x0 = [0.3; 0.15*ones(nx,1); 0.1*ones(nxsep,1); -1e-6*ones(nxsep,1)];  
sol0 = fsolve(@core2,x0);

x0 = [sol0;cden0';rpaden0'];
xp0= [zeros(size(sol0));zeros(nx,1);zeros(nxcat,1)];

if tden > 6000
    ptout = round(tden/60);
elseif tden >100
    ptout = 100;
else
    ptout = tden;
end
tspand = linspace(0,tden,ptout);

[time_den,y_den] = ode15i(@model2,tspand,x0,xp0);
hcyc = 2*(ik-1);
plot_data;
end