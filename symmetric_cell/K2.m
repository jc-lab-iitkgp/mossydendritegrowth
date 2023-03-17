function g = K2(x)
% Calculates the electric conductivity (kappa) of the electrolyte
c0 =	1.0793e-4;
c1 =	6.7461e-3;
c2 =    -5.2245e-3;
c3 =	1.3605e-3;
c4 =    -1.1724e-4;

g = c0 + c1*x + c2*x.^2 + c3*x.^3 + c4*x.^4;
end