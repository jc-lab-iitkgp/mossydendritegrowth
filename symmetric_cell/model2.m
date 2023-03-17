% 	model() 	
%
function f = model2(time,y,yp)
% Provides the residuals of the discretized model equations for their
% integration by the ode15i solver

global D F I lcat L nx nxcat omega rmax t0 vl ;

nxsep = nx -nxcat;
[count,usig] = fngrowth(vl,lcat,time,nxcat);

%%
%% Variable partition
v_cell = y(1,1);
p2 = y(2:nx+1,1);
p11 = y(nx+2:nx+nxsep+1,1); 
jn_a = y(nx+nxsep+2:nx+2*nxsep+1,1);
ce   = y(nx+2*nxsep+2:2*nx+2*nxsep+1,1);
rsa  = y(2*nx+2*nxsep+2:2*nx+3*nxsep+1,1);

cet   = yp(nx+2*nxsep+2:2*nx+2*nxsep+1,1);
rsat  = yp(2*nx+2*nxsep+2:2*nx+3*nxsep+1);


epa  = 1-(pi/6)*(rsa/rmax).^3;
epat = -(0.5*pi/rmax^3)*rsa.^2.*rsat;

aa = (pi/2*rsa.^2/rmax^3).*usig';                                % Surface area cathode

a   = [aa' zeros(1,nxcat)];

ep  = [epa' ones(1,nxcat)]; 
ept = [epat'.*usig zeros(1,nxcat)];

jn  = [jn_a' zeros(1,nxcat)];
Deff = D*ep.^1.5;

%%
f0 = core2([v_cell;p2;p11;jn_a],ce',rsa',time);

cex = dss010(0,L,nxsep,ce);
cex(1) = 0;
cex(nx) = -I*(1-t0)/(F*D);
cexx = dss044(0,L,nx,ce,cex,2,2);

f1 = cet' - (Deff.*cexx + a*(1-t0).*jn - ce'.*ept)./ep;

f2 = rsat' + omega*jn_a'.*usig;

f = [f0' f1 f2]';
end
% 	end model() 	