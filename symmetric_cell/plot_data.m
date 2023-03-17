tplot2 = [1 25 50 75 100];
tplot3 = [30 50 70 80 100];
lcor = ['r','g','b','m','k'];

% v_all(ptout*(ik-1)+1:ptout*(ik)) = y_den(:,1);
v_all = y_den(:,1);

p2_den = y_den(:,2:nx+1);
p11_den = y_den(:,nx+2:nx+nxsep+1); 
jn_a_den = y_den(:,nx+nxsep+2:nx+2*nxsep+1);
ce_den   = y_den(:,nx+2*nxsep+2:2*nx+2*nxsep+1);
rsa_den  = y_den(:,2*nx+2*nxsep+2:2*nx+3*nxsep+1);

epa_den  = 1-(pi/6)*(rsa_den/rmax).^3; 
aa_den  = (0.5*pi*rsa_den.^2/rmax^3);

theta = zeros(ptout,1);
for i=1:ptout
    [count, usig] = fngrowth(vl,lsep,time_den(i),nxsep);
for j =1:nxsep
    theta(i) = aa_den(i,j)*usig(j) + theta(i);
end
end


i2_a_den=zeros(ptout,nxsep); 
i2_a_den(:,1) = I;
for i=1:ptout
for j=2:nxsep
    jntemp = 0.5*jn_a_den(i,j-1) + 0.5*jn_a_den(i,j);
    i2_a_den(i,j) = aa_den(i,j)*F*jntemp*dxcat + i2_a_den(i,j-1);
end
end

figure(3)
hold on;
plot(tpla+time_den+tden*hcyc,v_all,'-');
title('Cell Voltage (V)');


figure(1)
hold on;
for i =1:5
    plot(x,p2_den(tplot2(i),:),'o','Color',lcor(i));
end   
title('Potential in electrolyte (V)');


figure(2)
hold on;
for i =1:5
    plot(x,ce_den(tplot2(i),:),'o','Color',lcor(i));
end
title('Concentration (mol/m^3)');


figure(4)
hold on;
for i =1:5
    plot(x_a,jn_a_den(tplot2(i),:),'o-','Color',lcor(i));
end
xlim([0 100e-6]);
title('Pore wall flux in moss j_n (mol/m^2/s)');


figure(5)
hold on;
for i =1:5
    plot(x_a,epa_den(tplot2(i),:),'o-','Color',lcor(i));
end
xlim([0 50e-6]);
title('Porosity');

figure(6)
hold on;
for i =1:5
    plot(x_a,rsa_den(tplot2(i),:),'o-','Color',lcor(i));
end
xlim([0 100e-6]);
title('radii');

figure(7)
hold on;
for i =1:5
    plot(x_a,aa_den(tplot2(i),:),'o-','Color',lcor(i));
end
xlim([0 100e-6]);
title('Specific Surface Area (m^2/m^3)');


% figure(7)
% hold on;
% for i =1:5
%     plot(x_a,i2_a_den(tplot2(i),:),'o-','Color',lcor(i));
% end
% xlim([0 lsep]);
% title('Ionic current density (A/m^2)');

% figure(8)
% hold on;
% for i =1:5
%     plot(x_a,p11_den(tplot2(i),:),'o-','Color',lcor(i));
% end

