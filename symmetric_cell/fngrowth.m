function [cnt,y] = fngrowth(v,l,time,nxsep)
% Calculates the time evolving unit step function
% v is the velocity at which dendrite is growing
% l is the maximuim length a moss can grow 
% time (sec)
y = zeros(1,nxsep);
cnt = 0;
x = linspace(0,l,nxsep);
    for j=1:nxsep
    if x(j) < time*v + 25e-6
    y(j) = 1;
    cnt = j;
    else
        y(j) = 0;
    end
    end
end