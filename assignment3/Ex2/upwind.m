function [ U, Utrue ] = upwind( uFn, a, dx, dt, xmin, xmax, tmax )

%
% a = 1/2;
% dx = 1/100;
N = (1/a)*(1/dx);
% dt = 16/1000;
% xmin = -1;
% xmax = 1;
t = 0;
steps = tmax/dt;

x = xmin-dx:dx:xmax+dx;

u0 = uFn(x,0);
u = u0;
unext = u0;

idx = 2:N+2;

U = zeros(steps+1, length(x));
U(1,:) = u0;
Utrue = U;

for n=2:steps+1
    t = t+dt;
    
    % BCs (TODO: ask if correct? u doesn't seem to damp)
    u(1) = uFn(-1,t);
    u(N+3) = uFn(1,t);
    
    % Next u
    unext(idx) = u(idx) - (a*dt/dx)*(u(idx) - u(idx-1));
    u = unext;
    U(n,:) = u;
    
    % Calculate the true u
    Utrue(n,:) = uFn(x,t);
end


end

