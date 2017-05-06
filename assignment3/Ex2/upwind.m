function [ U, Utrue ] = upwind( uFn, a, dx, dt, xmin, xmax, tmax )
% uFn - used to calculate u(x,0)
% a, dx, dt - user passed, used calculate Courant number (method signitature can be changed to pass in Cr directly if needed)
% xmin, xmax - spatial domain
% tmax - run until this time

% For this problem
% a = 1/2;
% dx = 1/100;
% dt = 16/1000;
% xmin = -1;
% xmax = 1;

N = (1/a)*(1/dx);
t = 0;
steps = tmax/dt;
cr = (a*dt/dx);

x = linspace(xmin, xmax, N+1);
j = 2:N+1;

u0 = uFn(x,0);
u = u0;
unext = u0;

U = zeros(steps+1, length(x));
U(1,:) = u0;
U(2:end,1) = uFn(-1,dt:dt:tmax);
U(2:end,end) = uFn(1,dt:dt:tmax);

Utrue = U; % this will be later overwritten anyway

for n=2:steps+1
    % Keep track of time for calculating the true solution
    t = t+dt; 
    % Next u
    unext(j) = u(j) - cr*(u(j) - u(j-1));
    unext(1) = u(1) - cr*(u(1) - u(end-1));
    % Update the iterate and save
    u = unext;
    U(n,:) = unext;    
    % Calculate the true u
    Utrue(n,:) = uFn(x,t);
end


end

