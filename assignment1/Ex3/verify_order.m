close all

%% Butcher tableau
f = @(x) [  x(1)+x(2)+x(3)-1;
            1/4*x(2)+x(3)-1/2;
            1/16*x(2)+x(3)-1/3;
            1/4*x(3)*x(6)-1/6;
            x(4)-1/4;
            x(5)+x(6)-1;
         ];
x = fsolve(f,zeros(6,1));
b1 = x(1); b2 = x(2); b3 = x(3); a21 = x(4); a31 = x(5); a32 = x(6);

butcher.AT = [0 0 0; a21 0 0; a31 a32 0]';
butcher.b  = [b1 b2 b3]';
butcher.c  = [0 1/3 1]';
butcher.d  = [b1 b2 b3]' - [1/8 1/2 3/8]'; % b - bhat
butcher.bhat = [1/8 1/2 3/8]';
butcher.stages = 3;


%% TEST EQUATION 
addpath('../Ex1');
tspan = [0; 10];
x0 = 1;
lambda = -1;
hs = 10; % how many stepsizes to run
h = logspace(-1,-3,hs);
Lerr = zeros(hs,1);
LEerr = zeros(hs,1);

%% Plotting
for i=1:hs
    % Convert to number of steps
    n = ceil((tspan(2) - tspan(1))/h(i) + 1);
    % Analytical Solution
    T = linspace(tspan(1),tspan(2),n);
    X = exp(lambda*T);
    % ERK designed
    [T1,X1,Err1] = ExplicitRungeKutta(@TestEquation,tspan,x0,n,butcher,lambda);
    % Local error
    Lerr(i) = abs(X(2) - X1(2));
    LEerr(i) = abs(Err1(2));
end
figure
loglog(h,Lerr,'-b',h,h.^(3+1),'-.b',h,LEerr,'-r',h,h.^(2+1),'-.r')
legend('Designed Runge-Kutta', 'O(h^3) help line', 'Designed Runge-Kutta (error)', 'O(h^2) help line')
