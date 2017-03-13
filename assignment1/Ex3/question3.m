close all
clear all

% Calculate our method's Butcher's tableau
format rat
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
butcher.stages = 3;

% Test equation
addpath('../Ex1'); % if ran from Ex3 directory 
tspan = [0; 10];
x0 = 1;
lambda = -1;

% Plot method vs 'true' solution
ns = [10 15 25 50];
figure
format long
for k=1:length(ns)
    [T1,X1,Err1] = ExplicitRungeKutta(@TestEquation,tspan,x0,ns(k),butcher,lambda);
    % Analytical solution
    T = linspace(tspan(1),tspan(2),ns(k))';
    X = exp(lambda*T);

    subplot(2,length(ns)/2,k)
    plot(T1,X1(:),'-og','LineWidth',1.3);
    hold on
    plot(T,X,'--or','LineWidth',0.5);
    legend('Designed Runge-Kutta', 'Analytical solution')
    xlabel('time')
    title(sprintf('Test equation (n = %d, x_0 = %d, \\lambda = %d)', ns(k), x0, lambda))
end

