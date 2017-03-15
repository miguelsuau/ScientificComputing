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
hs = [0.1 0.01 0.001];
columns = 2;
figure
format long
for k=0:length(hs)-1
    rows = length(hs);
    n = ((tspan(2)-tspan(1))/hs(k+1))+1;
    [T1,X1,Err1] = ExplicitRungeKutta(@TestEquation,tspan,x0,n,butcher,lambda);
    % Analytical solution
    T = linspace(tspan(1),tspan(2),n)';
    X = exp(lambda*T);

    subplot(rows,columns,k*columns+1)
    plot(T1,X1(:),'-g','LineWidth',1.3);
    hold on
    plot(T,X,'--r','LineWidth',0.5);
    legend('Designed Runge-Kutta', 'Analytical solution')
    xlabel('time')
    title(sprintf('Test equation ($h = 10^{%0.0f}$, $x_0 = %d$, $\\lambda = %d$)', log10(hs(k+1)), x0, lambda), 'Interpreter', 'latex')
    
    subplot(rows,columns,k*columns+2)
    err = abs(X - X1);
    maxerr = max(err);
    plot(err)
    hold on
    xlim([0 length(err)-1]); title(sprintf('Error_{MAX} = %0.2e', maxerr));
    refline([0 maxerr]);
end

