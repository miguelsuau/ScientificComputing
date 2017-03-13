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
addpath('..'); % if ran from Ex3 directory 
tspan = [0; 50];
x0 = [2; 0];
mu = 3;

% Comparison
opts = odeset('Stats','on','Jacobian',@VanderPolJac);
sol = ode15s(@VanderPolFun,tspan,x0,opts,mu); % get stats from ode15s
T2 = sol.x'; X2 = sol.y';

% Plot
ns = [175 180 250 500 length(T2)];
hs = [0.1 0.01 0.001];
fig = figure;
format long

for k=0:length(hs)-1
    columns = 3;
    rows = length(hs);
    n = ((tspan(2)-tspan(1))/hs(k+1))+1;
    % ERK designed
    [T1,X1,Err1] = ExplicitRungeKutta(@VanderPolfunjac,tspan,x0,n,butcher,mu);
    err1 = abs(Err1(:,1));
    err2 = abs(Err1(:,2));

    subplot(rows,columns,k*columns+1)
    plot(T1,X1(:,1),'-g','LineWidth',1.3);
    xlim(tspan'); ylim([-5 5]); xlabel('time'); ylabel('x_1');
    hold on
    plot(T2,X2(:,1),'--r','LineWidth',0.5);
    legend('Designed Runge-Kutta', 'ode15s')
    title(...
        sprintf('Step size = %0.1e, \\mu = %d', hs(k+1), mu)...
    )
    
    subplot(rows,columns,k*columns+2)
    plot(T1,X1(:,2),'-g','LineWidth',1.3);
    xlim(tspan'); ylim([-5 5]); xlabel('time'); ylabel('x_2');
    hold on
    plot(T2,X2(:,2),'--r','LineWidth',0.5);
    legend('Designed Runge-Kutta', 'ode15s')
    
    % errors
    merr1 = max(err1); merr2 = max(err2);
    fprintf('\nStep size = %0.1e\nMax error = (%0.3e, %0.3e)\n', hs(k+1), merr1, merr2);
    subplot(rows,columns,k*columns+3)
    plot(err1, '-m'); 
    hold on
    xlim([0 length(err1)-1]); title(sprintf('Error_{MAX} = (%0.2e, %0.2e)', merr1, merr2));
    hline1 = refline([0 merr1]); hline1.Color = 'm'; 
    hline2 = refline([0 merr2]); hline2.Color = 'b'; 
    plot(err2, '-b');
end

set(fig, 'Position', get(0, 'ScreenSize'));

%print(fig, 'ex3_vaderpol_ode15s_fp','-dpdf','-fillpage')
