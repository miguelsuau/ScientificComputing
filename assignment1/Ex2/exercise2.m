%% Exercise 2  The Van der Pol system
addpath('../')
addpath('../Ex1')
tspan = [0; 50];
y0 = [2; 0];
n = 1000;
h = (tspan(2)-tspan(1))/(n-1);
mu = 3;

%% Explicit Euler
[T1,Y1] = ExplicitEuler(@VanderPolfunjac,tspan,n,y0,mu);

%% Implicit Euler
[T2,Y2] = ImplicitEuler(@VanderPolfunjac,tspan,n,y0,mu);

%% Trapezoidal
[T3,Y3] = Trapezoidal(@VanderPolfunjac,tspan,n,y0,mu);

%% Classical Runge-Kutta
[T4,Y4] = ClassicalRungeKutta(@VanderPolfunjac,tspan,n,y0,mu);


%% DOPRI54
% Call ERKSolverErrorEstimationParameters.m to get Butcher's Tableau

butcher = ERKSolverErrorEstimationParameters('DOPRI54');
[T5,Y5,Err5] = ExplicitRungeKutta(@VanderPolfunjac,tspan,y0,n,butcher,mu);


subplot(2,1,1)
l1 = plot(T1,Y1(:,1));
hold on
l2 = plot(T2,Y2(:,1));
l3 = plot(T3,Y3(:,1));
l4 = plot(T4,Y4(:,1));
l5 = plot(T5,Y5(:,1));
axis([0 50 -Inf Inf])
subplot(2,1,2)
l6 = plot(T1,Y1(:,2));
hold on
l7 = plot(T2,Y2(:,2));
l8 = plot(T3,Y3(:,2));
l9 = plot(T4,Y4(:,2));
l10 = plot(T5,Y5(:,2));
axis([0 50 -Inf Inf])
hL = legend([l1,l2,l3,l4,l5],{'Explicit Euler','Implicit Euler','Trapezoidal','Classic Runge Kutta','DOPRI54 (5)'});

%% MU = 100

tspan = [0; 250];
y0 = [2; 0];
n = 100000;
h = (tspan(2)-tspan(1))/(n-1);

mu = 100;
figure
%% Explicit Euler
[T1,Y1] = ExplicitEuler(@VanderPolfunjac,tspan,n,y0,mu);

%% Implicit Euler
[T2,Y2] = ImplicitEuler(@VanderPolfunjac,tspan,n,y0,mu);

%% Trapezoidal
[T3,Y3] = Trapezoidal(@VanderPolfunjac,tspan,n,y0,mu);

%% Classical Runge-Kutta
[T4,Y4] = ClassicalRungeKutta(@VanderPolfunjac,tspan,n,y0,mu);

%% DOPRI54
% Call ERKSolverErrorEstimationParameters.m to get Butcher's Tableau

butcher = ERKSolverErrorEstimationParameters('DOPRI54');
[T5,Y5,Err5] = ExplicitRungeKutta(@VanderPolfunjac,tspan,y0,n,butcher,mu);

subplot(2,1,1)
l1 = plot(T1,Y1(:,1));
hold on
l2 = plot(T2,Y2(:,1));
l3 = plot(T3,Y3(:,1));
l4 = plot(T4,Y4(:,1));
l5 = plot(T5,Y5(:,1));
axis([0 250 -Inf Inf])
subplot(2,1,2)
l6 = plot(T1,Y1(:,2));
hold on
l7 = plot(T2,Y2(:,2));
l8 = plot(T3,Y3(:,2));
l9 = plot(T4,Y4(:,2));
l10 = plot(T5,Y5(:,2));
axis([0 250 -Inf Inf])
hL = legend([l1,l2,l3,l4,l5],{'Explicit Euler','Implicit Euler','Trapezoidal','Classic Runge Kutta','DOPRI54 (5)'});

%% ERROR ESTIMATION

h = logspace(-1,-3,10);
Lerr = zeros(5,10);
tspan = [0; 250];
y0 = [2; 0];
mu = 100;
for i=1:10
    
    n = ceil((tspan(2) - tspan(1))/h(i) + 1);
    n2 = ceil((tspan(2) - tspan(1))/h(i) + 1)*2;

    % Explicit Euler
    [T1,X1] = ExplicitEuler(@VanderPolfunjac,tspan,n,y0,mu);
    
    [T1Double,X1Double] = ExplicitEuler(@VanderPolfunjac,tspan,n2,y0,mu);
    
    Lerr(1,i) = 2*abs(X1Double(3,1) - X1(2,1));
    
    % Implicit Euler
    [T2,X2] = ImplicitEuler(@VanderPolfunjac,tspan,n,y0,mu);
    [T2Double,X2Double] = ImplicitEuler(@VanderPolfunjac,tspan,n2,y0,mu);
    
    Lerr(2,i) = 2*abs(X2Double(3,1) - X2(2,1));
    
    % Trapezoidal
    [T3,X3] = Trapezoidal(@VanderPolfunjac,tspan,n,y0,mu);
    [T3Double,X3Double] = Trapezoidal(@VanderPolfunjac,tspan,n2,y0,mu);
    
    Lerr(3,i) = 4/3*abs(X3Double(3,1) - X3(2,1));
    
    % Classical Runge-Kutta
    [T4,X4] = ClassicalRungeKutta(@VanderPolfunjac,tspan,n,y0,mu);
    [T4Double,X4Double] = ClassicalRungeKutta( ...
                          @VanderPolfunjac,tspan,n2,y0,mu);
    
    Lerr(4,i) = 16/15*abs(X4Double(3,1) - X4(2,1));
    
    % DOPRI54
    % Call ERKSolverErrorEstimationParameters.m to get Butcher's Tableau

    butcher = ERKSolverErrorEstimationParameters('DOPRI54');
    [T5,X5,Err5] = ExplicitRungeKutta( ...
                   @VanderPolfunjac,tspan,y0,n,butcher,mu);
               
    Lerr(5,i) = abs(Err5(2,1));
    
end

% local error estimates
figure
loglog(h,Lerr,'LineWidth',1.5)
r = legend('Explicit Euler','Implicit Euler','Trapezoidal','Runge Kutta 4','DOPRI54 (5)')
print('VanderPolError','-dpng')