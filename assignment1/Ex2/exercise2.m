%% Exercise 2  The Van der Pol system
addpath('../')
addpath('../Ex1')
tspan = [0; 50];
n = 10000;
y0 = [2; 0];
mu = 10;

%% Explicit Euler
%[T,Y] = ExplicitEuler(@VanderPolfunjac,tspan,n,y0,mu);

%% Implicit Euler
[T1,Y1] = ImplicitEuler(@VanderPolfunjac,tspan,n,y0,mu);

%% Trapezoidal
[T3,Y3] = Trapezoidal(@VanderPolfunjac,tspan,n,y0,mu);

%% Classical Runge-Kutta
[T4,Y4] = ClassicalRungeKutta(@VanderPolfunjac,tspan,n,y0,mu);

plot(T3,Y3(:,2))
hold on
plot(T4,Y4(:,2))

%% DOPRI54
% Call ERKSolverErrorEstimationParameters.m to get Butcher's Tableau

butcher = ERKSolverErrorEstimationParameters('DOPRI54');
[T5,Y5,Err5] = ExplicitRungeKutta(@VanderPolfunjac,tspan,y0,n,butcher,mu);

% DIFFERENCE BETWEEN THE TWO DOPRI54
butcher.b = [5179/57600; 0; 7571/16695; 393/640; -92097/339200; 187/2100; 1/40];
[T6,Y6,Err6] = ExplicitRungeKutta(@VanderPolfunjac,tspan,y0,n,butcher,mu);

% Explicit methods won't work for mu large. That is, for very stiff
% problems GIVE REASON AND EXPLAIN WHAT STIFF MEANS.

% Look up embeded methods error estimation
