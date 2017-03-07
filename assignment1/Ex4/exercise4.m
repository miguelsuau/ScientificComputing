%% Exercise 4 Step Size Controller
addpath('../')
tspan = [0; 50];
n = 1000;
y0 = [2; 0];
mu = 10;

abstol = 10e-6;
reltol = 10e-8;

%% Explicit Euler Adaptive Step
[T,Y] = ExplicitEuler_AdaptiveStep(...
        @VanderPolfunjac,tspan,n,y0,abstol,reltol,mu);

%% Implicit Euler Adaptive Step
[T1,Y1] = ImplicitEuler_AdaptiveStep(...
          @VanderPolfunjac,tspan,n,y0,abstol,reltol,mu);
