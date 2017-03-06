%% Exercise 4 Step Size Controller

tspan = [0; 50];
n = 100;
y0 = [2 0];
mu = 3;

abstol = 10e-6;
reltol = 10e-8;

%% Explicit Euler Adaptive Step
[T,Y] = ExplicitEuler_AdaptiveStep(...
        @VanderPolfunjac,tspan,n,y0,abstol,reltol,mu);
 