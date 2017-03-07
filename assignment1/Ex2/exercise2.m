%% Exercise 2  The Van der Pol system
addpath('../')
addpath('../Ex1')
tspan = [0; 50];
n = 10000;
y0 = [2; 0];
mu = 10;

abstol = 10e-6;
reltol = 10e-8;

%% Explicit Euler
%[T,Y] = ExplicitEuler(@VanderPolfunjac,tspan,n,y0,mu);

%% Implicit Euler
[T1,Y1] = ImplicitEuler(@VanderPolfunjac,tspan,n,y0,mu);

%% Trapezoidal
[T3,Y3] = Trapezoidal(@VanderPolfunjac,tspan,n,y0,mu);

plot(T1,Y1(:,1))
hold on
plot(T3,Y3(:,1))
