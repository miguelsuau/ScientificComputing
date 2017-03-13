%% Exercise 4 Step Size Controller
close all
addpath('../')
tspan = [0; 50];
n = 1000;
y0 = [2; 0];
mu = 3;

abstol = 1e-3;
reltol = 1e-3;

%% Explicit Euler Adaptive Step
[T1,Y1,info1] = ExplicitEuler_AdaptiveStep(...
          @VanderPolfunjac,tspan,n,y0,abstol,reltol,'PI',mu);
AdaptiveStepPlot(T1,Y1,info1)

%% Implicit Euler Adaptive Step
[T2,Y2,info2] = ImplicitEuler_AdaptiveStep(...
          @VanderPolfunjac,tspan,n,y0,abstol,reltol,'asymptotic',mu);
AdaptiveStepPlot(T2,Y2,info2)

%% Trapezoidal Adaptive Step
[T3,Y3,info3] = Trapezoidal_AdaptiveStep(...
          @VanderPolfunjac,tspan,n,y0,abstol,reltol,'asymptotic',mu);
AdaptiveStepPlot(T3,Y3,info3)

%% Classical Runge-Kutta      
[T4,Y4,info4] = ClassicalRungeKutta_AdaptiveStep( ... 
          @VanderPolfunjac,tspan,n,y0,abstol,reltol,'asymptotic',mu);
AdaptiveStepPlot(T4,Y4,info4)

%% Classical Runge-Kutta      
butcher = ERKSolverErrorEstimationParameters('DOPRI54');
[T5,Y5,Err,info5] = ExplicitRungeKutta_AdaptiveStep( ... 
          @VanderPolfunjac,tspan,n,y0,abstol,reltol,'asymptotic',butcher,mu);
AdaptiveStepPlot(T5,Y5,info5)

%% Adaptive step size plot
function AdaptiveStepPlot(T,Y,info)

figure
subplot(2,2,1)
plot(T,Y(:,1))
subplot(2,2,3)
plot(T,Y(:,2))
subplot(2,2,2)
plot(info.rvec)
axis([0 length(info.hvec) 0.6 1])
subplot(2,2,4)
plot(info.hvec)
axis([0 length(info.hvec) 0 0.7])

end