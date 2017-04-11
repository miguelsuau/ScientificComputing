%% Exercise 4 Step Size Controller
close all
addpath('../')
tspan = [0; 50];
n = 1;
y0 = [2; 0];
mu = 3;

abstol = 1e-3;
reltol = 1e-3;

%% Explicit Euler Adaptive Step
[T1,Y1,info1] = ExplicitEuler_AdaptiveStep(...
          @VanderPolfunjac,tspan,n,y0,abstol,reltol,'PI',mu);
%AdaptiveStepPlot(T1,Y1,info1)

%% Implicit Euler Adaptive Step
[T2,Y2,info2] = ImplicitEuler_AdaptiveStep(...
          @VanderPolfunjac,tspan,n,y0,abstol,reltol,'PI',mu);
%AdaptiveStepPlot(T2,Y2,info2)

%% Trapezoidal Adaptive Step
[T3,Y3,info3] = Trapezoidal_AdaptiveStep(...
          @VanderPolfunjac,tspan,n,y0,abstol,reltol,'PI',mu);
AdaptiveStepPlot(T3,Y3,info3)

%% Classical Runge-Kutta 

[T4,Y4,info4] = ClassicalRungeKutta_AdaptiveStep2( ... 
          @VanderPolfunjac,tspan,n,y0,abstol,reltol,'PI',butcher,mu);
AdaptiveStepPlot(T4,Y4,info4)

%% DOPRI54     
butcher = ERKSolverErrorEstimationParameters('DOPRI54');
[T5,Y5,Err,info5] = ExplicitRungeKutta_AdaptiveStep( ... 
          @VanderPolfunjac,tspan,n,y0,abstol,reltol,'PI',butcher,mu);
AdaptiveStepPlot(T5,Y5,info5)
print('AdaptiveStepDOPRI54','-dpng')

%% Adaptive step size plot
function AdaptiveStepPlot(T,Y,info)

figure
subplot(2,2,1)
plot(T,Y(:,1))
r = title('$y(t)$');
set(r,'Interpreter','latex')
axis([0 50 -2.5 2.5])
subplot(2,2,3)
plot(T,Y(:,2))
r = title('$\dot{y(t)}$');
set(r,'Interpreter','latex')
axis([0 50 -6 6])
subplot(2,2,2)
plot(linspace(0,50,length(info.rvec)),info.rvec)
axis([0 50 0 1.1])
r = title('$r(t)$');
set(r,'Interpreter','latex')
subplot(2,2,4)
plot(linspace(0,50,length(info.hvec)),info.hvec)
axis([0 50 0 0.3])
r = title('$h(t)$');
set(r,'Interpreter','latex')

end