clear all
close all

%% HARMONIC OSCILLATOR 


tspan = [0; 10];
h = 0.5;
n = ceil((tspan(2) - tspan(1))/h + 1);
x0 = [1; 0];

%% Analytical Solution
T = linspace(tspan(1),tspan(2),n);
X = cos(T);

%% Explicit Euler
[T1,X1] = ExplicitEuler(@HarmonicOscillator,tspan,n,x0);

%% Implicit Euler
[T2,X2] = ImplicitEuler(@HarmonicOscillator,tspan,n,x0);

%% Trapezoidal
[T3,X3] = Trapezoidal(@HarmonicOscillator,tspan,n,x0);

%% Classical Runge-Kutta
[T4,X4] = ClassicalRungeKutta(@HarmonicOscillator,tspan,n,x0);

%% DOPRI54
% Call ERKSolverErrorEstimationParameters.m to get Butcher's Tableau

butcher = ERKSolverErrorEstimationParameters('DOPRI54');
[T5,X5,Err5] = ExplicitRungeKutta(@HarmonicOscillator,tspan,x0,n,butcher);

% DIFFERENCE BETWEEN THE TWO DOPRI54
butcher.b = [5179/57600; 0; 7571/16695; 393/640; -92097/339200; 187/2100; 1/40];
[T6,X6,Err6] = ExplicitRungeKutta(@HarmonicOscillator,tspan,x0,n,butcher);


figure
plot(T,X,'--k');
hold on
plot(T1,X1(:,1),'LineWidth',1.5);
plot(T2,X2(:,1),'LineWidth',1.5);

legend('Analytical solution','Explicit Euler','Implicit Euler','location','southwest')

%% ERROR ESTIMATION

h = logspace(-1,-2,10);
Lerr = zeros(6,10);
LerrHat = zeros(6,10);
Gerr = zeros(6,10);

for i=1:10
    tspan = [0; 10];
    n = ceil((tspan(2) - tspan(1))/h(i) + 1);
    n2 = ceil((tspan(2) - tspan(1))/h(i) + 1)*2;
    x0 = [1; 0];

    % Analytical Solution
    T = linspace(tspan(1),tspan(2),n);
    X = zeros(length(T),2);
    X(:,1) = cos(T)';
    X(:,2) = -sin(T)';
    
    % Explicit Euler
    [T1,X1] = ExplicitEuler(@HarmonicOscillator,tspan,n,x0);
    
    [T1Double,X1Double] = ExplicitEuler(@HarmonicOscillator,tspan,n2,x0);
    
    Lerr(1,i) = abs(X(2) - X1(2,1));
    LerrHat(1,i) = 2*abs(X1Double(3,1) - X1(2,1));
    Gerr(1,i) = abs(X(T == 10) - X1(abs(T1 - 10) <= 1e-3,1));
    
    % Implicit Euler
    [T2,X2] = ImplicitEuler(@HarmonicOscillator,tspan,n,x0);
    [T2Double,X2Double] = ImplicitEuler(@HarmonicOscillator,tspan,n2,x0);
    
    Lerr(2,i) = abs(X(2) - X2(2,1));
    LerrHat(2,i) = 2*abs(X2Double(3,1) - X2(2,1));
    Gerr(2,i) = abs(X(T == 10) - X2(abs(T2 - 10) <= 1e-3,1));
    
    % Trapezoidal
    [T3,X3] = Trapezoidal(@HarmonicOscillator,tspan,n,x0);
    [T3Double,X3Double] = Trapezoidal(@HarmonicOscillator,tspan,n2,x0);
    
    Lerr(3,i) = abs(X(2,1) - X3(2,1));
    LerrHat(3,i) = 4/3*abs(X3Double(3,1) - X3(2,1));
    Gerr(3,i) = abs(X(T == 10) - X3(abs(T3 - 10) <= 1e-3,1));
    
    % Classical Runge-Kutta
    [T4,X4] = ClassicalRungeKutta(@HarmonicOscillator,tspan,n,x0);
    [T4Double,X4Double] = ClassicalRungeKutta( ...
                          @HarmonicOscillator,tspan,n2,x0);
    
    Lerr(4,i) = abs(X(2) - X4(2,1));
    LerrHat(4,i) = 16/15*abs(X4Double(3,1) - X4(2,1));
    Gerr(4,i) = abs(X(T == 10) - X4(abs(T4 - 10) <= 1e-3,1));
    
    % DOPRI54
    % Call ERKSolverErrorEstimationParameters.m to get Butcher's Tableau

    butcher = ERKSolverErrorEstimationParameters('DOPRI54');
    [T5,X5,Err5] = ExplicitRungeKutta( ...
                   @HarmonicOscillator,tspan,x0,n,butcher);
               
    Lerr(5,i) = abs(X(2) - X5(2,1));
    LerrHat(5,i) = abs(Err5(2,1));
    Gerr(5,i) = abs(X(T == 10) - X5(abs(T5 - 10) <= 1e-3,1));
    
    % DOPRI54 order 4
    butcher.b = [5179/57600; 0; 7571/16695; 393/640;
                -92097/339200; 187/2100; 1/40];
    [T6,X6,Err6] = ExplicitRungeKutta( ... 
                   @HarmonicOscillator,tspan,x0,n,butcher);
    
    Lerr(6,i) = abs(X(2) - X6(2,1));
    LerrHat(6,i) = abs(Err6(2,1));
    Gerr(6,i) = abs(X(T == 10) - X6(abs(T6 - 10) <= 1e-3,1));
    
end
figure


% local error
subplot(1,2,1)
loglog(h,Lerr(1:3,:),'LineWidth',1.5)
hold on
ax = gca;
ax.ColorOrderIndex = 2;
loglog(h,h.^2,'--',h,h.^3,'--')
r = legend('Explicit Euler','Implicit Euler','Trapezoidal','$O(h^2)$','$O(h^3)$','location','southeast');
set(r,'Interpreter','latex')
subplot(1,2,2)
loglog(h,Lerr(4:6,:),'LineWidth',1.5)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
loglog(h,h.^5,'--',h,h.^6,'--')
%axis([10e-3 10e-2 10e-16 10e-4 ])
r = legend('Runge Kutta 4','DOPRI54 (5)','DOPRI54 (4)','$O(h^5)$','$O(h^6)$','location','northwest');
set(r,'Interpreter','latex')

print('LocalErrorTest2','-dpng')

% Global error t = 10
figure
loglog(h,Gerr(1:3,:),'LineWidth',1.5)
hold on
loglog(h,Gerr(4:6,:),'LineWidth',1.5)
legend('Explicit Euler','Implicit Euler','Trapezoidal','Classic Runge Kutta','DOPRI54 (5)','DOPRI54 (4)','location','southeast')

print('GlobalErrorTest2','-dpng')
% local error estimates
figure
subplot(1,2,1)
loglog(h,LerrHat(1:3,:),'LineWidth',1.5)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
loglog(h,Lerr(1:3,:),'--')
r = legend('Explicit Euler','Implicit Euler','Trapezoidal','True Explicit', ...
       'True Implicit','True Trapezoidal','location','southeast')
set(r,'fontsize',7)
subplot(1,2,2)
loglog(h,LerrHat(4:6,:),'LineWidth',1.5)
hold on
ax = gca;
ax.ColorOrderIndex = 1;
loglog(h,Lerr(4:6,:),'--')
r = legend('Runge Kutta 4','DOPRI54 (5)','DOPRI54 (4)', ...
       'Runge Kutta true','DOPRI54 (5) true', ...
       'DOPRI54 (4) true','location','southeast')
set(r,'fontsize',7)
print('LocalEstimateTest2','-dpng')