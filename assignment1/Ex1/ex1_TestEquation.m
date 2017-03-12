%% Exercise 1 

%% TEST EQUATION 

tspan = [0; 10];
n = 50;
x0 = 1;
lambda = -1;

%% Analytical solution
T = linspace(tspan(1),tspan(2),n)';
X = exp(lambda*T);

%% Explicit Euler
[T1,X1] = ExplicitEuler(@TestEquation,tspan,n,x0,lambda);

%% Implicit Euler
[T2,X2] = ImplicitEuler(@TestEquation,tspan,n,x0,lambda);

%% Trapezoidal
[T3,X3] = Trapezoidal(@TestEquation,tspan,n,x0,lambda);

%% Classical Runge-Kutta
[T4,X4] = ClassicalRungeKutta(@TestEquation,tspan,n,x0,lambda);

plot(T,X(:),'--k')
hold on
plot(T1,X1(:),'LineWidth',1.5)
plot(T2,X2(:),'LineWidth',1.5)

legend('Analytical solution','Explicit Euler','Implicit Euler')

%% DOPRI54
% Call ERKSolverErrorEstimationParameters.m to get Butcher's Tableau

butcher = ERKSolverErrorEstimationParameters('DOPRI54');
[T5,X5,Err5] = ExplicitRungeKutta(@TestEquation,tspan,x0,n,butcher,lambda);

% DIFFERENCE BETWEEN THE TWO DOPRI54
butcher.b = [5179/57600; 0; 7571/16695; 393/640; -92097/339200; 187/2100; 1/40];
[T6,X6,Err6] = ExplicitRungeKutta(@TestEquation,tspan,x0,n,butcher,lambda);

%% ERROR ESTIMATION

h = logspace(-1,-3,10);
Lerr = zeros(6,10);
LerrHat = zeros(6,10);
Gerr = zeros(6,10);

for i=1:10
    n = ceil((tspan(2) - tspan(1))/h(i) + 1);
    n2 = ceil((tspan(2) - tspan(1))/h(i)*2 + 1);
    % Analytical Solution
    T = linspace(tspan(1),tspan(2),n);
    X = exp(lambda*T);
    
    % Explicit Euler
    [T1,X1] = ExplicitEuler(@TestEquation,tspan,n,x0,lambda);
    
    [T1Double,X1Double] = ExplicitEuler(@TestEquation,tspan,n2,x0,lambda);
    
    Lerr(1,i) = abs(X(2) - X1(2));
    LerrHat(1,i) = abs(X1Double(2) - X1(2));
    Gerr(1,i) = abs(X(T == 10) - X1(abs(T1 - 10) <= 1e-3));
    
    % Implicit Euler
    [T2,X2] = ImplicitEuler(@TestEquation,tspan,n,x0,lambda);
    [T2Double,X2Double] = ImplicitEuler(@TestEquation,tspan,n2,x0,lambda);
    
    Lerr(2,i) = abs(X(2) - X2(2));
    LerrHat(2,i) = abs(X2Double(2) - X2(2));
    Gerr(2,i) = abs(X(T == 10) - X2(abs(T2 - 10) <= 1e-3));
    
    % Trapezoidal
    [T3,X3] = Trapezoidal(@TestEquation,tspan,n,x0,lambda);
    [T3Double,X3Double] = Trapezoidal(@TestEquation,tspan,n2,x0,lambda);
    
    Lerr(3,i) = abs(X(2) - X3(2));
    LerrHat(3,i) = abs(X3Double(2) - X3(2));
    Gerr(3,i) = abs(X(T == 10) - X3(abs(T3 - 10) <= 1e-3));
    
    % Classical Runge-Kutta
    [T4,X4] = ClassicalRungeKutta(@TestEquation,tspan,n,x0,lambda);
    [T4Double,X4Double] = ClassicalRungeKutta( ...
                          @TestEquation,tspan,n2,x0,lambda);
    
    Lerr(4,i) = abs(X(2) - X4(2));
    LerrHat(4,i) = abs(X4Double(2) - X4(2));
    Gerr(4,i) = abs(X(T == 10) - X4(abs(T4 - 10) <= 1e-3));
    
    % DOPRI54
    % Call ERKSolverErrorEstimationParameters.m to get Butcher's Tableau

    butcher = ERKSolverErrorEstimationParameters('DOPRI54');
    [T5,X5,Err5] = ExplicitRungeKutta( ...
                   @TestEquation,tspan,x0,n,butcher,lambda);
               
    Lerr(5,i) = abs(X(2) - X5(2));
    LerrHat(5,i) = abs(Err5(2));
    Gerr(5,i) = abs(X(T == 10) - X5(abs(T5 - 10) <= 1e-3));
    
    % DOPRI54 order 4
    butcher.b = [5179/57600; 0; 7571/16695; 393/640;
                -92097/339200; 187/2100; 1/40];
    [T6,X6,Err6] = ExplicitRungeKutta( ... 
                   @TestEquation,tspan,x0,n,butcher,lambda);
    
    Lerr(6,i) = abs(X(2) - X6(2));
    LerrHat(6,i) = abs(Err6(2));
    Gerr(6,i) = abs(X(T == 10) - X6(abs(T6 - 10) <= 1e-3));
    
end
figure

% local error
loglog(h,Lerr(1:3,:))
legend('Explicit Euler','Implicit Euler','Trapezoidal','location','southeast')
figure
loglog(h,Lerr(4:6,:))
legend('Classic Runge Kutta','DOPRI54','DOPRI542','location','southeast')

% Global error t = 10
figure
loglog(h,Gerr(1:3,:))
legend('Explicit Euler','Implicit Euler','Trapezoidal','location','southeast')
title('Global error')
figure
loglog(h,Gerr(4:6,:))
legend('Classic Runge Kutta','DOPRI54','DOPRI542','location','southeast')
title('Global error')

% local error estimates
figure
loglog(h,LerrHat(1:3,:))
legend('Explicit Euler','Implicit Euler','Trapezoidal','location','southeast')
figure
loglog(h,LerrHat(4:6,:))
legend('Classic Runge Kutta','DOPRI54','DOPRI542','location','southeast')