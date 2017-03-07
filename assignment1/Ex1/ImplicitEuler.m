function [T,Y] = ImplicitEuler(funJac,tspan,n,y0,varargin)

% Compute step size
h = (tspan(2)-tspan(1))/(n-1);

% Set tolerance and max number of iterations for the solver
tol = 1.0e-8;
maxit = 100;

% Implicit Euler's method
Y = y0;
T = tspan(1);
for i=1:(n-1)
    
    T(i+1) = T(i)+h;
    % Compute explicit Euler and use it as initial value
    yinit = Y(:,i) + h*feval(funJac,T(i),Y(:,i),varargin{:});
    
    % Solve implicit equation
    Y(:,i+1) = NewtonsImplicitEuler(funJac, ...
               T(i),Y(:,i),h,yinit,tol,maxit,varargin{:});
end
Y = Y';
end
