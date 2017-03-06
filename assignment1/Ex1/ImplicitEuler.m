function [T,Y] = ImplicitEuler(funJac,tspan,n,y0,varargin)

m = length(y0);

% Compute step size
h = (tspan(1)-tspan(2))/n;

% Allocate memory for the solution
Y = zeros(n,m);
T = zeros(n,1);

% Set tolerance and max number of iterations for the solver
tol = 1.0e-8;
maxit = 100;

% Implicit Euler's method
Y(1,:) = y0;
T(1) = tspan(1);
for i=1:n
    
    T(i+1) = T(i)+h;
 
    % Compute explicit Euler and use it as initial value
    yinit = Y(i,:) + h*feval(funJac,T(i,:),Y(i,:),varargin{:});
    
    % Solve implicit equation
    Y(i+1,:) = NewtonsMethod(ResidualFunJac, ...
               yinit,tol,maxit,varargin{:});
end

end

function Residual