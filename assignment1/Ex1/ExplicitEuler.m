function [T,Y] = ExplicitEuler(funJac,tspan,n,y0,varargin)

m = length(y0);

% Compute step size
h = (tspan(2)-tspan(1))/(n-1);

% Allocate memory for the solution
Y = zeros(m,n);
T = zeros(n,1);

% Explicit Euler's method
Y(:,1) = y0;
T(1) = tspan(1);

for i=1:(n-1)
    T(i+1) = T(i)+h;
    Y(:,i+1) = Y(:,i) + h*feval(funJac,T(i),Y(:,i),varargin{:});
end
Y = Y';
end