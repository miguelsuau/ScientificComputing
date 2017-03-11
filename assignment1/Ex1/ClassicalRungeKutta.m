function [Ts,Ys] = ClassicalRungeKutta(fun,tspan,n,y0,varargin)
% Compute step size
h = (tspan(2)-tspan(1))/(n-1);

% Number of stages
s = 4; 

% Allocate memory for the solution
m = size(y0,1);
Ts = zeros(n,1);
Ys = zeros(n,m);

% Allocate memory for the stages
T = zeros(s,1);
Y = zeros(m,s);

% Classical Runge-Kutta
Ys(1,:) = y0;
Ts(1) = tspan(1);

y = y0;
t = tspan(1);
for i=1:(n-1)
    
    % Stage 1
    T(1) = t;
    Y(:,1) = y;
    F(:,1) = feval(fun,T(1),Y(:,1),varargin{:});
    % Stage 2
    T(2) = t + 0.5*h;
    Y(:,2) = y + 0.5*h*F(:,1);
    F(:,2) = feval(fun,T(2),Y(:,2),varargin{:});
    % Stage 3
    T(3) = t + 0.5*h;
    Y(:,3) = y + 0.5*h*F(:,2);
    F(:,3) = feval(fun,T(3),Y(:,3),varargin{:});
    % Stage 4
    T(4) = t + h;
    Y(:,4) = y + h*F(:,3);
    F(:,4) = feval(fun,T(4),Y(:,4),varargin{:});
    
    % Update y and t
    t = t + h;
    y = y + h*(1/6*F(:,1) + 1/3*F(:,2) + 1/3*F(:,3) + 1/6*F(:,4));
    
    % Store solution
    Ts(i+1) = t;
    Ys(i+1,:) = y;
    
    % Since there is no embeded method the error must be calculated
    % by step doubling
    
end

end