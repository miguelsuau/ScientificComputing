function [Ts,Ys,Err,info] = ExplicitRungeKutta(fun,tspan,y0,n,butcher,varargin)
% Generic algorithm for explicit Runge-Kutta methods

% Extract Butcher Tableau
s = butcher.stages;
AT = butcher.AT;
b = butcher.b;
c = butcher.c;
d = butcher.d;

% Compute step size
h = (tspan(2)-tspan(1))/(n-1);

% Precalculated parameters
hAT = h*AT;
hb = h*b;
hc = h*c;
hd = h*d;

% Allocate memory for the stages
m = size(y0,1);
T = zeros(1,s);
Y = zeros(m,s);
F = zeros(m,s);

% Allocate memory for the solution
Ts = zeros(n,1);
Ys = zeros(n,m);
Err = zeros(n,m);

% Explicit Runge-Kutta methods
t = tspan(1);
y = y0;
Ts(1) = t;
Ys(1,:) = y';

% Info
funeval = 0;

for i = 1:(n-1)
    % First stage
    T(1) = t;
    Y(:,1) = y;
    F(:,1) = feval(fun,T(1),Y(:,1),varargin{:});
    funeval = funeval + 1;
    
    % Following stages
    for j = 2:s
        T(j) = t + hc(j);
        Y(:,j) = y + F(:,1:j-1)*hAT(1:j-1,j);
        F(:,j) = feval(fun,T(j),Y(:,j),varargin{:});
        funeval = funeval + 1;
    end
    
    % Update t and y and calulate the error with the embeded method
    t = t + h;
    y = y + F*hb;
    err = F*hd; 
    
    % Store solution
    Ts(i+1) = t;
    Ys(i+1,:) = y';
    Err(i+1,:) = err';
    
end

info.nFun = funeval;

end
