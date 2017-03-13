function [T,Y] = Trapezoidal_AdaptiveStep(...
                 funJac,tspan,n,y0,abstol,reltol,type,varargin)

switch type
    case 'asymptotic'
        a = 1/3; b = 0; c = 0;
    case 'PI'
        a = 1/6; b = 1/6; c = 1/2;
end
% Controller parameters
epstol = 0.8;
facmin = 0.1;
facmax = 5.0;

% Set tolerance and max number of iterations for the solver
tol = 1.0e-8;
maxit = 100;

% Trapezoidal method
Y = y0;
T(1) = tspan(1);

% Compute initial step size
h0 = (tspan(2)-tspan(1))/(n-1);
hn = h0;

while T(end) < tspan(2)
    
    % Size of last step
    if T(end)+hn > tspan(2)
        hn = tspan(2) - T(end);
    end
    
    f = feval(funJac,T(end),Y(:,end),varargin{:});
    
    hp = hn;
    rp = epstol;
    
    AcceptStep = false;
    while ~ AcceptStep
        h = hn;
        
        % Step size h
        
        % Compute explicit Euler and use it as initial value
        yinit = Y(:,end) + h*f;

        % Solve implicit equation
        Y1 = NewtonsTrapezoidal(funJac, ...
                     T(end),Y(:,end),h,yinit,tol,maxit,varargin{:});
        
        % Step size h/2
        
        hm = 0.5*h;
        Tm = T(end) + hm;
        
        % Compute first half step
        yinit = Y(:,end) + hm*f;
        Ym = NewtonsTrapezoidal(funJac,...
             T(end),Y(:,end),hm,yinit,tol,maxit,varargin{:});
        
        % Compute second half step
        fm = feval(funJac,Tm,Ym,varargin{:});
        yinit = Ym + hm*fm;
        Yhat = NewtonsTrapezoidal(funJac,...
               Tm,Ym,hm,yinit,tol,maxit,varargin{:});
        
        % Error estimation
        e = Y1 - Yhat;
        r = max(abs(e)./max(abstol,abs(Yhat).*reltol));
        
        % Check conditon
        AcceptStep = r <=1;
        
        % step size controller (Asymptotic or second order PI)
        hn = max(facmin,min((epstol/r)^a*(epstol/rp)^b*(h/hp)^-c,facmax))*h;
        rp = r;
        hp = h;
    end
    
    T(end+1) = T(end)+h;    
    Y(:,end+1) = Yhat;
    
end

Y = Y';
end