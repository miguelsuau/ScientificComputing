function [T,Y] = ImplicitEuler_AdaptiveStep(... 
                 funJac,tspan,n,y0,abstol,reltol,varargin)

% Controller parameters
epstol = 0.8;
facmin = 0.1;
facmax = 5.0;

% Set tolerance and max number of iterations for the solver
tol = 1.0e-8;
maxit = 100;

% Implicit Euler's method
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
    
    h = h0;
    AcceptStep = false;
    while ~ AcceptStep
        h = hn;    

        % Step size h

        % Compute explicit Euler and use it as initial value
        yinit = Y(:,end) + h*f;

        % Solve implicit equation
        Y1 = NewtonsImplicitEuler(funJac, ...
                   T(end),Y(:,end),h,yinit,tol,maxit,varargin{:});

        % Step size h/2
        hm = 0.5*h;
        Tm = T(end) + hm;
        yinit = Y(:,end) + hm*f;
        Ym = NewtonsImplicitEuler(funJac, ...
                   T(end),Y(:,end),hm,yinit,tol,maxit,varargin{:});
        fm = feval(funJac,Tm,Ym,varargin{:});
        yinit = Ym + hm*fm;
        Yhat = NewtonsImplicitEuler(funJac, ...
                   Tm,Ym,hm,yinit,tol,maxit,varargin{:});

        % Error estimation
        e = Y1 - Yhat;
        r = max(abs(e)./max(abstol,abs(Yhat).*reltol));

        % Check condition
        AcceptStep = r <=1;

        % Asymptotic step size controller
        hn = max(facmin,min(sqrt(epstol/r),facmax))*h;
    end
    
    T(end+1) = T(end)+h;
    Y(:,end+1) = Yhat;
    
end
Y = Y';
end
