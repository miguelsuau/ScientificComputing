function [T,Y] = ExplicitEuler_AdaptiveStep(...
                 fun,tspan,n,y0,abstol,reltol,varargin)
% Controller parameters
epstol = 0.8;
facmin = 0.1;
facmax = 5.0; 

% Explicit Euler's method
Y = y0;
T(1) = tspan(1);

% Compute initial step size
hn = (tspan(2)-tspan(1))/(n-1);

while T(end) < tspan(2)
    % hn = h0;
    if T(end)+hn > tspan(2)
        hn = tspan(2)-T(end);
    end
    
    f = feval(fun,T(end),Y(:,end),varargin{:});
       
    AcceptStep = false;
    while ~ AcceptStep
        h = hn;
        % Step size h
        Y1 = Y(:,end) + h*f;
    
        % Step size h/2
        hm = 0.5*h;
        Tm = T(end) + hm;
        Ym = Y(:,end) + hm*f;
        fm = feval(fun,Tm,Ym,varargin{:});
        Yhat = Ym + hm*fm;
    
        % Error estimation
        e = Y1 - Yhat;
        r = max(abs(e)./max(abstol,abs(Yhat).*reltol));
        
        % Check condition
        AcceptStep = r <= 1;
        
        % Asymptotic step size controller
        hn = max(facmin,min(sqrt(epstol/r),facmax))*h;
    end
    
    T(end+1) = T(end) + h;
    Y(:,end+1) = Yhat; 

end
Y = Y';
end