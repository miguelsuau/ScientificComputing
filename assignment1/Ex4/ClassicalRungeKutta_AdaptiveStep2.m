function [Ts,Ys,Err,info] = ClassicalRungeKutta_AdaptiveStep2( ...
                            fun,tspan,n,y0,abstol,reltol,type,butcher,varargin)

% Generic algorithm for explicit Runge-Kutta methods
switch type
    case 'asymptotic'
        ki = 1/5; kp = 0;
    case 'PI'
        ki = 0.4/5; kp = 0.3/5;
end

% Controller parameters
epstol = 0.8;
facmin = 0.1;
facmax = 5.0;

% Extract Butcher Tableau
s = butcher.stages;
AT = butcher.AT;
b = butcher.b;
c = butcher.c;

% Compute initial step size
hn = (tspan(2)-tspan(1))/(n-1);

% Allocate memory for the stages
m = size(y0,1);
T = zeros(1,s);
Y = zeros(m,s);
F = zeros(m,s);

% Explicit Runge-Kutta methods
t = tspan(1);
y = y0;
Ts(1) = t;
Ys(1,:) = y';
Err = zeros(1,m);
% Information for testing
funeval = 0;
hvec = [];
rvec = [];

while Ts(end) < tspan(2)
    
    % Size of last step
    if T(end)+hn > tspan(2)
        hn = tspan(2) - T(end);
    end
    rp = epstol;
    
    AcceptStep = false;
    while ~ AcceptStep
        h = hn;
        
        % STEP SIZE H
        % Precalculated parameters
        hAT = h*AT;
        hb = h*b;
        hc = h*c;
        
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
        y1 = y + F*hb;
        
        % STEP SIZE H/2
        hm = 0.5*h;
        % Precalculated parameters
        hmAT = hm*AT;
        hmb = hm*b;
        hmc = hm*c;
        
        % Following stages
        for j = 2:s
            T(j) = t + hmc(j);
            Y(:,j) = y + F(:,1:j-1)*hmAT(1:j-1,j);
            F(:,j) = feval(fun,T(j),Y(:,j),varargin{:});
            funeval = funeval + 1;
        end
        ym = y + F*hmb;
        tm = t + hm;
        
        % First stage
        T(1) = tm;
        Y(:,1) = ym;
        F(:,1) = feval(fun,T(1),Y(:,1),varargin{:});
        funeval = funeval + 1;
    
        % Following stages
        for j = 2:s
            T(j) = tm + hmc(j);
            Y(:,j) = ym + F(:,1:j-1)*hmAT(1:j-1,j);
            F(:,j) = feval(fun,T(j),Y(:,j),varargin{:});
            funeval = funeval + 1;
        end
        yhat = ym + F*hmb;
        
        
        % Error estimation
        e = abs(yhat - y1);
        r = max(abs(e)./max(abstol,abs(yhat).*reltol));
        
         % Check conditon
        AcceptStep = r <=1;
        
        % step size controller (Asymptotic or second order PI)
        hn = max(facmin,min((epstol/r)^ki*(rp/r)^kp,facmax))*h;
        rp = r;
        hvec(end+1) = h;
        rvec(end+1) = r;
    
    end
    % Update t and y
    t = t + h;
    y = yhat;
     
    % Store solution
    Ts(end+1) = t;
    Ys(end+1,:) = y';
    Err(end+1,:) = e';
    
end
info.funeval = funeval;
info.naccept = length(Ts);
info.nreject = length(hvec) - length(Ts) + 1;
info.hvec = hvec;
info.rvec = rvec;
end
