function [Ts,Ys,info] = ClassicalRungeKutta_AdaptiveStep( ... 
                   fun,tspan,n,y0,abstol,reltol,type,varargin)

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

% Compute initial step size
hn = (tspan(2)-tspan(1))/(n-1);

% Number of stages
s = 4; 

% Allocate memory for the stages
m = size(y0,1);
T = zeros(s,1);
Y = zeros(m,s);

% Classical Runge-Kutta
Ys(1,:) = y0;
Ts(1) = tspan(1);

y = y0;
t = tspan(1);
hvec = [];
rvec = [];

funeval = 0;
while Ts(end) < tspan(2)
    
    % Size of last step
    if Ts(end)+hn > tspan(2)
        hn = tspan(2) - Ts(end);
    end
    
    rp = epstol;
    % Stage 1
    T(1) = t;
    Y(:,1) = y;
    F(:,1) = feval(fun,T(1),Y(:,1),varargin{:});
    funeval = funeval+1;
    
    AcceptStep = false;
    while ~ AcceptStep
        h = hn;
         
        % STEP SIZE H
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
        
        % Update y1
        y1 = y + h*(1/6*F(:,1) + 1/3*F(:,2) + 1/3*F(:,3) + 1/6*F(:,4));
        
        % STEP SIZE H/2
        hm = 0.5*h;
        % Stage 1
        Tm(1) = t;
        Ym(:,1) = y;
        Fm(:,1) = F(:,1);
        % Stage 2
        Tm(2) = t + 0.5*hm;
        Ym(:,2) = y + 0.5*hm*Fm(:,1);
        Fm(:,2) = feval(fun,Tm(2),Ym(:,2),varargin{:});
        % Stage 3
        Tm(3) = t + 0.5*hm;
        Ym(:,3) = y + 0.5*h*Fm(:,2);
        Fm(:,3) = feval(fun,Tm(3),Ym(:,3),varargin{:});
        % Stage 4
        Tm(4) = t + hm;
        Ym(:,4) = y + hm*Fm(:,3);
        Fm(:,4) = feval(fun,Tm(4),Ym(:,4),varargin{:});
      
        % Update ym and tm
        tm = t + hm;
        ym = y + hm*(1/6*Fm(:,1) + 1/3*Fm(:,2) + 1/3*Fm(:,3) + 1/6*Fm(:,4));
        
        % Stage 1
        Tm(1) = tm;
        Ym(:,1) = ym;
        Fm(:,1) = feval(fun,Tm(1),Ym(:,1),varargin{:});
        % Stage 2
        Tm(2) = tm + 0.5*hm;
        Ym(:,2) = ym + 0.5*hm*Fm(:,1);
        Fm(:,2) = feval(fun,Tm(2),Ym(:,2),varargin{:});
        % Stage 3
        Tm(3) = tm + 0.5*hm;
        Ym(:,3) = ym + 0.5*hm*Fm(:,2);
        Fm(:,3) = feval(fun,Tm(3),Ym(:,3),varargin{:});
        % Stage 4
        Tm(4) = tm + hm;
        Ym(:,4) = ym + hm*Fm(:,3);
        Fm(:,4) = feval(fun,Tm(4),Ym(:,4),varargin{:});
        
        % Update yhat
        yhat = ym + hm*(1/6*Fm(:,1) + 1/3*Fm(:,2) + 1/3*Fm(:,3) + 1/6*Fm(:,4));
        
        % Error estimation
        e = y1 - yhat;
        r = max(abs(e)./max(abstol,abs(yhat).*reltol));
        
        % Check condition
        AcceptStep = r <=1;

        % step size controller (Asymptotic or second order PI)
        hn = max(facmin,min((epstol/r)^ki*(rp/r)^kp,facmax))*h;
        rp = r;
        hvec(end+1) = h;
        rvec(end+1) = r;
        funeval = funeval + 10;
    end
    % Store solution
    t = t + h;
    y = yhat;
    Ts(end+1) = t;
    Ys(end+1,:) = y;
    
    % Since there is no embeded method the error must be calculated
    % by step doubling
    info.funeval = funeval;
    info.naccept = length(Ts);
    info.nreject = length(hvec) - length(Ts) + 1;
    info.hvec = hvec;
    info.rvec = rvec;

end

end