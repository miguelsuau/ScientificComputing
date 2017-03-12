function [Tout,Xout,info,stats] = ESDIRK23(fun,jac,t0,tf,x0,h0,absTol,relTol,varargin)

% ESDIRK23 
% Modified version for Ex5 according to the tips given in the lecture.
%=========================================================================
% Runge-Kutta method parameters
gamma = 1-1/sqrt(2);
a31 = (1-gamma)/2;
AT = [0 gamma a31;0 gamma a31;0 0 gamma];
c  = [0; 2*gamma; 1];
b  = AT(:,3);
bhat = [    (6*gamma-1)/(12*gamma); ...
            1/(12*gamma*(1-2*gamma)); ...
            (1-3*gamma)/(3*(1-2*gamma))    ];
d  = b-bhat;
% p  = 2; % ex5
% phat = 3; % ex5
s = 3;


% error and convergence controller
epsilon = 0.8;
tau = 0.1*epsilon; %0.005*epsilon;
itermax = 20;
% ke0 = 1.0/phat;
% ke1 = 1.0/phat;
% ke2 = 1.0/phat;
% alpharef = 0.3;
% alphaJac = -0.2;
% alphaLU  = -0.2;
% hrmin = 0.01;
% hrmax = 10;
%========================================================================
tspan = [t0 tf]; % carsten
info = struct(...
            'nStage',    s,       ... % carsten
            'absTol',    absTol,  ... % carsten % ex5
            'relTol',    relTol,  ... % carsten % ex5
            'iterMax',   itermax, ... % carsten
            'tspan',     tspan,   ... % carsten
            'nFun',      0, ...
            'nJac',      0, ...
            'nLU',       0, ...
            'nBack',     0, ...
            'nStep',     0, ...
            'nAccept',   0, ...
            'nFail',     0, ...
            'nDiverge',  0, ...
            'nSlowConv', 0);


        
% Main ESDIRK Integrator
%========================================================================
nx = size(x0,1);
F = zeros(nx,s);
t = t0;
x = x0;
h = h0;
IG = eye(length(x0)); % ex5 replacement for g

[F(:,1),~]  = feval(fun,t,x,varargin{:}); % ex5 no need for g
info.nFun = info.nFun+1;
[dfdx,~] = feval(jac,t,x,varargin{:}); % ex5 no need for g
info.nJac = info.nJac+1;
%FreshJacobian = true; % ex5
if (t+h)>tf
    h = tf-t;
end
hgamma = h*gamma;
dRdx = IG - hgamma*dfdx;
[L,U,pivot] = lu(dRdx,'vector');
info.nLU = info.nLU+1;
%hLU = h; % ex5

%FirstStep = true; % ex5
%ConvergenceRestriction = false; % ex5
%PreviousReject = false; % ex5
iter = zeros(1,s);

% Output
chunk = 100;
Tout = zeros(chunk,1);
Xout = zeros(chunk,nx);
%Gout = zeros(chunk,nx); % ex5

Tout(1,1) = t;
Xout(1,:) = x.';
%Gout(1,:) = g.'; % ex5

while t<tf
    info.nStep = info.nStep+1;
    %=====================================================================
    % A step in the ESDIRK method
    i=1;   
    diverging = false;
    SlowConvergence = false; % carsten
    alpha = 0.0;
    Converged = true;
    while (i<s) && Converged
        % Stage i=2,...,s of the ESDIRK Method
        i=i+1;
        phi = x + F(:,1:i-1)*(h*AT(1:i-1,i)); % ex5

        % Initial guess for the state % ex5 removed duplicate code
        dt = c(i)*h;
        %G = g + dt*F(:,1);
        X = x + dt*F(:,1); % ex5
        T = t+dt;
            
        [F(:,i),~] = feval(fun,T,X,varargin{:}); % ex5
        info.nFun = info.nFun+1;
        R = X - hgamma*F(:,i) - phi; % ex5
        rNewton = norm(R./(absTol + abs(X).*relTol), inf); % ex5
        Converged = (rNewton < tau);
        
        % Newton Iterations
        while ~Converged && ~diverging && ~SlowConvergence
            iter(i) = iter(i)+1;
            dX = U\(L\(R(pivot,1)));
            info.nBack = info.nBack+1;
            X = X - dX;
            rNewtonOld = rNewton;
            [F(:,i),~] = feval(fun,T,X,varargin{:}); % ex5
            info.nFun = info.nFun+1;
            R = X - hgamma*F(:,i) - phi; % ex5
            rNewton = norm(R./(absTol + abs(X).*relTol), inf); % ex5
            alpha = max(alpha,rNewton/rNewtonOld);
            Converged = (rNewton < tau);
            diverging = (alpha >= 1);
            SlowConvergence = (iter(i) >= itermax); 
        end
        diverging = (alpha >= 1)*i; % carsten, recording which stage is diverging
    end
    %if diverging, i, iter, pause, end
    nstep = info.nStep;
    stats.t(nstep) = t;
    stats.h(nstep) = h;
    stats.r(nstep) = NaN;
    stats.iter(nstep,:) = iter;
    stats.Converged(nstep) = Converged;
    stats.Diverged(nstep)  = diverging;
    stats.AcceptStep(nstep) = false;
    stats.SlowConv(nstep)  = SlowConvergence*i; % carsten, recording which stage is converging to slow (reaching maximum no. of iterations)
    iter(:) = 0; % carsten
    %=====================================================================
    
    % Error estimation
    e = F*(h*d);
    r = norm(e./(absTol + abs(X).*relTol), inf); % ex5
    r = max(r,eps);
    stats.r(nstep) = r;
    t = T;
    x = X;
    F(:,1) = F(:,s);  

    % Jacobian Update Strategy
    [dfdx,~] = feval(jac,t,x,varargin{:}); % ex5
    info.nJac = info.nJac+1;
    hgamma = h*gamma;
    dRdx = IG - hgamma*dfdx; % ex5
    [L,U,pivot] = lu(dRdx,'vector');
    info.nLU = info.nLU+1;
    info.nFail = info.nFail + ~Converged; % ex5
    info.nDiverge = info.nDiverge + (~Converged && diverging); % ex5
    
    %=====================================================================
    % Storage of variables for output % ex5
    info.nAccept = info.nAccept + 1;
    nAccept = info.nAccept;
    if nAccept > length(Tout);
       Tout = [Tout; zeros(chunk,1)];
       Xout = [Xout; zeros(chunk,nx)];
    end
    Tout(nAccept,1) = t;
    Xout(nAccept,:) = x.';
end
info.nSlowConv = length(find(stats.SlowConv)); % carsten
Tout = Tout(1:nAccept,1);
Xout = Xout(1:nAccept,:);

