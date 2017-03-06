function x = NewtonsMethod(FunJac,tk,xk,dt,x0,tol,maxit,varargin)
k = 0;
x = x0;
t = tk + dt;
[f,J] = feval(FunJac,t,x,varargin{:});
R = x - dt*f - xk;
I = eye(size(xk));
while ((k < maxit) && (norm(R,'inf') > tol))
    k = k+1;
    dRdx = I - dt*J;
    dx = dRdx\R;
    x = x - dx;
    [f,J] = feval(FunJac,t,x,varargin{:});
    R = x - dt*f - xk;
end