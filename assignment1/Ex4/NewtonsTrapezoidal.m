function [y,funeval] = NewtonsTrapezoidal(FunJac,tk,yk,h,y0,tol,maxit,varargin)

k = 0;
t = tk + h;
y = y0;
f1 = feval(FunJac,tk,yk,varargin{:});
[f2,J2] = feval(FunJac,t,y,varargin{:});
funeval = 2;
R = y - 0.5*h*(f1+f2) - yk;
I = eye(size(yk,1));

while ((k < maxit) && (norm(R,'inf') > tol))
    k = k+1;
    dRdy = I - 0.5*h*J2;
    dy = dRdy\R;
    y = y - dy;
    f1 = feval(FunJac,t,yk,varargin{:});
    [f2,J2] = feval(FunJac,t,y,varargin{:});
    funeval = funeval + 2;
    R = y - 0.5*h*(f1+f2) - yk;
end

end