function [y,funeval] = NewtonsImplicitEuler(FunJac,tk,yk,h,y0,tol,maxit,varargin)

k = 0;
t = tk + h;
y = y0;
[f,J] = feval(FunJac,t,y,varargin{:});
R = y -h*f - yk;
I = eye(size(yk,1));
funeval = 1;
while ((k < maxit) && (norm(R,'inf') > tol))
    k = k+1;
    dRdy = I - h*J;
    dy = dRdy\R;
    y = y - dy;
    [f,J] = feval(FunJac,t,y,varargin{:});
    funeval = funeval + 1;
    R = y - h*f - yk;
end

end
    