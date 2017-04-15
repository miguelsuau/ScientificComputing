function F = form_rhs5(m,f,u)

h = 1/(m+1);

% form RHS for 5-point Laplacean

% use u to calculate the bounds
[X,Y] = meshgrid(linspace(0,1,m+2), linspace(0,1,m+2)); % account for outside of the grid
U = u(X,Y);

Iint = 2:m+1;
Jint = 2:m+1;
Xint = X(Iint,Jint);
Yint = Y(Iint,Jint);

F = f(Xint,Yint);

F(:,1) = F(:,1) - U(2:m+1,1)/h^2;
F(:,m) = F(:,m) - U(2:m+1,m+2)/h^2;
F(1,:) = F(1,:) - U(1,2:m+1)/h^2;
F(m,:) = F(m,:) - U(m+2,2:m+1)/h^2;

F = reshape(F,m*m,1);

end