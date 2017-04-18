function F = form_rhs5(m,f,u)

h = 1/(m+1);

% form RHS for 5-point Laplacean
[X,Y] = meshgrid(0:h:1, 0:h:1);
f = f(X,Y);
u = u(X,Y);
f(2,2:end-1) = f(2,2:end-1) - u(1,2:end-1)*(m+1)^2;
f(end-1,2:end-1) = f(end-1,2:end-1) - u(end,2:end-1)*(m+1)^2;
f(2:end-1,2) = f(2:end-1,2) - u(2:end-1,1)*(m+1)^2;
f(2:end-1,end-1) = f(2:end-1,end-1) - u(2:end-1,end)*(m+1)^2;
F = f(2:end-1,2:end-1);


% % use u to calculate the bounds
% [X,Y] = meshgrid(linspace(0,1,m+2), linspace(0,1,m+2)); % account for outside of the grid
% U = u(X,Y);
% 
% Iint = 2:m+1;
% Jint = 2:m+1;
% Xint = X(Iint,Jint);
% Yint = Y(Iint,Jint);
% 
% F = f(Xint,Yint);
% b = zeros(m,m);
% b(:,1) = b(:,1) + U(2:m+1,1);
% b(:,m) = b(:,m) + U(2:m+1,m+2);
% b(1,:) = b(1,:) + U(1,2:m+1);
% b(m,:) = b(m,:) + U(m+2,2:m+1);
% % F(:,1) = F(:,1) - U(2:m+1,1)*(1/h^2);
% % F(:,m) = F(:,m) - U(2:m+1,m+2)*(1/h^2);
% % F(1,:) = F(1,:) - U(1,2:m+1)*(1/h^2);
% % F(m,:) = F(m,:) - U(m+2,2:m+1)*(1/h^2);
% b = (1/h^2)*reshape(b,m*m,1);
% F = reshape(F,m*m,1);
% F = F-b;

end