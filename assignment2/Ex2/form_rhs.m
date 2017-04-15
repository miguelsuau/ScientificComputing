function [ F ] = form_rhs( m, f, u )
% form RHS for 9-point Laplacean
h = 1/(m+1);
hmult = 1/(6*h^2);

[X,Y] = meshgrid(linspace(0,1,m+2), linspace(0,1,m+2)); % account for outside of the grid
Iint = 2:m+1;
Jint = 2:m+1;
Xint = X(Iint,Jint);
Yint = Y(Iint,Jint);

F = f(Xint,Yint);
U = u(X,Y);

% account for BC
% common for all
F(:,1) = F(:,1) - hmult * 2 * U(2:m+1, 1    ); % left column
F(:,m) = F(:,m) - hmult * 2 * U(2:m+1, m+2  ); % right column
F(1,:) = F(1,:) - hmult * 2 * U(1,     2:m+1); % top row
F(m,:) = F(m,:) - hmult * 2 * U(m+2,   2:m+1); % bottom row
% one inner point shifted
F(2:(m-1), 1) = F(2:(m-1), 1) - hmult * 1 * U(3:m, 1  ); % left column
F(2:(m-1), m) = F(2:(m-1), m) - hmult * 1 * U(3:m, m+2); % right column
F(1, 2:(m-1)) = F(1, 2:(m-1)) - hmult * 1 * U(1,   3:m); % top row
F(m, 2:(m-1)) = F(m, 2:(m-1)) - hmult * 1 * U(m+2, 3:m); % bottom row
% corner points
F(1,1) = F(1,1) - hmult * U(1,1);
F(1,m) = F(1,m) - hmult * U(1,m+2);
F(m,1) = F(m,1) - hmult * U(m+2,1);
F(m,m) = F(m,m) - hmult * U(m+2,m+2);


% finalize according to (3.19)
F = reshape(F, m*m, 1);
F = F + ((1/12)*h^2)*poisson5(m)*F;


end

