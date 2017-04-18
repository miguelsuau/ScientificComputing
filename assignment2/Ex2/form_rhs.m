function [ F ] = form_rhs( m, f, u )
% form RHS for 9-point Laplacean
h = 1/(m+1);
hmult = 1/(6*h^2);
%hmult = 1; % used for debug with symbolic vars

[X,Y] = meshgrid(linspace(0,1,m+2), linspace(0,1,m+2)); % account for outside of the grid
Iint = 2:m+1;
Jint = 2:m+1;
Xint = X(Iint,Jint);
Yint = Y(Iint,Jint);

F = f(Xint,Yint);
% contribution from 3.19 in LeVeque
contrib = reshape( (((1/12)*h^2)*poisson5(m))*reshape(F,m*m,1) , m, m);
F = F + contrib;
% boundaries
U = u(X,Y);
U(2:m+1,2:m+1) = 0;
% Uncomment to test with symbolic variables to see if Matlab code matches the calculations on paper
% U = sym('U%d%d', [m+2 m+2]);
% U(2:m+1,2:m+1) = 0;
% F = sym('F%d%d', [m m]);

% common for all (similar to 5-point L RHS)
%  - x -
%  x - x
%  - x -
F(:,1) = F(:,1) - hmult * 4 * U(2:m+1, 1    ); % left column
F(:,m) = F(:,m) - hmult * 4 * U(2:m+1, m+2  ); % right column
F(1,:) = F(1,:) - hmult * 4 * U(1,     2:m+1); % top row
F(m,:) = F(m,:) - hmult * 4 * U(m+2,   2:m+1); % bottom row
% account for corner points of the 9-point L stencil 
%  x - x
%  - - -
%  x - x
for k=1:m
    F(k,1) = F(k,1) - hmult * U(k,1) - hmult * U(k+2,1);
    F(k,m) = F(k,m) - hmult * U(k,m+2) - hmult * U(k+2,m+2);
    
    F(1,k) = F(1,k) - hmult * U(1,k) - hmult * U(1,k+2);
    F(m,k) = F(m,k) - hmult * U(m+2,k) - hmult * U(m+2,k+2);
end
% account for the contribution from corners of the grid that was added
% twice in the previous step 
F(1,1) = F(1,1) + hmult * 1 * U(1,1);
F(m,m) = F(m,m) + hmult * 1 * U(m+2,m+2);
F(m,1) = F(m,1) + hmult * 1 * U(m+2,1);
F(1,m) = F(1,m) + hmult * 1 * U(1,m+2);
% finalize
% FF = f(Xint,Yint);
% FC = ((1/12)*h^2)*poisson5(m)*reshape(FF,m*m,1);
% F = FF(:) + FC + F(:);
F = reshape(F, m*m, 1);

end

