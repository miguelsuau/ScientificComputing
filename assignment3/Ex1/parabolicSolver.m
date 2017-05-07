function U = parabolicSolver(boundaryFun,h,k,theta,mu)
% The following function solves the heat diffusion equation using the
% theta-scheme. It takes as input the two step sizes h and k and the
% parameters theta and mu and returns the solution at every grid point
% Define the grid size
M = ceil(2/h); 
N = ceil(1/k); 
% Build A0
A0 = diag(-2*ones(M+1,1)) + diag(ones(M,1),-1) + diag(ones(M,1),1);
A0(1,:) = zeros(1,M+1);
A0(end,:) = zeros(1,M+1);
I = eye(M+1);
x = linspace(-1,1,M+1);
t = linspace(0,1,N+1);
% Obtain initial and boundary conditions
U(:,1) = boundaryFun(x,0);
g = [boundaryFun(-1,t); zeros(M-1,N+1); boundaryFun(1,t)];
% Solve the ODEs system
for k=2:N+1
    U(:,k) = (I-theta*mu*A0)\((I+(1-theta)*mu*A0)*U(:,k-1)+g(:,k)-g(:,k-1));
end
end