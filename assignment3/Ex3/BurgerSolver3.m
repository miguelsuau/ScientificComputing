function U = BurgerSolver3(boundaryFun,h,k,epsilon,tmax)
M = ceil(2/h);
N = ceil(tmax/k);
x = linspace(-1,1,M+1);
t = linspace(0,tmax,N+1);
U = zeros(M+1,N+1);
% Build A0
A0 = diag(-2*ones(M+1,1)) + diag(ones(M,1),-1) + diag(ones(M,1),1);
A0(1,:) = zeros(1,M+1);
A0(end,:) = zeros(1,M+1);
% Build B0
B0 = diag(-1*ones(M+1,1)) + diag(ones(M,1),1);
B0(1,:) = zeros(1,M+1);
B0(end,:) = zeros(1,M+1);
I = eye(M+1);
% Obtain initial and boundary conditions
U(:,1) = boundaryFun(x,0,epsilon);
g = [boundaryFun(-1,t,epsilon); zeros(M-1,N+1); boundaryFun(1,t,epsilon)];
for k=2:N+1
    U(:,k) = (I-0.5*epsilon*k/h^2*A0)\((I+0.5*epsilon*k/h^2*A0)*U(:,k-1)-U(:,k-1).*(k/h*B0*U(:,k-1))+g(:,k)-g(:,k-1));
end
[T,X] = meshgrid(t,x);
mesh(X,T,U)
end