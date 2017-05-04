function U = BurgerSolver(boundaryFun,h,k,epsilon)
M = ceil(2/h);
N = ceil(1/k);
A0 = diag(-2*ones(M+1,1)) + diag(ones(M,1),-1) + diag(ones(M,1),1);
A0(1,:) = zeros(1,M+1);
A0(end,:) = zeros(1,M+1);
I = eye(M+1);
x = linspace(-1,1,M+1);
t = linspace(0,1,N+1);
U(:,1) = boundaryFun(x,0);
g(:,1) = [boundaryFun(-1,0); zeros(M-1,1); boundaryFun(1,0)];

for k=2:N
    g(:,k) = [boundaryFun(-1,t(k)); zeros(M-1,1); boundaryFun(1,t(k))];
    U(:,k) = (I-theta*mu*A0)\((I+(1-theta)*mu*A0)*U(:,k-1)+g(:,k)-g(:,k-1));
end


end