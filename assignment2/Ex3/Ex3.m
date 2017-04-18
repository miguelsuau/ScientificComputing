

%% 3.2 Under/Over-relaxed Jacobi Smoothing

% eigenvalues expression 4.90 REVIEW
feig = @(p,q,h,omega) (1-omega)+omega.*2/h^2*(cos(p.*pi*h)-1+cos(q.*pi*h)-1);

m = 1000;
h = 1/(m+1);

% obtain all combinations of p and q
p = m/2:m;
q = m/2:m;
[P,Q] = meshgrid(p,q);

% memory allocation
maxeig = zeros(100,1);

% for each value of omega find the maximum eigenvalue
omega = linspace(0,2,100);
for i = 1:100
    eig = feig(P,Q,h,omega(i));
    maxeig(i) = max(abs(eig));
end
% choose omega that makes maxeig minimum
[~,idx] = min(maxeig);
omegaopt = omega(idx);

%%
% Jacobi Smoothing using the optimal omega
% U initial guess
U = zeros(m);
F = zeros(m);
omegaopt = 2/3;
maxit = 100;
for i = 1:maxit
   U = smooth(U,omegaopt,m,F); 
end
