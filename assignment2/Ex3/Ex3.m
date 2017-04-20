%% 3.2 Under/Over-relaxed Jacobi Smoothing

% eigenvalues expression 4.90 REVIEW
feig = @(p,q,h,omega) (1-omega)+omega.*(cos(p.*pi*h)+cos(q.*pi*h));

m = 1000;
h = 1/(m+1);

% obtain all combinations of p and q
p = m/2:m;
q = m/2:m;
[P,Q] = meshgrid(p,q);

% memory allocation
maxeig = zeros(100,1);

% for each value of omega find the maximum eigenvalue
omega = linspace(0,2,1000);
for i = 1:1000
    eig = feig(P,Q,h,omega(i));
    maxeig(i) = max(max(abs(eig)));
end
% choose omega that makes maxeig minimum
[minmax,idx] = min(maxeig);
omegaopt = omega(idx);
plot(omega,maxeig,'linewidth',1.6)
hold on
plot(omegaopt,minmax,'ro','linewidth',1.6)
legend('maximum eigenvalue','optimal \omega','location','northwest')
xlabel('\omega')
ylabel('max(\rho_{p,q}(\omega))','interpreter','tex')
print('omegaopt','-dpng')