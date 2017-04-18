

%% 3.2 Under/Over-relaxed Jacobi Smoothing
feig = @(p,q,h,omega) (1-omega)+omega.*2/h^2*(cos(p.*pi*h)-1+cos(q.*pi*h)-1);
m = 10000;
omega = linspace(0,2,100);
h = 1/(m+1);
p = m/2:m;
q = m/2:m;
[P,Q] = meshgrid(p,q);
maxeig = zeros(100,1);

for i = 1:100
    eig = feig(P,Q,h,omega(i));
    maxeig(i) = max(max(abs(eig)));
end