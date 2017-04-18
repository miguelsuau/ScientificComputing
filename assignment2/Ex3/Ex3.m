

%% 3.2 Under/Over-relaxed Jacobi Smoothing
feig = @(p,q,h,omega) (1-omega) + omega*2/h^2*(cos(p*pi*h)-1+cos(q*pi*h)-1);
m = 100;
omega = linspace(0,2,100);
h = 1/(m+1);
eig = zeros(m/2-m,m/2-m);
maxeig = zeros(100,1);
for i = 1:100
    for p = m/2:m
        for q = m/2:m
            eig(p,q) = feig(p,q,h,omega(i));
        end
    end
    maxeig(i) = max(max(abs(eig)));
end