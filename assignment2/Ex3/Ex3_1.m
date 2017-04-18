m = 30;
h = 1/(m+1);
[X,Y] = meshgrid(0:h:1,0:h:1);

% Change to -Au=-F (`pcg` requires positive definite system)
f0 = @(X,Y) -16 * (pi ^ 2) * (...
    cos((4 * pi * X .* Y)) .* (X .^ 2)...
    + cos((4 * pi * X .* Y)) .* (Y .^ 2)...
    + 2 * sin((4 * pi * (X + Y)))...
);
b = -f0(X,Y);
Afun = @(u) reshape(-Amult(u,m+2), (m+2)*(m+2), 1);

u = pcg(Afun,b(:),1e-6,10);
% max(max(abs(Amult(u,m+2) - reshape(poisson5(m+2)*u,m+2,m+2))))

figure(1)
surf(X,Y,reshape(u,m+2,m+2));