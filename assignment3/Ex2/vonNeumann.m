
v = 0.8; % Courant number
% g(xi)
% a circle, centered at 1-cr with radius cr
% 1-v-v*(cos(phi) - 1i*sin(phi))
g = @(phi) 1 - v - v*exp(1i*phi); % phi = xi*dx

xi = 2*pi/1; % 2*pi/L, L = 1
phi = xi*0.01; % xi*dx, dx = 1/100

G = g(phi);

time_steps = 80/(16/1000);

fprintf(' Von Neumann analysis\n%s\n------\n', 'g(xi) = 1-v-v*(cos(xi*dx) - 1i*sin(xi*dx)); phi = xi*dx');
fprintf('Courant number = %f\n', v);
fprintf('xi = %f\n', xi);
fprintf('phi = %f\n', phi);
fprintf('g(2*pi): Re = %f, Im = %f\n', real(G), imag(G));
fprintf('Abs(g(xi)) = %f\n', abs(G));
fprintf('Angle(g(xi)) = %f radians, %f degrees\n', angle(G), rad2deg(angle(G)));

fprintf('------\nTime steps = %d\n', time_steps)
fprintf('Diffusion = %f\n', abs(G^time_steps));
fprintf('y = %f radians, %f degrees\n', angle(G), rad2deg(angle(G)));
fprintf('y_exact = -Cr*phi = %f\n', -v*phi);

