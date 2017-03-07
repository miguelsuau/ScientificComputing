function [f,J] = VanderPolfunjac(t,y,mu)

% Evaluate function
f = [y(2); mu*(1-y(1)^2)*y(2)-y(1)];
% Evaluate Jacobian
J = [0 1; -y(2)*mu*2*y(1)+1 mu*(1-y(1)^2)];

end