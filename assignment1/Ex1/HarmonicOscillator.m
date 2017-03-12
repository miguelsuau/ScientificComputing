function [f,J] = HarmonicOscillator(t,x)

% Evaluate function
f = [x(2); -x(1)];
% Evaluate Jacobian
J = [0 1; 1 0];

end