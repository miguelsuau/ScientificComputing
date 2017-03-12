function [f,J] = TestEquation(t,x,lambda)

% Evaluate function
f = lambda*x;
% Evaluate Jacobian
J = lambda;

end