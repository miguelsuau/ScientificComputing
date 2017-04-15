function [G,J] = FunJac(U,h)
% This function builds the vector G and its Jacobian to solve the BVP using
% Newton's method.

G = 0.01/h^2*( U(1:end-2) - 2*U(2:end-1) + U(3:end)) ... 
    + U(2:end-1).*((U(3:end) - U(1:end-2))/(2*h) - 1);

J = 1/h^2 * (diag(-2*0.01 + 0.5*h*(U(3:end) - U(1:end-2))-h^2) ...
    + diag(0.01 - 0.5*h*U(3:end-1),-1) + diag(0.01 + 0.5*h*U(2:end-2),1));

end