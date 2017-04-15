function A = poisson9(m)

% This approach has problem with overlapping corners, but it can be used
% TODO: If there is time try this one also, just modify the S matrix to
% account for the corners in the first and last row
% e = ones(m,1);
% S = spdiags([e 16*e -30*e 16*e e], [-2 -1 0 1 2], m, m);
% I = speye(m);
% A = kron(I,S)+kron(S,I);
h = 1/(m+1);
e = ones(m,1);
S = -spdiags([e 10*e e], [-1 0 1], m, m);
I = spdiags([-1/2*e e -1/2*e], [-1 0 1], m, m);
A = kron(I,S)+kron(S,I);
A = A/(6*h^2);
end