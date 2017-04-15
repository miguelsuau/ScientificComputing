function A = poisson5(m)

h = 1/(m+1);
e = ones(m,1);
S = spdiags([e -2*e e], [-1 0 1], m, m);
I = speye(m);
A = kron(I,S)+kron(S,I);
A = A/h^2;

end
