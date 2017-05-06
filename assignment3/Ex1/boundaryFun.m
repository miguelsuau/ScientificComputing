function u = boundaryFun(x,t)

alpha = [1 4 16];
u = 0;
for i = 1:3
    u = u + cos(alpha(i)*x)*exp(-alpha(i)^2*t);
end

end