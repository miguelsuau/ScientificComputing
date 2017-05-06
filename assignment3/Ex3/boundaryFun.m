function u = boundaryFun(x,t,epsilon)
% This function allows to evaluate the solution of the viscous Burger's
% equation at the boundaries
    u = -tanh((x+0.5-t)/(2*epsilon)) + 1;
end