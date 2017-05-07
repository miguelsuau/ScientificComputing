function [ combined ] = nonuniformgrid( points )
% Generate non-uniform grid between -1 and 1, where most of the points are
% centered around 0. 
% Number of points has to be even

if mod(points,2) == 1
    points = points - 1;
end
n = points/2;
r = 0.8;

lg = [0, (r - 1)/(r^(n-1) - 1)*cumsum(r.^(0:(n-2)))];
rg = 1+[0, (r - 1)/(r^(n-1) - 1)*cumsum(r.^((n-2):-1:0))];

if mod(points,2) == 1
    combined = [lg(:); 0; rg(:)]'; % combine left and right vectors
else 
    combined = [lg(:); rg(:)]'; % combine left and right vectors
end
combined = combined - 1;    % center around 0

% Uncomment to visualize the spacing
%figure(1);
%plot(combined, zeros(length(combined)), 'o-');

end