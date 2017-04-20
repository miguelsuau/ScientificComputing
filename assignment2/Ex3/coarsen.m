function [ Rc ] = coarsen( R, m )
% Restriction
% see:
% https://cn.inside.dtu.dk/cnnet/filesharing/download/10824fda-2212-4f3c-8baf-54cc461f365f
% slide 9

mc = (m-1)/2;
R = reshape(R,m+2,m+2);
% idx = 2:2:length(R);
% Rc = 1/4*(R(idx+1,idx-1) + R(idx+1,idx+1) + R(idx-1, idx-1) + R(idx-1, idx+1));

Rc = R(1:2:end,1:2:end);

% R = reshape(R,m,m);
% ind = 2:2:length(R);
%  
% Rc = ( 4*R(ind,ind) + ...
%                     2*(R(ind-1,ind)+R(ind+1,ind)+ ...
%                        R(ind,ind-1)+R(ind,ind+1)) + ...
%                       (R(ind-1,ind-1)+R(ind-1,ind+1)+ ...
%                        R(ind+1,ind-1)+R(ind+1,ind+1)) )/16;

end

