function [  ] = implot( alpha, beta, val, titleStr )
%IMPLOT Summary of this function goes here
%   Detailed explanation goes here
imagesc(alpha,beta,val,[0 1]);
grid on
colorbar
axis image
axis xy
xlabel('real','fontsize',14);
ylabel('imag','fontsize',14);
title(titleStr,'fontsize',14)

end

