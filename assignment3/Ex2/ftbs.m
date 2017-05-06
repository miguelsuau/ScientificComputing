function [ Uout ] = ftbs( u, cr, wave_periods, n_wave_length )
%FTBS Summary of this function goes here
%   Detailed explanation goes here

a = 1/2;
innerpts = 100; % todo pass in a, \Delta x = L / 100
x = linspace(-1,1,innerpts+2);
t = linspace(0,a^(-1)*wave_periods,a^(-1)*n_wave_length*wave_periods);

[X,T] = meshgrid(x,t);
U = u(X,T);

m = length(t);

ux0 = u(x,0);
ux0(abs(ux0) < 1e-3) = 0;

u1t = u(1,t);
u1t(abs(u1t) < 1e-3) = 0;

[um, un] = size(U);
Uout = zeros(size(U));
for k=1:(2*(a)^1 * n_wave_length):length(t) % ?
    Uout(k,:) = ux0; % periodic BCs
end
%Uout(1,:) = u(x,0);
Uout(:,1) = u1t; % at fixed x
Uout(:,end) = u1t; % at fixed x

% for n=2:m
%     if mod(n, (a^(-1))*n_wave_length) == 1 % skip the period boundary at 1, 101, ...
%         %fprintf('Skipping at idx %d\n', n);
%         continue
%     end
%     
%     Uprev = Uout(n-1,:); % previous U
%     Un = Uprev;
%     Un(2:end-1) = 0; % keep the BC for fixed x untouched
%     for j=2:length(Un)-1
%         Un(j) = Uprev(j) - cr*(Uprev(j) - Uprev(j-1));
%     end
%     Uout(n,:) = Un; % add current result to output
% end

% for n=2:um
%     Uprev = Uout(n-1,:); % previous U
%     Un = Uprev;
%     %Un(2:end-1) = 0; % keep the BC for fixed x untouched
%     for j=2:length(Un)
%         Un(j) = Uprev(j) - cr*(Uprev(j) - Uprev(j-1));
%     end
%     
%     % Show the moving wave
%     if (n/n_wave_length <= 1.8)
%         figure(3); clf; 
%         plot(x,Un,'bo-',x,U(n,:),'r'); 
%         axis([-1 1 -2 2]); 
%         title({[sprintf('Time = %2.2f s', n/n_wave_length)]});
%         ylabel('u'); xlabel('x');
%         pause(0.0001);
%     end
%     
%     Uout(n,:) = Un; % add current result to output
% end

x = linspace(-1,1,innerpts+1); % innerpts = 100
t = linspace(0,a^(-1)*wave_periods,a^(-1)*n_wave_length*wave_periods);

[X,T] = meshgrid(x,t);
U = u(X,T);
Uout = U;

u0 = u(x,0);
ut = u(1,t); % u(1,t) == u(-1,t)
u = u0; % initialize the first BC u(x,0)
idx = 2:length(x);

for tt=2:um % skip 1st index, go all the way to the end
    u(1) = ut(tt); % BC u(1,t)
    u(idx) = u(idx) - cr*(u(idx) - u(idx-1));
    
    % Show the moving wave
    if (tt/n_wave_length <= 1.8)
        figure(3); clf;
        plot(x,u,'bo-',x,U(tt,:),'r'); 
        axis([-1 1 -2 2]); 
        title({[sprintf('Time = %2.2f s', tt/n_wave_length)]});
        ylabel('u'); xlabel('x');
        pause(0.0001);
    end
    Uout(tt,:) = u;
end

end

