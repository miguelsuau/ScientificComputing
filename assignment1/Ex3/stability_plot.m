close all
clear all

alpha = -4:0.01:4;
beta = -4:0.01:4;
A = [   
        0               1/4              -7/5;     
        0               0                12/5;     
        0               0                0  
]';
b = [   
        -1/6;     
        8/9;     
        5/18;    
];
c = [0; 1/3; 1];
d = [ -7/24;    
       7/18;    
      -7/72 ];
nreal = length(alpha);
nimag = length(beta);
I = eye(size(A));
e = ones(size(A,1),1);

for kreal = 1:nreal
    for kimag = 1:nimag
        z = alpha(kreal) + 1i*beta(kimag);
        tmp = (I-z*A)\e;
        R = 1 + z*b'*tmp;
        Ehat = z*d'*tmp;
        f = exp(z);
        E = R-f;
        EhatmE = Ehat-E;
        absR(kimag,kreal) = abs(R);
        absEhatmE(kimag,kreal) = abs(EhatmE);
        absEhat(kimag,kreal) = abs(Ehat);
        absE(kimag,kreal) = abs(E);
        absF(kimag,kreal) = abs(f);
    end
end

figure
fs = 14;
subplot(1,2,1)
imagesc(alpha,beta,absR,[0 1]);
grid on
colorbar
axis image
axis xy
xlabel('real','fontsize',fs);
ylabel('imag','fontsize',fs);
title('3rd order |R(z)|','fontsize',fs)

% embedded
b = b - d; % bhat
for kreal = 1:nreal
    for kimag = 1:nimag
        z = alpha(kreal) + 1i*beta(kimag);
        tmp = (I-z*A)\e;
        R = 1 + z*b'*tmp;
        Ehat = z*d'*tmp;
        f = exp(z);
        E = R-f;
        EhatmE = Ehat-E;
        absR(kimag,kreal) = abs(R);
        absEhatmE(kimag,kreal) = abs(EhatmE);
        absEhat(kimag,kreal) = abs(Ehat);
        absE(kimag,kreal) = abs(E);
        absF(kimag,kreal) = abs(f);
    end
end

subplot(1,2,2)
imagesc(alpha,beta,absR,[0 1]);
grid on
colorbar
axis image
axis xy
xlabel('real','fontsize',fs);
ylabel('imag','fontsize',fs);
title('Embedded |R(z)|','fontsize',fs)