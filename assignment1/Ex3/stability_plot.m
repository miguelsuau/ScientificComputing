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
d = [ 7/24;    
      -7/18;    
      7/72 ];
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
        absR(kimag,kreal) = abs(R); % stability plot R(z)
        absEhatmE(kimag,kreal) = abs(EhatmE); % err's est err E(z) - (R(z) - exp(z))
        absEhat(kimag,kreal) = abs(Ehat); % est err E(z)
        absE(kimag,kreal) = abs(E); % true err R(z) - exp(z)
        absF(kimag,kreal) = abs(f); % exp(z)
    end
end

figure
subplot(4,2,1)
implot(alpha, beta, absR, '3rd order |R(z)| (stability plot)')
subplot(4,2,3)
implot(alpha, beta, absEhat, '3rd order |E(z)| (est. error)')
subplot(4,2,5)
implot(alpha, beta, absE, '3rd order |R(z) - exp(z)| (true error)')
subplot(4,2,7)
implot(alpha, beta, absEhatmE, '3rd order |E(z) - R(z) + exp(z)|')

% embedded
b = [1/8; 1/2; 3/8]; % bhat
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

subplot(4,2,2)
implot(alpha, beta, absR, 'Embedded |R(z)| (stability plot)')
subplot(4,2,4)
implot(alpha, beta, absEhat, 'Embedded |E(z)| (est. error)')
subplot(4,2,6)
implot(alpha, beta, absE, 'Embedded |R(z) - exp(z)| (true error)')
subplot(4,2,8)
implot(alpha, beta, absEhatmE, 'Embedded |E(z) - R(z) + exp(z)|')