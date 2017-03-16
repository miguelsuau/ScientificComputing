close all; clear all;
addpath('../Ex3');

%% Butcher tableau
gamma = 1-1/sqrt(2);
a31 = (1-gamma)/2;
AT = [0 gamma a31;0 gamma a31;0 0 gamma]; A = AT';
c  = [0; 2*gamma; 1];
b  = AT(:,3);
bhat = [    (6*gamma-1)/(12*gamma); ...
            1/(12*gamma*(1-2*gamma)); ...
            (1-3*gamma)/(3*(1-2*gamma))    ];
d  = bhat-b;
s = 3;

%% Stability related calculations
alpha = -10:0.01:10;
beta = -10:0.01:10;
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

%% Stability plot etc.
figure
subplot(2,2,1)
implot(alpha, beta, absR, '|R(z)| (stability plot)')
subplot(2,2,2)
implot(alpha, beta, absEhat, '|E(z)| (est. error)')
subplot(2,2,3)
implot(alpha, beta, absE, '|R(z) - exp(z)| (true error)')
subplot(2,2,4)
implot(alpha, beta, absEhatmE, '|E(z) - R(z) + exp(z)|')

%% 
syms z gamma; 
assume(gamma,'real');
a31 = (1-gamma)/2;
AT = [0 gamma a31;0 gamma a31;0 0 gamma]; A = AT';
b  = AT(:,3);

nR = det(eye(3)-z*A+z*ones(3,1)*b')
dR = det(eye(3)-z*A)

limit(abs(nR/dR),z,-inf)