function K = elemKVEM(EX,mu,coor)
Nv = size(coor,1);
% plane stress
C = EX/(1-mu^2)*[1,mu,0;mu,1,0;0,0,(1-mu)/2];

[dM,Pis,Pi,area] = calculatePi(coor);

Pis = blkdiag(Pis,Pis);
Pi = blkdiag(Pi,Pi);

A1 = [1,0;0,0;0,1];
A2 = [0,0;0,1;1,0];
A = [A1,A2];

dM1 = blkdiag(dM',dM');

G0 = area*dM1'*A'*C*A*dM1;

Kc = Pis'*G0*Pis;

II = eye(size(Pi));
% alpha = 1/10*sum(volume)*trace(Kc);
alpha = 1/Nv*trace(Kc);
Ks = alpha * ((II - Pi)' * (II - Pi));

AK = Kc+Ks;
K = AK;
