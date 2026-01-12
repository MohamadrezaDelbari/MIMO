%% Problem 2
Pt=10^-3;
Pn=10^-4./[1 0.25 0.09];
P1 = waterfill(Pt,Pn)

%% Problem 3
H=[0.2 -0.2 0.2;0.1 0.4 -1;1 0.2 1];
s=svd(H);
Pt=10;
Pn=1./s'.^2;
P2 = waterfill(Pt,Pn)
