%% Resilient Communication Systems (RCS)
% Mohamadreza Delbari
% This code is for channel estimation

%% Problem 1
clear variables
clc
M=2; % # of Tx antenna
N=3; % # of Rx antenna
Pt=2;
SNR=10000;
H=[0.1+0.1*1i	-0.3-0.2*1i;0.2+0.5*1i	0.4+0.1*1i;1i	0.4];
n=sqrt(Pt/SNR)*randn(N,1)+1i*sqrt(Pt/SNR)*randn(N,1);
n=0*n;
x1=sqrt(1/(2))*[1 1i;1i 1];
x2=sqrt(1/(4))*[1 1i;1i 1];
x3=sqrt(1)*[0.99 0.2i;0 0.99];
x4=sqrt(1)*[0.1 1i;0 0.99];

y1=H*x1+n;
y2=H*x2+n;
y3=H*x3+n;
y4=H*x4+n;

H1=(M/Pt)*y1*x1';
H2=(M/Pt)*y2*x2';
H3=(M/Pt)*y3*x3';
H4=(M/Pt)*y4*x4';

%% Problem 2
clear variables
clc
SNRdB=0:5:30;
SNR=db2pow(SNRdB);
M=1; N=3;
J1=M^2*N./SNR;
M=2; N=3;
J2=M^2*N./SNR;
AA=[SNRdB' J1' J2'];
semilogy(SNRdB,J1,SNRdB,J2)
