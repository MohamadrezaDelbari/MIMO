%% Resilient Communication Systems (RCS)
% Mohamadreza Delbari
% This code compare the performance of different linear MIMO detectors
%% Problem 2
clear variables
clc
M=3; % # of Tx antenna
N=4; % # of Rx antenna
Aset=[-1 1]; % BPSK modulation alphabet
t=1; % SNR counter
I=1000000; % Number of iterations
Ptx=1; % Normalized transmitter power
for SNRdb=-10:10:20
   S_MRC=0; % Sum of error bits
   S_ZF=0; % Sum of error bits
   S_MMSE=0; % Sum of error bits
   SNR=db2pow(SNRdb); % Changing dB to power
   sigma2_n=Ptx/SNR; % Noise power
   for i=1:I
   H=(randn(N,M)+1i*randn(N,M))*db2pow(0); % Rayleigh channel.
   x=(floor(rand(M,1)*2)-0.5)*2; % From Aset
   n=randn(N,1)*sqrt(sigma2_n); % Noise generation
   y=H*x+n; % Received signal
   %% Linear detection (soft detection)
   x_soft_MRC=inv(diag(diag(H'*H)))*H'*y;
   x_soft_ZF=inv(H'*H)*H'*y;
   x_soft_MMSE=inv(H'*H+SNR^-1*eye(M))*H'*y;
   %% Hard detection
   x_hard_MRC=(x_soft_MRC>0)-(x_soft_MRC<0);
   x_hard_ZF=(x_soft_ZF>0)-(x_soft_ZF<0);
   x_hard_MMSE=(x_soft_MMSE>0)-(x_soft_MMSE<0);
   %% Sum of error
   S_MRC=S_MRC+sum(x_hard_MRC~=x);
   S_ZF=S_ZF+sum(x_hard_ZF~=x);
   S_MMSE=S_MMSE+sum(x_hard_MMSE~=x);
   end
   Pr_MRC(t)=S_MRC/(M*I);
   Pr_ZF(t)=S_ZF/(M*I);
   Pr_MMSE(t)=S_MMSE/(M*I);
   t=t+1;
   SNRdb
end
t=(-10:10:20);
semilogy(t,Pr_MRC,t,Pr_ZF,t,Pr_MMSE)
A=[t' Pr_MRC' Pr_ZF' Pr_MMSE'];
%% Problem 1
clear variables
clc
M=2; % # of Tx antenna
N=3; % # of Rx antenna
H=[0.2	-0.2;0.1	0.4;1	0.2];
SNR=10;
W_MRC=inv(diag(diag(H'*H)))*H';
W_ZF=inv(H'*H)*H';
W_MMSE=inv(H'*H+SNR^-1*eye(M))*H';
y=[0.6;-0.1;1.3];
x_MRC=W_MRC*y;
x_ZF=W_ZF*y;
x_MMSE=W_MMSE*y;

