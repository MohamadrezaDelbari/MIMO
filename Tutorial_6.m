%% Resilient Communication Systems (RCS)
% Mohamadreza Delbari
% This code compare the performance of different non-linear MIMO detectors
clear variables
clc
%% Problem 1
clear variables
clc
M=2; % # of Tx antenna
N=3; % # of Rx antenna
H=[0.2	-0.2;0.1	0.4;1	0.2];
SNR=10;
L_ZF=[1 0;0.83 1];
D_ZF=[0.88 0; 0 0.24];
W_FF_ZF=inv(D_ZF)*inv(L_ZF')*H';
W_FB_ZF=L_ZF-eye(M);
L_MMSE=[1 0;0.59 1];
D_MMSE=[1.03 0; 0 0.34];
W_FF_MMSE=inv(D_MMSE)*inv(L_MMSE')*H';
W_FB_MMSE=L_MMSE-eye(M);
y=[0.6;-0.1;1.3];
r_ZF=W_FF_ZF*y;
r_MMSE=W_FF_MMSE*y;
%% Problem 2
M=2; % # of Tx antenna
N=3; % # of Rx antenna
Aset=[-1-1i 1-1i 1+1i -1+1i]; % BPSK modulation alphabet
t=1; % SNR counter
I=100000; % Number of iterations
Ptx=1; % Normalized transmitter power
for SNRdb=-10:5:20
   S_ZF=0; % Sum of error bits
   S_MMSE=0; % Sum of error bits
   S_ML=0;
   S_SP=0;
   SNR=db2pow(SNRdb); % Changing dB to power
   sigma2_n=Ptx/SNR; % Noise power
   if SNRdb==20
       I=200000;
   end
   for ii=1:I
   H=(randn(N,M)+1i*randn(N,M))*db2pow(0); % Rayleigh channel.
   x=(floor(rand(M,1)*2)-0.5)*2+1i*(floor(rand(M,1)*2)-0.5)*2; % From Aset
   n=randn(N,1)*sqrt(sigma2_n)+1i*randn(N,1)*sqrt(sigma2_n); % Noise generation
   y=H*x+n; % Received signal
   %% Non-linear detection (soft detection)
   tic
   [L_ZF,D_ZF]=ldl(H'*H); %Cholesky decomposition of H'*H for ZF-SIC
   L_ZF=L_ZF';
   r_ZF=inv(D_ZF)*inv(L_ZF)'*H'*y; %Signal after feedforward ZF-SIC
   W_FB_ZF=L_ZF-eye(M); %Feedback Filter
   T_ZF(ii)=toc;
   tic
   [L_MMSE,D_MMSE]=ldl(H'*H+SNR^-1*eye(M));%Cholesky decomposition for MMSE-SIC
   L_MMSE=L_MMSE';
   r_MMSE=inv(D_MMSE)*inv(L_MMSE)'*H'*y; %Signal after feedforward ZF-SIC
   W_FB_MMSE=L_MMSE-eye(M); %Feedback Filter
   T_MMSE(ii)=toc;
   if norm(diag(diag(L_MMSE))-eye(M))==0 && norm(diag(diag(L_ZF))-eye(M))==0
   x_hat_ZF=zeros(M,1);
   x_hat_MMSE=zeros(M,1);
   for i=M:-1:1
       tic
       d_ZF(i)=r_ZF(i)-W_FB_ZF(i,:)*x_hat_ZF;
       x_hat_ZF(i)=(real(d_ZF(i))>0)-(real(d_ZF(i))<0)+1i*((imag(d_ZF(i))>0)-(imag(d_ZF(i))<0));
       T_ZF(ii)=T_ZF(ii)+toc;
       tic
       d_MMSE(i)=r_MMSE(i)-W_FB_MMSE(i,:)*x_hat_MMSE;
       x_hat_MMSE(i)=(real(d_MMSE(i))>0)-(real(d_MMSE(i))<0)+1i*((imag(d_MMSE(i))>0)-(imag(d_MMSE(i))<0));
       T_MMSE(ii)=T_MMSE(ii)+toc;
   end
   %% ML and Sphere decoding
   minML=100000;
   tic
   for i_1=1:length(Aset)
   for i_2=1:length(Aset)
   %for i_3=1:length(Aset)
       %x_ML=[Aset(i_1);Aset(i_2);Aset(i_3)];
       x_ML=[Aset(i_1);Aset(i_2)];
       if norm(y-H*x_ML)<minML
           minML=norm(y-H*x_ML);
           x_ML_final=x_ML;
       end
   %end
   end
   end
   T_ML(ii)=toc;
   % minSP=100000;
   % tic
   % [Q_SP,R_SP]=qr(H);
   % yhat=Q_SP'*y;
   % for i_1=1:length(Aset)
   % f1=norm(yhat(3)-R_SP(3,3)*Aset(i_1),2)^2;
   % if f1<minSP
   % for i_2=1:length(Aset)
   %     f2=norm(yhat(2)-R_SP(2,2:3)*[Aset(i_2);Aset(i_1)],2)^2;
   %     if f1+f2<minSP
   %         for i_3=1:length(Aset)
   %             f3=norm(yhat(1)-R_SP*[Aset(i_3);Aset(i_2);Aset(i_1)],2)^2;
   %             if f1+f2+f3<minSP
   %                 minSP=f3+f2+f1;
   %                 x_SP=[Aset(i_3);Aset(i_2);Aset(i_1)];
   %             end
   %         end
   %     end
   % end
   % end
   % end
   %
   % T_SP(ii)=toc;
   %% Sum of error
   S_ZF=S_ZF+sum(x_hat_ZF~=x);
   S_MMSE=S_MMSE+sum(x_hat_MMSE~=x);
   S_ML=S_ML+sum(x_ML_final~=x);
   %S_SP=S_SP+sum(x_SP~=x);
   else
       ii=ii-1;
   end
   end
   Pr_ZF(t)=S_ZF/(2*M*I);
   Pr_MMSE(t)=S_MMSE/(2*M*I);
   Pr_ML(t)=S_ML/(2*M*I);
   Pr_SP(t)=S_SP/(2*M*I);
   t=t+1;
   SNRdb
end
t=(-10:5:20);
semilogy(t,Pr_ML,t,Pr_ZF,t,Pr_MMSE)
A=[t' Pr_ML' Pr_ZF' Pr_MMSE'];

