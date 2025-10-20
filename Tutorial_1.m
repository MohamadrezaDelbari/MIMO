clear variables
clc
%_____________Problem 1 (Part a)________________
syms x
A=[x -x;-x x]
B=[-x 2*x;2*x x]
eig(A)
det(A)
eig(B)
det(B)
det(A*B)
trace(A+B)
%_____________Problem 2_________________________
A=[1 0.5 -0.5;-1 1 -1]
B=[-1 1 -1;-0.5 0.5 1]
[V,D,U]=svd(A)
R1=fixed.qlessQR(A)