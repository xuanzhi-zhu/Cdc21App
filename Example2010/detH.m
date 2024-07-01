clear all

A = [0 1;-1 0]; %invertible
B = [0;1]; %
E = B;
H = eye(2); %full state feedback

dim_n = size(B,1);
dim_m = size(B,2);
dim_p = size(H,1);
Ad = [A B;zeros(dim_m,dim_n) zeros(dim_m,dim_m)];

Hd = [H zeros(dim_p,dim_m)];

L1d=transpose(lqr(Ad',Hd',eye(3),eye(2))); %1st Luenberger: real eig(Ad-L1d*Hd)=-1.018 and -1.018 and -0.623
L2d=transpose(lqr(Ad',Hd',2*eye(3),eye(2))); %2nd Luenberger: real eig(Ad-L2d*Hd)=-1.395 and -1.395 and -0.768

F1d=Ad-L1d*Hd;
F2d=Ad-L2d*Hd;

%%%%
delta=0:0.01:100;

% H2d=@(x) eye(3)-expm(F2d*x)*expm(-F1d*x);
% 
% plot(delta,arrayfun(@(x) det(H2d(x)),delta))

%%%
rankM=@(x) inv(A)*(expm(A*x)-eye(2))*E;
plot(delta,arrayfun(@(x) rank(rankM(x)),delta))
