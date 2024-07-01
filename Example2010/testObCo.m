close all
clear all

A = [0 1;0 0];
B = [0;1];
C = [1 0;0 1];

A = [0 1;0 0];
B = [0;1];
C = [1 0];


dim_n = size(B,1);
dim_m = size(B,2);
dim_p = size(C,1);


Az = [A B;zeros(dim_m,dim_n) zeros(dim_m,dim_m)];
Bz = [B;zeros(dim_m,dim_m)];
Cz = [C zeros(dim_p,dim_m)];

ObAC = obsv(A,C) 
unobAC = length(A)-rank(ObAC) % Number of unobservable states
CoAB = ctrb(A,B)
uncoAB = length(A) - rank(CoAB) % Number of uncontrollable states

ObACz = obsv(Az,Cz) 
unobACz = length(Az)-rank(ObACz) % Number of unobservable states
CoABz = ctrb(Az,Bz)
uncoABz = length(Az) - rank(CoABz) % Number of uncontrollable states