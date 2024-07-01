function [t,j,xi] = main_multijump()
    rng('shuffle'); %reseed random number generator
    global delta sigma A B E H K L1d L2d d
    delta = 1; %half of desired convergence instant
    sigma = 1; %SOD of ETM on sensor-controller channel
    A = [0 1;-1 0]; %invertible
    B = [0;1]; %
    E = B;
    H = eye(2); %full state feedback
    
    d = 1; %constant disturbance
    
%     rank(inv(A)*(expm(A*2*pi)-eye(2))*E) %full rank = min{size(B)} = 1
    
    dim_n = size(B,1);
    dim_m = size(B,2);
    dim_p = size(H,1);
    Ad = [A B;zeros(dim_m,dim_n) zeros(dim_m,dim_m)];
    Bd = [B;zeros(dim_m,dim_m)];
    Hd = [H zeros(dim_p,dim_m)];
    % ObAH = obsv(A,H) 
    % unobAH = length(A)-rank(ObAH) % Number of unobservable states
    % CoAB = ctrb(A,B)
    % uncoAB = length(A) - rank(CoAB) % Number of uncontrollable states
%     ObAHd = obsv(Ad,Hd) 
%     unobAHd = length(Ad)-rank(ObAHd) % Number of unobservable states
    % CoABd = ctrb(Ad,Bd)
    % uncoABd = length(Ad) - rank(CoABd) % Number of uncontrollable states


    K = lqr(A,B,eye(2),1); %LQR on controller side
%     F = A-B*K;
%     eigs(F) %Re(F)<0
%     Kd = [K eye(1)];
    
    L1d=transpose(lqr(Ad',Hd',eye(3),eye(2))); %1st Luenberger: real eig(Ad-L1d*Hd)=-1.018 and -1.018 and -0.623
    L2d=transpose(lqr(Ad',Hd',2*eye(3),eye(2))); %2nd Luenberger: real eig(Ad-L2d*Hd)=-1.395 and -1.395 and -0.768
%     F1d=Ad-L1d*Hd;
%     F2d=Ad-L2d*Hd;
%     eigs(F1d) %Re(F1d)<0
%     eigs(F2d) %Re(F2d)<0

%     H1d=(eye(3)-expm(F1d*delta)*expm(-F2d*delta))^(-1);
%     H2d=(eye(3)-expm(F2d*delta)*expm(-F1d*delta))^(-1);
%     det(H1d) %\neq 0
%     det(H2d) %\neq 0
    
    TSPAN = [0 10];
    JSPAN = [0 1000];
    rule  = 2;
    x_0 = rand(2,1)*10-5;
    x_0 = [1;-1];
    xs_0 = [0;0];
    tau_0 = 0;
    hd_0 = 0;
    q_0 = 0;
    
    xi_0 = [x_0;xs_0;tau_0;hd_0;q_0];
    [t,j,xi] = HyEQsolver(@F,@G,@C,@D,xi_0,TSPAN,JSPAN,rule);
    
    figure(1)
    subplot(2,2,1) %1st component of: x, xs
    plot(t,[xi(:,1) xi(:,3)])
    legend('x_1','xs_1')
    subplot(2,2,2) %2nd component of: x, xs
    plot(t,[xi(:,2) xi(:,4)])
    legend('x_2','xs_2')
    subplot(2,2,3)
    plot(t,xi(:,6))
    legend('hd')
    subplot(2,2,4)
    plot(t,j)
    legend('events')
   
    figure(2)
    out_sc = zeros(size(t)); 
    for i=1:1:length(t)
        x = xi(i,1:2);
        xs = xi(i,3:4);
        out_sc(i) = ETM_sc.D(x',xs');
    end
    plot(t,out_sc)
    legend('channel events')
    
    save exam2010.mat t j xi out_sc
end    

function dxi = F(xi)
    x = xi(1:2);
    xs = xi(3:4);
    tau =  xi(5);
    hd = xi(6);
    q = xi(7);
    
    u = controller.output(xs); %use predictor state
%     u = 0; %null
    dx = plant.dynamics(x,u);
    y = plant.output(x);
    dxs = controller.dynamics(xs,hd);
    dtau = 1;
    dhd = sensor.dynamics_flow(hd,tau,x,xs);
    dq = 0;
    
    dxi = [dx;dxs;dtau;dhd;dq];
end

function next_xi = G(xi)
    x = xi(1:2);
    xs = xi(3:4);
    tau =  xi(5);
    hd = xi(6);
    q = xi(7);
    
    out_i = (q==0); %intial time update xs to x
    out_sc = ETM_sc.D(x,xs); %trigger
    
    next_x = x;
    next_xs = x;
    next_q = 1;

    if (out_i == 1) && (out_sc == 0)
        next_tau = tau;
        next_hd = hd;
    elseif (out_i == 0) && (out_sc == 1)
        next_tau = 0;
        next_hd = sensor.dynamics_jump(hd,tau,x,xs);
    else
        next_tau = tau;
        next_hd = hd;
    end

    next_xi = [next_x;next_xs;next_tau;next_hd;next_q];
end

function out = C(xi)
    x = xi(1:2);
    xs = xi(3:4);
    tau =  xi(5);
    hd = xi(6);
    q = xi(7);
    
    out_i = (q~=0);
    out_sc = ETM_sc.C(x,xs); %channel event
    
    out = out_i & out_sc;
end

function out = D(xi)
    x = xi(1:2);
    xs = xi(3:4);
    tau =  xi(5);
    hd = xi(6);
    q = xi(7);
    
    out_i = (q==0);
    out_sc = ETM_sc.D(x,xs); %channel event
   
    out = out_i | out_sc;
end

