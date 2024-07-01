function [td,jd,xid] = stateFeedback()
    rng('shuffle'); %reseed random number generator
    global delta sigma Ad Bd Hd Kd L1d L2d H1d H2d d
    delta = 1; %half of desired convergence instant
    sigma = 1; %SOD of ETM on sensor-controller channel
    A = [0 1;-1 0]; %invertible
    B = [0;1]; %
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

    Kd = [K eye(1)];
    
    L1d = transpose(lqr(Ad',Hd',eye(3),eye(2))); %1st Luenberger: real eig(Ad-L1d*Hd)=-1.018 and -1.018 and -0.623
    L2d = transpose(lqr(Ad',Hd',2*eye(3),eye(2))); %2nd Luenberger: real eig(Ad-L2d*Hd)=-1.395 and -1.395 and -0.768
    F1d=Ad-L1d*Hd;
    F2d=Ad-L2d*Hd;
%     eigs(F1d) %Re(F1d)<0
%     eigs(F2d) %Re(F2d)<0

    H1d = (eye(3)-expm(F1d*delta)*expm(-F2d*delta))^(-1);
    H2d = (eye(3)-expm(F2d*delta)*expm(-F1d*delta))^(-1);
%     det(H1d) %\neq 0
%     det(H2d) %\neq 0
    
    TSPAN = [0 10];
    JSPAN = [0 1000];
    rule  = 2;
    xd_0 = [rand(2,1)*10-5;d]; % initial of d has to be d
    xd_0 = [1;-1;d]; % initial of d has to be d
    hxd1_0 = zeros(3,1);
    hxd2_0 = zeros(3,1);
    xcd_0 = zeros(3,1);
    tau_0 = 0;
    
    xid_0 = [xd_0;hxd1_0;hxd2_0;xcd_0;tau_0];
    [td,jd,xid] = HyEQsolver(@F,@G,@C,@D,xid_0,TSPAN,JSPAN,rule);
    
    figure(1)
    subplot(2,2,1) %1st component of: state, obs1, obs2, pr
    plot(td,[xid(:,1) xid(:,4) xid(:,7) xid(:,10)])
    legend('x_1','hx1_1','hx2_1','x_c1')
    subplot(2,2,2) %2nd component of: state, obs1, obs2, pr
    plot(td,[xid(:,2) xid(:,5) xid(:,8) xid(:,11)])
    legend('x_2','hx1_2','hx2_2','x_c2')
    subplot(2,2,3) %disturbance: real, obs1, obs2, pr
    plot(td,[xid(:,3) xid(:,6) xid(:,9) xid(:,12)])
    legend('d','hd1','hd2','d_c')
    subplot(2,2,4)
    plot(td,jd)
    legend('events')
    
    figure(2)
    out_scd = zeros(size(td)); 
    for i=1:1:length(td)
        hxd1 = xid(i,4:6);
        hxd2 = xid(i,7:9);
        xcd = xid(i,10:12);
        out_scd(i) = ETM_sc.D(hxd1',hxd2',xcd');
    end
    plot(td,out_scd)
    legend('channel events')
    
    save exam2021.mat td jd xid out_scd
end    

function dxid = F(xid)
    xd = xid(1:3);
    hxd1 = xid(4:6);
    hxd2 = xid(7:9);
    xcd = xid(10:12);
    
    u = controller.output(xcd); %use predictor state
%     u = 0; %null
    dxd = plant.dynamics(xd,u);
    y = plant.output(xd);
    [dhxd1,dhxd2] = sensor.dynamics_flow(hxd1,hxd2,y,u);
    dxcd = predictor.dynamics_flow(xcd,u);
    dtau = 1;
    
    dxid = [dxd;dhxd1;dhxd2;dxcd;dtau];
end

function next_xi = G(xid)
    xd = xid(1:3);
    hxd1 = xid(4:6);
    hxd2 = xid(7:9);
    xcd = xid(10:12);
    tau = xid(13);
    
    out_s = sensor.D(tau); %observer event
    out_sc = ETM_sc.D(hxd1,hxd2,xcd); %channel event
    if (out_s == 1) && (out_sc == 0)
        next_xd = xd;
        [next_hxd1,next_hxd2] = sensor.dynamics_jump(hxd1,hxd2);
        next_xcd = xcd;
        next_tau = 0;
    elseif (out_s == 0) && (out_sc == 1)
        next_xd = xd;
        next_hxd1 = hxd1;
        next_hxd2 = hxd2;
        next_xcd = predictor.dynamics_jump(hxd1,hxd2);
        next_tau = tau;
    else
        next_xd = xd;
        
        index = round(rand(1)); %generate 0 or 1 randomly
        [next_hx1d_p1,next_hx2d_p1] = sensor.dynamics_jump(hxd1,hxd2);
        next_hx1d_p2 = hxd1;
        next_hx2d_p2 = hxd2;
        next_hx1d_p12 = [next_hx1d_p1 next_hx1d_p2];
        next_hx2d_p12 = [next_hx2d_p1 next_hx2d_p2];
        next_xcd_p12 = [xcd predictor.dynamics_jump(hxd1,hxd2)];
        next_tau_p12 = [0 tau];
        
        next_hxd1 = next_hx1d_p12(:,index+1); %choose either one
        next_hxd2 = next_hx2d_p12(:,index+1);
        next_xcd = next_xcd_p12(:,index+1);
        next_tau = next_tau_p12(:,index+1);
    end
    next_xi = [next_xd;next_hxd1;next_hxd2;next_xcd;next_tau];
end

function out = C(xid)
    xd = xid(1:3);
    hxd1 = xid(4:6);
    hxd2 = xid(7:9);
    xcd = xid(10:12);
    tau = xid(13);
    
    out_s = sensor.C(tau);
    out_sc = ETM_sc.C(hxd1,hxd2,xcd);
    
    out = out_s & out_sc;
end

function out = D(xid)
    xd = xid(1:3);
    hxd1 = xid(4:6);
    hxd2 = xid(7:9);
    xcd = xid(10:12);
    tau = xid(13);
    
    out_s = sensor.D(tau);
    out_sc = ETM_sc.D(hxd1,hxd2,xcd);
    
    out = out_s | out_sc;
end

