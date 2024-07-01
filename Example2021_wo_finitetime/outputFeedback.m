function [t,j,xi] = outputFeedback()
    global delta sigma A B H K L1 L2 %(A,H) observable; F1 and F2 are Hurwitz; eye(2)-expm(F2*delta)*expm(-F1*delta) is invertible
    delta = 5; %desired convergence instant
    sigma = 2; %SOD of ETM on sensor-controller channel
    A = [0 1;0 0]; %plant
    B = [0;1];
    H = [1 0]; %full state feedback
%     Ob = obsv(A,H) 
%     unob = length(A)-rank(Ob) % Number of unobservable states
%     Co = ctrb(A,B)
%     unco = length(A) - rank(Co) % Number of uncontrollable states

    K = lqr(A,B,eye(2),1); %LQR on controller side

    %output feedback
    L1 = [1;1]; %1st Luenberger: real eig(A-L1*H)=-0.5 and -0.5
%     L2 = [2;1]; %2nd Luenberger: real eig(A-L2*H)=-1 and -1
    L2 = [-1;1]; %unstable: real eig(A-L2*H)=0.5 and 0.5
    % % check assumption eye(2)-expm(F2*delta)*expm(-F1*delta) is invertible
%     F1=A-L1*H;
%     F2=A-L2*H;
%     % % 
%     (eye(2)-expm(F2*delta)*expm(-F1*delta))^(-1)
    % 
    % L=transpose(lqr(A',H',eye(2),eye(2)));%optimal: minimum square estimator

    TSPAN = [0 250];
    JSPAN = [0 1000];
    rule  = 1;
    x_0 = rand(2,1)*10-5;
    hx1_0 = rand(2,1)*10-5;
    hx2_0 = rand(2,1)*10-5; %not necessary the same
    xc_0 = rand(2,1)*10-5;
%     hx2_0 = hx1_0;
%     xc_0 = hx1_0;
    tau_0 = 0;
    xi_0 = [x_0;hx1_0;hx2_0;xc_0;tau_0];
    [t,j,xi] = HyEQsolver(@F,@G,@C,@D,xi_0,TSPAN,JSPAN,rule);
    
    figure(1)
    subplot(2,2,1) %1st component of: state, obs1, obs2, pr
    plot(t,[xi(:,1) xi(:,3) xi(:,5) xi(:,7)])
    legend('x_1','hx1_1','hx2_1','x_c')
    subplot(2,2,2) %2nd component of: state, obs1, obs2, pr
    plot(t,[xi(:,2) xi(:,4) xi(:,6) xi(:,8)])
    legend('x_2','hx1_2','hx2_2','x_c')
    subplot(2,2,3)
    plot(t,xi(:,9))
    legend('timer')
    subplot(2,2,4)
    plot(t,j)
    legend('events')
    
    figure(2)
    out_sc = zeros(size(t)); 
    for i=1:1:length(t)
        hx1 = xi(i,3:4);
        hx2 = xi(i,5:6);
        xc = xi(i,7:8);
        out_sc(i) = ETM_sc.D(hx1',hx2',xc');
    end
    plot(t,out_sc)
    legend('channel events')
end    

function dxi = F(xi)
    x = xi(1:2);
    hx1 = xi(3:4);
    hx2 = xi(5:6);
    xc = xi(7:8);
    
    u = controller.output(xc); %use predictor state
%     u = 0; %null
    dx = plant.dynamics(x,u);
    y = plant.output(x);
    [dhx1,dhx2] = sensor.dynamics_flow(hx1,hx2,y,u);
    hxc = predictor.dynamics_flow(xc,u);
    dtau = 1;
    
    dxi = [dx;dhx1;dhx2;hxc;dtau];
end

function next_xi = G(xi)
    x = xi(1:2);
    hx1 = xi(3:4);
    hx2 = xi(5:6);
    xc = xi(7:8);
    tau = xi(9);
    
    out_s = sensor.D(tau); %observer event
    out_sc = ETM_sc.D(hx1,hx2,xc); %channel event
    if (out_s == 1) && (out_sc == 0)
        next_x = x;
        [next_hx1,next_hx2] = sensor.dynamics_jump(hx1,hx2);
        next_xc = xc;
        next_tau = 0; %if next_t=t, then simulation will stop at t=delta, here is a looping timer with interval of delta
    elseif (out_s == 0) && (out_sc == 1)
        next_x = x;
        next_hx1 = hx1;
        next_hx2 = hx2;
        next_xc = predictor.dynamics_jump(hx1,hx2);
        next_tau = tau;
    else
        next_x = x;
        
        index = round(rand(1)); %generate 0 or 1 randomly
        [next_hx1_p1,next_hx2_p1] = sensor.dynamics_jump(hx1,hx2);
        next_hx1_p2 = hx1;
        next_hx2_p2 = hx2;
        next_hx1_p12 = [next_hx1_p1 next_hx1_p2];
        next_hx2_p12 = [next_hx2_p1 next_hx2_p2];
        next_xc_p12 = [xc predictor.dynamics_jump(hx1,hx2)];
        next_tau_p12 = [0 tau];
        
        next_hx1 = next_hx1_p12(:,index+1); %choose either one
        next_hx2 = next_hx2_p12(:,index+1);
        next_xc = next_xc_p12(:,index+1);
        next_tau = next_tau_p12(:,index+1);
    end
    next_xi = [next_x;next_hx1;next_hx2;next_xc;next_tau];
end

function out = C(xi)
    hx1 = xi(3:4);
    hx2 = xi(5:6);
    xc = xi(7:8);
    tau = xi(9);
    
    out_s = sensor.C(tau);
    out_sc = ETM_sc.C(hx1,hx2,xc);
    
    out = out_s & out_sc;
end

function out = D(xi)
    hx1 = xi(3:4);
    hx2 = xi(5:6);
    xc = xi(7:8);
    tau = xi(9);
    
    out_s = sensor.D(tau);
    out_sc = ETM_sc.D(hx1,hx2,xc);
    
    out = out_s | out_sc;
end


