classdef sensor
    methods (Static)
        function dhd = dynamics_flow(hd,tau,x,xs)
            dhd = 0;
        end
        function next_hd = dynamics_jump(hd,tau,x,xs)
            global A E
            invA = [0 -1;1 0];
            M = invA*(expm(A*tau)-eye(2))*E;
            next_hd = hd + (M'*M)\(M'*(x-xs));
        end
    end
end