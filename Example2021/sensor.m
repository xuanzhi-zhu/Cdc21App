classdef sensor
    methods (Static)
        function [dhx1,dhx2] = dynamics_flow(hx1,hx2,y,u)
            global Ad Bd Hd L1d L2d
            dhx1 = Ad*hx1+Bd*u+L1d*(y-Hd*hx1);
            dhx2 = Ad*hx2+Bd*u+L2d*(y-Hd*hx2);
        end
        function [next_hx1,next_hx2] = dynamics_jump(hx1,hx2)
            global H1d H2d
            next_hx1 = H1d*hx1+H2d*hx2;
            next_hx2 = H1d*hx1+H2d*hx2;
        end
        function out = C(tau)
            global delta
            if tau <= delta
                out = 1;
            else
                out = 0;
            end
        end
        function out = D(tau)
            global delta
            if tau >= delta
                out = 1;
            else
                out = 0;
            end
        end
    end
end