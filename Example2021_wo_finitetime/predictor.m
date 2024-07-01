classdef predictor
    methods (Static)
        function dxc = dynamics_flow(xc,u)
            global Ad Bd
            dxc = Ad*xc+Bd*u;
        end
        function next_xc = dynamics_jump(hx1)
            next_xc = hx1;
        end
    end
end