classdef predictor
    methods (Static)
        function dxc = dynamics_flow(xc,u)
            global A B
            dxc = A*xc+B*u;
        end
        function next_xc = dynamics_jump(hx1,hx2)
            next_xc = hx1;
        end
    end
end