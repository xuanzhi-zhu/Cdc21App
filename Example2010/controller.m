classdef controller
    methods (Static)
        function dxs = dynamics(xs,hd)
            global A B K E
            dxs = (A-B*K)*xs + E*hd;
        end
        function u = output(x)
            global K
            u = -K*x;
        end          
    end
end