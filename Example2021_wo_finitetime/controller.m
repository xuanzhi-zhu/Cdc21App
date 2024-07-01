classdef controller
    methods (Static)
        function dxc = dynamics()
            dxc=0;
        end
        function u = output(x)
            global Kd
            u = -Kd*x;
        end          
    end
end