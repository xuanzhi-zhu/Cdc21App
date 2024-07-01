classdef plant
    methods (Static)
        function dx = dynamics(x,u)
            global Ad Bd
            dx = Ad*x + Bd*u;
        end
        function y = output(x)
            global Hd
            y = Hd*x;
        end
    end
end

