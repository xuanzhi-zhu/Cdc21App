classdef plant
    methods (Static)
        function dx = dynamics(x,u)
            global A B E d
            dx = A*x + B*u +E*d;
        end
        function y = output(x)
            global H
            y = H*x;
        end
    end
end

