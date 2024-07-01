classdef ETM_sc
    methods (Static)
        function out = C(x,xs)
            global sigma
            if norm(x-xs)<=sigma
                out = 1;
            else
                out = 0;
            end
        end
        function out = D(x,xs)
            global sigma
            if norm(x-xs)>=sigma
                out = 1;
            else
                out = 0;
            end
        end
    end
end

