classdef ETM_sc
    methods (Static)
        function out = C(hx1,xc)
            global sigma
            if norm(xc-hx1)<=sigma
                out = 1;
            else
                out = 0;
            end
        end
        function out = D(hx1,xc)
            global sigma
            if norm(xc-hx1)>=sigma
                out = 1;
            else
                out = 0;
            end
        end
    end
end

