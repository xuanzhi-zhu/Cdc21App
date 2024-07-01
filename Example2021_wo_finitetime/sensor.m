classdef sensor
    methods (Static)
        function dhx = dynamics_flow(hx1,y,u)
            global Ad Bd Hd L1d
            dhx = Ad*hx1+Bd*u+L1d*(y-Hd*hx1);
        end
    end
end