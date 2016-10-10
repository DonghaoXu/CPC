classdef copper < Property
    properties (SetAccess=private)
        name='copper';
    end
    methods 
%         function y=cp(~,T)
%             y=[];
%         end
        
%         function y=rho(~,T)
%             y=[];
%         end
        
%         function y=mu(~,T)
%             y=[];
%         end
        
        function y=k(~,T)
            % Temperature(K)  Thermal Conductivity(kW/m-K)
            % 200             0.413
            % 250             0.404
            % 273             0.401
            % 300             0.398
            % 350             0.394
            % 400             0.392
            % 500             0.388
            % 600             0.383
            % 700             0.377
            y=-5.834e-7*T.^3+8.746e-4*T.^2-4.679e-1*T+475.8;  %W/m-K
        end
    end
end