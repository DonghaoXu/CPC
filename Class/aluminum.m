classdef aluminum < Property
    properties (SetAccess=private)
        name='aluminum';
    end
    methods  
        function y=cp(~,T)
            % Temperature(K)  Specific heat(J/kg-K)
            % 200             790.5
            % 250             855.4
            % 298.15          897
            % 350             930.6
            % 400             955.5
            % 500             994.8
            % 600             1034
            y=1.199e-5*T.^3-1.516e-2*T.^2+6.512*T-6.285;    %J/kg-K
        end
        
%         function y=rho(~,T)
%             y=[];
%         end
        
%         function y=mu(~,T)
%             y=[];
%         end
        
        function y=k(~,T)
            % Temperature(K)  Thermal Conductivity(kW/m-K)
            % 250             0.235
            % 273             0.236
            % 300             0.237
            % 350             0.240
            % 400             0.240
            % 500             0.237
            % 600             0.232
            % 700             0.226
            y=-1.592e-4*T.^2+1.289e-1*T+213.1;  %W/m-K
        end
    end
end