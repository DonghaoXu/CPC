classdef glass < Property
    properties (SetAccess=private)
        name='glass';
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
            % 250             0.00128
            % 273             0.00133
            % 300             0.00138
            % 350             0.00145
            % 400             0.00151
            % 500             0.00162
            % 600             0.00175
            % 700             0.00192
            y=1.344e-3*T+9.632e-1;  %W/m-K
        end
    end
end