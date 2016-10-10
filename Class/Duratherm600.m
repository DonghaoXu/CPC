classdef Duratherm600 < Property
    properties (SetAccess=private)
        name='Duratherm600';
    end
    methods 
        function y=cp(~,T)
            y=3.429*T+1084; %J/kg-K
        end
        
        function y=rho(~,T)
            y=-0.6801*T+1012.6; %kg/m3
        end
        
        function y=mu(~,T)
            % Temperature(K)  Dynamic viscosity(1e-3 kg/m-s)
            % 273             33.82
            % 283             20.04
            % 293             12.84
            % 303             8.760
            % 313             6.290
            % 323             4.770
            % 333             3.730
            % 343             3.000
            % 353             2.460
            % 363             2.060
            % 373             1.750
            % 383             1.510
            % 393             1.310
            % 403             1.160
            % 413             1.030
            % 423             0.920
            % 433             0.840
            % 443             0.760
            % 453             0.700
            % 463             0.640
            % 473             0.590
            % 483             0.550
            % 493             0.520
            % 503             0.480
            % 513             0.460
            % 523             0.430
            y=1.405e+2./(T-235.2).^2.294;   %kg/m-s
        end
        
        function y=k(~,T)
            y=-8.602e-5*T+0.1662;   %W/m-K
        end
        
    end
end