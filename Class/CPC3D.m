classdef CPC3D
    properties
        AccHalfAng = pi/3;
        CR = 1.1748;
        RadReceiver = 0.0277;
        RadGlass = 0.0326;
        ThickGlass = 0.0016;
        Gap= 0.0048;
        RadUtube = 0.006;
        ThickReceiver = 0.002;
        ThickReflector = 0.002;
        L = 1.738;
        beta = deg2rad(25); 
        VacPres = 1e-5; %bar
    end
    properties (Dependent)
        Width;RCR;Height;
    end
    
    methods
        function r = CPC3D(varargin)
            if nargin == 0
                return
            end
            for i = 1: 2: length(varargin)
                if ismember(varargin{i}, properties(CPC3D))
                    r.(varargin{i}) = varargin{i + 1};
                else
                    fprintf('You cannot set %s.', varargin{i})
                end
            end
        end
        function self = set.RCR(self, rcr)
            self.RCR = rcr;
%             [w1,w2]=obj.minmaxwidth;
            self.CR= rcr * (w2 / w1 - 1) + 1;
        end
    end
    
    methods (Static)
%         function [f, g] = rcrfun(acc, cr, cr_max, rcr, ra, rg, gap, w, w_max, h, h_max)
%             % acc in RAD
%             w_max = CPC3D.wmaxfun(acc, ra, rg, gap);
%             cr_max = w_max / (2 * pi * ra);
%             rcr_cal = (cr - 1) / (cr_max - 1);
%             f = rcr - rcr_cal;
%         end
%         function f = wmaxfun(acc, cr, cr_max, rcr, ra, rg, gap, w, w_max, h, h_max)
%             f = (2 * sqrt((rg + gap).^2 - ra.^2) + (pi + 2 * asin(ra ./ (rg + gap))) * ra) ./ sin(acc) - w_max;
%             g = [-(2 * sqrt((rg + gap).^2 - ra.^2) + (pi + 2 * asin(ra ./ (rg + gap))) * ra) .* cos(acc) ./ (sin(acc)).^2;... % acc
%                 0;...%cr
%                 0;...%cr_max
%                 0;...%rcr
%                 -2 * ra ./ sqrt((rg + gap).^2 - ra.^2) + ];
%         end
%         function cpc
    end
end