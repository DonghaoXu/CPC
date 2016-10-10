classdef Optical3D
    properties(SetAccess = private)
        incidentAngle; rays;
    end
    properties
        nRays; geometry;
    end
    
    methods
        function r = Optical3D(varargin)
            if nargin == 0
                r.nRays = 1000;
                r.rays = repmat(ray(), r.nRays, 1);
                r.geometry = CPC3D();
                [theta_t, theta_l] = meshgrid(0: 2: 89, 0: 2: 89);
                r.incidentAngle = [theta_t(:), theta_l(:)];
                return
            end
            for i = 1: 2: length(varargin)
                if ismember(varargin{i}, {'nRays', 'geometry'})
                    r.(varargin{i}) = varargin{i + 1};
                else
                    fprintf('You cannot set %s.', varargin{i})
                end
            end
        end
        function self = set.incidentAngle(self, ia)
            ia = double(ia);
            assert(size(ia, 2) == 2,...
                'Incident angle is specified by TRANSVERSAL and LONGITUDINAL angles.')
            self.incidentAngle = ia;
        end
    end
end