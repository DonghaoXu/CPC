classdef ray < handle
    properties
        lock;
    end
    properties(GetAccess = public, SetAccess = protected)
        x, y, index;
    end
    properties(SetAccess = immutable)
        maxn;
    end
    properties(SetAccess = private)
        cur = uint8(0);
    end
    properties(Dependent)
        ia;
    end
    
    methods
        function r = ray(varargin)
            r.maxn = uint8(10);
            r.x = zeros(r.maxn, 1);
            r.y = zeros(r.maxn, 1);
            r.index = cell(r.maxn, 1);
            r.lock = true;
            if nargin == 0
                return
            end
            for i = 1: 2: length(varargin)
                if ismember(varargin{i}, {'maxn', 'lock'})
                    r.(varargin{i}) = varargin{i + 1};
                else
                    fprintf('You cannot set %s.', varargin{i})
                end
            end
            r.lock = any(r.lock);
            r.maxn = uint8(r.maxn);
            assert(r.cur <= r.maxn,...
                'cur should be <= maxn.');
        end
        function newobj = copy(oldobj)
            newobj = ray('maxn', oldobj.maxn);
            p = properties(oldobj);
            for i = 1: numel(p)
                if ismember(p{i}, {'maxn', 'ia'})
                    continue
                end
                newobj.(p{i}) = oldobj.(p{i});
            end
        end
        function ia = get.ia(self)
            ia = zeros(self.cur, 1);
            p = [self.x, self.y];
            pvec = diff(p);
            pvec1 = pvec(1: end - 1, :);
            pvec2 = pvec(2: end, :);
            ia(2: end -1) = pi - acos(dot(pvec1, pvec2, 2)...
                ./ sqrt(sum(pvec1.^2, 2))...
                ./ sqrt(sum(pvec2.^2, 2)));
        end
    end
    
    methods
        function [exflag, msg] = add(self, x, y, index)
            % add x, y, ia, index to ray
            % exflag: -1 Success
            %          0 Instance is locked
            %          1 Error occurs
            %          2 Maximum number of incidence is reached
            if ~self.lock
                self.cur = self.cur + 1;
                if self.cur <= self.maxn
                    try
                        self.x(self.cur) = x;
                        self.y(self.cur) = y;
                        self.index{self.cur} = index;
                        exflag = -1;
                        msg = 'Success';
                    catch ex
                        exflag = 2;
                        msg = ex.message;
                        fprintf(msg)
                    end
                else
                    exflag = 1;
                    msg = 'Maximum number of incidence is reached';
                end
            else
                exflag = 0;
                msg = 'Instance is locked';
            end
        end
        function [exflag, msg] = drop(self)
            % drop last point
            % exflag: -1 Success
            %          0 Instance is locked
            %          1 Current cursor is at 0
            if self.lock
                exflag = 0;
                msg = 'Instance is locked';
                return
            end
            if self.cur > 0
                self.x(self.cur) = 0;
                self.y(self.cur) = 0;
                self.index{self.cur} = [];
                self.cur = self.cur - 1;
                exflag = -1;
                msg = 'Success';
            else
                exflag = 1;
                msg = 'Current cursor is at 0.';
            end
            
        end
    end
end