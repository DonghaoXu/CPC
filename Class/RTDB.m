classdef RTDB < handle
    properties
        db;
    end
    properties (Constant)
        colname = {'acc', 'RCR', 'rrg', 'rgap', 'rtg',...
            'interfact', 'opteff'};
    end
    properties (Access = private)
        colnames = {'Acceptance half angle',...
            'RCR',...
            'Relative glass radius',...
            'Relative gap',...
            'Relative glass thickness',...
            'Intercept factor',...
            'Optical efficiency'};
        units = {'{\circ}', '', '', '', '', '', ''};
    end
    
    methods
        function self = RTDB(celldata)
            self.db = table([], [], [], [], [], [], [],...
                'variablenames', self.colname);
            if nargin > 0
                self.insert(celldata);
            end
        end
        
        function insert(self, celldata)
            dbnew = cell2table(celldata,...
                'variablenames', self.colname);
            self.db = [self.db; dbnew];
        end
        
        function delete(self, val)
            % row number provided
            if isa(val, 'numeric')
                self.db(val, :) = [];
                return
            end
            % bool array provided
            if isa(val, 'logical')
                self.db(val, :) = [];
                return
            end
        end
        
        function disp(self)
            self.db
        end
        
        function m = length(self)
            m = size(self.db, 1);
        end
        
        function flag = is_valid(self)
            % complex number checking
            flag = isreal(table2array(self.db(:, 1:5)));
            flag = flag && isreal(cell2mat(table2array(self.db(:, 6: 7))));
            if ~flag
                fprintf('Complex number check fail')
                return
            end
            
            % range check
            flag = flag && all(self.db.acc <= 90) && all(self.db.acc > 0);
            flag = flag && all(self.db.RCR <= 1) && all(self.db.RCR >=0);
            flag = flag && all(self.db.rrg > 1);
            flag = flag && all(self.db.rgap >= 0);
            flag = flag && all(self.db.rtg > 0) && all(self.db.rtg < 1);
            if ~flag
                fprintf('Range check fail for configurations')
                return
            end
            flag = flag && all(all(cell2mat(self.db.interfact) >= 0));
            flag = flag && all(all(cell2mat(self.db.opteff) <= 1)) ...
                && all(all(cell2mat(self.db.opteff) >= 0));
            if ~flag
                fprintf('Range check fail for interfact or opteff')
                return
            end
            
        end
    end
    
    methods
        function sortby(self, name)
            if isa(name, 'cell')
                assert(all(cellfun(@(x) isa(x, 'char'), name)), ...
                    'Input should be a cell of name strings')
            elseif isa(name, 'numeric')
                name = self.colname(name);
            end
            self.db = sortrows(self.db, name);
        end
        
        function unique(self)
            [~, i] = unique(self.db(:, 1: 5));
            self.db = self.db(i, :);
        end
        
    end
    
    methods (Access = private, Static)
        function index = encode(num_array)
            [m, n] = size(num_array);
            coded_array = num_array .* repmat(10.^(2 * (n - 1): -2: 0), m, 1);
            index = sum(coded_array, 2);
        end
    end
end