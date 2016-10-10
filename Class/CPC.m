classdef CPC 
    properties
        AccHalfAng;CR;RadReceiver;RadGlass;ThickGlass;Gap;RadUtube;
        ThickReceiver=0.002;ThickReflector=0.002;L=1.5;beta;VacPres=1e-5;
    end
    properties (SetAccess=private)
        Width;RCR;Height;
    end
    methods % Constructor
        function r=CPC(varargin)
            % Accept inputs as RayTraceData/Struct class/Numeric vectors
            if nargin>0
                if isa(varargin{1},'RayTraceData')
                    r=varargin{1}.geometry;
                    return
                end
                if isa(varargin{1},'NumericThermal') || isa(varargin{1},'SimpleThermal')
                    r=varargin{1}.optical.geometry;
                    return
                end
                if isa(varargin{1},'CPC')
                    r=varargin{1};
                    return
                end
                if isa(varargin{1},'struct')
                    r.Width=varargin{1}.Width;
                    r.RadReceiver=varargin{1}.RadReceiver;
                    r.RadGlass=varargin{1}.RadGlass;
                    r.ThickGlass=varargin{1}.ThickGlass;
                    r.Gap=varargin{1}.Gap;
                    r.CR=varargin{1}.CR;
                    r.AccHalfAng=varargin{1}.AccHalfAng;
                    return
                end
                f=@(x) isa(x,'double');
                if all(cellfun(f,varargin))
                    in=[varargin(:);cell(12,1)];
                    in=in(1:12);
                    r.AccHalfAng=in{1};
                    r.CR=in{2};
                    r.RadReceiver=in{3};
                    r.RadGlass=in{4};
                    r.ThickGlass=in{5};
                    r.Gap=in{6};
                    r.RadUtube=in{7};
                    r.ThickReceiver=in{8};
                    r.ThickReflector=in{9};
                    r.L=in{10};
                    r.beta=in{11};
                    r.VacPres=in{12};
                    return
                end
                error('Expected input includes Class RayTraceData, Class NumericThermal, Class SimpleThermal, Struct history, or double')
            end 
        end
        
        function self = toy(self)
            self.AccHalfAng = randi(90);
            self.RadReceiver = 1;
            self.RadGlass = rand() * (2 - 1.1) + 1.1;
            self.ThickGlass = (self.RadGlass - self.RadReceiver) * (rand() * 0.6 + 0.2);
            self.Gap = rand() * 0.2;
            self = self.setRCR(rand());
        end
        
        function obj=setRCR(obj,rcr)
            obj.RCR=rcr;
            [w1,w2]=obj.minmaxwidth;
            obj.CR=rcr*(w2/w1-1)+1;
        end
        
        function r=check(obj)
            r=0;
            if isempty(obj.AccHalfAng) || ~isa(obj.AccHalfAng,'numeric')
                error('My angle! My angle! It''s a single positive number less than 90! .AccHalfAng=')
            end
            if isempty(obj.CR) || ~isa(obj.CR,'numeric')
                error('My concentration ratio! My concentration ratio! It''s a single positive number! .CR=')
            end
            if isempty(obj.RadReceiver) || ~isa(obj.RadReceiver,'numeric')
                error('My receiver radius! My receiver radius! It''s a single positive number! .RadReceiver=')
            end
            if isempty(obj.RadGlass) || ~isa(obj.RadGlass,'numeric')
                error('My glass radius! My glass radius! It''s a single positive number! .RadGlass=')
            end
            if isempty(obj.ThickGlass) || ~isa(obj.ThickGlass,'numeric')
                error('My glass thickness! My glass thickness! It''s a single positive number! .ThickGlass=')
            end
            if isempty(obj.Gap) || ~isa(obj.Gap,'numeric')
                error('My gap! My gap! It''s a single positive number! .Gap=')
            end
            [wmin,wmax]=obj.minmaxwidth;
            if obj.CR>wmax/wmin
                fprintf('Are you kidding me? The concentration ratio you gave is even greater than the maximum: %f',wmax/wmin)
            end
            if obj.ThickGlass>obj.RadGlass-obj.RadReceiver
                fprintf('Are you kidding me? The glass thickness you gave is even greater than space: %f',obj.RadGlass-obj.RadReceiver)
            end
        end
        
        function r=checkoptical(obj)
            r=obj.check;
        end
        
        function r=checknumeric(obj)
            r=obj.checkoptical;
            if isempty(obj.RadUtube) || ~isa(obj.RadUtube,'double')
                error('My u-tube radius! My u-tube radius! It''s a single positive number! .RadUtube=')
            end
            if isempty(obj.ThickReceiver) || ~isa(obj.ThickReceiver,'double')
                error('My receiver thickness! My receiver thickness! It''s a single positive number! .ThickReceiver=')
            end
            if isempty(obj.ThickReflector) || ~isa(obj.ThickReflector,'double')
                error('My reflector thickness! My reflector thickness! It''s a single positive number! .ThickReflector=')
            end
            if isempty(obj.L) || ~isa(obj.L,'double')
                error('My length! My length! It''s a single positive number! .L=')
            end
            if isempty(obj.beta) || ~isa(obj.beta,'double')
                error('My tilted angle! My tilted angle! It''s a single positive number! .beta=')
            end
            if isempty(obj.VacPres) || ~isa(obj.VacPres,'double')
                error('My vacuum pressure! My vacuum pressure! It''s a small number in bar. .VacPres=')
            end
            if obj.RadUtube>obj.RadReceiver/2
                error('Are you kidding me? The u-tube radius you gave is even greater than half receiver radius: %f',obj.RadReceiver/2)
            end
        end
        
        function r=checksimple(obj)
            r=obj.checkoptical;
            if isempty(obj.RadUtube) || ~isa(obj.RadUtube,'double')
                error('My u-tube radius! My u-tube radius! It''s a single positive number! .RadUtube=')
            end
            if isempty(obj.L) || ~isa(obj.L,'double')
                error('My length! My length! It''s a single positive number! .L=')
            end
            if isempty(obj.beta) || ~isa(obj.beta,'double')
                error('My tilted angle! My tilted angle! It''s a single positive number! .beta=')
            end
            if isempty(obj.VacPres) || ~isa(obj.VacPres,'double')
                error('My vacuum pressure! My vacuum pressure! It''s a small number in bar. .VacPres=')
            end
            if obj.RadUtube>obj.RadReceiver/2
                error('Are you kidding me? The u-tube radius you gave is even greater than half receiver radius: %f',obj.RadReceiver/2)
            end
        end
    end
    
    methods % Geometry calculation
        function obj=update(obj)
            obj.check;
            obj.Width=2*pi*obj.RadReceiver*obj.CR;
            [w_min,w_max]=obj.minmaxwidth;
            obj.RCR=(obj.Width-w_min)/(w_max-w_min);
            [yl1,yl2]=obj.profile; yl=[yl1;yl2];
            h=max(yl(:,2))-min(yl(:,2));
            obj.Height=h;
        end
        
        function [w_min,w_max]=minmaxwidth(obj)
            theta_deg=obj.AccHalfAng;
            theta=theta_deg*pi/180;
            r_a=obj.RadReceiver;
            r_g=obj.RadGlass;
            gap=obj.Gap;
            w_max=(2*sqrt((r_g+gap)^2-r_a^2)+(pi+2*asin(r_a/(r_g+gap)))*r_a)/sin(theta);
            w_min=2*pi*r_a;
        end
             
        function [spr,spl]=jointpoint(obj)
            % The joint point between parabolic and involution
            [~,spr,spl]=sprefl(obj.RadReceiver,...
                obj.RadGlass,obj.Gap,...
                obj.AccHalfAng*pi/180);
        end
        
        function [yl1,yl2,yr1,yr2,yc,yg_o,yg_i,ya]=profile(obj)
            [~,spl]=obj.jointpoint();
            [yl1,yl2,yr1,yr2,yc,yg_o,yg_i,ya]=profileAnaXCPC(...
                obj.AccHalfAng*pi/180,...
                obj.RadReceiver,...
                obj.CR*2*pi*obj.RadReceiver,...
                obj.RadGlass,...
                obj.ThickGlass,...
                obj.Gap,spl);
        end
    end
    
    methods % Visualization    
        function h=show(obj)
%             figure;
            [yl1,yl2,yr1,yr2,yc,yg_o,yg_i,~]=obj.profile();
            plot(yl1(:,1),yl1(:,2),'color',[38 188 213]/255)
            hold on
            axis equal
            plot(yl2(:,1),yl2(:,2),'color',[38 188 213]/255)
            plot(yr1(:,1),yr1(:,2),'color',[38 188 213]/255)
            plot(yr2(:,1),yr2(:,2),'color',[38 188 213]/255)
            plot(yc(:,1),yc(:,2),'color',[254 67 101]/255)
            plot(yg_o(:,1),yg_o(:,2),'color',[29 191 151]/255)
            plot(yg_i(:,1),yg_i(:,2),'color',[29 191 151]/255)
            if ~isempty(obj.RadUtube)
                a=0:pi/20:2*pi;
                xut=obj.RadReceiver-obj.RadUtube;
                plot(obj.RadUtube*cos(a)+xut,obj.RadUtube*sin(a),'color',[254 67 101]/255)
                plot(obj.RadUtube*cos(a)-xut,obj.RadUtube*sin(a),'color',[254 67 101]/255)
            end
            xl=xlim;
            yl=ylim;
            xlim(xl*1.05);
            ylim(yl*1.05);
            set(gcf,'color',[1 1 1])
            h=gcf;
        end
    end
end