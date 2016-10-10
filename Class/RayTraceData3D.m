classdef RayTraceData3D
    
    properties (SetAccess=private)
        IncidentAngles;
    end
    properties
        NumberOfRays;geometry;
        x;y;z;index;ia;
    end
    
    methods % Constructor
        function r=RayTraceData3D(varargin)
            % Accept inputs as RayTraceData/Struct class/Numeric vectors
            if nargin>0
                if isa(varargin{1},'RayTraceData')
                    r=varargin{1};
                    return
                end
                if isa(varargin{1},'CPC')
                    r.geometry=varargin{1};
                    if nargin>1 && isa(varargin{2},'double')
                        r.NumberOfRays=varargin{2};
                    else
                        r.NumberOfRays=1000;
                    end
                    return
                end
                if isa(varargin{1},'NumericThermal') || isa(varargin{1},'SimpleThermal')
                    r=varargin{1}.optical;
                    return
                end
                if isa(varargin{1},'struct')
                    r.x=varargin{1}.x;
                    r.y=varargin{1}.y;
                    r.index=varargin{1}.index;
                    r.ia=varargin{1}.ia;
                    r.geometry=CPC(varargin{1});
                    r.NumberOfRays=size(varargin{1}.x,1);
                    return
                end
                f=@(x) isa(x,'double');
                if ~any(~cellfun(f,varargin))
                    in=[varargin(:);cell(13,1)];
                    in=in(1:8);
                    r.geometry=CPC(varargin);
                    r.NumberOfRays=in{13};
                    return
                end
                error('Expected input includes Class RayTraceData, Struct history, or double')
            end
        end
        
        function r=check(obj)
            r=obj.geometry.checkoptical;
            % check accordance
            if size(obj.x,1)~=obj.NumberOfRays
                fprintf('The number of rows in x coordinates doesn''t match the number of rays.');
            end
            if size(obj.y,1)~=obj.NumberOfRays
                fprintf('The number of rows in y coordinates doesn''t match the number of rays.');
            end
            if size(obj.index,1)~=obj.NumberOfRays
                fprintf('The number of rows in indices doesn''t match the number of rays.');
            end
            if size(obj.ia,1)~=obj.NumberOfRays
                fprintf('The number of rows in incident angles doesn''t match the number of rays.');
            end
            if size(obj.x,3)~=numel(obj.IncidentAngles)
                fprintf('The number of third dimension in x coordinates doesn''t match the number of rays.');
            end
            if size(obj.y,3)~=numel(obj.IncidentAngles)
                fprintf('The number of third dimension in y coordinates doesn''t match the number of rays.');
            end
            if size(obj.index,3)~=numel(obj.IncidentAngles)
                fprintf('The number of third dimension in indices doesn''t match the number of rays.');
            end
            if size(obj.ia,3)~=numel(obj.IncidentAngles)
                fprintf('The number of third dimension in incident angles doesn''t match the number of rays.');
            end
            % check settings
            if isempty(obj.NumberOfRays)
                fprintf('NumberOfRays is missing. .NumberOfRays=');
            end
        end
        
        function obj=setIA(obj, ia)
            obj.IncidentAngles = ia;
        end
    end
    
    methods % Simulation
        function obj=go(obj,varargin)
            % doesn't update obj
            % geometry:
            % 1. AccHalfAng
            % 2. CR
            % 3. RadReceiver
            % 4. RadGlass
            % 5. ThickGlass
            % 6. Gap
            % Simulation settings
            % 7. IncidentAngles
            % 8. NumberOfRays
            % full version: update(obj,iai,tol)
            obj.check;
            in=[varargin(:);cell(2,1)]; in=in(1:2);
            if isempty(in{2})
                tol=1e-5;
            else
                tol=varargin{2};
            end
            if isempty(in{1})
                iai=obj.IncidentAngles;
            else
                iai=in{1};
            end
            obj.IncidentAngles=iai;
            theta_deg=obj.geometry.AccHalfAng;
            theta=theta_deg*pi/180;
            r_a=obj.geometry.RadReceiver;
            r_g=obj.geometry.RadGlass;
            t_g=obj.geometry.ThickGlass;
            gap=obj.geometry.Gap;
            l=obj.NumberOfRays+1;
            [yl1,yl2,yr1,yr2,~,~,~,~]=obj.geometry.profile();
            history=RayTracing(yl2,yl1,yr2,yr1,r_g,t_g,r_a,gap,l,iai,tol,theta);
            obj.x=history.x;
            obj.y=history.y;
            obj.index=history.index;
            obj.ia=history.ia;
        end
        
        function obj = go3D(obj, sp, angles)
            % sp is a n x 3 matrix: (x, y, z)
            % angles is a n x 2 matrix: [transversal, longitudinal]
            assert(size(sp,1) == size(angles, 1))
            obj.IncidentAngles = angles;
            tol = 1e-5;
            theta_deg=obj.geometry.AccHalfAng;
            theta=theta_deg*pi/180;
            r_a=obj.geometry.RadReceiver;
            r_g=obj.geometry.RadGlass;
            t_g=obj.geometry.ThickGlass;
            gap=obj.geometry.Gap;
            [yl1,yl2,yr1,yr2,~,~,~,~]=obj.geometry.profile();
            nray=size(sp, 1);
            np = 10;
            
            xt = zeros(nray, np);
            yt = zeros(nray, np);
            zt = zeros(nray, np);
            it = zeros(nray, np);
            at = zeros(nray, np);
            for i = 1: nray
                oi = NaN;
                k = 1;
                agl = angles(i, :);
                theta_t = agl(i, 1);
                theta_l = agl(i, 2);
                vra = [tan(deg2rad(theta_t)), -1, tan(deg2rad(theta_l))];
                vra = normv(vra);
                in = sp(i, :);
                ia3D = NaN;
                
                tx = zeros(1, np);
                ty = zeros(1, np);
                tz = zeros(1, np);
                ti = zeros(1, np);
                tia = zeros(1, np);
                
                while oi~=0 && k<=np
                    tx(k) = in(1);
                    ty(k) = in(2);
                    tz(k) = in(3);
                    ti(k) = oi;
                    tia(k) = real(ia3D);
                    if k ~= np
                        [vra,in,ia3D,oi]=refl3D(vra,in,yl2,yl1,yr2,yr1,r_g,t_g,r_a,gap,tol,theta);
                    end
                    if in(3) > obj.geometry.L || in(3) < 0
                        break
                    end
                    k=k+1;
                end
                xt(i, :) = tx;
                yt(i, :) = ty;
                zt(i, :) = tz;
                it(i, :) = ti;
                at(i, :) = tia;
                
            end
            obj.x = xt;
            obj.y = yt;
            obj.z = zt;
            obj.index = it;
            obj.ia = at;
        end
        
        function obj=update(obj)
            obj.geometry=obj.geometry.update;
            obj=obj.go;
        end
        
        function sp = initrays3D(obj, angles, l)
            L = obj.geometry.Width * 100;
            [~, yl2] = obj.geometry.profile();
            left = min(yl2(:, 1));
            top = obj.geometry.RadGlass;
            sp_base = [linspace(left, -left, l)', repmat(top, l, 1), zeros(l, 1)];
            theta_t = angles(:, 1);
            theta_l = angles(:, 2);
            vra = [tan(deg2rad(theta_t)), -ones(size(theta_t)), tan(deg2rad(theta_l))];
            vra = normv(vra);
            sp = repmat(sp_base, 1, 1, size(angles, 1)) - ...
                repmat(reshape(L * vra', 1, 3, size(angles, 1)), l, 1);
            
        end
    end
    
    methods % Optical calculation
        function [q,mesh]=getefficiency(obj,type,varargin)
            % type: 'local' or 'total'
            % full version1: q=getefficiency(obj,'total',surfaceproperty,ctrl,angle)
            % q: 1.receiver 2.reflector 3.glass
            % full version2: getefficiency(obj,'local',surfaceproperty,ctrl,angle,mesh)
            % ctrl: ctrl(1) absorber ctrl(2) reflector ctrl(3) glass
            % ctrl: 0 for ignore the angular dependence, other for consider
            obj.check;
            in=[varargin(:);cell(4,1)]; in=in(1:4);
            if isempty(in{1})
                in{1}='SurfaceProperty.mat';
            end
            if isempty(in{2})
                in{2}=[1 1 1];
            end
            if isempty(in{3})
                a=obj.IncidentAngles;
                in{3}=[-a(end:-1:2);a];
            end
            if isempty(in{4})
                [q,mesh]=OpticalEff(obj,type,in{1},in{2},in{3});
            else
                [q,mesh]=OpticalEff(obj,type,in{1},in{2},in{3},in{4});
            end
        end
        
        function [intf,mesh]=getintercept(obj,type,varargin)
            % type: 'local' or 'total'
            % full version1: getintercept(obj,'local',angle,mesh)
            % full version2: getintercept(obj,'total',angle)
            obj.check;
            in=[varargin(:);cell(2,1)]; in=in(1:2);
            if isempty(in{1})
                a=obj.IncidentAngles;
                in{1}=[-a(end:-1:1);a];
            end
            if isempty(in{2})
                [intf,mesh]=InterFactor(obj,type,in{1});
            else
                [intf,mesh]=InterFactor(obj,type,in{1},in{2});
            end
        end
        
        function w_aper=aperture(obj,angle)
            obj.geometry=obj.geometry.update;
            r_a=obj.geometry.RadReceiver;
            w=obj.geometry.Width;
            h=obj.geometry.Height;
            w_abs=(r_a-h*sin((abs(angle(:)))/180*pi))./cos((abs(angle(:)))/180*pi)+w/2;
            w_aper=max(w,w_abs);
        end
        
        function [q,mesh]=getefficiency3D(obj,type,varargin)
            % type: 'local' or 'total'
            % full version1: q=getefficiency(obj,'total',surfaceproperty,ctrl,angle)
            % q: 1.receiver 2.reflector 3.glass
            % full version2: getefficiency(obj,'local',surfaceproperty,ctrl,angle,mesh)
            % ctrl: ctrl(1) absorber ctrl(2) reflector ctrl(3) glass
            % ctrl: 0 for ignore the angular dependence, other for consider
            obj.check;
            in=[varargin(:);cell(4,1)]; in=in(1:4);
            if isempty(in{1})
                in{1}='SurfaceProperty.mat';
            end
            if isempty(in{2})
                in{2}=[1 1 1];
            end
            if isempty(in{3})
                a=obj.IncidentAngles;
                in{3}=[-a(end:-1:2);a];
            end
            if isempty(in{4})
                [q,mesh]=OpticalEff3D(obj,type,in{1},in{2},in{3});
            else
                [q,mesh]=OpticalEff3D(obj,type,in{1},in{2},in{3},in{4});
            end
        end
    end
    
    methods % Visualization
        function [r,intf,mesh]=showintercept(obj,type,varargin)
            % type: 'local' or 'total'
            % r=showintercept(obj,'local',angle,intf,mesh)
            % r=showintercept(obj,'total',angle,intf)
            in=[varargin(:);cell(3,1)]; in=in(1:3);
            
            if strcmpi(type,'local')
                if isempty(in{1})
                    in{1}=obj.IncidentAngles;
                end
                if isempty(in{2}) || isempty(in{3})
                    fprintf('No data. No worry. I''ll get you covered.\n')
                    [intf,mesh]=obj.getintercept(type,in{1},in{3});
                elseif numel(in{1})~=size(in{2},2)
                    fprintf('Angle number is different from columns of intercept factor.\nNo worry. I''ll get you covered.\n')
                    [intf,mesh]=obj.getintercept(type,in{1},in{3});
                else
                    intf=in{2};
                    mesh=in{3};
                end
                tag=mesh.tag.tag;
                pmid=mesh.pmid;
                pmr=pmid; pmr(:,1)=pmid(:,2); pmr(:,2)=-pmid(:,1);
                loc=-atan2(pmr(:,2),pmr(:,1))*180/pi;
                id_abs=tag==2 | tag==-3 | tag==-4;
                id_glass=tag==1;
                intf_abs=intf(id_abs,:);
                intf_glass=intf(id_glass,:);
                ymax_abs=full(max(intf_abs(:))*1.1);
                ymax_glass=full(max(intf_glass(:))*1.1);
                loc_abs=loc(id_abs); [~,i_abs]=sort(loc_abs);
                loc_glass=loc(id_glass); [~,i_glass]=sort(loc_glass);
                figure, hold on
                ax_abs=subplot(1,2,1);
                ax_glass=subplot(1,2,2);
                f=gcf;
                f.Position=[200,250,1000,400];
                for i=1:numel(in{1})
                    habs=plot(ax_abs,loc_abs(i_abs),intf_abs(i_abs,i));
                    hglass=plot(ax_glass,loc_glass(i_glass),intf_glass(i_glass,i));
                    ax_abs.YLim=[0,ymax_abs];
                    ax_abs.XLim=[-180,180];
                    ax_glass.YLim=[0,ymax_glass];
                    ax_glass.XLim=[-180,180];
                    xlabel(ax_abs,'Absorber local angle ({\circ})')
                    ylabel(ax_abs,'Intercept factor')
                    xlabel(ax_glass,'Glass local angle ({\circ})')
                    ylabel(ax_glass,'Intercept factor')
                    title(ax_abs,sprintf('Incidence angle = %.1f Deg.',in{1}(i)),'fontweight','normal')
                    title(ax_glass,sprintf('Incidence angle = %.1f Deg.',in{1}(i)),'fontweight','normal')
                    FigureFormat(12,4,8);
                    pause(0.01)
                    if i==numel(in{1})
                        break
                    end
                    delete(habs)
                    delete(hglass)
                end
            else
                if isempty(in{1})
                    in{1}=obj.IncidentAngles;
                end
                if isempty(in{2})
                    fprintf('No data. No worry. I''ll get you covered.\n')
                    [intf,mesh]=obj.getintercept(type,in{1});
                elseif numel(in{1})~=size(in{2},2)
                    fprintf('Angle number is different from columns of intercept factor.\nNo worry. I''ll get you covered.\n')
                    [intf,mesh]=obj.getintercept(type,in{1});
                else
                    intf=in{2};
                    mesh='hehe';
                end
                figure
                angle=in{1}; [angle,i]=sort(angle);
                plot(angle,intf(:,i))
                xlabel('Incident angle ({\circ})')
                ylabel('Intercept factor')
                xlim([0,90])
                legend('Absorber','Reflector','Glass')
                f=gcf;
                f.Position=[400,250,500,400];
                FigureFormat(6,4,8)
            end
            r=gcf;
        end
        
        function [r,q,mesh]=showefficiency(obj,type,varargin)
            % type: 'local' or 'total'
            % r=showintercept(obj,'local',angle,q,mesh)
            % r=showintercept(obj,'total',angle,q)
            in=[varargin(:);cell(3,1)]; in=in(1:3);
            
            if strcmpi(type,'local')
                if isempty(in{1})
                    in{1}=obj.IncidentAngles;
                end
                if isempty(in{2}) || isempty(in{3})
                    fprintf('No data. No worry. I''ll get you covered.\n')
                    [q,mesh]=obj.getefficiency(type,'SurfaceProperty.mat',[1 1 1],in{1},in{3});
                elseif numel(in{1})~=size(in{2},2)
                    fprintf('Angle number is different from columns of intercept factor.\nNo worry. I''ll get you covered.\n')
                    [q,mesh]=obj.getefficiency(type,'SurfaceProperty.mat',[1 1 1],in{1},in{3});
                else
                    q=in{2};
                    mesh=in{3};
                end
                tag=mesh.tag.tag;
                pmid=mesh.pmid;
                pmr=pmid; pmr(:,1)=pmid(:,2); pmr(:,2)=-pmid(:,1);
                loc=-atan2(pmr(:,2),pmr(:,1))*180/pi;
                id_abs=tag==2 | tag==-3 | tag==-4;
                id_glass=tag==1;
                q_abs=q(id_abs,:);
                q_glass=q(id_glass,:);
                ymax_abs=full(max(q_abs(:))*1.1);
                ymax_glass=full(max(q_glass(:))*1.1);
                loc_abs=loc(id_abs); [~,i_abs]=sort(loc_abs);
                loc_glass=loc(id_glass); [~,i_glass]=sort(loc_glass);
                figure, hold on
                ax_abs=subplot(1,2,1);
                ax_glass=subplot(1,2,2);
                f=gcf;
                f.Position=[200,250,1000,400];
                for i=1:numel(in{1})
                    habs=plot(ax_abs,loc_abs(i_abs),q_abs(i_abs,i));
                    hglass=plot(ax_glass,loc_glass(i_glass),q_glass(i_glass,i));
                    ax_abs.YLim=[0,ymax_abs];
                    ax_abs.XLim=[-180,180];
                    ax_glass.YLim=[0,ymax_glass];
                    ax_glass.XLim=[-180,180];
                    xlabel(ax_abs,'Absorber local angle ({\circ})')
                    ylabel(ax_abs,'Efficiency')
                    xlabel(ax_glass,'Glass local angle ({\circ})')
                    ylabel(ax_glass,'Efficiency')
                    title(ax_abs,sprintf('Incidence angle = %.1f Deg.',in{1}(i)),'fontweight','normal')
                    title(ax_glass,sprintf('Incidence angle = %.1f Deg.',in{1}(i)),'fontweight','normal')
                    FigureFormat(12,4,8);
                    pause(0.01)
                    if i==numel(in{1})
                        break
                    end
                    delete(habs)
                    delete(hglass)
                end
            else
                if isempty(in{1})
                    in{1}=obj.IncidentAngles;
                end
                if isempty(in{2})
                    fprintf('No data. No worry. I''ll get you covered.\n')
                    [q,mesh]=obj.getefficiency(type,'SurfaceProperty.mat',[1 1 1],in{1});
                elseif numel(in{1})~=size(in{2},2)
                    fprintf('Angle number is different from columns of intercept factor.\nNo worry. I''ll get you covered.\n')
                    [q,mesh]=obj.getefficiency(type,'SurfaceProperty.mat',[1 1 1],in{1});
                else
                    q=in{2};
                    mesh='hehe';
                end
                figure
                angle=in{1}; [angle,i]=sort(angle);
                plot(angle,q(:,i))
                xlabel('Incident angle ({\circ})')
                ylabel('Efficiency')
                xlim([0,90])
                ylim([0,1])
                legend('Absorber','Reflector','Glass')
                f=gcf;
                f.Position=[400,250,500,400];
                FigureFormat(6,4,8)
            end
            r=gcf;
        end
        
        function h=showraytrace(obj,varargin)
            % varargin
            % 1. Incident angles to display
            % 2. Density, 0~1
            % 3. Display control
            %       0: Animate and delete all incident angles
            %       1: Wait for keyboard
            %       other: Keep all the rays
            % Default:	varargin{1}: no rays
            %           varargin{2}: min(1,100/obj.NumberOfRays)
            %           varargin{3}: 0
            % full version: show(obj,angles,density,dctrl)
            obj.geometry.show;
            if ~isempty(varargin)
                in=[varargin(:);cell(3,1)];
                in=in(1:3);
                if ~isempty(in{1})
                    % Generate the incident angle arrays
                    iav=in{1};
                    % Generate the indices of rays to be displayed
                    den=in{2};
                    if isempty(den)
                        den=min(1,100/obj.NumberOfRays);
                    end
                    totalray=den*obj.NumberOfRays;
                    rays=linspace(1,obj.NumberOfRays,totalray);
                    rays=unique(round(rays));
                    % Generate the control
                    dctrl=in{3};
                    if isempty(dctrl)
                        dctrl=0;
                    end
                    % Make figures
                    for i=1:numel(iav)
                        inang=iav(i);
                        i_ia=find(abs(abs(inang)-obj.IncidentAngles)<1e-6);
                        if numel(i_ia)==1
                            s=((inang>0)-0.5)*2;
                            xray=obj.x(rays,:,i_ia)*s;
                            yray=obj.y(rays,:,i_ia);
                            rh=plot(xray',yray','k');
                            title(sprintf('Incidence angle = %.1f Deg.', inang))
                            pause(0.01)
                            if dctrl==1
                                keyboard
                                delete(rh)
                            end
                            if dctrl==0
                                delete(rh)
                            end
                        else
                            fprintf('The ray tracing data doesn''t contain the specified incident angle: %.2f\n',inang)
                        end
                    end
                end
            end
            h=gcf;
        end
    end
    
    methods % not general
        function obj = fullangle(obj)
            [theta_t, theta_l] = meshgrid(linspace(0, 89, 41), linspace(0, 89, 41));
            angles = [theta_t(:), theta_l(:)];
            l = 500;
            sp = obj.initrays3D(angles, l);
            sp = reshape(permute(sp, [1 3 2]), [], 3); 
            % sp is now (l*size(angles, 1)) x 3
            angles = reshape(repmat(angles', l, 1), 2, [])';
            obj = obj.go3D(sp, angles);
        end
    end
    
%     methods
%         function history = go_big(obj, l)
%             tol = 1e-5;
%             iai = rand(l, 1) * 89;
%             theta_deg = obj.geometry.AccHalfAng;
%             theta = deg2rad(theta_deg);
%             r_a = obj.geometry.RadReceiver;
%             r_g = obj.geometry.RadGlass;
%             t_g = obj.geometry.ThickGlass;
%             gap = obj.geometry.Gap;
%             [yl1, yl2, yr1, yr2, ~, ~, ~, ~] = obj.geometry.profile();
%             history = RayTraceData.RayTracing_big(yl2, yl1, yr2, yr1, r_g, t_g, r_a, gap, l, iai, tol, theta);
%         end
%     end
%     
%     methods (Static, Access = private)
%         function history = RayTracing_big(yl2, yl1, yr2, yr1, r_g, t_g, r_a, gap, l, ias, tol, theta)
%             assert(l == numel(ias))
%             ias = deg2rad(ias(:));
%             nray = 10;
%             top = yl1(1, :);
%             if ~isnan(yl2(1, 1))
%                 top = yl2(end, :);
%             end
%             pleft = top;
%             L = 100;
%             hit = zeros(l, 1);
%             sp = zeros(l, 2);
%             
%             parfor i = 1: l
%                 lambda = rand();
%                 iai = ias(i);
%                 v1 = r_a * [cos(iai), sin(iai)];
%                 % v(2) = top(2). So ll should be as below
%                 ll = (top(2) - r_a * sin(iai)) / sin(-pi / 2 + iai);
%                 v2 = ll * [cos(-pi / 2 + iai), sin(-pi / 2 + iai)];
%                 v = v1 + v2;
%                 if v(1) > -top(1)
%                     pright = v;
%                 else
%                     pright = [-top(1), top(2)];
%                 end
%                 spi = pleft * lambda + pright * (1 - lambda) - v2 / ll * L;
%                 sp(i, :) = spi;
%                 fprintf('Incident angle: %.2f, Coordinates: [%.2f, %.2f]\n',...
%                     rad2deg(iai), spi(1), spi(2))
%                 % Start ray tracing
%                 k = 1;
%                 vra= v2 / ll;
%                 in = sp(i, :);
%                 while k<=nray
%                     [vra, in, ~, oi]=refl(vra,in,yl2,yl1,yr2,yr1,r_g,t_g,r_a,gap,tol,theta);
%                     if oi == -1
%                         hit(i) = 1;
%                         break
%                     elseif oi == 0
%                         break
%                     end
%                     k = k + 1;
%                 end
%             end
%             
%             history.hit = hit;
%             history.sp = sp;
%             history.ia = ias;
%         end
%     end
end

