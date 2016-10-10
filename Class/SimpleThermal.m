classdef SimpleThermal
    properties
        optical;outdoor;operation;
        HTF=Duratherm600();
        Reflector=aluminum();
        Receiver=copper();
        Glass=glass();
        rho_grd=0.2;
    end
    properties (SetAccess=protected)
        Area;VFA;solarflux;
        temperature;Ambient=air();
        efficiency;
    end
    properties (Access=private)
        other;
    end
    methods % Constructor
        function r=SimpleThermal(varargin)
            if nargin>0
                if isa(varargin{1},'OptimalDesign')
                    x=varargin{1};
                    r.optical=x.optical;
                    r.outdoor=x.outdoor;
                    r.operation=x.operation;
                    r.Area=x.Area;
                    r.VFA=x.VFA;
                    r.solarflux=x.solarflux;
                    r.efficiency=x.efficiency;
                    return
                end
                if isa(varargin{1},'NumericThermal')
                    r=varargin{1};
                    return
                end
                if isa(varargin{1},'RayTraceData')
                    r.optical=varargin{1};
                    return
                end
                if isa(varargin{1},'CPC')
                    r.optical=RayTraceData(varargin{1});
                    return
                end
                if isa(varargin{1},'struct')
                    r.optical=RayTraceData(varargin{1});
                    return
                end
                if isa(varargin{1},'numeric') % Redefine
                    r.optical.geometry=CPCgeo(varargin{1});
                    return
                end
            end
        end
        
        function obj=set.outdoor(obj,weather)
            if isa(weather,'tmy3')
                obj.outdoor=weather;
            else
                error('The .outdoor should be a tmy3 class.')
            end
        end
        
        function obj=set.operation(obj,op)
            if ~isa(op,'struct')
                error('Operation is a struct with fields ''T_in'' and ''MassFR''.')
            elseif numel(op)~=1
                error('I can only take one operation condition at a time.')
            else
                op.T_in;op.MassFR;
                obj.operation=op;
            end
        end
               
        function obj=set.rho_grd(obj,r)
            if r<=0 || r>=1
                error('Did you ever see a ground reflectance being negative or greater than 1? Hehe')
            else
                obj.rho_grd=r;
            end
        end
        
        function obj=set.Reflector(obj,r)
            if isa(r,'Property')
                obj.Reflector=r;
            else
                error('Property! Property Class!')
            end
        end
        
        function obj=set.Receiver(obj,r)
            if isa(r,'Property')
                obj.Receiver=r;
            else
                error('Property! Property Class!')
            end
        end
        
        function obj=set.Glass(obj,r)
            if isa(r,'Property')
                obj.Glass=r;
            else
                error('Property! Property Class!')
            end
        end
        
        function obj=set.HTF(obj,r)
            if isa(r,'Property')
                obj.HTF=r;
            else
                error('Property! Property Class!')
            end
        end
        
        function obj=setsf(obj,q_sf)
            obj.solarflux=q_sf;
        end
        
        function r=getsf(obj)
            r=obj.solarflux;
        end
        
        function obj=setambient(obj,amb)
            obj.Ambient=amb;
        end
        
        function r=check(obj)
            obj.optical.check;
            r=obj.optical.geometry.checksimple;
            if isempty(obj.optical.geometry.Height) || isempty(obj.optical.geometry.Width)
                error('Wow, wow, wow. Hold on, pal. Height or Width are missing. Update it first.')
            end
            if isempty(obj.rho_grd)
                error('Wow, wow, wow. Hold on, pal. Set the ground reflectivity first: .rho_grd=')
            end
            if isempty(obj.outdoor)
                fprintf('Just a suggestion. Set the outdoor condition: .outdoor=')
            end
        end
    end
    
    methods % Simulation
        function [r,obj]=go(obj,varargin)
            % go(obj): simulate the conditions of obj.solarflux
            % go(obj,angle,beam,diffuse,T_amb,Vel,op): simulate the condition where
            % incident angle, beam and diffuse are as specified
            in=[varargin(:);cell(6,1)]; in=in(1:6);
            obj.check;
            switch nargin
                case 1
                    if isempty(obj.solarflux)
                        [~,ang]=obj.outdoor.angles(obj.optical.geometry.beta);
                        [obj.solarflux,Beam,Diffuse]=obj.getsolarflux;
                    else
                        ang=1:size(obj.solarflux,2);
                        Beam=sum(obj.solarflux)*2;
                        Diffuse=0;
                    end
                    T_amb=obj.outdoor.DryBulb;
                    op=obj.operation;
                otherwise
                    try
                        horzcat(in{1:5});
                        vertcat(in{1:5});
                    catch
                        error('numel(Diffuse)~=numel(Beam) || numel(Diffuse)~=numel(angle) || numel(Diffuse)~=numel(T_amb) || numel(Diffuse)~=numel(Vel). I''m tired of talking.')
                    end
                    obj.other.outdoor=obj.outdoor;
                    obj.outdoor=tmy3();
                    ang=in{1};
                    Beam=in{2};
                    Diffuse=in{3};
                    T_amb=in{4};
                    Vel=in{5};
                    op=in{6};
                    obj.outdoor=obj.outdoor.dummy(numel(ang));
                    obj.outdoor.Wspd=Vel(:);
                    obj.outdoor.DryBulb=T_amb(:);
                    obj.solarflux=obj.getsolarflux([],ang,Beam,Diffuse);
            end
            if isempty(obj.VFA)
                obj=obj.vfa;
            end
            obj.Ambient=air();
            obj.Ambient=obj.Ambient.setdefault(obj.outdoor.DryBulb);
            G=Beam+Diffuse;
            T_result=cell(1,numel(ang));
            options = optimoptions('fsolve','Jacobian','on','Display','iter-detailed','MaxIter',40);
            dt=[5,50,500,5000];
            r_g=obj.optical.geometry.RadGlass;
            t_g=obj.optical.geometry.ThickGlass;
            L=obj.optical.geometry.L;
            r_a=obj.optical.geometry.RadReceiver;
            t_a=r_a/14;
            r_u=obj.optical.geometry.RadUtube;
            for j=1:numel(ang)
                obj.operation=op(j);
                T0=obj.initemp(G(j),T_amb(j));
                [T,~,exitflag]=fsolve(@(T) obj.heatflux(T,j),T0,options);
                i=0;
                if exitflag~=1
                    T=T0;
                    while exitflag~=1 && i<numel(dt)
                        i=i+1;
                        fprintf('Constructing a new guess value...\n')
                        m=zeros(size(T));
                        c=zeros(size(T));
                        m(1)=2*pi*r_g*t_g*L*2230; % borosilicate glass: 2230kg/m3
                        m(3)=2*pi*r_a*t_a*L*8960; % copper: 8960kg/m3
                        m(2)=obj.Area.A_reflector*t_a*L*2700; % Aluminum: 2700kg/m3
                        m(4)=2*pi*r_u^2*L*800; % HTF: 800kg/m3
                        c(1)=750; % J/kg-K
                        c(3)=385; % J/kg-K
                        c(2)=obj.Reflector.cp(T(2));
                        c(4)=obj.HTF.cp(T(4));
                        q=obj.heatflux(T,j);
                        dT=q(:)./m(:)./c(:);
                        T=T(:)+dT*dt(i);
                        [T,~,exitflag]=fsolve(@(T) obj.heatflux(T,j),T,options);
                    end
                end
                if exitflag~=1
                    fprintf('Equations not solved. Update guess values.\n')
                end
                T_result{j}=T;
            end
            r=T_result;
            obj.temperature=r;
        end
        
        function obj=update(obj)
            [~,obj]=obj.go;
        end
        
        function [grid,obj]=meshing(obj,varargin)
            % varargin{1}: smallest element size
            if isempty(obj.optical.geometry.RadUtube)
                error('You forgot to give a u-tube radius, bro!')
            end
            [pt,t,tag,normvec,lface]=CPCmesh1d(obj.optical.geometry,varargin);
            description=cell(6,2);
            description(:,1)=num2cell([1;2;-3;-4;3;4]);
            description(:,2)={'Glass';...
                'Absorber (no shared part)';...
                'The face shared by absorber and u-tube inlet';...
                'The face shared by absorber and u-tube outlet';...
                'U-tube inlet (no shared part)';...
                'U-tube outlet (no shared part)'};
            tagsum=struct('description',1,...
                'tag',tag);
            tagsum.description=description;
            [yl1,yl2]=obj.optical.geometry.profile; yl=[yl1;yl2];
            yl=yl(~isnan(yl(:,1)),:); [~,I]=sort(yl(:,1)); yl=yl(I,:);
            dyl=diff(yl); lcpc=sum(sqrt(sum(dyl.^2,2)))*2;
            [row,col]=mkt2t_DX(t);
            grid=struct('points',pt,...
                'connectivity',t,...
                'tag',tagsum,...
                'normvec',normvec,...
                'lface',lface,...
                'pmid',(pt(t(:,1),:)+pt(t(:,2),:))/2,...
                'lcpc',lcpc,...
                'La',sum(lface(tag==2 | tag==-3 | tag==-4)),...
                'Lg',sum(lface(tag==1)),...
                'row',row,...
                'col',col);
        end
        
        function [FA,obj]=vf(obj)
            meshsimple=obj.meshing;
            tag=meshsimple.tag.tag;
            lface=meshsimple.lface;
            obj.Area.A_glass=meshsimple.Lg;
            obj.Area.A_absorber=sum(lface(tag==2 | tag==-3 | tag==-4));
            obj.Area.A_reflector=meshsimple.lcpc;
            vftemp=viewfactor(meshsimple,obj.optical.geometry);
            obj.VFA.FA_glass_sky=sum(vftemp(tag==1,end).*lface(tag==1));
            obj.VFA.FA_glass_reflector=obj.Area.A_glass-obj.VFA.FA_glass_sky;
            obj.VFA.FA_absorber_glass=obj.Area.A_absorber;
            obj.VFA.FA_reflector_sky=obj.Area.A_reflector-obj.VFA.FA_glass_reflector;
            FA=obj.VFA;
        end
        
        function [FA,obj]=vfa(obj)
            [FA,obj]=obj.vf;
        end
        
        function [q_sf,Beam,Diffuse,obj]=getsolarflux(obj,varargin)
            % getsolarflux(obj): calculate the solar flux for all outdoor
            % conditions in obj.outdoor
            % getsolarflux(obj,surp): Uses surfaceproperty from surp. If
            % surp is not specified, use the default
            % getsolarflux(obj,surp,angle,beam,diffuse): calculate the solar
            % flux for incident angle, b/d specified by angle, beam, and
            % diffuse
            in=[varargin(:);cell(4,1)]; in=in(1:4);
            obj.check;
            if numel(varargin)<2
                [~,ang]=obj.outdoor.angles(obj.optical.geometry.beta);
                [Beam,dsky,dgrd]=obj.outdoor.tilt(obj.rho_grd,obj.optical.geometry.beta);
                Diffuse=dsky+dgrd;
            elseif numel(in{2})~=numel(in{3}) || numel(in{2})~=numel(in{4})
                error('numel(Diffuse)~=numel(Beam) || numel(Diffuse)~=numel(angle). I''m tired of talking.')
            else
                ang=in{2};
                Beam=in{3};
                Diffuse=in{4};
            end
            q=obj.optical.getefficiency('total',in{1},[1 1 1],[]); %W/m2
            lns=obj.optical.geometry.L;
            ia=obj.optical.IncidentAngles;
            ii=interp1([-ia(end:-1:2);ia],1:2*numel(ia)-1,ang);
            ii=round(ii);
            q_sf=[q(:,ii)*lns*diag(Beam)+mean(q,2)*lns*Diffuse(:)';zeros(1,numel(ang))]; %W/m
            w_aper=obj.optical.aperture(ang);
            q_sf=q_sf*diag(w_aper); %W
            obj.solarflux=q_sf;
            obj.efficiency=q;
        end
        
        function T0=initemp(obj,varargin)
            % initemp(obj,GHI,T_amb)
            in=[varargin(:);cell(2,1)]; in=in(1:2);
            if isempty(in{1})
                GHI=obj.outdoor.GHI;
            else
                GHI=in{1};
            end
            if isempty(in{2})
                in{2}=obj.outdoor.DryBulb;
            end
            T_amb=in{2}; T_amb=T_amb(:)';
            T_in=obj.operation.T_in; T_in=repmat(T_in,1,numel(T_amb));
            A=obj.optical.geometry.Width*obj.optical.geometry.L;
            T_out=T_in+GHI(:)'*A*0.4/obj.HTF.cp(T_in)/obj.operation.MassFR;
            T0=[T_amb+5;T_out+10;T_amb+2;T_out];
        end
        
        function [q,Jac,varargout]=heatflux(obj,varargin)
            % heatflux(obj,T,q_sf,outdoor)
            % heatflux(obj,T,j)
            % Old form:
            % [q,Jac,qd]=heatflux(Tc,r_u,q_sf,MassFR,T_in,T_amb,T_sky,Vel)
            % global tag lface lns thick lcpc VFA epsilon sgm row col r_a r_g La Lg h w
            
            % Initial check
            in=[varargin(:);cell(4,1)]; in=in(1:3);
            if isempty(in{1})
                if isempty(obj.temperature)
                    error('No temperature, no heat flux.')
                else
                    Tc=obj.temperature{1};
                end
            else
                Tc=in{1};
            end
            if isempty(in{2})
                in{2}=1;
            end
            switch numel(in{2})
                case 1
                    i=in{2};
                    q_sf=obj.solarflux(:,i);
                    in{3}=obj.outdoor.getweather(i);
                    ambient=obj.Ambient.getdefault(i);
                otherwise
                    q_sf=in{2}(:,1);
                    if isempty(in{3})
                        error('You already gave a solar flux. Just keep giving a weather, OK?')
                    else
                        ambient=air();
                        ambient.setdefault(in{3}.DryBulb(1));
                    end
            end
            
            
            % Substitute the old inputs
            T_amb=in{3}.DryBulb(1);
            T_sky=T_amb-6;
            Vel=in{3}.Wspd(1);
            T_in=obj.operation.T_in;
            MassFR=obj.operation.MassFR;
            r_a=obj.optical.geometry.RadReceiver;
            r_g=obj.optical.geometry.RadGlass;
            t_g=obj.optical.geometry.ThickGlass;
            r_u=obj.optical.geometry.RadUtube;
            w=obj.optical.geometry.Width;
            L=obj.optical.geometry.L;
            h=obj.optical.geometry.Height;
            T_g=Tc(1); T_refl=Tc(2); T_a=Tc(3); T_htf=Tc(4);
            A_glass=obj.Area.A_glass;
            A_reflector=obj.Area.A_reflector;
            A_absorber=obj.Area.A_absorber;
            A=[A_glass;A_reflector;A_absorber];
            A_ut=2*pi*r_u*L;
            FA_glass_sky=obj.VFA.FA_glass_sky;
            FA_glass_reflector=obj.VFA.FA_glass_reflector;
            FA_absorber_glass=obj.VFA.FA_absorber_glass;
            FA_reflector_sky=obj.VFA.FA_reflector_sky;
            lcpc=obj.Area.A_reflector;
            p_vac=obj.optical.geometry.VacPres; %bar
            sgm=5.67e-8; %W/m2-K4
            
            % Convection
            % glass and ambient
            [~,h_glass_air]=convection.Qglass2air(Vel,T_amb,T_g,r_g,h,w,ambient);
            % reflector and ambient
            [~,h_refl_air_cavity,~,h_refl_air_back]=convection.Qrefl2air(Vel,T_refl,T_amb,h,w,lcpc,ambient);
            h_refl_air=h_refl_air_cavity+h_refl_air_back;
            % absorber and HTF
            h_absorber_htf=convection.htube2htf_simple(obj.HTF,(T_htf+T_in)/2,MassFR,r_u,2*L);
            % glass and absorber
            [~,h_absorber_glass]=convection.Annulus(T_a,T_g,r_a,r_g-t_g,p_vac); % Q in W/m h in W/m-K
            A_conv=[-h_glass_air*A_glass-h_absorber_glass*L,0,h_absorber_glass*L,0;...
                0,-h_refl_air*A_reflector,0,0;...
                h_absorber_glass*L,0,-h_absorber_glass*L-h_absorber_htf*A_ut,h_absorber_htf*A_ut/2;...
                0,0,h_absorber_htf*A_ut,-h_absorber_htf*A_ut/2-MassFR*cp_HTF(T_htf)];
            S_conv=[h_glass_air*A_glass*T_amb;...
                h_refl_air*A_reflector*T_amb;...
                h_absorber_htf*A_ut/2*T_in;...
                -h_absorber_htf*A_ut/2*T_in+MassFR*cp_HTF(T_in)*T_in];
            q_conv=A_conv*Tc+S_conv;
            Jac_conv=A_conv; Jac_conv(4,4)=Jac_conv(4,4)-MassFR*3.429*T_htf;
            
            % Radiation
            epsilon=[0.91;0.1;0.04];
            EA=diag(epsilon.*A./(1-epsilon));
            A_rad1=[0,FA_glass_reflector,0;...
                FA_glass_reflector,0,0;...
                0,0,0];
            A_rad1=EA+diag([A_glass;A_reflector;0])-A_rad1;
            S_rad1=EA*sgm*Tc(1:3,1).^4+...
                [FA_glass_sky;FA_reflector_sky;0]*sgm*T_sky^4;
            J1=A_rad1\S_rad1;
            q_rad1=EA*(J1-sgm.*Tc(1:3,1).^4);
            Jac_rad1=EA*(A_rad1\EA-eye(3))*4*sgm*diag(Tc(1:3,1).^3);
            
            A_rad2=[0,0,FA_absorber_glass;...
                0,0,0;...
                FA_absorber_glass,0,0];
            A_rad2=EA+diag(sum(A_rad2,2))-A_rad2;
            S_rad2=EA*sgm*Tc(1:3,1).^4;
            J2=A_rad2\S_rad2;
            q_rad2=EA*(J2-sgm.*Tc(1:3,1).^4);
            Jac_rad2=EA*(A_rad2\EA-eye(3))*4*sgm*diag(Tc(1:3,1).^3);
            
            q_rad=q_rad1+q_rad2;
            q_rad=[q_rad;0];
            Jac_rad=Jac_rad1+Jac_rad2;
            % Jac_rad=blkdiag(Jac_rad,0);
            Jac_rad=[Jac_rad,zeros(3,1);zeros(1,4)];
            
            %heat balance: f=[glass;reflector;absorber;fluid]
            %sf [absorber;reflector;glass;0] --> [glass;reflector;absorber;0]
            q_sf=q_sf([3,2,1,4]);
            q=q_conv+q_rad+q_sf;
            Jac=Jac_conv+Jac_rad;
            if nargout>2
                qd.q_UF=h_absorber_htf*A_ut*((T_in+Tc(4))/2-Tc(3));
                qd.q_FF=MassFR*(cp_HTF(T_in)*T_in-cp_HTF(Tc(4))*Tc(4));
                qd.q_AG=h_absorber_glass*L*(Tc(1)-Tc(3));
                qd.q_GA=h_glass_air*A_glass*(T_amb-Tc(1));
                qd.q_RA=h_refl_air*A_reflector*(T_amb-Tc(2));
                qd.q_RGS_Grad=q_rad1(1);
                qd.q_RGS_Rrad=q_rad1(2);
                qd.q_GGA_Grad=q_rad2(1);
                qd.q_GGA_Arad=q_rad2(3);
                varargout{1}=qd;
            end
        end
    end
    
    methods % Visualization
        function r=dispsolarflux(obj)
            if isempty(obj.solarflux)
                error('You can''t ask for what I don''t have.')
            else
                r.receiver=obj.solarflux(1,:);
                r.reflector=obj.solarflux(2,:);
                r.glass=obj.solarflux(3,:);
                r.units='W';
            end
        end
        
        function r=disptemperature(obj)
            if isempty(obj.temperature)
                error('You can''t ask for what I don''t have.')
            else
                t=cell2mat(obj.temperature);
                r.receiver=t(1,:);
                r.reflector=t(2,:);
                r.glass=t(3,:);
                r.htf=t(4,:);
                r.units='K';
            end
        end
    end
end