classdef NumericThermal
    properties
        optical;outdoor;operation;
        Nl=20;
        rho_grd=0.2;
        Reflector=aluminum();
        Receiver=copper();
        Glass=glass();
        HTF=Duratherm600();
        mesh;
    end
    
    properties (SetAccess=private)
        VF;VFA;solarflux;other;
        temperature;Ambient=air();
    end
    
    properties (Dependent)
        
    end
    
    methods % Constructor
        function r=NumericThermal(varargin)
            if nargin>0
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
                    r.optical.geometry.update;
                    return
                end
                if isa(varargin{1},'struct')
                    r.optical=RayTraceData(varargin{1});
                    r.optical.geometry=r.optical.geometry.update;
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
        
        function obj=set.Nl(obj,N)
            if isa(N,'numeric')
                obj.Nl=ceil(abs(N(1)));
            else
                error('I need a positive integer.')
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
        
        function obj=set.mesh(obj,m)
            if isa(m,'struct')
                m.points;m.connectivity;m.tag.tag;
                m.tag.description;m.normvec;m.lface;
                m.pmid;m.lcpc;m.La;m.Lg;m.row;m.col;
                obj.mesh=m;
            else
                error('Mesh is pretty complex. Why not use .meshing to generate one?')
            end
        end
        
        function obj=setsf(obj,q_sf)
            obj.solarflux=q_sf;
        end
        
        function obj=setambient(obj,amb)
            obj.Ambient=amb;
        end
        
        function r=check(obj)
            obj.optical.check;
            r=obj.optical.geometry.checknumeric;
            if isempty(obj.optical.geometry.Height) || isempty(obj.optical.geometry.Width)
                error('Wow, wow, wow. Hold on, pal. Height or Width are missing. Update it first.')
            end
            if isempty(obj.rho_grd)
                error('Wow, wow, wow. Hold on, pal. Set the ground reflectivity first: .rho_grd=')
            end
            if isempty(obj.Nl)
                error('Wow, wow, wow. Hold on, pal. Set the trough grid# first: .Nl=')
            end
            if isempty(obj.outdoor)
                fprintf('Just a suggestion. Set the outdoor condition: .outdoor=')
            end
        end
    end
    
    methods % Simulation
        function [r,obj]=go(obj,varargin)
            % go(obj): simulate the conditions of obj.solarflux
            % go(obj,angle,beam,diffuse,T_amb,Vel,op)
            in=[varargin(:);cell(6,1)]; in=in(1:6);
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
                obj.VFA=obj.vfa;
            end
            if isempty(op)
                op=repmat(obj.operation,numel(ang),1);
            end
            obj.Ambient=air();
            obj.Ambient=obj.Ambient.setdefault(obj.outdoor.DryBulb);
            G=Beam+Diffuse;
            T_result=cell(1,numel(ang));
            options = optimoptions('fsolve','Jacobian','on','Display','iter-detailed','MaxIter',20);
            dt=[5,50,500,5000];
            tag=obj.mesh.tag.tag;
            lface=obj.mesh.lface;
            lns=obj.optical.geometry.L/obj.Nl;
            t_g=obj.optical.geometry.ThickGlass;
            t_a=obj.optical.geometry.ThickReceiver;
            lcpc=obj.mesh.lcpc;
            t_cpc=obj.optical.geometry.ThickReflector;
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
                        m(tag==1,:)=repmat(lface(tag==1)*t_g*lns*2230,1,obj.Nl); % borosilicate glass: 2230kg/m3
                        m(tag==2 | abs(tag)==3 | abs(tag)==4,:)=repmat(lface(tag==2 | abs(tag)==3 | abs(tag)==4)*t_a*lns*8960,1,obj.Nl); % copper: 8960kg/m3
                        m(end-2,:)=lcpc*t_cpc*lns*2700; % Aluminum: 2700kg/m3
                        m(end-1:end,:)=pi*r_u^2*lns*800; % HTF: 800kg/m3
                        c(tag==1,:)=750; % J/kg-K
                        c(tag==2 | abs(tag)==3 | abs(tag)==4,:)=385; % J/kg-K
                        c(end-2,:)=obj.Reflector.cp(T(end-2,:));
                        c(end-1:end,:)=obj.HTF.cp(T(end-1:end,:));
                        q=obj.heatflux(T,j);
                        dT=q./m(:)./c(:);
                        T=reshape(T(:)+dT*dt(i),[],obj.Nl);
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
            [pt,t,tag,normvec,lface]=CPCmesh1d(obj.optical.geometry,varargin{:});
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
            obj.mesh=grid;
        end
        
        function [vfactor,obj]=vf(obj)
            if isempty(obj.mesh)
                [~,obj]=obj.meshing;
            end
            obj.VF=viewfactor(obj.mesh,obj.optical.geometry);
            vfactor=obj.VF;
        end
        
        function [FA,obj]=vfa(obj)
            if isempty(obj.VF)
                [~,obj]=obj.vf;
            end
            [r,c,vftemp]=find(obj.VF);
            lall=[obj.mesh.lface;obj.mesh.lcpc];
            vfatemp=vftemp.*lall(r);
            obj.VFA=sparse(r,c,vfatemp);
            FA=obj.VFA;
        end
        
        function [q_sf,Beam,Diffuse,obj]=getsolarflux(obj,varargin)
            % getsolarflux(obj): calculate the solar flux for all outdoor
            % conditions in obj.outdoor
            % getsolarflux(obj,surp): Uses surfaceproperty from surp. If
            % surp is not specified, use the default
            % getsolarflux(obj,surp,angle): calculate the solar flux for
            % incident angle specified by angle and beam/diffuse in the
            % first obj.outdoor
            % getsolarflux(obj,surp,angle,beam,diffuse): calculate the solar
            % flux for incident angle, b/d specified by angle, beam, and
            % diffuse
            % getsolarflux(obj,surp,[],beam,diffuse): calculate the solar flux
            % for the first incident angle in obj.outdoor and b/d specifed
            % by beam/diffuse
            if isempty(obj.mesh)
                obj.mesh=obj.meshing;
            end
            in=[varargin(:);cell(4,1)]; in=in(1:4);
            obj.check;
            if numel(varargin)<2
                [~,ang]=obj.outdoor.angles(obj.optical.geometry.beta);
                [Beam,dsky,dgrd]=obj.outdoor.tilt(obj.rho_grd,obj.optical.geometry.beta);
                Diffuse=dsky+dgrd;
            else
                if isempty(in{3}) || isempty(in{4})
                    [Beam,dsky,dgrd]=obj.outdoor.tilt(obj.rho_grd,obj.optical.geometry.beta,1);
                    Diffuse=dsky+dgrd;
                else
                    Beam=in{3};
                    Diffuse=in{4};
                end
                if isempty(in{2})
                    [~,ang]=obj.outdoor.angles(obj.optical.geometry.beta,1);
                    ang=repmat(ang,1,numel(Beam));
                else
                    ang=in{2};
                    if numel(Diffuse)~=numel(Beam)
                        error('numel(Diffuse)~=numel(Beam). I''m tired of talking.')
                    elseif numel(Diffuse)==1
                        Beam=repmat(Beam,1,numel(ang));
                        Diffuse=repmat(Diffuse,1,numel(ang));
                    end
                end
            end
            q=obj.optical.getefficiency('local',in{1},[1 1 1],[],obj.mesh); %W/m2
            lns=obj.optical.geometry.L/obj.Nl;
            ia=obj.optical.IncidentAngles;
            ii=interp1([-ia(end:-1:2);ia],1:2*numel(ia)-1,ang);
            ii=round(ii);
            q_sf=[q(:,ii)*lns*diag(Beam)+mean(q,2)*lns*Diffuse(:)';zeros(2,numel(ang))]; %W/m
            w_aper=obj.optical.aperture(ang);
            q_sf=q_sf*diag(w_aper); %W
            obj.solarflux=q_sf;
        end
        
        function T0=initemp(obj,varargin)
            % initemp(obj,GHI,outdoor)
            % 1 to size(t,1) is temperature for glass, absorber, u-tube
            % size(t,1)+1 is reflector temperature
            % size(t,1)+2 is inlet HTF
            % size(t,1)+3 is outlet HTF
            in=[varargin(:);cell(2,1)]; in=in(1:2);
            if isempty(in{1})
                GHI=obj.outdoor.GHI;
            else
                GHI=in{1};
            end
            if isempty(in{2})
                in{2}=obj.outdoor.DryBulb;
            end
            
            T_amb=in{2};
            tag=obj.mesh.tag.tag;
            N=obj.Nl;
            T_in=obj.operation.T_in;
            A=obj.optical.geometry.Width*obj.optical.geometry.L;
            T_out=T_in+GHI(1)*A*0.4/obj.HTF.cp(T_in)/obj.operation.MassFR;
            
            T=T_amb*ones(size(tag,1)+1,N);
            T(tag==1,:)=T_amb+5;
            T(tag==2 | abs(tag)==3 | abs(tag)==4,:)=T_out+10;
            T(end,:)=T_amb+2;
            lambda1=(0:N-1)/(2*N-1);
            lambda2=(2*N-1:-1:N)/(2*N-1);
            T_htf1=(1-lambda1)*T_in+lambda1*T_out;
            T_htf2=(1-lambda2)*T_in+lambda2*T_out;
            T0=[T;T_htf1;T_htf2];
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
            tag=obj.mesh.tag.tag;
            lface=obj.mesh.lface;
            lcpc=obj.mesh.lcpc;
            lns=obj.optical.geometry.L/obj.Nl;
            sgm=5.67e-8; %W/m2-K4
            row=obj.mesh.row;
            col=obj.mesh.col;
            r_a=obj.optical.geometry.RadReceiver;
            r_g=obj.optical.geometry.RadGlass;
            r_u=obj.optical.geometry.RadUtube;
            La=obj.mesh.La;
            Lg=obj.mesh.Lg;
            h=obj.optical.geometry.Height;
            w=obj.optical.geometry.Width;
            MassFR=obj.operation.MassFR;
            T_in=obj.operation.T_in;
            T_amb=in{3}.DryBulb(1);
            T_sky=T_amb-6;
            Vel=in{3}.Wspd(1);
            t_g=obj.optical.geometry.ThickGlass;
            t_a=obj.optical.geometry.ThickReceiver;
            t_cpc=obj.optical.geometry.ThickReflector;
            thick=[t_g;t_a;t_a;t_a;t_cpc];
            epsilon_g=0.91; epsilon_abs=0.04; epsilon_ref=0.1;
            epsilon=[epsilon_g;epsilon_abs;epsilon_abs;epsilon_abs;epsilon_ref;1];
            p_vac=obj.optical.geometry.VacPres;
            % Let's go
            [Ncell,N]=size(Tc);
            tagall=[tag;5;7];
            % 1 glass
            % 2 absorber
            % 3 u-tube inlet  -3 u-tube inlet shared with absorber
            % 4 u-tube outlet  -4 u-tube outlet shared with absorber
            % 5 reflector
            % 7 sky
            lall=[lface;lcpc];
            
            % Conduction E-W
            Row_EW_cond=repmat(row,1,N)+repmat((0:N-1)*Ncell,size(row,1),1); Row_EW_cond=Row_EW_cond(:);
            Col_EW_cond=repmat(col,1,N)+repmat((0:N-1)*Ncell,size(col,1),1); Col_EW_cond=Col_EW_cond(:);
            k=zeros(Ncell,N);
            k(Ncell-1:Ncell,:)=obj.HTF.k(Tc(Ncell-1:Ncell,:)); %htf
            k(Ncell-2,:)=obj.Reflector.k(Tc(Ncell-2,:)); %reflector
            k(tag==1,:)=obj.Glass.k(Tc(tag==1,:)); %W/m-K
            k(tag==2 | abs(tag)==3 | abs(tag)==4,:)=obj.Receiver.k(Tc(tag==2 | abs(tag)==3 | abs(tag)==4,:)); %W/m-K
            D_ew_per_area=2./(repmat(lface(row),N,1)./k(Row_EW_cond)+repmat(lface(col),N,1)./k(Col_EW_cond));   %W/m2-K
            D_ew=D_ew_per_area.*repmat(thick(abs(tag(row))),N,1)*lns;   %W/K
            
            % Conductivity N-S
            Row_NS_cond=repmat((1:Ncell-3)',1,N-1)+repmat((0:N-2)*Ncell,Ncell-3,1); Row_NS_cond=Row_NS_cond(:);
            Col_NS_cond=Row_NS_cond+Ncell;
            D_ns_per_area=2*k(Col_NS_cond)./(lns+lns*k(Col_NS_cond)./k(Row_NS_cond));
            D_ns=D_ns_per_area.*repmat(thick(abs(tag)).*lface,N-1,1);
            Row_refl_cond=(Ncell-2:Ncell:Ncell*(N-1)-2)';
            Col_refl_cond=Row_refl_cond+Ncell;
            D_refl=2*k(Row_refl_cond)./(lns+lns*k(Row_refl_cond)./k(Col_refl_cond))*thick(5)*lcpc;
            Row_htf=[(Ncell-1:Ncell:Ncell*(N-1)-1)';...
                Ncell*N-1;...
                (Ncell*N:-Ncell:2*Ncell)'];
            Col_htf=[(Ncell*2-1:Ncell:Ncell*N-1)';...
                Ncell*N;...
                (Ncell*(N-1):-Ncell:Ncell)'];
            D_htf=2./(lns./k(Row_htf)+lns./k(Col_htf))*pi*r_u^2; %W/K
            
            % Conduction all
            cond_P=sparse([Row_EW_cond;Col_EW_cond;Row_NS_cond;Col_NS_cond;Row_refl_cond;Col_refl_cond;Row_htf;Col_htf],...
                [Col_EW_cond;Row_EW_cond;Col_NS_cond;Row_NS_cond;Col_refl_cond;Row_refl_cond;Col_htf;Row_htf],...
                [D_ew;D_ew;D_ns;D_ns;D_refl;D_refl;D_htf;D_htf],Ncell*N,Ncell*N);
            cond_P=sum(cond_P,2); cond_P(Ncell-1)=cond_P(Ncell-1)+2/lns*k_HTF(T_in)*pi*r_u^2;
            A_cond=sparse([Row_EW_cond;Col_EW_cond;Row_NS_cond;Col_NS_cond;Row_refl_cond;Col_refl_cond;Row_htf;Col_htf;(1:Ncell*N)'],...
                [Col_EW_cond;Row_EW_cond;Col_NS_cond;Row_NS_cond;Col_refl_cond;Row_refl_cond;Col_htf;Row_htf;(1:Ncell*N)'],...
                [D_ew;D_ew;D_ns;D_ns;D_refl;D_refl;D_htf;D_htf;-cond_P]);
            S_cond=sparse(Ncell-1,1,2/lns*k_HTF(T_in)*pi*r_u^2*T_in,Ncell*N,1);
            q_cond=A_cond*Tc(:)+S_cond;
            Jac_cond=A_cond;
            
            % Convection
            S_conv=zeros(Ncell,N);
            % U-tube and HTF
            iu_i=abs(tag)==3; iu_o=abs(tag)==4;
            Row_UF_conv=repmat([find(iu_i);find(iu_o)],1,N)+repmat((0:N-1)*Ncell,sum(iu_i+iu_o),1);
            Row_UF_conv=Row_UF_conv(:);
            Col_UF_conv=repmat([(Ncell-1)*ones(sum(iu_i),1);Ncell*ones(sum(iu_o),1)],1,N)...
                +repmat((0:N-1)*Ncell,sum(iu_i+iu_o),1);
            Col_UF_conv=Col_UF_conv(:);
            HA_UF_conv=[lface(iu_i)*convection.htube2htf(obj.HTF,Tc(Ncell-1,:),MassFR,r_u,lns*(0.5:N-0.5))*lns;...
                lface(iu_o)*convection.htube2htf(obj.HTF,Tc(Ncell,:),MassFR,r_u,lns*(2*N-0.5:-1:N+0.5))*lns]; %W/K
            HA_UF_conv=HA_UF_conv(:);
            A_UF_conv=sparse([Row_UF_conv;Col_UF_conv],[Col_UF_conv;Row_UF_conv],...
                [HA_UF_conv;HA_UF_conv],Ncell*N,Ncell*N);
            P_UF_conv=sum(A_UF_conv,2);
            A_UF_conv=sparse([Row_UF_conv;Col_UF_conv;(1:Ncell*N)'],[Col_UF_conv;Row_UF_conv;(1:Ncell*N)'],...
                [HA_UF_conv;HA_UF_conv;-P_UF_conv],Ncell*N,Ncell*N);
            
            % HTF and HTF
            HA_FF_conv1=MassFR*cp_HTF(Tc(Row_htf)); % corresponding to [Col_htf,Row_htf]
            HA_FF_conv2=-MassFR*cp_HTF(Tc([Row_htf;Ncell])); % corresponding to [Row_htf;Ncell],[Row_htf;Ncell]
            S_conv(Ncell-1,1)=cp_HTF(T_in)*MassFR*T_in;
            A_FF_conv=sparse([Col_htf;Row_htf;Ncell],[Row_htf;Row_htf;Ncell],...
                [HA_FF_conv1;HA_FF_conv2],Ncell*N,Ncell*N);
            S_conv_FF=zeros(Ncell,N);
            S_conv_FF(Ncell-1,1)=S_conv(Ncell-1,1);
            
            % Between absorber and glass
            ig=tag==1; ia=tag==2 | tag==-3 | tag==-4;
            T_g_avg=mean(Tc(ig,:));
            T_a_avg=mean(Tc(ia,:));
            [~,h_abs_glass]=convection.Annulus(T_a_avg,T_g_avg,r_a,r_g-t_g,p_vac); % Q in W/m h in W/m-K
            Col_AG_conv=repmat([find(ig);find(ia)],1,sum(ig+ia)); Row_AG_conv=Col_AG_conv';
            HA_AG_conv=[-ones(sum(ig),1)/sum(ig);ones(sum(ia),1)/sum(ia)]*[lface(ig)/Lg;-lface(ia)/La]'; % (sum(T_g)/ng-sum(T_a)/na)*Lgi/Lg
            Col_AG_conv=repmat(Col_AG_conv(:),1,N)+repmat(0:Ncell:Ncell*(N-1),numel(Col_AG_conv),1);
            Row_AG_conv=repmat(Row_AG_conv(:),1,N)+repmat(0:Ncell:Ncell*(N-1),numel(Row_AG_conv),1);
            HA_AG_conv=HA_AG_conv(:)*h_abs_glass*lns;
            A_AG_conv=sparse(Row_AG_conv,Col_AG_conv,HA_AG_conv,Ncell*N,Ncell*N);
            
            % Between glass and ambient
            [~,h_glass_air]=convection.Qglass2air(Vel,T_amb,T_g_avg,r_g,h,w,ambient);
            Col_GA_conv=repmat(find(ig),1,sum(ig)); Row_GA_conv=Col_GA_conv';
            HA_GA_conv=-ones(sum(ig),1)/sum(ig)*lface(ig)'/Lg;
            Col_GA_conv=repmat(Col_GA_conv(:),1,N)+repmat(0:Ncell:Ncell*(N-1),numel(Col_GA_conv),1);
            Row_GA_conv=repmat(Row_GA_conv(:),1,N)+repmat(0:Ncell:Ncell*(N-1),numel(Row_GA_conv),1);
            HA_GA_conv=HA_GA_conv(:)*h_glass_air*lns*2*pi*r_g;
            A_GA_conv=sparse(Row_GA_conv,Col_GA_conv,HA_GA_conv,Ncell*N,Ncell*N);
            S_conv(ig,:)=lface(ig)/Lg*h_glass_air*lns*2*pi*r_g*T_amb;
            S_conv_GA=zeros(Ncell,N);
            S_conv_GA(ig,:)=S_conv(ig,:);
            
            % Between reflector and ambient
            [~,h_refl_air_cavity,~,h_refl_air_back]=convection.Qrefl2air(Vel,Tc(Ncell-2,:),T_amb,h,w,lcpc,ambient);
            h_refl_air=h_refl_air_cavity+h_refl_air_back;
            Row_RA_conv=(Ncell-2:Ncell:Ncell*N-2)';
            HA_RA_conv=h_refl_air'*lcpc*lns;
            A_RA_conv=sparse(Row_RA_conv,Row_RA_conv,-HA_RA_conv,Ncell*N,Ncell*N);
            S_conv(Ncell-2,:)=HA_RA_conv*T_amb;
            S_conv_RA=zeros(Ncell,N);
            S_conv_RA(Ncell-2,:)=S_conv(Ncell-2,:);
            
            % All convection
            A_conv=A_UF_conv+A_FF_conv+A_AG_conv+A_GA_conv+A_RA_conv;
            q_conv=A_conv*Tc(:)+S_conv(:);
            Jac_conv=A_conv+sparse([Col_htf;Row_htf;Ncell],[Row_htf;Row_htf;Ncell],...
                [MassFR*3.429*Tc(Row_htf);-MassFR*3.429*Tc([Row_htf;Ncell])],Ncell*N,Ncell*N);
            
            % Solve radiative flux
            [ra,ca,vfaa]=find(obj.VFA);
            Aepsilon=epsilon(abs(tagall(1:end-1)))./(1-epsilon(abs(tagall(1:end-1)))).*lall*lns;
            K=sparse(1:Ncell-2,1:Ncell-2,Aepsilon);
            B=4*sgm*sparse(1:Ncell*N,1:Ncell*N,Tc(:).^3);
            S_rad=repmat(Aepsilon,1,N)*sgm.*Tc(1:Ncell-2,:).^4;
            
            % Among reflector, glass, and sky
            ir1=tagall(ra)==1;
            ic1=tagall(ca)==5 | tagall(ca)==7;
            i1=ir1 & ic1; % all entries from glass to sky, reflector
            ir2=tagall(ra)==5;
            ic2=tagall(ca)==1 | tagall(ca)==7;
            i2=ir2 & ic2; % all entries from reflector to sky, glass
            i=i1 | i2;
            VFA1=sparse(ra(i),ca(i),vfaa(i),Ncell-2,Ncell-1);
            rad_P=sum(VFA1,2)+Aepsilon;
            [r,c,vfa]=find(VFA1(:,1:end-1));
            A_rad=sparse([r;(1:Ncell-2)'],[c;(1:Ncell-2)'],[-vfa;rad_P]);
            S_r=repmat(VFA1(:,end),1,N)*sgm*T_sky^4+S_rad;
            J_rad=A_rad\S_r;
            q_rad1=(J_rad-sgm*Tc(1:Ncell-2,:).^4).*repmat(Aepsilon,1,N);
            BLK=blkdiag(K*(A_rad\K-speye(Ncell-2)),sparse(2,2));
            BLK=repmat({BLK},1,N);
            BLK=blkdiag(BLK{:});
            Jac_rad_1=BLK*B;
            
            % Among glass, glass and absorber
            ir=tagall(ra)==1 | tagall(ra)==2 | tagall(ra)==-3 | tagall(ra)==-4;
            ic=tagall(ca)==1 | tagall(ca)==2 | tagall(ca)==-3 | tagall(ca)==-4;
            i=ir & ic;
            VFA2=sparse(ra(i),ca(i),vfaa(i),Ncell-2,Ncell-2);
            rad_P=sum(VFA2,2)+Aepsilon;
            [r,c,vfa]=find(VFA2);
            A_rad=sparse([r;(1:Ncell-2)'],[c;(1:Ncell-2)'],[-vfa;rad_P]);
            J_rad=A_rad\S_rad;
            q_rad2=(J_rad-sgm*Tc(1:Ncell-2,:).^4).*repmat(Aepsilon,1,N);
            BLK=blkdiag(K*(A_rad\K-speye(Ncell-2)),sparse(2,2));
            BLK=repmat({BLK},1,N);
            BLK=blkdiag(BLK{:});
            Jac_rad_2=BLK*B;
            
            % Among absorber and absorber
            ir=tagall(ra)==2 | tagall(ra)==3 | tagall(ra)==4;
            ic=tagall(ca)==2 | tagall(ca)==3 | tagall(ca)==4;
            i=ir & ic;
            VFA3=sparse(ra(i),ca(i),vfaa(i),Ncell-2,Ncell-2);
            rad_P=sum(VFA3,2)+Aepsilon;
            [r,c,vfa]=find(VFA3);
            A_rad=sparse([r;(1:Ncell-2)'],[c;(1:Ncell-2)'],[-vfa;rad_P]);
            J_rad=A_rad\S_rad;
            q_rad3=(J_rad-sgm*Tc(1:Ncell-2,:).^4).*repmat(Aepsilon,1,N);
            BLK=blkdiag(K*(A_rad\K-speye(Ncell-2)),sparse(2,2));
            BLK=repmat({BLK},1,N);
            BLK=blkdiag(BLK{:});
            Jac_rad_3=BLK*B;
            
            % all radiation
            q_rad=[q_rad1+q_rad2+q_rad3;sparse(2,N)];
            Jac_rad=Jac_rad_1+Jac_rad_2+Jac_rad_3;
            
            q=q_cond+q_conv+q_rad(:)+repmat(q_sf,N,1);
            if nargout>2
                qd.q_cond=reshape(q_cond,[],N);
                qd.q_UF=reshape(A_UF_conv*Tc(:),[],N);
                qd.q_FF=reshape(A_FF_conv*Tc(:)+S_conv_FF(:),[],N);
                qd.q_AG=reshape(A_AG_conv*Tc(:),[],N);
                qd.q_GA=reshape(A_GA_conv*Tc(:)+S_conv_GA(:),[],N);
                qd.q_RA=reshape(A_RA_conv*Tc(:)+S_conv_RA(:),[],N);
                qd.q_RGS=[q_rad1;zeros(2,N)];
                qd.q_GGA=[q_rad2;zeros(2,N)];
                qd.q_AA=[q_rad3;zeros(2,N)];
                varargout{1}=qd;
            end
            Jac=Jac_cond+Jac_conv+Jac_rad;
        end
        
        function [T,q]=avgtq(obj,T0,qd)
            tag=obj.mesh.tag.tag;
            T.abs=mean(mean(T0(tag==2 | abs(tag)==3 | abs(tag)==4,:),2));
            T.glass=mean(mean(T0(tag==1,:),2));
            T.refl=mean(T0(end-2,:));
            T.htf=T0(end,1);
            q.q_UF=sum(sum(qd.q_UF(tag==2 | abs(tag)==3 | abs(tag)==4,:)));
            q.q_FF=sum(sum(qd.q_FF(end-1:end,:)));
            q.q_AG=sum(sum(qd.q_AG(tag==2 | abs(tag)==3 | abs(tag)==4,:)));
            q.q_GA=sum(sum(qd.q_GA(tag==1,:)));
            q.q_RA=sum(sum(qd.q_RA(end-2,:)));
            q.q_RGS_Grad=sum(sum(qd.q_RGS(tag==1,:)));
            q.q_RGS_Rrad=sum(sum(qd.q_RGS(end-2,:)));
            q.q_GGA_Grad=sum(sum(qd.q_GGA(tag==1,:)));
            q.q_GGA_Arad=sum(sum(qd.q_GGA(tag==2 | tag==-3 | tag==-4,:)));
        end
        
    end
    
    methods % Visualization
        function r=showmesh2d(obj)
            r=checktypequiver(obj.mesh);
        end
        
        function r=showmesh3d(obj)
            tag=obj.mesh.tag.tag;
            t=obj.mesh.connectivity;
            pt=obj.mesh.points;
            N=obj.Nl; Np=size(pt,1); 
            L=obj.optical.geometry.L;
            lns=L/N;
            ptx=repmat(pt(:,1),N+1,1);
            pty=ones(Np,1)*lns*(0:N); pty=pty(:);
            ptz=repmat(pt(:,2),N+1,1);
            pn=[ptx,pty,ptz];
            t=repmat(t,N,1)+reshape(repmat([(0:N-1),(0:N-1)]*Np,size(t,1),1),[],2);
            tn=[t,t(:,2)+Np,t(:,1)+Np];
            ax=cell(1,2);
            ax{1}=subplot(1,2,1);
            ax{2}=subplot(1,2,2);
            pcol=[1 1 1];
            ii=repmat(tag==1,N,1);
            patch('vertices',pn,'faces',tn(ii,:),...
                'edgecolor','k','facecolor','flat',...
                'facevertexcdata',pcol,'parent',ax{1});
            ii=repmat(tag==2 | abs(tag)==3 | abs(tag)==4,N,1);
            patch('vertices',pn,'faces',tn(ii,:),...
                'edgecolor','k','facecolor','flat',...
                'facevertexcdata',pcol,'parent',ax{2});
            axis(ax{1},'off')
            axis(ax{2},'off')
            view(ax{1},30,30)
            view(ax{2},30,30)
            set(gcf,'color',[1 1 1])
            r=gcf;
            r.Position=[240,200,800,400];
        end
        
        function r=showtemperature(obj,varargin)
            % y=VisualCPC(T,pt,t,L)
            ind=ismember({'receiver','glass','htf'},lower(varargin));
            m=sum(ind);
            if m==0
                m=3;
                ind=[true true true];
            end
            id=find(ind);
            tag=obj.mesh.tag.tag;
            T=obj.temperature;
            t=obj.mesh.connectivity;
            pt=obj.mesh.points;
            L=obj.optical.geometry.L;
            N=size(T{1},2);
            part3d=[tag==2 | abs(tag)==3 | abs(tag)==4,...
                tag==1,...
                abs(tag)==3 | abs(tag)==4];
%             part3d=repmat(part3d,N,1);
            partmesh=[tag==2 | tag==-3 | tag==-4,...
                tag==1];
            ttl={'Receiver','Glass','HTF'};
            ax=cell(2,m);
            figure
            y=gcf;
            y.Position=[40,100,1200,600];
            for j=1:m
                ax{1,j}=subplot(2,m,j);
%                 ax{1,j}.Position=[(j-1)/m+1/m*0.2,0.55,1/m*0.6,0.4];
                ax{2,j}=subplot(2,m,j+m);
%                 ax{2,j}.Position=[(j-1)/m+1/m*0.2,0.05,1/m*0.6,0.4];
            end
            h=cell(2,m);
            x=L/N/2:L/N:L;
            for i=1:numel(T)
                for j=1:m
                    i3d=part3d(:,id(j));
                    Ttemp=T{i};
                    if id(j)==3
                        Ttemp(abs(tag)==3,:)=repmat(Ttemp(end-1,:),sum(abs(tag)==3),1);
                        Ttemp(abs(tag)==4,:)=repmat(Ttemp(end,:),sum(abs(tag)==4),1);
                    end
                    h{1,j}=VisualCPC(ax{1,j},Ttemp(i3d,:),pt,t(i3d,:),L);
                    title(ax{1,j},ttl{id(j)})
                    
                    switch id(j)
                        case 3
                            Ttemp=[T{i}(end-1,:),T{i}(end,end:-1:1)];
                            xm=[x,x(end:-1:1)];
                            h{2,j}=plot(ax{2,j},xm,Ttemp);
                            xlabel(ax{2,j},'Local distance (m)')
                            ylabel(ax{2,j},'Fluid Temperature ({\circ}C)')
                        otherwise
                            Ttemp=T{i};
                            imesh=partmesh(:,id(j));
                            Ttemp=Ttemp(imesh,:);
                            pmid=(pt(t(imesh,1),:)+pt(t(imesh,2),:))/2;
                            a=-atan2(-pmid(:,1),pmid(:,2))*180/pi;
                            [a,isort]=sort(a);
                            [xm,am]=meshgrid(x,a);
                            h{2,j}=mesh(ax{2,j},am,xm,Ttemp(isort,:)); %#ok<CPROP>
                            xlabel(ax{2,j},'Local angle ({\circ})')
                            ylabel(ax{2,j},'Local distance (m)')
                            zlabel(ax{2,j},'Temperature ({\circ}C)')
                            view(ax{2,j},40,50)
                            colormap default
                    end
                end
                pause(0.01)
                if i==numel(T)
                    break
                end
                cellfun(@delete,h)
            end
            r=gcf;
        end
    end
    
    methods % Utility
        function r=description(obj)
            r=obj;
        end
    end
end