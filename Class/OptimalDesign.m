classdef OptimalDesign < SimpleThermal
    properties (SetAccess=private)
        outdoor_bin; energy; hours;
        intensity;Coef;
%         efficiency; 
    end
    properties (Access=private)
        other;
    end
    properties (Dependent)
    end
    methods % Constructor
        function r=OptimalDesign(val)
            if nargin>0
                if isa(val,'SimpleThermal')
                    r.optical=val.optical;
                    r.outdoor=val.outdoor;
                    r.operation=val.operation;
                    r.HTF=val.HTF;
                    r.Reflector=val.Reflector;
                    r.Receiver=val.Receiver;
                    r.Glass=val.Glass;
                    r.rho_grd=val.rho_grd;
                    return
                end
                if isa(val,'RayTraceData')
                    r.optical=val;
                    return
                end
                if isa(val,'CPC')
                    r.optical=RayTraceData(val);
                    return
                end
                if isa(val,'struct')
                    r.optical=RayTraceData(val);
                    return
                end
            end
        end
        
        function [c,obj]=getCoef(obj)
            bl=cellfun(@isempty,{obj.optical.x,obj.optical.y,obj.optical.index,obj.optical.ia});
            if any(bl) % do interpolation
                obj.optical.geometry=obj.optical.geometry.update;
                q=obj.interpRT;
            else
                q=obj.optical.getefficiency('total',[],[1 1 1],obj.optical.IncidentAngles); %W/m2
            end
            obj.efficiency=q;
            temp=obj; 
            Vel=mean(temp.outdoor.Wspd)*ones(4,1);
            T_amb=mean(temp.outdoor.DryBulb)*ones(4,1);
            ang=[0;0;round(temp.optical.geometry.AccHalfAng/2);round((90+temp.optical.geometry.AccHalfAng)/2)];
            B=[1000;0;800;200]; d=[0;1000;200;300];
            [T,temp]=temp.go(ang,B,d,T_amb,Vel);
            T=cell2mat(T); T_o=T(4,:);
            T_in=temp.operation.T_in;
            q=temp.operation.MassFR*(temp.HTF.cp(T_o).*T_o...
                -temp.HTF.cp(T_in)*T_in);
            A=[temp.solarflux(1,:)',T_amb-T_in];
            c=(A'*A)\A'*q(:);
            obj.Coef=c;
        end
        
        function [q,obj]=getefficiency(obj)
            bl=cellfun(@isempty,{obj.optical.x,obj.optical.y,obj.optical.index,obj.optical.ia});
            obj.optical.geometry=obj.optical.geometry.update;
            if any(bl) % do interpolation
                q=obj.interpRT;
            else
                q=obj.optical.getefficiency('total',[],[1 1 1],obj.optical.IncidentAngles); %W/m2
            end
            obj.efficiency=q;
        end
    end
    methods % Simulation
        function [DV,obj,exitflag]=maxannual(obj,varargin)
            % maxannual(obj,DV0,lb,ub)
            % DesignVariable
            % 1. theta_deg:    acceptance half angle
            % 2. rcr:          relative concentration ratio
            % 3. rr_a:         relative radius of absorber
            % 4. rr_g:         relative radius of glass envelope
            % 5. rt_g:         relative thickness of glass envelope
            % 6. rgap:         relative gap distance between glass and reflector
            % 7. rr_u:         relative radius of u-tube
            % 8. L:            relative length of the trough
            % 9. beta:         tilted angle
            in=[varargin(:);cell(3,1)]; in=in(1:3);
            % Default constants: r_a=0.05, rtg=0.2, rgap=0, L=L
            ub=[80,1,0.05,2,0.2,0,0.4,obj.optical.geometry.L,obj.outdoor.latitude+25];
            lb=[10,0,0.05,1.1,0.2,0,0.1,obj.optical.geometry.L,0];
            id=[1,2,4,5,6];
            if ~isempty(in{2})
                in{2}=in{2}(:)';
                if any(lb(id)-in{2}(id)>0)
                    fprintf('Default\tSupplied\n')
                    fprintf('%.2f\t%.2f\n',[lb(id)',in{2}(id)'])
                    error('ID: 1,2,4,5,6 of lower bound is even lower than the limit')
                else
                    lb=in{2};
                end
            end
            if ~isempty(in{3})
                in{3}=in{3}(:)';
                if any(in{3}(id)-ub(id)>0)
                    fprintf('Default\tSupplied\n')
                    fprintf('%.2f\t%.2f\n',[ub(id)',in{3}(id)'])
                    error('ID: 1,2,4,5,6 of upper bound is even greater than the limit')
                else
                    ub=in{3};
                end
            end
            if isempty(in{1})
                DV0=[60,0.5,0.04,1.1,... %1~4
                    0.2,0.1,0.2,... %5~7
                    obj.optical.geometry.L,25];
%                 DV0=obj.initdv(lb,ub);
            else
                DV0=varargin{1};
            end
            options = optimoptions('fmincon','Display','iter-detailed','maxiter',20);
            [DV,~,exitflag]=fmincon(@(DV) -obj.annualobfun(DV),DV0,[],[],[],[],lb,ub,[],options);
            dv=DV;
            cpcgeo=CPC(dv(1),... % acceptance half angle
                [],... % CR
                dv(3),... % receiver radius
                dv(3)*dv(4),... % glass radius
                dv(3)*(dv(4)-1)*dv(5),... % glass thickness
                dv(3)*dv(6),... % gap
                dv(3)*dv(7));  % u-tube radius
            cpcgeo.L=dv(8);
            cpcgeo.beta=dv(9);
            [w1,w2]=cpcgeo.minmaxwidth;
            cpcgeo.CR=dv(2)*(w2/w1-1)+1;
            cpcgeo.VacPres=obj.optical.geometry.VacPres;
            [~,obj]=obj.annual(cpcgeo);
        end
        
        function [qdot,obj]=annualobfun(obj,DesignVariable)
            % DesignVariable
            % 1. theta_deg:    acceptance half angle
            % 2. rcr:          relative concentration ratio
            % 3. rr_a:         relative radius of absorber
            % 4. rr_g:         relative radius of glass envelope
            % 5. rt_g:         relative thickness of glass envelope
            % 6. rgap:         relative gap distance between glass and reflector
            % 7. rr_u:         relative radius of u-tube
            % 8. L:            relative length of the trough
            % 9. beta:         tilted angle
            
            % Construct an OptimalDesign object
            dv=DesignVariable;
%             cpcgeo=CPC(dv(1),... % acceptance half angle
%                 [],... % CR
%                 dv(3),... % receiver radius
%                 dv(3)*dv(4),... % glass radius
%                 dv(3)*(dv(4)-1)*dv(5),... % glass thickness
%                 dv(3)*dv(6),... % gap
%                 dv(3)*dv(7));  % u-tube radius
%             cpcgeo.L=dv(8);
%             cpcgeo.beta=dv(9);
%             [w1,w2]=cpcgeo.minmaxwidth;
%             cpcgeo.CR=dv(2)*(w2/w1-1)+1;
%             cpcgeo.VacPres=1e-5;
            cpcgeo=OptimalDesign.dv2cpc(dv);
%             cpcgeo=cpcgeo.update;
            [qdot,obj]=obj.annual(cpcgeo,1);
        end
        
        function [qdot,obj]=annual(obj,varargin)
            % annual(obj,a class to specify a new geometry)
            in=[varargin(:);cell(2,1)]; in=in(1:2);
            if isempty(in{1})
                rt=obj.optical;
            elseif isa(in{1},'CPC')
                rt=RayTraceData(in{1});
            elseif isa(in{1},'RayTraceData')
                rt=in{1};
            else
                error('You can take a rest sitting there giving nothing. Otherwise, you just need to give me a CPC class, better with ray trace data (RayTraceData class).')
            end
            rt.geometry.checksimple;
            obj.optical=rt;
%             bl=cellfun(@isempty,{obj.optical.x,obj.optical.y,obj.optical.index,obj.optical.ia});
%             obj.optical.geometry=obj.optical.geometry.update;
%             if any(bl) % do interpolation
%                 q=obj.interpRT;
%             else
%                 q=obj.optical.getefficiency('total',[],[1 1 1],obj.optical.IncidentAngles); %W/m2
%             end
            [~,obj]=obj.getefficiency;
            obj.other.weather=obj.outdoor; % store the original outdoor condition
            % Remember to update geometry before this line
            Aper=obj.optical.geometry.Width*obj.optical.geometry.L;
            obj.operation.MassFR=Aper/20;
            % Bin method
            iabin=0:89;
            [wth,beam,d1,d2,obj.hours]=obj.binangle(iabin);
            diffuse=d1+d2;
            obj.other.beam=beam;
            obj.other.diffuse=diffuse;
            % Assign a modified weather data to outdoor
            obj.outdoor=wth;
            % gogogo
            id=obj.hours~=0;
            Tr=obj.go(iabin(id),beam(id),diffuse(id),wth.DryBulb(id),wth.Wspd(id));
            T=cell(1,90);
            T(id)=Tr;
            T_in=obj.operation.T_in;
            T(~id)=mat2cell(repmat([300;300;T_in;T_in],1,sum(~id)),4,ones(1,sum(~id)));
            % Restore the outdoor
            obj.outdoor_bin=wth;
            obj.outdoor=obj.other.weather;
            obj.temperature=T;
            % Energy
            T_o=cell2mat(T); T_o=T_o(4,:);
            E=obj.operation.MassFR*(obj.HTF.cp(T_o).*T_o...
                -obj.HTF.cp(T_in)*T_in).*obj.hours';
            obj.energy=E;
            qdot=sum(E)/Aper;
            obj.intensity=qdot;
            if isempty(in{2})
                if ispc
                    sep=' ';
                else
                    sep='';
                end
                p=['/Users/Donghao',sep,'Xu/Dropbox/DonghaoCode/Data/Result'];
                try
                    obj.secretdb(fullfile(p,[date,'.txt']));
                catch
                    p='/Users/Donghao/Dropbox/DonghaoCode/Data/Result';
                    obj.secretdb(fullfile(p,[date,'.txt']));
                end
            end
        end
        
        function [qdot,obj]=annualutd(obj,varargin)
            in=[varargin(:);cell(2,1)]; in=in(1:2);
            if isempty(in{1})
                rt=obj.optical;
            elseif isa(in{1},'CPC')
                rt=RayTraceData(in{1});
            elseif isa(in{1},'RayTraceData')
                rt=in{1};
            else
                error('You can take a rest sitting there giving nothing. Otherwise, you just need to give me a CPC class, better with ray trace data (RayTraceData class).')
            end
            rt.geometry.checksimple;
            obj.optical=rt;
            [~,obj]=obj.getefficiency;
            Aper=obj.optical.geometry.Width*obj.optical.geometry.L;
            obj.operation.MassFR=Aper/20;
            [~,ang]=obj.outdoor.angles(obj.optical.geometry.beta);
            [beam,d1,d2]=obj.outdoor.tilt(obj.rho_grd,obj.optical.geometry.beta);
            d=d1+d2;
            tamb=obj.outdoor.DryBulb;
            vel=obj.outdoor.Wspd;
            id=abs(ang)<=89;
            T=obj.go(ang(id),beam(id),d(id),tamb(id),vel(id));
            T_o=cell2mat(T); T_o=T_o(4,:); T_in=obj.operation.T_in;
            E=obj.operation.MassFR*(obj.HTF.cp(T_o).*T_o...
                -obj.HTF.cp(T_in)*T_in);
            E=E(E>0);
            obj.energy=E;
            qdot=sum(E)/Aper;
            obj.intensity=qdot;
        end
        
        function [w,beam,diffuse_sky,diffuse_grd,hours]=binangle(obj,varargin)
            in=[varargin(:);cell(1,1)]; in=in(1);
            if isempty(in{1})
                iabin=0:89;
            else
                iabin=in{1};
            end
            nbin=numel(iabin);
            w=tmy3();
            w.code=obj.outdoor.code; w.station=obj.outdoor.station; w.state=obj.outdoor.state;
            w.latitude=obj.outdoor.latitude; w.longitude=obj.outdoor.longitude; w.timezone=obj.outdoor.timezone;
            w.doy=zeros(nbin,1);
            w.hod=zeros(nbin,1);
            w.GHI=zeros(nbin,1);
            w.DNI=zeros(nbin,1);
            w.dHI=zeros(nbin,1);
            w.DryBulb=zeros(nbin,1);
            w.Wspd=zeros(nbin,1);
            % Use regression to fit the model
            % q=a1*(B*eta_B_g(theta)+d*eta_d_g)+a2*(B*eta_B_abs(theta)+d*eta_d_abs)+a3*dT
            % Use the three conditions below:
            % 1. B=1000 theta=0 d=0 dT=T_in-mean(T_amb)
            % 2. B=0 theta=0(does not matter) d=1000 dT=T_in-mean(T_amb)
            % 3. B=800 theta=theta_c/2 d=200 dT=T_in-mean(T_amb)
            % 4. B=200 theta=(theta_c+90)/2 d=300 dT=T_in-mean(T_amb)
            temp=obj; 
            Vel=mean(temp.outdoor.Wspd)*ones(4,1);
            T_amb=mean(temp.outdoor.DryBulb)*ones(4,1);
            ang=[0;0;round(temp.optical.geometry.AccHalfAng/2);round((90+temp.optical.geometry.AccHalfAng)/2)];
            B=[1000;0;800;200]; d=[0;1000;200;300];
            [T,temp]=temp.go(ang,B,d,T_amb,Vel);
            T=cell2mat(T); T_o=T(4,:);
            T_in=temp.operation.T_in;
            q=temp.operation.MassFR*(temp.HTF.cp(T_o).*T_o...
                -temp.HTF.cp(T_in)*T_in);
            A=[temp.solarflux(1,:)',T_amb-T_in];
            c=(A'*A)\A'*q(:);
%             c=obj.getCoef;
            % Start to check each bin
            [~,ia2D]=obj.outdoor.angles(obj.optical.geometry.beta);
            [B,d1,d2]=obj.outdoor.tilt(obj.rho_grd,obj.optical.geometry.beta);
            Vel=obj.outdoor.Wspd; T_amb=obj.outdoor.DryBulb;
            eta=obj.efficiency;
            beam=zeros(nbin,1); 
            diffuse_sky=zeros(nbin,1);
            diffuse_grd=zeros(nbin,1);
            hours=zeros(nbin,1);
            for i=1:nbin
                if i==1
                    id=abs(ia2D)-iabin(i)<(iabin(i+1)-iabin(i))/2;
                elseif i==nbin
                    id=abs(ia2D)-iabin(i)>(iabin(i-1)-iabin(i))/2 & abs(ia2D)-iabin(i)<90-iabin(i);
                else
                    id=abs(ia2D)-iabin(i)<(iabin(i+1)-iabin(i))/2 & abs(ia2D)-iabin(i)>(iabin(i-1)-iabin(i))/2;
                end
                if any(id)
                    Bi=B(id); dskyi=d1(id); dgrdi=d2(id); di=dskyi+dgrdi;
                    veli=Vel(id); tambi=T_amb(id); ia=ia2D(id);
                    w_aper=obj.optical.aperture(ia); L=obj.optical.geometry.L;
                    a_aper=w_aper*L;
                    A=[(Bi.*eta(1,floor(abs(ia))+1)'+di*mean(eta(1,:))).*a_aper(:),...
                        tambi-T_in];
                    idcrit=A*c>0;
                    hours(i)=sum(idcrit);
                    if any(hours(i))
                        beam(i)=mean(Bi(idcrit));
                        diffuse_sky(i)=mean(dskyi(idcrit));
                        diffuse_grd(i)=mean(dgrdi(idcrit));
                        w.Wspd(i)=mean(veli(idcrit));
                        w.DryBulb(i)=mean(tambi(idcrit));
                    end
                end
            end
        end
        
        function [r,obj]=go(obj,ang,Beam,Diffuse,T_amb,Vel)
            % go(obj,angle,beam,diffuse,T_amb,Vel): simulate the condition where
            % incident angle, beam and diffuse are as specified
            obj.outdoor.Wspd=Vel(:);
            obj.outdoor.DryBulb=T_amb(:);
            if isempty(obj.VFA)
                [~,obj]=obj.vfa;
            end
            % Calculate solarflux
            q=obj.efficiency;
            ii=round(abs(ang))+1;
            L=obj.optical.geometry.L;
            q_sf=[q(:,ii)*L*diag(Beam)+mean(q,2)*L*Diffuse(:)';zeros(1,numel(ang))]; %W/m
            w_aper=obj.optical.aperture(ang);
            q_sf=q_sf*diag(w_aper); %W
            obj.solarflux=q_sf;
            % 
            obj.Ambient=air();
            obj.Ambient=obj.Ambient.setdefault(obj.outdoor.DryBulb);
            G=Beam+Diffuse;
            T_result=cell(1,numel(ang));
            options = optimoptions('fsolve','Jacobian','on','Display','off','MaxIter',20);
            dt=[5,50,500,5000];
            r_g=obj.optical.geometry.RadGlass;
            t_g=obj.optical.geometry.ThickGlass;
            r_a=obj.optical.geometry.RadReceiver;
            t_a=r_a/14;
            r_u=obj.optical.geometry.RadUtube;
            for j=1:numel(ang)
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
                    if exitflag~=1
                        fprintf('Equations not solved. Update guess values.\n')
                    else
                        fprintf('Equation solved.\n')
                    end
                end
                T_result{j}=T;
            end
            r=T_result;
            obj.temperature=r;
        end
        
        function [dv0,acc,rcr,q]=initdv(obj,lb,ub,varargin)
%             dv0=zeros(9,1);
            % fixed
            % receiver radius
            dv0(3)=ub(3);
            % glass radius
            dv0(4)=(ub(4)+lb(4))/2;
            % glass thickness
            dv0(5)=lb(5);
            % relative gap
            dv0(6)=lb(6);
            % u-tube radius
            dv0(7)=(ub(7)+lb(7))/2;
            % Length
            dv0(8)=ub(8);
            % Tilted angle
            dv0(9)=(lb(9)+ub(9))/2;
            % varied
            nacc=min(15,ceil((ub(1)-lb(1))/5)+1);
            acc=linspace(lb(1),ub(1),nacc);
            nrcr=min(11,ceil((ub(2)-lb(2))/0.1)+1);
            rcr=linspace(lb(2),ub(2),nrcr);
            q=zeros(nacc,nrcr);
            if isempty(varargin)
                nparfor=0;
            else
                nparfor=varargin{1};
            end
            parfor (i=1:nacc,nparfor)
                qtemp=zeros(1,nrcr);
                dv=dv0;
                dv(1)=acc(i);
                for j=1:nrcr
                    dv(2)=rcr(j);
                    try
                        qtemp(j)=obj.annualobfun(dv);
                    catch
                        warning('Possible errors occur at:')
                        dv
                    end
                end
                q(i,:)=qtemp;
            end
            [qmax,iacc]=max(q);
            [~,ircr]=max(qmax);
            iacc=iacc(ircr);
            dv0(1)=acc(iacc); dv0(2)=rcr(ircr);
        end
        
        function q=interpRT(obj,varargin)
            % interpRT(obj,'intercept' or 'efficiency',CPC class)
            in=[varargin(:);cell(2,1)]; in=in(1:2);
            if isempty(in{1})
                qi=9; % efficiency
            elseif strcmpi(in{1}(1:3),'int')
                qi=8; % intercept factor
            else
                qi=9;
            end
            if isempty(in{2})
                cpcgeo=obj.optical.geometry;
            else
                cpcgeo=CPC(in{2});
            end
            Data=load('RayTraceDataBase.mat');
            Data=Data.Data;
            a=size(Data,1); tol=1e-5;
            acc=cpcgeo.AccHalfAng;
            rcr=cpcgeo.RCR;
            rrg=cpcgeo.RadGlass/cpcgeo.RadReceiver;
            rtg=cpcgeo.ThickGlass/(rrg-1)/cpcgeo.RadReceiver;
            rgap=cpcgeo.Gap/cpcgeo.RadReceiver;
            dacc=cell2mat(Data(2:a,1)); udacc=sort(unique(dacc));
            drcr=cell2mat(Data(2:a,2)); udrcr=sort(unique(drcr));
            drrg=cell2mat(Data(2:a,4)); udrrg=sort(unique(drrg));
            drtg=cell2mat(Data(2:a,5)); udrtg=sort(unique(drtg));
            drgap=cell2mat(Data(2:a,7)); udrgap=sort(unique(drgap));
            if acc<udacc(1)-tol || acc>udacc(end)+tol
                error('Acceptance half angle out of bound.')
            end
            if rcr<udrcr(1)-tol || rcr>udrcr(end)+tol
                error('Relative concentration ratio out of bound.')
            end
            if rrg<udrrg(1)-tol || rrg>udrrg(end)+tol
                error('Glass radius out of bound.')
            end
            if rtg<udrtg(1)-tol || rtg>udrtg(end)+tol
                error('Relative glass thickness out of bound.')
            end
            if rgap<udrgap(1)-tol || rgap>udrgap(end)+tol
                error('Relative gap out of bound.')
            end
            coef=ones(a-1,1);
            % find each par
            temp=cell(5,3);
            temp(:,1)={acc;rcr;rrg;rtg;rgap};
            temp(:,2)={udacc;udrcr;udrrg;udrtg;udrgap};
            temp(:,3)={dacc;drcr;drrg;drtg;drgap};
            for i=1:5
                xq=temp{i,1};
                x=temp{i,2};
                ind=xq-x<=tol;
                ii=find(ind,1,'first');
                xub=x(ii);
                iub=abs(temp{i,3}-xub)<1e-10;
                if ii==1
                    coef(~iub)=0;
                else
                    xlb=x(ii-1);
                    cub=(xq-xlb)/(xub-xlb);
                    coef(iub)=coef(iub)*cub;
                    ilb=abs(temp{i,3}-xlb)<tol;
                    coef(ilb)=coef(ilb)*(1-cub);
                    coef(~ilb & ~iub)=0;
                end
            end
            if abs(sum(coef)-1)>tol
                fprintf('Super warning: the coefficients don''t sum to 1. Hmmmm, suspicous....\n')
                keyboard
            end
            Data=Data(2:a,:);
            Q=Data(coef>0,qi); Q=cell2mat(Q);
            n=sum(coef>0);
            row=repmat([1;2;3],1,n);
            col=1:3*n;
            cmat=repmat(coef(coef>0)',3,1);
            q=full(sparse(row,col,cmat)*Q);
        end
        
        function p=secretdb(obj,varargin)
            % make table
            % 'location','AccHalfAng','RCR','CR','r_a','rr_g','rr_ut','rt_g','rgap','length','beta','T_in','intensity','datetime(now)'
            ctgr={'location','AccHalfAng__DEG','RCR','CR','r_a__m','rr_g','rr_ut','rt_g','rgap','length__m','beta__DEG','T_in__K','intensity__Wh__m2','datetime__now'};
            nobj=numel(obj); ncat=numel(ctgr);
            result=cell(nobj,ncat);
            w=[obj.outdoor]';
            result(:,1)={w.station}';
            result(:,ncat)=repmat({now},nobj,1);
            geo=[obj.optical]'; geo=[geo.geometry]';
            data=zeros(nobj,ncat-2);
            data(:,1)=[geo.AccHalfAng]';
            data(:,2)=[geo.RCR]';
            data(:,3)=[geo.CR]';
            data(:,4)=[geo.RadReceiver]';
            data(:,5)=[geo.RadGlass]'./[geo.RadReceiver]';
            data(:,6)=[geo.RadUtube]'./[geo.RadReceiver]';
            data(:,7)=[geo.ThickGlass]'./([geo.RadGlass]'-[geo.RadReceiver]');
            data(:,8)=[geo.Gap]'./[geo.RadReceiver]';
            data(:,9)=[geo.L]';
            data(:,10)=[geo.beta]';
            op=[obj.operation]';
            data(:,11)=[op.T_in]';
            data(:,12)=[obj.intensity]';
            result(:,2:ncat-1)=num2cell(data);
            result=cell2table(result,'variablenames',ctgr);
            if isempty(varargin)
                secret=[num2str(now),'.txt'];
            else
                secret=varargin{1};
            end
            if exist(secret,'file')
                temp=readtable(secret);
                tempin=temp.intensity__Wh__m2;
                if iscell(tempin)
                    tempin=cellfun(@str2num,tempin);
                    temp.intensity__Wh__m2=tempin;
                end
                result=[result;temp];
            end
            writetable(result,secret)
            p=secret;
        end
        
        function h=binsummary(obj)
            f=figure;
            f.Position=[140,100,1000,600];
            % Cumulative total hour over ia
            ax1=subplot(2,2,1);
            cumhours=cumsum(obj.hours);
            [~,ang]=obj.outdoor.angles(obj.optical.geometry.beta);
            id=abs(ang)<90; ang=abs(ang(id));
            hcdf=cdfplot(ang);
            xx=get(hcdf,'xdata'); yy=get(hcdf,'ydata')*4380;
            delete(hcdf)
            plot(ax1,cumhours)
            hold on
            plot(ax1,xx,yy,'--')
            xlabel(ax1,'Incident angle ({\circ})')
            ylabel(ax1,'Total hours (hrs)')
            legend(ax1,'Useful hours','Available hours',...
                'location','northwest')
            xlim(ax1,[0,90])
            
            % Cumulative total beam and diffuse over ia
            ax2=subplot(2,2,2);
            [b,d1,d2]=obj.outdoor.tilt(obj.rho_grd,obj.optical.geometry.beta,id);
            d=d1+d2; [ang,i]=sort(ang); b=b(i); d=d(i);
            cumb=cumsum(b); cumd=cumsum(d);
            Beam=obj.other.beam.*obj.hours;
            cumbeam=cumsum(Beam);
            Diffuse=obj.other.diffuse.*obj.hours;
            cumdiffuse=cumsum(Diffuse);
            plot(ax2,cumbeam)
            hold on
            plot(ax2,cumdiffuse)
            plot(ax2,ang,cumb,'--')
            plot(ax2,ang,cumd,'--')
            xlabel(ax2,'Incident angle ({\circ})')
            ylabel(ax2,'Irradiance (Wh/m2)')
            legend(ax2,'Useful beam','Useful diffuse',...
                'Availabe beam','Available diffuse',...
                'location','northwest')
            xlim(ax2,[0,90])
            
            % Cumulative total absorbed solarflux over ia
            ax3=subplot(2,2,3);
            q=obj.efficiency;
            q_sf=q*diag(Beam)+mean(q,2)*Diffuse(:)'; %Wh/m2
            cumrec=cumsum(q_sf(1,:));
            cumref=cumsum(q_sf(2,:));
            cumglass=cumsum(q_sf(3,:));
            plot(ax3,cumrec)
            hold on
            plot(ax3,cumref)
            plot(ax3,cumglass)
            xlabel(ax3,'Incident angle ({\circ})')
            ylabel(ax3,'Solar flux (Wh/m2)')
            legend(ax3,'Solar flux absorbed by receiver',...
                'Solar flux absorbed by reflector',...
                'Solar flux absorbed by glass',...
                'location','northwest')
            xlim(ax3,[0,90])
            
            % Efficiency over ia
            ax4=subplot(2,2,4);
            plot(ax4,q')
            xlabel(ax4,'Incident angle ({\circ})')
            ylabel(ax4,'Optical efficiency')  
            legend('Receiver','Reflector','Glass')
            xlim(ax4,[0,90])
            ylim(ax4,[0,1])
            
            h=gcf;
            set(h,'color',[1 1 1])
        end
        
        function [varspace,q,h]=dvanalysis(obj,varargin)
            in=[varargin(:);cell(1,1)]; in=in(1);
            ub=[80,1,0.1,2,0.8,0.2,0.4,2,60];
            lb=[10,0,0.01,1.1,0.2,0,0.1,1,0];
            step=[5,0.1,0.01,0.1,0.1,0.05,0.1,0.25,5];
            vars={'Acceptance half angle ({\circ})',...
                'RCR','Receiver radius (m)','Relative glass radius',...
                'Relative glass thickness','Relative gap',...
                'Relative u-tube radius','Length (m)','Tilted angle ({\circ})'};
            dv0=OptimalDesign.cpc2dv(obj.optical.geometry);
            q0=obj.annualobfun(dv0);
            nvar=numel(ub);
            varspace=cell(nvar,1);
            ax=cell(nvar,1);
            q=cell(nvar,1);
            a=3; b=3;
            f=figure;
            f.Position=[140,0,1000,800];
            if isempty(in{1})
                npar=4;
            else
                npar=in{1};
            end
            parfor (i=1:nvar,npar)
                od=obj;
                dv=dv0;
                vartemp=[lb(i):step(i):ub(i),dv0(i)]';
                qtemp=zeros(size(vartemp)); qtemp(end)=q0;
                for j=1:numel(vartemp)-1
                    fprintf('%s: %d/%d\n',vars{i},j,numel(vartemp)-1)
                    dv(i)=vartemp(j);
                    qtemp(j)=real(od.annualobfun(dv));
                end
                [vartemp,id]=unique(vartemp);
                qtemp=qtemp(id);
                q{i}=qtemp;
                varspace{i}=vartemp;
            end
            qmax=max(cell2mat(q)); qmin=min(cell2mat(q));
            for i=1:nvar
                ax{i}=subplot(a,b,i);
                plot(ax{i},varspace{i},q{i},'o-','markerfacecolor',[1 1 1])
                hold on
                plot(ax{i},dv0(i),q0,'r^','markerfacecolor','r')
                xlabel(ax{i},vars{i})
                ylabel(ax{i},'Energy per aperture area (Wh/m^2)')
                ylim(ax{i},[qmin*0.9,qmax*1.1])
            end
            h=gcf;
        end
    end
    
    methods (Static) % Utility
        function fd=opticaldb(surp,sourcepath,targetpath)
            % surp: surface property
            % sourcepath: the path that stores the ray trace data
            % targetpath: the path where the database is saved
            if isempty(surp)
                surp='SurfaceProperty.mat';
            end
            list=dir(fullfile(sourcepath,'*.mat'));
            a=size(list,1);
            Data=cell(a+1,9);
            Data{1,1}='Acceptance half angle';
            Data{1,2}='Relative concentration ratio';
            Data{1,3}='Geometric concentration ratio';
            Data{1,4}='Radius of glass';
            Data{1,5}='Relative thickness of glass';
            Data{1,6}='Thickness of glass';
            Data{1,7}='Extra gap';
            Data{1,8}='Intercept factor';
            Data{1,9}='Optical efficiency';
            for i=2:a+1
                fprintf('%d/%d\n',i-1,a)
                temp=cell(1,9);
                tload=load(fullfile(sourcepath,list(i-1).name));
                rt=RayTraceData(tload.history);
                rt.geometry.CR=rt.geometry.CR-1e-4;
                rt.geometry=rt.geometry.update;
                temp{1}=round(rt.geometry.AccHalfAng);
                temp{2}=rt.geometry.RCR;
                temp{3}=rt.geometry.CR;
                temp{4}=rt.geometry.RadGlass;
                temp{5}=rt.geometry.ThickGlass/(rt.geometry.RadGlass-1);
                temp{6}=rt.geometry.ThickGlass;
                temp{7}=rt.geometry.Gap;
                temp{8}=rt.getintercept('total',0:89);
                temp{9}=rt.getefficiency('total',surp,[1 1 1],0:89);
                Data(i,:)=temp;
            end
            fd=fullfile(targetpath,'RayTraceDataBase.mat');
            save(fd,'Data')
        end
        
        function dv=cpc2dv(cpc)
            if isa(cpc,'CPC')
                cpc=cpc.update;
                dv(1)=cpc.AccHalfAng;
                dv(2)=cpc.RCR;
                dv(3)=cpc.RadReceiver;
                dv(4)=cpc.RadGlass/dv(3);
                dv(5)=cpc.ThickGlass/(cpc.RadGlass-dv(3));
                dv(6)=cpc.Gap/dv(3);
                dv(7)=cpc.RadUtube/dv(3);
                dv(8)=cpc.L;
                dv(9)=cpc.beta;
            else
                error('Good morning, my name is {\bold}CPC{\bold}2dv!')
            end
        end
        
        function cpc=dv2cpc(dv)
            if isa(dv,'numeric')
                cpc=CPC(dv(1),... % acceptance half angle
                    [],... % CR
                    dv(3),... % receiver radius
                    dv(3)*dv(4),... % glass radius
                    dv(3)*(dv(4)-1)*dv(5),... % glass thickness
                    dv(3)*dv(6),... % gap
                    dv(3)*dv(7));  % u-tube radius
                cpc.L=dv(8);
                cpc.beta=dv(9);
                [w1,w2]=cpc.minmaxwidth;
                cpc.CR=dv(2)*(w2/w1-1)+1;
                cpc.VacPres=1e-5;
            else
                error('Good morning, my name is cpc2{\bold}DV{\bold}!')
            end
        end
        
        function dv=table2dv(tbl)
            if isa(tbl,'table')
                dv=zeros(size(tbl,1),9);
                dv(:,1)=tbl.AccHalfAng__DEG;
                dv(:,2)=tbl.RCR;
                dv(:,3)=tbl.r_a__m;
                dv(:,4)=tbl.rr_g;
                dv(:,5)=tbl.rt_g;
                dv(:,6)=tbl.rgap;
                dv(:,7)=tbl.rr_ut;
                dv(:,8)=tbl.length__m;
                dv(:,9)=tbl.beta__DEG;
            else
                error('Good morning, my name is {\bold}TABLE{\bold}2dv!')
            end
        end
        
        function cpc=table2cpc(tbl)
            cpc=repmat(CPC(),size(tbl,1),1);
            dv=OptimalDesign.table2dv(tbl);
            for i=1:size(tbl,1)
                cpc(i)=OptimalDesign.dv2cpc(dv(i,:));
            end
        end
        
        function fd=intdb(sourcefilelist,tarpath)
            fd=tarpath;
            if isempty(sourcefilelist)
                return
            end
            r=cell(numel(sourcefilelist),1);
            for i=1:numel(sourcefilelist)
                r{i}=readtable(sourcefilelist(i).name);
                tempin=r{i}.intensity__Wh__m2;
                if iscell(tempin)
                    tempin=cellfun(@str2num,tempin);
                    r{i}.intensity__Wh__m2=tempin;
                end
            end
            r=vertcat(r{:});
            writetable(r,tarpath)
        end
        
        function grptbl=grouptable(tbl,varargin)
            % grouptable(tbl,'tin','length')
            r=size(tbl,1);
            loc=tbl.location;
            in=varargin;
            if isempty(in)
                in{1}='I AssUrE yOU nObOdY wOUld tYpE In thEse wOrds';
            end
            id=ismember({'tin','length'},in);
            scn=loc;n=2;
            if id(1)
                tin=num2str(tbl.T_in__K);
                tin=[repmat('T_in=',r,1),tin];
                tin=cellstr(tin);
                scn=cellfun(@strcat,scn,tin,'uniformoutput',false);
                n=n+1;
            end
            if id(2)
                len=num2str(tbl.length__m);
                len=[repmat('Length=',r,1),len];
                len=cellstr(len);
                scn=cellfun(@strcat,scn,len,'uniformoutput',false);
                n=n+1;
            end
            [scnuni,iall,iuni]=unique(scn);
            nscn=numel(scnuni);
            grptbl=cell(nscn,n);
            grptbl(:,1)=loc(iall);
            n=2;
            if id(1)
                grptbl(:,n)=tin(iall);
                n=n+1;
            end
            if id(2)
                grptbl(:,n)=len(iall);
                n=n+1;
            end
            for i=1:nscn
                id=iuni==i;
                temptbl=tbl(id,:);
                [~,isort]=sort(temptbl.intensity__Wh__m2,'descend');
                grptbl{i,n}=temptbl(isort,:);
            end
        end
        
        function optbl=findoptimal(tbl,varargin)
            grptbl=OptimalDesign.grouptable(tbl,'tin','length');
            grptbl=grptbl(:,4);
            f1=@(x) x(1,:);
            optbl=cellfun(f1,grptbl,'uniformoutput',false);
            optbl=vertcat(optbl{:});
            optbl=OptimalDesign.grouptable(optbl,varargin{:});
        end
    end
    
    methods (Access=protected, Hidden)
        function [qdot,obj]=annual_notsave(obj,varargin)
            % annual(obj,a class to specify a new geometry)
            in=[varargin(:);cell(1,1)]; in=in(1);
            if isempty(in{1})
                rt=obj.optical;
            elseif isa(in{1},'CPC')
                rt=RayTraceData(in{1});
            elseif isa(in{1},'RayTraceData')
                rt=in{1};
            else
                error('You can take a rest sitting there giving nothing. Otherwise, you just need to give me a CPC class, better with ray trace data (RayTraceData class).')
            end
            rt.geometry.checksimple;
            obj.optical=rt;
            bl=cellfun(@isempty,{obj.optical.x,obj.optical.y,obj.optical.index,obj.optical.ia});
            if any(bl) % do interpolation
                obj.optical.geometry=obj.optical.geometry.update;
                q=obj.interpRT;
            else
                q=obj.optical.getefficiency('total',[],[1 1 1],obj.optical.IncidentAngles); %W/m2
            end
            obj.efficiency=q;
            obj.other.weather=obj.outdoor; % store the original outdoor condition
            % Remember to update geometry before this line
            Aper=obj.optical.geometry.Width*obj.optical.geometry.L;
            obj.operation.MassFR=Aper/20;
            % Bin method
            iabin=0:89;
            [wth,beam,d1,d2,obj.hours]=obj.binangle(iabin);
            diffuse=d1+d2;
            % Assign a modified weather data to outdoor
            obj.outdoor=wth;
            % gogogo
            id=obj.hours~=0;
            Tr=obj.go(iabin(id),beam(id),diffuse(id),wth.DryBulb(id),wth.Wspd(id));
            T=cell(1,90);
            T(id)=Tr;
            T_in=obj.operation.T_in;
            T(~id)=mat2cell(repmat([300;300;T_in;T_in],1,sum(~id)),4,ones(1,sum(~id)));
            % Restore the outdoor
            obj.outdoor_bin=wth;
            obj.outdoor=obj.other.weather;
            obj.temperature=T;
            % Energy
            T_o=cell2mat(T); T_o=T_o(4,:);
            E=obj.operation.MassFR*(obj.HTF.cp(T_o).*T_o...
                -obj.HTF.cp(T_in)*T_in);
            obj.energy=E;
            qdot=sum(E)/Aper;
            obj.intensity=qdot;
        end
    end
end