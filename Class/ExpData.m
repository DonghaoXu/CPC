classdef ExpData
    properties
        Date;Time;Diffuse;Beam;Ambient;Wspd;
        TS1;TS2;TS3;TS4;MassFR;NoArrays;
    end
    
    properties (Dependent)
        Energy;AvailableSolar;Efficiency;
        dT;dTI;fB;ia3D;ia2D;ia2D_Long;
        AVGBeam;AVGDiffuse;AVGWspd;AVGEnergy;
        AVGdT;AVGdTI;AVGEff;AVGfB;AVGMassFR;
        AVGia3D;AVGia2D;AVGia2D_Long;AVGTS2;AVGTS3;effstats;
        AVGAmbient;AVGOperation;n;
    end
    
    properties (Access=protected)
        HTF=Duratherm600; cpcgeo;
        loc;
    end
    
    methods
        function y=ExpData(val)
%             par=[60,1.1748,0.027699,0.032572,0.001596,0.004780,...
%                 0.006,0.002,0.002,1.74,25,1e-5];
            y.cpcgeo=CPC(60,1.1748,0.027699,0.032572,...
                0.001596,0.004780,...
                0.006,0.002,0.002,1.738,25,1e-5);
            y.cpcgeo=y.cpcgeo.update;
            if nargin>0
                t=datenum(val.Date,'mmddyyyy');
                y.Date=datestr(t);
                y.Time=val.Time;
                y.Diffuse=val.DiffuseSolar;
                y.Beam=val.DirectSolar;
                y.Ambient=val.AmbientTemp+273;
                y.Wspd=val.WindVelocity*0.44704; %m/s
                y.TS1=val.TS1+273;
                y.TS2=val.TS2+273;
                y.TS3=val.TS3+273;
                y.TS4=val.TS4+273;
                y.NoArrays=val.NoArrays;
                y.MassFR=val.MassFlowRate/10/y.NoArrays; %kg/s per tube
                y.loc=tmy3;
                y.loc.latitude=40.4419;
                y.loc.longitude=-86.9125;
                y.loc.timezone=-75;
                y.loc.doy=repmat(datenum([y.Date(1:7),'0000']),numel(y.Time),1);
                y.loc.hod=y.Time;
            end
        end
        
        function y=get.n(obj)
            y=numel(obj.Time);
        end
        
        function y=get.ia3D(obj)
            y=obj.loc.angles(obj.cpcgeo.beta);
        end
        
        function y=get.ia2D(obj)
            [~,y]=obj.loc.angles(obj.cpcgeo.beta);
        end
        
        function y = get.ia2D_Long(obj)
            y = rad2deg(atan(sqrt(tan(deg2rad(obj.ia3D)).^2 -...
                tan(deg2rad(obj.ia2D)).^2)));
        end
        
        function y=get.AVGTS2(obj)
            y=mean(obj.TS2);
        end
        
        function y=get.AVGTS3(obj)
            y=mean(obj.TS3);
        end
        
        function y=get.Energy(obj)
            cp=(obj.HTF.cp(obj.TS2)+obj.HTF.cp(obj.TS2))/2;
            y=obj.MassFR.*cp.*(obj.TS3-obj.TS2);
        end
        
        function y=get.AVGEnergy(obj)
            y=mean(obj.Energy);
        end
        
        function y=get.AvailableSolar(obj)
            aper=obj.cpcgeo.Width*obj.cpcgeo.L;
            y=(obj.Diffuse+obj.Beam)*aper;
        end
        
        function y=get.Efficiency(obj)
            y=obj.Energy./obj.AvailableSolar;
        end
        
        function y=get.effstats(obj)
            y1=(max(obj.Efficiency)-min(obj.Efficiency))/2;
            y2=mean(obj.Efficiency);
            y3=std(obj.Efficiency);
            y=[y2;y3;y1];
        end
        
        function y=get.AVGEff(obj)
            y=mean(obj.Energy)/mean(obj.AvailableSolar);
        end
        
        function y=get.dT(obj)
            y=(obj.TS2+obj.TS3)/2-obj.Ambient;
        end
        
        function y=get.AVGdT(obj)
            y=mean(obj.dT);
        end
        
        function y=get.dTI(obj)
            y=obj.dT./(obj.Diffuse+obj.Beam);
        end
        
        function y=get.AVGdTI(obj)
            y=obj.AVGdT/(obj.AVGDiffuse+obj.AVGBeam);
        end
        
        function y=get.AVGBeam(obj)
            y=mean(obj.Beam);
        end
        
        function y=get.AVGDiffuse(obj)
            y=mean(obj.Diffuse);
        end
        
        function y=get.fB(obj)
            y=obj.Beam./(obj.Beam+obj.Diffuse);
        end
        
        function y=get.AVGfB(obj)
            y=obj.AVGBeam./(obj.AVGBeam+obj.AVGDiffuse);
        end
        
        function y=get.AVGWspd(obj)
            y=mean(obj.Wspd);
        end
        
        function y=get.AVGMassFR(obj)
            y=mean(obj.MassFR);
        end
        
        function y=get.AVGia3D(obj)
            y=mean(obj.ia3D);
        end
        
        function y=get.AVGia2D(obj)
            y=mean(obj.ia2D);
        end
        
        function y = get.AVGia2D_Long(obj)
            y = mean(obj.ia2D_Long);
        end
        
        function y=get.AVGAmbient(obj)
            y=mean(obj.Ambient);
        end
        
        function y=get.AVGOperation(obj)
            y.MassFR=obj.AVGMassFR;
            y.T_in=mean(obj.TS2);
        end
        
        function [bl,ind,obj]=singlecheck(obj)
            ind=[obj.Diffuse]'>0;
            ind=ind & [obj.Beam]'>0;
            ind=ind & [obj.Ambient]'>0;
            ind=ind & [obj.TS2]'>0;
            ind=ind & [obj.TS3]'>0;
            ind=ind & [obj.MassFR]'>0;
            bl=~any(~ind);
        end
        
        function c=getcpc(obj)
            c=obj.cpcgeo;
        end
    end
    methods
        function [obj,ind]=arraycheck(obj)
            f=@(x) x.singlecheck;
            ind=arrayfun(f,obj);
            obj=obj(ind);
        end
        
        function [mdl,tbl]=regdata(obj,varargin)
            in=[varargin(:);cell(1,1)]; in=in(1);
            if isempty(in{1})
                in{1}='1';
            end
            varname={'fB','fd','dTI1','dTI2','Eff'};
            B=[obj.AVGBeam]'; d=[obj.AVGDiffuse]';
            tbl=table(B./(B+d),d./(B+d),[obj.AVGdTI]',...
                [obj.AVGdT]'.^2./(B+d),...
                [obj.AVGEff]',...
                'variablenames',varname);
            switch in{1}
                case {'second','2'}
                    mdl=fitlm(tbl,'Eff~fB+fd+dTI1+dTI2','intercept',false,'RobustOpts','on');
                otherwise
                    mdl=fitlm(tbl,'Eff~fB+fd+dTI1','intercept',false,'RobustOpts','on');
            end
            ExpData.plotmdl(mdl)
        end
        
        function y=filter(obj,varargin)
            ind=true(numel(obj),1);
            if nargin>0
                for i=1:numel(varargin)
                    ind=ind & varargin{i};
                end
            end
            y=obj(ind);
            y.histsum;
        end
        
        function y=plot(obj,xvar,yvar,varargin)
            unc=obj.uct_measure;
            switch lower(xvar)
                case 'ia'
                    x=abs([obj.AVGia2D]');
                    xlb=0;
                    xub=0;
                case 'fb'
                    x=[obj.AVGfB]';
                    xlb=x-unc.Beam(:,1)./(unc.Beam(:,1)+unc.Diffuse(:,3));
                    xub=unc.Beam(:,3)./(unc.Beam(:,3)+unc.Diffuse(:,1))-x;
                case 'dti'
                    x=[obj.AVGdTI]';
                    xlb=x-unc.dTI(:,1);
                    xub=unc.dTI(:,3)-x;
                case 'dt'
                    x=[obj.AVGdT]';
                    xlb=x-unc.dT(:,1);
                    xub=unc.dT(:,3)-x;
            end
            switch lower(yvar)
                case 'eff'
                    y=[obj.AVGEff]';
                    ylb=y-unc.Efficiency(:,1);
                    yub=unc.Efficiency(:,3)-y;
                case 'energy'
                    y=[obj.AVGEnergy]';
                    ylb=0;
                    yub=0;
            end
%             figure
            
            % y uncertainty
            if any(yub) && nargin>3
                y = errorbar(x,y,ylb,yub,'o');
            else
                y = plot(x,y,'o');
            end
            % x uncertainty
%             if any(xu)
%                 h2=plot(xu(:,[1,1]),y)
%             end
        end
        
        function [sm,f]=histsum(obj)
            varname={'f_B','Global solar radiation (W/m^2)',...
                '{\Delta}T ({\circ}C)','Incident angle ({\circ})'};
            nvar=numel(varname);
            data={[obj.AVGfB]',...
                [obj.AVGBeam]'+[obj.AVGDiffuse]',...
                [obj.AVGdT]',...
                [obj.AVGia2D]',[obj.AVGdTI]',[obj.AVGTS3]'-273};
            ax=cell(nvar,1);
            a=2; b=2;
            f=figure;
            f.Position=[190,100,900,600];
            for i=1:nvar
                ax{i}=subplot(a,b,i);
                histogram(ax{i},data{i},10)
                xlabel(ax{i},varname{i})
            end
            f=gcf;
            mm=[cellfun(@min,data);cellfun(@max,data)];
            varname={'f_B','I_g','dT','ia','dti','T'};
            mm=array2table(mm,'variablenames',varname);
            mt=table({'min';'max'},'variablenames',{'Bound'});
            sm=[mt,mm];
        end
        
        function [mdl,tbl,robj]=bestreg(obj)
            [mdl,tbl]=obj.regdata;
            robj=obj;
%             p=mdl.NumPredictors;
            while 1
%                 n=mdl.NumObservations;
                ind=abs(mdl.Residuals.Studentized)<3;
%                 ind=mdl.Diagnostics.Leverage<p*2/n;
%                 ind=abs(mdl.Residuals.Raw)<0.05;
                if ~any(~ind)
                    break
                end
                robj=robj(ind);
                tbl=tbl(ind,:);
                mdl=fitlm(tbl,'Eff~fB+fd+dTI1','intercept',false,'RobustOpts','on');
            end
            ExpData.plotmdl(mdl);
            robj.histsum;
        end
        
        function [mdl,tbl,robj]=bestangreg(obj,varargin)
            robj=obj;
%             p=mdl.NumPredictors;
            while 1
%                 n=mdl.NumObservations;
                [mdl,tbl,robj]=robj.angregdata(varargin{:});
                ind=abs(mdl.Residuals.Studentized)<3;
%                 ind=mdl.Diagnostics.Leverage<p*2/n;
%                 ind=abs(mdl.Residuals.Raw)<0.05;
                if ~any(~ind)
                    break
                end
                robj=robj(ind);
            end
            ExpData.plotmdl(mdl);
            robj.histsum;
        end
        
        function [mdl,tbl,robj,sind]=angregdata(obj,varargin)
            in=[varargin(:);cell(1,1)]; in=in(1);
            if isempty(in{1})
                in{1}='1';
            end
            ia=abs([obj.AVGia2D]');
            fb=[obj.AVGfB]';
            dti=[obj.AVGdTI]';
            dti2=[obj.AVGdT]'.^2./([obj.AVGBeam]'+[obj.AVGDiffuse]');
            eff=[obj.AVGEff]';
            edges=[0:5:60,68,75,90];
            figure;
            h=histogram(ia,edges);
            xlabel('Incident angle ({\circ})')
            id=find(h.Values>=10);
            nid=numel(id);
%             data=cell(nid,1);
            varname=cell(nid+4,1);
            sind=false(numel(obj),1);
            temp=zeros(numel(obj),nid+4);
            for i=1:nid
                ind=ia>=edges(i) & ia<=edges(i+1);
                sind=sind | ind;
                temp(ind,i)=fb(ind);
                temp(ind,nid+1)=1-fb(ind);
                temp(ind,nid+2)=dti(ind);
                temp(ind,nid+3)=dti2(ind);
                temp(ind,nid+4)=eff(ind);
                varname{i}=sprintf('fb_%d_%d',edges(id(i)),edges(id(i)+1));
%                 data{i}=temp;
            end
            robj=obj(sind);
            temp=temp(temp(:,end)~=0,:);
%             data=vertcat(data{:});
            varname(nid+1:nid+4)={'fd';'dTI1';'dTI2';'Eff'};
            tbl=array2table(temp,'variablenames',varname);
            switch in{1}
                case {'2','second'}
                    mdl=fitlm(tbl,'linear','intercept',false,'robustopts','on');
                otherwise
                    mdl=fitlm(tbl(:,[1:nid+2,nid+4]),'linear','intercept',false,'robustopts','on');
            end
            ExpData.plotmdl(mdl);
        end
        
        function ind=condition(obj,var1,var2,varargin)
            in=[varargin(:);cell(3,1)]; in=in(1:3);
            if isempty(in{1})
                nbin1=10;
            else
                nbin1=in{1};
            end
            if isempty(in{2})
                nbin2=10;
            else
                nbin2=in{2};
            end
            if isempty(in{3})
                k=1;
            else
                k=in{3};
            end
            switch var1
                case 'dti'
                    x1=[obj.AVGdTI]';
                case 'dt'
                    x1=[obj.AVGdT]';
                case 'Ig'
                    x1=[obj.AVGBeam]'+[obj.AVGDiffuse]';
                case 'beam'
                    x1=[obj.AVGBeam]';
                case 'diffuse'
                    x1=[obj.AVGDiffuse]';
                case 'fb'
                    x1=[obj.AVGfB]';
                case 'ia'
                    x1=[obj.AVGia2D]';
                otherwise
                    error('Filter by whwat?')
            end
            switch var2
                case 'dti'
                    x2=[obj.AVGdTI]';
                case 'dt'
                    x2=[obj.AVGdT]';
                case 'Ig'
                    x2=[obj.AVGBeam]'+[obj.AVGDiffuse]';
                case 'beam'
                    x2=[obj.AVGBeam]';
                case 'diffuse'
                    x2=[obj.AVGDiffuse]';
                case 'fb'
                    x2=[obj.AVGfB]';
                case 'ia'
                    x2=[obj.AVGia2D]';
            end
            x=[x1,x2];
            h=hist3(x,[nbin1,nbin2]);
            [~,id]=sort(h(:),'descend');
            k=min(k,numel(id));
            id=id(k);
            i2=ceil(id/nbin2); i1=id-(i2-1)*nbin1;
            b1=min(x1)+(i1-[1,0])/nbin1*(max(x1)-min(x1));
            b2=min(x2)+(i2-[1,0])/nbin2*(max(x2)-min(x2));
            ind=x1>=b1(1) & x1<=b1(2) & x2>=b2(1) & x2<=b2(2);
        end
        
        function y=uct_measure(obj)
            uct_massfr=0.001;
            uct_temp=0.3; %K
            uct_solar=0.08;
            sqrtn=sqrt([obj.n]');
            effub=[obj.AVGEff]'.*(1+uct_massfr./sqrtn).*(1+uct_temp./sqrtn./([obj.AVGTS3]'-[obj.AVGTS2]'))./(1-uct_solar./sqrtn);
            efflb=[obj.AVGEff]'.*(1-uct_massfr./sqrtn).*(1-uct_temp./sqrtn./([obj.AVGTS3]'-[obj.AVGTS2]'))./(1+uct_solar./sqrtn);
            meff=[obj.AVGEff]';
            dtiub=[obj.AVGdTI]'.*(1+uct_temp./sqrtn./[obj.AVGdT]')./(1-uct_solar./sqrtn);
            dtilb=[obj.AVGdTI]'.*(1-uct_temp./sqrtn./[obj.AVGdT]')./(1+uct_solar./sqrtn);
            mdti=[obj.AVGdTI]';
            mdt=[obj.AVGdT]';
            dtub=[obj.AVGdT]'.*(1+uct_temp./sqrtn./[obj.AVGdT]');
            dtlb=[obj.AVGdT]'.*(1-uct_temp./sqrtn./[obj.AVGdT]');
            mtinc=[obj.AVGTS3]'-[obj.AVGTS2]';
            tincub=mtinc.*(1+uct_temp./sqrtn./mtinc);
            tinclb=mtinc.*(1-uct_temp./sqrtn./mtinc);
            mb=[obj.AVGBeam]';
            bub=mb.*(1+uct_solar./sqrtn);
            blb=mb.*(1-uct_solar./sqrtn);
            md=[obj.AVGDiffuse]';
            dub=md.*(1+uct_solar./sqrtn);
            dlb=md.*(1-uct_solar./sqrtn);
            y=table([efflb,meff,effub],[dtilb,mdti,dtiub],[dtlb,mdt,dtub],...
                [tinclb,mtinc,tincub],[blb,mb,bub],[dlb,md,dub],...
                'variablenames',{'Efficiency','dTI','dT','TempInc','Beam','Diffuse'});
        end
        
        function r=uct_method(obj,mdl)
            stats=[obj.effstats];
            err_method=stats(2,:);
            sst=mdl.SST;
            r=sum(err_method.^2)/sst;
            histogram(err_method,10)
        end
        
        function [r,dtexp,dtmdl]=validate(obj,varargin)
            in=[varargin(:);cell(1,1)]; in=in(1);
            c=obj.getcpc;
            if any(strcmpi(in{1},{'Numeric','num'}))
                thm=NumericThermal(c);
            else
                thm=SimpleThermal(c);
            end
            thm.optical=thm.optical.update;
            [~,thm]=thm.vfa;
            ang=[obj.AVGia2D]';
            beam=[obj.AVGBeam]';
            diffuse=[obj.AVGDiffuse]';
            T_amb=[obj.AVGAmbient]';
            vel=[obj.AVGWspd]';
            op=[obj.AVGOperation]';
            r=thm.go(ang,beam,diffuse,T_amb,vel,op);
            f=@(x) x(end,1);
            to=cellfun(f,r);
            ti=[obj.AVGTS2];
            dtexp=[obj.AVGTS3]-ti;
            dtmdl=to-ti;
        end
        
        function [d,r,dtexp,dtmdl]=findvacpres(obj,p,varargin)
            % varargin: contains a string to determine the use of numerical
            % or simple
            in=[varargin(:);cell(1,1)]; in=in(1);
            np=numel(p);
            r=cell(np,1);
            dtexp=cell(np,1);
            dtmdl=cell(np,1);
            d=zeros(np,1);
            ang=[obj.AVGia2D]';
            beam=[obj.AVGBeam]';
            diffuse=[obj.AVGDiffuse]';
            T_amb=[obj.AVGAmbient]';
            vel=[obj.AVGWspd]';
            op=[obj.AVGOperation]';
            ti=[obj.AVGTS2];
            dtexpi=[obj.AVGTS3]-ti;
            f=@(x) x(end,1);
            c=obj.getcpc;
            if any(strcmpi(in{1},{'Numeric','num'}))
                thm=NumericThermal(c);
            else
                thm=SimpleThermal(c);
            end
            thm.optical=thm.optical.update;
            [~,thm]=thm.vfa;
            parfor i=1:np
                temp=thm;
                temp.optical.geometry.VacPres=p(i);
                ri=temp.go(ang,beam,diffuse,T_amb,vel,op);
                to=cellfun(f,ri);
                dtmdli=to-ti;
                r{i}=ri;
                dtexp{i}=dtexpi;
                dtmdl{i}=dtmdli;
                d(i)=sum((dtexpi-dtmdli).^2);
            end
        end
    end
    methods (Static)
        function plotmdl(mdl)
            figure;
            ax=axes;
            cmap=mycmap();
            expc=cmap{2}(2,:);
            plot(ax,mdl.Variables.dTI1,mdl.Variables.Eff,'ko','markerfacecolor',expc)
            hold on
            regc=cmap{4}(3,:);
            plot(ax,mdl.Variables.dTI1,mdl.Fitted,'ks','markerfacecolor',regc)
            xlabel('{\Delta}T/I_g (m^2-{\circ}C/W)')
            ylabel('Overall efficiency')
            legend('Experiment data','Fitting data')
            figure;
%             ax1=subplot(1,2,1);
%             copyobj(allchild(ax),ax1);
%             xlabel('{\Delta}T/I_g (m^2-{\circ}C/W)')
%             ylabel('Overall efficiency')
%             legend('Experiment data','Fitting data')
% %             ax1.Position=[0.05,0.1,0.4,0.8];
%             ax2=subplot(1,2,2);
            ax2=axes;
            plot(mdl.Variables.Eff,mdl.Fitted,'o')
            hold on
            lmt=[ax2.XLim;ax2.YLim];
            rg=[min(lmt(:,1)),max(lmt(:,2))];
            plot(rg,rg)
            axis equal
            xlabel('Experimental data')
            ylabel('Fitting data')
%             ax2.Position=[0.55,0.1,0.4,0.8];
        end
        
        function plotresid(mdl)
            ax1=subplot(1,2,1);
            mdl.plotResiduals('fitted');
            ax2=subplot(1,2,2);
            mdl.plotResiduals('probability')
        end
        
        function Data_raw=extractdata(fname)
            raw=xlsread(fname,2);
            Data_raw.Time=raw(:,1);
            Data_raw.OAT=raw(:,2);
            Data_raw.OAH=raw(:,3);
            Data_raw.Global=raw(:,4);
            Data_raw.Wspd=raw(:,5);
            Data_raw.FR=raw(:,6);
            Data_raw.TS1=raw(:,8);
            Data_raw.TS2=raw(:,9);
            Data_raw.TS3=raw(:,10);
            Data_raw.TS4=raw(:,11);
            Data_raw.Total=raw(:,14);
            Data_raw.Diffuse=raw(:,15);
            if size(raw,2)>24
                Data_raw.Surf1=raw(:,25);
                Data_raw.Surf2=raw(:,26);
                Data_raw.Surf3=raw(:,27);
            end
        end
        
        function sep232013
            fname='report_09232013_Analysis.xlsx';
            Data_raw=ExpData.extractdata(fname);
            Data_raw.Beam=Data_raw.Total-Data_raw.Diffuse;
            dtbl=struct2table(Data_raw);
            ind=dtbl.FR>0;
            dtbl=dtbl(ind,:);
            n=datenum('Sep-23-0000');
            beta=25; gama=0;
            [ia3D,~,alpha,phi,~,~,dtbl.SolarTime]=IncAngle(n,dtbl.Time*24,40.44,-86.91,-75,beta,gama,1);
            ia2D=IncAng2D(alpha,phi,beta);
            dtbl.SolarTime=dtbl.SolarTime/24;
            load('Exp2013')
            d=d2013(datenum(vertcat(d2013.Date))==datenum('Sep-23-2013'));
            tm=vertcat(d.Time)/24;
            startend=arrayfun(@(x) numel(x.Time),d);
            startend=[1,cumsum(startend)];
            idsteady=dtbl.Time>=tm(1) & dtbl.Time<=tm(end);
            c=d.getcpc;
            Area=c.Width*c.L;
            htf=Duratherm600;
            cp_htf=htf.cp((dtbl.TS3+dtbl.TS2)/2+273);
            rho_htf=htf.rho((dtbl.TS3+dtbl.TS2)/2+273);
            dtbl.eff=dtbl.FR/4/10*3.78/60.*rho_htf.*cp_htf.*(dtbl.TS3-dtbl.TS2)...
                ./(dtbl.Total*Area*1000);
            cmap=mycmap();
            % Figure 1
            figure; hold on
            ax1=gca;
            plot(ax1,dtbl.Time,dtbl.Total,'--','color',[1 1 1]*0.3)
            h(1)=plot(ax1,dtbl.Time(idsteady),dtbl.Total(idsteady),...
                'o-','color',cmap{2}(3,:));
            plot(ax1,dtbl.Time,dtbl.Beam,'--','color',[1 1 1]*0.3);
            h(2)=plot(ax1,dtbl.Time(idsteady),dtbl.Beam(idsteady),...
                '^-','color',cmap{3}(3,:));
            ax1.Position(4)=0.8;
            ax1_pos=ax1.Position;
            ax2=axes('Position',ax1_pos,...
                'YAxisLocation','right',...
                'XAxisLocation','top',...
                'color','none');
            ax2.NextPlot='add';
            plot(ax2,dtbl.SolarTime,dtbl.TS2,'--','color',[1 1 1]*0.3);
            h(3)=plot(ax2,dtbl.SolarTime(idsteady),dtbl.TS2(idsteady),...
                's-','color',cmap{5}(3,:));
            plot(ax2,dtbl.SolarTime,dtbl.OAT,'--','color',[1 1 1]*0.3);
            h(4)=plot(ax2,dtbl.SolarTime(idsteady),dtbl.OAT(idsteady),...
                'd-','color',cmap{1}(3,:));
            ax2.YLim=ax2.YLim*2;
            datetick(ax1,'x',16)
            [~,~,~,~,~,~,ast]=IncAngle(n,ax1.XTick*24,40.44,-86.91,-75,[],[],1);
            ax2.XTick=ast/24;
            ax2.XLim=ast([1,end])/24;
            ax2.XTickLabel=datestr(ax2.XTick);
            ylabel(ax1,'Solar irradiance (W/m^2)')
            ylabel(ax2,'Temperature ({\circ}C)')
            plot(ax1,[tm(startend),tm(startend)],ax1.YLim,':','color',[1 1 1]*0.3)
            t1=text(tm(1),420,'Quasi-steady states\rightarrow',...
                'HorizontalAlignment','right','parent',ax1);
            xlabel(ax1,'Local time')
            xlabel(ax2,'Solar time')
            l=legend(h,'Global radiation','Beam radiation',...
                'Inlet temperature','Ambient temperature',...
                'location','northwest');
            FigureFormat(5,10/3,7)
            % Figure 2
            clear h
            figure; hold on
            ax1=gca;
            h(1)=plot(ax1,dtbl.Time(idsteady),ia3D(idsteady),...
                '-o','color',cmap{1}(3,:));
            h(2)=plot(ax1,dtbl.Time(idsteady),ia2D(idsteady),...
                '-s','color',cmap{2}(3,:));
            ax1.Position(4)=0.8;
            ax1_pos=ax1.Position;
            ax2=axes('Position',ax1_pos,...
                'YAxisLocation','right',...
                'XAxisLocation','top',...
                'color','none');
            ax2.NextPlot='add';
            h(3)=plot(ax2,dtbl.SolarTime(idsteady),dtbl.Beam(idsteady)./dtbl.Total(idsteady),...
                '-^','color',cmap{3}(3,:));
            ax2.YLim=[0,1];
            ax1.YLim(2)=ax1.YLim(2)+10;
            ax1.XTickLabel=datestr(ax1.XTick);
            [~,~,~,~,~,~,ast]=IncAngle(n,ax1.XTick*24,40.44,-86.91,-75,[],[],1);
            ax2.XTick=ast/24;
            ax2.XLim=ast([1,end])/24;
            ax2.XTickLabel=datestr(ax2.XTick);
            plot(ax1,[tm(startend),tm(startend)],ax1.YLim,':','color',[1 1 1]*0.3)
            xlabel(ax1,'Local time')
            xlabel(ax2,'Solar time')
            ylabel(ax1,'Angle ({\circ})')
            ylabel(ax2,'Beam fraction')
            legend(h,'Incident angle','Traversal projection','Beam fraction',...
                'location','southeast')
            FigureFormat(5,10/3,7)
            % Figure 3
            clear h
            figure; hold on
            ax1=gca;
            h(1)=plot(ax1,dtbl.Time(idsteady),dtbl.TS3(idsteady)-dtbl.TS2(idsteady),...
                '-o','color',cmap{1}(3,:));
            h(2)=plot(ax1,dtbl.Time(idsteady),(dtbl.TS3(idsteady)+dtbl.TS2(idsteady))/2-dtbl.OAT(idsteady)-130,...
                '-s','color',cmap{2}(3,:));
            ax1.Position(4)=0.8;
            ax1_pos=ax1.Position;
            ax2=axes('Position',ax1_pos,...
                'YAxisLocation','right',...
                'XAxisLocation','top',...
                'color','none');
            ax2.NextPlot='add';
            h(3)=plot(ax2,dtbl.SolarTime(idsteady),dtbl.eff(idsteady),...
                '^','color',cmap{3}(3,:));
            ax2.YLim=[0,1];
            ax1.YLim=[-15,20];
            ax1.YTick=-15:5:20;
            ax1.YTickLabel={'-15','-10','-5','0','(135)5','140','145','150'};
            ax1.XTickLabel=datestr(ax1.XTick);
            [~,~,~,~,~,~,ast]=IncAngle(n,ax1.XTick*24,40.44,-86.91,-75,[],[],1);
            ax2.XTick=ast/24;
            ax2.XLim=ast([1,end])/24;
            ax2.XTickLabel=datestr(ax2.XTick);
            plot(ax1,[tm(startend),tm(startend)],ax1.YLim,':','color',[1 1 1]*0.3)
            xlabel(ax1,'Local time')
            xlabel(ax2,'Solar time')
            ylabel(ax1,'Temperature ({\circ}C)')
            ylabel(ax2,'Overall efficiency')
            legend(h,'T_{out}-T_{in}','{\Delta}T=(T_{in}+T_{out})/2-T_{amb}','{\eta}',...
                'location','southeast')
            FigureFormat(5,10/3,7)
            
            
        end
    end
    methods
        % Specific functions NOT general
        function [h,sm]=effvsdt(obj)
            h=cell(4,1);
            c=mycmap;
            h{1}=obj.plot('dt','eff'); hold on
            ind=obj.condition('Ig','fb');
            dt=obj(ind);
            ind=dt.condition('ia','ia',4,4);
            dt=dt(ind);
            dt=dt([dt.AVGdT]>15);
            h{2}=dt.plot('dt','eff');
            ylim([0,1])
            xlabel('{\Delta}T ({\circ}C)')
            ylabel('Overall efficiency')
            x=h{2}.XData(:); y=h{2}.YData(:);
            xlim([min(x)-5,max(x)+5])
            [x,i]=sort(x); y=y(i);
            mdl=fitlm(x,y,'robustopts','on');
            h{3}=plot(x([1,end]),mdl.Fitted([1,end]),'-','color',c{5}(3,:));
            mdl=fitlm(x,y,'quadratic','robustopts','on');
            xquad=linspace(min(x),max(x),30);
            h{4}=plot(xquad,predict(mdl,xquad(:)),'--','color',c{4}(3,:));
            legend('Unselected data points',...
                'Selected data points',...
                'Linear fitting','Quadratic fitting')
            f=gcf;
            sm=dt.histsum;
            close(gcf)
            str={'f_B: 0.84~0.91',...
                'I_g: 1000~1094W/m^2',...
                'Incident angle: 16~26{\circ}'};
            annotation(f,'textbox',[0.15,0.8,0.1,0.1],...
                'String',str);
            FigureFormat(5,10/3,7)
            h{1}.Marker='.';
            h{1}.MarkerEdgeColor=[1 1 1]*0.7;
            h{2}.MarkerEdgeColor='k';%c{1}(4,:);
            h{2}.MarkerFaceColor=c{1}(2,:);
            h{2}.MarkerSize=4;
            print('3.4vsdT.tif','-dtiff')
        end
        
        function [h,sm]=effvsfb(obj)
            h=cell(4,1);
            c=mycmap;
            h{1}=obj.plot('fb','eff'); hold on
            ind=obj.condition('dti','ia',[],[],5);
            dt=obj(ind);
            h{2}=dt.plot('fb','eff');
            ylim([0,1])
            xlabel('f_B')
            ylabel('Overall efficiency')
            x=h{2}.XData(:); y=h{2}.YData(:);
            xlim([min(x)-(max(x)-min(x))*0.05,max(x)+(max(x)-min(x))*0.05])
            [x,i]=sort(x); y=y(i);
            mdl=fitlm(x,y,'robustopts','on');
            h{3}=plot(x([1,end]),mdl.Fitted([1,end]),'-','color',c{5}(3,:));
            mdl=fitlm(x,y,'quadratic','robustopts','on');
            xquad=linspace(min(x),max(x),30);
            h{4}=plot(xquad,predict(mdl,xquad(:)),'--','color',c{4}(3,:));
            legend('Unselected data points',...
                'Selected data points',...
                'Linear fitting','Quadratic fitting')
            f=gcf;
            sm=dt.histsum;
            close(gcf)
            str={'{\Delta}T/I_g: 0.06~0.07{\circ}C-m^2/W',...
                'Incident angle: 27~36{\circ}'};
            annotation(f,'textbox',[0.15,0.8,0.1,0.1],...
                'String',str);
            FigureFormat(5,10/3,7)
            h{1}.Marker='.';
            h{1}.MarkerEdgeColor=[1 1 1]*0.7;
            h{2}.MarkerEdgeColor='k';%c{1}(4,:);
            h{2}.MarkerFaceColor=c{1}(2,:);
            h{2}.MarkerSize=4;
            print('3.4vsfB.tif','-dtiff')
        end
        
        function [h,sm]=effvsia(obj)
            h=cell(4,1);
            c=mycmap;
            h{1}=obj.plot('dt','eff'); hold on
            ind=obj.condition('fb','dti',[],[],1);
            dt=obj(ind);
            h{2}=dt.plot('ia','eff');
            ylim([0,1])
            xlabel('Incident angle ({\circ})')
            ylabel('Overall efficiency')
            x=h{2}.XData(:); y=h{2}.YData(:);
            xlim([min(x)-(max(x)-min(x))*0.05,max(x)+(max(x)-min(x))*0.05])
            [x,i]=sort(x); y=y(i);
            mdl=fitlm(x,y,'robustopts','on');
            h{3}=plot(x([1,end]),mdl.Fitted([1,end]),'-','color',c{5}(3,:));
            mdl=fitlm(x,y,'quadratic','robustopts','on');
            xquad=linspace(min(x),max(x),30);
            h{4}=plot(xquad,predict(mdl,xquad(:)),'--','color',c{4}(3,:));
            legend('Unselected data points',...
                'Selected data points',...
                'Linear fitting','Quadratic fitting')
            f=gcf;
            sm=dt.histsum;
            close(gcf)
            str={'{\Delta}T/I_g: 0.120~0.143{\circ}C-m^2/W',...
                'f_B: 0.83~0.92'};
            annotation(f,'textbox',[0.15,0.8,0.1,0.1],...
                'String',str);
            FigureFormat(5,10/3,7)
            h{1}.Marker='.';
            h{1}.MarkerEdgeColor=[1 1 1]*0.7;
            h{2}.MarkerEdgeColor='k';%c{1}(4,:);
            h{2}.MarkerFaceColor=c{1}(2,:);
            h{2}.MarkerSize=4;
            print('3.4vsia.tif','-dtiff')
        end
        
        function obj=findind(obj)
            y=obj.uct_measure;
            unc=y.Efficiency(:,3)-y.Efficiency(:,1);
            ind=unc<0.1;
            india=abs([obj.AVGia2D]')<60;
            ind=ind & india;
            obj=obj(ind);
            [~,~,obj]=obj.bestangreg;
        end
        
        function y=surftemp(obj)
            datenumber=datenum(obj.Date);
            datestring=datestr(datenumber,'mmddyyyy');
            fname=['report_',datestring,'_Analysis.xlsx'];
            data=xlsread(fname,2);
            ms_data=datevec(data(:,1));
            ms_data=ms_data(:,4:5);
            ms_data=ms_data(:,1)*100+ms_data(:,2);
            intp=floor(obj.Time([1,end]));
            decp=round((obj.Time([1,end])-intp)*60);
            ms_obj=intp*100+decp;
            ind=ms_data>=ms_obj(1) & ms_data<=ms_obj(2);
            y=data(ind,25:27);
        end
        
        function [y,st]=surftable(obj)
            dt=vertcat(obj.Date);
            f=@(x) x.Time([1,end])';
            tm=arrayfun(f,obj,'uniformoutput',false);
            tm=cell2mat(tm(:)); 
            f=@(x) x.surftemp;
            st=arrayfun(f,obj,'uniformoutput',false);
            mst=cellfun(@mean,st(:),'uniformoutput',false);
            mst=cell2mat(mst);
%             dst=cellfun(@std,st(:),'uniformoutput',false);
%             dst=cell2mat(dst);
            f=@(x) range([x.TS2,x.TS3,x.Beam,x.Diffuse]);
            datarange=arrayfun(f,obj,'uniformoutput',false);
            datarange=cell2mat(datarange(:));
            uctscore=std(datarange./repmat([2,2,64,64],size(datarange,1),1),0,2);
            y=table(dt,tm,[obj.AVGTS2]'-273,mst(:,1),mst(:,2),mst(:,3),uctscore,...
                'variablenames',{'Date','Time','TS2','Bot_RB','Mid_LB','Top_MT','UncertaintyScore'});
        end
        
        function [y,st,obj]=surftable_cor(obj)
            [y,st]=obj.surftable;
            ind=y{:,4:6}<20;
            ind=~any(ind,2);
            y=y(ind,:);
            st=st(ind); st=st(:);
            obj=obj(ind);
            ilong=datenum(y.Date)<=datenum('24-Jul-2014');
            tarx=y{:,4:6};
            [x,i]=sort(tarx,2);
            [~,irev]=sort(i,2);
            m=size(x,1);
            x2=[ones(m,1),x(:,2)];
            b2=(x2(ilong,:)'*x2(ilong,:))\(x2(ilong,:)'*x(ilong,1));
            x3=[ones(m,1),x(:,3)];
            b3=(x3(ilong,:)'*x3(ilong,:))\(x3(ilong,:)'*x(ilong,1));
            xcor=[x(:,1),x2*b2,x3*b3];
            inds=sub2ind(size(x),repmat((1:m)',1,3),irev);
            tarxcor=xcor(inds);
            y{:,4:6}=tarxcor;
            n=datenum(y.Date)-datenum('Dec-31-2013');
            tm=mean(y.Time,2);
            [~,~,~,~,~,~,y.SolarTime]=IncAngle(n,tm,40.44,-86.91,-75,[],[],1);
            y.ia2D=[obj.AVGia2D]';
            y.Total=[obj.AVGBeam]'+[obj.AVGDiffuse]';
            y.Diffuse=[obj.AVGDiffuse]';
            y.Beam=[obj.AVGBeam]';
            y.fb=[obj.AVGfB]';
            y.OAT=[obj.AVGAmbient]'-273;
            ii=mat2cell(i,ones(size(i,1),1));
            stcor=cellfun(@(x,y) x(:,y),st,ii,'uniformoutput',false);
            dimrow=cellfun(@(x) size(x,1),stcor);
            stcor=cell2mat(stcor);
            m=size(stcor,1);
            stcor(:,2:3)=[ones(m,1),stcor(:,2),ones(m,1),stcor(:,3)]*blkdiag(b2,b3);
            stcor=mat2cell(stcor,dimrow);
            iir=mat2cell(irev,ones(size(irev,1),1));
            st=cellfun(@(x,y) x(:,y),stcor,iir,'uniformoutput',false);
        end
        
        function surffigure(obj)
            [~,st,obj]=obj.surftable_cor;
            stlong=st(15:23);
            stlong=cell2mat(stlong(:));
            objlong=obj(15:23);
            timelong=vertcat(objlong.Time)/24;
            n=datenum('Jul-22-0000');
            [~,~,~,~,~,~,solartimelong]=IncAngle(n,timelong*24,40.44,-86.91,-75,[],[],1);
            solartimelong=solartimelong/24;
            cmap=mycmap;
            % Figure 1
            figure; hold on
            h(1)=plot(timelong,vertcat(objlong.Beam)+vertcat(objlong.Diffuse),...
                'o','color',cmap{3}(3,:));
            h(2)=plot(timelong,vertcat(objlong.Beam),...
                '^','color',cmap{3}(2,:));
            ax1=gca;
            ax1_pos=ax1.Position;
            ax2=axes('Position',ax1_pos,...
                'YAxisLocation','right',...
                'XAxisLocation','top',...
                'color','none');
            ax2.NextPlot='add';
            h(3)=plot(ax2,solartimelong,stlong(:,1),...
                'd','color',cmap{2}(1,:));
            h(4)=plot(ax2,solartimelong,stlong(:,2),...
                '+','color',cmap{2}(2,:));
            h(5)=plot(ax2,solartimelong,stlong(:,3),...
                'x','color',cmap{2}(3,:));
            h(6)=plot(ax2,solartimelong,vertcat(objlong.Ambient)-273,...
                '*','color',cmap{1}(3,:));
            h(7)=plot(ax2,solartimelong,vertcat(objlong.TS2)-273,...
                '<','color',cmap{5}(3,:));
            ax1.YLim=[0,1300];
            ax2.YLim=[25,90];
            datetick(ax1,'x',16)
            [~,~,~,~,~,~,ast]=IncAngle(n,ax1.XTick*24,40.44,-86.91,-75,[],[],1);
            ast=ast/24;
            ax2.XLim=ast([1,end]);
            ax2.XTick=ast;
            ax2.XTickLabel=datestr(ast);
            ylabel(ax1,'Solar irradiance (W/m^2)')
            ylabel(ax2,'Temperature ({\circ}C)')
            xlabel(ax1,'Local time')
            xlabel(ax2,'Solar time')
            legend(h,'Global radiation','Beam radiation',...
                'TS1','TS2','TS3',...
                'Ambient temperature','Inlet temperature')
            FigureFormat(6,4,7)
            
            sttran=st(29:36);
            sttran=cell2mat(sttran(:));
            objtran=obj(29:36);
            timetran=vertcat(objtran.Time)/24;
            n=datenum('Jul-22-0000');
            [~,~,~,~,~,~,solartimetran]=IncAngle(n,timetran*24,40.44,-86.91,-75,[],[],1);
            solartimetran=solartimetran/24;
            % Figure 2
            clear h
            figure; hold on
            h(1)=plot(timetran,vertcat(objtran.Beam)+vertcat(objtran.Diffuse),...
                'o','color',cmap{3}(3,:));
            h(2)=plot(timetran,vertcat(objtran.Beam),...
                '^','color',cmap{3}(2,:));
            ax1=gca;
            ax1_pos=ax1.Position;
            ax2=axes('Position',ax1_pos,...
                'YAxisLocation','right',...
                'XAxisLocation','top',...
                'color','none');
            ax2.NextPlot='add';
            h(3)=plot(ax2,solartimetran,sttran(:,1),...
                'd','color',cmap{2}(1,:));
            h(4)=plot(ax2,solartimetran,sttran(:,2),...
                '+','color',cmap{2}(2,:));
            h(5)=plot(ax2,solartimetran,sttran(:,3),...
                'x','color',cmap{2}(3,:));
            h(6)=plot(ax2,solartimetran,vertcat(objtran.Ambient)-273,...
                '*','color',cmap{1}(3,:));
            h(7)=plot(ax2,solartimetran,vertcat(objtran.TS2)-273,...
                '<','color',cmap{5}(3,:));
            ax1.YLim=[0,1300];
            ax2.YLim=[20,90];
            datetick(ax1,'x',16)
            [~,~,~,~,~,~,ast]=IncAngle(n,ax1.XTick*24,40.44,-86.91,-75,[],[],1);
            ast=ast/24;
            ax2.XLim=ast([1,end]);
            ax2.XTick=ast;
            ax2.XTickLabel=datestr(ast);
            ylabel(ax1,'Solar irradiance (W/m^2)')
            ylabel(ax2,'Temperature ({\circ}C)')
            xlabel(ax1,'Local time')
            xlabel(ax2,'Solar time')
            legend(h,'Global radiation','Beam radiation',...
                'TS1','TS2','TS3',...
                'Ambient temperature','Inlet temperature',...
                'location','northwest')
            FigureFormat(6,4,7)
        end
    end
end