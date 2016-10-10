classdef tmy3
    properties
        code;station;state;
        latitude;longitude;timezone;
        doy;hod;
        GHI;DNI;dHI;
        DryBulb;DewPoint;RH;Wspd;
    end
    
    methods
        function r=tmy3(val)
            % read .csv file
            if nargin>0
                if strcmpi(val(end-2:end),'csv') 
                    if strcmpi(val(1:7),'http://')
                        options=weboptions('Timeout',30);
                        filename=websave([tempname,'.csv'],val,options);
                    else
                        filename=val;
                    end
                    fid = fopen(filename);
                    out = textscan(fid,'%s','delimiter',',');
                    fclose(fid);
                    out=out{1,1};
                    loc=out(1:7,1);
                    r.latitude=str2double(loc{5,1});
                    r.longitude=str2double(loc{6,1});
                    r.timezone=str2double(loc{4,1})*15;
                    r.station=loc{2,1}(2:end-1);
                    r.state=loc{3,1};
                    r.code=loc{1,1};
                    data=reshape(out(8:end,1),71,8761)';
                    r.doy=reshape(repmat(1:365,24,1),[],1);
                    r.hod=reshape(repmat((1:24)',1,365),[],1);
                    r.GHI=cellfun(@str2num,data(2:end,5));
                    r.DNI=cellfun(@str2num,data(2:end,8));
                    r.dHI=cellfun(@str2num,data(2:end,11));
                    r.DryBulb=cellfun(@str2num,data(2:end,32))+273;
                    r.DewPoint=cellfun(@str2num,data(2:end,35));
                    r.RH=cellfun(@str2num,data(2:end,38))/100;
                    r.Wspd=cellfun(@str2num,data(2:end,47));
                    if strcmpi(val(1:7),'http://')
                        delete(filename)
                    end
                elseif isa(val,'NumericThermal') || isa(val,'SimpleThermal')
                    r=val.outdoor;
                end
            end
        end
        
        function [ia3D,ia2D,alpha,phi,delta,ha]=angles(obj,beta,varargin)
            % Input:
            % beta: tilted angle (DEG)
            % gama: surface azimuth (DEG)
            % Output:
            % ia3D: 3D incident angle (DEG)
            % ia2D: 2D incident angle projected to the cross-section of
            % certain trough (DEG)
            % alpha: solar altitude (DEG)
            % phi: solar azimuth (DEG)
            % delta: solar declination angle (DEG)
            % ha: hour angle (DEG)
            % varargin: specify the date, eg. '3/23'.
            in=[varargin(:);cell(1,1)]; in=in(1);
            if isempty(in{1})
                ind=1:numel(obj.hod);
            elseif isa(in{1},'numeric') || isa(in{1},'logical')
                ind=in{1};
            elseif isa(in{1},'char')
                ind=datenum(['0/',in{1}]);
            else
                in{1}=cellfun(@(x) strcat('0/',x),in{1},'uniformoutput',false);
                ind=datenum(in{1});
            end
            [ia3D,~,alpha,phi,delta,ha]=IncAngle(obj.doy(ind),obj.hod(ind),...
                obj.latitude,obj.longitude,obj.timezone,...
                beta,0);
            ia2D=IncAng2D(alpha,phi,beta);
        end
        
        function [beam,diffuse_sky,diffuse_grd]=tilt(obj,rho_grd,beta,varargin)
            % rho_grd: ground reflectance
            % beta: tilted angle
            % varargin: specify the date, eg. '3/23'.
            in=[varargin(:);cell(1,1)]; in=in(1);
            if isempty(in{1})
                ind=1:numel(obj.GHI);
            elseif isa(in{1},'numeric') || isa(in{1},'logical')
                ind=in{1};
            elseif isa(in{1},'char')
                ind=datenum(['0/',in{1}]);
            else
                in{1}=cellfun(@(x) strcat('0/',x),in{1},'uniformoutput',false);
                ind=datenum(in{1});
            end
            ia3D=obj.angles(beta,ind)*pi/180;
            beam=max(0,obj.DNI(ind).*cos(ia3D));
            diffuse_sky=obj.dHI(ind)*(1+cos(beta*pi/180))/2;
            diffuse_grd=obj.GHI(ind)*rho_grd*(1-cos(beta*pi/180))/2;
        end
        
        function r=getweather(obj,i)
            r=repmat(obj,1,numel(i));
            for j=1:numel(i)
                r(j).doy=obj.doy(i(j));
                r(j).hod=obj.hod(i(j));
                r(j).GHI=obj.GHI(i(j));
                r(j).DNI=obj.DNI(i(j));
                r(j).dHI=obj.dHI(i(j));
                r(j).DryBulb=obj.DryBulb(i(j));
                r(j).Wspd=obj.Wspd(i(j));
            end
        end
        
        function r=dummy(obj,n)
            r=obj;
            r.doy=sparse(n,1);
            r.hod=sparse(n,1);
            r.GHI=sparse(n,1);
            r.DNI=sparse(n,1);
            r.dHI=sparse(n,1);
            r.DryBulb=sparse(n,1);
            r.Wspd=sparse(n,1);
        end
        
        function [b,ax1,l12,ax2,l21,l22]=monthly(obj,ind)
            cmap=mycmap();
            d=ceil((1:8760)'/24);
            m=month(d(ind));
            bm=accumarray(m,obj.GHI(ind)-obj.dHI(ind));
            dfs=accumarray(m,obj.dHI(ind));
            ta=accumarray(m,obj.DryBulb(ind),[],@mean)-273;
            vel=accumarray(m,obj.Wspd(ind),[],@mean);
            [ax1,l11,l12]=plotyy(1:12,bm+dfs,1:12,bm./(bm+dfs),'bar','plot');
            delete(l11)
            hold on
            b=bar(ax1(1),[bm,dfs],'stacked');
            b(1).FaceColor=cmap{3}(3,:);
            b(2).FaceColor=cmap{2}(3,:);
            l12.Marker='o';
            l12.MarkerFaceColor=[1 1 1];
            l12.MarkerEdgeColor=cmap{1}(3,:);
            l12.Color=cmap{1}(3,:);
            ax1(2).YColor='k';
            ax1(1).XLim=[0,13];
            ax1(2).XLim=[0,13];
            y2=max(bm+dfs);
            y2=num2str(y2);
            ax1(1).YLim(2)=10^(numel(y2)-1)*(str2double(y2(1))+1+0.5*(str2double(y2(2))>=0.5));
            ax1(2).YLim=[0,1];
            ax1(1).XTickLabel={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
            ax1(1).YTick=linspace(ax1(1).YLim(1),ax1(1).YLim(2),6);
            ax1(2).YTick=linspace(0,1,6);
%             ax(2).YTickLabel={'0',}
            ylabel(ax1(1),'Monthly total radiation (Wh/m^2)')
            ylabel(ax1(2),'Monthly fraction of beam radiation')
            FigureFormat(5,10/3,7)
            l12.MarkerSize=5;
            l12.LineWidth=0.5;
            legend('Horizontal beam radiation','Horizontal diffuse radiation','Fraction of beam radiation',...
                'location','northwest')
            
            figure;
            [ax2,l21,l22]=plotyy(1:12,ta,1:12,vel);
            ax2(1).XTick=1:12;
            ax2(1).XLim=[0,13];
            ax2(2).XLim=[0,13];
            ax2(1).XTickLabel={'Jan','Feb','Mar','Apr','May','Jun','Jul','Aug','Sep','Oct','Nov','Dec'};
            ax2(1).YColor='k';
            ax2(2).YColor='k';
            ax2(2).YLim=[0,10];
            ax2(1).YTick=linspace(ax2(1).YLim(1),ax2(1).YLim(2),6);
            ax2(2).YTick=linspace(ax2(2).YLim(1),ax2(2).YLim(2),6);
            l21.Marker='o';
            l21.Color=cmap{5}(2,:);
            l22.Marker='^';
            l22.Color=cmap{1}(3,:);
            ylabel(ax2(1),'Monthly Average ambient temperature ({\circ}C)')
            ylabel(ax2(2),'Monthly Average wind speed (m/s)')
            legend('Ambient temperature','Wind speed')
            FigureFormat(5,10/3,7)
            l21.MarkerSize=5;
            l21.MarkerFaceColor=[1 1 1];
            l22.MarkerSize=5;
            l22.MarkerFaceColor=[1 1 1];
        end
    end
    
    methods (Static)
        function fd=downloadtmy3(city,varargin)
            url='http://rredc.nrel.gov/solar/old_data/nsrdb/1991-2005/data/tmy3/';
            if ismac
                sep='';
            elseif ispc
                sep=' ';
            end
            if isempty(varargin)
                fpath=['\Users\Donghao',sep,'Xu\Dropbox\DonghaoCode\Data\Weather\'];
            else
                fpath=varargin{1};
            end
            usafn=tmy3.findcity(city);
            link=strcat(url,usafn{1},'TYA.csv');
            options=weboptions('Timeout',30);
            fd=websave([fpath,usafn{1},'TYA.csv'],link,options);
        end
        
        function [usafn,fname]=findcity(city)
            d=load('USAF.mat');
            USAF=d.USAF;
            citylist=USAF(:,2);
            aa=strfind(citylist,city);
            id=~cellfun(@isempty,aa);
            usafn=USAF(id,:);
            fname=usafn(:,1);
            apdx=repmat({'TYA.csv'},numel(fname),1);
            fname=cellfun(@horzcat,fname,apdx,'uniformoutput',false);
        end
    end
end