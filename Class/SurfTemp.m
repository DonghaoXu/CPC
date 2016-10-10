classdef SurfTemp < ExpData
    properties
        
    end
    
    methods (Static)
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
            Data_raw.Surf1=raw(:,25);
            Data_raw.Surf2=raw(:,26);
            Data_raw.Surf3=raw(:,27);
        end
        
        function crt=samplecriteria()
            crt.GlobalMin=0;
            crt.GlobalMax=inf;
            crt.BeamMin=0;
            crt.BeamMax=inf;
            crt.DiffuseMin=0;
            crt.DiffuseMax=inf;
            crt.BeamFractionMin=0;
            crt.BeamFractionMax=inf;
            crt.WspdMin=0;
            crt.WspdMax=4.5; %m/s
            crt.NoData=10;
            crt.VarGlobal=64;
            crt.VarTemp=2;
            crt.VarFR=0.4;
        end
        
        function y=filterdata(Data_raw,crt)
            % prepare the data
            Data_raw.Global = Data_raw.Total;
            Data_raw.Beam = Data_raw.Total - Data_raw.Diffuse;
            Data_raw.BeamFraction = Data_raw.Beam ./ Data_raw.Total;
            
            % Boundary criteria
            ind = (Data_raw.Global > crt.GlobalMin) & ...
                (Data_raw.Global < crt.GlobalMax) & ...
                (Data_raw.Beam > crt.BeamMin) & ...
                (Data_raw.Beam < crt.BeamMax) & ...
                (Data_raw.Diffuse > crt.DiffuseMin) & ...
                (Data_raw.Diffuse < crt.DiffuseMax) & ...
                (Data_raw.BeamFraction > crt.BeamFractionMin) & ...
                (Data_raw.BeamFraction < crt.BeamFractionMax) & ...
                (Data_raw.Wspd > crt.WspdMin) & ...
                (Data_raw.Wspd < crt.WspdMax);
            id = find(ind);

            % Quasi-steady state criteria
            Data_tbl = struct2table(Data_raw);
            Data = Data_tbl(ind, {'Beam','Diffuse','TS2','TS3','FR'});
            Data_dif = diff(table2array(Data));
            
            
            % # of data criteria
            did = diff(id);
            find(did ~= 1);
            
            
            %filter by wind speed and solar irradiance
            i_first=find(Data_raw(:,4)>crt.globalmin,1,'first');
            i_last=find(Data_raw(:,4)>Solar_min,1,'last');
            Data=Data_raw(i_first:i_last,:);
            m=size(Data,1);
            
            %difference of two consecutive data points
            solardiff=[Data(2:end,4)-Data(1:end-1,4);0];
            T2diff=[Data(2:end,9)-Data(1:end-1,9);0];
            T3diff=[Data(2:end,10)-Data(1:end-1,10);0];
            FRdiff=[Data(2:end,6)-Data(1:end-1,6);0];
            
            %locate the points that do not satisfy
            ind_range=(Data(:,5)>wind_max | Data(:,4)<Solar_min | Data(:,6)<1)*2;
            ind_diff=(abs(solardiff)>64 | abs(T2diff)>2 | abs(T3diff)>2 | abs(FRdiff)>0.4);
            ind_all=ind_range+ind_diff;
            ind_data=find(ind_all~=0);
            ind_data=[0;ind_data;m];
            bool_range=(ind_all(ind_all~=0)>1);
            bool_range=[bool_range;0];
            %find how long it lasts between two unsatisfying points
            Duration=ind_data(2:end)-ind_data(1:end-1)-bool_range;
            %the nth period lasts more than 10 min, stored in ind_dura
            ind_dura=find(Duration>9);
            %there are m periods
            m=size(ind_dura,1);
            if m==0
                fy=[];
                return
            end
            dataset_raw=repmat(struct('data',[]),m,1);
            for i=1:m
                a=ind_data(ind_dura(i))+ind_range(ind_dura(i))/2+1;
                b=ind_data(ind_dura(i)+1)-ind_range(ind_dura(i)+1)/2;
                dataset_raw(i,1).data=Data(a:b,:);
            end
            
            %check the range of solar irradiance and temperatures
            numperiod=zeros(m,1);
            ftime=Data_raw(1,1);
            parfor i=1:m
                solartemp=dataset_raw(i,1).data;
                dataset_fil=datafilter(solartemp,ftime);
                dataset(i,1).row=dataset_fil;
                numperiod(i,1)=size(dataset_fil,1);
            end
            numperiod=[0;numperiod];
            y=zeros(sum(numperiod),2);
            for i=1:m
                y(sum(numperiod(1:i,1))+1:sum(numperiod(1:i+1,1)),:)=dataset(i,1).row;
            end
            y=y(y(:,1)~=0,:);
            m=size(y,1);
            if m~=0
                parfor i=1:m
                    fy(i,1)=dataproc(filename,Data_raw,Property,y(i,:),x);
                end
            else
                fy=[];
            end
            
        end
        
        
        function y=datafilter(data,time)
            i=0;
            m=size(data,1);
            y=zeros(floor(m/10),2);
            datatemp=data;
            while size(datatemp,1)>9
                i=i+1;
                [y(i,:),datatemp2]=datafilter_single(datatemp,time);
                datatemp=datatemp2;
            end
        end
    end
end