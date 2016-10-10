function [intf,varargout]=InterFactor(data,type,angle,varargin)
% local optical efficiency
% varargin{1}: mesh struct (optional)
% x: ctrl for absorber, reflector, glass. 0 for ignore angular dependence,
% other for consider
if isa(data,'char')
    data=load(data,'history');
    history=RayTraceData(data.history);
else
    history=RayTraceData(data);
end

[~,id]=ismember(abs(angle),history.IncidentAngles);
if any(~id)
    f=@(x) fprintf('I don''t have ray tracing data for angle: %.1f',x);
    cellfun(f,mat2cell(angle(~id),sum(~id),1))
    error('Bang!')
end
X=history.x(:,:,id);
X(:,:,angle<0)=-X(:,:,angle<0);
Y=history.y(:,:,id);
ID=history.index(:,:,id);
[a,b,c]=size(X);

% reflector
ID_refl=ID<=4 & ID>=1;
% absorber
ID_abs=(ID==-1);
% glass
ID_g_e=(ID==7 | ID==8);
% ID_g_i=(ID==-7 | ID==-8);

if strcmpi(type,'local')
    if isempty(varargin)
        thm=NumericThermal(history.geometry);
        mesh=thm.meshing;
%         mesh=thm.mesh;
    else
        mesh=varargin{1};
    end
    t=mesh.connectivity;
    pt=mesh.points;
    tag=mesh.tag.tag;
    varargout{1}=mesh;
    
    intf=zeros(size(tag,1),c);
    % Absorber
    [i,j]=find(ID_abs);
    x_abs=X((j-1)*a+i); y_abs=Y((j-1)*a+i);
    col_abs=ceil(j/b);
    t_abs=t(tag==2 | tag==-3 | tag==-4,:);
    bmin=atan2(pt(t_abs(:,1),2),pt(t_abs(:,1),1));
    bmax=atan2(pt(t_abs(:,2),2),pt(t_abs(:,2),1));
    bin=[bmin,bmax]; bin=sort(bin,2); Nb=size(bin,1);
    [~,iend]=max(bin(:,2)-bin(:,1));
    loc_abs=atan2(y_abs,x_abs); Nloc=size(loc_abs,1);
    bmin=bin(:,1); bmax=bin(:,2)-1e-8;
    imin=repmat(loc_abs',Nb,1)>repmat(bmin,1,Nloc);
    imax=repmat(loc_abs',Nb,1)<repmat(bmax,1,Nloc);
    ibin=imin & imax;
    ibin(iend,:)=~any(ibin);
    [rs,cs]=find(ibin);
    intf_abs=sparse(rs,col_abs(cs),1,size(t_abs,1),c);
    
    % Reflector
    [~,j]=find(ID_refl);
    col_refl=ceil(j/b);
    intf_refl=sparse(ones(size(j)),col_refl,1,1,c);
    
    % Glass
    [i,j]=find(ID_g_e);
    x_glass=X((j-1)*a+i); y_glass=Y((j-1)*a+i);
    col_glass=ceil(j/b);
    t_glass=t(tag==1,:);
    bmin=atan2(pt(t_glass(:,1),2),pt(t_glass(:,1),1));
    bmax=atan2(pt(t_glass(:,2),2),pt(t_glass(:,2),1));
    bin=[bmin,bmax]; bin=sort(bin,2); Nb=size(bin,1);
    [~,iend]=max(bin(:,2)-bin(:,1));
    loc_glass=atan2(y_glass,x_glass); Nloc=size(loc_glass,1);
    bmin=bin(:,1); bmax=bin(:,2)-1e-8;
    imin=repmat(loc_glass',Nb,1)>repmat(bmin,1,Nloc);
    imax=repmat(loc_glass',Nb,1)<repmat(bmax,1,Nloc);
    ibin=imin & imax;
    ibin(iend,:)=~any(ibin);
    [rs,cs]=find(ibin);
    intf_glass=sparse(rs,col_glass(cs),1,size(t_glass,1),c);
     
    % All
    intf(tag==2 | tag==-3 | tag==-4,:)=intf_abs;
    intf(tag==1,:)=intf_glass;
    intf=[intf;intf_refl]/a;
else
    intf=zeros(3,c);
    % Absorber
    [i,j]=find(ID_abs);
    intf_abs=sparse(i,j,1,a,b*c);
    intf_abs=sum(intf_abs); intf_abs=full(reshape(intf_abs,b,c));
    intf(1,:)=sum(intf_abs);
    % Reflector
    [i,j]=find(ID_refl);
    intf_refl=sparse(i,j,1,a,b*c);
    intf_refl=sum(intf_refl); intf_refl=full(reshape(intf_refl,b,c));
    intf(2,:)=sum(intf_refl);
    % Glass
    [i,j]=find(ID==7);
    intf_glass=sparse(i,j,1,a,b*c);
    intf_glass=sum(intf_glass); intf_glass=full(reshape(intf_glass,b,c));
    intf(3,:)=sum(intf_glass);
    intf=intf/a;
    varargout{1}={'You'' calling angular/total whatever efficiency. It''s none of my business! Let go of me!'};
end
end