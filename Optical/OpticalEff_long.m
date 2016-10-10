function [q,varargout]=OpticalEff_long(data,type,surp,x,angle,varargin)
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

if isa(surp,'char')
    load(surp);
else
    R_g_e=surp{1};
    R_g_i=surp{2};
    R_refl=surp{3};
    R_abs=surp{4};
    T_g=surp{5};
end

assert(numel(angle) == 1);
assert(angle < 90 & angle >=0);

id = history.IncidentAngles == 0;
assert(sum(id) == 1);
X=history.x(:,:,id);
Y=history.y(:,:,id);
IA=history.ia(:,:,id);
ID=history.index(:,:,id);
[a,b,c]=size(X);

% preprocess IA
IA = rad2deg(acos(cos(deg2rad(IA)) * cos(deg2rad(angle))));

RatioLeave=ones(a,b,c);
% reflector
ID_refl=ID<=4 & ID>=1;
if x(2)
    RatioLeave(ID_refl)=interp1(0:90,R_refl,IA(ID_refl));
else
    RatioLeave(ID_refl)=R_refl(1);
end
% absorber
ID_abs=(ID==-1);
if x(1)
    RatioLeave(ID_abs)=interp1(0:90,R_abs,IA(ID_abs));
else
    RatioLeave(ID_abs)=R_abs(1);
end
% glass
ID_g_e=(ID==7 | ID==8);
ID_g_i=(ID==-7 | ID==-8);
if x(3)
    RatioLeave(ID_g_e)=(1-interp1(0:90,R_g_e,IA(ID_g_e))).*interp1(0:90,T_g,IA(ID_g_e));
    RatioLeave(ID_g_i)=1-interp1(0:90,R_g_i,IA(ID_g_i));
else
    RatioLeave(ID_g_e)=(1-R_g_e(1))*T_g(1);
    RatioLeave(ID_g_i)=1-R_g_i(1);
end

for j=3:b
    RatioLeave(:,j,:)=RatioLeave(:,j-1,:).*RatioLeave(:,j,:);
end
Incoming=RatioLeave;
Incoming(:,2:b,:)=RatioLeave(:,1:b-1,:);
SolarFlux=Incoming-RatioLeave;

if strcmpi(type,'local')
    if isempty(varargin)
        thm=NumericThermal(history.geometry);
        mesh=thm.meshing;
    else
        mesh=varargin{1};
    end
    t=mesh.connectivity;
    pt=mesh.points;
    tag=mesh.tag.tag;
    varargout{1}=mesh;
    
    q=zeros(size(tag,1),c);
    % Absorber
    [i,j]=find(ID_abs);
    x_abs=X((j-1)*a+i); y_abs=Y((j-1)*a+i); sf_abs=SolarFlux((j-1)*a+i);
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
    q_abs=sparse(rs,col_abs(cs),sf_abs(cs),size(t_abs,1),c);
    
    % Reflector
    [i,j]=find(ID_refl);
    col_refl=ceil(j/b);
    sf_refl=SolarFlux((j-1)*a+i);
    q_refl=sparse(ones(size(j)),col_refl,sf_refl,1,c);
    
    % Glass
    [i,j]=find(ID_g_e);
    x_glass=X((j-1)*a+i); y_glass=Y((j-1)*a+i);
    col_glass=ceil(j/b);
    Re=interp1(0:90,R_g_e,IA((j-1)*a+i));
    T=interp1(0:90,T_g,IA((j-1)*a+i));
    Ri=interp1(0:90,R_g_i,IA((j-1)*a+i));
    sf_glass=Incoming((j-1)*a+i).*(1-Re).*(1-T).*(1+T.*Ri)./(1-Re.*Ri.*T.^2);
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
    q_glass=sparse(rs,col_glass(cs),sf_glass(cs),size(t_glass,1),c);
    
    % All
    q(tag==2 | tag==-3 | tag==-4,:)=q_abs;
    q(tag==1,:)=q_glass;
    q=[q;q_refl]/a;
%     r_a=history.geometry.RadReceiver;
%     w=history.geometry.Width;
%     [yl1,yl2]=history.geometry.profile; yl=[yl1;yl2];
%     h=max(yl(:,2))-min(yl(:,2));
%     w_abs=(r_a-h*sin((abs(angle))/180*pi))./cos((abs(angle))/180*pi)+w/2;
%     w_aper=max(w,w_abs);
%     q=q.*repmat(w_aper',size(q,1),1); %m
else
    q=zeros(3,c);
    % Absorber
    [i,j]=find(ID_abs);
    q_abs=sparse(i,j,SolarFlux((j-1)*a+i),a,b*c);
    q_abs=sum(q_abs); q_abs=full(reshape(q_abs,b,c));
    q(1,:)=sum(q_abs);
    % Reflector
    [i,j]=find(ID_refl);
    q_refl=sparse(i,j,SolarFlux((j-1)*a+i),a,b*c);
    q_refl=sum(q_refl); q_refl=full(reshape(q_refl,b,c));
    q(2,:)=sum(q_refl);
    % Glass
    [i,j]=find(ID_g_e);
    T=interp1(0:90,T_g,IA((j-1)*a+i));
    q_glass=sparse(i,j,RatioLeave((j-1)*a+i).*(1./T-1),a,b*c);
    q_glass=sum(q_glass); q_glass=full(reshape(q_glass,b,c));
    q(3,:)=sum(q_glass);
    q=q/a;
    varargout{1}={'You'' calling angular/total whatever efficiency. It''s none of my business! Let go of me!'};
end
end