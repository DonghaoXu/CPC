function [p,t,tag,n,AB]=CPCmesh1d(varargin)
%   r_a:     Radius of absorber
%   r_g:     Radius of glass
%   gap:    Extra gap between reflector and glass
%   acc:    Acceptance half angle
%   cr:     Geometric concentration ratio
%   N:      Number of CPCs in series
if nargin==0
    error('Get me some inputs, dude! I need a geometry class input, or r_a, r_g, r_u, gap, theta, and cr!')
end
if ~isa(varargin{1},'CPC')
    error('Oh, please, dude. Read about CPC class and get me one!')
end

if numel(varargin)<=2
    r_a=varargin{1}.RadReceiver;
    r_g=varargin{1}.RadGlass;
    r_u=varargin{1}.RadUtube;
    gap=varargin{1}.Gap;
    theta=varargin{1}.AccHalfAng;
    cr=varargin{1}.CR;
    elesize=[];
    if numel(varargin)==2 && isa(varargin{2},'numeric')
        elesize=varargin{2};
    end
elseif numel(varargin)>=6
    r_a=varargin{1};
    r_g=varargin{2};
    r_u=varargin{3};
    gap=varargin{4};
    theta=varargin{5};
    cr=varargin{6};
    elesize=[];
    if numel(varargin)>=7 && isa(varargin{7},'numeric')
        elesize=varargin{7};
    end
end
tol=r_a*1e-6;
h0=min([r_a*pi/60,r_u/2,elesize]);
% d=h0*5;
% t_g=0;

%0. Auto-check if theta is in Rad or Degree
if theta>pi/2     %In Degree
    theta=theta*pi/180;
end

%1. Check if geometric concentration ratio is out of bound
w_max=(2*sqrt((r_g+gap)^2-r_a^2)+(pi+2*asin(r_a/(r_g+gap)))*r_a)/sin(theta);
w=2*pi*r_a*cr;
if w-w_max>tol
    error('Geometric concentration ratio is out of bound. The maximum concentration ratio is %.3f.',w_max/2/pi/r_a)
else
    w=min(w,w_max);
end

%2. CPC reflector
% [~,~,spl]=sprefl(r_a,r_g,gap,theta);
% [yl1,yl2,~,~,~,~,~,~]=profileAnaXCPC(theta,r_a,w,r_g,t_g,gap,spl);
% yl=[yl1;yl2];
% [~,I]=sort(yl(:,1));
% yl=yl(I,:);
% yl(end,1)=0;
% yl1_sampled=sample(yl(end:-1:1,:),d);
% yall1_sampled=[yl1_sampled;-yl1_sampled(:,1),yl1_sampled(:,2)];
% p_cpc=unique(yall1_sampled,'rows');
% Ncpc=size(p_cpc,1);
% t_cpc=[1:Ncpc-1;2:Ncpc]';
% tag_cpc=zeros(size(t_cpc,1),1);
p_cpc=[];
Ncpc=0;
t_cpc=[];
tag_cpc=[];

%3. Glass envelope
step_glass=h0/r_g;
step_glass=pi/ceil(pi/step_glass);
cgx=r_g*(cos(step_glass/2:step_glass:2*pi))';
cgy=r_g*(sin(step_glass/2:step_glass:2*pi))';
p_glass=[cgx,cgy];
Ng=size(p_glass,1);
t_glass=[(1:Ng)',[2:Ng,1]'];
tag_glass=ones(size(t_glass,1),1);
ab_g=p_glass(t_glass(:,2),:)-p_glass(t_glass(:,1),:);
n_g=[ab_g(:,2),-ab_g(:,1)];
AB_g=sqrt(sum(n_g.^2,2));
n_g=n_g./[AB_g,AB_g];

%4. Absorber
step_abs=h0/r_a;
step_abs=pi/ceil(pi/step_abs);
cax=r_a*(cos(step_abs/2:step_abs:2*pi))';
cay=r_a*(sin(step_abs/2:step_abs:2*pi))';
p1x=cax(cay<r_u & cax<0 & cay>-r_u);
p1x1=p1x(1);
p1x2=p1x(end);
p1y=cay(cay<r_u & cax<0 & cay>-r_u);
p1y1=p1y(1);
p1y2=p1y(end);
p2x1=cax(cay<r_u & cax>0 & cay>0);
p2x1=p2x1(end);
p2y1=cay(cay<r_u & cax>0 & cay>0);
p2y1=p2y1(end);
p2x2=cax(cay<0 & cax>0 & cay>-r_u);
p2x2=p2x2(1);
p2y2=cay(cay<0 & cax>0 & cay>-r_u);
p2y2=p2y2(1);
Nabs=size(cax,1);
p_abs=[cax,cay];
t_abs=[(1:Nabs)',[2:Nabs,1]'];
tag_abs=2*ones(size(t_abs,1),1);
iut=max(abs([p_abs(t_abs(:,1),2),p_abs(t_abs(:,2),2)]),[],2)<r_u;
tag_abs(iut & p_abs(t_abs(:,1),1)>0)=-4;
tag_abs(iut & p_abs(t_abs(:,1),1)<0)=-3;
ab_abs=p_abs(t_abs(:,2),:)-p_abs(t_abs(:,1),:);
n_abs=[ab_abs(:,2),-ab_abs(:,1)];
AB_abs=sqrt(sum(n_abs.^2,2));
n_abs=n_abs./[AB_abs,AB_abs];
% nn=floor(asin(r_u/r_a)/step_abs-0.5);
% i_out=t_abs(:,1);
% i_out=i_out([1:nn,Nabs-nn:Nabs]);
% i_out=i_out(:);
% tag_abs(max(abs([p_abs(t_abs(:,1),2),p_abs(t_abs(:,2),2)]),[],2)<r_u)=3;
% tag_abs(i_out)=4;

%5. u-tube
step_utube=h0/r_u;
step_utube=pi/ceil(pi/step_utube);
cux=r_u*(cos(pi/2:-step_utube:-pi/2))'-(r_a-r_u);
cuy=r_u*(sin(pi/2:-step_utube:-pi/2))';
cux1=[p1x1;cux;p1x2];
cux2=[p2x1;-cux;p2x2];
cuy1=[p1y1;cuy;p1y2];
cuy2=[p2y1;cuy;p2y2];
p_ut1=[cux1,cuy1];
Nut=size(p_ut1,1);
t_ut11=(1:Nut-1)';
t_ut12=(2:Nut)';
t_ut1=[t_ut11(:),t_ut12(:)];
ab_ut1=p_ut1(t_ut1(:,1),:)-p_ut1(t_ut1(:,2),:);
n_ut1=[ab_ut1(:,2),-ab_ut1(:,1)];
AB_ut1=sqrt(sum(n_ut1.^2,2));
n_ut1=n_ut1./[AB_ut1,AB_ut1];
p_ut2=[cux2,cuy2];
t_ut2=t_ut1;
ab_ut2=p_ut2(t_ut2(:,2),:)-p_ut2(t_ut2(:,1),:);
n_ut2=[ab_ut2(:,2),-ab_ut2(:,1)];
AB_ut2=sqrt(sum(n_ut2.^2,2));
n_ut2=n_ut2./[AB_ut2,AB_ut2];
tag_ut1=3*ones(size(t_ut1,1),1);
tag_ut2=4*ones(size(t_ut2,1),1);

%6. Assemble a CPC
p=[p_cpc;p_glass;p_abs;p_ut1;p_ut2];
t=[t_cpc;t_glass+Ncpc;t_abs+Ncpc+Ng;t_ut1+Ncpc+Ng+Nabs;t_ut2+Ncpc+Ng+Nabs+Nut];
n=[n_g;n_abs;n_ut1;n_ut2];
AB=[AB_g;AB_abs;AB_ut1;AB_ut2];
% 1 glass
% 2 absorber
% 3 u-tube inlet  -3 u-tube inlet shared with absorber
% 4 u-tube outlet  -4 u-tube outlet shared with absorber
% 5 reflector
tag=sparse([tag_cpc;tag_glass;tag_abs;tag_ut1;tag_ut2]);

%7. Extra topology


%8. Replicate
N=1;
Np=size(p,1);
px=repmat(p(:,1),1,N)+repmat(w*((1:N)-1)-w*(N-1)/2,Np,1);
py=repmat(p(:,2),1,N);
p=[px(:),py(:)];
Nt=size(t,1);
t1=repmat(t(:,1),1,N)+repmat(Np*(0:N-1),Nt,1);
t2=repmat(t(:,2),1,N)+repmat(Np*(0:N-1),Nt,1);
t=[t1(:),t2(:)];
tag=repmat(tag,N,1);
[t,ia]=unique(t,'rows');
tag=tag(ia);
[p,~,ic]=unique(p,'rows');
t=ic(t(:));
t=reshape(t,[],2);

%9. Counterclockwise regulation
p1=p(t(:,1),:); p2=p(t(:,2),:);
n_12=p2-p1;
c=cross([n_12,zeros(size(n_12,1),1)],[n,zeros(size(n_12,1),1)]);
t(c(:,3)<0,:)=[t(c(:,3)<0,2),t(c(:,3)<0,1)];
end