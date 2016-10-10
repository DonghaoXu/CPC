function VF=viewfactor(varargin)
% Input was:
% mesh struct: tag,normvec,pt,t,pmid,lface,lcpc
% CPCgeo class: top_ref,r_a,r_g,r_u
% VF is the view factor matrix
% VF is (N+1)x(N+2) if t is Nx2
% VF(:,N+1) represents the view factor to reflector
% VF(:,N+2) represents the view factor to sky
% VF(i,j) is the view factor from i to j
% int(R^2/(R^2+(x-y)^2)^2,x,y)=-atan((x-y)/R)*(x-y)/(2R)
% Between any two surfaces i and j, assuming to be infinitely long
% VF(i,j)=abs(R11+R22-R12-R21)/2/L_i
%         L_j
%      --------
%      |\    /|
%      | \ R12|
%   R11|  \/  |R22
%      |  /\  |
%      | / R21|
%      |/    \|
%      --------
%         L_i
% global tag normvec pt t pmid lface lcpc

if nargin==0
    error('...I''m silenced by the voidness in the input fields...')
end

if nargin==2 && isa(varargin{1},'struct') && isa(varargin{2},'CPC')
    tag=varargin{1}.tag.tag;
    normvec=varargin{1}.normvec;
    pt=varargin{1}.points;
    t=varargin{1}.connectivity;
    pmid=varargin{1}.pmid;
    lface=varargin{1}.lface;
    lcpc=varargin{1}.lcpc;
    [yl1,yl2]=varargin{2}.profile; yl=[yl1;yl2]; 
    yl=yl(~isnan(yl(:,1)),:); [~,I]=sort(yl(:,1)); yl=yl(I,:);
    top_ref=yl(1,:);
    r_a=varargin{2}.RadReceiver;
    r_g=varargin{2}.RadGlass;
    r_u=varargin{2}.RadUtube;
elseif nargin==8 && isa(varargin{8},'CPC')
    tag=varargin{1};
    normvec=varargin{2};
    pt=varargin{3};
    t=varargin{4};
    pmid=varargin{5};
    lface=varargin{6};
    lcpc=varargin{7};
    [yl1,yl2]=varargin{2}.profile; yl=[yl1;yl2]; 
    yl=yl(~isnan(yl(:,1)),:); [~,I]=sort(yl(:,1)); yl=yl(I,:);
    top_ref=yl(1,:);
    r_a=varargin{8}.RadReceiver;
    r_g=varargin{8}.RadGlass;
    r_u=varargin{8}.RadUtube;
elseif nargin==5 && isa(varargin{1},'struct')
    tag=varargin{1}.tag.tag;
    normvec=varargin{1}.normvec;
    pt=varargin{1}.points;
    t=varargin{1}.connectivity;
    pmid=varargin{1}.pmid;
    lface=varargin{1}.lface;
    lcpc=varargin{1}.lcpc;
    top_ref=varargin{2};
    r_a=varargin{3};
    r_g=varargin{4};
    r_u=varargin{5}; 
elseif nargin==11
    tag=varargin{1};
    normvec=varargin{2};
    pt=varargin{3};
    t=varargin{4};
    pmid=varargin{5};
    lface=varargin{6};
    lcpc=varargin{7};
    top_ref=varargin{8};
    r_a=varargin{9};
    r_g=varargin{10};
    r_u=varargin{11};    
else
    error('The smartest way for input is a mesh struct and a CPC class. Otherwise, balabala... Why not come in and have a drink?')
end



N=size(t,1);

% 2. ID (just consider i to j where j>i) row->col
id=repmat((1:N)',1,N); ihalf=(id<repmat(1:N,N,1));
[row,col]=find(ihalf);
VF=zeros(size(row));

%1. glass to sky
igr=tag==1; row_g2sky=id(igr,1); col_g2sky=(N+2)*ones(size(row_g2sky));
n_P=normvec(igr,:);
n_PA=repmat(top_ref,size(row_g2sky,1),1)-pmid(igr,:);
PA=sqrt(sum(n_PA.^2,2)); n_PA=n_PA./[PA,PA];
n_PB=repmat([-top_ref(1),top_ref(2)],size(row_g2sky,1),1)-pmid(igr,:);
PB=sqrt(sum(n_PB.^2,2)); n_PB=n_PB./[PB,PB];
iA=sum(n_PA.*n_P,2); iB=sum(n_PB.*n_P,2);
R11=sqrt(sum((pt(t(igr,1),:)-repmat(top_ref,size(row_g2sky,1),1)).^2,2));
R22=sqrt(sum((pt(t(igr,2),:)-repmat([-top_ref(1),top_ref(2)],size(row_g2sky,1),1)).^2,2));
R21=sqrt(sum((pt(t(igr,2),:)-repmat(top_ref,size(row_g2sky,1),1)).^2,2));
R12=sqrt(sum((pt(t(igr,1),:)-repmat([-top_ref(1),top_ref(2)],size(row_g2sky,1),1)).^2,2));
VF_g2sky=(R12+R21-R11-R22)/2./lface(row_g2sky); % default: glass sees A and B
%1.1 glass doesn't see A nor B and faces down
VF_g2sky(iA<0 & iB<0 & n_P(:,2)<0)=0;
%1.2 glass doesn't see A nor B and faces up
VF_g2sky(iA<0 & iB<0 & n_P(:,2)>0)=1;
%1.3 glass sees A but not B
VF_g2sky(iA>0 & iB<0)=((R21(iA>0 & iB<0)-R11(iA>0 & iB<0))./lface(row_g2sky((iA>0 & iB<0)))+1)/2;
%1.3 glass sees B but not A
VF_g2sky(iA<0 & iB>0)=((R12(iA<0 & iB>0)-R22(iA<0 & iB>0))./lface(row_g2sky((iA<0 & iB>0)))+1)/2;

% n_PA=[n_PA,zeros(size(row_g2sky))];
% sina=cross(n_PA,n_P); sina=sina(:,3);
% cosa=sum(n_PA.*n_P,2);
% alpha=atan2(sina,cosa);
% sina(alpha<-pi/2 | alpha>alpha_c)=-1; sina(alpha>pi/2 & alpha<alpha_c)=1;
% 
% n_PB=[n_PB,zeros(size(row_g2sky))];
% sinb=cross(n_PB,n_P); sinb=sinb(:,3);
% cosb=sum(n_PB.*n_P,2);
% beta=atan2(sinb,cosb);
% sinb(beta<-pi/2 & beta>-alpha_c)=-1; sinb(beta>pi/2 | beta<-alpha_c)=1;
% 
% VF_g2sky=(sinb-sina)/2;
% VF_g2sky=max(0,VF_g2sky);

% glass to reflector
row_g2ref=[row_g2sky;(N+1)*ones(size(row_g2sky));]; 
col_g2ref=[(N+1)*ones(size(row_g2sky));row_g2sky];
VF_g2ref=1-VF_g2sky;
VF_g2ref=[VF_g2ref;VF_g2ref.*lface(row_g2sky)/lcpc];

% glass to glass
igr=tag(row)==1; igc=tag(col)==1; igg=igr & igc; id_gg=find(igg);
n_PA=pmid(col(igg),:)-pmid(row(igg),:); PA=sqrt(sum(n_PA.^2,2)); n_PA=n_PA./[PA,PA];
n_P=-normvec(row(igg),:); cosP=sum(n_P.*n_PA,2);
iggs=(cosP<sqrt(r_g^2-r_a^2)/r_g);
n_A=-normvec(col(id_gg(iggs)),:); cosA=-sum(n_A.*n_PA(iggs,:),2);
if any(cosA<0), error('check'); end
% VF(id_gg(iggs))=cosP(iggs).*cosA.*lface(col(id_gg(iggs))).*atan(lns./PA(iggs))./PA(iggs);
% VF(id_gg(iggs))=cosP(iggs).*cosA.*lface(col(id_gg(iggs))).*pi/2./PA(iggs);
R11=sqrt(sum((pt(t(row(id_gg(iggs)),1),:)-pt(t(col(id_gg(iggs)),1),:)).^2,2));
R22=sqrt(sum((pt(t(row(id_gg(iggs)),2),:)-pt(t(col(id_gg(iggs)),2),:)).^2,2));
R21=sqrt(sum((pt(t(row(id_gg(iggs)),2),:)-pt(t(col(id_gg(iggs)),1),:)).^2,2));
R12=sqrt(sum((pt(t(row(id_gg(iggs)),1),:)-pt(t(col(id_gg(iggs)),2),:)).^2,2));
VF(id_gg(iggs))=abs(R11+R22-R12-R21)/2./lface(row(id_gg(iggs)));

% glass to absorber
igr=tag(row)==1; iac=(tag(col)==2 | tag(col)==-3 | tag(col)==-4);
iga=igr & iac; id_ga=find(iga);
n_PA=pmid(col(iga),:)-pmid(row(iga),:); PA=sqrt(sum(n_PA.^2,2)); n_PA=n_PA./[PA,PA];
n_A=normvec(col(iga),:); cosA=-sum(n_A.*n_PA,2); igas=cosA>0; % cosA>0 indicates seeing
% n_P=-normvec(row(id_ga(igas)),:); cosP=sum(n_P.*n_PA(igas,:),2); 
% VF(id_ga(igas))=cosP.*cosA(igas).*lface(col(id_ga(igas))).*atan(lns./PA(igas))./PA(igas);
% VF(id_ga(igas))=cosP.*cosA(igas).*lface(col(id_ga(igas))).*pi/2./PA(igas);
R11=sqrt(sum((pt(t(row(id_ga(igas)),1),:)-pt(t(col(id_ga(igas)),1),:)).^2,2));
R22=sqrt(sum((pt(t(row(id_ga(igas)),2),:)-pt(t(col(id_ga(igas)),2),:)).^2,2));
R21=sqrt(sum((pt(t(row(id_ga(igas)),2),:)-pt(t(col(id_ga(igas)),1),:)).^2,2));
R12=sqrt(sum((pt(t(row(id_ga(igas)),1),:)-pt(t(col(id_ga(igas)),2),:)).^2,2));
VF(id_ga(igas))=abs(R11+R22-R12-R21)/2./lface(row(id_ga(igas)));

% absorber to absorber
iar=tag(row)==2; iac=tag(col)==2; iaa=iar & iac; id_aa=find(iaa);
n_PA=pmid(col(iaa),:)-pmid(row(iaa),:); PA=sqrt(sum(n_PA.^2,2)); n_PA=n_PA./[PA,PA];
n_PC1=repmat([-(r_a-r_u),0],sum(iaa),1)-pmid(row(iaa),:); PC1=sqrt(sum(n_PC1.^2,2));
n_PC2=repmat([(r_a-r_u),0],sum(iaa),1)-pmid(row(iaa),:); PC2=sqrt(sum(n_PC2.^2,2));
PC1PA=sum(n_PA.*n_PC1,2); PC2PA=sum(n_PA.*n_PC2,2);
dC1=sqrt(PC1.^2-PC1PA.^2); dC2=sqrt(PC2.^2-PC2PA.^2);
iaas=(dC1>r_u & dC2>r_u);
n_P=-normvec(row(id_aa(iaas)),:); cosP=sum(n_P.*n_PA(iaas,:),2);
n_A=-normvec(col(id_aa(iaas)),:); cosA=-sum(n_A.*n_PA(iaas,:),2);
if any(cosP<0 | cosA<0), error('check'); end
% VF(id_aa(iaas))=cosP.*cosA.*lface(col(id_aa(iaas))).*atan(lns./PA(iaas))./PA(iaas);
% VF(id_aa(iaas))=cosP.*cosA.*lface(col(id_aa(iaas))).*pi/2./PA(iaas);
R11=sqrt(sum((pt(t(row(id_aa(iaas)),1),:)-pt(t(col(id_aa(iaas)),1),:)).^2,2));
R22=sqrt(sum((pt(t(row(id_aa(iaas)),2),:)-pt(t(col(id_aa(iaas)),2),:)).^2,2));
R21=sqrt(sum((pt(t(row(id_aa(iaas)),2),:)-pt(t(col(id_aa(iaas)),1),:)).^2,2));
R12=sqrt(sum((pt(t(row(id_aa(iaas)),1),:)-pt(t(col(id_aa(iaas)),2),:)).^2,2));
VF(id_aa(iaas))=abs(R11+R22-R12-R21)/2./lface(row(id_aa(iaas)));

% absorber to glass
iar=(tag(row)==2 | tag(row)==-3 | tag(row)==-4); igc=tag(col)==1; 
iag=iar & igc; id_ag=find(iag);
n_PA=pmid(col(iag),:)-pmid(row(iag),:); PA=sqrt(sum(n_PA.^2,2)); n_PA=n_PA./[PA,PA];
n_P=normvec(row(iag),:); cosP=sum(n_P.*n_PA,2); iags=cosP>0; % cosA>0 indicates seeing
% n_A=-normvec(col(id_ag(iags)),:); cosA=-sum(n_A.*n_PA(iags,:),2); 
% VF(id_ag(iags))=cosP(iags).*cosA.*lface(col(id_ag(iags))).*atan(lns./PA(iags))./PA(iags);
% VF(id_ag(iags))=cosP(iags).*cosA.*lface(col(id_ag(iags))).*pi/2./PA(iags);
R11=sqrt(sum((pt(t(row(id_ag(iags)),1),:)-pt(t(col(id_ag(iags)),1),:)).^2,2));
R22=sqrt(sum((pt(t(row(id_ag(iags)),2),:)-pt(t(col(id_ag(iags)),2),:)).^2,2));
R21=sqrt(sum((pt(t(row(id_ag(iags)),2),:)-pt(t(col(id_ag(iags)),1),:)).^2,2));
R12=sqrt(sum((pt(t(row(id_ag(iags)),1),:)-pt(t(col(id_ag(iags)),2),:)).^2,2));
VF(id_ag(iags))=abs(R11+R22-R12-R21)/2./lface(row(id_ag(iags)));

% absorber to u-tube (tag==3)
iar=tag(row)==2; iuc=tag(col)==3;
iau=iar & iuc; id_au=find(iau);
n_PA=pmid(col(iau),:)-pmid(row(iau),:); PA=sqrt(sum(n_PA.^2,2)); n_PA=n_PA./[PA,PA];
n_PC3=repmat([-(r_a-r_u),0],sum(iau),1)-pmid(row(iau),:);
n_PC4=repmat([(r_a-r_u),0],sum(iau),1)-pmid(row(iau),:); PC4=sqrt(sum(n_PC4.^2,2));
PC3PA=sum(n_PA.*n_PC3,2); 
PC4PA=sum(n_PA.*n_PC4,2); dC4=sqrt(PC4.^2-PC4PA.^2);
iaus=(PC3PA>PA) & (~(dC4<r_u & PC4PA<PA));
n_P=-normvec(row(id_au(iaus)),:); cosP=sum(n_P.*n_PA(iaus,:),2);
n_A=normvec(col(id_au(iaus)),:); cosA=-sum(n_A.*n_PA(iaus,:),2);
if any(cosP<0 | cosA<0), error('check'); end
% VF(id_au(iaus))=cosP.*cosA.*lface(col(id_au(iaus))).*atan(lns./PA(iaus))./PA(iaus);
% VF(id_au(iaus))=cosP.*cosA.*lface(col(id_au(iaus))).*pi/2./PA(iaus);
R11=sqrt(sum((pt(t(row(id_au(iaus)),1),:)-pt(t(col(id_au(iaus)),1),:)).^2,2));
R22=sqrt(sum((pt(t(row(id_au(iaus)),2),:)-pt(t(col(id_au(iaus)),2),:)).^2,2));
R21=sqrt(sum((pt(t(row(id_au(iaus)),2),:)-pt(t(col(id_au(iaus)),1),:)).^2,2));
R12=sqrt(sum((pt(t(row(id_au(iaus)),1),:)-pt(t(col(id_au(iaus)),2),:)).^2,2));
VF(id_au(iaus))=abs(R11+R22-R12-R21)/2./lface(row(id_au(iaus)));

% absorber to u-tube (tag==4)
iar=tag(row)==2; iuc=tag(col)==4;
iau=iar & iuc; id_au=find(iau);
n_PA=pmid(col(iau),:)-pmid(row(iau),:); PA=sqrt(sum(n_PA.^2,2)); n_PA=n_PA./[PA,PA];
n_PC3=repmat([-(r_a-r_u),0],sum(iau),1)-pmid(row(iau),:); PC3=sqrt(sum(n_PC3.^2,2));
n_PC4=repmat([(r_a-r_u),0],sum(iau),1)-pmid(row(iau),:);
PC3PA=sum(n_PA.*n_PC3,2); dC3=sqrt(PC3.^2-PC3PA.^2);
PC4PA=sum(n_PA.*n_PC4,2); 
iaus=(PC4PA>PA) & (~(dC3<r_u & PC3PA<PA));
n_P=-normvec(row(id_au(iaus)),:); cosP=sum(n_P.*n_PA(iaus,:),2);
n_A=normvec(col(id_au(iaus)),:); cosA=-sum(n_A.*n_PA(iaus,:),2);
if any(cosP<0 | cosA<0), error('check'); end
% VF(id_au(iaus))=cosP.*cosA.*lface(col(id_au(iaus))).*atan(lns./PA(iaus))./PA(iaus);
% VF(id_au(iaus))=cosP.*cosA.*lface(col(id_au(iaus))).*pi/2./PA(iaus);
R11=sqrt(sum((pt(t(row(id_au(iaus)),1),:)-pt(t(col(id_au(iaus)),1),:)).^2,2));
R22=sqrt(sum((pt(t(row(id_au(iaus)),2),:)-pt(t(col(id_au(iaus)),2),:)).^2,2));
R21=sqrt(sum((pt(t(row(id_au(iaus)),2),:)-pt(t(col(id_au(iaus)),1),:)).^2,2));
R12=sqrt(sum((pt(t(row(id_au(iaus)),1),:)-pt(t(col(id_au(iaus)),2),:)).^2,2));
VF(id_au(iaus))=abs(R11+R22-R12-R21)/2./lface(row(id_au(iaus)));

% u-tube to u-tube
iur=(tag(row)==3 | tag(row)==4); iuc=(tag(col)==3 | tag(col)==4);
iuu=iur & iuc; id_uu=find(iuu);
n_PA=pmid(col(iuu),:)-pmid(row(iuu),:); PA=sqrt(sum(n_PA.^2,2)); n_PA=n_PA./[PA,PA];
n_P=normvec(row(iuu),:); cosP=sum(n_P.*n_PA,2);
n_A=normvec(col(iuu),:); cosA=-sum(n_A.*n_PA,2); 
iuus=cosP>0 & cosA>0;
% VF(id_uu(iuus))=cosP(iuus).*cosA(iuus).*lface(col(id_uu(iuus))).*atan(lns./PA(iuus))./PA(iuus);
% VF(id_uu(iuus))=cosP(iuus).*cosA(iuus).*lface(col(id_uu(iuus))).*pi/2./PA(iuus);
R11=sqrt(sum((pt(t(row(id_uu(iuus)),1),:)-pt(t(col(id_uu(iuus)),1),:)).^2,2));
R22=sqrt(sum((pt(t(row(id_uu(iuus)),2),:)-pt(t(col(id_uu(iuus)),2),:)).^2,2));
R21=sqrt(sum((pt(t(row(id_uu(iuus)),2),:)-pt(t(col(id_uu(iuus)),1),:)).^2,2));
R12=sqrt(sum((pt(t(row(id_uu(iuus)),1),:)-pt(t(col(id_uu(iuus)),2),:)).^2,2));
VF(id_uu(iuus))=abs(R11+R22-R12-R21)/2./lface(row(id_uu(iuus)));

% u-tube (tag==3) to absorber
iur=tag(row)==3; iac=tag(col)==2;
iua=iur & iac; id_ua=find(iua);
n_PA=pmid(col(iua),:)-pmid(row(iua),:); PA=sqrt(sum(n_PA.^2,2)); n_PA=n_PA./[PA,PA];
n_PC=repmat([(r_a-r_u),0],sum(iua),1)-pmid(row(iua),:); PC=sqrt(sum(n_PC.^2,2));
PCPA=sum(n_PA.*n_PC,2); dC=sqrt(PC.^2-PCPA.^2);
n_P=normvec(row(iua),:); cosP=sum(n_PA.*n_P,2);
iuas=cosP>0 & dC>r_u;
n_A=normvec(col(id_ua(iuas)),:); cosA=-sum(n_A.*n_PA(iuas,:),2);
if any(cosA<0), error('check'); end
% VF(id_ua(iuas))=cosP(iuas).*cosA.*lface(col(id_ua(iuas))).*atan(lns./PA(iuas))./PA(iuas);
% VF(id_ua(iuas))=cosP(iuas).*cosA.*lface(col(id_ua(iuas))).*pi/2./PA(iuas);
R11=sqrt(sum((pt(t(row(id_ua(iuas)),1),:)-pt(t(col(id_ua(iuas)),1),:)).^2,2));
R22=sqrt(sum((pt(t(row(id_ua(iuas)),2),:)-pt(t(col(id_ua(iuas)),2),:)).^2,2));
R21=sqrt(sum((pt(t(row(id_ua(iuas)),2),:)-pt(t(col(id_ua(iuas)),1),:)).^2,2));
R12=sqrt(sum((pt(t(row(id_ua(iuas)),1),:)-pt(t(col(id_ua(iuas)),2),:)).^2,2));
VF(id_ua(iuas))=abs(R11+R22-R12-R21)/2./lface(row(id_ua(iuas)));

% u-tube (tag==4) to absorber
iur=tag(row)==4; iac=tag(col)==2;
iua=iur & iac; id_ua=find(iua);
n_PA=pmid(col(iua),:)-pmid(row(iua),:); PA=sqrt(sum(n_PA.^2,2)); n_PA=n_PA./[PA,PA];
n_PC=repmat([-(r_a-r_u),0],sum(iua),1)-pmid(row(iua),:); PC=sqrt(sum(n_PC.^2,2));
PCPA=sum(n_PA.*n_PC,2); dC=sqrt(PC.^2-PCPA.^2);
n_P=normvec(row(iua),:); cosP=sum(n_PA.*n_P,2);
iuas=cosP>0 & dC>r_u;
n_A=normvec(col(id_ua(iuas)),:); cosA=-sum(n_A.*n_PA(iuas,:),2);
if any(cosA<0), error('check'); end
% VF(id_ua(iuas))=cosP(iuas).*cosA.*lface(col(id_ua(iuas))).*atan(lns./PA(iuas))./PA(iuas);
% VF(id_ua(iuas))=cosP(iuas).*cosA.*lface(col(id_ua(iuas))).*pi/2./PA(iuas);
R11=sqrt(sum((pt(t(row(id_ua(iuas)),1),:)-pt(t(col(id_ua(iuas)),1),:)).^2,2));
R22=sqrt(sum((pt(t(row(id_ua(iuas)),2),:)-pt(t(col(id_ua(iuas)),2),:)).^2,2));
R21=sqrt(sum((pt(t(row(id_ua(iuas)),2),:)-pt(t(col(id_ua(iuas)),1),:)).^2,2));
R12=sqrt(sum((pt(t(row(id_ua(iuas)),1),:)-pt(t(col(id_ua(iuas)),2),:)).^2,2));
VF(id_ua(iuas))=abs(R11+R22-R12-R21)/2./lface(row(id_ua(iuas)));

% Integrate
row=row(VF~=0); col=col(VF~=0); VF=VF(VF~=0);
VF=sparse([row;col],[col;row],[VF;VF.*lface(row)./lface(col)],N,N);
[row,col,VF]=find(VF);
% % glass normalization
% ig=tag(row)==1; rowt=row(ig); colt=col(ig); VFt=VF(ig);
% VFt=sparse(rowt,colt,VFt,N,N); VFtot=sum(VFt,2);
% VF(ig)=VF(ig)./VFtot(rowt);
% % absorber and utube normalization
% iar=(tag(row)==2 | tag(row)==-3 | tag(row)==-4); igc=tag(col)==1; iag=iar & igc;
% rowt=row(iag); colt=col(iag); VFt=VF(iag);
% VFt=sparse(rowt,colt,VFt,N,N); VFtot=sum(VFt,2);
% VF(iag)=VF(iag)./VFtot(rowt);
% iar=(tag(row)==2 | tag(row)==3 | tag(row)==4); 
% iac=(tag(col)==2 | tag(col)==3 | tag(col)==4); iaa=iar & iac;
% rowt=row(iaa); colt=col(iaa); VFt=VF(iaa);
% VFt=sparse(rowt,colt,VFt,N,N); VFtot=sum(VFt,2);
% VF(iaa)=VF(iaa)./VFtot(rowt);
VF=sparse([row;row_g2ref;row_g2sky],[col;col_g2ref;col_g2sky],[VF;VF_g2ref;VF_g2sky],N+1,N+2);
VF(N+1,N+2)=1-sum(VF(N+1,1:N+1)); % reflector to sky
end