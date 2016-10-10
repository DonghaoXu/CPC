function y=VisualCPC(ax,T,pt,t,L)
%check
if size(T,1)~=size(t,1)
    error('Temperature and connectivity have different sizes'); 
end

if isempty(ax)
    figure
    ax=gca;
end

N=size(T,2); Np=size(pt,1); lns=L/N;
Tmin=min(min(T)); Tmax=max(max(T));
set(gcf,'color',[1 1 1]);
ptx=repmat(pt(:,1),N+1,1);
pty=ones(Np,1)*lns*(0:N); pty=pty(:);
ptz=repmat(pt(:,2),N+1,1);
pn=[ptx,pty,ptz];
t=repmat(t,N,1)+reshape(repmat([(0:N-1),(0:N-1)]*Np,size(t,1),1),[],2);
tn=[t,t(:,2)+Np,t(:,1)+Np];
lam=(T(:)-Tmin)/(Tmax-Tmin);
% pcol=[max(0,2*lam-1),1-abs(2*lam-1),max(0,1-2*lam)];
pcol = lam * [1 1 1];
y=patch('vertices',pn,'faces',tn,...
    'edgecolor','none','facecolor','flat',...
    'facevertexcdata',pcol,'parent',ax);
view(ax,40,50)
colorbarlabel=(Tmin:(Tmax-Tmin)/10:Tmax)';
colorbarlabel=cellstr(num2str(colorbarlabel,'%.2f\n'));
colorbar(ax,'yticklabel',colorbarlabel)
lam=(0:0.001:1)';
% map=[max(0,2*lam-1),1-abs(2*lam-1),max(0,1-2*lam)];
map = lam * [1 1 1];
colormap(ax,map)
axis(ax,'off')
end