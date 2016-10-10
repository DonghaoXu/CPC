function [history,time]=RayTracing(yl2,yl1,yr2,yr1,r_g,t_g,r_a,gap,l,iai,tol,theta)
%l=size(sp,1)+1;
nray=10;
iai=iai(:);
na=size(iai,1);
history.x=zeros(l-1,nray,na);
history.y=zeros(l-1,nray,na);
history.x(:,:,:)=NaN;
history.y(:,:,:)=NaN;
history.index=zeros(l-1,nray,na);
history.ia=zeros(l-1,nray,na);
history.ia(:,:,:)=NaN;
time=zeros(na,1);
topm=[yl2(end,:);yl1(1,:)];
[~,ind]=max(topm(:,2));
top=topm(ind,:);
L=norm(top)+1;
lam=(0:1/l:1)';
lam=lam(2:end-1);

%matlabpool
for j=1:na
    A=[-sin(iai(j)/180*pi),cos(iai(j)/180*pi);cos(iai(j)/180*pi),sin(iai(j)/180*pi)];
    b_rec=[L;r_g];
    b_rl=[L;top(1)*cos(iai(j)/180*pi)+top(2)*sin(iai(j)/180*pi)];
    b_rr=[L;-top(1)*cos(iai(j)/180*pi)+top(2)*sin(iai(j)/180*pi)];
    %b_rl=[L;yl2(end,1)*cos(iai(j)/180*pi)+yl2(end,2)*sin(iai(j)/180*pi)];
    %b_rr=[L;yr2(end,1)*cos(iai(j)/180*pi)+yr2(end,2)*sin(iai(j)/180*pi)];
    coor_rec=A\b_rec;
    coor_rl=A\b_rl;
    coor_rr=A\b_rr;
    sp1=coor_rl;
    if coor_rec(1)>coor_rr(1)
        sp2=coor_rec;
    else
        sp2=coor_rr;
    end
    spa=lam*sp1'+(1-lam)*sp2';
    inix=spa(:,1);
    iniy=spa(:,2);
    va=[sin(iai(j)/180*pi),-cos(iai(j)/180*pi)];
    history.x(:,1,j)=spa(:,1);
    history.y(:,1,j)=spa(:,2);
    history.ia(:,1,j)=iai(j);
    tempx=history.x(:,:,j);
    tempy=history.y(:,:,j);
    tempi=history.index(:,:,j);
    tempia=history.ia(:,:,j);
    ia1=iai(j);
    fprintf('Incident angle: %.1f, %d/%d... ',iai(j),j,na)
    tic;
    parfor i=1:l-1
        %i
        oi=NaN;
        k=1;
        vra=va;
        in=[inix(i),iniy(i)];
        tx=zeros(1,nray);
        ty=zeros(1,nray);
        ti=zeros(1,nray);
        tia=zeros(1,nray);
        tx(:)=NaN;
        ty(:)=NaN;
        ia=ia1;
        while oi~=0 && k<=nray
            %k
            tx(k)=in(1);
            ty(k)=in(2);
            ti(k)=oi;
            tia(k)=real(ia);
            %oi=abs(oi);
            if k~=nray
                [vra,in,ia,oi]=refl(vra,in,yl2,yl1,yr2,yr1,r_g,t_g,r_a,gap,tol,theta);
            end    
            k=k+1;
        end
        tempx(i,:)=tx;
        tempy(i,:)=ty;
        tempi(i,:)=ti;
        tempia(i,:)=tia;
    end
    tcom=toc;
    time(j)=toc;
    fprintf('Time elapsed: %.4fs\n',tcom)
    history.x(:,:,j)=tempx;
    history.y(:,:,j)=tempy;
    tempi(:,1)=0;
    history.index(:,:,j)=tempi;
    history.ia(:,:,j)=tempia;
end
%matlabpool close
end