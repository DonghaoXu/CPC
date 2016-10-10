function y=checktype(pt,t,tag,normvec)
color=[251,178,23;
    1,77,103;
    237,222,139;
    96,143,159;
    254,67,101;
    252,157,154;
    161,23,21;
    118,77,57;
    3,101,100;
    107,194,53];
color=color/255;
figure;
hold on
ut=unique(tag);
pmid=(pt(t(:,1),:)+pt(t(:,2),:))/2;
for j=1:size(ut,1)
    i=find(tag==ut(j));
    for k=1:size(i,1)
        x=[pt(t(i(k),1),1),pt(t(i(k),2),1)];
        y=[pt(t(i(k),1),2),pt(t(i(k),2),2)];
        plot(x,y,'color',color(j,:));
    end
    quiver(pmid(i,1),pmid(i,2),normvec(i,1),normvec(i,2),'color',color(j,:));
end
axis equal
axis off
set(gcf,'color',[1 1 1])
y=gcf;
% fprintf('Press any key to continue...\n')
% pause
% close
end