function y=checktypequiver(varargin)
if nargin==0
    error('Get me some inputs, dude! I need points, connectivity, normal vector, and tags, or a mesh struct!')
end
if isa(varargin{1},'struct')
    pt=varargin{1}.points;
    t=varargin{1}.connectivity;
    tag=varargin{1}.tag.tag;
    normvec=varargin{1}.normvec;
elseif isa(varargin{1},'numeric') && numel(varargin)==4
    pt=varargin{1};
    t=varargin{2};
    tag=varargin{3};
    normvec=varargin{4};
end
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
    quiver(pt(t(i,1),1),pt(t(i,1),2),...
        pt(t(i,2),1)-pt(t(i,1),1),pt(t(i,2),2)-pt(t(i,1),2),...
        'color',color(j,:));
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