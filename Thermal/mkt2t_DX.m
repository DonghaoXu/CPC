function [row,col]=mkt2t_DX(t)
% Smart
N=size(t,1);
tt=repmat(t(:),1,N);
id=repmat((1:N)',2,N);

t1=repmat(t(:,1)',2*N,1);
t2=repmat(t(:,2)',2*N,1);

ihalf=(id>repmat(1:N,2*N,1));

id1=(tt==t1 & ihalf);
id2=(tt==t2 & ihalf);

[~,row1]=find(id1);
col1=id(id1);

[~,row2]=find(id2);
col2=id(id2);

row=[row1;row2];
col=[col1;col2];
end