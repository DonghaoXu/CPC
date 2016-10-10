function cmap=mycmap(varargin)
cmap=cell(5,1);
cmap(1)={[200,212,236;...
    137,171,217;...
    1,104,179;...
    0,90,160;...
    1,73,132]/256};
cmap(2)={[213,233,206;...
    164,213,148;...
    18,170,89;...
    16,148,49;...
    6,123,42]/256};
cmap(3)={[254,240,203;...
    254,223,143;...
    250,191,1;...
    215,165,4;...
    192,140,5]/256};
cmap(4)={[255,233,201;...
    249,205,134;...
    243,152,1;...
    211,130,0;...
    174,105,1]/256};
cmap(5)={[248,208,200;...
    243,154,136;...
    230,1,44;...
    200,0,37;...
    164,0,27]/256};
if ~isempty(varargin)
    temp=vertcat(cmap{:});
    colormap(temp)
    t=pi/4:pi/2:2*pi+pi/4;
    figure; axis equal; hold on
    for i=1:size(temp,1)
        x=2*ceil(i/size(cmap{1,1},1))-1+cos(t);
        y=2*mod(i-1,size(cmap{1,1},1))+sin(t);
        fill(x,y,temp(i,:))
    end
end
end