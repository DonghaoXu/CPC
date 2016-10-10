function FigureFormat(width,height,Fontsize)
leftinch=0.1;
bottominch=0.1;

set(gcf,'paperposition',[leftinch,bottominch,width,height])
set(gcf,'papersize',[width+leftinch*2,height+bottominch*2])
set(findall(gcf,'-property','FontSize'),...
    'FontSize',Fontsize)
set(findall(gcf,'-property','FontName'),...
    'FontName','Cambria')
set(findall(gcf,'-property','linewidth'),...
    'linewidth',0.1)
set(findall(gcf,'-property','MarkerSize'),...
    'MarkerSize',3)