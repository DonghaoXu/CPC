function [it1,it2]=getInter(p)
%return the intersection with y=0 for a curve identified by vector x
% global s R r

[~,n]=size(p);
if n~=2
    fprintf('\nWrong dimensions in profile vector\n');
    return;
end

x=p(:,1); y=p(:,2);
ysign=sign(y);
ind=find(diff(ysign));

it1=[inf,0];
it2=it1;
if ~isempty(ind)
    id=x(ind)>0 & x(ind+1)>0;
    if any(id)
        x=[x(ind(id)),x(ind(id)+1)];
        y=[y(ind(id)),y(ind(id)+1)];
        [~,i]=min(x(:,1));
        x=x(i,:);
        y=y(i,:);
        it1=[x(1),y(1)];
        it2=[x(2),y(2)];
    end
end


% if m==0
%     y=[inf,0];
%     z=y;
%     return;
% end

% ind=1;
% xp=[x(ind:end,1),sign(x(ind,2))*x(ind:end,2)];
% ind=find(xp(:,2)<=0,1,'first');
% k=1;
% sind=0;
% y=[inf,0;inf,0;inf,0];
% z=y;
% while ~isempty(ind)
%     y(k,:)=x(ind+sind-1,:);
%     z(k,:)=x(ind+sind,:);
%     k=k+1;
%     sind=ind+sind-1;
%     xp=[xp(ind:end,1),sign(xp(ind,2))*xp(ind:end,2)];
%     xp(1,2)=xp(1,2)+10e-5;
%     ind=find(xp(:,2)<=0,1,'first');    
% end
% 
% ind=y(:,1)>0 & z(:,1)>0;
% if any(ind)
%     y=y(ind,:);
%     z=z(ind,:);
%     [~,ind]=min(y(:,1));
%     y=y(ind,:);
%     z=z(ind,:);
% else
%     y=[inf,0];
%     z=y;
% end




% xp=[x(:,1),sign(x(1,2))*x(:,2)];
% ind1=find(xp(:,2)<=0, 1, 'first');
% 
% if isempty(ind1)
%     y=[inf,0];
%     z=y;
% elseif id==flag
%     distance=sum(xp.^2,2);
%     [~,ind2]=min(distance);
%     if abs(ind2-ind1)<2
%         ind1=find(xp(ind1:end,2)>0,1,'first');
%         if isempty(ind1)
%             y=[inf,0];
%             z=y;
%         else
%             y=x(ind1-1,:);
%             z=x(ind1,:);      
%         end
%     else
%         y=x(ind1-1,:);
%         z=x(ind1,:); 
%     end  
% else
%     y=x(ind1-1,:);
% 	z=x(ind1,:);
% %    switch(flag)
% %        case 5
% %            y=[x(ind1,1)-(x(ind1,1)-x(ind1-1,1))/(x(ind1,2)-x(ind1-1,2))*x(ind1,2),0];
% %            z=y;
% %         case 6
% %             rs=-s*R;
% %             y=[rs(1)-sqrt(r^2-rs(2)^2),0];
% %             z=y;
% %        otherwise
% %            y=x(ind1-1,:);
% %            z=x(ind1,:);
% %    end
%     %y1=x(ind1,1);
%     %-(x(ind1,1)-x(ind1-1,1))/(x(ind1,2)-x(ind1-1,2))*x(ind1,2);
%     %y2=x(ind2,1);
%     %-(x(ind2,1)-x(ind2+1,1))/(x(ind2,2)-x(ind2+1,2))*x(ind2,2);
%     %y=min(y1,y2);
% end

end