function deg=get_pre_or_post_turning_direction(vector,dis,i,j,idx0,dat_JB,pre,post)
if isequal(pre,'y')

    x1=dat_JB.x{i,1}(idx0(j));
    y1=dat_JB.y{i,1}(idx0(j));
    if idx0(j)>dis
        x2=dat_JB.x{i,1}(idx0(j)-dis+1);
        y2=dat_JB.y{i,1}(idx0(j)-dis+1);
    else
        x2=dat_JB.x{i,1}(1);
        y2=dat_JB.y{i,1}(1);
    end
    x3=x1-x2;
    y3=y1-y2;

    deg=atan2d(y3.*(vector(1))-vector(2).*x3,x3.*(vector(1))+y3.*vector(2));
    deg=0-deg;
elseif isequal(post,'y')
    len=length(dat_JB.x{i,1});
    x2=dat_JB.x{i,1}(idx0(j));
    y2=dat_JB.y{i,1}(idx0(j));
    if  idx0(j)<len-dis+1
        x1=dat_JB.x{i,1}(idx0(j)+dis-1);
        y1=dat_JB.y{i,1}(idx0(j)+dis-1);
    else
        x1=dat_JB.x{i,1}(len);
        y1=dat_JB.y{i,1}(len);
    end
    x3=x1-x2;
    y3=y1-y2;

    deg=atan2d(y3.*(vector(1))-vector(2).*x3,x3.*(vector(1))+y3.*vector(2));
    deg=0-deg;
end
end
