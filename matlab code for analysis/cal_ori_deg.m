function deg=cal_ori_deg(x1,x2,y1,y2,vector)
%vector is the odor position
x3=x1-x2;
y3=y1-y2;
if height(vector)==1
    deg=atan2d(y3.*(vector(1))-vector(2).*x3,x3.*(vector(1))+y3.*vector(2));
else
    deg=atan2d(y3.*(vector(:,1))-vector(:,2).*x3,x3.*(vector(:,1))+y3.*vector(:,2));
end

end