function [t_abs,t_p,x]=calculate_relative_p(dat,fields,bins)
%this function will calculate the relative p for each larva
%may delete some larvae with t_sum larger than ns
%we set bins=20; bins is the range for x position 
tsum_pos=find(matches(fields,'tsum'));
theta_pos=find(matches(fields,'heading direction'));
et_pos=find(matches(fields,'et'));
% x just means the x axis for heading direction which ranges from 0-360 degrees

x=-180:bins:180;


for j=1:length(dat{theta_pos,1})
    for i=1:length(x)-1
        ind1=find(dat{theta_pos,1}{j,1}>=x(i)&dat{theta_pos,1}{j,1}<x(i+1));
        if isempty(ind1)
            t_abs(i,j)=0;
            t_p(i,j)=0;
            continue 
        elseif ind1(1)==1
            ind1(1)=[];
        end 
        t_abs_larva=sum(dat{et_pos,1}{j,1}(ind1)-dat{et_pos,1}{j,1}(ind1-1));
        t_abs(i,j)=t_abs_larva;
        t_p(i,j)=t_abs_larva/dat{tsum_pos,1}(j);
        clear t_abs_larva ind1
    end
    
end