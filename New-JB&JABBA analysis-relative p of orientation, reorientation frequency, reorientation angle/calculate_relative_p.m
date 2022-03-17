function dat_JB=calculate_relative_p(dat_JB,bins)
%this function will calculate the relative p of orientation for each larva
%we set bins=20; bins is the range for x position
theta=-180:bins:180;
for j=1:length(dat_JB.orientation)
    for i=1:length(theta)-1
        ind1=find(dat_JB.orientation{j,1}>=theta(i)&dat_JB.orientation{j,1}<theta(i+1));
        if isempty(ind1)
            dat_JB.p_each_larva(i,j)=0;
            continue
        elseif ind1(1)==1
            ind1(1)=[];
        end
        t_abs_larva=sum(dat_JB.et{j,1}(ind1)-dat_JB.et{j,1}(ind1-1));
        dat_JB.p_each_larva(i,j)=t_abs_larva/dat_JB.tsum(j);
        clear t_abs_larva ind1
    end
end
    dat_JB.p(:,1)=mean(dat_JB.p_each_larva,2);
    dat_JB.p(:,2)=std(dat_JB.p_each_larva,[],2)./sqrt(length(dat_JB.p_each_larva));
end

