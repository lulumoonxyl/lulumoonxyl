function [p_timeseries,PI]=cal_navindex_vs_t(data,et,series,edge)
%this function will calculate the preferential index: (number of larvae
%close to the odor-number of larvae away from the odor)/total number of
%larvare for the timeseries [0:900]
%edge is how many sections we want to divide the x pos
len=length(series)-1;
all_larva=zeros(len,1);larva_region=zeros(len,length(edge)-1);
for i=1:len

    for j=1:length(et)
        if i==len
            idx=find(et{j,1}>=series(i)&et{j,1}<=series(i+1));
        else
            idx=find(et{j,1}>=series(i)&et{j,1}<series(i+1));
            if isempty(idx)
                continue
            end
            % we divide the plan into three regions (based on the x values)
            % namly 1(far away),2,3(closet to odor);2 is sorta a gray region

            N=histcounts(data{j,1}(idx),edge);
            [M,I]=max(N);
            all_larva(i,1)=all_larva(i,1)+1;
            larva_region(i,I)=larva_region(i,I)+1;
            clear N M I
        end
    end

    p_timeseries=larva_region./all_larva;
    p=larva_region(:,3)./(larva_region(:,1)+larva_region(:,3));
    %calculate the preferential index using # in region 1-# in region
    %2/region1+2
    PI(:,1)=(larva_region(:,3)-larva_region(:,1))./(larva_region(:,1)+larva_region(:,3));
    PI(:,2)=1.96.*sqrt(p.*(1-p)./(larva_region(:,1)+larva_region(:,3)));
end