function [dat_grouped,PI_timeseries]=calculate_preferential_index_vs_time(dat_grouped,ii,time,xmax,PI_timeseries)
%this function will calculate the preferential index: (number of larvae
%close to the odor-number of larvae away from the odor)/total number of
%larvare for the timeseries [0:900]

for j=1:length(dat_grouped.et)
    for i=1:length(time)-1
        idx=find(dat_grouped.et{j,1}>=time(i)&dat_grouped.et{j,1}<time(i+1));
        idx1=find(dat_grouped.x{j,1}(idx)>xmax/2);
        if ~isempty(idx)
            all_larvae_in_t(i,j)=1;
        else
            all_larvae_in_t(i,j)=0;
        end
        if ~isempty(idx1)
            larvae_close_to_odor(i,j)=1;
        else
            larvae_close_to_odor(i,j)=0;
        end
    end
end
total=sum(all_larvae_in_t,2);
close=sum(larvae_close_to_odor,2);
further=total-close;
PI=(close-further)./total;
%get the 95% CI for plot later
p=close./total;
CI=1.96.*sqrt(p.*(1-p)./total);
PI_timeseries{ii,1}=PI;
PI_timeseries{ii,2}=CI;
end