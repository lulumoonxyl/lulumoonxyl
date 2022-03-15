function dat_grouped=eliminate_tt(dat_grouped,n)
%this function will eliminate trajectory with tracking time less than n=10s
if isempty(n)==0
fn=fieldnames(dat_grouped);
for i =1:length(dat_grouped.et)
    %calculate the tracking time for each trajectory 
    tsum(i,1)=dat_grouped.et{i,1}(end)-dat_grouped.et{i,1}(1);
end
%find the trajectory with tracking time less than n seconds
idx=find(tsum<n);
if isempty(idx)==0
    for j=1:length(fn)
        %eliminate the trajectory with tracking time less than n seconds
        dat_grouped.(fn{j})(idx)=[];
    end 
end 
end 
 
end 