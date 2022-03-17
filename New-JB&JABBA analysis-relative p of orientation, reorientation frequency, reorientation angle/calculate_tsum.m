function [dat_JB,dat_JAABA]=calculate_tsum(dat_JB,dat_JAABA,n)
%this function will get the tsum for each larvae and remove the cells that
%have tracking time less than n s

for i=1:length(dat_JB.et)
    dat_JB.tsum(i,1)=dat_JB.et{i,1}(end)-dat_JB.et{i,1}(1);
end
if ~isempty(n)
    idx=find(dat_JB.tsum<n);
    fn=fieldnames(dat_JB);
    fn_JAABA=fieldnames(dat_JAABA);
    for j=1:length(fn)
        dat_JB.(fn{j})(idx)=[];
    end
     for j=1:length(fn_JAABA)
        dat_JAABA.(fn_JAABA{j})(idx)=[];
    end
end

end