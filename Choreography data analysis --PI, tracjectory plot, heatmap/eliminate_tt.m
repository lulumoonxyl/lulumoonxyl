function [fields_grouped,dat_grouped]=eliminate_tt(fields_grouped,dat_grouped,n)
if isempty(n)==0
t=find(matches(fields_grouped,'et'));
for i =1:length(dat_grouped{t,1})
    tsum(i,1)=dat_grouped{t,1}{i,1}(end)-dat_grouped{t,1}{i,1}(1);
end

idx=find(tsum<n);
if isempty(idx)==0
    for j=1:length(dat_grouped)
        dat_grouped{j,1}(idx)=[]
    end 
end 
end 
 
end 