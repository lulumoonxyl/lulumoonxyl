function [dat_grouped,fields_grouped]=group_all_data_based_on_AN(dat,fields)
ind=find(matches(fields,'AN'));
idx=find(diff(dat{ind,1})~=0);
idx=[0;idx;length(dat{ind,1})];
for j=1:length(dat)
    for i=1:length(idx)-1
        dat_grouped{j,1}{i,1}=dat{j,1}(idx(i)+1:idx(i+1));
    end
end
fields_grouped=fields;
clear idx i j
end