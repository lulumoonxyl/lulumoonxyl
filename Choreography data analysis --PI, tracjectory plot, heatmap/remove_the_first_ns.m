function [dat_remove_n_s,dat_grouped_remove_n_s,fields_remove_ns,fields_grouped_remove_ns]=remove_the_first_ns(dat,dat_grouped,fields,fields_grouped,n)
%this function will remove the first n seconds for the dat and dat_grouped
et_pos=find(matches(fields,'et'));
idx=find(dat{et_pos,1}>=n);


for i=1:length(dat)
    dat_remove_n_s{i,1}=dat{i,1}(idx);
end

clear idx
et_pos1=find(matches(fields_grouped,'et'));
w=1;
for i=1:length(dat_grouped{et_pos1,1})
    idx=find(dat_grouped{et_pos1,1}{i,1}>=n);
    if isempty(idx)
        clear idx
        continue
    else
        for j=1:length(dat_grouped)
            dat_grouped_remove_n_s{j,1}{w,1}=dat_grouped{j,1}{i,1}(idx);
            
        end
        w=w+1;
    end
end
fields_remove_ns=fields;
fields_grouped_remove_ns=fields_grouped;
end