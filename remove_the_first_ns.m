function [dat_remove_n_s,fields_remove_ns]=remove_the_first_ns(dat,fields,n)
%this function will remove the first n seconds for the dat and dat_grouped

et=find(matches(fields,'et'));
AN=find(matches(fields,'AN'));
xs=find(matches(fields,'xspine'));
ys=find(matches(fields,'yspine'));

w=1;
for i=1:length(dat{et,1})
    idx=find(dat{et,1}{i,1}>=n);
    if isempty(idx)
        clear idx
        continue
    else
        for j=1:length(dat)
            if j==AN;
                dat_remove_n_s{j,1}{w,1}=dat{j,1}{i,1};
            elseif j==xs|j==ys
                dat_remove_n_s{j,1}{w,1}=dat{j,1}{i,1}(idx,:);
            else 
                dat_remove_n_s{j,1}{w,1}=dat{j,1}{i,1}(idx);
            end
        end
        w=w+1;
    end
end

fields_remove_ns=fields;
end