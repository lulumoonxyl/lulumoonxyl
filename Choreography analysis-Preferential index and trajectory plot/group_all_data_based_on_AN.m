function dat_grouped=group_all_data_based_on_AN(dat)
%this function will groupe the data based on their animal number 
fn=fieldnames(dat);
idx=find(diff(dat.AN)~=0);
idx=[0;idx;length(dat.AN)];
for j=1:length(fn)
    for i=1:length(idx)-1
        dat_grouped.(fn{j}){i,1}=dat.(fn{j})(idx(i)+1:idx(i+1));
    end
end
clear idx i j
end