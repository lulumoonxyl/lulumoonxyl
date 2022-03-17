function data_split_q=get_data_for_each_quadrant(data_split,str)
fn=fieldnames(data_split);
idx=find(contains(fn,str));
fn_q=fn(idx);
for i=1:length(fn_q)
    data_split_q.(fn_q{i})=data_split.(fn_q{i});
end 
end 