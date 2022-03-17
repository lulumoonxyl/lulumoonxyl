function dat_JAABA=arrange_JAABA_data(dat_JAABA)
%this function will rearrange the JAABA data array for the first several n
%cells
fn=fieldnames(dat_JAABA);
for i=1:length(fn)
    dat_JAABA.(fn{i})=vertcat(dat_JAABA.(fn{i}){:});
for j=1:length(dat_JAABA.(fn{i}))
    dat_JAABA.(fn{i}){j,1}= dat_JAABA.(fn{i}){j,1}';
end

end
clear i j
end