% if the function is not recognized, make sure the corresponding .m file is
% in the directory
function [dat,xmax,ymax]=vertcat_all_data (dat)
fn=fieldnames(dat);
% input is the cell array extracted from the dat file
% this function will just make them in to a long array
for i =1:length(fn)
    dat.(fn{i})=vertcat(dat.(fn{i}){:});
end

xmax=max(dat.x);
ymax=max(dat.y);
%if odor is in the left side, we need to flip it left and right in order to
%plot the odor in the right side
dat.x=xmax-dat.x;

end

