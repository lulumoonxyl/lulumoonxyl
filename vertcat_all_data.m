% if the function is not recognized, make sure the corresponding .m file is
% in the directory
function [dat,xmax,ymax]=vertcat_all_data (dat,fields)

% input is the cell array extracted from the dat file
% this function will just make them in to a long array
for i =1:length(dat)
    dat{i,1}=vertcat(dat{i,1}{:});
end
x_pos=find(matches(fields,'x'));
y_pos=find(matches(fields,'y'));
xmax=max(dat{x_pos,1});
ymax=max(dat{y_pos,1});
%if odor is in the left side
dat{x_pos,1}=xmax-dat{x_pos,1};

end

