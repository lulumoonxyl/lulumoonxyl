function [dat_edge_removed,dat_grouped_edge_removed,fields_edge_removed,fields_grouped_edge_removed]=eliminate_edge(dat, dat_grouped,fields,fields_grouped,edge_x_y)
%this will remove the edge of all the data
%Note: dat_edge_removed can only be used to draw heatmap or trajectory
% dat_grouped_edge_removed can be used to calculate PI and tt
%edge_x_y=[xmin,xmax;ymin,ymax]-->[5,222;25,250]
if isempty(dat)~=1&isempty(fields)~=1
x1=find(matches(fields,'x'));
y1=find(matches(fields,'y'));

idx=find(dat{y1,1}>edge_x_y(2,1)&dat{y1,1}<edge_x_y(2,2)&dat{x1,1}>edge_x_y(1,1)&dat{x1,1}<edge_x_y(1,2));
for i=1:length(dat)
    dat_edge_removed{i,1}=dat{i,1}(idx);
end
clear idx
fields_edge_removed=fields;
else 
    dat_edge_removed={};
    fields_edge_removed={};
end 

x2=find(matches(fields_grouped,'x'));
y2=find(matches(fields_grouped,'y'));
w=1;
dat_grouped_edge_removed={};

for i=1:length(dat_grouped{x2,1})
    idx=find(dat_grouped{y2,1}{i,1}>edge_x_y(2,1)&dat_grouped{y2,1}{i,1}<edge_x_y(2,2)&dat_grouped{x2,1}{i,1}>edge_x_y(1,1)&dat_grouped{x2,1}{i,1}<edge_x_y(1,2));
    if isempty(idx)
        clear idx
        continue
        
    elseif length(idx)==length(dat_grouped{y2,1}{i,1})
        for j=1:length(dat_grouped)
            dat_grouped_edge_removed{j,1}{w,1}=dat_grouped{j,1}{i,1};
        end
        w=w+1;
    else
        ind=find(diff(idx)>1);
        ind=[0;ind;length(idx)];
        for b=1:length(ind)-1
            for j=1:length(dat_grouped)
                dat_grouped_edge_removed{j,1}{w,1}=dat_grouped{j,1}{i,1}(idx(ind(b)+1):idx(ind(b+1)));
            end
            w=w+1;
        end    
    end
    
    
    fields_grouped_edge_removed=fields_grouped;
end
end

