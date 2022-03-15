function [dat_edge,dat_grouped_edge]=eliminate_edge(dat, dat_grouped,edge_x_y)
%this will remove the edge of all the data
%edge_x_y=[xmin,xmax;ymin,ymax]-->[5,222;25,250]
fn=fieldnames(dat);
idx=find(dat.y>edge_x_y(2,1)&dat.y<edge_x_y(2,2)&dat.x>edge_x_y(1,1)&dat.x<edge_x_y(1,2));
for i=1:length(fn)
    dat_edge.(fn{i})=dat.(fn{i})(idx);
end
clear idx

fn_grouped=fieldnames(dat_grouped);
w=1;
dat_grouped_edge={};

for i=1:length(dat_grouped.x)
    idx=find(dat_grouped.y{i,1}>edge_x_y(2,1)&dat_grouped.y{i,1}<edge_x_y(2,2)&dat_grouped.x{i,1}>edge_x_y(1,1)&dat_grouped.x{i,1}<edge_x_y(1,2));
    if isempty(idx)
        clear idx
        continue

    elseif length(idx)==length(dat_grouped.y{i,1})
        for j=1:length(fn_grouped)
            dat_grouped_edge.(fn_grouped{j}){w,1}=dat_grouped.(fn_grouped{j}){i,1};
        end
        w=w+1;
    else
        ind=find(diff(idx)>1);
        ind=[0;ind;length(idx)];
        for b=1:length(ind)-1
            for j=1:length(fn_grouped)
                dat_grouped_edge.(fn_grouped{j}){w,1}=dat_grouped.(fn_grouped{j}){i,1}(idx(ind(b)+1):idx(ind(b+1)));
            end
            w=w+1;
        end
    end

end
end

