function [dat_edge_removed,fields_edge_removed]=eliminate_edge_JB(dat,fields)

%dat_grouped_edge_removed can be used to calculate PI and tt


x2=find(matches(fields,'x'));
y2=find(matches(fields,'y'));
AN=find(matches(fields,'AN'));
xs=find(matches(fields,'xspine'));
ys=find(matches(fields,'yspine'));

w=1;

for i=1:length(dat{x2,1})
    idx=find(dat{y2,1}{i,1}>25&dat{y2,1}{i,1}<250&dat{x2,1}{i,1}>5&dat{x2,1}{i,1}<222);
    if isempty(idx)
        clear idx
        continue
        
    elseif length(idx)==length(dat{y2,1}{i,1})
        for j=1:length(dat)
            dat_edge_removed{j,1}{w,1}=dat{j,1}{i,1};
        end
        w=w+1;
    else
        ind=find(diff(idx)>1);
        ind=[0;ind;length(idx)];
        for b=1:length(ind)-1
            for j=1:length(dat)
                if j==AN
                    dat_edge_removed{j,1}{w,1}=dat{j,1}{i,1};
                elseif j==xs|j==ys
                    dat_edge_removed{j,1}{w,1}=dat{j,1}{i,1}(idx(ind(b)+1):idx(ind(b+1)),:);
                else 
                   dat_edge_removed{j,1}{w,1}=dat{j,1}{i,1}(idx(ind(b)+1):idx(ind(b+1))); 
                end
               
            end
             w=w+1;
        end
        fields_edge_removed=fields;
    end
end
end 



