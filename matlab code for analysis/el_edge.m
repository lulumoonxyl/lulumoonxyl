function dat_edge=el_edge(dat,edge)
%this function will eliminate the x and y value that are in the edge also
%elimiante other datapoints have the same ind
name=dat.Properties.VariableNames;
ind_t=find(matches(name,'t0s')|matches(name,'t1s')|matches(name,"AN"));
name(ind_t)=[];
ind_spine=find(contains(name,"spine"));
w=1;

for i=1:length(dat.AN)
    %find the ind where x and y value are in the region of interest instead
    %of in the edge
    idx=find(dat.x{i,1}>edge(1,1)&dat.x{i,1}<edge(1,2)&dat.y{i,1}>edge(2,1)&dat.y{i,1}<edge(2,2));

    if isempty(idx)||length(idx)==1
        clear idx
        continue
    else
        % if larvae going back and forth in the edge, we need to divide it
        % into different trajectories
        ind=find(diff(idx)>1);
        ind=[0;ind;length(idx)];
        for b=1:length(ind)-1
            
            dat_edge.AN(w,1)=dat.AN(i,1);
            for k=1:length(name)
                if ismember(k,ind_spine)
                    dat_edge.(name{k}){w,1}=dat.(name{k}){i,1}(idx(ind(b)+1):idx(ind(b+1)),:);
                else
                    dat_edge.(name{k}){w,1}=dat.(name{k}){i,1}(idx(ind(b)+1):idx(ind(b+1)));
                end
            end
            %eliminate the edge for the t0 t1 idx
            ind1=idx(ind(b)+1);
            ind2=idx(ind(b+1));
            ind_array=(ind1:1:ind2)';

            dat0=dat.t0s{i,1};
            dat1=dat.t1s{i,1};

            idx0=find(ismember(dat0,dat.et{i,1}(ind_array)));
            idx1=find(ismember(dat1,dat.et{i,1}(ind_array)));

            if isempty(idx0)||isempty(idx1)
                %if no data falls into the region of your interest
                dat_edge.t0s{w,1}={};
                dat_edge.t1s{w,1}={};
            else
                
                idx2=intersect(idx0,idx1);
                %get the t0s and t1s
                
                 dat_edge.t0s{w,1}=dat.t0s{i,1}(idx2);
                dat_edge.t1s{w,1}=dat.t1s{i,1}(idx2);

            end

            w=w+1;
        end


    end

end

dat_edge=struct2table(dat_edge);
end