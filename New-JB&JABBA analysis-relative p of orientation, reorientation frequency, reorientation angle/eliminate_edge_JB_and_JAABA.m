function [dat_JB_edge,dat_JAABA_edge]=eliminate_edge_JB_and_JAABA(dat_JB_ns,edge,dat_JAABA_ns)
w=1;
fn=fieldnames(dat_JB_ns);
fn_JAABA=fieldnames(dat_JAABA_ns);
for i=1:length(dat_JB_ns.x)
    idx=find(dat_JB_ns.y{i,1}>edge(2,1)&dat_JB_ns.y{i,1}<edge(2,2)&dat_JB_ns.x{i,1}>edge(1,1)&dat_JB_ns.x{i,1}<edge(1,2));
    if isempty(idx)
        clear idx
        continue

    elseif length(idx)==length(dat_JB_ns.y{i,1})
        for j=1:length(fn)
            dat_JB_edge.(fn{j}){w,1}=dat_JB_ns.(fn{j}){i,1};
        end

        for j=1:length(fn_JAABA)
            dat_JAABA_edge.(fn_JAABA{j}){w,1}=dat_JAABA_ns.(fn_JAABA{j}){i,1};
        end
        w=w+1;

    else
        %what if there are multiple sessions that fall into the region of
        %interest
        ind=find(diff(idx)>1);
        ind=[0;ind;length(idx)];
        for b=1:length(ind)-1
            for j=1:length(fn)
                if isequal(fn{j},'AN')
                    dat_JB_edge.(fn{j}){w,1}=dat_JB_ns.(fn{j}){i,1};
                elseif isequal(fn{j},'xspine')|isequal(fn{j},'yspine')
                    dat_JB_edge.(fn{j}){w,1}=dat_JB_ns.(fn{j}){i,1}(idx(ind(b)+1):idx(ind(b+1)),:);
                else
                    dat_JB_edge.(fn{j}){w,1}=dat_JB_ns.(fn{j}){i,1}(idx(ind(b)+1):idx(ind(b+1)));
                end

            end

            % eliminate edge for JAABA data
            dat_JAABA_edge.AN{w,1}=dat_JAABA_ns.AN{i,1};

            ind1=idx(ind(b)+1);
            ind2=idx(ind(b+1));
            ind_array=(ind1:1:ind2)';

            dat0=dat_JAABA_ns.t0_idx{i,1};
            dat1=dat_JAABA_ns.t1_idx{i,1};

            idx0=find(ismember(dat0,ind_array));
            idx1=find(ismember(dat1,ind_array));

            if isempty(idx0)||isempty(idx1)
                %if no data falls into the region of your interest
                dat_JAABA_edge.t0_idx{w,1}={};
                dat_JAABA_edge.t1_idx{w,1}={};
                dat_JAABA_edge.t0s{w,1}={};
                dat_JAABA_edge.t1s{w,1}={};

            else
                idx2=intersect(idx0,idx1);
                %get the correct index for the t0 and t1
                t0_edge=dat0(idx2)-ind1+1;
                t1_edge=dat1(idx2)-ind1+1;
                %assign them to the new structure
                dat_JAABA_edge.t0_idx{w,1}=t0_edge;
                dat_JAABA_edge.t0_idx{w,1}=t1_edge;
                dat_JAABA_edge.t0s{w,1}=dat_JAABA_ns.t0s{i,1}(idx2);
                dat_JAABA_edge.t1s{w,1}=dat_JAABA_ns.t1s{i,1}(idx2);
                
            end
            
            w=w+1;
        end
    end
    
end

end 





