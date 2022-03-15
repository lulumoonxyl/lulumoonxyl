function [dat_ns,dat_grouped_ns]=remove_the_first_ns(dat,dat_grouped,ns)
%% this function will remove the first n seconds for the dat
    fn=fieldnames(dat);
    idx=find(dat.et>=ns);
    for i=1:length(fn)
        dat_ns.(fn{i})=dat.(fn{i})(idx);
    end
    clear idx
%% remove the first n seconds for the dat_grouped
    fn_grouped=fieldnames(dat_grouped);
    w=1;
    for i=1:length(dat_grouped.et)
        idx=find(dat_grouped.et{i,1}>=ns);
        if isempty(idx)
            clear idx
            continue
        else
            for j=1:length(fn_grouped)
                dat_grouped_ns.(fn_grouped{j}){w,1}=dat_grouped.(fn_grouped{j}){i,1}(idx);
            end
            w=w+1;
        end
    end

end