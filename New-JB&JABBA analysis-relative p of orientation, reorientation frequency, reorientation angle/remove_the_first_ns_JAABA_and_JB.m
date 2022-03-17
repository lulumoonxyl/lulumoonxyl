function [dat_JB_ns,dat_JAABA_ns]=remove_the_first_ns_JAABA_and_JB(dat_JAABA,ns,dat_JB)
%this function will remove data in the first n s; and this function should
%be placed after we align the JB and JABBA data arrays
fn=fieldnames(dat_JB);
fn_JAABA=fieldnames(dat_JAABA);
w=1;
for i=1:length(dat_JB.et)
    idx=find(dat_JB.et{i,1}>=ns);
    if isempty(idx)
        clear idx
        continue
    else
        for j=1:length(fn)
            if isequal((fn{j}),'AN')
                dat_JB_ns.AN{w,1}=dat_JB.AN{i,1};
            elseif isequal((fn{j}),'xspine')||isequal((fn{j}),'yspine')
                dat_JB_ns.(fn{j}){w,1}=dat_JB.(fn{j}){i,1}(idx,:);
            else
                dat_JB_ns.(fn{j}){w,1}=dat_JB.(fn{j}){i,1}(idx);
            end
        end
        %remove the first ns from JAABA data

        for j=1:length(fn_JAABA)
            dat_JAABA_ns.(fn_JAABA{j}){w,1}=dat_JAABA.(fn_JAABA{j}){i,1};
        end

        if length(idx)~=length(dat_JB.et{i,1})
            ind=idx(1);
            idx1=find(dat_JAABA_ns.t0_idx{w,1}<ind|dat_JAABA_ns.t1_idx{w,1}<ind);
            dat_JAABA_ns.t0_idx{w,1}(idx1)=[];
            dat_JAABA_ns.t1_idx{w,1}(idx1)=[];
            dat_JAABA_ns.t0_idx{w,1}=dat_JAABA_ns.t0_idx{w,1}-ind+1;
            dat_JAABA_ns.t1_idx{w,1}=dat_JAABA_ns.t1_idx{w,1}-ind+1;
            clear ind idx1

            idx2=find(dat_JAABA_ns.t0s{w,1}<ns|dat_JAABA_ns.t1s{w,1}<ns);
            dat_JAABA_ns.t0s{w,1}(idx2)=[];
            dat_JAABA_ns.t1s{w,1}(idx2)=[];
        end
        w=w+1;

    end
end

clearvars -except dat*
end