function dat_JAABA=align_JB_JAABA_data(dat_JB,dat_JAABA)
%this function will eliminate some data from the JAABA array as it contains
%slightly more cells than JB array
fn=fieldnames(dat_JAABA);
w=1;
for i=1:length(dat_JB.AN)
    while dat_JB.AN{i,1}~=dat_JAABA.AN{w,1}
        for j=1:length(fn)
            dat_JAABA.(fn{j})(w)=[];
        end
    end
    w=w+1;
end

if length(dat_JB.AN)<length(dat_JAABA.AN)
    idx=length(dat_JAABA.AN):-1:length(dat_JB.AN)+1;
    for j=1:length(fn)
        dat_JAABA.(fn{j})(idx)=[];
    end
end

% this function will also find the corresponding index for the t0s and t1s
for i=1:length(dat_JB.et)
    dat0=dat_JAABA.t0s{i,1};
    dat1=dat_JAABA.t1s{i,1};
    time=dat_JB.et{i,1};
    
    ind0=find(~ismember(dat0,time))';
    ind1=find(~ismember(dat1,time))';
    ind=[ind0,ind1];
    dat0(ind)=[];
    dat1(ind)=[];
    
    idx0=find(ismember(time,dat0));
    idx1=find(ismember(time,dat1));
    dat_JAABA.t0_idx{i,1}=idx0;
    dat_JAABA.t1_idx{i,1}=idx1;
end

clearvars -except dat_JAABA 
end

