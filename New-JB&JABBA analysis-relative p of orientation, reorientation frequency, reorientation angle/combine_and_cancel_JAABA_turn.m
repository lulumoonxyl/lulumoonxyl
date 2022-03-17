function dat_JAABA=combine_and_cancel_JAABA_turn(dat_JAABA)
for i=1:length(dat_JAABA.t1_idx)
    if ~isempty(dat_JAABA.t0_idx{i,1})&&~isempty(dat_JAABA.t1_idx{i,1})
        if length(dat_JAABA.t0_idx{i,1})>1
            turn_interval=dat_JAABA.t0_idx{i,1}(2:end)-dat_JAABA.t1_idx{i,1}(1:end-1);
            %if the difference between two turning events is less than 5
            %sections, we combine them 
            ind=find(turn_interval<=5);
            dat_JAABA.t0_idx{i,1}(ind+1)=[];
            dat_JAABA.t1_idx{i,1}(ind)=[];
            dat_JAABA.t0s{i,1}(ind+1)=[];
            dat_JAABA.t1s{i,1}(ind)=[];
        end
        
        turn_diff=dat_JAABA.t1_idx{i,1}-dat_JAABA.t0_idx{i,1};
        %after combination if the turning event takes less than 4
        %timepoints to turn, delete it 
        idx=find(turn_diff<=4);
        dat_JAABA.t1_idx{i,1}(idx)=[];
        dat_JAABA.t0_idx{i,1}(idx)=[];
        dat_JAABA.t1s{i,1}(idx)=[];
        dat_JAABA.t0s{i,1}(idx)=[];
    end
end
end