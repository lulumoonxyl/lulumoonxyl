function dat_JB=change_in_theta_for_turning_event(dat_JB,dat_JAABA,dis,vector)
%this function will calculate the x where each turning event occurs and the
%change in heading direction of each turning event
%it will also calculate the pre-turning and  post_turning heading direction

for i=1:length(dat_JAABA.t0_idx)
    idx0=dat_JAABA.t0_idx{i,1};
    idx1=dat_JAABA.t1_idx{i,1};
    
    if isempty(idx0)||isempty(idx1)
        dat_JB.reorientation_angle{i,1}=[];dat_JB.turning_x{i,1}=[];
        dat_JB.pre_turning_angle{i,1}=[];dat_JB.post_turning_angle{i,1}=[];
        
        clear idx0 idx1
        continue;
    end
    
    if idx1(end)==length(dat_JB.x{i,1})
        idx0(end)=[];idx1(end)=[];
    end
    
    if isempty(idx0)||isempty(idx1)
        dat_JB.reorientation_angle{i,1}=[];dat_JB.turning_x{i,1}=[];
        dat_JB.pre_turning_angle{i,1}=[];dat_JB.post_turning_angle{i,1}=[];
        clear idx0 idx1
        continue;
    end
    
    if idx0(1)==1
        idx0(1)=[]; idx1(1)=[];
    end
    
    if isempty(idx0)||isempty(idx1)
        dat_JB.reorientation_angle{i,1}=[];dat_JB.turning_x{i,1}=[];
        dat_JB.pre_turning_angle{i,1}=[];dat_JB.post_turning_angle{i,1}=[];
        
        clear idx0 idx1
        continue;
    end
    %get the x position where turning occurs
    x0=dat_JB.x{i,1}(idx0);
    x1=dat_JB.x{i,1}(idx1);
    dat_JB.turning_x{i,1}=x0+(x1-x0)./2;
    
    for j=1:length(idx0)
        %calculate the pre-turning heading direction based on
        %interval-->how many timepoints
        pre_deg=get_pre_or_post_turning_direction(vector,dis,i,j,idx0,dat_JB,'y','n');
        %calculate the post-turning heading direction
        post_deg=get_pre_or_post_turning_direction(vector,dis,i,j,idx1,dat_JB,'n','y');
        dat_JB.pre_turning_angle{i,1}(j,1)=pre_deg;
        dat_JB.post_turning_angle{i,1}(j,1)=post_deg;
        
        turn_deg=dat_JB.orientation{i,1}(idx0(j):idx1(j));
        deg=[pre_deg;post_deg];
        [q,N]=quadrant_check(deg);
        [q_turn,N_turn]=quadrant_check(turn_deg);
        %get the absolute value of reorientation angle, somehow we need to clarify
        %whether it will be clockwise and counterclockwise, also angle in the
        %quadrant [2;3] needs more attention.
        
        if (~isequal(q,[2;3]))&&(~isequal(q,[3;2]))
            if length(N_turn)==4
                turn_theta=360-abs(post_deg-pre_deg);
            else
                turn_theta=abs(post_deg-pre_deg);
            end
        else
            if length(N_turn)==4
                turn_theta=abs(post_deg)+abs(pre_deg);
            else
                turn_theta=360-abs(post_deg)-abs(pre_deg);
            end
        end
        dat_JB.reorientation_angle{i,1}(j,1)=turn_theta;
        
    end
end
clearvars -except dat*
end

