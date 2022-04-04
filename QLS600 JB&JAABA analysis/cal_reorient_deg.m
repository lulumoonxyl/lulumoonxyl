function reorient_deg=cal_reorient_deg(pre_deg,post_deg)
%calculate the absolute value of reorientation angle using pre/post-turning
%degree 
for i=1:length(pre_deg)
    if isempty(pre_deg{i,1})||isempty(post_deg{i,1})
        reorient_deg{i,1}=[];
        continue
    end 
    %if the pre/post angles are in second or third quadrant ([90 180] and [-180 190])
    %we the reorientation angle will be 360-abs(pre)-abs(post)
    %the rest will just be post-pre
    reorient_deg{i,1}=abs(post_deg{i,1}-pre_deg{i,1});
    [q_pre,N_pre]=quadrant_check(pre_deg{i,1});
    [q_post,N_post]=quadrant_check(post_deg{i,1});
    
    idx=find((q_pre==2&q_post==3)|(q_pre==3&q_post==2));
    reorient_deg{i,1}(idx)=360-reorient_deg{i,1}(idx);
end 
end