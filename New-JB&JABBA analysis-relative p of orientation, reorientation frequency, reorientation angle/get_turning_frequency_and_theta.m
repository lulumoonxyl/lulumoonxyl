function [result,result_split]=get_turning_frequency_and_theta(dat_JB,deg,split_or_not,bins,large_or_small)
%this function will get the turning frequency (turning event/tracking
%time for each larvae
%based on the pre-turning heading direction (4 groups: [-45 45],[-180
%-135]&[135 180],[45,135],[-135,-45]; and turning frequency for all
%positions if split_or_not is not empty
%a large turning event is delta theta>deg

result.x=[0:bins:220]';
if isequal(split_or_not,'y')
    result_split.x=result.x;
end
%% calculate the turning frequency and large turning event/all turning ratio ,and delta theta even
%for each larva
for i=1:length(result.x)-1
    w=1;u=1;q=[1 1 1 1];qq=[1 1 1 1];
    for j=1:length(dat_JB.x)
        %find how many datapoints this larva has in this xbin
        idx=find(dat_JB.x{j,1}>=result.x(i)&dat_JB.x{j,1}<result.x(i+1));
        
        if isempty(idx)
            continue
        else
            if idx(1)==1
                idx(1)=[];
            end
            %find the number of turning events that this larva has in this
            %xbin
            if isempty(idx)
            continue
            end 
            str='';
            idx2=find(dat_JB.turning_x{j,1}>=result.x(i)&dat_JB.turning_x{j,1}<result.x(i+1));
            [result,w,u]=calculate_tf_reorientation_angle(dat_JB,result,idx,idx2,i,j,w,u,deg,str,large_or_small);
            
        end
        %% get the tf,large tf and delta theta for different heading directions
        if isequal(split_or_not,'y')
            h_dir=dat_JB.orientation{j,1}(idx);
            dir=[-180 -135 -45 45 135 180];
            for d=1:length(dir)-2
                if d==1
                    idx_quadrant{d,1}=find((h_dir>=dir(d)&h_dir<dir(d+1))|(h_dir<dir(end)&h_dir>=dir(end-1)));
                else
                    idx_quadrant{d,1}=find(h_dir<dir(d+1)&h_dir>=dir(d));
                end
            end
            
            %find the number of turn happens in this xbin with heading
            %direction in the specific quadrant
            pre_turn_deg=dat_JB.pre_turning_angle{j,1}(idx2);
            for d=1:length(dir)-2
                if d==1
                    idx_pre{d,1}=find((pre_turn_deg>=dir(d)&pre_turn_deg<dir(d+1))|(pre_turn_deg<dir(end)&pre_turn_deg>=dir(end-1)));
                else
                    idx_pre{d,1}=find(pre_turn_deg<dir(d+1)&pre_turn_deg>=dir(d));
                end
            end
            
            %%
            
            str_quadrant={'_q1','_q2','_q3','_q4'};
            for d=1:length(idx_quadrant)
                if ~isempty(idx_quadrant{d,1})
                    [result_split,q(d),qq(d)]=calculate_tf_reorientation_angle(dat_JB,result_split,...
                        idx(idx_quadrant{d,1}),idx2(idx_pre{d,1}),i,j,q(d),qq(d),deg,str_quadrant(d),large_or_small);
                else
                    continue
                end
            end
            
        end
    end
    
end

end


