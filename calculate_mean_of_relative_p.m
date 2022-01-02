function [t_p_mean,t_p_sem,t_abs_p]=calculate_mean_of_relative_p(dat_JB,fields_JB,t_p,t_abs)
%this function will calculate the relative probability of theta ith 2
%methods
%1)calculate the p for each larva and get the mean and sem
%2)calculate the sum of all tracking time and get the relative p by
%dividing it by the tsum 
tsum_pos=find(matches(fields_JB,'tsum'));
t_p_mean=mean(t_p,2);
t_p_sem=std(t_p,[],2)./length(t_p);

t_abs_sum=sum(t_abs,2);
tsum=sum(dat_JB{tsum_pos,1});

t_abs_p=t_abs_sum./tsum;
end 