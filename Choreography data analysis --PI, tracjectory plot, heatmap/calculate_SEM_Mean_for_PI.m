function [PI_mean,PI_t_mean,PI_sem,PI_t_sem]=calculate_SEM_Mean_for_PI(dat_grouped,fields_grouped,ii,PI_mean,PI_t_mean,PI_sem,PI_t_sem)
%This function will just organize the data and get the mean and sem for the
%navigation index which will be plotted later

%this function will also save the PI for each larvae which will be used for
%stat test ANOVA 
pos=find(matches(fields_grouped,'PI mean'));
%find out the nan pos and eliminate those

Total_PI=dat_grouped{pos,1};
idx=find(isnan(Total_PI));
if isempty(idx)~=1
Total_PI(idx)=[];
end 
clear idx

PI_mean(ii,1)=mean(Total_PI);
PI_sem(ii,1)=std(Total_PI)/sqrt(length(Total_PI));


pos1=find(matches(fields_grouped,'PI t mean'));
if isempty(pos1)
    PI_t_mean=[];
else 
    [row,col]=size(dat_grouped{pos1,1});
    for i=1:col
        PI=dat_grouped{pos1,1}(:,i);
        idx=find(isnan(PI));
        PI(idx)=[];
        PI_t_mean(ii,i)=mean(PI);
        PI_t_sem(ii,i)=std(PI)/sqrt(length(PI));
        clear idx 
    end 
end

end 