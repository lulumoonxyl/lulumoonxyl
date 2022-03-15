function [dat_grouped,PI_mean,PI_sem,PI_t_mean,PI_t_sem]=calculate_navigational_index_velocity_vs_speed(dat_grouped,ii,t,PI_mean,PI_sem,PI_t_mean,PI_t_sem)
%PI = v_x/speed in each timepoint for each trajectory
%the input t will allow you to separate the recording into different
%duration

fn=fieldnames(dat_grouped);

for i=1:length(dat_grouped.sp)
    x_diff=diff(dat_grouped.x{i,1});
    t_diff=diff(dat_grouped.et{i,1});
    speed=dat_grouped.sp{i,1}(1:end-1);
    PI=(x_diff./t_diff)./speed;
    
    idx=find(speed==0);
    idx0=find(t_diff==0);
    idx3=[idx;idx0];
    PI(idx3)=[];
    dat_grouped.PI{i,1}=PI;
    dat_grouped.PI_mean(i,1)=mean(PI);
    dat_grouped.PI_sample_size=length(dat_grouped.PI_mean);
    clear PI idx x_diff t_diff speed idx0 idx3
end
%calculate the overall mean of PI for each condition
PI_mean(ii,1)=mean(dat_grouped.PI_mean);
PI_sem(ii,1)=std(dat_grouped.PI_mean)/sqrt(length(dat_grouped.PI_mean));

if isempty(t)~=1
    time=0:t:900;
    for j=1:length(time)-1
        w=1;
        for i=1:length(dat_grouped.sp)
            idx=find(dat_grouped.et{i,1}>=time(j)&dat_grouped.et{i,1}<time(j+1));
            if isempty(idx)
                clear idx
                continue
            end
            x_diff=diff(dat_grouped.x{i,1}(idx));
            t_diff=diff(dat_grouped.et{i,1}(idx));
            speed=dat_grouped.sp{i,1}(idx);
            speed(end)=[];
            PI_t=(x_diff./t_diff)./speed;
            
            idx1=find(speed==0);
            idx2=find(t_diff==0);
            idx=[idx1;idx2];
            PI_t(idx)=[];
            dat_grouped.PI_t{w,j}=PI_t;
            dat_grouped.PI_mean_t(w,j)=mean(PI_t);
            w=w+1;
            clear PI_t idx x_diff t_diff speed idx1 idx2 idx
        end
        idx=find(cellfun(@isempty,dat_grouped.PI_t(:,j)));
        PI=dat_grouped.PI_mean_t(:,j);
        PI(idx)=[];
        PI_t_mean(ii,j)=mean(PI);
        PI_t_sem(ii,j)=std(PI)/sqrt(length(PI));
         dat_grouped.PI_sample_size_t(j,1)=length(PI);
        clear idx
    end
end


end