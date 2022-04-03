function [dat_grouped,PI_x_mean,PI_x_sem]=calculate_PI_vx_sp_xpos(dat_grouped,xpos)
%dat_grouped is the data from choreography
%this function will get the PI based on the x position
%you can place this function after the calculate_PI_vx_sp function
%xpos will just be 0:xbins:220; if xbins=20, xpos=0:20:220
l=length(xpos)-1;
dat_grouped.PI_x_larva=cell(l,1);
for i=1:l
    w=1;
    for j=1:length(dat_grouped.x)
        if i==l
            %find how many points that the larva spend in this xbin
            idx=find(dat_grouped.x{j,1}>=xpos(i)&dat_grouped.x{j,1}<=xpos(i+1));
        else
            idx=find(dat_grouped.x{j,1}>=xpos(i)&dat_grouped.x{j,1}<xpos(i+1));
        end

        if isempty(idx)||length(idx)==1
            clear idx
            continue
        end
        x_diff=diff(dat_grouped.x{j,1}(idx));
        y_diff=diff(dat_grouped.y{j,1}(idx));
        t_diff=diff(dat_grouped.et{j,1}(idx));

        idx0=find(t_diff==0);
        idx1=find(x_diff==0&y_diff==0);
        idx2=[idx0;idx1];
        x_diff(idx2)=[];
        y_diff(idx2)=[];
        t_diff(idx2)=[];

        vx=x_diff./t_diff;
        vy=y_diff./t_diff;
        PI=vx./sqrt(vx.^2+vy.^2);% get the PI for each timepoint in the xbin, sp=sqrt(vx.^2+vy.^2)

        PI_larva=mean(PI);
        dat_grouped.PI_x_larva{i,1}(w,1)=PI_larva;
        w=w+1;
        clear PI *diff  idx*
    end
    if ~isempty(dat_grouped.PI_x_larva{i,1})
        PI_x_mean(1,i)=mean(dat_grouped.PI_x_larva{i,1});
        PI_x_sem(1,i)=std(dat_grouped.PI_x_larva{i,1})/sqrt(length(dat_grouped.PI_x_larva{i,1}));
        dat_grouped.PI_sample_size_x(1,i)=length(dat_grouped.PI_x_larva{i,1});
    else
        PI_x_mean(1,i)=0;  PI_x_sem(1,i)=0; dat_grouped.PI_sample_size_x(1,i)=0;
    end
    clear idx
end

end
