function [dat_array,dat_grouped,fields_grouped]=calculate_navigational_index(dat_grouped,fields_grouped,ii,t,dat_array)
%may also calculate the PI based on different periods
%t will usually be 300s or []

sp_pos=find(matches(fields_grouped,'sp'));
et_pos=find(matches(fields_grouped,'et'));
x_pos=find(matches(fields_grouped,'x'));
l=length(fields_grouped);
fields_grouped{l+1}='PI';
fields_grouped{l+2}='PI mean';


for i=1:length(dat_grouped{sp_pos,1})
    x_diff=diff(dat_grouped{x_pos,1}{i,1});
    t_diff=diff(dat_grouped{et_pos,1}{i,1});
    speed=dat_grouped{sp_pos,1}{i,1}(1:end-1);
    PI=(x_diff./t_diff)./speed;
    
    idx=find(speed==0);
    PI(idx)=[];
    dat_grouped{l+1}{i,1}=PI;
    dat_grouped{l+2}(i,1)=mean(PI);
    clear PI idx x_diff t_diff speed
end
PI_mean(ii,1)=mean(dat_grouped{l+2,1});
PI_SEM(ii,1)=std(dat_grouped{l+2,1})/length(dat_grouped{l+2,1});
dat_array{ii,1}=dat_grouped{l+2};

if isempty(t)~=1
    time=0:t:900;
    fields_grouped{l+3}='PI t';
    fields_grouped{l+4}='PI t mean';
    for i=1:length(dat_grouped{sp_pos,1})
        for j=1:length(time)-1
            idx=find(dat_grouped{et_pos}{i,1}>=time(j)&dat_grouped{et_pos}{i,1}<time(j+1));
            if isempty(idx)
                clear idx
                dat_grouped{l+3}{i,j}=nan;
                dat_grouped{l+4}(i,j)=nan;
                continue
            end 
            x_diff=diff(dat_grouped{x_pos,1}{i,1}(idx));
            t_diff=diff(dat_grouped{et_pos,1}{i,1}(idx));
            speed=dat_grouped{sp_pos,1}{i,1}(idx);
            speed(end)=[];
            PI_t=(x_diff./t_diff)./speed;
            
            idx1=find(speed==0);
            PI_t(idx1)=[];
            dat_grouped{l+3}{i,j}=PI_t;
            dat_grouped{l+4}(i,j)=mean(PI_t);
            clear PI_t idx x_diff t_diff speed idx1
        end
    end
    
    
end
end