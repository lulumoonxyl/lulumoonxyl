function [t,p1,ttl_t,tsum]=bin_data_sum_t(data,et,series,varargin)
%this function will bin the data based on the series, and it will sum all
%the time that the larvae spend after binning with series
%varargins will be 1)t0_idx: the index where turning starts 2)t1_idx:
%index where turning events stop-->using these input, we will eliminate the
%turning event when compute the time; if t0_idx and t1_idx are empty, we
%will not eliminate the time when turning event occurs
l=length(series)-1;
ld=length(data);

t0_idx={};t1_idx={};
for i=1:2:length(varargin)

    if strcmp(varargin{i},"t0_idx")
        t0_idx=varargin{i+1};
    elseif strcmp(varargin{i},"t1_idx")
        t1_idx=varargin{i+1};

    end
end
%create empty arry for the t for each larvae and each bin

t=NaN(ld,l);
p=NaN(ld,l);

ttl_t=zeros(ld,1);

for j=1:ld
    %find all the index where turning events occur
    idx=[];
    ttl_t(j,1)=et{j,1}(end)-et{j,1}(1);

    if ~isempty(t0_idx)&&~isempty(t1_idx)
        if ~isempty(t0_idx{j,1})
            for m=1:length(t0_idx{j,1})
                idx=vertcat(idx,(t0_idx{j,1}(m):t1_idx{j,1}(m))');
            end
        end

        %find the total time that turning occurs, subtrct it from the tsum of
        %this larvae
        if ~isempty(idx)
            turn_t=sum(et{j,1}(idx)-et{j,1}(idx-1));
            ttl_t(j,1)=ttl_t(j,1)-turn_t;
            clear turn_t 
        end
        
    end
    %bin the data according to the series
    for i=1:l
        if i==l
            ind1=find(data{j,1}>=series(i)&data{j,1}<=series(i+1));
        else

            ind1=find(data{j,1}>=series(i)&data{j,1}<series(i+1));
        end

        if ~isempty(ind1)
            if ind1(1)==1
                ind1(1)=[];
            end

            if ~isempty(ind1)
                ind2=find(ismember(ind1,idx));
                if ~isempty(ind2)
                    ind1(ind2)=[];
                end
                if ~isempty(ind1)
                    t(j,i)=sum(et{j,1}(ind1)-et{j,1}(ind1-1));
                    clear ind*
                end
            end
        end
    end

end

%sum all the time in the t matrix to confirm we do not count something
%twice

tsum=sum(t,2,"omitnan");% tracking time of each animal, tsum has a size of (ld,1)
%% get the p1 by dividing the t by all tsum 
p1=sum(t,1,'omitnan')./sum(t,'all','omitnan');
p1=p1';
%% since each larvae has different tracking time, we need to weight them to
%%get the mean of probability and std
% p=t./tsum;
% w=tsum./sum(tsum);
% pw=p.*w;
% p_mean(:,1)=sum(pw,'omitnan')';
% for i=1:width(p)
%     idx=find(isnan(p(:,i)));
%     len(i,1)=length(p(:,1))-length(idx);
% end
% p_mean(:,2)=(std(p,w,'omitnan')')./sqrt(len);

end

