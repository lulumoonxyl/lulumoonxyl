function [t,p,p_bin,t_series2,tsum,ttl_t]=bin_data_sum_t2(data1,data2,series1,series2,et,varargin)
%this function will bin the data based on two series and then count all the
%time that it spend in each x and hd bins etc
%data1 and series1 is the orientation
%data2 and series2 is the x or y or t
%ttl_t is the total tracking time of each animal =tsum
%t_series2 is the tracking time for each object spend in each bin in
%series2

l=length(series1)-1;
ld=length(data1);
l2=length(series2)-1;
t0_idx={};t1_idx={};
for i=1:2:length(varargin)

    if strcmp(varargin{i},"t0_idx")
        t0_idx=varargin{i+1};
    elseif strcmp(varargin{i},"t1_idx")
        t1_idx=varargin{i+1};
    end
end
%create empty arry for the t for each larvae and each bin
t=NaN(l,l2);
t_series2=zeros(ld,l2);
% p=NaN(l,l2);
% p_obj=cell(l,l2);
% w=ones(l,l2);
% 
% p_mean=zeros(l,l2);
% p_sem=zeros(l,l2);

for j=1:ld
    %find all the index where turning events occur
    ttl_t(j,1)=et{j,1}(end)-et{j,1}(1);
    idx=[];
    if ~isempty(t0_idx)&&~isempty(t1_idx)
        if ~isempty(t0_idx{j,1})
            for m=1:length(t0_idx{j,1})
                idx=vertcat(idx,(t0_idx{j,1}(m):t1_idx{j,1}(m))');
            end

            %find the total time that turning occurs, substract it from the tsum of
            %this larvae
            turn_t=sum(et{j,1}(idx)-et{j,1}(idx-1));
            ttl_t(j,1)=et{j,1}(end)-et{j,1}(1)-turn_t;
            clear turn_t
        end
    end
    t_an=NaN(l,l2);

    %bin the data1 according to the series1, and data2 according to series2
    for i=1:l

        for z=1:l2
            if z==l2 && i==l
                ind1=find(data1{j,1}>=series1(i)&data1{j,1}<=series1(i+1)&data2{j,1}>=series2(z)&data2{j,1}<=series2(z+1));
            elseif z~=l2 && i~=l
                ind1=find(data1{j,1}>=series1(i)&data1{j,1}<series1(i+1)&data2{j,1}>=series2(z)&data2{j,1}<series2(z+1));
            elseif z==l2
                ind1=find(data1{j,1}>=series1(i)&data1{j,1}<series1(i+1)&data2{j,1}>=series2(z)&data2{j,1}<=series2(z+1));
            elseif i==l
                ind1=find(data1{j,1}>=series1(i)&data1{j,1}<=series1(i+1)&data2{j,1}>=series2(z)&data2{j,1}<series2(z+1));
            end

            if ~isempty(ind1)
                if ind1(1)==1
                    ind1(1)=[];
                end

                if ~isempty(ind1)
                    ind2=find(ismember(ind1,idx));
                    ind1(ind2)=[];
                    if ~isempty(ind1)
                        t_an(i,z)=sum(et{j,1}(ind1)-et{j,1}(ind1-1));

                        % this will generate a matrix of the time each object spend in each bin combination
                        clear ind*
                    end
                end
            end
        end
    end

    %% compute the probability of orientation for each object if it is not nan
%     t_an_sum=sum(t_an,1,"omitnan");
%     for i=1:l
%         for z=1:l2
%             if ~isnan(t_an(i,z))
%                 p_obj{i,z}(w(i,z),1)=t_an(i,z)/t_an_sum(1,z);
%                 w(i,z)=w(i,z)+1;
%             end
%         end
%     end
    %sum all the time of one object
    tsum(j,1)=sum(t_an,'all',"omitnan");
    %sum the time of one object spent in the binning of series2
    t_series2(j,:)=sum(t_an,1,'omitnan');
    %add the time of this object to previos object:variable t
    T=cat(3,t_an,t);
    t=sum(T,3,'omitnan');
    clear T t_an t_an_sum

end
%% Compute the p of each animal to series 2
% for i=1:l
%     for z=1:l2
%         p_mean(i,z)=mean(p_obj{i,z},'omitnan');
%         p_sem(i,z)=std(p_obj{i,z},'omitnan')/length(p_obj{i,z});
%     end
% end
%% sum all the time in the t matrix to confirm we do not count something
%twice

p=t./sum(tsum,'omitnan');%divide by total time
p_bin=t./sum(t,1,'omitnan');
% p_mean(:,1)=(mean(p,1,'omitnan'))';
% for i=1:width(p)
%     idx=find(isnan(p(:,i)));
%     len(i,1)=length(p(:,1))-length(idx);
% end
% p_mean(:,2)=(std(p,0,1,'omitnan')')./len;

end