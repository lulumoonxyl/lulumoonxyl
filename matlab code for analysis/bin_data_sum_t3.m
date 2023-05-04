function [t,p,p_bin,tsum,ttl_t]=bin_data_sum_t3(data1,data2,data3,series1,series2,series3,et,varargin)
%this function will bin the data based on two series and then count all the
%time that it spend in each x and hd bins etc
%data1 and series1 is the orientation
%data2/3 and series2/3 is the x or y or t
%ttl_t is the total tracking time of each animal =tsum
%t_series2 is the tracking time for each object spend in each bin in
%series2

l=length(series1)-1;
ld=length(data1);
l2=length(series2)-1;
l3=length(series3)-1;
t0_idx={};t1_idx={}; %if this is not [], it will eliminate all the time the object spends on turning

for i=1:2:length(varargin)

    if strcmp(varargin{i},"t0_idx")
        t0_idx=varargin{i+1};
    elseif strcmp(varargin{i},"t1_idx")
        t1_idx=varargin{i+1};
    end
end
%create empty arry for the t for each larvae and each bin
t=NaN(l,l2,l3);
% t_series3=zeros(ld,l2,l3);
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
    t_an=NaN(l,l2,l3);

    %bin the data1 according to the series1, and data2/3 according to
    %series2/3

    for i=1:l
        for z=1:l2
            for k=1:l3
                if z==l2 && i==l && k==l3
                    ind1=find(data1{j,1}>=series1(i)&data1{j,1}<=series1(i+1)&data2{j,1}>=series2(z)&data2{j,1}<=series2(z+1)&data3{j,1}>=series3(k)&data3{j,1}<=series3(k+1));
                elseif z~=l2 && i~=l && k~=l3
                    ind1=find(data1{j,1}>=series1(i)&data1{j,1}<series1(i+1)&data2{j,1}>=series2(z)&data2{j,1}<series2(z+1)&data3{j,1}>=series3(k)&data3{j,1}<series3(k+1));
                elseif z==l2 && k==l3
                    ind1=find(data1{j,1}>=series1(i)&data1{j,1}<series1(i+1)&data2{j,1}>=series2(z)&data2{j,1}<=series2(z+1)&data3{j,1}>=series3(k)&data3{j,1}<=series3(k+1));
                elseif i==l && k==l3
                    ind1=find(data1{j,1}>=series1(i)&data1{j,1}<=series1(i+1)&data2{j,1}>=series2(z)&data2{j,1}<series2(z+1)&data3{j,1}>=series3(k)&data3{j,1}<=series3(k+1));
                elseif i==l && z==l2
                    ind1=find(data1{j,1}>=series1(i)&data1{j,1}<=series1(i+1)&data2{j,1}>=series2(z)&data2{j,1}<=series2(z+1)&data3{j,1}>=series3(k)&data3{j,1}<series3(k+1));
                elseif z==l2
                    ind1=find(data1{j,1}>=series1(i)&data1{j,1}<series1(i+1)&data2{j,1}>=series2(z)&data2{j,1}<=series2(z+1)&data3{j,1}>=series3(k)&data3{j,1}<series3(k+1));
                elseif i==l
                    ind1=find(data1{j,1}>=series1(i)&data1{j,1}<=series1(i+1)&data2{j,1}>=series2(z)&data2{j,1}<series2(z+1)&data3{j,1}>=series3(k)&data3{j,1}<series3(k+1));
                elseif k==l3
                    ind1=find(data1{j,1}>=series1(i)&data1{j,1}<series1(i+1)&data2{j,1}>=series2(z)&data2{j,1}<series2(z+1)&data3{j,1}>=series3(k)&data3{j,1}<=series3(k+1));
                end
                if ~isempty(ind1)
                    if ind1(1)==1
                        ind1(1)=[];
                    end

                    if ~isempty(ind1)
                        ind2=find(ismember(ind1,idx));
                        ind1(ind2)=[];
                        if ~isempty(ind1)
                            t_an(i,z,k)=sum(et{j,1}(ind1)-et{j,1}(ind1-1));

                            % this will generate a matrix of the time each object spend in each bin combination
                            clear ind*
                        end
                    end
                end
            end
        end
    end

    %% sum all the time of one object
    tsum(j,1)=sum(t_an,'all',"omitnan");
    %% sum the time of one object spent in the binning of series2
%     t_series3(j,:,:)=sum(t_an,1,'omitnan');
    %% add the time of this object to previos object:variable t
    T=cat(4,t_an,t);
    t=sum(T,4,'omitnan');
    clear T t_an t_an_sum

end

%% sum all the time in the t matrix to confirm we do not count something
%twice

p=t./sum(tsum,'omitnan');%divide by total time
p_bin=t./sum(t,1,'omitnan');
end