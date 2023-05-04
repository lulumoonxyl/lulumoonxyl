function [r_mean,pos]=bin_data_mean(data1,series1,data)
%this function will bin data1 based on series1 and calculate the mean value
%from the corresponding results using data

l=length(series1)-1;
r_mean=nan(l,2);
sum=0;
pos=cell(l,1);
for i=1:l

    if i==l
        ind1=find(data1>=series1(i)&data1<=series1(i+1));
        sum=sum+length(ind1);
    else
        ind1=find(data1>=series1(i)&data1<series1(i+1));
        sum=sum+length(ind1);
    end
    %if there are turning happens in the bin
    if ~isempty(ind1)
         r_mean(i,1)=mean(data(ind1),'omitnan');
         idx=find(isnan(data(ind1)));
         r_mean(i,2)=std(data(ind1),'omitnan')/sqrt(length(ind1)-length(idx));
         clear idx
    end

    pos{i,1}=data(ind1);%this gathers all the data point for each bin
    clear ind1

end
end