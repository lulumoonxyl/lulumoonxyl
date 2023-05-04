function [c]=bin_data_count_tracked(data,series)
%this function will bin data according to series  and count the number of
%data in each bin combination

l=length(series)-1;
c=zeros(l,1);

for i=1:l
    for j=1:length(data)
        if i==l
            ind1=find(data{j,1}>=series(i)&data{j,1}<=series(i+1));
        else
            ind1=find(data{j,1}>=series(i)&data{j,1}<series(i+1));
        end
    
    %simply add one to each bin if the ind1 is not 0
    if ~isempty(ind1)
        c(i,1)=c(i,1)+1;
    end
    clear ind1
    end 
end
end
