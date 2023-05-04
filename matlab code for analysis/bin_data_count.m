function [c,pos]=bin_data_count(data,series)
%this function will bin data according to series  and count the number of
%data in each bin combination

l=length(series)-1;
c=nan(l,1);
pos=cell(l,1);

for i=1:l
    if i==l
        ind1=find(data>=series(i)&data<=series(i+1));
    else
        ind1=find(data>=series(i)&data<series(i+1));
    end 
    %simply add the number of turning event into the big matrix
    c(i,1)=length(ind1);
    pos{i,1}=ind1; %save the index of the data happens in each bin
    clear ind1

end
end
