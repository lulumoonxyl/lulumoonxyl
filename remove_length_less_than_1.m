function dat=remove_length_less_than_1(dat,fields)
%this function is before calculating the tsum 
%it will remove any subarray with length less than 1
et=find(matches(fields,'et'));
for i=length(dat{et,1}):-1:1
    if length(dat{et,1}{i,1})==1
        for j=1:length(dat)
            dat{j,1}(i)=[];
        end 
    end 
end 
end 