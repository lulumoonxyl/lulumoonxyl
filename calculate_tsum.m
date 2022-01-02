function [dat,fields]=calculate_tsum(dat,fields,n)
%this function will get the tsum for each larvae and remove the cells that
%have tracking time less than n s
l=length(fields);
fields{l+1,1}='tsum';

et_pos=find(matches(fields,'et'));

for i=1:length(dat{et_pos,1})
    dat{l+1,1}(i,1)=dat{et_pos,1}{i,1}(end)-dat{et_pos,1}{i,1}(1);
end 
if isempty(n)==0
    idx=find(dat{l+1,1}<n);
    for j=1:length(dat)
        dat{j,1}(idx)=[];
    end 
    clear idx
end 

end