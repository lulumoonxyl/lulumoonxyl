function [PI_mean,PI_sem,PI_larva,PI_size,group]=bin_PI_cell(data1,series1,data,name)
%% this function will bin each cell of data based on the binning of series1
%and dat1

l=length(series1)-1;
PI_mean=zeros(1,l);
PI_sem=zeros(1,l);
%% group, PI_larva, PI_size are used for krusal-wallis test later
group=cell(l,1);
PI_larva=cell(l,1);
PI_size=zeros(l,1);

for i=1:l
    w=1;
    for j=1:length(data1)
        %loop through each cell of the data list and get the corresponding PI
        %based on the bin of series and data
        if i==l
            idx=find(data1{j,1}>=series1(i)&data1{j,1}<=series1(i+1));
        else
            idx=find(data1{j,1}>=series1(i)&data1{j,1}<series1(i+1));
        end
        
        if ~isempty(idx)

            PI_larvae=mean(data{j,1}(idx));
            PI_larva{i,1}(w,1)=PI_larvae;
            w=w+1;
        end
        clear idx
    end
    PI_size(i,1)=length(PI_larva{i,1});
    group{i,1}=repmat({append(name,"_",num2str(series1(i)),"-",num2str(series1(i+1)))},PI_size(i,1),1);
    PI_mean(1,i)=mean(PI_larva{i,1},'omitnan');
    idx=find(isnan(PI_larva{i,1}));
    PI_sem(1,i)=std(PI_larva{i,1},'omitnan')/sqrt(length(PI_larva{i,1})-length(idx));
    clear PI_larvae idx
end
end