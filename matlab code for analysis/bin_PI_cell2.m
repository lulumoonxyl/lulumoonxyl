function [PI_mean,PI_sem]=bin_PI_cell2(data1,series1,data2,series2,data)
%% this function will bin each cell of data based on the binning of series1
%and dat1

height=length(series1)-1;
width=length(series2)-1;
PI_mean=zeros(height,width);
PI_sem=zeros(height,width);

for i=1:height
    for k=1:width
        PI_larva=[];
        w=1;
        for j=1:length(data1)
            %loop through each cell of the data list and get the corresponding PI
            %based on the bin of series and data
            if i==height && k==width
                ind=find(data1{j,1}>=series1(i)&data1{j,1}<=series1(i+1)&data2{j,1}>=series2(k)&data2{j,1}<=series2(k+1));
            elseif i==height
                ind=find(data1{j,1}>=series1(i)&data1{j,1}<=series1(i+1)&data2{j,1}>=series2(k)&data2{j,1}<series2(k+1));
            elseif k==width
                ind=find(data1{j,1}>=series1(i)&data1{j,1}<series1(i+1)&data2{j,1}>=series2(k)&data2{j,1}<=series2(k+1));
            else
                ind=find(data1{j,1}>=series1(i)&data1{j,1}<series1(i+1)&data2{j,1}>=series2(k)&data2{j,1}<series2(k+1));
            end

            if ~isempty(ind)
                PI_larvae=mean(data{j,1}(ind));
                PI_larva(w,1)=PI_larvae;
                w=w+1;
                clear ind
            end
        end
        if ~isempty(PI_larva)

            PI_mean(i,k)=mean(PI_larva,'omitnan');
            idx=find(isnan(PI_larva));
            PI_sem(i,k)=std(PI_larva,'omitnan')/sqrt(length(PI_larva)-length(idx));
            clear PI_larva idx
        end
    end
end
end