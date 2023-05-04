function [deg_all,deg_mean,deg_sem]=bin_data_mean2(data1,series1,data2,series2,data)
%this function will bin data1/2 according to series 1/2 and calculate the
%mean of the data in each bin combination
%if data1 is orientation, series1 needs to be -180:10:180 for the heatmap
%plot
l=length(series1)-1;
l2=length(series2)-1;
deg_all=cell(l,l2);
deg_mean=nan(l,l2);
%len=0;

for i=1:l
    for z=1:l2
        %create an array to store the data for each HD and x/y/t bin
        %combination

        if z==l2 && i==l
            ind1=find(data1>=series1(i)&data1<=series1(i+1)&data2>=series2(z)&data2<=series2(z+1));
        elseif z~=l2 && i~=l
            ind1=find(data1>=series1(i)&data1<series1(i+1)&data2>=series2(z)&data2<series2(z+1));
        elseif z==l2
            ind1=find(data1>=series1(i)&data1<series1(i+1)&data2>=series2(z)&data2<=series2(z+1));
        elseif i==l
            ind1=find(data1>=series1(i)&data1<=series1(i+1)&data2>=series2(z)&data2<series2(z+1));
        end
        %if there are turning happens in the bin
        if ~isempty(ind1)
            deg_all{i,z}=data(ind1);
            deg_mean(i,z)=mean(data(ind1));
            deg_sem(i,z)=std(data(ind1))/sqrt(length(data(ind1)));
        end
        clear ind1
    end
end

end