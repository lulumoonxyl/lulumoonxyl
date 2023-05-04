function [c]=bin_data_count2(data1,series1,data2,series2)
%this function will bin data1/2 according to series 1/2 and count the number of
%data in each bin combination
%if data1 is orientation, series1 needs to be -180:10:180 for the heatmap
%plot
%data1 and data2 should be an array not a cell array
l=length(series1)-1;

l2=length(series2)-1;
c=nan(l,l2);

for i=1:l
    for z=1:l2
        if z==l2 && i==l
            ind1=find(data1>=series1(i)&data1<=series1(i+1)&data2>=series2(z)&data2<=series2(z+1));
        elseif z~=l2 && i~=l
            ind1=find(data1>=series1(i)&data1<series1(i+1)&data2>=series2(z)&data2<series2(z+1));
        elseif z==l2
            ind1=find(data1>=series1(i)&data1<series1(i+1)&data2>=series2(z)&data2<=series2(z+1));
        elseif i==l
            ind1=find(data1>=series1(i)&data1<=series1(i+1)&data2>=series2(z)&data2<series2(z+1));
        end
        %simply add the number of turning event into the big matrix
        c(i,z)=length(ind1);
        clear ind1
    end
end


end