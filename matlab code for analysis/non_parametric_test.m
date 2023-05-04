function [p_mwu,group_comp,p_k]=non_parametric_test(data,group,size,displayopt)
%data is the vertcat array of all data points
%group is the string array of the groups for all datapoints
%size is the number of datapoints for each group
len=length(size);
p_k=nan;

for i=1:len
    if i==1
        start(i,1)=1;
        stop(i,1)=size(i,1);
    else

        start(i,1)=start(i-1,1)+size(i-1,1);
        stop(i,1)=stop(i-1,1)+size(i,1);
    end
end
start
stop
size
if len>2
    if displayopt==1
        p_k=kruskalwallis(data,group); %first,use this test to check whether there is a difference between multiple group
    else
        p_k=kruskalwallis(data,group,'off');
    end
    if p_k<=0.05
        w=1;
        for i=1:len-1
            for j=i+1:len
                ind=isnan(data(start(i,1):stop(i,1)));
                ind1=isnan(data(start(j,1):stop(j,1)));
                disp(append('NAN length in the data for ',group(start(i,1)),num2str(len(ind))));
                disp(append('NAN length in the data for ',group(start(j,1)),num2str(len(ind1))));
                p_mwu(w,1)=ranksum(data(start(i,1):stop(i,1)),data(start(j,1):stop(j,1))); %use mann-whitney u-test to find which two groups have significant difference
                group_comp(w,1)=append(group(start(i,1)),"-",group(start(j,1)));
                w=w+1;

            end
        end
    else
        a=unique(group)
        p_mwu=nan;
        group_comp=" ";
        for i=1:length(a)
            group_comp =append(a(i),"_",group_comp);
        end
        warning (append("There is no significant difference between the data of different groups with Kruskal-Wallis test, p=", num2str(p_k)));
    end
elseif len==2
    %when there are less than 3 groups, perform the mann-whitney test
    %directly
    p_mwu=ranksum(data(start(1,1):stop(1,1)),data(start(2,1):stop(2,1))); %use mann-whitney u-test to find which two groups have significant difference
    group_comp=append(group(start(1,1)),"-",group(start(2,1)));
    if p_mwu>0.05
        warning (append("There is no significant difference between the preferential index of 2 groups with Mann-Whitney U-test, p=", num2str(p_mwu)));

    end
else
    disp("There is only one group for comparison. Cannot perform pairwise test!");
    p_mwu=nan; 
    group_comp=nan;
end
end