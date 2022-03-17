function [result,w,u]=calculate_tf_reorientation_angle(dat_JB,result,idx,idx2,i,j,w,u,deg,str,large)
if isequal(large,'large')
    name={'turning_frequency','turning_frequency_large','large_over_all','theta'};
else
    name={'turning_frequency','turning_frequency_small','large_over_all','theta'};
end
name=append(name,str);
%calculate the total amount of time that this larva spends in this xbin
tt=sum(dat_JB.et{j,1}(idx)-dat_JB.et{j,1}(idx-1));

%calculate the turning frequency for each larvae is each x bins
tf_all=length(idx2)/tt;
%calculate the large turning frequency (theta>=deg)
if isequal(large,'large')
    %get the number of large turning events from idx 2
    idx3=find(dat_JB.reorientation_angle{j,1}(idx2)>=deg);
else
    idx3=find(dat_JB.reorientation_angle{j,1}(idx2)<deg);
end
tf_large=length(idx3)/tt;

result.(name{1}){i,1}(w,1)=tf_all;
result.(name{2}){i,1}(w,1)=tf_large;
w=w+1;

if ~isempty(idx2)
    result.(name{3}){i,1}(u,1)=length(idx3)/length(idx2);
    result.(name{4}){i,1}(u,1)=mean(dat_JB.reorientation_angle{j,1}(idx2));
    u=u+1;
end
end