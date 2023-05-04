function [arry]=assign_num_3quadrants(input)
if iscell(input) && width(input)==1
    arry{1,1}=vertcat(input{4,1},input{5,1});
    arry{2,1}=vertcat(input{6,1},input{7,1},input{2,1},input{3,1});
    arry{3,1}=vertcat(input{1,1},input{8,1});
elseif iscell(input) && width(input)==2
    for i=1:2
        arry{1,i}=vertcat(input{4,i},input{5,i});
        arry{2,i}=vertcat(input{6,i},input{7,i},input{2,i},input{3,i});
        arry{3,i}=vertcat(input{1,i},input{8,i});
    end
else
    arry(1,:)=input(4,:)+input(5,:);
    arry(2,:)=input(6,:)+input(7,:)+input(2,:)+input(3,:);
    arry(3,:)=input(1,:)+input(8,:);
end
end