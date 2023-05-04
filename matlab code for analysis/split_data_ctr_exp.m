function [ctr,exp,ctr_sem,exp_sem,ctr_name,exp_name,ctr_color,exp_color]=split_data_ctr_exp(data1,idx,name,color,varargin)
%there are two conditions:
% 1) data1 is a matrix containing the same info--> we
%then need a data2 to represent sem if you do calculate it
% the shape of data 1/2 should be n row*ii_num col
% 2) data1 is a cell, each cell contains two columns for mean and sem
%the shape of data1 should be {ii_num,1}-->which has n rows and 2 columns
% 3) if there is only 1 ctr group, the ctr,exp,ctr_sem,exp_sem output will be
% the same format as data1/2
%4) if there is more than 1 ctr groups, all these output will be cell arrays
data2=[];

for i=1:2:length(varargin)
    if strcmp(varargin{i},'sem')
        data2=varargin{i+1};
    elseif strcmp(varargin{i},'input_cond')
        input_cond=varargin{i+1};
    end
end

%idx is the index of the control group
if length(idx)==1
    %name1 is actually the names of the exp groups not the ctr group; the
    %same applies to the color1 variable
    ctr_sem={};exp_sem={};

    color1=color(idx);
    name1=name(idx);

    name(idx)=[];
    color(idx,:)=[];

    ctr_name=strings(length(name),1);
    exp_name=strings(length(name),1);
    ctr=cell(length(name),1);
    exp=cell(length(name),1);
    ctr_color=zeros(length(name),3);
    exp_color=zeros(length(name),3);
    %%  for condition 2) data1 is a cell
    if iscell(data1)
        ctr1=data1{idx,1};
        data1(idx)=[];
        
        for i=1:length(data1)
            ctr_name(i,1)=name1;
            exp_name(i,1)=name(i);

            ctr_color(i,:)=color1;
            exp_color(i,:)=color(i,:);

            exp{i,1}=data1{i,1};
            ctr{i,1}=ctr1;
        end

    elseif ismatrix(data1)
        %% for condition 1) data1/2 is a matrix
        ctr1=data1(:,idx);
        data1(:,idx)=[];
        for i=1:width(data1)

            ctr_name(i,1)=name1;
            exp_name(i,1)=name(i);

            ctr_color(i,:)=color1;
            exp_color(i,:)=color(i,:);

            ctr{i,1}=ctr1;
            exp{i,1}=data1(:,i);
        end
        % if there is an data2 input
        if ~isempty(data2)
            ctr_sem1=data2(:,idx);
            data2(:,idx)=[];
            for i=1:width(data2)
                exp_sem{i,1}=data2(:,i);
                ctr_sem{i,1}=ctr_sem1;
            end
        end
    end
else

    ctr_sem={};exp_sem={};
    l=length(name)-length(idx);
    ctr_name=strings(l,1);
    exp_name=strings(l,1);
    ctr=cell(l,1);
    exp=cell(l,1);
    ctr_color=zeros(l,3);
    exp_color=zeros(l,3);

    w=1;
    for i=1:length(idx)
        num2str(i)
        ind=find(contains(input_cond,num2str(i)));
        ind1=find(~ismember(ind,idx(i)));
        
        for j=1:length(ind1)

            if iscell(data1)
                ctr_name(w,1)=name(idx(i));
                exp_name(w,1)=name(ind(ind1(j)));
                ctr_color(w,:)=color(idx(i),:);
                exp_color(w,:)=color(ind(ind1(j)),:);

                ctr{w,1}=data1{idx(i),1};
                exp{w,1}=data1{ind(ind1(j)),1};
                name1(w,1)=name(ind(ind1(j)));
                color1(w,:)=color(ind(ind1(j)),:);
                w=w+1;
            elseif ismatrix(data1)
                ctr_name(w,1)=name(idx(i));
                exp_name(w,1)=name(ind(ind1(j)));
                ctr_color(w,:)=color(idx(i),:);
                exp_color(w,:)=color(ind(ind1(j)),:);

                ctr{w,1}=data1(:,idx(i));
                exp{w,1}=data1(:,ind(ind1(j)));

                name1(w,1)=name(ind(ind1(j)));
                color1(w,:)=color(ind(ind1(j)),:);

                if  ~isempty(data2)
                    ctr_sem{w,1}=data2(:,idx(i));
                    exp_sem{w,1}=data2(:,ind(ind1(j)));
                end
                w=w+1;
            end
        end

    end
end
end