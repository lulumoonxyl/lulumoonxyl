%% stat test for odor concentration
%load excel
T=readtable('C:\Users\Fei Huang\OneDrive - McGill University\matlab\Analysis based on larvae\olfactory_output_pdf\odor concentration.xlsx','UseExcel',true,'Sheet','trial 9(data for analysis)');
%% 1-1) compute the diffference in concentration for each timepoint across groups
cond=["H2O";"c10n3GA";"c10n2GA";"c10n1GA";"c10n5EA";"c10n3EA";"c10n2EA"];
time=["0min";"3min";"6min";"9min";"12min";"15min"];
pos=["row1";"row3";"row7"];
name=T.Properties.VariableNames;
p_mwu=cell(42,1); group_comp=cell(42,1); p_k=cell(42,1);
for i=1:7
    for j=1:6
        a=append(cond(i),'_',time(j));
        idx=find(name==a)
        b=T.(a);
        group=convertCharsToStrings(b);
        %remove empty string
        group(cellfun('isempty',group))=[];

        data=table2array(T(:,idx-1));
        data(find(isnan(data)))=[];

        size=zeros(3,1);
        r=["row1";"row3";"row7"];
        for z=1:3
            size(z,1)=length(find(contains(group,r(z))));
        end
        [p_mwu{6*i-6+j,1},group_comp{6*i-6+j,1},p_k{6*i-6+j,1}]=non_parametric_test(data,group,size,0);
    end
end


group_pos=cell(6,3);
dat_pos=cell(6,3);
size_pos=zeros(6,3);

p_mwu_pos=cell(42,1); group_comp_pos=cell(42,1); p_k_pos=cell(42,1);
for i=1:6
    a=time(i);
    idx=find(contains(name,a));
    cond_name=name(idx);

    ind_str=find(rem(idx,2)==0)
    ind_val=find(rem(idx,2)==1)
    for j=1:length(ind_str)
        b=T.(cond_name{1,ind_str(j)});
        data=T.(cond_name{1,ind_val(j)});
        data(find(isnan(data)))=[];
        % get the string for each condition
        group=convertCharsToStrings(b);
        group(cellfun('isempty',group))=[];
        % select the one with either 'row1','row3','row7'
        % then group them differently
        % also group the data

        for k=1:length(pos)
            ind_pos=find(contains(group,pos(k)));
            group_pos{i,k}=[group(ind_pos);group_pos{i,k}];
            dat_pos{i,k}=[data(ind_pos);dat_pos{i,k}];
            size_pos(i,k)=length(dat_pos{i,k});
       end
    end
    
end

