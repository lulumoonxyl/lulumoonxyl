for ii=1:5

% choredir='D:\choreography-result\CS@CS_cap_test\CS@CS_cap test_exp2'
choredir='D:\choreography-result\CS@CS_cap_test\CS@CS_EA_gen0 to gen2'

%% get filelist

full = dir(choredir);
n = 1;
sear = fullfile(choredir,'*');
d = 1;
while length(d) > 0
    sear = fullfile(sear,'*');
    d = dir(sear);
    full = vertcat(full,d);
    n = n+1;
end
%
filt = [full.isdir];
full = full(filt);
names = {full.name}';

expr = '^\d\d\d\d\d\d\d\d_\d\d\d\d\d\d';
filt = regexp(names,expr);
filt = cellfun(@(x) ~isempty(x), filt);
d = full(filt);


%% group genotypes

names = {d.folder}';
splits = cellfun(@(x) split(x,'\'), names, 'UniformOutput', false);
%cellfun: apply a function to the array
splits = cellfun(@(x) [x{end-1},'\',x{end}], splits, 'UniformOutput', false);
[uname,na,nb] = unique(splits);

%% plot loop

% import specs
delimiter = ' ';
startRow = 0;
formatSpec = '%s%f%f%f%[^\n\r]';

%% -------------------------------------------------------------------%

idx = find(nb == ii);
% import chore data for first folder (100uW stimulation)
et = {};
dat_x= {};
animal_num={};
dat_y ={};
speed ={};


for jj = idx'
    dirname = fullfile(d(jj).folder,d(jj).name);
    fname = dir(fullfile(dirname,['*x.dat']));
    f = contains({fname.name},['.x.']);
    try
        fname = fullfile(fname(f).folder,fname(f).name);
    catch
        fprintf(['empty folder ' dirname]);
        continue
    end
    
    fname_y = dir(fullfile(dirname,['*y.dat']));
    f_y = contains({fname_y.name},['.y.']);
    try
        fname_y = fullfile(fname_y(f_y).folder,fname_y(f_y).name);
    catch
        fprintf(['empty folder ' dirname]);
        continue
    end
    
    fname_sp = dir(fullfile(dirname,['*speed.dat']));
    f_sp  = contains({fname_sp .name},['.speed.']);
    try
        fname_sp  = fullfile(fname_sp (f_sp ).folder,fname_sp (f_sp ).name);
    catch
        fprintf(['empty folder ' dirname]);
        continue
    end
    
    fileID = fopen(fname,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
    fclose(fileID);
    
    fileID_y = fopen(fname_y,'r');
    dataArray_y = textscan(fileID_y, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
    fclose(fileID_y);
    
    fileID_sp = fopen(fname_sp,'r');
    dataArray_sp = textscan(fileID_sp, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
    fclose(fileID_sp);
    
    animal_num=vertcat(animal_num, dataArray{:,2}); %add animal num
    et = vertcat(et,dataArray{:,3});
    speed =vertcat(speed,dataArray_sp{:,4})
    dat_x = vertcat(dat_x,dataArray{:,4});
    dat_y = vertcat(dat_y,dataArray_y{:,4});
    clear dataArray dataArray_y dataArray_sp
    
    
end

%delete the cell where x and y do not have the same length
% for i=1:length(dat_x)
%     lx=length(dat_x{i,1});
%     ly=length(dat_y{i,1});
%     lsp=length(speed{i,1});
%     if lx==ly&&lx==lsp
%         continue 
%     else 
%         dat_x{i,1}={}
%         dat_y{i,1}={}
%         speed{i,1}={}
%         et{i,1}={}
%         animal_num{i,1}={}
%     end 
% end 

x1=vertcat(dat_x{:});

y1=vertcat(dat_y{:});
t1=vertcat(et{:});
% % % % 
% % % %  ind=find(t>=600 &t<900);
figure(3);
hold on;
 subplot(1,2,2);
plot(x1,y1,'b.','MarkerSize',2);
title(uname(ii));
axis tight
axis square;
hold off

   
    
for i=1:2
    figure(1)
    subplot(1,2,i)
    xlim([-200 200]);
    ylim([-200 200]);
end 
    

    
%     for i=1:length(dat_x)
%         ind=find(dat_x{i,1}<=20|dat_x{i,1}>=240|dat_y{i,1}>=220|dat_y{i,1}<=5);
%         if isempty(ind)
%             continue
%         else
%             dat_x{i,1}(ind)=[];
%             dat_y{i,1}(ind)=[];
%             speed{i,1}(ind)=[];
%             animal_num{i,1}(ind)=[];
%             et{i,1}(ind)=[];
%         end
%         clear ind
%         
%     end
    
    dat_x_correct={};
    dat_y_correct={};
 for i=1: length(animal_num) 
%      if isempty(animal_num{i,1})
%          continue
%      end 
     [N{i},na1{i},nb1{i}] = unique(animal_num{i});
 
 for j=1:length(N{i})
     x= dat_x{i}(na1{i}(j));
     y= dat_y{i}(na1{i}(j));
     
     if j==length(N{i})
      l=length(dat_x{i});
         dat_x_correct{i,1}(na1{i}(j):l,1)=dat_x{i}(na1{i}(j):l)-x;
         dat_y_correct{i,1}(na1{i}(j):l,1)=dat_y{i}(na1{i}(j):l)-y;
         
     else
         
         
         dat_x_correct{i,1}(na1{i}(j):na1{i}(j+1)-1,1)=dat_x{i}(na1{i}(j):na1{i}(j+1)-1)-x;
         dat_y_correct{i,1}(na1{i}(j):na1{i}(j+1)-1,1)=dat_y{i}(na1{i}(j):na1{i}(j+1)-1)-y;
     end
     clear x y 
 end 
 end 

% %% plot the trajectory for the whole period

x2=vertcat(dat_x_correct{:});
y2=vertcat(dat_y_correct{:});
figure(1)
hold on;

subplot(1,2,ii)
plot(x3,y2,'b.','MarkerSize', 1);
title(uname(ii));
axis tight
xline(0,'--r');
axis square;
         



% for j=1:4
% subplot(2,3,6)
%  xlim([0 225]);
% axis tight
% plot(dat_x{j},dat_y{j},'b.','MarkerSize', 1);
%          title('Trajectories for 10^-^3 GA,white');
% end        
%% calculate the x velocity
% for i=1:length(dat_x)
% dat_x{i}=-dat_x{i};
% end 
velocity={};
velocity_x_half={};
velocity_t_three={};
for i= 1:length(animal_num)
    
    for j=1:length(N{i})
        if j==length(N{i})
            x=dat_x{i}(na1{i}(j):end);
            t=et{i}(na1{i}(j):end);
            sp = speed{i}(na1{i}(j)+1:end);
        else
            x=dat_x{i}(na1{i}(j):na1{i}(j+1)-1);
            t=et{i}(na1{i}(j):na1{i}(j+1)-1);
            sp = speed{i}(na1{i}(j)+1:na1{i}(j+1)-1);
        end
        %calculate the index in the x direction
        x_diff=diff(x);
        t_diff=diff(t);
        l=length(x);
        x1=x(2:l);
       
        t1=t(2:l);
        
        ind0=find(x1>=125)%larvae is on the odor side
        ind2=find(x1<125)% when the x position is smaller than -113, the larvae is on the side opposite to the odor
        
        ind3=find(t1>=0&t1<300);
        ind4=find(t1>=300&t1<600);
        ind5=find(t1>=600&t1<900);
        
        v=(x_diff./t_diff)./sp;
        ind1= find(sp==0);
        v(ind1)=nan;
        velocity{i,1}{j,1}=v;
        
        velocity_x_half{i,1}{j,1}=v(ind0);%the first column is on the odor side
        velocity_x_half{i,2}{j,1}=v(ind2)%opposite side of the odor
        
        velocity_t_three{i,1}{j,1}=v(ind3);
        velocity_t_three{i,2}{j,1}=v(ind4);
        velocity_t_three{i,3}{j,1}=v(ind5);
        clear x t x_diff t_diff sp ind1 v ind0 ind2 x1 ind3 ind4 ind5
    end
    
end

if ii==4
    ind0=1
elseif ii==3
    ind0=2
elseif ii==2
    ind0=3
else ii==1
    ind0=4

end
v1=vertcat(velocity{:});
v2=vertcat(v1{:});% this is for calculating the index for each group
ind=find(isnan(v2));
v2(ind)=[];
index(ind0)=mean(v2);

st(ind0)=std(v2);
SEM(ind0)=st(ind0)/sqrt(length(v2))

clear v1 v2 st ind
% we then need to calculate the index for each exp and plot them 
for i=1:length(velocity)
    v=vertcat(velocity{i,1}{:});
    ind=find(isnan(v));
    v(ind)=[];
    index_exp{ind0,1}(i,1)=mean(v);
    SEM_exp{ind0,1}(i,1)=(std(v))/sqrt(length(v));
    clear v ind
end 

% calculate the overall index for larvae in different position 
v3=vertcat(velocity_x_half{:,1});
v4=vertcat(velocity_x_half{:,2});
v5{:,1}=vertcat(v3{:});
v5{:,2}=vertcat(v4{:});

for i=1:2
    ind=find(isnan(v5{1,i}));
    v5{1,i}(ind)=[];
    index_x_half(ind0,i)=mean(v5{1,i});
    SEM_x_half(ind0,i)=std(v5{1,i})/sqrt(length(v5{1,i}));
    clear ind
end 



clear v3 v4 v5
for i=1:length(velocity)
    for j=1:2
        v=vertcat(velocity_x_half{i,j}{:});
        ind=find(isnan(v));
        v(ind)=[];
        index_x_half_exp{ind0,j}(i,1)=mean(v);
        SEM_x_half_exp{ind0,j}(i,1)=(std(v))/sqrt(length(v));
        clear v ind
    end
    
end

% calculate the overall index for larvae in different time
v6=vertcat(velocity_t_three{:,1});
v7=vertcat(velocity_t_three{:,2});
v8=vertcat(velocity_t_three{:,3});
v9{:,1}=vertcat(v6{:});
v9{:,2}=vertcat(v7{:});
v9{:,3}=vertcat(v8{:});

for i=1:3
    ind=find(isnan(v9{1,i}));
    v9{1,i}(ind)=[];
    index_t_three(ind0,i)=mean(v9{1,i});
    SEM_t_three(ind0,i)=std(v9{1,i})/sqrt(length(v9{1,i}));
    clear ind;
end 


clear v6 v7 v8 v9

for i=1:length(velocity)
    for j=1:3
        v=vertcat(velocity_t_three{i,j}{:});
        ind=find(isnan(v));
        v(ind)=[];
        index_t_three_exp{ii,j}(i,1)=mean(v);
        SEM_t_three_exp{ii,j}(i,1)=(std(v))/sqrt(length(v));
        clear v ind
    end
    
end

clearvars -except SEM index index_exp SEM_exp SEM_x_half index_x_half index_t_three_exp SEM_t_three_exp index_x_half_exp SEM_x_half_exp...
    index_t_three SEM_t_three

end 

for i=1:5
    figure(1)
    hold on;
    subplot(2,3,5)
    title('two')
    xlim([-200 200]);
    ylim([-200 200]);
end 
%% arrange the plot for the index and SEM of each exp 
figure(2)
hold on 
   b=bar(index2)
   errorbar(index2,SEM2,'.')
   b.FaceColor='flat';
%     b.CData(2,:)=[1 1 1];
for i=1:numel(index_exp2)
    l=length(index_exp2{1,i});
    plot(ones(l,1)*i,index_exp2{1,i},'k.','MarkerSize',10);
%     errorbar(ones(l,1)*i,index_exp1{i,1},SEM_exp1{i,1},'b.');
end 

xlim([0 5]);
xticks([0:1:5]);
% yline(0,'r--');
xticklabels({' ',' H2O','EA','H2O','GA'})

ylabel('v_x/<s>');

axis square
xtickangle(30);

figure(2)
title('eliminate edge')

%% plot the bar for the different time
figure(3)

hold on;
hBar = bar(index_t_three2,1);  
% Return �bar� Handle
for k1 = 1:size(index_t_three2,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: �XOffset� Is An Undocumented Feature; This Selects The �bar� Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
end

errorbar(ctr, ydt,SEM_t_three2' , '.r')  
legend('0-300s','300-600s','600-900s')

[rol,col]=size(index_t_three_exp2);
for i=1:col
    figure(3)
   
    hold on
    for j=1:rol
        l=length(index_t_three_exp2{j,i});
        plot(ones(l,1)*ctr(i,j),index_t_three_exp2{j,i},'b.','MarkerSize',10);
%         errorbar(ones(l,1)*j,index_t_three_exp1{j,i},SEM_t_three_exp1{j,i},'b.');
    end


end
    xlim([0 5]);
    xticks([0:1:5]);
    ylabel('v_x/<s> ');
    xticklabels({' ','3 caps, H2O','3 caps, 10^-^2 GA','5 caps, H2O','5 caps, 10^-^2 GA'})
    xticklabels({' ','exp 1, EA','exp 1, H2O','exp 2, EA','exp 2, H2O','exp 3, EA','exp 3, H2O'})
    xtickangle(30);
 clear ctr ydt k1 hBar

%% plot the index and sem for position 
figure(1)

hold on 
hBar = bar(index_x_half2,1);  

for k1 = 1:size(index_x_half2,2)
    ctr(k1,:) = bsxfun(@plus, hBar(1).XData, hBar(k1).XOffset');    % Note: �XOffset� Is An Undocumented Feature; This Selects The �bar� Centres
    ydt(k1,:) = hBar(k1).YData;                                     % Individual Bar Heights
end
hold on
errorbar(ctr, ydt,SEM_x_half2','.r') 
legend('Odor side','Opposite side')% Return �bar� Handle

[rol,col]=size(index_x_half_exp2);
for i=1:col
    figure(1)
    
    hold on
    for j=1:rol
        l=length(index_x_half_exp2{j,i});
        plot(ones(l,1)*ctr(i,j),index_x_half_exp2{j,i},'b.','MarkerSize',10);
%         errorbar(ones(l,1)*j,index_x_half_exp{j,i},SEM_x_half_exp{j,i},'b.');
    end

   

end
    xlim([0 3]);
     xticks([0:1:3]);
    ylabel('v_x/<s> for control, nothing');
    xticklabels({' ','gen0, 5*10^-^3 GA','gen0, 5*10^-^3 Oil','gen1, 10^-^3 GA','gen1, 10^-^3 Oil','gen2, 10^-^3 GA','gen2, 10^-^3 Oil'})
    xticklabels({' ','3 caps, H2O','3 caps, 10^-^2 GA','5 caps, H2O','5 caps, 10^-^2 GA'});
    xticklabels({'','gen0,10^-5 EA','gen0, H2O','gen1, 10^-^5 EA','gen1, H2O','gen2, 10^-^5','gen2, H2O'})
    xtickangle(30);
 clear ctr ydt k1 hBar
%% arrange the plot
%white: Oil, GA (low to high), H2O, EA
%grape:Oil, GA (low to high), H2O, EA
%gen 0?Oil, GA (low to high), H2O, EA
 index1(1,5)=index(1);
 SEM1(1,5)=SEM(1);
 a(3:end)=SEM1(1,2:end)
 SEM1=a
 clear index st SEM 
x=categorical ({'one','two','three','four','five'});
x=reordercats(x,{'one','two','three','four','five'});
 x=categorical({'White, 5*10^-^3 Oil','White,10^-^3 GA','White, 5*10^-^3 GA','White, 10^-^2 GA','White, H2O','White, 10^-^5 EA',...
    'Grape, 5*10^-^3 Oil','Grape, 10^-^3 GA','Grape, 5*10^-^3 GA',...
    'Gen0, Brown, 5*10^-^3 Oil','Gen0, Brown, 10^-^3 GA','Gen0, Brown, 5*10^-^3 GA','Gen0, Brown, H2O','Gen0, Brown, 10^-^5 EA',...
    'Gen1, Brown, 10^-^3 Oil','Gen1, Brown, 10^-^3 GA','Gen1, Brown, H2O','Gen1, Brown, 10^-^5 EA','Gen2, Brown, 10^-^3 Oil','Gen2, Brown, 10^-^3 GA'
    });
x=reordercats(x,{'White, 5*10^-^3 Oil','White,10^-^3 GA','White, 5*10^-^3 GA','White, 10^-^2 GA','White, H2O','White, 10^-^5 EA',...
    'Grape, 5*10^-^3 Oil','Grape, 10^-^3 GA','Grape, 5*10^-^3 GA',...
    'Gen0, Brown, 5*10^-^3 Oil','Gen0, Brown, 10^-^3 GA','Gen0, Brown, 5*10^-^3 GA','Gen0, Brown, H2O','Gen0, Brown, 10^-^5 EA',...
    'Gen1, Brown, 10^-^3 Oil','Gen1, Brown, 10^-^3 GA','Gen1, Brown, H2O','Gen1, Brown, 10^-^5 EA','Gen2, Brown, 10^-^3 Oil','Gen2, Brown, 10^-^3 GA'});
x=categorical({ 'White, 5*10^-^3 Oil','White,10^-^3 GA','White, H2O','White, 10^-^5 EA'})
x=reordercats(x,{ 'White, 5*10^-^3 Oil','White,10^-^3 GA','White, H2O','White, 10^-^5 EA'});

ind=[1 2 5 6];
x=categorical({'Gen2, Brown, 10^-^3 GA','Gen2, Brown, 10^-^3 Oil','White,10^-^3 GA'});

x=reordercats(x,{'Gen2, Brown, 10^-^3 GA','Gen2, Brown, 10^-^3 Oil','White,10^-^3 GA'});
figure(2);
 hold on;
%  b=bar(x,index1(2,:))
%  errorbar(index1(2,:),SEM1(2,:),'.');
   b=bar(index)
   errorbar(index,SEM,'.')
%    axis square
 ylabel('v_x/<s> for 5*10^-^3 GA');
 b.FaceColor='flat';
  b.CData(1,:)=[1 1 1];
  b.CData(1,:)=[0 0 0];
  b.CData(2,:)=[1 0 0];
  axis square
  ylim([-0.05 0.17])
  xlim([0 6])
 
 l=[5 13 17];% h2o is black histogram [0 0 0]
 l=[6 14 18 ];%EA is the red one [1 0 0]
 l=[1 7 10 15 19];%oil is white [1 1 1]
 
for i=1:length(l)
 b.CData(l(i),:)=[1 1 1];
end 
%% get the number of larvae tracked per area as an index
figure(1)
for i=1:2
 subplot(2,1,i)
hold on;
h=findobj(gca,'Type','Line')
x=h.XData;
y=h.YData;

PI{i,1}=y(1790:1800)'
index(i)=mean(PI{i,1});
SEM(i)=std(PI{i,1})/sqrt(11);
clear x y h 
end 

index1(2,13)=index(1)
SEM1(2,13)=SEM(1)

clear index SEM PI x b i
%% plot the index of each exp 
%% calculate index for different time and different area
%filter shorter one --may not need to do this 
%% calculate the probability of tracking time 