% this code will get the x,y,midline and AN and group all the timestamps
% together
function group_choreo_data(tracker_num,genotype,condition,odor_side,varargin)
%% odor side=='left', then we will flip the x postion to make the odor in the right side just to make it easier to intepret plots
% varargin input: edge should have the strcut like [xmin, xmax;ymin ymax];
% You need to check where the locotion of odorants is using the chore.jar
% keep in mind that this code does no take into account that the odor is on
% the top or down
%% identify the input and output folder
choredir=fullfile("/project/6010970/screen/choreography-results",tracker_num,genotype,condition);
outdir=fullfile("/project/6010970/screen/olfactory_output/choreography",tracker_num,genotype,condition)
% choredir=fullfile("D:\choreography-results",tracker_num,genotype,condition);
% outdir=fullfile("D:\olfactory_output\choreo\",genotype,condition);
edge=[];
for i=1:2:length(varargin)
    if strcmp(varargin{i},'edge')
        edge=varargin{i+1};
    end 
end 
%% load data
%get raw choreography data (x,y,midline) from the .txt files
dat=get_data_from_dat_file(choredir);
%% concatenate all the data
%later
fn=fieldnames(dat);
for i =1:length(fn)
    dat.(fn{i})=vertcat(dat.(fn{i}){:});
end
clear fn
xmax=max(dat.x);
ymax=max(dat.y);

%if odor is in the left side, we need to flip it left and right in order to
%plot the odor in the right side
if strcmp(odor_side,"left")% left is defined whaen draw with the map function
    dat.x=xmax-dat.x;
end
%move the data to (0,0)
ymin=min(dat.y);
xmin=min(dat.x);
dat.x=dat.x-xmin;
dat.y=dat.y-ymin;
%% grouped data based on their aninmal number (AN)
idx=find(diff(dat.AN)~=0);

fn={"AN","et","x","y","mid"};
idx=[0;idx;length(dat.AN)];
for j=1:length(fn)
    for i=1:length(idx)-1
        if strcmp(fn{j},"mid")||strcmp(fn{j},"AN")
            dat_grouped.(fn{j})(i,1)=dat.(fn{j})(idx(i)+3);
        else
            dat_grouped.(fn{j}){i,1}=dat.(fn{j})(idx(i)+1:idx(i+1));
        end
    end
end
dat_grouped=struct2table(dat_grouped);
clear idx dat
%% eliminate the edge data if there is this varargin input
w=1;
if ~isempty(edge)
    for i=1:height(dat_grouped)

        idx=find(dat_grouped.x{i,1}>edge(1,1)&dat_grouped.x{i,1}<edge(1,2)&dat_grouped.y{i,1}>edge(2,1)&dat_grouped.y{i,1}<edge(2,2));
        if isempty(idx)||length(idx)==1
            clear idx
            continue
        else
            % if larvae going back and forth in the edge, we need to divide it
            % into different trajectories
            ind=find(diff(idx)>1);
            ind=[0;ind;length(idx)];
            for b=1:length(ind)-1
                dat_edge.AN(w,1)=dat_grouped.AN(i,1);
                dat_edge.mid(w,1)=dat_grouped.mid(i,1);
                dat_edge.x{w,1}=dat_grouped.x{i,1}(idx(ind(b)+1):idx(ind(b+1)));
                dat_edge.y{w,1}=dat_grouped.y{i,1}(idx(ind(b)+1):idx(ind(b+1)));
                dat_edge.et{w,1}=dat_grouped.et{i,1}(idx(ind(b)+1):idx(ind(b+1)));
                w=w+1;
            end
    		clear idx
        end
    end
    clear dat_grouped;
    dat_grouped=struct2table(dat_edge);
    clear dat_edge;
end

%% move data to the (0,0)
ymin=min(vertcat(dat_grouped.y{:}));
xmin=min(vertcat(dat_grouped.x{:}));
for i=1:length(dat_grouped.x)
	dat_grouped.y{i,1}=dat_grouped.y{i,1}-ymin;
	dat_grouped.x{i,1}=dat_grouped.x{i,1}-xmin;
end 
%% add two new variables xcentered and ycentered which can be used for
%the centered trajectory later (all trajectory will have a start point at [0,0])
for i=1:length(dat_grouped.x)
    x=dat_grouped.x{i,1}(1,1);
    y=dat_grouped.y{i,1}(1,1);
    dat_grouped.xcentered{i,1}=dat_grouped.x{i,1}-x;
    dat_grouped.ycentered{i,1}=dat_grouped.y{i,1}-y;
end


%% calculate tsum
for i=1:length(dat_grouped.x)
    tsum(i,1)=dat_grouped.et{i,1}(end)-dat_grouped.et{i,1}(1);
end
dat_grouped.tsum=tsum;
clear tsum
%% eliminate tt less than 10s
idx=find(dat_grouped.tsum<10);
dat_grouped(idx,:)=[];
%% eliminate data have 99% within 3L*of its midline
for i=1:length(dat_grouped.x)
    x=dat_grouped.x{i,1}(1,1);
    y=dat_grouped.y{i,1}(1,1);
    L=dat_grouped.mid(i,1);
    l=length(dat_grouped.x{i,1});
    dis=sqrt((dat_grouped.x{i,1}-x).^2+(dat_grouped.y{i,1}-y).^2);
    ind=find(dis>=3*L);
    p_out_L(i,1)=length(ind)/l;
    clear x y ind dis L l dis

end
%% plot the prob vs frequency
%     p_bin=0:0.01:1;
%     for i=1:length(p_bin)-1
%         if i==length(p_bin)
%             idx=find(p_out_L>=p_bin(i)&p_out_L<=p_bin(i+1));
%             
%         else
%             idx=find(p_out_L>=p_bin(i)&p_out_L<p_bin(i+1));            
%         end
%         pro(i,1)=length(idx)/length(p_out_L);
% clear idx
%     end
%     prob=[0+0.01/2:0.01:1-0.01/2]';
%     mlt_subplt(pro,prob,100,1,1,1,append('Frequency of traj within 3*L circle@',genotype,'@',condition),'Probability of datapoint within 3L','Frequency',[0 0 0],append(genotype,'@',condition),'bar','ii_num',1);

idx=find(p_out_L<0.01);
dat_grouped(idx,:)=[];

%% plot and save the trajectory 
x=vertcat(dat_grouped.x{:});
y=vertcat(dat_grouped.y{:});
color=[0 0 0];
mlt_subplt(y,x,1,1,1,1,append('Trajectory with elimination,',genotype,'@',condition),'X(mm)','Y(mm)',color,append(genotype,'@',condition),'scatter');
save_all_figures(outdir)
close all;
clearvars -except dat_grouped outdir list odor_side

%save data

filename=fullfile(outdir,'data.mat');
if ~isfolder(outdir)
    mkdir(outdir)
end

if isfile(filename)
    delete(filename);
end
save(filename,'dat_grouped');
dat_height=height(dat_grouped)
clear
end

