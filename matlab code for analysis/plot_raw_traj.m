function []=plot_raw_traj(tracker_num,genotype,condition,odor_side)
%% plot the raw trajectory of all timestamps and decide whether we want to eliminate the edge
%now only set odor side to be'left' and do not consider other situations 
%if odor_side=='left', we will flip the x postion to make the odor in the
%right side 
%% identify the input and output folder
choredir=fullfile("/project/6010970/screen/choreography-results",tracker_num,genotype,condition);
outdir=fullfile("/project/6010970/screen/olfactory_output/choreography",tracker_num,genotype,condition);
%% load data
%get raw choreography data (x,y,midline) from the .txt files
dat=get_data_from_dat_file(choredir);
%% concatenate all the data

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
%% plot the traj
color=[0 0 0];
mlt_subplt(dat.y,dat.x,1,1,1,1,append('Raw Trajectory@',genotype,'@',condition),'X(mm)','Y(mm)',color,append(genotype,'@',condition),'scatter');
save_all_figures(outdir)
close all;
end 