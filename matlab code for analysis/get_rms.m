function []=get_rms(tracker_num,genotype,condition,filename)

%% get the directory for all the data
w=1;

for i=1:length(genotype)
    for j=1:length(condition)
        cond(w,1)=fullfile("/project/6010970/screen/olfactory_output/JB_JAABA",tracker_num,genotype{i,1},condition{j,1},filename);
        name(w,1)=append(genotype{i,1},"@",condition{j,1});
        w=w+1;
    end
end

ii_num=length(cond);

for ii=1:ii_num
    cond(ii)
    load(cond(ii));
    %% 1-1) Compute the turning frequency-->count of turning event/time
    % get the data ready for binning
    %     pre=vertcat(dat.pre_deg{:});
    %     post=vertcat(dat.post_deg{:});
    %     turn_x=vertcat(dat.turn_x{:});
    %     turn_y=vertcat(dat.turn_y{:});
    reori_deg=vertcat(dat.reorient_deg_abs{:});
    %     t0s1=vertcat(dat.t0s{:});t0s=vertcat(t0s1{:});
    
    rms_cond=sqrt((1/length(reori_deg))*sum(reori_deg.^2))
end
end