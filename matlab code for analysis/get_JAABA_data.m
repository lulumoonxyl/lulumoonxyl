function [dat,uname]=get_JAABA_data(choredir)

t0s={};
t1s={};
AN={};

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

filt = [full.isdir];
full = full(filt);
names = {full.name}';

expr = '^\d\d\d\d\d\d\d\d_\d\d\d\d\d\d';
filt = regexp(names,expr);
filt = cellfun(@(x) ~isempty(x), filt);
d = full(filt);

names = {d.folder}';
splits = cellfun(@(x) split(x,'\'), names, 'UniformOutput', false);
%cellfun: apply a function to the array
splits = cellfun(@(x) [x{end-1},'\',x{end}], splits, 'UniformOutput', false);
[uname,na,nb] = unique(splits);

if length(uname)>1
    warning('The conditions of timestamp folders are not consistant');
    quit;
end


for jj = 1:length(nb)
    dirname = fullfile(d(jj).folder,d(jj).name);
    fname = dir(fullfile(dirname,['*scores_Turn_ga_updated.mat']));
    f = contains({fname.name},['scores_Turn_ga_updated.']);

    if length(f)==1
        fname = fullfile(fname(f).folder,fname(f).name);
        fileID = open(fname);
        t0s = vertcat (t0s, {fileID.allScores.allScores.t0sSeconds'});
        t1s = vertcat (t1s, {fileID.allScores.allScores.t1sSeconds'});
    else
        disp('This folder does not contain the score_Turn_ga_updated.mat. Therefore, use score_Turn_ga.mat instead');
        fname = dir(fullfile(dirname,['*scores_Turn_ga.mat']));
        f = contains({fname.name},['scores_Turn_ga.']);

        fname_et = dir(fullfile(dirname,['perframe\timestamps.mat']));
        f_et = contains({fname_et.name},['timestamps']);
        if length(f)==1
        fname = fullfile(fname(f).folder,fname(f).name);
        fileID = open(fname);

        fname_et=fullfile(fname_et(f_et).folder,fname_et(f_et).name);
        fileID_et=open(fname_et);
        % sort the time 
        et=cell2mat(fileID_et.data);
        et=unique(et);
        et=sort(et);
        % get the time in second where turning happens
        % this does not modify the turning event
        % short event and events occur right after each other will not be
        % fixed
        % This issue is fixed when we run the group_JB_JAABA_data code
        t0_idx={fileID.allScores.t0s'}; t0_idx=vertcat(t0_idx{:});
        t1_idx={fileID.allScores.t1s'}; t1_idx=vertcat(t1_idx{:});

        t0sSeconds=cell(length(t0_idx),1);
        t1sSeconds=cell(length(t1_idx),1);

        for k=1:length(t0_idx)
           t0sSeconds{k,1}=et(t0_idx{k,1});
           t1sSeconds{k,1}=et(t1_idx{k,1});
        end 
        t0s = vertcat (t0s, {fileID.allScores.allScores.t0sSeconds'});
        t1s = vertcat (t1s, {fileID.allScores.allScores.t1sSeconds'});

        else
            t=vertcat(t0_idx{:});
            disp("This folder does not have any JAABA score! Reprocess and double check the data!")
            quit
        end 
    end

    fname_trx = dir(fullfile(dirname,['*trx.mat']));
    f_trx = contains({fname_trx.name},['trx.']);

    try
        fname_trx = fullfile(fname_trx(f).folder,fname_trx(f).name);
    catch
        fprintf(['empty folder ' dirname]);
        continue
    end

    fileID_trx = open(fname_trx);
    AN =vertcat (AN, {fileID_trx.trx.id}');
    clear fileID fileID_trx fname fname_trx f_trx f
end
t0s=vertcat(t0s{:});t1s=vertcat(t1s{:});AN=vertcat(AN{:});
dat=table(t0s,t1s,AN);
clearvars -except uname dat
end