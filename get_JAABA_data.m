function [dat1,fields1,uname]=get_JAABA_data(folder_name,ii)
choredir=folder_name;
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

names = {d.folder}';
splits = cellfun(@(x) split(x,'\'), names, 'UniformOutput', false);
%cellfun: apply a function to the array
splits = cellfun(@(x) [x{end-1},'\',x{end}], splits, 'UniformOutput', false);
[uname,na,nb] = unique(splits);

delimiter = ' ';
startRow = 0;
formatSpec = '%s%f%f%f%[^\n\r]';
idx = find(nb == ii);
% dat.tstart={};
% dat.tend={};
dat.t0process={};
dat.t1process={};
dat.AN={};
dat.et={};

for jj = idx'
    dirname = fullfile(d(jj).folder,d(jj).name);
    fname = dir(fullfile(dirname,['*scores_Turn_ga_updated.mat']));
    f = contains({fname.name},['scores_Turn_ga_updated.']);
    fname = fullfile(fname(f).folder,fname(f).name);
    
    fileID = open(fname);
%     dat.tstart= vertcat(dat.tstart, {fileID.allScores.allScores.tStart}');
%     dat.tend = vertcat (dat.tend, {fileID.allScores.allScores.tEnd}');
    dat.t0process = vertcat (dat.t0process, {fileID.allScores.allScores.t0sProcessed}');
    dat.t1process = vertcat (dat.t1process, {fileID.allScores.allScores.t1sProcessed}');
    
    
    clear fileID fname f
    
    fname = dir(fullfile(dirname,['*trx.mat']));
    f = contains({fname.name},['trx.']);
    fname = fullfile(fname(f).folder,fname(f).name);
    fileID = open(fname);
    dat.AN =vertcat (dat.AN, {fileID.trx.id}');
    dat.et =vertcat (dat.et, {fileID.trx.timestamps}');
    clear fileID fname f
end
dat1= struct2cell(dat);
fields1=fieldnames(dat);

end