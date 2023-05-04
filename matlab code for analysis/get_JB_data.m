function [dat,uname]=get_JB_data(choredir)
% get the data from JB analysis
%create empty strcture where we store the data

xspine = {};
yspine ={};
AN={};
et={};

%get all the subfolder name
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
%cellfun: apply a function to the array to get the last section of the name of the file
splits = cellfun(@(x) [x{end-1},'\',x{end}], splits, 'UniformOutput', false);
[uname,na,nb] = unique(splits);%find how many conditions we have based on the filename

if length(uname)>1
    warning('The conditions of timestamp folders are not consistant');
    quit;
end


for jj = 1:length(nb)
    dirname = fullfile(d(jj).folder,d(jj).name);
    fname = dir(fullfile(dirname,['*trx.mat']));%open the trx file in each folder
    f = contains({fname.name},['trx.']);
    try
        fname = fullfile(fname(f).folder,fname(f).name);
    catch
        fprintf(['empty folder ' dirname]);
        continue
    end

    fileID = open(fname);
    et= vertcat(et, {fileID.trx.t}');
    xspine = vertcat (xspine, {fileID.trx.x_spine}');%save the specific data into each field of the structure
    yspine = vertcat (yspine, {fileID.trx.y_spine}');
    AN= vertcat(AN, {fileID.trx.numero_larva_num}');
    clear fileID

end
AN=vertcat(AN{:});
dat=table(xspine,yspine,AN,et);
end