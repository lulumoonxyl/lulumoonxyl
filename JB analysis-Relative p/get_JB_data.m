function [dat1,fields1,uname]=get_JB_data(folder_name,ii)
% get the data from JB analysis

    choredir =folder_name;
    

        dat.xspine = {};
        dat.yspine ={};
        dat.AN={};
        dat.et={};
    
    
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
    
    delimiter = ' ';
    startRow = 0;
    formatSpec = '%s%f%f%f%[^\n\r]';
    idx = find(nb == ii);
    for jj = idx'
        dirname = fullfile(d(jj).folder,d(jj).name);
        fname = dir(fullfile(dirname,['*trx_reduced.mat']));
        f = contains({fname.name},['trx_reduced.']);
        fname = fullfile(fname(f).folder,fname(f).name);
        fileID = open(fname);
        dat.et= vertcat(dat.et, {fileID.trx.t}');
        dat.xspine = vertcat (dat.xspine, {fileID.trx.x_spine}');
        dat.yspine = vertcat (dat.yspine, {fileID.trx.y_spine}');
        dat.AN= vertcat(dat.AN, {fileID.trx.numero_larva_num}');
        clear fileID
        
    end

dat1= struct2cell(dat);
fields1=fieldnames(dat);
end