function [dat1,fields1,uname]=get_data_from_dat_file(name_of_folder,ii)
% remember to check whether choreography results have the same length;
% if not try to run it in the student room
%this function will combine data with the same ii in different folder 

        dat.et = {};
        dat.x= {};
        dat.AN={};
        dat.y ={};
        dat.sp ={};
 
    choredir=name_of_folder;
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
    idx = find(nb == ii);
    % import chore data for first folder ?
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
        
        dat.AN=vertcat(dat.AN, dataArray{:,2}); %add animal num
        dat.et = vertcat(dat.et,dataArray{:,3});
        dat.sp =vertcat(dat.sp,dataArray_sp{:,4})
        dat.x = vertcat(dat.x,dataArray{:,4});
        dat.y = vertcat(dat.y,dataArray_y{:,4});
        clear dataArray dataArray_y dataArray_sp
    end
    clearvars -except ii dat uname i name_of_folder

dat1= struct2cell(dat);
fields1=fieldnames(dat);
end
