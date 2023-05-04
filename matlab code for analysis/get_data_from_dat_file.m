function dat=get_data_from_dat_file(choredir)
%this function will get the 's','y', and 'midline' dat files from all the
%timestamps and group them into one table.

dat.et = {};
dat.x= {};
dat.AN={};
%dat.AN_y={};
%dat.AN_mid={};
dat.y ={};
dat.mid={};
%dat.sp ={};
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

%%get all the timestamp folder inside one condition
%Important: change '\' to '/' based on your system
names = {d.folder}';
splits = cellfun(@(x) split(x,'/'), names, 'UniformOutput', false);
%cellfun: apply a function to the array
splits = cellfun(@(x) [x{end-1},'/',x{end}], splits, 'UniformOutput', false);
[uname,na,nb] = unique(splits);
if length(uname)>1
    warning('The conditions of timestamp folders are not consistant');
    quit;
end
delimiter = ' ';
startRow = 0;
formatSpec = '%s%f%f%f%[^\n\r]';


%% import the data that are from the same condition
for jj = 1:length(nb)

    dirname = fullfile(d(jj).folder,d(jj).name);
    fname = dir(fullfile(dirname,['*x.dat']));%open txt file and save the data
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

    fname_mid = dir(fullfile(dirname,['*midline.dat']));
    f_mid = contains({fname_mid.name},['.midline.']);
    try
        fname_mid = fullfile(fname_mid(f_mid).folder,fname_mid(f_mid).name);
    catch
        fprintf(['empty folder ' dirname]);
        continue
    end

    %         fname_sp = dir(fullfile(dirname,['*speed.dat']));
    %         f_sp  = contains({fname_sp .name},['.speed.']);
    %         try
    %             fname_sp  = fullfile(fname_sp (f_sp ).folder,fname_sp (f_sp ).name);
    %         catch
    %             fprintf(['empty folder ' dirname]);
    %             continue
    %         end

    fileID = fopen(fname,'r');
    dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
    fclose(fileID);

    fileID_y = fopen(fname_y,'r');
    dataArray_y = textscan(fileID_y, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
    fclose(fileID_y);

    fileID_mid = fopen(fname_mid,'r');
    dataArray_mid = textscan(fileID_mid, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
    fclose(fileID_mid);

    %         fileID_sp = fopen(fname_sp,'r');
    %         dataArray_sp = textscan(fileID_sp, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
    %         fclose(fileID_sp);
    %
    dat.AN=vertcat(dat.AN, dataArray{:,2}); %add animal num
    dat.et = vertcat(dat.et,dataArray{:,3});
    %         dat.sp =vertcat(dat.sp,dataArray_sp{:,4});
    dat.x = vertcat(dat.x,dataArray{:,4});
    dat.y = vertcat(dat.y,dataArray_y{:,4});
    %dat.AN_y=vertcat(dat.AN_y, dataArray_y{:,2});
    dat.mid = vertcat(dat.mid,dataArray_mid{:,4});
    %dat.AN_mid=vertcat(dat.AN_mid, dataArray_mid{:,2});
    clear dataArray dataArray_y dataArray_mid dirname
end
clearvars -except ii dat uname i 

end
