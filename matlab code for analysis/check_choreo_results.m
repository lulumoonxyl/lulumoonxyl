function []=check_choreo_results(tracker_num,genotype,condition)
choredir=fullfile("/project/6010970/screen/choreography-results",tracker_num,genotype,condition);
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
fn={"x","y","midline","speed085","speed","cast","curve","crabspeed","kink","bias","area","dir","morpwidth"};
%% import the data that are from the same condition
for jj = 1:length(nb)

    dirname = fullfile(d(jj).folder,d(jj).name);
    for k=1:length(fn)
        %open txt file and save the data
        fname= dir(fullfile(dirname,[append('*',(fn{k}),'.dat')]));
        f= contains({fname.name},[append('.',(fn{k}),'.')]);
        try
            fname= fullfile(fname(f).folder,fname(f).name);
        catch
            disp(append('empty ',(fn{k}),'.dat file in', dirname));
            continue
        end
        fileID = fopen(fname,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
        fclose(fileID);

        %check wheather the timestamps column of all dat files are the same
        timestamp=dataArray{1,1};
        idx{k,1}=find(strcmp(timestamp{1,1},timestamp)==0);
        len_file(1,k)=length(timestamp);
        if ~isempty(idx{k,1})
            disp(append('Abnormal animal number in ',(fn{k}),'.dat file,',num2str(length(idx{k,1})),'rows,' ,'In folder:',dirname));
        end
	if strcmp(fn{k},"x")
	xmax=max(dataArray{1,4})
	xmin=min(dataArray{1,4})
	elseif strcmp(fn{k},"y")
	ymax=max(dataArray{1,4})
	ymin=min(dataArray{1,4})
	end
        clear f fname fileID timestamp dataArray
    end

    an=length(unique(len_file));
    if an~=1
        disp(append('Abnormal length of dat file in folder:',dirname));
        for k=1:length(len_file)
            disp(append((fn{k}),'.dat file has ',num2str(len_file(1,k)), ' rows '));
        end
    else
        disp(append('All files (except outline and spine) have the same length for folder ',dirname));
    end
    clear dirname idx len_file an

end
end