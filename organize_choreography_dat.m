%% renaming choreography data and move them to some other folders if we want
%basically delete '_week1/2/3' in the .dat name

%choredir='D:\choreography-result\CS@CS_replications\*\';
%outdir='C:\Users\feihu\OneDrive - McGill University\matlab\demo_data\choreography-result';
%expression='_week+[1-3]_'-->the part of naming that you want to replace 
%replace='_';
function cond=organize_choreography_dat(choredir,outdir,cond,expression,replace)
for ii=1:cond
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
    %this gets the naming of the file, for example splits='o_five_10n1GA_week1_0s1x900s0s#n#n#n@40'
    
    splits = cellfun(@(x) [x{end}], splits, 'UniformOutput', false);
    %this will delete week1/2/3 from the naming so that we can group them
    %together
    week=["_week1","_week2","_week3"];
    splits=erase(splits,week);
    [uname,na,nb] = unique(splits);
    
    delimiter = ' ';
    startRow = 0;
    formatSpec = '%s%f%f%f%[^\n\r]';
    idx = find(nb == ii);
  
    file_destination=fullfile(outdir,uname(ii));
    
    
    for jj = idx'
        
        dirname = fullfile(d(jj).folder,d(jj).name);
        fname = dir(fullfile(dirname,['*.dat']));
        oldnames={fname.name};
        newnames=regexprep(oldnames,expression,replace)
        
        file_dir=fullfile(file_destination,d(jj).name);
        file_dir=cell2mat(file_dir);
        if ~isdir(file_dir)
            mkdir(file_dir);
        end
        
        for k=1:length(oldnames)
            if ~strcmp(oldnames{k},newnames{k})
                movefile(fullfile(dirname,oldnames{k}),fullfile(dirname,newnames{k}))
            end
            
            if ~isempty(outdir)
            file_name_new=fullfile(dirname,newnames{k});
            copyfile(file_name_new, file_dir)
            end 
        end
        
        clear oldnames newnames dirname fname k file_dir
    end
    clearvars -except ii choredir outdir
end
end