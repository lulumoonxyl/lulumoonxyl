
choredir='D:\jb-results\t88\CS@CS\replication'
outdir="D:\organized data qls600\jb_result";
%% reduce jb data
for ii=1:5
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
        fname = dir(fullfile(dirname,['*trx_reduced.mat']));
        f = contains({fname.name},['trx_reduced.']);
        fname = fullfile(fname(f).folder,fname(f).name);
        fileID = open(fname);
        trx=fileID.trx;
        f=fieldnames(trx);
        %ewmovw 50 fields from the trx file. 
        %fields_removed={f{1};f{3};{f{14:15}}';{f{26:46}}';{f{50:55}}';{f{57:58}}';{f{66:82}}'};
        %fields_removed={f(11:50)};
        fields_removed={{f{4:6}}';{f{12:13}}';f{21};f{5};{f{32:44}}'};
        fields_removed=vertcat(fields_removed{:});
        trx=rmfield(trx,fields_removed);
        trx_t=struct2table(trx);
        file_dir=fullfile(file_destination,d(jj).name);
        
       
        if ~isdir(file_dir)
            mkdir(file_dir);
        end
        outname='trx.mat';
        filename=fullfile(file_dir,outname);
        
        if isfile(filename)
            delete(filename);
        end
        
        if ~isfile(filename)
            save(filename,'trx_t','-v7.3');
        end
        
        clear fileID filename file_dir trx trx_t
      
    end
    clearvars -except outdir choredir ii
end
