%% get jb data
ii1=1;
ii2=1;

choredir ='D:\jb-results\t88\CS@CS'
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
idx = find(nb == ii1);
dat.x = {};
dat.y ={};
dat.xspine = {};
dat.yspine ={};

dat.AN={};
dat.et={};

for jj = idx'
    dirname = fullfile(d(jj).folder,d(jj).name);
    fname = dir(fullfile(dirname,['*trx.mat']));
    f = contains({fname.name},['trx.']);
    
    fname = fullfile(fname(f).folder,fname(f).name);
    
    fileID = open(fname);
    dat.et= vertcat(dat.et, {fileID.trx.t}');
    dat.x = vertcat (dat.x, {fileID.trx.x_center}');
    dat.y = vertcat (dat.y, {fileID.trx.y_center}');
    dat.xspine = vertcat (dat.xspine, {fileID.trx.x_spine}');
    dat.yspine = vertcat (dat.yspine, {fileID.trx.y_spine}');
    dat.AN= vertcat(dat.AN, {fileID.trx.numero_larva_num}'); 
    clear fileID

end 

clearvars -except dat ii2
%% get data from the jaaba result
choredir ='D:\JAABA-results\t88\CS@CS'
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
idx = find(nb == ii2);
dat.tstart={};
dat.tend={};
dat.t0process={};
dat.t1process={};
AN_jaaba={}
et_jaaba={};

for jj = idx'
    dirname = fullfile(d(jj).folder,d(jj).name);
    fname = dir(fullfile(dirname,['*scores_Turn_ga_updated.mat']));
    f = contains({fname.name},['scores_Turn_ga_updated.']);
    fname = fullfile(fname(f).folder,fname(f).name);
    
    fileID = open(fname);
    dat.tstart= vertcat(dat.tstart, {fileID.allScores.allScores.tStart}');
    dat.tend = vertcat (dat.tend, {fileID.allScores.allScores.tEnd}');
    dat.t0process = vertcat (dat.t0process, {fileID.allScores.allScores.t0sProcessed}');
    dat.t1process = vertcat (dat.t1process, {fileID.allScores.allScores.t1sProcessed}');
    
    
    clear fileID fname f 
  
    fname = dir(fullfile(dirname,['*trx.mat']));
    f = contains({fname.name},['trx.']);
    fname = fullfile(fname(f).folder,fname(f).name);
    fileID = open(fname);
    AN_jaaba =vertcat (AN_jaaba, {fileID.trx.id}');
    et_jaaba =vertcat (et_jaaba, {fileID.trx.timestamps}');
    clear fileID fname f 
end 
clearvars -except AN_jaaba et_jaaba dat

%% calculate the relative probability of heading direction from the jb results
%use xspine and yspine for the heading direction 
dat.dir={};
for i=1:length(dat.xspine)
    x1=dat.xspine{i,1}(:,1);
    x2=dat.xspine{i,1}(:,7);
    x3=x1-x2;
    
    y1=dat.yspine{i,1}(:,1);
    y2=dat.yspine{i,1}(:,7);
    y3=y1-y2;
    
    degree=atan2d(y3.*1-0.*x3,x3.*1+y3.*0);
    dat.dir{i,1}=degree;
    clear degree x1 x2 x3 y1 y2 y3 
end 

% delete some data from jaaba 

dat.AN =vertcat(dat.AN{:});
AN_jaaba=vertcat(AN_jaaba{:});

AN={}; x={}; y={}; xspine={}; yspine={}; et={}; dir={};

ind=find(diff(dat.AN)<0);
for i=1:length(dat.tstart)
    if i==1
        AN{i,1}=dat.AN(1:ind(i));
        x{i,1}={dat.x{1:ind(i)}};
        y{i,1}={dat.y{1:ind(i)}};
        xspine{i,1}={dat.xspine{1:ind(i)}};
        yspine{i,1}={dat.yspine{1:ind(i)}};
        et{i,1}={dat.et{1:ind(i)}};
        dir{i,1}={dat.dir{1:ind(i)}};
    
    elseif i==length(dat.tstart)
        AN{i,1}=dat.AN(ind(i-1)+1:end);
        x{i,1}={dat.x{ind(i-1)+1:end}};
        y{i,1}={dat.y{ind(i-1)+1:end}};
        xspine{i,1}={dat.xspine{ind(i-1)+1:end}};
        yspine{i,1}={dat.yspine{ind(i-1)+1:end}};
        et{i,1}={dat.et{ind(i-1)+1:end}}; 
        dir{i,1}={dat.dir{ind(i-1)+1:end}};
    else 
        AN{i,1}=dat.AN(ind(i-1)+1:ind(i));
        x{i,1}={dat.x{ind(i-1)+1:ind(i)}};
        y{i,1}={dat.y{ind(i-1)+1:ind(i)}};
        xspine{i,1}={dat.xspine{ind(i-1)+1:ind(i)}};
        yspine{i,1}={dat.yspine{ind(i-1)+1:ind(i)}};
        et{i,1}={dat.et{ind(i-1)+1:ind(i)}}; 
        dir{i,1}={dat.dir{ind(i-1)+1:ind(i)}};
    end 
end 

dat.AN=AN; dat.x=x; dat.y=y; dat.xspine=xspine; dat.yspine=yspine; dat.et=et; dat.dir=dir;

ind1=1; AN_jaaba1={}; et_jaaba1={};
for i=1:length(dat.tstart)
    
    AN_jaaba1{i,1}=AN_jaaba(ind1:ind1+length(dat.tstart{i,1})-1);
    et_jaaba1{i,1}=et_jaaba(ind1:ind1+length(dat.tstart{i,1})-1);
    ind1=ind1+length(dat.tstart{i,1});

end 
clearvars -except dat AN_jaaba1 et_jaaba1


 % get the relative probability for each exp 
 dat.head_dir=[]; degree2=-180:20:180; dat.tsum=[]; dat.p=[];
 
 for i=1:length(dat.dir)
     tsum=0;
     for j=1:length(dat.dir{i,1})
         tsum=tsum+dat.et{i,1}{1,j}(end)-dat.et{i,1}{1,j}(1);
     end
     dat.tsum(1,i)=tsum;
 end
 
 for i=1:length(dat.dir)
     for l=1:length(degree2)-1
      t1=0;   
         for j=1:length(dat.dir{i,1})
             ind=find(dat.dir{i,1}{1,j}>=degree2(l)&dat.dir{i,1}{1,j}<degree2(l+1));
             
             ind1=find(ind==1);
             ind(ind1)=[];
             t=dat.et{i,1}{1,j}(ind)-dat.et{i,1}{1,j}(ind-1);
             t1=t1+sum(t);
             clear ind1 ind t
         end
         dat.head_dir(l,i)=t1;
         
     end
 end
 
 degree1=-170:20:170;
 
 dat.p=dat.head_dir./dat.tsum;
 dat.p_mean=mean(dat.p,2)
 dat.p_SEM=std(dat.p,0,2)./sqrt(length(dat.x));
 
 figure();
 hold on ;
 
 
 subplot(2,2,1)
 patch([degree1(1:3)';flipud(degree1(1:3)')],[dat.p_mean(1:3)-dat.p_SEM(1:3);flipud(dat.p_mean(1:3)+dat.p_SEM(1:3))],[1 0.41 0.16],'EdgeColor',[1 1 1])
 patch([degree1(3:7)';flipud(degree1(3:7)')],[dat.p_mean(3:7)-dat.p_SEM(3:7);flipud(dat.p_mean(3:7)+dat.p_SEM(3:7))],[0.6 0.7 0.8],'EdgeColor',[1 1 1])
 patch([degree1(7:12)';flipud(degree1(7:12)')],[dat.p_mean(7:12)-dat.p_SEM(7:12);flipud(dat.p_mean(7:12)+dat.p_SEM(7:12))],[0.8 0.7 0.8],'EdgeColor',[1 1 1])
 patch([degree1(12:16)';flipud(degree1(12:16)')],[dat.p_mean(12:16)-dat.p_SEM(12:16);flipud(dat.p_mean(12:16)+dat.p_SEM(12:16))],[0.6 1 0.8],'EdgeColor',[1 1 1])
 patch([degree1(16:18)';flipud(degree1(16:18)')],[dat.p_mean(16:18)-dat.p_SEM(16:18);flipud(dat.p_mean(16:18)+dat.p_SEM(16:18))],[1 0.41 0.16],'EdgeColor',[1 1 1])
 axis sqaure 
 
 hold on;
 plot(degree1',dat.p_mean,'r-','LineWidth',2);
 
 clearvars -except dat degree1 degree2 AN_jaaba1 et_jaaba1
 
 %% cancel some data in the jaaba results 
 for i=1:length(dat.AN)
     ind=find(ismember(AN_jaaba1{i,1},dat.AN{i,1})==0);
     for j=1:length(ind)
         dat.tstart{i,1}(1,ind(j))=nan;
         dat.tend{i,1}(1,ind(j))=nan;
         dat.t0process{i,1}{1,ind(j)}=nan;
         dat.t1process{i,1}{1,ind(j)}=nan;
     end 
     clear ind
 end 
 
 %% match the timepoints in jaaba and jb results 
 for i=1:length(dat.et)
     w=1;
     for j=1:length(dat.et{i,1})
         while isnan(dat.tstart{i,1}(1,w));
             w=w+1;
         end 
         ind1=find(ismember(et_jaaba1{i,1}{w,1},dat.et{i,1}{1,j})==0);
         ind2=find(diff(ind1)>1);
         ind3=find(dat.t0process{i,1}{1,w}<=ind2);
         dat.t0process{i,1}{1,w}(ind3)=[];
         dat.t1process{i,1}{1,w}(ind3)=[];
         
         dat.t0process{i,1}{1,w}=dat.t0process{i,1}{1,w}(:)-ind2;
         dat.t1process{i,1}{1,w}=dat.t1process{i,1}{1,w}(:)-ind2;
         
         for t=length(ind1):-1:1
             et_jaaba1{i,1}{w,1}(ind1(t))=[];
         end 
        w=w+1;
         
     end 
 end 
 
 %% delete and combine some turning event
 
dat.run0process={}; dat.run1process={};
for i=1:length(dat.x)
    w=1;
    for j=1:length(dat.x{i,1})
         while isnan(dat.tstart{i,1}(1,w))
             w=w+1;
         end 
        
        if isempty(dat.t0process{i,1}{1,w})
            dat.run0process{i,1}{1,j}=1;
            dat.run1process{i,1}{1,j}=length(et_jaaba1{i,1}{w,1});
            w=w+1;
            continue;
        end 
        
        %if it's smaller than 4 combines them
        if length(dat.t0process{i,1}{1,w})~=1
        turn_diff=dat.t0process{i,1}{1,w}(2:end)-dat.t1process{i,1}{1,w}(1:end-1);
        ind=find(turn_diff<4);
        dat.t0process{i,1}{1,w}(ind+1)=[];
        dat.t1process{i,1}{1,w}(ind)=[];
        end
        
        turn_diff_single=dat.t1process{i,1}{1,w}-dat.t0process{i,1}{1,w};
        ind1=find(turn_diff_single<3);
        dat.t0process{i,1}{1,w}(ind1)=[];
        dat.t1process{i,1}{1,w}(ind1)=[];
        
        l=length(et_jaaba1{i,1}{w,1});
        
        ind2=(find(dat.t0process{i,1}{1,w}>=l));
        dat.t1process{i,1}{1,w}(ind2)=[];
        dat.t0process{i,1}{1,w}(ind2)=[];
        
        ind3=(find(dat.t1process{i,1}{1,w}>l));
        dat.t1process{i,1}{1,w}(ind3)=l;
      
        clear ind1 ind2 turn_diff_single turn_diff ind3 

         if isempty(dat.t0process{i,1}{1,w})
            dat.run0process{i,1}{1,j}=1;
            dat.run1process{i,1}{1,j}=l;
            w=w+1;
            continue;
          end
        
        if dat.t0process{i,1}{1,w}(1)<=2 & dat.t1process{i,1}{1,w}(end)>= l-1
            dat.run0process{i,1}{1,j}=dat.t1process{i,1}{1,w}(1:end-1)+1;
            dat.run1process{i,1}{1,j}=dat.t0process{i,1}{1,w}(2:end)-1;
            
        elseif dat.t0process{i,1}{1,w}(1)<=2 & dat.t1process{i,1}{1,w}(end)< l-1
            dat.run0process{i,1}{1,j}=dat.t1process{i,1}{1,w}+1;
            dat.run1process{i,1}{1,j}=dat.t0process{i,1}{1,w}(2:end)-1;
            
            dat.run1process{i,1}{1,j}(end+1)=l;
            
        elseif dat.t0process{i,1}{1,w}(1)>2 & dat.t1process{i,1}{1,w}(end)>= l-1
            
            dat.run0process{i,1}{1,j}=dat.t1process{i,1}{1,w}(1:end-1)+1;
            dat.run1process{i,1}{1,j}=dat.t0process{i,1}{1,w}-1;
            
            dat.run0process{i,1}{1,j}=[1;dat.run0process{i,1}{1,j}];
        else 
            dat.run0process{i,1}{1,j}=dat.t1process{i,1}{1,w}+1;
            dat.run1process{i,1}{1,j}=dat.t0process{i,1}{1,w}-1;
            
            dat.run1process{i,1}{1,j}(end+1)=l;
            dat.run0process{i,1}{1,j}=[1;dat.run0process{i,1}{1,j}];
        end

        w=w+1;
        
    end 
    
end 
clearvars -except dat degree1 degree2 AN_jaaba1 et_jaaba1
%% plot one example

figure()
hold on;
plot(dat.x{i,1}{1,j},dat.y{i,1}{1,j},'g.','MarkerSize',7);
ind0=dat.t0process{i,1}{1,w};
ind1=dat.t1process{i,1}{1,w};
for l=1:length(ind0)
    ind2{l,1}(:,1)=ind0(l)+1:ind1(l)-1
end 
ind2=vertcat(ind2{:});
plot(dat.x{i,1}{1,j}(ind2),dat.y{i,1}{1,j}(ind2),'y.','MarkerSize',9);
plot(dat.x{i,1}{1,j}(ind0),dat.y{i,1}{1,j}(ind0),'r.','MarkerSize',9);
plot(dat.x{i,1}{1,j}(ind1),dat.y{i,1}{1,j}(ind1),'b.','MarkerSize',9);

xlabel('x(mm)');
ylabel('y(mm)');
legend('Running','Turning','Turning Starts','Turning Stops');
text(dat.x{i,1}{1,j}(1),dat.y{i,1}{1,j}(1),'Trajectory start point');
axis square;
clear ind0 ind1
%% get the probability of reorientation from running duration and running direction
dat.rundur={};
dat.rundir={};
dat.reorientation_rate={};
for i=1:length(dat.run0process)
    for j=1:length(dat.run0process{i,1})
        
        for z=1:length(dat.run0process{i,1}{1,j})
            start=dat.run0process{i,1}{1,j}(z);
            stop=dat.run1process{i,1}{1,j}(z);
            
            x1=dat.x{i,1}{1,j}(stop);
            x2=dat.x{i,1}{1,j}(start);
            x3=x1-x2;
            
            y1=dat.y{i,1}{1,j}(stop);
            y2=dat.y{i,1}{1,j}(start);
            y3=y1-y2;
            
            degree=atan2d(y3*1-0*x3,x3*1+y3*0);
            
            dat.rundir{i,1}{j,1}(z,1)=degree;
            
            dat.rundur{i,1}{j,1}(z,1)=dat.et{i,1}{1,j}(stop)-dat.et{i,1}{1,j}(start);
            
            dat.reorientation_rate{i,1}{j,1}(z,1)=1/dat.rundur{i,1}{j,1}(z,1);
            clear degree x1 x2 x3 y1 y2 y3 start stop
        end 
    end 
end 
clearvars -except dat degree1 degree2 AN_jaaba1 et_jaaba1
dat.rundir=vertcat(dat.rundir{:});
dat.rundir=vertcat(dat.rundir{:});

dat.reorientation_rate=vertcat(dat.reorientation_rate{:});
dat.reorientation_rate=vertcat(dat.reorientation_rate{:});

dat.reorientation_p=[];
dat.reorientation_sem=[];

    for l=1:length(degree2)-1
        ind=find(dat.rundir>=degree2(l)&dat.rundir<degree2(l+1));
        p=mean(dat.reorientation_rate(ind));
        sem=(std(dat.reorientation_rate(ind)))./sqrt(length(ind));
        dat.reorientation_p(l,1)=p;
        dat.reorientation_sem(l,1)=sem
    end 
clearvars -except dat degree1 degree2 AN_jaaba1 et_jaaba1

figure();
 hold on ;
subplot(2,2,1)
 patch([degree1(1:3)';flipud(degree1(1:3)')],[ dat.reorientation_p(1:3)-dat.reorientation_sem(1:3);flipud( dat.reorientation_p(1:3)+dat.reorientation_sem(1:3))],[1 0.41 0.16],'EdgeColor',[1 1 1])
 patch([degree1(3:7)';flipud(degree1(3:7)')],[ dat.reorientation_p(3:7)-dat.reorientation_sem(3:7);flipud( dat.reorientation_p(3:7)+dat.reorientation_sem(3:7))],[0.6 0.7 0.8],'EdgeColor',[1 1 1])
 patch([degree1(7:12)';flipud(degree1(7:12)')],[ dat.reorientation_p(7:12)-dat.reorientation_sem(7:12);flipud( dat.reorientation_p(7:12)+dat.reorientation_sem(7:12))],[0.8 0.7 0.8],'EdgeColor',[1 1 1])
 patch([degree1(12:16)';flipud(degree1(12:16)')],[ dat.reorientation_p(12:16)-dat.reorientation_sem(12:16);flipud( dat.reorientation_p(12:16)+dat.reorientation_sem(12:16))],[0.6 1 0.8],'EdgeColor',[1 1 1])
 patch([degree1(16:18)';flipud(degree1(16:18)')],[ dat.reorientation_p(16:18)-dat.reorientation_sem(16:18);flipud( dat.reorientation_p(16:18)+dat.reorientation_sem(16:18))],[1 0.41 0.16],'EdgeColor',[1 1 1])
 hold on;
 plot(degree1', dat.reorientation_p,'r-','LineWidth',2); 
 ylabel('Reorientation Rate');
 xlabel('Heading Direction (degrees)');
%% probability of orientation after turning
dat.turndir_before={};
dat.turndir_after={};

for i=1:length(dat.et)
    w=1;
    for j=1:length(dat.et{i,1})
         while isnan(dat.tstart{i,1}(1,w))
             w=w+1;
         end 
        
         if isempty(dat.t0process{i,1}{1,w})
             dat.turndir_before{i,1}{j,1}=[];
             dat.turndir_after{i,1}{j,1}=[];
             w=w+1;
             continue
         end 
             ind1=dat.t0process{i,1}{1,w};
             ind2=dat.t1process{i,1}{1,w};
             dat.turndir_before{i,1}{j,1}=dat.dir{i,1}{1,j}(ind1);
             dat.turndir_after{i,1}{j,1}=dat.dir{i,1}{1,j}(ind2);
        w=w+1;
        clear ind1 ind2
    end 
end 
for i=1:length(dat.turndir_before)
    dat.turndir_before{i,1}=vertcat(dat.turndir_before{i,1}{:});
    dat.turndir_after{i,1}=vertcat(dat.turndir_after{i,1}{:});
end
clearvars -except dat degree1 degree2 AN_jaaba1 et_jaaba1
for i=1:length(dat.turndir_before)
    ind=find(abs(dat.turndir_before{i,1})>45&abs(dat.turndir_before{i,1})<135);
    dat.turndir_before_135_45{i,1}=dat.turndir_before{i,1}(ind);
    dat.turndir_after_135_45{i,1}=dat.turndir_after{i,1}(ind);
end 
clearvars -except dat degree1 degree2 AN_jaaba1 et_jaaba1
dat.accept=[];
dat.non_accept=[];
for i=1:length(dat.turndir_before)
    ind=find(abs(dat.turndir_before_135_45{i,1})>abs(dat.turndir_after_135_45{i,1}));
    dat.accept(i,1)=length(ind)/length(dat.turndir_before_135_45{i,1});
    dat.non_accept(i,1)=(length(dat.turndir_before_135_45{i,1})-length(ind))/length(dat.turndir_before_135_45{i,1});
    clear ind;
    
end 

dat.accept_mean(1,1)=mean(dat.accept);
dat.accept_mean(1,2)=mean(dat.non_accept);
dat.accept_sem(1,1)=std(dat.accept)/sqrt(length(dat.accept));
dat.accept_sem(1,2)=std(dat.non_accept)/sqrt(length(dat.non_accept));

x=categorical ({'up gradient','down gradient'});
x=reordercats(x,{'up gradient','down gradient'});
figure();
hold on;
b=bar(x,dat.accept_mean);
errorbar(dat.accept_mean,dat.accept_sem,'.')

ylabel('Heading direction after turning');

