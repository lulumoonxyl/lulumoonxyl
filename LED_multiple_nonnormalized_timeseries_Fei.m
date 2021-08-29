%% timeserires code - Adam Schwendt

% plot timeseries choreography result data (.dat files!)
% applicable to all choreograpgy data - just change the 'chore' value (could set this up as a function where 'chore' is an input variable)

%***YOU MUST MANUALLY EDIT THE FOLLOWING (IN ORDER)*****
    %SPECIFICATION OF CHOREOGRAPHY RESULTS
    %CHORE TYPE
    %idx LOOP DEPENDENT ON HOW MANY GENOTYPES TO PLOT TOGETHER
    %NORMALIZATION INTERVAL
    %TITLE
    %LEGEND 
%% Set parameters
for ii=1
% choredir = 'C:/Users/Unite 58/Desktop/AdamData/chore/choreography-results';
choredir = 'C:\Users\feihu\OneDrive\Desktop\CS@CS';
%choredir = 'F:\CS@CS' %uigetdir % specify the "choreography-results" (manual select)
% --> choose the folder with all of the genoytpes and their respective time stamps

% outdir = fullfile(pwd,'figures'); (Where do you want the figures to print?)
%outdir = fullfile('G:\Feidata\Summer 2021--Behavior algorithm to aversive odor GA and menthol\matlab'); 

%fileTypes = {'crabspeed','curve','kink','x','y','bias','speed','midline'};
%chore = {'x','y'}; % SET THE CHORE MANUALLY

%*****IMPORTANT*** Specify the TITLE BELOW; only the proper chore will be included


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
% 
filt = [full.isdir];
full = full(filt);
names = {full.name}';

expr = '^\d\d\d\d\d\d\d\d_\d\d\d\d\d\d';
filt = regexp(names,expr);
filt = cellfun(@(x) ~isempty(x), filt);
d = full(filt);

%% compare to figurelist (optional)
% I recommend you add something here so as not to remake all the figures
% every time

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

 %% -------------------------------------------------------------------%
% Get the x value 
%Previously you specified a folder that contains all of the genotypes for
%your choreography data --> now, this loop goes throught the folders in
%ALPHABETICAL ORDER! and calculates / NORMALIZES each series
  

 %declare ii for the loop 

    % determine files to import
    idx = find(nb == ii);
    % import chore data for first folder (100uW stimulation) 
    et = {};
    dat_x= {};
    animal_num={};
    dat_y ={};
    speed={};

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
        
        fileID = fopen(fname,'r');
        dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
        fclose(fileID);
        
       fileID_y = fopen(fname_y,'r');
        dataArray_y = textscan(fileID_y, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
        fclose(fileID_y);
%         
         animal_num=vertcat(animal_num, dataArray{:,2}); %add animal num
        et = vertcat(et,dataArray{:,3});
%         speed =vertcat(speed,dataArray{:,4}) 
        dat_x = vertcat(dat_x,dataArray{:,4});
        dat_y = vertcat(dat_y,dataArray_y{:,4});
        clear dataArray
        clear dataArray_y
    end
   

  result{ii,1}=animal_num; %store all the result in the same file 
    
   result{ii,2}=et; %store all the result in the same file 
    result{ii,3}=dat_x;
    result{ii,4}=dat_y;
    
    et = vertcat(et{:});
    dat_x = vertcat(dat_x{:});
    dat_y=vertcat(dat_y{:});
    speed=vertcat(speed{:});
    
  
    
    result{ii,6}=et; %store all the result in the same file 
    result{ii,7}=dat_x;
    result{ii,8}=dat_y;
    clearvars -except result
end  
    %plot map function
    time = 0:300:900; %plot based on per 5 min 
    figure(1);
    hold on;
    for i=1:length(time)-1;
        ind= find (result{ii,6}>=time(i)&result{ii,6}<time(i+1));
         subplot(3,3,i+6);
        plot(result{ii,7}(ind),result{ii,8}(ind),'.', 'MarkerSize', 1);
        title(['Trajectory of 10^-^3 GA, white, from t = ', num2str(time(i)),' to ', num2str(time(i+1)), 's']);
        xline(57,'r','LineWidth',1.5)
        xline(114,'r','LineWidth',1.5)
        xline(171,'r','LineWidth',1.5)
       
        
        axis tight
    end     
    hold off;
      
    figure();
    hold on;
    plot(dat_x,dat_y,'r.','MarkerSize',1);
     xlabel('x(mm)');
     ylabel('y(mm)');
    title('White, 10^-^5 EA');
    
%   for i=1:10  
%   result{i,10}=animal_num{i,1};  
%   end    
%   
%     for i=1: length(dat)
%     if (dat(i)<47)
%         dat(i)=-1;
%     elseif(dat(i)>57)
%         dat(i)=1;
%     else
%         dat(i)=0;
%     end 
%     
%     end  

    % calculate timeseries %
    bins = 0:0.5:ceil(max(et));
    nanarr = nan(1,length(bins));
    Y = discretize(et,bins);
    seri = accumarray(Y,speed,[],@mean);
    nanarr(1:length(seri)) = seri;
    seri = nanarr;

    figure(1);
    hold on;
    subplot(2,1,1)
    hold on;
    plot(bins, seri, '-');
    plot(bins2, seri2, '-');
%     plot(bins3, seri3, '-');
    xlabel('time (sec)');
    ylabel ('speed')
    legend('2*10^-4 GA','5*10^-^3 GA','5*10^-^3 Oil','10^-^2 GA');
    %% -------------------------------------------------------------------%
    
    idx2 = find(nb == 6);
    % import chore data for the second folder (200uW stimulation)
    et2 = {};
    dat2_x = {};
    dat2_y={};
    animal_num2={};
  speed2 = {};
    for jj = idx'
        dirname2 = fullfile(d(jj).folder,d(jj).name);
        fname2 = dir(fullfile(dirname2,['*speed.dat']));
        f = contains({fname2.name},['.speed.']);
        try
            fname2 = fullfile(fname2(f).folder,fname2(f).name);
        catch
            fprintf(['empty folder ' dirname2]);
            continue
        end
        
        
%         fname2_y = dir(fullfile(dirname2,['*y.dat']));
%         f_y = contains({fname2_y.name},['.y.']);
%         try
%             fname2_y = fullfile(fname2_y(f_y).folder,fname2_y(f_y).name);
%         catch
%             fprintf(['empty folder ' dirname2]);
%             continue
%         end
        
        fileID2 = fopen(fname2,'r');
        dataArray2 = textscan(fileID2, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
        fclose(fileID2);
        
        
%         fileID2_y = fopen(fname2_y,'r');
%         dataArray2_y = textscan(fileID2_y, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
%         fclose(fileID2_y);
%         
%         animal_num2 = vertcat(animal_num2, dataArray2{:,2});
 et2={};    
et2 = vertcat(et2,dataArray2{:,3});
%         dat2_x = vertcat(dat2_x,dataArray2{:,4});
%         dat2_y = vertcat(dat2_y,dataArray2_y{:,4});
       
speed2= vertcat(speed2,dataArray2{:,4});
        clear dataArray2
%         clear dataArray2_y
    end
    result{2,1}=animal_num
    result{2,2}=et
    result{2,3}=dat_x
    result{2,4}=dat_y
    
    et2 = vertcat(et2{:});
    dat2_x = vertcat(dat2_x{:});
    dat2_y=vertcat(dat2_y{:});
    speed2 =  vertcat (speed2 {:});
    
    result{2,6}=et; %store all the result in the same file
    result{2,7}=dat_x;
    result{2,8}=dat_y;
    
   %plot map function
   %plot based on per 5 min 
  
   figure(1)
   hold on;
    for i=1:length(time)-1;
        ind= find (result{2,6}>=time(i)&result{2,6}<time(i+1));
        subplot(2,3,i+3);
        plot(result{2,7}(ind),result{2,8}(ind),'.', 'MarkerSize', 1);
        title(['Trajectory of H2O from t = ', num2str(time(i)),' to ', num2str(time(i+1)), 's']);
        xline(57,'r','LineWidth',1.5)
        xline(114,'r','LineWidth',1.5)
        xline(171,'r','LineWidth',1.5)
       
        axis tight
    end     
      
%   for i=1:10
%   result{i+10,10}=animal_num2{i,1};  
%   end    
%   
%     
%     for i=1: length(dat2)
%     if (dat2(i)<47)
%         dat2(i)=-1;
%     elseif(dat2(i)>57)
%         dat2(i)=1;
%     else 
%         dat2(i)=0;
%         
%     end 
%     
%     end  

   % calculate timeseries %
    bins2 = 0:0.5:ceil(max(et2));
    nanarr2 = nan(1,length(bins2));
    Y2 = discretize(et2,bins2);
    seri2 = accumarray(Y2,speed2,[],@mean);
    nanarr2(1:length(seri2)) = seri2;
    seri2 = nanarr2;
    
    %-------------------------------------------------------------------%
    
     idx3 = find(nb == ii+2);
     %     % import chore data for the third folder (300uW stimulation)
     et3 = {};
%      dat3_y = {};
%      dat3_x={};
%      animal_num3={};
     speed3 ={};
     for jj = idx3'
         dirname3 = fullfile(d(jj).folder,d(jj).name);
         fname3 = dir(fullfile(dirname3,['*speed.dat']));
         f = contains({fname3.name},['.speed.']);
         try
             fname3 = fullfile(fname3(f).folder,fname3(f).name);
         catch
             fprintf(['empty folder ' dirname3]);
             continue
         end
         
%            fname3_y = dir(fullfile(dirname3,['*y.dat']));
%          f_y = contains({fname3_y.name},['.y.']);
%          try
%              fname3_y = fullfile(fname3_y(f_y).folder,fname3_y(f_y).name);
%          catch
%              fprintf(['empty folder ' dirname3]);
%              continue
%          end
%          
         
         fileID3 = fopen(fname3,'r');
         dataArray3 = textscan(fileID3, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
         fclose(fileID3);
         
%           fileID3_y = fopen(fname3_y,'r');
%          dataArray3_y = textscan(fileID3_y, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
%          fclose(fileID3_y);
         
         et3 = vertcat(et3,dataArray3{:,3});
%          animal_num3 =vertcat(animal_num3,dataArray3{:,2});
%          dat3_x = vertcat(dat3_x,dataArray3{:,4});
%          dat3_y = vertcat(dat3_y,dataArray3_y{:,4});
speed3 = vertcat(speed3,dataArray3{:,4});
         clear dataArray3
%          clear dataArray3_y
     end
     
     result{3,1}=animal_num3;
     result{3,2}=et3;
     result{3,3}=dat3_x;
     result{3,4}=dat3_y;
     
     et3 = vertcat(et3{:});
%      dat3_x = vertcat(dat3_x{:});
%      dat3_y = vertcat(dat3_y{:});
     
     speed3 = vertcat (speed3{:});
     
     result{3,6}=et3;
     result{3,7}=dat3_x;
     result{3,8}=dat3_y;
     
       
   figure(1)
   hold on;
    for i=1:length(time)-1;
        ind= find (result{3,6}>=time(i)&result{3,6}<time(i+1));
        subplot(3,3,i+6);
        plot(result{3,7}(ind),result{3,8}(ind),'.', 'MarkerSize', 1);
        title(['Trajectory of 5*10^-^3 Oil from t = ', num2str(time(i)),' to ', num2str(time(i+1)), 's']);
        xline(57,'r','LineWidth',1.5)
        xline(114,'r','LineWidth',1.5)
        xline(171,'r','LineWidth',1.5)
       
        axis tight
    end     
      
     
%      % calculate timeseries %
    bins3 = 0:0.5:ceil(max(et3));
    nanarr3 = nan(1,length(bins3));
    Y3 = discretize(et3,bins3);
    seri3 = accumarray(Y3,speed3,[],@mean);
    nanarr3(1:length(seri3)) = seri3;
    seri3 = nanarr3;
%     
%     %-------------------------------------------------------------------%
%     
%     idx4 = find(nb == ii+3);
%     % import chore data for the fourth folder (400uW stimulation)
%     et4 = {};
%     dat4 = {};
%     for jj = idx4'
%         dirname4 = fullfile(d(jj).folder,d(jj).name);
%         fname4 = dir(fullfile(dirname4,['*' chore{:} '.dat']));
%         f = contains({fname4.name},['.' chore{:} '.']);
%         try
%             fname4 = fullfile(fname4(f).folder,fname4(f).name);
%         catch
%             fprintf(['empty folder ' dirname4]);
%             continue
%         end
%         
%         fileID4 = fopen(fname4,'r');
%         dataArray4 = textscan(fileID4, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
%         fclose(fileID4);
% 
%         et4 = vertcat(et4,dataArray4{:,3});
%         dat4 = vertcat(dat4,dataArray4{:,4});
%         clear dataArray4
%     end
%     et4 = vertcat(et4{:});
%     dat4 = vertcat(dat4{:});
%     
%      % calculate timeseries %
%     bins4 = 0:0.5:ceil(max(et4));
%     nanarr4 = nan(1,length(bins4));
%     Y4 = discretize(et4,bins4);
%     seri4 = accumarray(Y4,dat4,[],@mean);
%     nanarr4(1:length(seri4)) = seri4;
%     seri4 = nanarr4;
    
     %-------------------------------------------------------------------%
%      
%     idx5 = find(nb == ii+4); 
%     %import chore data for fifth folder (No stimulation, Control) 
%     et5 = {};
%     dat5 = {};
%    for jj = idx5'
%       dirname5 = fullfile(d(jj).folder,d(jj).name);
%        fname5 = dir(fullfile(dirname5,['*' chore{:} '.dat']));
%        f = contains({fname5.name},['.' chore{:} '.']);
%        try
%            fname5 = fullfile(fname5(f).folder,fname5(f).name);
%         catch
%            fprintf(['empty folder ' dirname5]);
%            continue
%         end
%         
%        fileID5 = fopen(fname5,'r');
%        dataArray5 = textscan(fileID5, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
%        fclose(fileID5);
% 
%        et5 = vertcat(et5,dataArray5{:,3});
%        dat5 = vertcat(dat5,dataArray5{:,4});
%         clear dataArray5
%     end
%     et5 = vertcat(et5{:});
%     dat5 = vertcat(dat5{:});
%     
%      % calculate timeseries %
%     bins5 = 0:0.5:ceil(max(et5));
%     nanarr5 = nan(1,length(bins5));
%     Y5 = discretize(et5,bins5);
%     seri5 = accumarray(Y5,dat5,[],@mean);
%     nanarr5(1:length(seri5)) = seri5;
%     seri5 = nanarr5;
%   
%   %-------------------------------------------------------------------%
%     idx6 = find(nb == ii+5); 
%     % import chore data for the second folder (600uW stimulation) 
%     et6 = {};
%     dat6 = {};
%     for jj = idx6'
%         dirname6 = fullfile(d(jj).folder,d(jj).name);
%         fname6 = dir(fullfile(dirname6,['*' chore{:} '.dat']));
%         f = contains({fname6.name},['.' chore{:} '.']);
%         try
%             fname6 = fullfile(fname6(f).folder,fname6(f).name);
%         catch
%             fprintf(['empty folder ' dirname6]);
%             continue
%         end
%         
%         fileID6 = fopen(fname6,'r');
%         dataArray6 = textscan(fileID6, formatSpec, 'Delimiter', delimiter,'MultipleDelimsAsOne', true, 'HeaderLines' ,startRow, 'ReturnOnError', false);
%         fclose(fileID6);
%  
%         et6 = vertcat(et6,dataArray6{:,3});
%         dat6 = vertcat(dat6,dataArray6{:,4});
%         clear dataArray6
%     end
%     et6 = vertcat(et6{:});
%     dat6 = vertcat(dat6{:});
%     
%      % calculate timeseries %
%     bins6 = 0:0.5:ceil(max(et6));
%     nanarr6 = nan(1,length(bins6));
%     Y6 = discretize(et6,bins6);
%     seri6 = accumarray(Y6,dat6,[],@mean);
%     nanarr6(1:length(seri6)) = seri6;
%     seri6 = nanarr6;

  %-------------------------------------------------------------------%
    
    % PLOT TIMESERIES %
    %colour1 = [0, 0.4470, 0.7410];
    %colour2 = [0.8500, 0.3250, 0.0980];
    %colour3 = [0.9290, 0.6940, 0.1250];
    %colour4 = [0.4660, 0.6740, 0.1880];
    %colour5 = [0.6350, 0.0780, 0.1840];
    
    %Insert a rectangle for LED or audio stimulation
    %rectangle('Position', [44.7 0.1 30 1.0], 'FaceColor', [0.92, 0.92, 0.92], 'EdgeColor', [0.92, 0.92, 0.92]);
    %hold on
    
    p = plot(result{3,15},result{3,16}, 'Color', [0, 0, 0], 'Linewidth', 1 ); %First Folder - 100 uW
    hold on 
    p2 = plot(result{4,15},result{4,16}, 'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 1); %Second folder - 200uW
    hold on
%     p3 = plot(bins3,seri3, 'Color', [0.9290, 0.6940, 0.1250], 'Linewidth', 1); %Third folder - 300uW
%     hold on
%     p4 = plot(bins4,seri4, 'Color', [0, 0.4470, 0.7410], 'Linewidth', 1,'LineStyle','--'); %Fourth folder - 400uW
%     hold on
%     p3 = plot(bins3,seri3, 'Color', [0.6350, 0.0780, 0.1840], 'Linewidth', 1); %Fifth folder - Control
%     hold on
%     p4 = plot(bins4,seri4, 'Color', [0, 0, 0], 'Linewidth', 1); %Fifth folder - Control
%     hold on
    %SET THE X LIMIT *****
    xlim([0,600]);
    %ylim([2.5,4]);
    
    
   
  
    
 %-------------------------------------------------------------------%
 
    % AXIS LABELS %

 fileName=strrep(uname{ii},'/',' ');
    ax = gca;
    ax.YLabel.String =  ['X position (mm)']; %variable
    ax.XLabel.String = 'Time (seconds)'; 
    ch = chore; 
    ax.Title.String = ['X position']; %**MANUALLY EDIT** ' Audio Stimulation' or LED Stimulation and ?CHORE
 
    % Legend with sequential names, make sure it corresponds with the idx
    % loop!
    hleg = legend('42a-C','42a-X'); %Audio
    htitle= get(hleg,'Title');
    set(htitle,'String', 'Genotypes');
    
    %legend(dirname, dirname2, dirname3, dirname4, dirname5);
    
%-------------------------------------------------------------------%
    
    % SAVE TO OUTDIR %
    
    if ~isdir(outdir)
        mkdir(outdir);
    end
   
    outname = 'x_42a'; %% SET THIS YOURSELF
    
    %saved in matlab
    filepath = fullfile(outdir,outname);
    
    %orient the paper so that it is landscape
    orient landscape
    
    %saved in pdf in landscape
    print(filepath,'-dpdf','-painters','-fillpage');
    
    %closes figure
    hold off
    


